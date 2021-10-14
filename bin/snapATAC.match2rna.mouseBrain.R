#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input", required=TRUE, help="load snap RData")
parser$add_argument("-s", "--seurat", required=TRUE, help="load seurat RData")
#parser$add_argument("--cpu", default = 4, help="# of cpus [default %(default)s]")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

suppressPackageStartupMessages(library("SnapATAC"))
suppressPackageStartupMessages(library("Seurat"))
suppressPackageStartupMessages(library("GenomicRanges"));
suppressPackageStartupMessages(library("RColorBrewer"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("plyr"))
#library("future")
library("Matrix")
library("tictoc")

snapF = args$input
seuratF = args$seurat
#cpus = as.numeric(args$cpu)
outF = args$output

#----------------------
# load rna seurat
rna.se <- readRDS(seuratF)

#-------------------
# downsample to 200
idx.ls <- list()
for(l in unique(rna.se$ClusterName)){
  idx <- which(rna.se$ClusterName == l)
  if(length(idx)>200){
    set.seed(2020)
    idx.dn <- sample(idx, 200)
  }else{
    print(length(idx))
    set.seed(2020)
    idx.dn <- idx
  }
  idx.ls[[l]] <- idx.dn
}
idx.ls <- sort(unlist(idx.ls))

rna.se <- rna.se[, idx.ls]

rna.se <- NormalizeData(rna.se)
rna.se <- FindVariableFeatures(rna.se)
all.genes <- rownames(rna.se)
rna.se <- ScaleData(rna.se, features = all.genes)
rna.se <- RunPCA(rna.se)
#ElbowPlot(rna.subset.se)
rna.se <- FindNeighbors(rna.se, dims = 1:20)
rna.se <- FindClusters(rna.se, resolution = 0.5)
rna.se <- RunUMAP(rna.se, dims = 1:20)

#DimPlot(rna.subset.se, reduction = "umap")
#DimPlot(rna.subset.se, reduction = "umap", group.by = "ClusterName")
fwrite(rna.se[[]], paste(outF, "rna.seObj.meta.tsv", sep="."), sep="\t", quote=F, col.names=T, row.names=F)
saveRDS(rna.se, file=paste(outF, "rna.seObj.rds", sep="."))


#-------------------------
# load atac and convert to Seurat Obj
x.sp <- readRDS(snapF)
x.sp@metaData$cluster <- x.sp@metaData$L3cluster

#-------------------
# downsample to 200
idx.ls <- list()
for(l in unique(x.sp@metaData$cluster)){
  idx <- which(x.sp@metaData$cluster == l)
  if(length(idx)>200){
    set.seed(2020)
    idx.dn <- sample(idx, 200)
  }else{
    print(length(idx))
    set.seed(2020)
    idx.dn <- idx
  }
  idx.ls[[l]] <- idx.dn
}
idx.ls <- sort(unlist(idx.ls))

x.sp <- x.sp[idx.ls, ]


pca.use = x.sp@smat@dmat;
eigs.dims <- 1:ncol(pca.use)
pca.use = pca.use[,eigs.dims]

metaData.use = x.sp@metaData
data.use = t(x.sp@pmat)
peak.use = as.data.frame(x.sp@peak);
gmat.use = t(x.sp@gmat);

rownames(x = data.use) = peak.use$name;
colnames(x = data.use) = paste(x.sp@sample, x.sp@barcode, sep=".");
colnames(x = gmat.use) = paste(x.sp@sample, x.sp@barcode, sep=".");
rownames(x = pca.use)  = paste(x.sp@sample, x.sp@barcode, sep=".");
rownames(metaData.use) = paste(x.sp@sample, x.sp@barcode, sep=".");
colnames(x = pca.use) <- paste0("pca_", eigs.dims);

atac.se <- CreateSeuratObject(counts = data.use, assay = "ATAC");
atac.se[["ACTIVITY"]] <- CreateAssayObject(counts = gmat.use);
atac.se <- AddMetaData(atac.se, metadata = metaData.use);
atac.se$tech <- "atac"
DefaultAssay(atac.se) <- "ACTIVITY"

atac.se[["pca"]] <- new(Class = "DimReduc", cell.embeddings = pca.use,
           feature.loadings = matrix(0,0,0), feature.loadings.projected = matrix(0,0,0),
           assay.used ="ATAC", stdev = rep(1,length(eigs.dims)), 
           key ="pca_", jackstraw = new(Class = "JackStrawData"), misc = list()) 

atac.se <- NormalizeData(atac.se)
atac.se <- FindVariableFeatures(object = atac.se, assay="ACTIVITY");

all.genes <- rownames(atac.se)
atac.se <- ScaleData(atac.se, features = all.genes)      
#atac.se <- RunPCA(atac.se, verbose = FALSE)
atac.se <- RunUMAP(atac.se, dims = 1:20)
#atac.se <- FindNeighbors(atac.se, dims = 1:20)
#atac.se <- FindClusters(atac.se, resolution = 0.5)

#DimPlot(atac.se, reduction = "umap")
#DimPlot(atac.se, reduction = "umap", group.by="L3cluster")


# set parallel
#plan("multiprocess", workers = cpus)

#---------------------------
# identify variable genes from RNA
variable.genes = VariableFeatures(object = rna.se);

#--------------------------
# transfer anchors
rna2atac.transfer.anchors <- FindTransferAnchors(
    reference = rna.se, 
    query = atac.se, 
    features = variable.genes, 
    reference.assay = "RNA", 
    query.assay = "ACTIVITY", 
    reduction = "cca"
  );

#------------------------------------
# transfer label from RNA to ATAC
rna2atac.predictions <- TransferData(
    anchorset = rna2atac.transfer.anchors, 
    refdata = rna.se$ClusterName,
    weight.reduction = atac.se[["pca"]],
    dims = 1:20
  )

rna2atac.predictions.df <- data.frame(rna2atac.predicted.id = rna2atac.predictions$predicted.id, rna2atac.predict.max.score = apply(rna2atac.predictions[,-1], 1, max))
rownames(rna2atac.predictions.df) <- colnames(atac.se) 
atac.se <- AddMetaData(atac.se, metadata = rna2atac.predictions.df)

outfname = paste(outF,"rna2atac.predictScore.hist.pdf", sep=".")
pdf(outfname)
hist(
    rna2atac.predictions.df$rna2atac.predict.max.score,
    xlab="prediction score",
    col="lightblue",
    xlim=c(0, 1),
    main="rna2atac prediction"
  );
abline(v=0.5, col="red", lwd=2, lty=2);
dev.off()

refdata <- GetAssayData(
    object = rna.se, 
    assay = "RNA", 
    slot = "data"
  );

rna2atac.imputation <- TransferData(
    anchorset = rna2atac.transfer.anchors, 
    refdata = refdata, 
    weight.reduction = atac.se[["pca"]], 
    dims = 1:20
  );

# add imputate gmat to x.sp
atac.se[["RNA"]] <- CreateAssayObject(counts = rna2atac.imputation@data);

fwrite(atac.se[[]], paste(outF, "atac.seObj.meta.tsv", sep="."), sep="\t", quote=F, col.names=T, row.names=F)
saveRDS(atac.se, file=paste(outF, "atac.seObj.rds", sep="."))

#------------------------------------------
# plot original embed with transfered label
pdf(paste(outF,"origEmbed.transferLabel.pdf", sep="."))
# transfer RNA label to ATAC
DimPlot(atac.se, group.by = "rna2atac.predicted.id", label = TRUE, repel = TRUE, reduction="umap") + ggtitle("scRNA-seq cells: match to snATAC-seq") + NoLegend()
DimPlot(atac.se, group.by = "rna2atac.predicted.id", label = FALSE, repel = FALSE, reduction="umap") + ggtitle("scRNA-seq cells: match to snATAC-seq") + NoLegend()
DimPlot(atac.se, reduction = "umap", group.by="L3cluster", label = TRUE, repel = TRUE) + ggtitle("snATAC-seq cells") + NoLegend()
DimPlot(atac.se, reduction = "umap", group.by="L3cluster", label = FALSE, repel = FALSE) + ggtitle("snATAC-seq cells") + NoLegend()
DimPlot(rna.se, group.by = "ClusterName", label = TRUE, repel = TRUE, reduction="umap") + ggtitle("scRNA-seq cells") + NoLegend() 
DimPlot(rna.se, group.by = "ClusterName", label = FALSE, repel = FALSE, reduction="umap") + ggtitle("scRNA-seq cells") + NoLegend()
dev.off()

##################################
## co-embedding
##################################
coembed <- merge(x = rna.se, y = atac.se)
rm(rna.se)
rm(atac.se)

# Finally, we run PCA and UMAP on this combined object, to visualize the co-embedding of both
# datasets
coembed <- ScaleData(coembed, features = variable.genes, do.scale = FALSE)
coembed <- RunPCA(coembed, features = variable.genes, verbose = FALSE)
coembed <- RunUMAP(coembed, dims = 1:20)
coembed <- FindNeighbors(coembed, dims = 1:20)
coembed <- FindClusters(coembed, resolution = 0.5)

coembedF <- coembed[[]]
idents=Idents(coembed)
coembedF <- data.frame(coembedF, coembed.idents=idents)
fwrite(coembedF, paste(outF, "coembed.seObj.meta.tsv", sep="."), sep="\t", quote=F, col.names=T, row.names=T)
saveRDS(coembed, file=paste(outF,"coembed.rds", sep="."))

coor_umap <- coembed@reductions$umap@cell.embeddings
meta_coembed <- coembed[[c("ClusterName","cluster", "L3Color")]]
n <- length(unique(meta_coembed$ClusterName[!is.na(meta_coembed$ClusterName)]))
cluster_color <- colorRampPalette(brewer.pal(8, "Set1"))(n)
color_map <- setNames(cluster_color, unique(meta_coembed$ClusterName[!is.na(meta_coembed$ClusterName)]))
meta_coembed$cluster_color <- color_map[match(meta_coembed$ClusterName, names(color_map))]

meta_coembed <- cbind(meta_coembed, coor_umap)
meta_coembed[is.na(meta_coembed$L3Color), "L3Color"] <- "#e6e6e6"
meta_coembed[is.na(meta_coembed$cluster_color), "cluster_color"] <- "#e6e6e6"

pdf(paste(outF, "coembed.pdf", sep="."), width=6, height = 6)
DimPlot(coembed, group.by = "tech", cols=c("red", "black"))

cols_map <- setNames(unique(meta_coembed$cluster_color), unique(meta_coembed$ClusterName))
DimPlot(coembed, group.by = "ClusterName", label = FALSE, repel = FALSE, reduction="umap", cols=cols_map) + ggtitle("coembed: scRNA") + NoLegend()
DimPlot(rev(coembed), group.by = "ClusterName", label = TRUE, repel = TRUE, reduction="umap", cols=cols_map) + ggtitle("coembed: snRNA") + NoLegend()

cols_map <- setNames(unique(meta_coembed$L3Color), unique(meta_coembed$cluster))
DimPlot(coembed, group.by = "cluster", label = FALSE, repel = FALSE, reduction="umap", cols=cols_map) + ggtitle("coembed: snATAC") + NoLegend()
DimPlot(rev(coembed), group.by = "cluster", label = TRUE, repel = TRUE, reduction="umap", cols=cols_map) + ggtitle("coembed: snATAC") + NoLegend()

plot(rev(meta_coembed$UMAP_1), rev(meta_coembed$UMAP_2), pch=19, cex=0.2, col=rev(meta_coembed$cluster_color))
plot(meta_coembed$UMAP_1, meta_coembed$UMAP_2, pch=19, cex=0.2, col=meta_coembed$L3Color)
dev.off()


