#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-r", "--ref", required=TRUE, help="load rna RDS")
parser$add_argument("-q", "--quy", required=TRUE, help="load atac RDS")
parser$add_argument("--refAttr", required=TRUE, help="used ref attr")
parser$add_argument("--quyAttr", required=TRUE, help="used quy attr")
parser$add_argument("-m", "--minN", default = 50, help="# of min cells [default %(default)s]")
parser$add_argument("-d", "--downN", default = 100, help="# of sampled cells [default %(default)s]")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

suppressPackageStartupMessages(library("SnapATAC"))
suppressPackageStartupMessages(library("Seurat"))
suppressPackageStartupMessages(library("GenomicRanges"));
suppressPackageStartupMessages(library("RColorBrewer"))
suppressPackageStartupMessages(library("viridis"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("reshape2"))
suppressPackageStartupMessages(library("pheatmap"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("plyr"))
library("tictoc")

refF <- args$ref
quyF <- args$quy
refAttr <- args$refAttr
quyAttr <- args$quyAttr
minN <- as.numeric(args$minN)
dnN <- as.numeric(args$downN)
outF = args$output

#-------------------
# load data
ref.se <- readRDS(refF)
quy.se <- readRDS(quyF)

#-------------------
# downsample
refAttr.cnt <- table(ref.se@meta.data[, refAttr])
quyAttr.cnt <- table(quy.se@meta.data[, quyAttr])


refAttr.lst <- names(which(refAttr.cnt >= minN))
quyAttr.lst <- names(which(quyAttr.cnt >= minN))


# get ref index
ref.idx <- c()

for(a in refAttr.lst){
#  print(a)
  idx <- which(ref.se@meta.data[, refAttr] == a)
  if(length(idx) >= dnN){
    set.seed(2021)
    ref.used <- sample(idx, dnN)
  }else{
    ref.used <- idx
  }
  ref.idx <- c(ref.idx, ref.used)
}

ref.se <- ref.se[, ref.idx]


# get quy index
quy.idx <- c()

for(a in quyAttr.lst){
#  print(a)
  idx <- which(quy.se@meta.data[, quyAttr] == a)
  if(length(idx) >= dnN){
    set.seed(2021)
    quy.used <- sample(idx, dnN)
  }else{
    quy.used <- idx
  }
  quy.idx <- c(quy.idx, quy.used)
}

quy.se <- quy.se[, quy.idx]

#------------------------------
# reference-based clustering
quy.se@meta.data$queryCluster <- quy.se@meta.data[, quyAttr]
ref.se@meta.data$refCluster <- ref.se@meta.data[, refAttr]

quy.se@meta.data$set <- "query"
ref.se@meta.data$set <- "ref"

se.lst <- list(ref=ref.se, query=quy.se)

# SCTransform
for (i in names(se.lst)) {
    se.lst[[i]] <- SCTransform(se.lst[[i]], verbose = TRUE, return.only.var.genes = F)
}

se.features <- SelectIntegrationFeatures(object.list = se.lst, nfeatures = 3000)
se.lst <- PrepSCTIntegration(object.list = se.lst, anchor.features = se.features, verbose = T)

# This command returns dataset 5.  We can also specify multiple refs. (i.e. c(5,6))
# reference_dataset <- which(names(se.lst) == "ref")

se.anchors <- FindIntegrationAnchors(object.list = se.lst, normalization.method = "SCT",
    anchor.features = se.features)  # , reference = reference_dataset)
se.integrated <- IntegrateData(anchorset = se.anchors, normalization.method = "SCT")

se.integrated <- RunPCA(object = se.integrated, verbose = FALSE)
ElbowPlot(se.integrated, ndims = 50)

se.integrated <- RunUMAP(object = se.integrated, dims = 1:30)
se.integrated <- FindNeighbors(se.integrated, dims = 1:30)
se.integrated <- FindClusters(se.integrated, resolution = 0.5, algorithm = 2)

#DimPlot(se.integrated, reduction = "umap", group.by = c("set", "cluster"))
#DimPlot(se.integrated, reduction = "umap", group.by = c("queryCluster", "refCluster"))

saveRDS(se.integrated, file=paste(outF, "coembed.rds", sep="."))

#-------------------------
# plot & output results
pdf(paste(outF, "coembed.umap.pdf", sep="."), width = 18, height=6)
DimPlot(se.integrated, reduction = "umap", group.by = c("set", "seurat_clusters"))
DimPlot(se.integrated, group.by = c("seurat_clusters", "tech"))
DimPlot(se.integrated, group.by = c("refCluster", "queryCluster")) + 
  theme(legend.position = "top") + guides(color = guide_legend(nrow = 4, byrow = TRUE,
    override.aes = list(size = 2.5)))
DimPlot(se.integrated, group.by = c("refCluster", "queryCluster"), label = T) + 
  theme(legend.position = "top") + guides(color = guide_legend(nrow = 4, byrow = TRUE,
    override.aes = list(size = 2.5)))
dev.off()


out <- se.integrated@meta.data
umap <- se.integrated@reductions$umap@cell.embeddings
out <- cbind(out, umap)
write.table(out, paste(outF, ".coembed.meta.txt",sep=""), sep="\t", col.names = T, row.names = T, quote = F)


pdf(paste(outF, ".coembed.query.umap.pdf",sep=""), width=6, height = 4)
ggplot(out, aes(x=UMAP_1, y=UMAP_2, color=as.factor(queryCluster))) +
  geom_point(size=0.2) +
  theme_classic()
dev.off()

pdf(paste(outF, ".coembed.ref.umap.pdf",sep=""), width=6, height = 4)
ggplot(out[which(!is.na(out$refCluster)), ], aes(x=UMAP_1, y=UMAP_2, color=as.factor(refCluster))) +
  geom_point(size=0.2) +
  theme_classic()
dev.off()

#------------------------
# overlap score
metaout <- out

#-----------------------------------------------------
# t1: table with 2 columns: coembed labels, raw labels
# t2: table with 2 columns: coembed labels, raw labels
cal_ovlpScore <- function(t1, t2){
  t1.table <- table(t1)
  t2.table <- table(t2)
  t1.pct <- apply(t1.table, 2, function(x){x/sum(x)})
  t2.pct <- apply(t2.table, 2, function(x){x/sum(x)})
  t1.labels <- colnames(t1.pct)
  t2.labels <- colnames(t2.pct)
  ovlpScore.df <- data.frame(anno1=as.character(), anno2=as.character(), ovlpScore=as.numeric())
  for(t1.label in t1.labels){
    for(t2.label in t2.labels){
      t1.pct.df <- data.frame(t1.pct[,t1.label])
      colnames(t1.pct.df) <- "t1"
      t1.pct.df$ident <- rownames(t1.pct.df)
      t2.pct.df <- data.frame(t2.pct[,t2.label])
      colnames(t2.pct.df) <- "t2"
      t2.pct.df$ident <- rownames(t2.pct.df)
      comp.df <- join(t1.pct.df, t2.pct.df, by="ident", type="full")
      comp.df[is.na(comp.df)] <- 0
      comp.df$ident <- NULL
      comp.df <- t(comp.df)
      ovlpScore <- sum(apply(comp.df, 2, min))
      out <- data.frame(anno1=t1.label, anno2=t2.label, ovlpScore=ovlpScore)
      ovlpScore.df <- rbind(ovlpScore.df, out)
    }
  }
  return(ovlpScore.df)
}


# calculate overlap
ident2ref <- data.frame(idents=metaout$seurat_clusters, rna_label=metaout$refCluster)
ident2ref <- ident2ref[complete.cases(ident2ref), ]

ident2query <- data.frame(idents=metaout$seurat_clusters, atac_label=metaout$queryCluster)
ident2query <- ident2query[complete.cases(ident2query), ]

ovlpScore.df <- cal_ovlpScore(ident2ref, ident2query)
ovlpScore.df <- as.data.frame(ovlpScore.df)
write.table(ovlpScore.df, paste(outF,"coembed.overlapScore.tsv", sep="."), sep="\t", col.names = T, row.names = F, quote = F)

# ovlpScore.df.sel <- subset(ovlpScore.df, ovlpScore.df$ovlpScore>=0.2)
ovlpScore.mx <- reshape2::dcast(ovlpScore.df, anno1~anno2, value.var="ovlpScore", fill = 0)
ovlpScore.mx <- as.data.frame(ovlpScore.mx)
write.table(ovlpScore.mx, paste(outF,"coembed.overlapScore.mx.tsv", sep="."), sep="\t", col.names = T, row.names = F, quote = F)

ovlpScore.plot <- ovlpScore.mx
rownames(ovlpScore.plot) <- ovlpScore.plot$anno1
ovlpScore.plot$anno1 <- NULL
ovlpScore.plot <- as.matrix(ovlpScore.plot)

pdf(paste(outF, ".coembed.ovlpScore.pdf",sep=""))
pheatmap(ovlpScore.plot, main = "overlap score: row: drop-seq ref; col: atac")
dev.off()




