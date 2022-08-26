#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input", required=TRUE, help="load atac RData with gmat")
parser$add_argument("-a", "--attr", required=TRUE, help="used attr for UMAP plots: majorType")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

suppressPackageStartupMessages(library("SnapATAC"))
suppressPackageStartupMessages(library("Seurat"))
suppressPackageStartupMessages(library("GenomicRanges"));
suppressPackageStartupMessages(library("data.table"))
library("tictoc")

RDataf = args$input
attr = args$attr
atacOUT = args$output

# load atac
x.sp <- get(load(RDataf))

# conver to se obj
pca.use = x.sp@smat@dmat;
eigs.dims <- 1:ncol(pca.use)
pca.use = pca.use[,eigs.dims]

umap.use = x.sp@umap

metaData.use = x.sp@metaData
gmat.use = t(x.sp@gmat);

colnames(x = gmat.use) = paste(x.sp@sample, x.sp@barcode, sep=".");
rownames(x = pca.use)  = paste(x.sp@sample, x.sp@barcode, sep=".");
rownames(x = umap.use)  = paste(x.sp@sample, x.sp@barcode, sep=".");
rownames(metaData.use) = paste(x.sp@sample, x.sp@barcode, sep=".");
colnames(x = pca.use) <- paste0("pca_", eigs.dims);
colnames(x = umap.use) <- paste0("umap_", 1:2);

atac.se <- CreateSeuratObject(counts = gmat.use);
atac.se <- AddMetaData(atac.se, metadata = metaData.use);
atac.se$tech <- "atac"

atac.se[["pca"]] <- new(Class = "DimReduc", cell.embeddings = pca.use,
           feature.loadings = matrix(0,0,0), feature.loadings.projected = matrix(0,0,0),
           assay.used ="RNA", stdev = rep(1,length(eigs.dims)),
           key ="pca_", jackstraw = new(Class = "JackStrawData"), misc = list())

atac.se[["umap"]] <- new(Class = "DimReduc", cell.embeddings = umap.use,
           feature.loadings = matrix(0,0,0), feature.loadings.projected = matrix(0,0,0),
           assay.used ="RNA", stdev = rep(1,length(1:2)),
           key ="umap_", jackstraw = new(Class = "JackStrawData"), misc = list())


pdf(paste(atacOUT, "DimPlot.umap.pdf", sep="."))
DimPlot(atac.se, reduction = "umap")
DimPlot(atac.se, reduction = "umap", group.by = attr)
DimPlot(atac.se, reduction = "umap", group.by = attr, label = T)
dev.off()


saveRDS(atac.se, file = paste(atacOUT, ".cluster.seObj.rds", sep="") )

# NonN
idx <- which(atac.se$Class == "NonN")
atac.NonN <- subset(atac.se, cells=colnames(x = atac.se)[idx])
saveRDS(atac.NonN, file=paste(atacOUT, ".NonN.cluster.seObj.rds", sep=""))

# GABA
idx <- which(atac.se$Class == "GABA")
atac.GABA <- subset(atac.se, cells=colnames(x = atac.se)[idx])
saveRDS(atac.GABA, file=paste(atacOUT, ".GABA.cluster.seObj.rds", sep=""))

# GLUT
idx <- which(atac.se$Class == "GLUT")
atac.GLUT <- subset(atac.se, cells=colnames(x = atac.se)[idx])
saveRDS(atac.GLUT, file=paste(atacOUT, ".GLUT.cluster.seObj.rds", sep=""))


