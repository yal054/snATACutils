#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input", required=TRUE, help="load norm matrix")
parser$add_argument("-m", "--meta", required=TRUE, help="load meta table")
parser$add_argument("-a", "--attr", required=TRUE, help="used attr for UMAP plots: final_subclass")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

suppressPackageStartupMessages(library("SnapATAC"))
suppressPackageStartupMessages(library("Seurat"))
suppressPackageStartupMessages(library("GenomicRanges"));
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("plyr"))
library("Matrix")
library("tictoc")


normMx = args$input
metaMx = args$meta
attr = args$attr
rnaOUT = args$output

# load data
norm.mx <- readRDS(normMx)
meta.mx <- readRDS(metaMx)

# create seurat obj
colnames(norm.mx) <- rownames(meta.mx)
rna.se <- CreateSeuratObject(counts = norm.mx, assay = "RNA", project = "MTG")
rna.se <- AddMetaData(rna.se, metadata = meta.mx)
rna.se$tech <- "rna"
DefaultAssay(rna.se) <- "RNA"
save(rna.se, file=paste(rnaOUT, ".raw.seObj.RData", sep=""))


#-------------
# NonN
idx <- which(rna.se$class == "Non-Neuronal")
rna.NonN <- subset(rna.se, cells=colnames(x = rna.se)[idx])
save(rna.NonN, file=paste(rnaOUT, ".NonN.raw.seObj.RData", sep=""))

# normalized
rna.NonN <- NormalizeData(rna.NonN, normalization.method = "LogNormalize", scale.factor = 10000)
# Identification of highly variable features (feature selection)
rna.NonN <- FindVariableFeatures(rna.NonN, selection.method = "vst", nfeatures = 2000)
# Scaling the data
all.genes <- rownames(rna.NonN)
rna.NonN <- ScaleData(rna.NonN, features = all.genes)
# Perform linear dimensional reduction
rna.NonN <- RunPCA(rna.NonN, features = VariableFeatures(object = rna.NonN))
# Determine the ‘dimensionality’ of the dataset
pdf(paste(rnaOUT, "NonN.ElbowPlot.pdf", sep="."))
ElbowPlot(rna.NonN, ndims = 50)
dev.off()
# Cluster the cells
rna.NonN <- FindNeighbors(rna.NonN, dims = 1:20)
rna.NonN <- FindClusters(rna.NonN, resolution = 0.5)
# Run non-linear dimensional reduction (UMAP/tSNE)
rna.NonN <- RunUMAP(rna.NonN, dims = 1:20)
pdf(paste(rnaOUT, "NonN.DimPlot.umap.pdf", sep="."))
DimPlot(rna.NonN, reduction = "umap")
DimPlot(rna.NonN, reduction = "umap", group.by = attr)
DimPlot(rna.NonN, reduction = "umap", group.by = attr, label = T)
dev.off()

saveRDS(rna.NonN, file = paste(rnaOUT, ".NonN.cluster.seObj.rds", sep="") )



#-------------
# GABA
idx <- which(rna.se$class == "GABAergic")
rna.GABA <- subset(rna.se, cells=colnames(x = rna.se)[idx])
save(rna.GABA, file=paste(rnaOUT, ".GABA.raw.seObj.RData", sep=""))

# normalized
rna.GABA <- NormalizeData(rna.GABA, normalization.method = "LogNormalize", scale.factor = 10000)
# Identification of highly variable features (feature selection)
rna.GABA <- FindVariableFeatures(rna.GABA, selection.method = "vst", nfeatures = 2000)
# Scaling the data
all.genes <- rownames(rna.GABA)
rna.GABA <- ScaleData(rna.GABA, features = all.genes)
# Perform linear dimensional reduction
rna.GABA <- RunPCA(rna.GABA, features = VariableFeatures(object = rna.GABA))
# Determine the ‘dimensionality’ of the dataset
pdf(paste(rnaOUT, "GABA.ElbowPlot.pdf", sep="."))
ElbowPlot(rna.GABA, ndims = 50)
dev.off()
# Cluster the cells
rna.GABA <- FindNeighbors(rna.GABA, dims = 1:20)
rna.GABA <- FindClusters(rna.GABA, resolution = 0.5)
# Run non-linear dimensional reduction (UMAP/tSNE)
rna.GABA <- RunUMAP(rna.GABA, dims = 1:20)
pdf(paste(rnaOUT, "GABA.DimPlot.umap.pdf", sep="."))
DimPlot(rna.GABA, reduction = "umap")
DimPlot(rna.GABA, reduction = "umap", group.by = attr)
DimPlot(rna.GABA, reduction = "umap", group.by = attr, label = T)
dev.off()

saveRDS(rna.GABA, file = paste(rnaOUT, ".GABA.cluster.seObj.rds", sep="") )



#-------------
# GLUT
idx <- which(rna.se$class == "Glutamatergic")
rna.GLUT <- subset(rna.se, cells=colnames(x = rna.se)[idx])
save(rna.GLUT, file=paste(rnaOUT, ".GLUT.raw.seObj.RData", sep=""))

# normalized
rna.GLUT <- NormalizeData(rna.GLUT, normalization.method = "LogNormalize", scale.factor = 10000)
# Identification of highly variable features (feature selection)
rna.GLUT <- FindVariableFeatures(rna.GLUT, selection.method = "vst", nfeatures = 2000)
# Scaling the data
all.genes <- rownames(rna.GLUT)
rna.GLUT <- ScaleData(rna.GLUT, features = all.genes)
# Perform linear dimensional reduction
rna.GLUT <- RunPCA(rna.GLUT, features = VariableFeatures(object = rna.GLUT))
# Determine the ‘dimensionality’ of the dataset
pdf(paste(rnaOUT, "GLUT.ElbowPlot.pdf", sep="."))
ElbowPlot(rna.GLUT, ndims = 50)
dev.off()
# Cluster the cells
rna.GLUT <- FindNeighbors(rna.GLUT, dims = 1:20)
rna.GLUT <- FindClusters(rna.GLUT, resolution = 0.5)
# Run non-linear dimensional reduction (UMAP/tSNE)
rna.GLUT <- RunUMAP(rna.GLUT, dims = 1:20)
pdf(paste(rnaOUT, "GLUT.DimPlot.umap.pdf", sep="."))
DimPlot(rna.GLUT, reduction = "umap")
DimPlot(rna.GLUT, reduction = "umap", group.by = attr)
DimPlot(rna.GLUT, reduction = "umap", group.by = attr, label = T)
dev.off()

saveRDS(rna.GLUT, file = paste(rnaOUT, ".GLUT.cluster.seObj.rds", sep="") )



#--------------------
# get cluster
# clustering RNA-seq using Seurt 
# normalized 
rna.se <- NormalizeData(rna.se, normalization.method = "LogNormalize", scale.factor = 10000)

# Identification of highly variable features (feature selection)
rna.se <- FindVariableFeatures(rna.se, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(rna.se), 10)

# plot variable features with and without labels
pdf(paste(rnaOUT, "var.scatter.pdf", sep="."), width=8, height = 6)
plot1 <- VariableFeaturePlot(rna.se)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
dev.off()

# Scaling the data
all.genes <- rownames(rna.se)
rna.se <- ScaleData(rna.se, features = all.genes)

# Perform linear dimensional reduction
rna.se <- RunPCA(rna.se, features = VariableFeatures(object = rna.se))

# Determine the ‘dimensionality’ of the dataset
#rna.se <- JackStraw(rna.se, num.replicate = 100)
#rna.se <- ScoreJackStraw(rna.se, dims = 1:20)
pdf(paste(rnaOUT, "ElbowPlot.pdf", sep="."))
ElbowPlot(rna.se, ndims = 50)
dev.off()

# Cluster the cells
rna.se <- FindNeighbors(rna.se, dims = 1:20)
rna.se <- FindClusters(rna.se, resolution = 0.5)

head(Idents(rna.se), 5)

# Run non-linear dimensional reduction (UMAP/tSNE)
rna.se <- RunUMAP(rna.se, dims = 1:20)
pdf(paste(rnaOUT, "DimPlot.umap.pdf", sep="."))
DimPlot(rna.se, reduction = "umap")
DimPlot(rna.se, reduction = "umap", group.by = attr)
DimPlot(rna.se, reduction = "umap", group.by = attr, label = T)
dev.off()

saveRDS(rna.se, file = paste(rnaOUT, ".cluster.seObj.rds", sep="") )



