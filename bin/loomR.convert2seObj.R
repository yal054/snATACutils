#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input", required=TRUE, help="load loom matrix")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

suppressPackageStartupMessages(library("Seurat"))
library("tictoc")

loomF = args$input
outF = args$output

loom.ls <- readRDS(loomF)
raw.matrix <- loom.ls[['raw.matrix']]
gene.attrs <- loom.ls[['gene.attrs']]
cell.attrs <- loom.ls[['cell.attrs']]
knn <- loom.ls[['knn']]
mknn <- loom.ls[['mknn']]

# filter out cortex, olfactory, cerebrel nuclei, hippocampus
#regions.ls <- unique(cell.attrs[, "Region"])
#regions.sel <- c("CNS", "Cortex", "Hippocampus", "Hippocampus,Cortex", "Olfactory bulb", "Striatum dorsal", "Striatum ventral", "Dentate gyrus", "Striatum dorsal,Striatum ventra", "Striatum dorsal, Striatum ventral, Dentate gyrus", "Pallidum", "Striatum dorsal, Striatum ventral,Amygdala", "Striatum dorsal, Striatum ventral", "Telencephalon", "Brain")

#regions.sel.idx <- which(cell.attrs$Region %in% regions.sel)
#matrix.sel <- raw.matrix[, regions.sel.idx]
#cell.attrs.sel <- cell.attrs[regions.sel.idx, ]
#knn.sel <- knn[regions.sel.idx, regions.sel.idx]
#mknn.sel <- mknn[regions.sel.idx, regions.sel.idx]

# create seurat obj
rna.se <- CreateSeuratObject(counts = raw.matrix, project = "mouseBrain", assay = "RNA")
colnames(cell.attrs) <- make.names(colnames(cell.attrs))
rna.se <- AddMetaData(rna.se, metadata = cell.attrs);

rna.se[["orig.ident"]] <- rna.se[["ClusterName"]]
rna.se@graphs$RNA_nn <- knn
rna.se@graphs$RNA_snn <- mknn

umapF <- as.matrix(cell.attrs[,c("X_X", "X_Y")])
colnames(umapF) <- c("UMAP_1", "UMAP_2")
rna.se[["umap"]] <- new(Class = "DimReduc", cell.embeddings = umapF,
           feature.loadings = matrix(0,0,0), feature.loadings.projected = matrix(0,0,0),
           assay.used ="RNA", stdev = rep(1,ncol(umapF)),
           key ="UMAP_", jackstraw = new(Class = "JackStrawData"), misc = list())

outFname <- paste(outF, "seObj.rds", sep=".")
saveRDS(rna.se, file = outFname)




