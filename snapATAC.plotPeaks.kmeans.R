#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input", required=TRUE, help="load preprocessed RData")
parser$add_argument("-p", "--peak", required=TRUE, help="list of diffpeaks in txt format")
parser$add_argument("--npeaks", default = 2000, help="downsample rate for peaks [default %(default)s]")
parser$add_argument("-k", "--kmeansk", default = 15, help="k for kmeans [default %(default)s]")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

suppressPackageStartupMessages(library("SnapATAC"))
suppressPackageStartupMessages(library("data.table"))

RDataF = args$input
diffpeakF = args$peak
npeaks = as.numeric(args$npeaks)
kmeansk = as.numeric(args$kmeansk)
outF = args$output

load(RDataF)
x.sp <- makeBinary(x.sp, mat="pmat");
diffpeaks.ls <- fread(diffpeakF, sep="\t",header=T)

selected_diffpeaks <- diffpeaks.ls$name

allpeaks <- as.data.frame(x.sp@peak)$name
selected_diffpeaks_idx <- which(allpeaks %in% selected_diffpeaks)

clustName <- c()
agg.ls <- list()
clustN <- unique(x.sp@cluster)
for(clust in clustN){
    selected_barcodes_idx <- which(x.sp@cluster == clust)
    nBarcodes <- length(selected_barcodes_idx)
    selected_pmat <- x.sp@pmat[selected_barcodes_idx, selected_diffpeaks_idx]
    aggVal <- (Matrix::colSums(selected_pmat) + 1 )/ nBarcodes
#    aggVal.z = (aggVal - mean(aggVal)) / sd(aggVal);
    clustName <- paste("cluster", clust, sep="")
    agg.ls[[clustName]] <- aggVal
}

agg.df <- do.call(cbind, agg.ls)
rownames(agg.df) <- selected_diffpeaks

set.seed(10)

agg.df.z <- t(apply(agg.df, 1, scale))
colnames(agg.df.z) <- colnames(agg.df)

#outfname = paste(outF, ".aggVal.txt",sep="")
#write.table(agg.df, outfname, sep="\t", quote=F, col.names=T, row.names=T)

library(pheatmap)
library(RColorBrewer)
library(viridis)
library(dendsort)

sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))


# k-means
set.seed(10)
kclusters <- kmeans(agg.df.z, kmeansk)
anno_row <- as.data.frame(factor(kclusters$cluster))
colnames(anno_row) <- "anno_row"
tmp <- merge(anno_row, agg.df.z, by="row.names")
tmp <- tmp[order(tmp$anno_row),]

plotdat <- tmp[,3:ncol(tmp)]
rownames(plotdat) <- tmp$Row.names
row_lable <- tmp[,1:2]
rownames(row_lable) <- tmp$Row.names
row_lable$Row.names <- NULL

selected_tmp <- tmp[sample(rownames(tmp), npeaks),]
selected_tmp <- selected_tmp[order(selected_tmp$anno_row),]
selected_plotdat <- selected_tmp[,3:ncol(selected_tmp)]
rownames(selected_plotdat) <- selected_tmp$Row.names
selected_row_lable <- selected_tmp[,1:2]
rownames(selected_row_lable) <- selected_tmp$Row.names
selected_row_lable$Row.names <- NULL

hclust_cols <- sort_hclust(hclust(dist(t(selected_plotdat))))

outpdfname = paste(outF, ".diffpeaks.kmeans.pdf",sep="")
pdf(outpdfname, width=6, height=8)

pheatmap(
  mat               = selected_plotdat,
  scale             = 'none',
  color             = viridis(30),
  border_color      = NA,
  cluster_cols      = hclust_cols,
#  cluster_rows      = hclust_rows,
#  cluster_cols      = T,
  cluster_rows      = F,

  annotation_row    = selected_row_lable,
  show_colnames     = TRUE,
  show_rownames     = FALSE,
  drop_levels       = TRUE,
  fontsize          = 14,
  main              = "heatmap on diffpeaks"
)

dev.off()












