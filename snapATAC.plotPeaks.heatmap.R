#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input", required=TRUE, help="load preprocessed RData")
parser$add_argument("-p", "--peak", required=TRUE, help="list of diffpeaks in txt format")
parser$add_argument("--npeaks", default = 2000, help="downsample rate for peaks [default %(default)s]")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

suppressPackageStartupMessages(library("SnapATAC"))
suppressPackageStartupMessages(library("data.table"))

RDataF = args$input
diffpeakF = args$peak
npeaks = as.numeric(args$npeaks)
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
selected_agg.df <- agg.df[sample(rownames(agg.df), npeaks),]

agg.df.logfc <- log2(agg.df/rowMeans(agg.df))
selected_agg.df.logfc <- agg.df.logfc[sample(rownames(agg.df.logfc), npeaks),]

agg.df.z <- t(apply(agg.df, 1, scale))
colnames(agg.df.z) <- colnames(agg.df)
selected_agg.df.z <- agg.df.z[sample(rownames(agg.df.z), npeaks),]

#outfname = paste(outF, ".aggVal.txt",sep="")
#write.table(agg.df, outfname, sep="\t", quote=F, col.names=T, row.names=T)

library(pheatmap)
library(RColorBrewer)
library(viridis)
library(dendsort)

sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))

outpdfname = paste(outF, ".diffpeaks.heatmap.pdf",sep="")
pdf(outpdfname, width=6, height=8)
hclust_rows <- sort_hclust(hclust(dist(selected_agg.df.z)))
hclust_cols <- sort_hclust(hclust(dist(t(selected_agg.df.z))))
pheatmap(
  mat               = selected_agg.df.z,
  scale             = 'none',
  color             = viridis(50),
  border_color      = NA,
  cluster_cols      = hclust_cols,
  cluster_rows      = hclust_rows,
#  cluster_cols      = T,
#  cluster_rows      = T,
  show_colnames     = TRUE,
  show_rownames     = FALSE,
  drop_levels       = TRUE,
  fontsize          = 14,
  main              = "heatmap on diffpeaks"
)


hclust_rows <- sort_hclust(hclust(dist(selected_agg.df.logfc)))
hclust_cols <- sort_hclust(hclust(dist(t(selected_agg.df.logfc))))
pheatmap(
  mat               = selected_agg.df.logfc,
  scale             = 'none',
  color             = viridis(50),
  border_color      = NA,
  cluster_cols      = hclust_cols,
  cluster_rows      = hclust_rows,
#  cluster_cols      = T,
#  cluster_rows      = T,
  show_colnames     = TRUE,
  show_rownames     = FALSE,
  drop_levels       = TRUE,
  fontsize          = 14,
  main              = "heatmap on diffpeaks"
)


hclust_rows <- sort_hclust(hclust(dist(selected_agg.df)))
hclust_cols <- sort_hclust(hclust(dist(t(selected_agg.df))))
pheatmap(
  mat               = selected_agg.df,
  scale             = 'row',
  color             = viridis(50),
  border_color      = NA,
  cluster_cols      = hclust_cols,
  cluster_rows      = hclust_rows,
#  cluster_cols      = T,
#  cluster_rows      = T,
  show_colnames     = TRUE,
  show_rownames     = FALSE,
  drop_levels       = TRUE,
  fontsize          = 14,
  main              = "heatmap on diffpeaks"
)

dev.off()















