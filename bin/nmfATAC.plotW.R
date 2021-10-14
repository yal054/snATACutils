#!/usr/bin/env Rscript

# read results from ksklearn

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input", required=TRUE, help="input matrix")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()


library(data.table)
dataW <- fread(args$input,sep="\t")

library(pheatmap)
library(RColorBrewer)
library(viridis)
library(dendsort)

# scale by column
#tmp <- apply(dataW,2,scale)

normUnity <- function(x){
  sum <- sum(x)
  out <- x / sum(x)
}

tmp <- apply(dataW,1,normUnity)
tmp <- t(tmp)
mx <- tmp[sample(nrow(tmp), 5000), ]

sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))
mat_cluster_rows_W <- sort_hclust(hclust(dist(mx)))
#mat_cluster_cols_W <- sort_hclust(hclust(dist(t(mx))))

quantile_breaks <- function(xs, n = 30) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

#mat_breaks_W <- quantile_breaks(t(mx), n = 30)

pdf(paste(args$output,".W.pdf",sep=''))
pheatmap(
  mat               = mx,
  scale             = 'none',
  color             = viridis(30),
#  color             = viridis(length(mat_breaks_W) - 1),
#  breaks            = mat_breaks_W,
  border_color      = NA,
#  cluster_cols      = mat_cluster_cols_W,
  cluster_cols      = F,
  cluster_rows      = mat_cluster_rows_W,
  show_colnames     = FALSE,
  show_rownames     = FALSE,
  drop_levels       = TRUE,
  fontsize          = 14,
  main              = "decomp W"
)
dev.off()


norm01 <- function(x){
  min <- min(x)
  max <- max(x)
  out <- (x - min) / (max - min)
}




