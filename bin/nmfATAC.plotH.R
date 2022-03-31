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
dataH <- fread(args$input,sep="\t")

library(pheatmap)
library(RColorBrewer)
library(viridis)
library(dendsort)

# scale by column
#mx <- apply(dataH,2,scale)

normUnity <- function(x){
  sum <- sum(x)
  out <- x / sum(x)
}

mx <- apply(dataH,2,normUnity)

sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))
#mat_cluster_rows_H <- sort_hclust(hclust(dist(dataH)))
mat_cluster_cols_H <- sort_hclust(hclust(dist(t(mx))))

quantile_breaks <- function(xs, n = 30) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

# mat_breaks_H <- quantile_breaks(t(mx), n = 30)

pdf(paste(args$output,".H.pdf",sep=''))
pheatmap(
  mat               = mx,
  scale             = 'none',
  color             = viridis(30),
#  color             = viridis(length(mat_breaks_H) - 1),
#  breaks            = mat_breaks_H,
  border_color      = NA,
  cluster_cols      = mat_cluster_cols_H,
  cluster_rows      = F,
#  cluster_rows      = mat_cluster_rows_H,
  show_colnames     = TRUE,
  show_rownames     = FALSE,
  drop_levels       = TRUE,
  fontsize          = 14,
  main              = "decomp H"
)
dev.off()


nor01 <- function(x){
  min <- min(x)
  max <- max(x)
  out <- (x - min) / (max - min)
}

