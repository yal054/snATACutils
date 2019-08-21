#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--RData", required=TRUE, help="load clustered RData")
parser$add_argument("-m", "--marker", required=TRUE, help="load marker genes in txt format file")
parser$add_argument("--clustLabel", required=TRUE, help="cluster labels used for subset, 2 column: label, cluster #")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

suppressPackageStartupMessages(library("SnapATAC"))
library("ggplot2")

RDataF = as.character(args$RData)
geneF = as.character(args$marker)
clustLabelsF = args$clustLabel
outF = as.character(args$output)

load(RDataF)
clustLabels.df <- read.table(clustLabelsF)
colnames(clustLabels.df) <- c("label", "cluster")

col.anno <- clustLabels.df
rownames(col.anno) <- col.anno$V2
col.anno$cluster <- NULL

## Visulization

markermat <- read.table(geneF,sep="\t",head=F)
marker.genes <- markermat$V1

x.sp = scaleCountMatrix(
    x.sp,
    cov=x.sp@metaData$UQ,
    mat="gmat",
    method = "RPM"
    );

dat <- data.frame()
for(g in marker.genes){
  print(g)
  val <- as.numeric(x.sp@gmat[,which(x.sp@gmat@Dimnames[[2]] == g)])
  qVal <- pmin(1, val/quantile(val[which(val > 0)], 0.99));
  qVal[qVal < 0] = 0;
  qVal[qVal > 0.8] = 1;
  df <- data.frame(sample=x.sp@sample, barcode=x.sp@barcode, cluster=x.sp@cluster, Val=val, qVal=qVal, gene=g)
  dat <- rbind(dat, df)
}

aggMean.dat <- aggregate(dat[, 4:5], list(dat$cluster, dat$gene), mean)
colnames(aggMean.dat) <- c("Clusters", "Markers", "Val", "qVal")

#ggplot(aggMean.dat, aes(Clusters, Markers)) + geom_tile(aes(fill = Val), colour = "white") + scale_fill_gradient(low = "white",high = "steelblue")
#ggplot(aggMean.dat, aes(Clusters, Markers)) + geom_tile(aes(fill = Val)) + scale_fill_viridis_c() + theme_bw()

#ggplot(aggMean.dat, aes(Clusters, Markers)) + geom_tile(aes(fill = qVal), colour = "white") + scale_fill_gradient(low = "white",high = "steelblue")
#ggplot(aggMean.dat, aes(Clusters, Markers)) + geom_tile(aes(fill = qVal)) + scale_fill_viridis_c() + theme_bw()

library(pheatmap)
library(RColorBrewer)
library(viridis)
library(dendsort)
library(reshape2)

sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))

aggMean.dat.d <- dcast(aggMean.dat, Markers ~ Clusters)
rownames(aggMean.dat.d) <- aggMean.dat.d$Markers
aggMean.dat.d$Markers <- NULL

#hclust_rows <- sort_hclust(hclust(dist(aggMean.dat.d)))
hclust_cols <- sort_hclust(hclust(dist(t(aggMean.dat.d))))

outpdfname = paste(outF, ".plotGene.heatmap.pdf",sep="")
pdf(outpdfname, width=8, height=4)

pheatmap(
  mat               = aggMean.dat.d,
  scale             = 'none',
  color             = viridis(50),
  border_color      = NA,
  cluster_cols      = hclust_cols,
  cluster_rows      = F,
  annotation_col    = col.anno,
#  cluster_cols      = T,
#  cluster_rows      = T,
  show_colnames     = TRUE,
  show_rownames     = TRUE,
  drop_levels       = TRUE,
  fontsize          = 8,
  main              = "heatmap of marker genes"
)

dev.off()


