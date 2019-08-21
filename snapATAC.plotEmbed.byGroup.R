#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input", required=TRUE, help="load clustered RData")
parser$add_argument("--clustLabel", required=TRUE, help="cluster labels used for subset, 2 column: label, cluster #")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

suppressPackageStartupMessages(library("SnapATAC"))
library("tictoc")
library("umap")
library("Rtsne")

RDataF = args$input
clustLabelsF = args$clustLabel
outF = args$output

load(RDataF)
clustLabels.df <- read.table(clustLabelsF)

# map cluster to group id
library(plyr);
current.cluster.ids <- clustLabels.df[,2]
new.cluster.ids <- as.character(clustLabels.df[,1])

x.sp@cluster <- plyr::mapvalues(
    x = x.sp@cluster,
    from = current.cluster.ids,
    to = new.cluster.ids
    );

pdf(paste(outF, ".cluster.byGroup.pdf",sep=""))
#plotViz(
#    obj=x.sp,
#    method="tsne",
#    point.size=0.5,
#    point.shape=19,
#    point.alpha=0.8,
#    point.color="cluster",
#    text.add=FALSE,
#    down.sample=10000,
#    pdf.file.name=NULL,
#    pdf.width=7,
#    pdf.height=7
#    );

#plotViz(
#    obj=x.sp, 
#    method="tsne", 
#    point.size=0.5, 
#    point.shape=19, 
#    point.alpha=0.8, 
#    point.color="cluster",
#    text.add=TRUE,
#    text.size=1,
#    text.color="black",
#    text.halo.add=TRUE,
#    text.halo.color="white",
#    text.halo.width=0.2,
#    down.sample=10000,
#    pdf.file.name=NULL,
#    pdf.width=7, 
#    pdf.height=7
#    );

plotViz(
    obj=x.sp,
    method="umap",
    point.size=0.5,
    point.shape=19,
    point.alpha=0.8,
    point.color="cluster",
    text.add=FALSE,
    down.sample=10000,
    pdf.file.name=NULL,
    pdf.width=7,
    pdf.height=7
    );

plotViz(
    obj=x.sp, 
    method="umap", 
    point.size=0.5, 
    point.shape=19, 
    point.alpha=0.8, 
    point.color="cluster",
    text.add=TRUE,
    text.size=1,
    text.color="black",
    text.halo.add=TRUE,
    text.halo.color="white",
    text.halo.width=0.2,
    down.sample=10000,
    pdf.file.name=NULL,
    pdf.width=7, 
    pdf.height=7
    );


## Heretical clustering of the clusters
# calculate the ensemble signals for each cluster
#ensemble.ls = lapply(split(seq(length(x.sp@cluster)), x.sp@cluster), function(x){
#    colMeans(x.sp@bmat[x,])
#    })

# cluster using 1-cor as distance
#hc = hclust(as.dist(1 - cor(t(do.call(rbind, ensemble.ls)))), method="ward.D2");
#par(mfrow=c(1,1))
#plot(hc, hang=-1, xlab="");
dev.off()

