#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input", required=TRUE, help="load preprocessed RData")
parser$add_argument("-p", "--partition", required=TRUE, help="partition for python-leiden")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

suppressPackageStartupMessages(library("SnapATAC"))
library("tictoc")
library("umap")

RDataF = args$input
partitionF = args$partition
outF = args$output

load(RDataF)
partition <- read.table(partitionF, sep="\t", head=F)
x.sp@cluster <- as.factor(partition[,1])

pdf(paste(outF, ".refineCluster.pdf", sep=""))
plotViz(
    obj= x.sp,
    method="umap",
    main="Cluster",
    point.color=x.sp@cluster,
    point.size=0.2,
    point.shape=19,
    text.add=FALSE,
    text.size=1,
    text.color="black",
    down.sample=10000,
    legend.add=FALSE
  );

plotViz(
    obj= x.sp,
    method="umap",
    main="Cluster",
    point.color=x.sp@cluster,
    point.size=0.2,
    point.shape=19,
    text.add=TRUE,
    text.size=1,
    text.color="black",
    down.sample=10000,
    legend.add=FALSE
  );
dev.off()
## Heretical clustering of the clusters
# calculate the ensemble signals for each cluster
#ensemble.ls = lapply(split(seq(length(x.sp@cluster)), x.sp@cluster), function(x){
#    colMeans(x.sp@bmat[x,])
#    })

# cluster using 1-cor as distance
#hc = hclust(as.dist(1 - cor(t(do.call(rbind, ensemble.ls)))), method="ward.D2");
#par(mfrow=c(1,1))
#plot(hc, hang=-1, xlab="");
#dev.off()

## Create chromatin lanscape and identify cis-elements for each cluster seperately.

outmetaf <- paste(outF, ".refineCluster.meta.txt", sep="")
outmetamx <- data.frame(cellID=paste(x.sp@sample,x.sp@barcode,sep="."), x.sp@metaData, cluster=x.sp@cluster, x.sp@umap)
write.table(outmetamx, outmetaf, row.names=F, col.names=T, sep="\t", quote=F)

#outfname = paste(outF, ".refineCluster.RData", sep="")
#save(x.sp, file=outfname)

