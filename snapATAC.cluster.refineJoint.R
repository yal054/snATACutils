#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input", required=TRUE, help="load clustered RData")
parser$add_argument("--pc_num", default = 50, help="num of PCA components [default %(default)s]")
parser$add_argument("-d", "--dim", default = 20, help="num of PCA dims used for clustering [default %(default)s]")
parser$add_argument("--cpu", default = 1, help="# of cpus [default %(default)s]")
parser$add_argument("--seed", default = 2019, help="random seed [default %(default)s]")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

suppressPackageStartupMessages(library("SnapATAC"))
library("tictoc")
library("umap")
library("Rtsne")
library("data.table")

RDataF = args$input
dims = as.numeric(args$dim)
pc_num = as.numeric(args$pc_num)
cpus = args$cpu
seed = args$seed
outF = args$output

load(RDataF)

## Jaccard Matrix & Normlaization
tic("runJaccard")
x.sp = runJaccard(
    obj = x.sp,
    tmp.folder = tempdir(),
    mat = "bmat",
    max.var=2000,
    ncell.chunk=1000,
    do.par=TRUE,
    num.cores=cpus,
    seed.use=seed
    );
toc()

tic("runNormJaccard")
x.sp = runNormJaccard(
    obj = x.sp,
    tmp.folder = tempdir(),
    ncell.chunk=1000,
    method="normOVE",
    row.center=TRUE,
    row.scale=TRUE,
    low.threshold=-5,
    high.threshold=5,
    num.cores=cpus,
    seed.use=seed
    );
toc()

## PCA / svd
tic("runDimReduct")
x.sp = runDimReduct(
    x.sp,
    pc.num=pc_num,
    input.mat="jmat",
    method="svd",
    center=TRUE,
    scale=FALSE,
    seed.use=seed
    );
toc()


pdf(paste(outF, ".pre.plotDimResuce.pdf",sep=""))
plotDimReductElbow(
    obj=x.sp,
    point.size=1.5,
    point.shape=19,
    point.color="red",
    point.alpha=1,
    pdf.file.name=NULL,
    pdf.height=7,
    pdf.width=7,
    labs.title="PCA Elbow plot",
    labs.subtitle=NULL,
    );

plotDimReductPW(
    obj=x.sp,
    pca.dims=1:pc_num
    );
dev.off()


# KNN
tic("runKNN")
x.sp = runKNN(
    obj=x.sp,
    pca.dims=2:dims,
    weight.by.sd=FALSE,
    k=25,
    nn.eps=0.0,
#    snn=TRUE,
#    snn.prune=1/25,
    save.knn=FALSE,
    filename=NULL
    );
toc()

## Clustering
# R-igraph
tic("runCluster_R-igraph")
x.sp = runCluster(
    obj=x.sp,
    tmp.folder=tempdir(),
    louvain.lib="R-igraph",
    path.to.snaptools=NULL,
    seed.use=seed
    );
toc()

## Visulization
# tsne
tic("runViz_Rtsne")
x.sp = runViz(
    obj=x.sp,
    tmp.folder=tempdir(), 
    dims=2,
    pca.dims=2:dims, 
    weight.by.sd=FALSE,
    method="Rtsne",
    fast_tsne_path=NULL,
    Y.init=NULL,
    seed.use=seed,
    num.cores=cpus
    );
toc()

tic("runViz_umap")
x.sp = runViz(
    obj=x.sp, 
    tmp.folder=tempdir(),
    dims=2,
    pca.dims=2:dims, 
    weight.by.sd=FALSE,
    method="umap",
    fast_tsne_path=NULL,
    Y.init=NULL,
    seed.use=seed,
    num.cores=cpus
    );
toc()


pdf(paste(outF, ".cluster.pdf",sep=""))

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

plotViz(
    obj=x.sp,
    method="tsne",
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
    method="tsne", 
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

## Create chromatin lanscape and identify cis-elements for each cluster seperately.

outmetaf <- paste(outF, ".cluster.meta.txt", sep="")
outmetamx <- cbind(x.sp@sample, x.sp@metaData, x.sp@cluster, x.sp@tsne, x.sp@umap)
write.table(outmetamx, outmetaf, row.names=F, col.names=T, sep="\t", quote=F)

outfname = paste(outF, ".cluster.RData",sep="")
save(x.sp, file=outfname)

