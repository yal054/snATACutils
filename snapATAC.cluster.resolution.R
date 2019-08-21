#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input", required=TRUE, help="load preprocessed RData")
parser$add_argument("-d", "--pc_dim", default = 25, help="num of PCA dims used for clustering [default %(default)s]")
parser$add_argument("-r", "--resolution", default = 0.7, help="resolution of communities in Louvain clustering [default %(default)s]")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

suppressPackageStartupMessages(library("SnapATAC"))
library("tictoc")
library("umap")

RDataF = args$input
dims = args$pc_dim
resolution = args$resolution
outF = args$output

load(RDataF)

## Clustering

# python version
tic("runCluster_python-louvain")
x.sp = runCluster(
    obj=x.sp,
    louvain.lib="python-louvain",
    resolution=resolution,
    path.to.snaptools="~/apps/anaconda2/bin/snaptools",
    seed.use=10
    );
toc()

## Visulization
# tsne
tic("runViz_Rtsne")
x.sp = runViz(
    obj=x.sp, 
    dims=2,
    pca.dims=1:dims, 
    weight.by.sd=FALSE,
    method="Rtsne",
    fast_tsne_path=NULL,
    Y.init=NULL,
    seed.use=10,
    num.cores=5
    );
toc()

tic("runViz_umap")
x.sp = runViz(
    obj=x.sp, 
    dims=2,
    pca.dims=1:dims, 
    weight.by.sd=FALSE,
    method="umap",
    fast_tsne_path=NULL,
    Y.init=NULL,
    seed.use=10,
    num.cores=5
    );
toc()

#x.sp = runViz(
#    obj=x.sp, 
#    dims=2,
#    pca.dims=1:dims, 
#    weight.by.sd=FALSE,
#    method="fast_tsne",
#    fast_tsne_path="/projects/ps-renlab/yangli/scripts/FIt-SNE/bin/fast_tsne",
#    Y.init=NULL,
#    seed.use=10,
#    num.cores=5
#    );


pdf(paste(outF, ".cluster.resolution", resolution, ".pdf",sep=""))

plotViz(
    obj=x.sp,
    method="tsne",
    point.size=0.5,
    point.shape=19,
    point.alpha=0.8,
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
    method="umap",
    point.size=0.5,
    point.shape=19,
    point.alpha=0.8,
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
ensemble.ls = lapply(split(seq(length(x.sp@cluster)), x.sp@cluster), function(x){
    colMeans(x.sp@bmat[x,])
    })

# cluster using 1-cor as distance
hc = hclust(as.dist(1 - cor(t(do.call(rbind, ensemble.ls)))), method="ward.D2");
par(mfrow=c(1,1))
plot(hc, hang=-1, xlab="");
dev.off()


## Create chromatin lanscape and identify cis-elements for each cluster seperately.
outmetaf <- paste(outF, ".cluster.resolution", resolution, ".meta.txt", sep="")

outmetamx <- cbind(x.sp@metaData, x.sp@cluster, x.sp@tsne, x.sp@umap)
write.table(outmetamx, outmetaf, row.names=F, col.names=T, sep="\t", quote=F)

outfname = paste(outF, ".cluster.resolution", resolution, ".RData",sep="")
save(x.sp, file=outfname)

