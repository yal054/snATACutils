#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input", required=TRUE, help="load preprocessed RData")
parser$add_argument("--pc_num", default = 50, help="num of PCA components [default %(default)s]")
parser$add_argument("-d", "--pc_dim", default = 20, help="num of PCA dims used for clustering [default %(default)s]")
parser$add_argument("--cpu", default = 1, help="# of cpus [default %(default)s]")
parser$add_argument("--black_list", default="/projects/ps-renlab/yangli/genome/mm10/mm10.blacklist.bed.gz", help="black list file")
parser$add_argument("--bin_size", default = 5000, help="binSize to use [default %(default)s]")
parser$add_argument("--resolution", default = 0.5, help="resulution for python-louvain [default %(default)s]")
parser$add_argument("--path_to_snap", default = NULL, help="path to snaps [default %(default)s]")
#parser$add_argument("--path_to_snaptools", default = "/home/yangli1/apps/anaconda3/bin/snaptools", help="path to snaptools [default %(default)s]")
#parser$add_argument("--resolution", default = 0.5, help="resulution for python-louvain [default %(default)s]")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

suppressPackageStartupMessages(library("SnapATAC"))
suppressPackageStartupMessages(library("GenomicRanges"));
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library("foreach"))
suppressPackageStartupMessages(library("doParallel"))
library("tictoc")
library("umap")

RDataF = args$input
black_list = args$black_list
pc_num = as.numeric(args$pc_num)
dims = as.numeric(args$pc_dim)
bin_size = as.numeric(args$bin_size)
cpus = as.numeric(args$cpu)
resolution = as.numeric(args$resolution)
path_to_snap = args$path_to_snap
#path_to_snaptools = args$path_to_snaptools
outF = args$output

#--------------------------------
x.sp <- get(load(RDataF))

if(is.null(path_to_snap)){
print("using original path to snap")
}else{
print("changing the path to snap")
x.sp@file <- paste(path_to_snap, x.sp@sample, ".snap", sep="")
}


cnt2sample <- table(x.sp@sample)
sampleName <- names(which(cnt2sample<=10))
idx <- which(x.sp@sample %in% sampleName)
if(length(idx)>=1){
x.sp <- x.sp[-idx, ]
}

celln <- nrow(x.sp)

if (celln >= 20000){

## 0. sampling
sample.list <- unique(x.sp@sample)
sampleSize <- 20000

x.sp@metaData$log10UQ <- log10(x.sp@metaData$UQ + 1)

tic("sampling landmark")
row.covs.dens <- density(
    x = x.sp@metaData[,"log10UQ"],
    bw = 'nrd', adjust = 1
  );

sampling_prob <- 1 / (approx(x = row.covs.dens$x, y = row.covs.dens$y, xout = x.sp@metaData[,"log10UQ"])$y + .Machine$double.eps);
set.seed(2021);
idx.landmark.ds <- sort(sample(x = seq(nrow(x.sp)), size = sampleSize, prob = sampling_prob));
x.landmark.sp = x.sp[idx.landmark.ds,];
# remove single cell from sampels
cnt2sample <- table(x.landmark.sp@sample)
sampleName <- names(which(cnt2sample==1))
idx <- which(x.landmark.sp@sample %in% sampleName)
if(length(idx)>=1){
x.landmark.sp <- x.landmark.sp[-idx, ]
}

x.query.sp = x.sp[-idx.landmark.ds,];
# remove single cell from sampels
cnt2sample <- table(x.query.sp@sample)
sampleName <- names(which(cnt2sample==1))
idx <- which(x.query.sp@sample %in% sampleName)
if(length(idx)>=1){
x.query.sp <- x.query.sp[-idx, ]
}
toc()


## 1. identify usable features
tic("addBmatToSnap")
x.landmark.sp = addBmatToSnap(x.landmark.sp, bin.size=bin_size);
toc()

tic("makeBinary")
x.landmark.sp = makeBinary(x.landmark.sp, mat="bmat");
toc()

tic("filterBins")
black_list = read.table(black_list);
black_list.gr = GRanges(
    black_list[,1], 
    IRanges(black_list[,2], black_list[,3])
  );
idy = queryHits(
    findOverlaps(x.sp@feature, black_list.gr)
  );
if(length(idy) > 0){
    x.landmark.sp = x.landmark.sp[,-idy, mat="bmat"];
  };
x.landmark.sp

chr.exclude = seqlevels(x.landmark.sp@feature)[grep("random|chrM", seqlevels(x.landmark.sp@feature))];
idy = grep(paste(chr.exclude, collapse="|"), x.landmark.sp@feature);
if(length(idy) > 0){
    x.landmark.sp = x.landmark.sp[,-idy, mat="bmat"]
  };
x.landmark.sp

bin.cov = log10(Matrix::colSums(x.landmark.sp@bmat)+1);
pdf(paste(outF, ".BinCoverage.pdf",sep=""))
hist(
    bin.cov[bin.cov > 0], 
    xlab="log10(bin cov)", 
    main="log10(Bin Cov)", 
    col="lightblue", 
    xlim=c(0, 5)
  );
dev.off()

bin.cutoff = quantile(bin.cov[bin.cov > 0], 0.95);
idy = which(bin.cov <= bin.cutoff & bin.cov > 0);
x.landmark.sp = x.landmark.sp[, idy, mat="bmat"];
x.landmark.sp

#idx = which(Matrix::rowSums(x.landmark.sp@bmat) > 500);
#x.landmark.sp = x.landmark.sp[idx,];
#x.landmark.sp
toc()

#summarySnap(x.landmark.sp)
x.landmark.sp
outfname = paste(outF, ".landmark.RData",sep="")
save(x.landmark.sp, file=outfname)

#summarySnap(x.query.sp)
x.query.sp
outfname = paste(outF, ".query.RData",sep="")
save(x.query.sp, file=outfname)


# 2. embedding
tic("runDiffusionMaps")
x.landmark.sp = runDiffusionMaps(
    obj= x.landmark.sp,
    input.mat="bmat",
    num.eigs=50
  );
x.landmark.sp@metaData$landmark = 1;
toc()


# 3.extension
tic("extension")
num.files = length(sample.list);

registerDoParallel(cpus)
#x.query.ls <- foreach (i=1:num.files, .combine=snapRbind, .inorder=TRUE) %dopar% {
x.query.ls <- foreach (i=1:num.files, .inorder=TRUE) %dopar% {
print(sample.list[i])
x.query.sub <- x.query.sp[which(x.query.sp@sample == as.character(sample.list[i])), ]
x.query.sub = addBmatToSnap(x.query.sub, do.par=FALSE, num.cores=cpus, bin.size=bin_size);
x.query.sub = makeBinary(x.query.sub);

idy = unique(queryHits(findOverlaps(x.query.sub@feature, x.landmark.sp@feature)));
x.query.sub = x.query.sub[,idy, mat="bmat"];
x.query.sub = runDiffusionMapsExtension(
    obj1=x.landmark.sp,
    obj2=x.query.sub,
    input.mat="bmat"
  );
x.query.sub@metaData$landmark = 0;
x.query.sub = rmBmatFromSnap(x.query.sub)
return(x.query.sub)
}
closeAllConnections()

x.query.sp = Reduce(snapRbind, x.query.ls);
x.landmark.sp = rmBmatFromSnap(x.landmark.sp)
x.sp = snapRbind(x.landmark.sp, x.query.sp);
x.sp = x.sp[order(x.sp@metaData[,"sample"])];
toc()

rm(x.landmark.sp)
rm(x.query.sp)

outfname = paste(outF, ".pre.RData",sep="")
save(x.sp, file=outfname)

outmeta = paste(outF, ".pre.meta.txt",sep="")
write.table(x.sp@metaData, file=outmeta, sep="\t", quote=F, col.names=T, row.names=F)


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
    labs.subtitle=NULL
    );

plotDimReductPW(
    obj=x.sp,
    eigs.dims=1:pc_num
    );
dev.off()

}else{

tic("addBmatToSnap")
x.sp = addBmatToSnap(
    x.sp,
    do.par=TRUE,
    bin.size=bin_size,
    num.cores=cpus
    );
toc()
x.sp = makeBinary(x.sp);

# identify usable features
bin.cov = log(Matrix::colSums(x.sp@bmat)+1, 10)
bin.cutoff = quantile(bin.cov[bin.cov > 0], 0.95)
idy1 = which(bin.cov >= bin.cutoff | bin.cov == 0)

suppressPackageStartupMessages(library(GenomicRanges));
black_list = read.table(black_list);
black_list.gr = GRanges(
    black_list[,1],
    IRanges(black_list[,2], black_list[,3])
);

idy2 = queryHits(findOverlaps(x.sp@feature, black_list.gr));
idy = unique(c(idy1, idy2));
features.gr = x.sp@feature[-idy];
idy = unique(queryHits(findOverlaps(x.sp@feature, features.gr)));
x.sp = x.sp[,idy, mat="bmat"];

tic("runDiffusionMaps")
x.sp = runDiffusionMaps(
    obj= x.sp,
    input.mat="bmat",
    num.eigs=pc_num
  );
toc()

x.sp <- rmBmatFromSnap(x.sp)

outfname = paste(outF, ".pre.RData",sep="")
save(x.sp, file=outfname)

outmeta = paste(outF, ".pre.meta.txt",sep="")
write.table(x.sp@metaData, file=outmeta, sep="\t", quote=F, col.names=T, row.names=F)


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
    labs.subtitle=NULL
    );

plotDimReductPW(
    obj=x.sp,
    eigs.dims=1:pc_num
    );
dev.off()

}

# 4. UMAP and Cluster

tic("runKNN")
x.sp = runKNN(
    obj=x.sp,
    eigs.dims=1:dims,
    k=50
    );
toc()

## R-igraph
tic("runCluster")
x.sp = runCluster(obj=x.sp,
    tmp.folder=tempdir(),
    louvain.lib="R-igraph",
    seed.use=10
    );
toc()

#library("leiden")
#tic("runCluster_refineResolution")
#x.sp = runCluster(
#    obj=x.sp,
#    tmp.folder=tempdir(),
#    louvain.lib="leiden",
#    resolution = resolution,
#    seed.use=10
#    );
#toc()


tic("runViz_umap")
x.sp = runViz(
    obj=x.sp,
    tmp.folder=tempdir(),
    dims=2,
    eigs.dims=1:dims,
    method="umap",
    seed.use=10
  );
toc()


pdf(paste(outF, ".cluster.pdf",sep=""))
plotViz(
    obj=x.sp,
    method="umap",
    point.size=0.5,
    point.shape=19,
    point.alpha=0.8,
    point.color=x.sp@cluster,
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
    point.color=x.sp@cluster,
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
dev.off()

#x.sp <- rmBmatFromSnap(x.sp)

outmetaf <- paste(outF, ".cluster.meta.txt", sep="")
outmetamx <- cbind(x.sp@metaData, x.sp@cluster, x.sp@umap, res=resolution)
write.table(outmetamx, outmetaf, row.names=F, col.names=T, sep="\t", quote=F)

outfname = paste(outF, ".cluster.RData",sep="")
save(x.sp, file=outfname)

graph.knn <- x.sp@graph@mat
outfname = paste(outF, ".knn.mmtx", sep="")
writeMM(graph.knn,file=outfname)

