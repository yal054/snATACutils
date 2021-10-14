#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input", required=TRUE, help="load preprocessed RData")
parser$add_argument("--pc_num", default = 50, help="num of PCA components [default %(default)s]")
parser$add_argument("-d", "--pc_dim", default = 20, help="num of PCA dims used for clustering [default %(default)s]")
parser$add_argument("--cpu", default = 4, help="# of cpus [default %(default)s]")
parser$add_argument("--black_list", default="/projects/ps-renlab/yangli/genome/mm10/mm10.blacklist.bed.gz", help="black list file")
parser$add_argument("--bin_size", default = 5000, help="binSize to use [default %(default)s]")
parser$add_argument("--resolution", default = 0.5, help="resulution for python-louvain [default %(default)s]")
#parser$add_argument("--path_to_snaptools", default = "/home/yangli1/apps/anaconda3/bin/snaptools", help="path to snaptools [default %(default)s]")
#parser$add_argument("--resolution", default = 0.5, help="resulution for python-louvain [default %(default)s]")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

suppressPackageStartupMessages(library("SnapATAC"))
library("tictoc")
library("umap")

RDataF = args$input
black_list = args$black_list
pc_num = as.numeric(args$pc_num)
dims = as.numeric(args$pc_dim)
bin_size = as.numeric(args$bin_size)
cpus = as.numeric(args$cpu)
resolution = as.numeric(args$resolution)
#path_to_snaptools = args$path_to_snaptools
outF = args$output

#--------------------------------
x.sp <- get(load(RDataF))

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

tic("runKNN")
x.sp = runKNN(
    obj=x.sp,
    eigs.dims=1:dims,
    k=15
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

x.sp <- rmBmatFromSnap(x.sp)

outmetaf <- paste(outF, ".cluster.meta.txt", sep="")
outmetamx <- cbind(x.sp@metaData, x.sp@cluster, x.sp@umap, res=resolution)
write.table(outmetamx, outmetaf, row.names=F, col.names=T, sep="\t", quote=F)

outfname = paste(outF, ".cluster.RData",sep="")
save(x.sp, file=outfname)

graph.knn <- x.sp@graph@mat
outfname = paste(outF, ".knn.mmtx", sep="")
writeMM(graph.knn,file=outfname)

