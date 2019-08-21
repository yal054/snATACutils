#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("--list", required=TRUE, help="input txt with the information of path of snap file and name")
parser$add_argument("--project2ref", required=TRUE, help="inferred cluster labels for every cell")
parser$add_argument("--clustLabel", required=TRUE, help="cluster labels used for subset, 2 column: label, cluster #")
parser$add_argument("--path_to_snap", required=TRUE, help="path to dir of snap")
parser$add_argument("--cpu", default = 4, help="# of cpus [default %(default)s]")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

suppressPackageStartupMessages(library("SnapATAC"))
library(tictoc)

listF = args$list
project2refF = args$project2ref
clustLabelsF = args$clustLabel
path_to_snap = args$path_to_snap
cpus = as.numeric(args$cpu)
outF = args$output

# create snap obj based on file list
listf <- read.table(listF,sep="\t",header=F)
fnames <- as.character(listf[,2])
pnames <- as.character(listf[,1])

tic("createSnap")
message("create snap obj");
x.sp = createSnap(
    file=fnames,
    sample=pnames,
    num.cores=cpus
    );
toc()

file.ls <- x.sp@file
new.file.ls <- sapply(strsplit(x.sp@file, "/", fixed=TRUE), tail, 1)
new.file.ls <- paste(path_to_snap, new.file.ls, sep="")
x.sp@file <- new.file.ls

IDs <- paste(x.sp@sample, x.sp@barcode, sep=".")

# load cluster label into list
message("load cluster labels");
clustLabels.ls <- list()
clustLabels <- strsplit(readLines(clustLabelsF),"\t")
for(clustLabel in clustLabels){
label <- clustLabel[1]
cluster <- as.numeric(clustLabel[2])
if(label %in% names(clustLabels.ls)){
clustLabels.ls[[label]] <- c(clustLabels.ls[[label]], cluster)
}else{
clustLabels.ls[[label]] <- c(cluster)
}
}

# load projected / inferred cluster label
message("load inferred cluster ids");
project2ref <- read.table(project2refF, sep="\t", header=F)
colnames(project2ref) <- c("ids", "project2ref")

# extract snap obj by cluster and barocdes
labelNames <- names(clustLabels.ls)

subset_snap <- function(labelName, x.sp, clustLabels.ls, project2ref, IDs){
cluster.ls <- clustLabels.ls[[labelName]]
ids.sub <- project2ref[which(project2ref$project2ref %in% cluster.ls), "ids"]
IDs.sub.idx <- which(IDs %in% ids.sub)
x.sp.sub <- x.sp[IDs.sub.idx,];
outfname = paste(outF, ".", labelName, ".sub.RData",sep="")
save(x.sp.sub, file=outfname)
message(paste("subset snap obj to: ", labelName, sep=""));
}

# make parallel
library(parallel)
message("subset in parallel");
mclapply(labelNames, subset_snap, x.sp = x.sp, clustLabels.ls = clustLabels.ls, project2ref = project2ref, IDs = IDs, mc.cores = cpus)

