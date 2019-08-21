#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input", required=TRUE, help="load preprocessed RData")
parser$add_argument("--cpu", default = 12, help="# of cpus [default %(default)s]")
parser$add_argument("-n", "--num", default = 50, help="number of cells within one cluster [default %(default)s]")
parser$add_argument("--clustLabel", required=TRUE, help="cluster labels used for subset, 2 column: label, cluster #")
parser$add_argument("--path_to_out", required=TRUE, help="path to output dir [default %(default)s]")
parser$add_argument("--path_to_macs2", default = "/home/r3fang/anaconda2/bin/macs2", help="path to macs2 [default %(default)s]")
parser$add_argument("--path_to_snaptools", default = "/home/r3fang/anaconda2/bin/snaptools", help="path ot snaptools [default %(default)s]")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

suppressPackageStartupMessages(library("SnapATAC"))
library("tictoc")

RDataF = args$input
cpus = as.numeric(args$cpu)
Num = args$num 
clustLabelsF = args$clustLabel
path_to_out = args$path_to_out
path_to_macs2 = args$path_to_macs2
path_to_snaptools = args$path_to_snaptools
outF = args$output

load(RDataF)

file.ls <- x.sp@file
new.file.ls <- sapply(strsplit(x.sp@file, "/", fixed=TRUE), tail, 1)
new.file.ls <- paste(path_to_out, new.file.ls, sep="")
x.sp@file <- new.file.ls

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

peaks.gr = runMACSForAll(
    obj=x.sp,
    path.to.snaptools=path_to_snaptools,
    path.to.macs=path_to_macs2,
    output.prefix=outF,
    num.cores=cpus,
    min.cells=Num,
    gsize="mm", 
    buffer.size=500, 
    macs.options="--nomodel --shift 37 --ext 73 --qval 1e-2 -B --SPMR --call-summits",
    tmp.folder=tempdir()
    ); 

outfname = paste(outF, ".peaks.byGroup.RData",sep="")
save(peaks.gr, file=outfname)

runSnapAddPmat(
    obj=x.sp, 
    peak=peaks.gr, 
    path.to.snaptools=path_to_snaptools, 
    buffer.size=500,
    num.cores=cpus,
    tmp.folder=tempdir()
    );

x.sp = addPmatToSnap(
    obj=x.sp, 
    num.cores=cpus
    );

outfname = paste(outF, ".runMACS.byGroup.RData",sep="")
save(x.sp, file=outfname)


