#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input", required=TRUE, help="load clustered RData")
parser$add_argument("-p", "--peaksGR", required=TRUE, help="load peaksGR RData")
parser$add_argument("--cpu", default = 5, help="# of cpus [default %(default)s]")
parser$add_argument("--path_to_out", required=TRUE, help="path to output dir [default %(default)s]")
parser$add_argument("--path_to_snaptools", default = "/home/r3fang/anaconda2/bin/snaptools", help="path ot snaptools [default %(default)s]")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

suppressPackageStartupMessages(library("SnapATAC"))
library("tictoc")

RDataF = args$input
peaksGR = args$peaksGR
cpus = as.numeric(args$cpu)
path_to_out = args$path_to_out
path_to_macs2 = args$path_to_macs2
path_to_snaptools = args$path_to_snaptools
outF = args$output

load(RDataF)
peaks.gr <- get(load(peaksGR))

file.ls <- x.sp@file
new.file.ls <- sapply(strsplit(x.sp@file, "/", fixed=TRUE), tail, 1)
new.file.ls <- paste(path_to_out, new.file.ls, sep="")
x.sp@file <- new.file.ls

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

outfname = paste(outF, ".pmat.RData",sep="")
save(x.sp, file=outfname)

