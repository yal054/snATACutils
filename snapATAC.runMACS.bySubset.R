#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input", required=TRUE, help="load preprocessed RData")
parser$add_argument("--cpu", default = 12, help="# of cpus [default %(default)s]")
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
path_to_macs2 = args$path_to_macs2
path_to_snaptools = args$path_to_snaptools
outF = args$output

x.sp <- get(load(RDataF))

peaks.gr = runMACS(
    obj=x.sp,
    path.to.snaptools=path_to_snaptools,
    path.to.macs=path_to_macs2,
    output.prefix=outF,
    num.cores=cpus,
    gsize="mm",
    buffer.size=500,
    macs.options="--nomodel --shift 37 --ext 73 --qval 1e-2 -B --SPMR --call-summits",
    tmp.folder=tempdir()
    );

outfname = paste(outF, ".peaks.RData",sep="")
save(peaks.gr, file=outfname)

