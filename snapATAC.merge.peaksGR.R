#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input", required=TRUE, help="list of peaks.RData")
parser$add_argument("--cpu", default = 6, help="# of cpus [default %(default)s]")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

suppressPackageStartupMessages(library("SnapATAC"))
suppressPackageStartupMessages(library("GenomicRanges"))

inF = args$input
cpus = as.numeric(args$cpu)
outF = args$output

peak.files <- read.table(inF, sep="\t",header=F)

peaks.ls <- parallel::mclapply(as.list(levels(peak.files[,1])), function(x){
peaks.df <- get(load(x))
if((x=nrow(peaks.df)) > 0L){
    peaks.gr = GenomicRanges::GRanges(peaks.df[,1], IRanges(peaks.df[,2], peaks.df[,3]));
}else{
    peaks.gr = GenomicRanges::GRanges();
}}, mc.cores=cpus)

peaks.gr = GenomicRanges::reduce(do.call(c, peaks.ls));

outfname = paste(outF, ".peaksGR.RData",sep="")
save(peaks.gr, file=outfname)

