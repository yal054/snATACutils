#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input", required=TRUE, help="peaks.bed")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

suppressPackageStartupMessages(library("GenomicRanges"))
suppressPackageStartupMessages(library("data.table"))

inF = args$input
outF = args$output

peaks.df <- fread(inF, header=F, sep="\t")
peaks.df <- as.data.frame(peaks.df)

if((x=nrow(peaks.df)) > 0L){
    peaks.gr = GenomicRanges::GRanges(peaks.df[,1], IRanges(peaks.df[,2], peaks.df[,3]));
}else{
    peaks.gr = GenomicRanges::GRanges();
}

outfname = paste(outF, ".peaksGR.RData",sep="")
save(peaks.gr, file=outfname)

