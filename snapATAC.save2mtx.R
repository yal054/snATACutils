#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input", required=TRUE, help="input snap format file")
parser$add_argument("-m", "--matrix", default = "gmat", help="matrix name to export")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

suppressPackageStartupMessages(library("SnapATAC"))
library(tictoc)

snapF = args$input
matName = args$matrix
outF = args$output

load(snapF)

if (matName == "gmat"){
mat.sp <- x.sp@gmat
mat.x <- paste(x.sp@sample,x.sp@barcode,sep=".")
mat.y <- colnames(x.sp@gmat)
}else if (matName == "pmat"){
mat.sp <- x.sp@pmat
mat.x <- paste(x.sp@sample,x.sp@barcode,sep=".")
mat.y <- paste(x.sp@peak@seqnames, x.sp@peak@range, sep=":")
}else if (matName == "bmat"){
mat.sp <- x.sp@bmat
mat.x <- paste(x.sp@sample,x.sp@barcode,sep=".")
mat.y <- paste(x.sp@feature@seqnames, x.sp@feature@ranges, sep=":")
} 


outMM = paste(outF, matName, "mtx",sep=".")
writeMM(mat.sp, file=outMM)

library("data.table")
outx = paste(outF, matName, "xgi", sep=".")
fwrite(as.data.frame(mat.x), file=outx, sep="\t", quote=F, col.names=F, row.names=F)

outy = paste(outF, matName, "ygi", sep=".")
fwrite(as.data.frame(mat.y), file=outy, sep="\t", quote=F, col.names=F, row.names=F)


