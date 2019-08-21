#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input", required=TRUE, help="load pmat RData")
parser$add_argument("-n", "--ncell", default = 5000, help="split column to ncell * peak matrix[default %(default)s]")
parser$add_argument("--gtf", required=TRUE, help="anno gtf [default %(default)s]")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

suppressPackageStartupMessages(library("rtracklayer"))
library("scales", lib.loc="/home/yangli1/R/x86_64-pc-linux-gnu-library/3.4")
suppressPackageStartupMessages(library("SnapATAC"))
.libPaths("/home/yangli1/R/x86_64-pc-linux-gnu-library/3.4")
suppressPackageStartupMessages(library("Seurat"))
suppressPackageStartupMessages(library("tictoc"))

RDataF = args$input
n = as.numeric(args$ncell)
gtfF = args$gtf
outF = args$output

load(RDataF)

pmat <- t(x.sp@pmat)
colnames(pmat) <- paste(x.sp@sample, x.sp@barcode, sep=".")
peakinfo=data.frame(x.sp@peak)[,1:3]
rownames(pmat) = paste(peakinfo[,1], ":", peakinfo[,2], "-", peakinfo[,3], sep = "")

rm(x.sp)

# create a gene activity matrix from the peak matrix and GTF, using chromosomes 1:22, X, and Y.
# Peaks that fall within gene bodies, or 2kb upstream of a gene, are considered

# split matrix
splitN <- ceiling(ncol(pmat) / n)
pmat.ls <- list()

for (i in 1:splitN){
  j = i-1
  if(j==0){
    pmat.ls[[i]] <- pmat[,1:(n*i)]
    cat("round ", i, " : from ", 1, " to ", i*n, "\n")
  }else{
    if(n*i <= ncol(pmat)){
      pmat.ls[[i]] <- pmat[,((j*n)+1):(n*i)]
      cat("round ", i, " : from ", (j*n+1), " to ", n*i, "\n")
    }else{
      pmat.ls[[i]] <- pmat[,((j*n)+1):ncol(pmat)]
      cat("round ", i, " : from ", j*n+1, " to ", ncol(pmat), "\n")
    }
  }
}
test.dat <- pmat[,1:2]
rm(pmat)

# pmat.t <- t(pmat)
# pmat.ls <- split(pmat.t, rep(1:splitN, each = n, length.out=ncol(pmat)))
test <- CreateGeneActivityMatrix(peak.matrix = test.dat,
                                            annotation.file = gtfF,
                                            seq.levels =  c(paste("chr", 1:19,sep = ""), "chrX", "chrY"),
                                            upstream = 2000, verbose = TRUE)

tic("tic")
activity.matrix.out <- Matrix(nrow=dim(test)[1], ncol=0)
rm(test)
for (i in 1:splitN){
activity.matrix <- CreateGeneActivityMatrix(peak.matrix = pmat.ls[[i]], 
                                            annotation.file = gtfF, 
                                            seq.levels =  c(paste("chr", 1:19,sep = ""), "chrX", "chrY"), 
                                            upstream = 2000, verbose = TRUE)
activity.matrix.out <- Matrix::cBind(activity.matrix.out, activity.matrix)
}
toc()

outfname = paste(outF, ".geneActivity.RData",sep="")
save(activity.matrix.out, file=outfname)

