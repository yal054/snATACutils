#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--RData", required=TRUE, help="load clustered RData")
parser$add_argument("-r", "--refcluster", required=TRUE, help="load reference cluser metadata")
parser$add_argument("-n", "--num", default = 20, help="num of cell used of best match search")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

suppressPackageStartupMessages(library("SnapATAC"))

RDataF = as.character(args$RData)
refF = as.character(args$refcluster)
num = as.numeric(args$num)
outF = as.character(args$output)

load(RDataF)

refcluster <- read.table(refF,sep="\t",head=T)
refbarcode.ls <- paste(refcluster[,1], refcluster[,2], sep=".")

refcluster.dat <- data.frame(cbind(refbarcode.ls, refcluster[,8]))
colnames(refcluster.dat) <- c("barcode", "clust")

graph.snn <- x.sp@graph@mat
barcode.ls <- paste(x.sp@sample, x.sp@barcode, sep=".")

dnsamp.idx <- which(barcode.ls %in% refcluster.dat$barcode)

graph.snn.dnsamp <- graph.snn[dnsamp.idx,]

rownames(graph.snn.dnsamp) <- barcode.ls[dnsamp.idx]
colnames(graph.snn.dnsamp) <- barcode.ls

find_topn_cluster <- function(x, topn, ref2clust){
  y <- x[order(x,decreasing = T)][1:topn]
  z <- table(ref2clust[which(ref2clust$barcode %in% names(y)),2])
  clust <- names(which(z==max(z)))[1]
  return(clust)
}

top_match <- apply(graph.snn.dnsamp,2,find_topn_cluster, topn=num, ref2clust=refcluster.dat)
top_match.df <- data.frame(barcodes=names(top_match), project2cluster=as.numeric(unlist(top_match)))

outfname = paste(outF, ".project2ref.txt",sep="")
write.table(top_match.df, outfname, quote=F, col.names=T, row.names=F, sep="\t")
