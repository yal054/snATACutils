#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input", required=TRUE, help="load clustered RData")
parser$add_argument("-c", "--clustLabel", required=TRUE, help="2 col: anno, cluster")
parser$add_argument("-t", "--thread", default = 4, help="# of cpus [default %(default)s]")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

suppressPackageStartupMessages(library("SnapATAC"))
library("tictoc")

RDataF = args$input
clustLabelsF = args$clustLabel
cpus = as.numeric(args$thread)
outF = args$output

#--------------------------------
x.sp <- get(load(RDataF))

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

# extract snap obj by cluster and barocdes
labelNames <- names(clustLabels.ls)

subset_obj <- function(labelName, x.sp, clustLabels.ls){
cluster.ls <- clustLabels.ls[[labelName]]
idx.sub <- which(x.sp@cluster %in% cluster.ls)
x.sp.sub <- x.sp[idx.sub, ]
outfname = paste(outF, ".", labelName, ".sub.RData",sep="")
save(x.sp.sub, file=outfname)
message(paste("subset snap obj to: ", labelName, sep=""));
}

# make parallel
library(parallel)
message("subset in parallel");
mclapply(labelNames, subset_obj, x.sp = x.sp, clustLabels.ls = clustLabels.ls, mc.cores = cpus)



