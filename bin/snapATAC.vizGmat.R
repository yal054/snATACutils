#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input", required=TRUE, help="load clustered RData")
parser$add_argument("-n", "--num", default=5000, help="# of cells used for smooth")
parser$add_argument("-g", "--gene", required=TRUE, help="gene marker table")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

# use new version
suppressPackageStartupMessages(library("SnapATAC"));
suppressPackageStartupMessages(library("GenomicRanges"));
suppressPackageStartupMessages(library("plyr"));
library("tictoc")

RDataF = args$input
num = args$num
markerF <- args$gene
outF = args$output

#-------------------------------
# load data
tic("loadSnap")
x.sp <- get(load(RDataF))
toc()

#--------------------------
# load gmat
tic("addGmat")
x.sp <- addGmatToSnap(x.sp)
toc()

save(x.sp, file=paste(outF, ".Gmat.RData", sep=""))

#----------------------------------
# norm and smooth
set.seed(2021)
if(nrow(x.sp) > num){
x.sp <- x.sp[sample(1:nrow(x.sp), num)]
}

tic("runKNN")
x.sp = runKNN(
    obj=x.sp,
    eigs.dims=1:20,
    k=50
)
toc()


x.sp = scaleCountMatrix(
    x.sp,
    cov=x.sp@metaData$UQ,
    mat="gmat",
    method = "RPM"
    );

save(x.sp, file=paste(outF, ".normGmat.RData", sep=""))

# smooth the cell-by-gene matrix
x.sp = runMagic(
    obj=x.sp,
    input.mat="gmat",
    step.size=3
  );

save(x.sp, file=paste(outF, ".smoothGmat.RData", sep=""))


#-----------------
# rna marker
#markerF <- "rna.markers.txt"
markermat <- read.table(markerF, sep="\t",head=F)
marker.genes <- unique(markermat$V1)
marker.genes <- marker.genes[which(marker.genes %in% colnames(x.sp@gmat))]

pdf(paste(outF, ".smooth.plotGene.pdf", sep=""))
par(mfrow = c(2,2));
for(i in 1:length(marker.genes)){
plotFeatureSingle(
        obj=x.sp,
        feature.value=x.sp@gmat[, as.character(marker.genes[i])],
        method="umap", 
        main=marker.genes[i],
        point.size=0.1, 
        point.shape=19, 
        down.sample=5000,
        quantiles=c(0.01, 0.99)
  )
}
par(mfrow=c(1,1));
dev.off()



