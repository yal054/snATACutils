#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input", required=TRUE, help="load subset RData")
parser$add_argument("--down_rate", nargs='+', help="cutoff of downsample rate [default %(default)s]")
parser$add_argument("--bin_size", default = 5000, help="binSize to use [default %(default)s]")
parser$add_argument("--black_list", required=TRUE, help="black list file")
parser$add_argument("--pc_num", default = 50, help="num of PCA components [default %(default)s]")
parser$add_argument("--cpu", default = 4, help="# of cpus [default %(default)s]")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

suppressPackageStartupMessages(library("SnapATAC"))
library(tictoc)

RDataF = args$input
down_rate = as.numeric(args$down_rate)
bin_size = as.numeric(args$bin_size)
black_list = args$black_list
pc_num = as.numeric(args$pc_num)
cpus = as.numeric(args$cpu)
outF = args$output

tic("loadSnap")
x.sp <- get(load(RDataF))
toc()


downSampleSnap <- function(obj, down_rate, seed.use) {
  UseMethod("downSampleSnap", obj);
}

downSampleSnap.default <- function(obj, down_rate=0.3, seed.use=10){
    nBarcode = length(obj@barcode);
    set.seed(seed.use);
    selected_barcode = sample(obj@barcode, ceiling(nBarcode * down_rate));
    selected_idx = which(obj@barcode %in% selected_barcode);
    obj@barcode = obj@barcode[selected_idx];
    obj@file = obj@file[selected_idx];
    obj@sample = obj@sample[selected_idx];
    obj@metaData = obj@metaData[selected_idx,];
    return(obj);
}

if(!is.null(down_rate)){
message("Downsample barcode by rate: ", down_rate);
x.sp <- downSampleSnap(x.sp, down_rate)
}

tic("addBmatToSnap")
x.sp = addBmatToSnap(
    x.sp,
    bin.size=bin_size,
    num.cores=cpus
    );
toc()

pdf(paste(outF, ".pre.sta.pdf",sep=""))
plotBarcode(x.sp,
              pdf.file.name=NULL,
              pdf.width=7,
              pdf.height=7,
              col="grey",
              border="grey",
              breaks=50
              );
dev.off()

summarySnap(x.sp)

## Matrix Binarization
tic("makeBinary")
x.sp = makeBinary(x.sp, mat="bmat");
toc()

tic("calBmatCor")
calBmatCor(x.sp);
toc()

## Feature Selection (filter bins according to black list and chrom)
tic("filterBins")
black_list = read.table(black_list);
suppressPackageStartupMessages(library(GenomicRanges));
black_list.gr = GRanges(
                          black_list[,1],
                          IRanges(black_list[,2], black_list[,3])
                         );
idy1 = queryHits(findOverlaps(x.sp@feature, black_list.gr));
idy2 = grep("chrM|random", x.sp@feature);
idy = unique(c(idy1, idy2));
x.sp = x.sp[,-idy, mat="bmat"];

pdf(paste(outF, ".pre.BinCoverage.pdf",sep=""))
plotBinCoverage(
    x.sp,
    pdf.file.name=NULL,
    col="grey",
    border="grey",
    breaks=10,
    xlim=c(-6,6)
    );
dev.off()

x.sp = filterBins(x.sp, low.threshold=-2, high.threshold=2, mat="bmat");
# filter empty rows
idx <- which(Matrix::rowSums(x.sp@bmat)!=0)
x.sp@metaData <- x.sp@metaData[idx,]
x.sp@barcode <- x.sp@barcode[idx]
x.sp@sample <- x.sp@sample[idx]
x.sp@file <- x.sp@file[idx]
x.sp@bmat <- x.sp@bmat[idx,]
toc()

## Jaccard Matrix & Normlaization
tic("runJaccard")
x.sp = runJaccard(
    obj = x.sp,
    tmp.folder = tempdir(),
    mat = "bmat",
    max.var=2000,
    ncell.chunk=1000,
    do.par=TRUE,
    num.cores=cpus,
    seed.use=10
    );
toc()

tic("runNormJaccard")
x.sp = runNormJaccard(
    obj = x.sp,
    tmp.folder = tempdir(),
    ncell.chunk=1000,
    method="normOVE",
    row.center=TRUE,
    row.scale=TRUE,
    low.threshold=-5,
    high.threshold=5,
    num.cores=cpus,
    seed.use=10
    );
toc()

# rmBmat
tic("rmBmat")
rmBmatFromSnap(x.sp)
toc()

## PCA / svd
tic("runDimReduct")
x.sp = runDimReduct(
    x.sp,
    pc.num=pc_num,
    input.mat="jmat",
    method="svd",
    center=TRUE,
    scale=FALSE,
    seed.use=10
    );
toc()

pdf(paste(outF, ".pre.plotDimResuce.pdf",sep=""))
plotDimReductElbow(
    obj=x.sp, 
    point.size=1.5,
    point.shape=19,
    point.color="red",
    point.alpha=1,
    pdf.file.name=NULL,
    pdf.height=7,
    pdf.width=7,
    labs.title="PCA Elbow plot",
    labs.subtitle=NULL,
    );

plotDimReductPW(
    obj=x.sp, 
    pca.dims=1:pc_num
    );
dev.off()

outfname = paste(outF, ".pre.RData",sep="")
save(x.sp, file=outfname)

outmeta = paste(outF, ".pre.meta.txt",sep="")
write.table(x.sp@metaData, file=outmeta, sep="\t", quote=F, col.names=T, row.names=F)

