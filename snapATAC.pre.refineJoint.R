#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input", required=TRUE, help="input txt with the information of path of snap file and name")
parser$add_argument("--fragment_num", default = 1000, help="cutoff of fragment.num [default %(default)s]")
parser$add_argument("--umi_num", default = 500, help="cutoff of UMI.num [default %(default)s]")
parser$add_argument("--mito_ratio", default = 1, help="cutoff of mito.ratio [default %(default)s]")
parser$add_argument("--dup_ratio", default = 1, help="cutoff of dup.ratio [default %(default)s]")
parser$add_argument("--umap_ratio", default = 0, help="cutoff of umap.ratio [default %(default)s]")
parser$add_argument("--pair_ratio", default = 0, help="cutoff of pair.ratio [default %(default)s]")
parser$add_argument("--down_rate", default = 1, help="cutoff of downsample rate [default %(default)s]")
parser$add_argument("--bin_size", default = 5000, help="binSize to use [default %(default)s]")
parser$add_argument("--doublets", nargs='+', help="input txt with doublets")
parser$add_argument("--black_list", required=TRUE, help="black list file")
parser$add_argument("--cpu", default = 4, help="# of cpus [default %(default)s]")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

suppressPackageStartupMessages(library("SnapATAC"))
library(tictoc)

inputF = args$input
fragment_num = as.numeric(args$fragment_num)
umi_num = as.numeric(args$umi_num)
mito_ratio = as.numeric(args$mito_ratio)
dup_ratio = as.numeric(args$dup_ratio)
umap_ratio = as.numeric(args$umap_ratio)
pair_ratio = as.numeric(args$pair_ratio)
down_rate = as.numeric(args$down_rate)
bin_size = as.numeric(args$bin_size)
doubletsF = args$doublets
black_list = args$black_list
cpus = as.numeric(args$cpu)
outF = args$output

inputf <- read.table(inputF,sep="\t",header=F)
fnames <- as.character(inputf[,2])
pnames <- as.character(inputf[,1])

# function for downsample datasets
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

filterDoublets <- function(obj, doublets) {
  UseMethod("filterDoublets", obj);
}

filterDoublets.default <- function(obj, doublets=doubletsF){
    doublets = read.table(doubletsF, header=F, sep="\t")
    selected_idx <- which(!paste(obj@sample, obj@barcode,sep=".") %in% doublets$V1)
    obj@barcode = obj@barcode[selected_idx];
    obj@file = obj@file[selected_idx];
    obj@sample = obj@sample[selected_idx];
    obj@metaData = obj@metaData[selected_idx,];
    return(obj);
}


# create snap obj
tic("createSnap")
x.sp = createSnap(
    file=fnames,
    sample=pnames,
    num.cores=cpus
    );
toc()

# Barcode Selection
# filter cells based on the following cutoffs
tic("filterCells")
x.sp = filterCells(x.sp,
                   subset.names=c("fragment.num",
                                  "UMI",
                                  "mito.ratio",
                                  "dup.ratio",
                                  "umap.ratio",
                                  "pair.ratio"),
                   low.thresholds=c(fragment_num, umi_num, 0, 0, umap_ratio, pair_ratio),
                   high.thresholds=c(Inf, Inf, mito_ratio, dup_ratio, 1, 1)
                  );
toc()

# filter doublets
if(!is.null(doubletsF)){
  tic("filterDoublets")
  message("Filter out potential doublets: ", doubletsF);
  x.sp <- filterDoublets(x.sp, doubletsF)
  toc()
}

# downsample 
if(!is.null(down_rate)){
  if(down_rate != 1){
  tic("downsample")
  message("Downsample barcode by rate: ", down_rate);
  x.sp <- downSampleSnap(x.sp, down_rate)
  toc()
  }
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

x.sp = filterBins(x.sp, low.threshold=-Inf, high.threshold=2, mat="bmat");
toc()


outfname = paste(outF, ".pre.RData",sep="")
save(x.sp, file=outfname)

outmeta = paste(outF, ".pre.meta.txt",sep="")
write.table(x.sp@metaData, file=outmeta, sep="\t", quote=F, col.names=T, row.names=F)

