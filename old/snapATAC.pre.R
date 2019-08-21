#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input", required=TRUE, help="input snap format file")
parser$add_argument("--fragment_num", default = 1000, help="cutoff of fragment.num [default %(default)s]")
parser$add_argument("--mito_ratio", default = 0.1, help="cutoff of mito.ratio [default %(default)s]")
parser$add_argument("--dup_ratio", default = 0.5, help="cutoff of dup.ratio [default %(default)s]")
parser$add_argument("--umap_ratio", default = 0.8, help="cutoff of umap.ratio [default %(default)s]")
parser$add_argument("--pair_ratio", default = 0.9, help="cutoff of pair.ratio [default %(default)s]")
parser$add_argument("--bin_size", default = 5000, help="binSize to use [default %(default)s]")
parser$add_argument("--black_list", required=TRUE, help="black list file")
parser$add_argument("--pc_num", default = 50, help="num of PCA components [default %(default)s]")
parser$add_argument("--cpu", default = 5, help="# of cpus [default %(default)s]")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

library("RColorBrewer")

snapF = args$input
fragment_num = as.numeric(args$fragment_num)
mito_ratio = as.numeric(args$mito_ratio)
dup_ratio = as.numeric(args$dup_ratio)
umap_ratio = as.numeric(args$umap_ratio)
pair_ratio = as.numeric(args$pair_ratio)
bin_size = as.numeric(args$bin_size)
black_list = args$black_list
pc_num = as.numeric(args$pc_num)
cpus = as.numeric(args$cpu)
outF = args$output

suppressPackageStartupMessages(library("SnapATAC"))
x.sp = createSnap(snapF, metaData=TRUE);

pdf(paste(outF, ".raw.sta.pdf",sep=""))
plotBarcode(x.sp);
dev.off()

## Barcode Selection
# filter cells based on the following cutoffs
x.sp = filterCells(x.sp,
                   subset.names=c("fragment.num",
                                  "mito.ratio",
                                  "dup.ratio",
                                  "umap.ratio",
                                  "pair.ratio"),
                   low.thresholds=c(fragment_num, 0, 0, umap_ratio, pair_ratio),
                   high.thresholds=c(Inf, mito_ratio, dup_ratio, 1, 1)
                  );

pdf(paste(outF, ".pre.sta.pdf",sep=""))
plotBarcode(x.sp);
dev.off()

# in most cases, we recommand for 5kb resolution
x.sp = addBmat(x.sp, snapF, binSize=bin_size);
x.sp = addGmat(x.sp, snapF);

## Matrix Binarization
x.sp = makeBinary(x.sp, mat="bmat");
checkBinSize(x.sp);

## Feature Selection (filter bins according to black list and chrom)
black_list = read.table(black_list);
black_list.gr = GRanges(
                          black_list[,1],
                          IRanges(black_list[,2], black_list[,3])
                         );
idy1 = queryHits(findOverlaps(x.sp@feature, black_list.gr));
idy2 = grep("chrM|random", x.sp@feature);
idy = unique(c(idy1, idy2));
x.sp = x.sp[,-idy, mat="bmat"];


plotBinCoverage <- function(obj,...) {
  UseMethod("plotBinCoverage", obj);
}

#' @export
plotBinCoverage.default <- function(obj, ...){
    # check the input
    if(!(class(obj)=="snap")){stop(paste("Error @plotBarcode: obj is not a snap object!", sep=""))};
    cov = Matrix::colSums(obj@bmat);
    idy = seq(ncol(colSums(obj@bmat)));
    cov = cov[which(cov > 0)];
    cov = log10(cov + 1);
    cov = (cov - mean(cov)) / sd(cov);
    hist(cov, ...);
}

pdf(paste(outF, ".pre.BinCoverage.pdf",sep=""))
plotBinCoverage(x.sp);
dev.off()

x.sp = filterBins(x.sp, low.threshold=-2, high.threshold=2, mat="bmat");

## Jaccard Matrix & Normlaization
x.sp = calJaccard(
    x.sp,
    mat = "bmat",
    ncell.chunk=1000,
    max.var=5000,
    seed.use=10,
    norm.method="normOVE",
    row.center=TRUE,
    row.scale=TRUE,
    low.threshold=-5,
    high.threshold=5,
    keep.jmat=TRUE,
    do.par=TRUE,
    num.cores=cpus)

## PCA
x.sp = runPCA(
    x.sp,
    pc.num=pc_num,
    input.mat="nmat",
    method="svd",
    weight.by.sd=TRUE,
    center=TRUE,
    scale=FALSE,
    seed.use=10
    )

pdf(paste(outF, ".pre.plotPCA.pdf",sep=""))
plotPCA(x.sp, method="elbow");
plotPCA(x.sp, method="pairwise");
dev.off()

outfname = paste(outF, ".pre.RData",sep="")
save(x.sp, file=outfname)

outmeta = paste(outF, ".pre.meta.txt",sep="")
write.table(x.sp@metaData, file=outmeta, sep="\t", quote=F, col.names=T, row.names=F)

