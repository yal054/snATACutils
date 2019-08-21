#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("--project2ref", required=TRUE, help="inferred cluster labels for every cell")
parser$add_argument("--clustLabel", required=TRUE, help="cluster labels used for subset, 2 column: label, cluster #")
parser$add_argument("--path_to_out", required=TRUE, help="path to output dir [default %(default)s]")

parser$add_argument("--pc_num", default = 50, help="num of PCA components [default %(default)s]")
parser$add_argument("--cpu", default = 4, help="# of cpus [default %(default)s]")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

suppressPackageStartupMessages(library("SnapATAC"))
library(tictoc)

RdataF = args$input
clustLabelsF = args$sub
path_to_out = args$path_to_out
outF = args$output

load(RdataF)
file.ls <- x.sp@file
new.file.ls <- sapply(strsplit(x.sp@file, "/", fixed=TRUE), tail, 1)
new.file.ls <- paste(path_to_out, new.file.ls, sep="")
x.sp@file <- new.file.ls

clustLabels <- read.table(clustLabelsF,sep="\t",header=F)
labelNames <- unique(as.character(clustLabels[,1]))
clustNames <- as.character(clustLabels[,2])

for (label in labelNames){
clustNames <- clustLabels[which(clustLabels[,1] == label),2]
print(c(label, ": ", clustNames))

x.sub.sp = x.sp[which(x.sp@cluster == clustNames),];
rmBmatFromSnap(x.sub.sp)

outfname = paste(outF, ".", label, ".sub.RData",sep="")
save(x.sub.sp, file=outfname)
}


refclustF <- "/projects/ps-renlab/yangli/projects/CEMBA/test/test.project2cluster.txt"
refclust <- read.table(refclustF, sep="\t", header=F)
colnames(refclust) <- c("ids", "refclust")

refclust$sample <- data.frame(matrix(unlist(strsplit(as.character(refclust$ids), split = "\\.")), ncol = 2, byrow = T))[,1]
refclust$barcode <- data.frame(matrix(unlist(strsplit(as.character(refclust$ids), split = "\\.")), ncol = 2, byrow = T))[,2]

cls <- 1
refclust.sub <- refclust[which(refclust$refclust == cls), ]

x.sp.ls <- list()
use.samples <- unique(refclust$sample)
for(use.sample in use.samples){
RData <- paste(use.sample, ".snap.pre.RData", sep="")
barcode.ls <- refclust.sub[which(refclust.sub$sample == use.sample), "barcode"]
load(RData)
x.sp <- rmBmatFromSnap(x.sp)
x.sp.sub <- x.sp[which(x.sp@barcode %in% barcode.ls),]
x.sp.ls[[use.sample]] <- x.sp.sub
}

file.sources = list.files("/mnt/silencer2/home/yangli/scripts/SnapATAC/R/", pattern="*.R")
sapply(paste("/mnt/silencer2/home/yangli/scripts/SnapATAC/R/", file.sources, sep=""), source,.GlobalEnv)

x.sp.merge = Reduce(snapRbind, x.sp.ls);

path_to_out <- "/projects/ps-renlab/yangli/projects/CEMBA/test/"
file.ls <- x.sp.merge@file
new.file.ls <- sapply(strsplit(x.sp.merge@file, "/", fixed=TRUE), tail, 1)
new.file.ls <- paste(path_to_out, new.file.ls, sep="")
x.sp.merge@file <- new.file.ls

bin_size <- 10000
cpus <- 2
x.sp.merge = addBmatToSnap(
    x.sp.merge,
    bin.size=bin_size,
    num.cores=cpus
    );







require(stringr)
pnames <- head(unlist(str_split(tail(unlist(str_split(snapF, "/")),1),".snap")),1)

tic("readSnap")
x.sp = createSnap(
    file=snapF,
    sample=pnames,
    num.cores=cpus
    );

x.sp = addBmatToSnap(
    x.sp,
    bin.size=bin_size,
    num.cores=cpus
    );
toc()

pdf(paste(outF, ".raw.sta.pdf",sep=""))
plotBarcode(x.sp, 
              pdf.file.name=NULL, 
              pdf.width=7, 
              pdf.height=7, 
              col="grey",
              border="grey",
              breaks=50
              );
dev.off()

## Barcode Selection
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

pdf(paste(outF, ".qc.sta.pdf",sep=""))
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
library(GenomicRanges);
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
    pca.dims=1:pc_num,
    );
dev.off()

# KNN
#tic("runKNN")
#x.sp = runKNN(
#    obj=x.sp,
#    pca.dims=1:(pc_num/2),
#    weight.by.sd=FALSE,
#    k=50,
#    nn.eps=0.0,
#    snn=TRUE,
#    snn.prune=1/15,
#    save.knn=FALSE,
#    filename=NULL
#    );
#toc()

outfname = paste(outF, ".pre.RData",sep="")
save(x.sp, file=outfname)

outmeta = paste(outF, ".pre.meta.txt",sep="")
write.table(x.sp@metaData, file=outmeta, sep="\t", quote=F, col.names=T, row.names=F)

