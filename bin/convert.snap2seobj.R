#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input", required=TRUE, help="load snap RData")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

suppressPackageStartupMessages(library("SnapATAC"))
suppressPackageStartupMessages(library("Seurat"))
library("tictoc")

snapF = args$input
seobjF = args$output

tic("snap2seobj")
#------------------------------------
# load atac and convert to Seurat Obj

message("Loading snap object...")
sample_id = sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(snapF))
x.sp <- createSnap(snapF, sample=sample_id)
x.sp@metaData$sample <- sample_id # if you are interested in merging multiple samples, it is better to keep tracking the sample id.

# load cell-by-gene matrix
x.sp <- addGmatToSnap(x.sp, do.par=T, num.cores=1)

message("Converting snap object to seurat object...")
metaData.use = x.sp@metaData
gmat.use = t(x.sp@gmat);

colnames(x = gmat.use) = paste(x.sp@sample, x.sp@barcode, sep=".");
rownames(metaData.use) = paste(x.sp@sample, x.sp@barcode, sep=".");

se <- CreateSeuratObject(counts = gmat.use, assay = "ACTIVITY");
se <- AddMetaData(se, metadata = metaData.use);
se$tech <- "atac"
DefaultAssay(se) <- "ACTIVITY"
toc()

message("Writing seurat object to rds...")
saveRDS(se, file=paste(seobjF, "rds", sep="."))

