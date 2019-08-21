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
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

suppressPackageStartupMessages(library("SnapATAC"))
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

