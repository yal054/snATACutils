#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input", required=TRUE, help="input snap format file")
parser$add_argument("--fragment_num", default = 1000, help="cutoff of fragment.num [default %(default)s]")
parser$add_argument("--umi_num", default = 500, help="cutoff of UMI.num [default %(default)s]")
parser$add_argument("--mito_ratio", default = 1, help="cutoff of mito.ratio [default %(default)s]")
parser$add_argument("--dup_ratio", default = 1, help="cutoff of dup.ratio [default %(default)s]")
parser$add_argument("--umap_ratio", default = 0, help="cutoff of umap.ratio [default %(default)s]")
parser$add_argument("--pair_ratio", default = 0, help="cutoff of pair.ratio [default %(default)s]")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

suppressPackageStartupMessages(library("SnapATAC"))
library(tictoc)

snapF = args$input
fragment_num = as.numeric(args$fragment_num)
umi_num = as.numeric(args$umi_num)
mito_ratio = as.numeric(args$mito_ratio)
dup_ratio = as.numeric(args$dup_ratio)
umap_ratio = as.numeric(args$umap_ratio)
pair_ratio = as.numeric(args$pair_ratio)
bin_size = as.numeric(args$bin_size)
black_list = args$black_list
pc_num = as.numeric(args$pc_num)
cpus = as.numeric(args$cpu)
outF = args$output

tic("readSnap")
x.sp = createSnap(snapF, metaData=TRUE, des=snapF);
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

outmeta = paste(outF, ".qc.meta.txt",sep="")
write.table(x.sp@metaData, file=outmeta, sep="\t", quote=F, col.names=T, row.names=F)

