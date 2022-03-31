#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input", required=TRUE, help="load bam")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

library("tictoc")
suppressPackageStartupMessages(library("ATACseqQC"))

bamfile <- args$input
outF <- args$output

bamfile.labels <- gsub(".bam", "", basename(bamfile))

## generate fragement size distribution
pdf(paste(outF, ".fragSizeDist.pdf",sep=""), height = 6, width = 8)
fragSize <- fragSizeDist(bamfile, bamfile.labels)
dev.off()

write.table(fragSize, paste(outF, ".fragSizeDist.txt",sep=""), sep="\t",quote = F)

