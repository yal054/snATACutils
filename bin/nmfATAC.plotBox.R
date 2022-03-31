#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input", required=TRUE, help="input matrix")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

mx <- read.table(args$input,sep="\t",header=T)
library(ggplot2)

pdf(paste(args$output,".pdf",sep=''))
ggplot(mx, aes(ranks, contributes)) + geom_boxplot(aes(group = cut_width(ranks, 1), colour = ranks)) + theme_bw()
ggplot(mx, aes(ranks, sparseness)) + geom_boxplot(aes(group = cut_width(ranks, 1), colour = ranks)) + theme_bw()
ggplot(mx, aes(ranks, entropy)) + geom_boxplot(aes(group = cut_width(ranks, 1), colour = ranks)) + theme_bw()
dev.off()
