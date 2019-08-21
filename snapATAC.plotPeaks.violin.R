#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input", required=TRUE, help="input z-zcore for diffpeaks")
parser$add_argument("-n", "--num", default = 5000, help="number of rows to plot [default %(default)s]")
parser$add_argument("--from", default = 0, help="z-score from [default %(default)s]")
parser$add_argument("--to", default = 10, help="z-score to [default %(default)s]")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

suppressPackageStartupMessages(library("SnapATAC"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("ggplot2"))

inF = args$input
N = as.numeric(args$num)
zfrom = as.numeric(args$from)
zto = as.numeric(args$to)
outF = args$output

diffpeaks.z <- fread(inF, sep="\t",header=T)
diffpeaks.z.sample <- diffpeaks.z[sample(nrow(diffpeaks.z), N), ]

outpdfname = paste(outF, ".diffpeaks.violin.pdf",sep="")
pdf(outpdfname, width=6, height=8)
p <- ggplot(diffpeaks.z.sample, aes(factor(clusterCell), y_z)) + geom_violin(scale = "width", aes(fill=factor(clusterCell))) + facet_grid(factor(clusterPeak) ~ .) + theme_classic() + ylim(zfrom,zto)
plot(p)
dev.off()

