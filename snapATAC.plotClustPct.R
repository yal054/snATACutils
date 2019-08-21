#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input", required=TRUE, help="load metaTable")
parser$add_argument("-s", "--summary", required=TRUE, help="load dataset summary")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

meta = args$input
summaryF = args$summary
outF = args$output

metaTable <- read.table(meta, sep="\t", header=T)
#colnames(metaTable) <- c("sample", "barcode", "TN", "UM", "PP", "UQ", "CM", "cluster", "tsne.1", "tsne.2", "umap.1", "umap.2")
colnames(metaTable) <- c("sample", "barcode", "TN", "UM", "PP", "UQ", "CM", "cluster", "umap.1", "umap.2")
summaryTable <- read.table(summaryF, sep="\t", header=F)
colnames(summaryTable) <- c("datasets", "dissections", "regions", "parts")

merged.dat <- merge(summaryTable[,1:4], metaTable, by.x="datasets", by.y="sample")

library("ggplot2")
library("RColorBrewer")
library("reshape2")

pdfname = paste(outF, ".clustpct.pdf",sep="")
pdf(pdfname, width = 8, height = 6)

ggplot(merged.dat, aes(x=factor(cluster))) + geom_bar(aes(fill = dissections)) + theme_bw() + labs(title="dissections") + scale_fill_manual(values = colorRampPalette(brewer.pal(8,"Dark2"))(length(levels(merged.dat$dissections))))

clustTable <- table(merged.dat[,c("cluster","dissections")])
clustTable.m <- melt(clustTable)
ggplot(clustTable.m, aes(factor(cluster), factor(dissections))) + geom_point(aes(size = value, colour=value)) + theme_bw() + scale_size_continuous(range = c(0,12))

ggplot(merged.dat, aes(x=factor(cluster))) + geom_bar(aes(fill = regions)) + theme_bw() + labs(title="regions") + scale_fill_manual(values = colorRampPalette(brewer.pal(8,"Dark2"))(length(levels(merged.dat$regions))))

clustTable <- table(merged.dat[,c("cluster","regions")])
clustTable.m <- melt(clustTable)
ggplot(clustTable.m, aes(factor(cluster), factor(regions))) + geom_point(aes(size = value, colour=value)) + theme_bw() + scale_size_continuous(range = c(0,15))

ggplot(merged.dat, aes(x=factor(cluster))) + geom_bar(aes(fill = parts)) + theme_bw() + labs(title="parts") + scale_fill_manual(values = colorRampPalette(brewer.pal(8,"Dark2"))(length(levels(merged.dat$parts))))

clustTable <- table(merged.dat[,c("cluster","parts")])
clustTable.m <- melt(clustTable)
ggplot(clustTable.m, aes(factor(cluster), factor(parts))) + geom_point(aes(size = value, colour=value)) + theme_bw() + scale_size_continuous(range = c(0,15))

dev.off()


