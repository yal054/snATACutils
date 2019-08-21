#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--RData", required=TRUE, help="load clustered RData")
parser$add_argument("-m", "--marker", required=TRUE, help="load marker genes in txt format file")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

suppressPackageStartupMessages(library("SnapATAC"))
library("ggplot2")

RDataF = as.character(args$RData)
geneF = as.character(args$marker)
outF = as.character(args$output)

load(RDataF)

## Visulization

markermat <- read.table(geneF,sep="\t",head=F)
marker.genes <- markermat$V1


dat <- data.frame()
for(g in marker.genes){
  print(g)
  val <- as.numeric(x.sp@gmat[,which(x.sp@gmat@Dimnames[[2]] == g)])
  qVal <- pmin(1, val/quantile(val[which(val > 0)], 0.99));
  qVal[qVal < 0] = 0;
  qVal[qVal > 0.95] = 1;
  zVal <- scale(val)
  df <- data.frame(sample=x.sp@sample, barcode=x.sp@barcode, cluster=x.sp@cluster, qVal=qVal, zVal=zVal, gene=g)
  dat <- rbind(dat, df)
}

plotdat <- dat[sample(rownames(dat), 2000),]

outpdfname = paste(outF, ".plotGene.violin.pdf",sep="")
pdf(outpdfname, width=6, height=8)
ggplot(plotdat, aes(factor(cluster), zVal)) + geom_violin(scale = "width", aes(fill=factor(gene))) + facet_grid(factor(gene) ~ .) + ylim(0,25) + theme_bw()
ggplot(plotdat, aes(factor(cluster), qVal)) + geom_violin(scale = "width", aes(fill=factor(gene))) + facet_grid(factor(gene) ~ .) + ylim(0,1) + theme_bw()
dev.off()

