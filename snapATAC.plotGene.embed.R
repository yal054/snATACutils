#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--RData", required=TRUE, help="load clustered gmat RData")
parser$add_argument("-m", "--marker", required=TRUE, help="load marker genes in txt format file")
parser$add_argument("--cpu", default = 4, help="# of cpus [default %(default)s]")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

suppressPackageStartupMessages(library("SnapATAC"))
library("tictoc")
library("umap")

RDataF = as.character(args$RData)
cpus = as.numeric(args$cpu)
geneF = as.character(args$marker)
outF = as.character(args$output)

load(RDataF)

## Visulization

markermat <- read.table(geneF,sep="\t",head=F)
marker.genes <- markermat$V1

pdf(paste(outF, ".plotGene.embed.pdf",sep=""))

plotGene(
    obj=x.sp, 
    gene.names=marker.genes,
    viz.method="tsne",
    point.size=0.3,
    point.color="red",
    point.shape=19,
    background.point=TRUE,
    background.point.color="grey",
    background.point.alpha=0.3,
    background.point.size=0.1,
    background.point.shape=19,
    low.value=0.0,
    high.value=0.95,
    down.sample=5000,
    seed.use=10,
    plot.nrow=4,
    plot.ncol=4,
    pdf.file.name=NULL, 
    pdf.height=7, 
    pdf.width=7
    );

plotGene(
    obj=x.sp,
    gene.names=marker.genes,
    viz.method="umap",
    point.size=0.3,
    point.color="red",
    point.shape=19,
    background.point=TRUE,
    background.point.color="grey",
    background.point.alpha=0.3,
    background.point.size=0.1,
    background.point.shape=19,
    low.value=0.0,
    high.value=0.95,
    down.sample=5000,
    seed.use=10,
    plot.nrow=4,
    plot.ncol=4,
    pdf.file.name=NULL,
    pdf.height=7,
    pdf.width=7
    );

dev.off()


