#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--bed", required=TRUE, help="input bed")
parser$add_argument("-s", "--species", default="mm10", help="genome hg19, mm10")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

bedF = args$bed
sp = args$species
outPrefix = args$output

suppressPackageStartupMessages(library(rGREAT))
suppressPackageStartupMessages(library(data.table))
library("tictoc")

# load bed
tic("loadBED")
bed_list = fread(bedF);
bed.gr = GRanges(
    bed_list[[1]],
    IRanges(bed_list[[2]], bed_list[[3]])
  );
toc()

# submit jobs
tic("submit jobs")
job = submitGreatJob(gr=bed.gr, species=sp)
toc()

# access result tables
tic("getTab")
tb = getEnrichmentTables(job)
toc()


# output
tic("output")

bp <- tb[["GO Biological Process"]]
fwrite(bp, paste(outPrefix, ".rGREAT.bp.tsv", sep=""), sep="\t", col.names = T, row.names = F, quote = F)

mf <- tb[["GO Molecular Function"]]
fwrite(mf, paste(outPrefix, ".rGREAT.mf.tsv", sep=""), sep="\t", col.names = T, row.names = F, quote = F)

cc <- tb[["GO Cellular Component"]]
fwrite(cc, paste(outPrefix, ".rGREAT.cc.tsv", sep=""), sep="\t", col.names = T, row.names = F, quote = F)

pdf(paste(outPrefix, ".rGREAT.dist.pdf", sep=""), width=12, height=4)
plotRegionGeneAssociationGraphs(job)
dev.off()

res = plotRegionGeneAssociationGraphs(job)
res <- as.data.frame(res)
fwrite(res, paste(outPrefix, ".rGREAT.res.tsv", sep=""), sep="\t", col.names = T, row.names = F, quote = F)
toc()

