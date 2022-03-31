#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-g", "--gene", required=TRUE, help="input gene symbol list")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

geneF = args$gene
outPrefix = args$output

suppressPackageStartupMessages(library(enrichR))
suppressPackageStartupMessages(library(data.table))

# load bed
geneF <- fread(geneF, sep="\t", header=F)
gene.lst <- geneF[[1]]

dbs <- listEnrichrDbs()
dbs <- c("GO_Molecular_Function_2018", "GO_Cellular_Component_2018", "GO_Biological_Process_2018" ,"KEGG_2019_Mouse")
enriched <- enrichr(gene.lst, dbs)

# output
bp <- enriched[["GO_Biological_Process_2018"]]
mf <- enriched[["GO_Molecular_Function_2018"]]
cc <- enriched[["GO_Cellular_Component_2018"]]
kegg <- enriched[["KEGG_2019_Mouse"]]

fwrite(bp, paste(outPrefix, ".enrichR.bp.tsv", sep=""), sep="\t", col.names = T, row.names = F, quote = F)
fwrite(mf, paste(outPrefix, ".enrichR.mf.tsv", sep=""), sep="\t", col.names = T, row.names = F, quote = F)
fwrite(cc, paste(outPrefix, ".enrichR.cc.tsv", sep=""), sep="\t", col.names = T, row.names = F, quote = F)
fwrite(kegg, paste(outPrefix, ".enrichR.res.tsv", sep=""), sep="\t", col.names = T, row.names = F, quote = F)

