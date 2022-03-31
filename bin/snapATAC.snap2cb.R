#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input", required=TRUE, help="load preprocessed RData")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()
suppressPackageStartupMessages(library("SnapATAC"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("tictoc"))

RDataF = args$input
outF = args$output

# export marker genes
# load snap
x.sp <- get(load(RDataF))

# export gmat to mtx format
gmat <- x.sp@gmat
rownames(gmat) <- paste(x.sp@sample, x.sp@barcode, sep=".")
gmat <- t(gmat)

# export
geneFname <- paste(outF, "gene.tsv", sep=".")
fwrite(as.data.frame(rownames(gmat)), geneFname, sep="\t", col.names = F, row.names = F, quote = F)
cellFname <- paste(outF, "cell.tsv", sep=".")
fwrite(as.data.frame(colnames(gmat)), cellFname, sep="\t", col.names = F, row.names = F, quote = F)
mtxFname <- paste(outF, "gmat.mtx", sep=".")
writeMM(gmat, file=mtxFname)

# export coords
barcodes.ls <- paste(x.sp@sample, x.sp@barcode, sep=".")
#tsne.coords <- data.frame(cellId=barcodes.ls, x=x.sp@tsne[,1], y=x.sp@tsne[,2])
umap.coords <- data.frame(cellId=barcodes.ls, x=x.sp@umap[,1], y=x.sp@umap[,2])

#tsneFname <- paste(outF, "tsne.coords.tsv", sep=".")
#fwrite(tsne.coords, tsneFname, sep="\t", col.names = T, row.names = F, quote = F)
umapFname <- paste(outF, "umap.coords.tsv", sep=".")
fwrite(umap.coords, umapFname, sep="\t", col.names = T, row.names = F, quote = F)

# export metatable
metaTable <- data.frame(cellId=paste(x.sp@sample, x.sp@barcode, sep="."), x.sp@metaData, cluster=x.sp@cluster)
metaFname <- paste(outF, "meta.tsv", sep=".")
fwrite(metaTable, metaFname, sep="\t", col.names = T, row.names = F, quote = F)

