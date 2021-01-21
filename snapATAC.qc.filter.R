#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input", required=TRUE, help="input snap format file")
parser$add_argument("--fragment_num", default = 1000, help="cutoff of fragment.num [default %(default)s]")
parser$add_argument("--tsse2depth", nargs='+', help="input txt with tsse vs. depth")
parser$add_argument("--tsse_cutoff", default = 10, help="cutoff of tsse [default %(default)s]")
parser$add_argument("--cpu", default = 4, help="# of cpus [default %(default)s]")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

suppressPackageStartupMessages(library("SnapATAC"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("plyr"))
library(tictoc)

snapF = args$input
fragment_num = as.numeric(args$fragment_num)
tsse2depthF = args$tsse2depth
tsse_cutoff = as.numeric(args$tsse_cutoff)
cpus = as.numeric(args$cpu)
outF = args$output

require(stringr)
pnames <- head(unlist(str_split(tail(unlist(str_split(snapF, "/")),1),".DNA.snap")),1)

tic("readSnap")
x.sp = createSnap(
    file=snapF,
    sample=pnames,
    num.cores=cpus
    );
toc()

x.sp@metaData$sample <- pnames

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
message("Filter by TSS enrochment and depth: ", tsse2depthF);

tic("filterCells")
tsse2depth <- fread(tsse2depthF, sep="\t", header=T)
metaF <- join(x.sp@metaData, tsse2depth, by="barcode")
metaF <- metaF[, c("sample","barcode", "TN","UM","PP","UQ","CM","dup_rate","TSSe")]
metaF$log10UQ <- log10(metaF$UQ + 1)

library(viridisLite);
library(ggplot2);
p = ggplot(
    metaF,
    aes(x= log10UQ, y= TSSe)) +
    geom_point(size=0.1, col="grey") +
    theme_bw() +
    ggtitle("scatter QC") +
    ylim(0, 50) + xlim(0, 6) +
    geom_hline(yintercept=tsse_cutoff) + geom_vline(xintercept=log10(fragment_num)) +
    labs(x = "log10(UQ + 1)", y="TSS enrichment")
pdf(paste(outF, ".tsse2depth.raw.pdf",sep=""))
plot(p)
smoothScatter(metaF$log10UQ, metaF$TSSe)
abline(v=log10(fragment_num), h=tsse_cutoff)
dev.off()

x.sp@metaData <- metaF
idx.sel <- which(metaF$UQ >= fragment_num & metaF$TSSe >= tsse_cutoff)
x.sp = x.sp[idx.sel,];
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

## save
metaF <- x.sp@metaData
outfname = paste(outF, ".qc.filter.meta.txt",sep="")
fwrite(metaF, outfname, sep="\t", col.names=T, row.names=F, quote=F)

tic("saveSnap")
outfname = paste(outF, ".qc.filter.RData",sep="")
save(x.sp, file=outfname)
toc()

