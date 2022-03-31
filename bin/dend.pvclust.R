#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-p", "--pmat", required=TRUE, help="load pmat in cnt, cpm, ra")
parser$add_argument("-n", "--nboot", default=1000, help="num of bootstrap")
parser$add_argument("-t", "--ncpu", required=TRUE, help="# of cpu")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

#.libPaths("/home/yangli1/apps/anaconda3/envs/r_env/lib/R/library")
isFALSE <- function (x) {is.logical(x) && length(x) == 1L && !is.na(x) && !x}
suppressPackageStartupMessages(library("pvclust"))
suppressPackageStartupMessages(library("data.table"))
library("tictoc")

pmatF = args$pmat
ncpu = as.integer(args$ncpu)
bootN = as.integer(args$nboot)
outF = args$output

# func
cal_cv <- function(x){
  cv <- sd(x) / mean(x)
  return(cv)
}


pmx <- fread(pmatF, header=T, sep="\t")
rownames(pmx) <- pmx$peak
peak.use <- pmx$peak
pmx$peak <- NULL

# calculate CV
peak_cv <- apply(pmx, 1, cal_cv)
pdf(paste(outF, ".peak_cv.density.pdf", sep=""))
plot(density(peak_cv))
dev.off()

# extract features
m <- mean(peak_cv)
sd <- sd(peak_cv)
lb <- m-sd # left boundary
rb <- m+sd # right boundary

idx <- which(peak_cv >= 1.5)
#idx <- which(peak_cv >= rb)
pmx.sel <- pmx[idx, ]
peak.use <- peak.use[idx]

#pmx.sel <- pmx
# pvclust
set.seed(2020)
pmx.sel.log <- log(pmx.sel + 1)
pvclust.res <- pvclust(pmx.sel.log, method.dist="cor", method.hclust="ward.D2", nboot=bootN, parallel=ncpu)
outdend <- paste(outF, "dend.RData", sep=".")
save(pvclust.res, file=outdend)

pdf(paste(outF, "dend.pdf",sep="."), height = 12, width = 30)
plot(pvclust.res)
plot(pvclust.res)
pvrect(pvclust.res, alpha=0.95)
seplot(pvclust.res)
dev.off()

