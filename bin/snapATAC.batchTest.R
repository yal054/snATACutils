#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input", required=TRUE, help="load subset RData")
parser$add_argument("-b", "--batch", default = "replicate", help="batch to test")
#parser$add_argument("-m", "--meta", default = "/projects/ps-renlab/yangli/projects/CEMBA/01.joint_dat/rs1cemba/rs1cemba.full.cellMeta.tsv", help="meta table for cells")
parser$add_argument("-n", "--ncell", default = 5000, help="subsampling # of cell")
parser$add_argument("-d", "--ndim", default = 20, help="# of dims")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

suppressPackageStartupMessages(library("SnapATAC"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("RColorBrewer"))
suppressPackageStartupMessages(library("kBET"))
suppressPackageStartupMessages(library("lisi"))
library("tictoc")


RDataF = args$input
batchName = args$batch
#metaF = args$meta
ncell = as.numeric(args$ncell)
dims = as.numeric(args$ndim)
outF = args$output

x.sp <- get(load(RDataF))
#meta <- fread(metaF, sep="\t", header=T)
#x.sp@metaData <- join(x.sp@metaData[,c("sample", "barcode")], meta, by=c("sample", "barcode"))

if(nrow(x.sp)>ncell){
  set.seed(2020)
  idx <- sample(1:nrow(x.sp), ncell)
  x.sp <- x.sp[idx, ]
}


batch <- x.sp@metaData[,batchName]

###########################################
#  linear model fit of principal components
egi.data <- list()
egi.data$x <- x.sp@smat@dmat[, 1:dims]
egi.data$sdev <- x.sp@smat@sdev

batch.egi <- pcRegression(egi.data, batch, n_top = 20, tol = 1e-16)
out <- batch.egi$r2
out <- as.data.frame(out)
out$eig <- paste("eig", 1:nrow(out), sep="")
out$cluster <- tail(unlist(strsplit(outF, "/")),1)

outfname = paste(outF, ".pcRegression.res.txt",sep="")
write.table(out, outfname, quote=F, col.names=T, row.names=F, sep="\t")


##############################################
# k-Nearest neighbor batch-effect test (kBET)
data <- egi.data$x
batch.estimate <- kBET(data, batch, plot=F, k0=x.sp@graph@k) 
out <- batch.estimate$summary
out$class <- rownames(out)
out$cluster <- tail(unlist(strsplit(outF, "/")),1)

outfname = paste(outF, ".kBET.res.txt",sep="")
write.table(out, outfname, quote=F, col.names=T, row.names=F, sep="\t")


# Plot kBET's rejection rate
plot.data <- data.frame(class=rep(c('observed', 'expected'), 
                                  each=length(batch.estimate$stats$kBET.observed)), 
                        data =  c(batch.estimate$stats$kBET.observed,
                                  batch.estimate$stats$kBET.expected))

outfname = paste(outF, ".kBET.box.pdf",sep="")
pdf(outfname)
ggplot(plot.data, aes(class, data)) + geom_boxplot() +
     labs(x='Test', y='Rejection rate',title='kBET test results') +
     theme_bw() +
     scale_y_continuous(limits=c(0,1))
dev.off()

##########################################
# Local Inverse Simpson’s Index (LISI)
X <- data
meta_data <- data.frame(batch)
colnames(meta_data) <- batchName
attr.list <- batchName

res <- compute_lisi(X, meta_data, attr.list)
nlen <- length(unique(batch))
out <- data.frame(sample=x.sp@metaData$sample, barcode=x.sp@metaData$barcode, LISI=res[,1], residual=nlen-res[,1])

#meta_data0 <- meta_data
#meta_data0[,1] <- sample(meta_data0[,1])
#res0 <- compute_lisi(X, meta_data0, attr.list)

outfname = paste(outF, ".LISI.res.txt",sep="")
write.table(out, outfname, quote=F, col.names=T, row.names=F, sep="\t")

outfname = paste(outF, ".LISI.density.pdf",sep="")
pdf(outfname, width=6, height=3)
plot(density(out$LISI), xlim=c(0,nlen), main=paste("iLISI", outF, sep=" "))
dev.off()





