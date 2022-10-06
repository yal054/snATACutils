#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input", required=TRUE, help="load clustered RData")
parser$add_argument("-p", "--pset", required=TRUE, help="column name for used peak set: MajorType | SubType")
parser$add_argument("--path_to_peak", default="/projects/ps-renlab/yangli/projects/CEMBA/01.joint_dat/rs1cemba/tracks/cCREs/", help="path to peak set for each cluster")
parser$add_argument("-a", "--attr", required=TRUE, help="attribute's column name used for regression: MajorRegion | SubRegion | DissectionRegion | SubType")
parser$add_argument("-l", "--lrTest", required=TRUE, help="LR test results")
parser$add_argument("-s", "--sig", default = 0.01, help="significant cutoff for LR test [default %(default)s]")
parser$add_argument("-e", "--eFDR", default = 0.01, help="empirical FDR cutoff [default %(default)s]")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

suppressPackageStartupMessages(library("SnapATAC"))
suppressPackageStartupMessages(library("GenomicRanges"));
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("philentropy"))
suppressPackageStartupMessages(library("fdrtool"))
#suppressPackageStartupMessages(library("fitdistrplus"))
library("tictoc")

RDataF = args$input
pset = args$pset
path_to_peak = args$path_to_peak
lrF = args$lrTest
attr = args$attr
sig = as.numeric(args$sig)
eFDR = as.numeric(args$eFDR)
outF = args$output

#---------
# init
x.sp <- get(load(RDataF))

lr.df <- fread(lrF, sep="\t", header=T)
colnames(lr.df) <- c("id", "ncell_acc", "peak", "p_value", "q_value") 
#lr.df <- subset(lr.df, lr.df$ncell_acc>0)

#----------------
# filter used peaks
type.lst <- unique(x.sp@metaData[,pset])

if(length(type.lst)==1){
inf = list.files(path = path_to_peak, pattern = paste("*", type.lst, "*.bed",sep=""), full.names=T )
peak.lst <- fread(inf, sep="\t", header=F);
}else{
peak.lst <- list()
for(i in type.lst){
  inf = list.files(path = path_to_peak, pattern = paste("*", i, "*.bed",sep=""), full.names=T )
  inp = fread(inf, sep="\t", header=F);
  peak.lst[[i]] <- inp
}
peak.lst <- do.call(rbind, peak.lst)
peak.lst <- unique(peak.lst)
}

peak.gr = GRanges(
    peak.lst[[1]],
    IRanges(peak.lst[[2]], peak.lst[[3]])
  );
peak.gr$peak <- peak.lst[[4]]

#-----------------------
# filter snap obj

idy = queryHits(
    findOverlaps(x.sp@peak, peak.gr)
  );
x.sp <- x.sp[, idy, "pmat"]


#-----------------------
# Jensen-Shannon divergence
# 1. proprotion for each cluster
attr.lst <- unique(x.sp@metaData[, attr])

tic("get/norm prop")
prop.lst <- list()
for(t in attr.lst){
#  print(t)
  idx <- which(x.sp@metaData[, attr] == t)
  n <- length(idx)
  prop <- Matrix::colSums(x.sp@pmat[idx, ]>0)
  prop <- prop / n
  prop.lst[[t]] <- prop
}

prop.df <- do.call(cbind, prop.lst)

# 2. scale factor
medacc.lst <- list()
for(t in attr.lst){
#  print(t)
  idx <- which(x.sp@metaData[, attr] == t)
  medacc <- Matrix::rowSums(x.sp@pmat[idx, ]>0)
  medacc <- median(medacc)
  medacc.lst[[t]] <- medacc
}

medacc.df <- do.call(cbind, medacc.lst)
scalef.df <- mean(log10(medacc.df[1,])) / log10(medacc.df)

# 3. correct/norm and scale
prop.norm.df <- sweep(prop.df, 2, scalef.df[1,], "*")
prop.scale.df <- t(apply(prop.norm.df, 1, function(x){x/sum(x)}))
prop.scale.df[is.na(prop.scale.df)] <- 0
toc()

# 4. Jensen–Shannon divergence & Speciﬁcity
calJSD <- function(p, q){
    P <- p
    Q <- q
    x <- rbind(P,Q)
    suppressMessages(JSD(x))
}

tic("JS calculation")

if(identical(lr.df$peak, as.character(x.sp@peak))){
jss.df <- data.frame(id=as.character(), 
                     ncell_acc=as.numeric(), 
                     peak=as.character(), 
                     p_value=as.numeric(), 
                     q_value=as.numeric(), 
                     JSD=as.numeric(), 
                     JSS=as.numeric(), 
                     z=as.numeric(), 
                     JSscore=as.numeric(), 
                     attr=as.character(),
                     group=as.character(),
                     zPval=as.numeric(),
                     zQval=as.numeric(),
                     zLocFDR=as.numeric(),
                     zCUT=as.numeric(),
                     JSSeFDRcut=as.numeric(),
                     passJSS=as.character(),
                     JSscoreeFDRcut=as.numeric(),
                     passJSscore=as.character())

for(i in 1:length(attr.lst)){
  attr.sel <- attr.lst[i]
  Q <- rep(0, length(attr.lst))
  Q[i] <- 1
  print(attr.sel)
  jsd <- apply(prop.scale.df, 1, calJSD, q=Q)
  jss <- 1 - sqrt(jsd)
  out <- lr.df
  out$JSD <- jsd
  out$JSS <- jss
  out$z <- ( jss - mean(jss) ) / sd(jss)
  out$JSscore <- jss^2 * prop.norm.df[,i]
  out$attr <- attr.sel
  out$group <- "bgd"
  idx <- which(out$q_value <= sig)
  out[idx, "group"] <- "diff"
  # fit
  JSSeFDRcut <- quantile(out[which(out$group=="bgd"), ][["JSS"]], 1-eFDR)
  JSscoreeFDRcut <- quantile(out[which(out$group=="bgd"), ][["JSscore"]], 1-eFDR)
  fdrOUT <- fdrtool(out$z, plot=FALSE)
  zCUT <- fdrOUT$param[,"cutoff"]
  out$zPval <- fdrOUT$pval
  out$zQval <- fdrOUT$qval
  out$zLocFDR <- fdrOUT$lfdr
  out$zCUT<- zCUT
  out$JSSeFDRcut <- JSSeFDRcut
  out$passJSS <- "no"
  out[which(out$JSS > JSSeFDRcut), "passJSS"] <- "yes"
  out$JSscoreeFDRcut <- JSscoreeFDRcut
  out$passJSscore <- "no"
  out[which(out$JSscore > JSscoreeFDRcut), "passJSscore"] <- "yes"  
  jss.df <- rbind(jss.df, out)
}
}

toc()

# 5. output 
fwrite(jss.df, file=paste(outF, "JSStest.out.tsv", sep="."), quote =F, col.names = T, row.names = F, sep="\t")

pdfname <- paste(outF, "JSStest.fit.pdf", sep=".")
pdf(pdfname)

for(i in 1:length(attr.lst)){
  attr.sel <- attr.lst[i]
  dat <- jss.df[which(jss.df$attr == attr.sel), ]
  zCUT <- unique(dat$zCUT)
  JSSeFDRcut <- unique(dat$JSSeFDRcut)

p1 <- ggplot(dat, aes(x=JSS)) +
  geom_density(aes(fill=group), alpha=0.5) +
  theme_bw() + ggtitle(paste(attr.sel, " (JSS)", sep="")) +
  geom_vline(xintercept = JSSeFDRcut, color="black") +
  theme(axis.text.x=element_text(colour = "black", size = 8),
        axis.text.y=element_text(colour = "black", size = 8),
        axis.title.x=element_text(colour = "black", size = 8),
        axis.title.y=element_text(colour = "black", size = 8))
  
p2 <- ggplot(dat, aes(x=z)) +
  geom_density(fill="grey", alpha=0.5) +
  theme_bw() + ggtitle(paste(attr.sel, " (z-score of JSS)", sep="")) +
  geom_vline(xintercept = zCUT) +
  theme(axis.text.x=element_text(colour = "black", size = 8),
        axis.text.y=element_text(colour = "black", size = 8),
        axis.title.x=element_text(colour = "black", size = 8),
        axis.title.y=element_text(colour = "black", size = 8))
print(p1)
print(p2)
}
dev.off()




