#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input", required=TRUE, help="load clustered RData")
parser$add_argument("-a", "--attr", required=TRUE, help="attribute's column name used for regression: MajorRegion | SubRegion | DissectionRegion | SubType")
parser$add_argument("-l", "--lrTest", required=TRUE, help="LR test results")
parser$add_argument("-s", "--sig", default = 0.001, help="significant cutoff for LR test [default %(default)s]")
parser$add_argument("-e", "--eFDR", default = 0.05, help="empirical FDR cutoff [default %(default)s]")
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
suppressPackageStartupMessages(library("mixtools"))
#suppressPackageStartupMessages(library("fitdistrplus"))
library("tictoc")

RDataF = args$input
lrF = args$lrTest
attr = args$attr
sig = as.numeric(args$sig)
eFDR = as.numeric(args$eFDR)
outF = args$output

#---------
# init
x.sp <- get(load(RDataF))

lr.df <- fread(lrF, sep="\t", header=T)

#-----------------------
# Jensen-Shannon divergence
#----------------
# fit mixture
find.cutoff <- function(model, proba=0.5, i=index.lower) {
    ## Cutoff such that Pr[drawn from bad component] == proba
    f <- function(x) {
        proba - (model$lambda[i]*dnorm(x, model$mu[i], model$sigma[i]) /
                     (model$lambda[1]*dnorm(x, model$mu[1], model$sigma[1]) + model$lambda[2]*dnorm(x, model$mu[2], model$sigma[2])))
        }
        return(uniroot(f=f, lower = min(model$mu[1],model$mu[2]), upper = max(model$mu[1], model$mu[2]))$root)  # Careful with division by zero if changing lower and upper
}

# 1. proprotion for each cluster
attr.lst <- unique(x.sp@metaData[, attr])

tic("get/norm prop")
prop.lst <- list()
for(t in attr.lst){
#  print(t)
  idx <- which(x.sp@metaData[, attr] == t)
#  n <- length(idx)
  prop <- Matrix::colSums(x.sp@gmat[idx, ])
  prop <- prop / sum(prop) * 10^4
  prop.lst[[t]] <- prop
}

prop.df <- do.call(cbind, prop.lst)

#gene.sel.lst <- c()
#for(t in attr.lst){
#  print(t)
#  log10prop <- log10(prop.df[,t])
#  log10prop <- log10prop[!is.infinite(log10prop)]  
#  mixmdl = normalmixEM(log10prop, k=2)
#  index.lower <- which.min(mixmdl$mu)
#  cutoffs <- find.cutoff(mixmdl, proba=0.5, i=index.lower)
#  idx.sel <- which(log10(prop.df[,t]) >= cutoffs)
#  gene.sel <- rownames(prop.df)[idx.sel]
#  gene.sel.lst <- c(gene.sel.lst, gene.sel)
#}
#gene.sel.lst <- unique(gene.sel.lst)


# fit cpm
cnt <- Matrix::colSums(x.sp@gmat)
cpm <- cnt / sum(x.sp@metaData$UQ) * 10^6
log10cpm <- log10(cpm+1)
#plot(density(log10cpm))
mixmdl <- normalmixEM(log10cpm, k=2)
index.lower <- which.min(mixmdl$mu)
cutoffs <- find.cutoff(mixmdl, proba=0.5, i=index.lower)
idx.sel <- which(log10(cpm+1) >= cutoffs)
gene.sel.lst <- rownames(prop.df)[idx.sel]


idx <- match(gene.sel.lst, lr.df$gene_short_name)
lr.df <- lr.df[idx, ]

# 2. scale factor
medacc.lst <- list()
for(t in attr.lst){
#  print(t)
  idx <- which(x.sp@metaData[, attr] == t)
  medacc <- Matrix::rowSums(x.sp@gmat[idx, ])
  medacc <- median(medacc)
  medacc.lst[[t]] <- medacc
}

medacc.df <- do.call(cbind, medacc.lst)
scalef.df <- mean(log10(medacc.df[1,])) / log10(medacc.df)

# 3. correct/norm and scale
prop.norm.df <- sweep(prop.df, 2, scalef.df[1,], "*")
idx <- match(gene.sel.lst, rownames(prop.norm.df))
prop.norm.df <- prop.norm.df[idx, ]

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

if(identical(lr.df$gene_short_name, as.character(rownames(prop.scale.df)))){
jss.df <- data.frame(gene_short_name=as.character(),
                     num_cells_expressed=as.numeric(),
                     gene_id=as.character(),
                     p_value=as.numeric(),
                     q_value=as.numeric(),
                     id = as.character(),
                     type = as.character(),
                     JSD=as.numeric(),
                     JSS=as.numeric(),
                     z=as.numeric(),
                     normCPM = as.numeric(),
                     JSscore=as.numeric(),
                     attr=as.character(),
                     group=as.character(),
                     zPval=as.numeric(),
                     zQval=as.numeric(),
                     zLocFDR=as.numeric(),
                     zCUT=as.numeric(),
                     JSSeFDRcut=as.numeric(),
                     passJSS=as.character())

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
  out$normCPM <- prop.norm.df[,i]
  out$JSscore <- jss^2 * prop.norm.df[,i]
  out$attr <- attr.sel
  out$group <- "bgd"
  idx <- which(out$q_value <= sig)
  out[idx, "group"] <- "diff"
  # fit
  fdrOUT <- fdrtool(out$z, plot=FALSE)
  zCUT <- fdrOUT$param[,"cutoff"]
  out$zPval <- fdrOUT$pval
  out$zQval <- fdrOUT$qval
  out$zLocFDR <- fdrOUT$lfdr
  out$zCUT<- zCUT
  JSSeFDRcut <- quantile(out[which(out$group=="bgd"), ][["JSS"]], 1-eFDR)
  out$JSSeFDRcut <- JSSeFDRcut
  out$passJSS <- "no"
  out[which(out$JSS > JSSeFDRcut), "passJSS"] <- "yes"
  jss.df <- rbind(jss.df, out)
}
}

toc()

# 5. output 
fwrite(jss.df, file=paste(outF, "JSStest.out.tsv", sep="."), quote =F, col.names = T, row.names = F, sep="\t")

pdfname <- paste(outF, "JSStest.fit.pdf", sep=".")
pdf(pdfname)
plot(density(log10cpm))
abline(v=cutoffs)
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




