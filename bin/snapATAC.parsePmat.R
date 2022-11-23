#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input", required=TRUE, help="load subset RData")
parser$add_argument("--cpu", default = 4, help="# of cpus [default %(default)s]")
parser$add_argument("-a", "--alpha", default = 1.5, help="alpha for beta distribution")
parser$add_argument("-b", "--beta", default = 2, help="beta for beta distribution")
parser$add_argument("-f", "--frac", default = 0.1, help="fraction of cell have accessiable peaks")
parser$add_argument("-p", "--pep", default = 0.05, help="the desired posterior error probability (PEP)")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

suppressPackageStartupMessages(library("SnapATAC"))
suppressPackageStartupMessages(library("data.table"))
library(tictoc)

RDataF = args$input
a = as.numeric(args$alpha)
b = as.numeric(args$beta)
f = as.numeric(args$frac)
pep = as.numeric(args$pep)
outF = args$output

x.sp <- get(load(RDataF))

CPMnorm <- function(obj){
  pmat <- obj@pmat
  O <- Matrix::colSums(pmat)
  U <- sum(obj@metaData$UQ)
  cpm <- O / U * 10^6
  return(cpm)
}

RAnorm <- function(obj){
  pmat <- obj@pmat
  O <- Matrix::colSums(pmat)
  C <- sum(O)
  RA <- O / C * 10^6
  return(RA)
}

calFrac <- function(obj){
  pmat <- makeBinary(obj, "pmat")@pmat
  F <- Matrix::colSums(pmat) / nrow(pmat)
  return(F)
}

# Return probability that at least half the cells express, if we have observed k of n cells expressing
# Args:
#        k (int):    Number of observed positive cells
#        n (int):    Total number of cells
# Remarks:
#        Probability that at least half the cells express, when we observe k positives among n cells is:
#            p|k,n = 1-(betainc(1+k, 1-k+n, 0.5)*gamma(2+n)/(gamma(1+k)*gamma(1-k+n))/beta(1+k, 1-k+n)
# Note:
#        The formula was derived in Mathematica by computing
#            Probability[x > f, {x \[Distributed] BetaDistribution[1 + k, 1 + n - k]}]
#        and then replacing f by 0.5

p_half <- function(k, n, f=0.5){
  a = 1.5
  b = 2
  incb = pbeta(f, a + k, b - k + n)
  if (incb == 0){
    p = 1.0
  }else{
    p = 1L - exp(log(incb) + lbeta(a + k, b - k + n) + lgamma(a + b + n) - lgamma(a + k) - lgamma(b - k + n))
  }
  return(p)
}

# Trinarize a vector, grouped by labels, using a beta binomial model
#    Args:
#        array (ndarray of ints):    The input vector of ints
#        labels (ndarray of ints):    Group labels 0, 1, 2, ....
#        pep (float):                The desired posterior error probability (PEP)
#    Returns:
#        ps (ndarray of float):        The posterior probability of expression in at least a fraction f
#        expr_by_label (ndarray of float): The trinarized expression pattern (one per label)
#    Remarks:
#        We calculate probability p that at least half the cells express (in each group),
#        and compare with pep, setting the binary pattern to 1 if p > pep,
#        -1 if p < (1 - pep) and 0 otherwise.
betabinomial_trinarize <- function(k, n, f, pep){
  ps = as.numeric(lapply(k, p_half, n, f))
  access_by_label <- rep(0,length(ps))
  access_by_label[which(ps > (1-pep))] <- 1
  access_by_label[which(ps < pep)] <- -1
  pval <- 1-ps
  out <- list(prob=ps, pval=pval, access_label=access_by_label)
  return(out)
}


tic("calCnt")
CT <- Matrix::colSums(x.sp@pmat)
CNT1 <- Matrix::colSums(x.sp@pmat[which(x.sp@metaData$replicate == 1), ]) 
CNT2 <- Matrix::colSums(x.sp@pmat[which(x.sp@metaData$replicate == 2), ])
tic()

tic("calFrac")
Frac <- calFrac(x.sp)
toc()

tic("CPMnorm")
CPM <- CPMnorm(x.sp)
toc()

tic("RAnorm")
RA <- RAnorm(x.sp.sub)
toc()

tic("Trinarize")
pmat <- makeBinary(x.sp.sub, "pmat")@pmat
k <- Matrix::colSums(pmat)
n <- nrow(pmat)
tri_out <- betabinomial_trinarize(k, n, f, pep)
toc()
accessFDR <- p.adjust(tri_out$pval, method="BH")

message(paste("write norm score for: ", RDataF, sep=""));
out <- data.frame(peak=x.sp@peak$name, CT=CT, Cnt1=CNT1, Cnt2=CNT2, Frac=Frac, CPMnorm=CPM, RAnorm=RA, accessProb=tri_out$prob, accessPval=tri_out$pval, accessFDR=accessFDR, accessLabel=tri_out$access_label)
outfname = paste(outF, ".parsePmat.txt",sep="")
write.table(out, outfname, quote=F, col.names=T, row.names=F, sep="\t")







