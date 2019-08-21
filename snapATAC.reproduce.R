#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-r", "--ref", required=TRUE, help="load reference for datasets")
parser$add_argument("--metaj", required=TRUE, help="load joint metatable")
parser$add_argument("--metar", required=TRUE, help="load replicate metatable")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()
refF <- args$ref
metajF <- args$metaj
metarF <- args$metar
outF <- args$output

.libPaths("/home/yangli1/R/x86_64-pc-linux-gnu-library/3.4")
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("magrittr")) # need to run every time you start R and want to use %>%
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(require(gridExtra))


# calculate Kullback-Leibler Divergence (Relative Entropy) between rep1 and rep2
datRef <- read.table(refF, header = T, sep="\t")
metaDat <- read.table(metajF, sep="\t", header = T)
names(metaDat)[names(metaDat) == 'x.sp.sample'] <- 'datasets'
names(metaDat)[names(metaDat) == 'x.sp.cluster'] <- 'cluster'
metaDat.m <- merge(datRef, metaDat, by = "datasets")

metaDat.m2 <- metaDat.m %>%
  group_by(cluster, replicate) %>%
  summarise(count=n()) %>%
  mutate(perc=count/sum(count))
outFname <- paste(outF, ".rep.perc.txt", sep="")
write.table(metaDat.m2, outFname, sep="\t", quote=F, col.names=T, row.names=F)

pdf(paste(outF, ".rep.barplot.pdf", sep=""))
p1 <- ggplot(metaDat.m, aes(cluster)) + geom_bar(aes(fill = factor(replicate))) + theme_minimal(base_size = 14)

brks <- c(0, 0.25, 0.5, 0.75, 1)
p2 <- ggplot(metaDat.m2, aes(x = factor(cluster), y = perc, fill = factor(replicate))) +
  geom_bar(stat="identity", width = 0.7) +
  scale_y_continuous(breaks = brks, labels = scales::percent(brks)) +
  labs(x = "cluster", y = NULL, fill = "replicate") +
  theme_minimal(base_size = 14)
grid.arrange(p1, p2, nrow=2)
dev.off()

# load entropy library
suppressPackageStartupMessages(library("entropy", quiet = TRUE))
# probabilities for two random variables
freqs1 = metaDat.m2[which(metaDat.m2$replicate==1), "perc"]$perc
freqs2 = metaDat.m2[which(metaDat.m2$replicate==2), "perc"]$perc

# KL divergence from X1 to X2
kld <- KL.plugin(freqs1, freqs2)


# calculate Rand Index / Adjusted Rand Index
suppressPackageStartupMessages(library("fossil", quiet = TRUE))

metaDat.r <- read.table(metarF, sep="\t", header = T)
names(metaDat.r)[names(metaDat.r) == 'x.sp.sample'] <- 'datasets'
names(metaDat.r)[names(metaDat.r) == 'x.sp.cluster'] <- 'cluster'

metaDat.r.m <- merge(metaDat[,c("datasets", "barcode", "cluster")], metaDat.r[,c("datasets", "barcode", "cluster")], by = c("datasets", "barcode"))

cls.r.x <- metaDat.r.m$cluster.x
cls.r.y <- metaDat.r.m$cluster.y
ari <- adj.rand.index(cls.r.x, cls.r.y)

suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(dendsort))
suppressPackageStartupMessages(library(grid))
#sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))

# confusion matrix and purity
confusion.r <- table(metaDat.r.m[,c("cluster.x", "cluster.y")])

row_max.r <- apply(confusion.r, 1, max)
col_max.r <- apply(confusion.r, 2, max)
purity.r.row <- sum(row_max.r) / sum(confusion.r)
purity.r.col <- sum(col_max.r) / sum(confusion.r)
cat("purity.r: ", purity.r.row, purity.r.col, "\n")

cal_pct <- function(x) {
    return (x  / sum(x))
}

confusion.r.pctrow <- apply(confusion.r, 1, cal_pct)
confusion.r.pctcol <- t(apply(confusion.r, 2, cal_pct))

outFname <- paste(outF, ".rep.confusion.raw.mx", sep="")
write.table(confusion.r, outFname, sep="\t", quote=F, col.names=T, row.names=T)
outFname <- paste(outF, ".rep.confusion.pctbyrow.mx", sep="")
write.table(confusion.r.pctrow, outFname, sep="\t", quote=F, col.names=T, row.names=T)
outFname <- paste(outF, ".rep.confusion.pctbycol.mx", sep="")
write.table(confusion.r.pctcol, outFname, sep="\t", quote=F, col.names=T, row.names=T)


pdf(paste(outF, ".rep.confusion.pdf", sep=""))

setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")
pheatmap(
  mat               = t(confusion.r),
  scale             = 'none',
  color             = viridis(50),
  border_color      = NA,
  cluster_cols      = T,
  cluster_rows      = T,
  show_colnames     = TRUE,
  show_rownames     = TRUE,
  drop_levels       = TRUE,
  fontsize          = 14,
  main              = "confusion matrix: raw"
)
setHook("grid.newpage", NULL, "replace")
grid.text("merged cluster", y=-0.07, gp=gpar(fontsize=14))
grid.text("replicate cluster", x=-0.07, rot=90, gp=gpar(fontsize=14))

setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")
pheatmap(
  mat               = confusion.r.pctrow,
  scale             = 'none',
  color             = viridis(50),
  border_color      = NA,
  cluster_cols      = T,
  cluster_rows      = T,
  show_colnames     = TRUE,
  show_rownames     = TRUE,
  drop_levels       = TRUE,
  fontsize          = 14,
  main              = paste("confusion matrix: normalized by joint cluster\npurity: ", purity.r.row, sep="")
)
setHook("grid.newpage", NULL, "replace")
grid.text("merged cluster", y=-0.07, gp=gpar(fontsize=14))
grid.text("replicate cluster", x=-0.07, rot=90, gp=gpar(fontsize=14))


setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")
pheatmap(
  mat               = confusion.r.pctcol,
  scale             = 'none',
  color             = viridis(50),
  border_color      = NA,
  cluster_cols      = T,
  cluster_rows      = T,
  show_colnames     = TRUE,
  show_rownames     = TRUE,
  drop_levels       = TRUE,
  fontsize          = 14,
  main              = paste("confusion matrix: normalized by replica cluster\npurity: ", purity.r.col, sep="")
)
setHook("grid.newpage", NULL, "replace")
grid.text("merged cluster", y=-0.07, gp=gpar(fontsize=14))
grid.text("replicate cluster", x=-0.07, rot=90, gp=gpar(fontsize=14))
dev.off()

out.sta <- c(outF, kld, ari, purity.r.row, purity.r.col)
names(out.sta) <- c("outF", "kld", "ari", "purity.row", "purity.col")
out.sta <- as.data.frame(t(out.sta))

outFname <- paste(outF, ".rep.sta.txt", sep="")
write.table(out.sta, outFname, sep="\t", quote=F, col.names=T, row.names=F)
