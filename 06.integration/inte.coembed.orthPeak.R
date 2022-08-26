#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-r", "--ref", required=TRUE, help="load rna RDS")
parser$add_argument("-q", "--quy", required=TRUE, help="load atac RDS")
parser$add_argument("--refAttr", required=TRUE, help="used ref attr")
parser$add_argument("--quyAttr", required=TRUE, help="used quy attr")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

suppressPackageStartupMessages(library("SnapATAC"))
suppressPackageStartupMessages(library("GenomicRanges"));
suppressPackageStartupMessages(library("RColorBrewer"))
suppressPackageStartupMessages(library("viridis"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("reshape2"))
suppressPackageStartupMessages(library("pheatmap"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("plyr"))
library("tictoc")

refF <- args$ref
quyF <- args$quy
refAttr <- args$refAttr
quyAttr <- args$quyAttr
outF = args$output

#-------------------
# load data
mm.sp <- readRDS(refF)
mm.sp@metaData$species <- "mm"
hg.sp <- readRDS(quyF)
hg.sp@metaData$species <- "hg"

# reassign bmat
identical(mm.sp@peak$name, hg.sp@peak$name)

mm.sp@bmat <- mm.sp@pmat
mm.sp@feature <- mm.sp@peak

hg.sp@bmat <- hg.sp@pmat
hg.sp@feature <- hg.sp@peak

# reassign metaData
mm.sp@metaData$ref <- mm.sp@metaData[, refAttr]
mm.sp@metaData$quy <- NA
hg.sp@metaData$ref <- NA
hg.sp@metaData$quy <- hg.sp@metaData[, quyAttr]


colnames(mm.sp@metaData)[which(colnames(mm.sp@metaData)=="tsse")] <- "TSSe"
colnames(hg.sp@metaData)[which(colnames(hg.sp@metaData)=="log10UQ")] <- "logUMI"
used.colnames <- intersect(colnames(mm.sp@metaData), colnames(hg.sp@metaData))

mm.sp@metaData <- mm.sp@metaData[,used.colnames] 
hg.sp@metaData <- hg.sp@metaData[,used.colnames] 

# merge
x.sp.ls <- c(mm.sp, hg.sp)
x.sp = Reduce(snapRbind, x.sp.ls);

# cluster on ortholog pmat 
x.sp <- makeBinary(x.sp, mat="bmat")

# remove 0 columns
#idy <- which(colSums(x.sp@bmat)==0)
#x.sp@bmat <- x.sp@bmat[, -idy]
#x.sp@feature <- x.sp@feature[-idy]


tic("runDiffusionMaps")
x.sp = runDiffusionMaps(
    obj= x.sp,
    input.mat="bmat",
    num.eigs=50
  );
toc()

pdf(paste(outF, ".pre.plotDimResuce.pdf",sep=""))
plotDimReductElbow(
    obj=x.sp,
    point.size=1.5,
    point.shape=19,
    point.color="red",
    point.alpha=1,
    pdf.file.name=NULL,
    pdf.height=7,
    pdf.width=7,
    labs.title="PCA Elbow plot",
    labs.subtitle=NULL
    );

plotDimReductPW(
    obj=x.sp,
    eigs.dims=1:50
    );
dev.off()

# correct species
iniDim <- 4
dims <- 20
#vars <- "species"
#x.sp = runHarmony(
#    obj=x.sp, 
#    eigs.dim=1:dims, 
#    meta_data=x.sp@metaData,
#    vars_use = vars
#  )

# KNN
tic("runKNN")
x.sp = runKNN(
    obj=x.sp,
    eigs.dims=iniDim:dims, # dim 1&2 are species difference
    k=50
    );
toc()


## Visulization
tic("runViz_umap")
x.sp = runViz(
    obj=x.sp,
    tmp.folder=tempdir(),
    dims=2,
    eigs.dims=iniDim:dims,
    method="umap",
    seed.use=10
  );
toc()

## R-igraph
tic("runCluster_R-igraph")
x.sp = runCluster(obj=x.sp,
    tmp.folder=tempdir(),
    louvain.lib="R-igraph",
    seed.use=10
    );
toc()


pdf(paste(outF, ".coembed.cluster.pdf",sep=""))
plotViz(
    obj= x.sp,
    method="umap",
    main="Cluster",
    point.color=x.sp@cluster,
    point.size=0.2,
    point.shape=19,
    text.add=FALSE,
    text.size=1,
    text.color="black",
    down.sample=10000,
    legend.add=FALSE
  );

plotViz(
    obj= x.sp,
    method="umap",
    main="Cluster",
    point.color=x.sp@cluster,
    point.size=0.2,
    point.shape=19,
    text.add=TRUE,
    text.size=1,
    text.color="black",
    down.sample=10000,
    legend.add=FALSE
  );

plotViz(
    obj= x.sp,
    method="umap",
    main="mouse",
    point.color=x.sp@metaData$ref,
    point.size=0.2,
    point.shape=19,
    text.add=TRUE,
    text.size=1,
    text.color="black",
    down.sample=10000,
    legend.add=FALSE
  );

plotViz(
    obj= x.sp,
    method="umap",
    main="human",
    point.color=x.sp@metaData$quy,
    point.size=0.2,
    point.shape=19,
    text.add=TRUE,
    text.size=1,
    text.color="black",
    down.sample=10000,
    legend.add=FALSE
  );

plotViz(
    obj= x.sp,
    method="umap",
    main="species",
    point.size=0.2,
    point.shape=19,
    point.color=x.sp@metaData[,"species"],
    text.add=FALSE,
    text.size=1.5,
    text.color="black",
    down.sample=10000,
    legend.add=TRUE
  );
dev.off()

x.sp@metaData$cluster <- x.sp@cluster
save(x.sp, file=paste(outF, "coembed.RData", sep="."))

metaout <- x.sp@metaData
fwrite(metaout, file=paste(outF, "coembed.meta.txt", sep="."), sep="\t", quote = F, col.names = T, row.names = F)


#------------------------
# overlap score
#metaout <- x.sp@metaData

#-----------------------------------------------------
# t1: table with 2 columns: coembed labels, raw labels
# t2: table with 2 columns: coembed labels, raw labels
cal_ovlpScore <- function(t1, t2){
  t1.table <- table(t1)
  t2.table <- table(t2)
  t1.pct <- apply(t1.table, 2, function(x){x/sum(x)})
  t2.pct <- apply(t2.table, 2, function(x){x/sum(x)})
  t1.labels <- colnames(t1.pct)
  t2.labels <- colnames(t2.pct)
  ovlpScore.df <- data.frame(anno1=as.character(), anno2=as.character(), ovlpScore=as.numeric())
  for(t1.label in t1.labels){
    for(t2.label in t2.labels){
      t1.pct.df <- data.frame(t1.pct[,t1.label])
      colnames(t1.pct.df) <- "t1"
      t1.pct.df$ident <- rownames(t1.pct.df)
      t2.pct.df <- data.frame(t2.pct[,t2.label])
      colnames(t2.pct.df) <- "t2"
      t2.pct.df$ident <- rownames(t2.pct.df)
      comp.df <- join(t1.pct.df, t2.pct.df, by="ident", type="full")
      comp.df[is.na(comp.df)] <- 0
      comp.df$ident <- NULL
      comp.df <- t(comp.df)
      ovlpScore <- sum(apply(comp.df, 2, min))
      out <- data.frame(anno1=t1.label, anno2=t2.label, ovlpScore=ovlpScore)
      ovlpScore.df <- rbind(ovlpScore.df, out)
    }
  }
  return(ovlpScore.df)
}


# calculate overlap
ident2ref <- data.frame(idents=metaout$cluster, rna_label=metaout$ref)
ident2ref <- ident2ref[complete.cases(ident2ref), ]

ident2query <- data.frame(idents=metaout$cluster, atac_label=metaout$quy)
ident2query <- ident2query[complete.cases(ident2query), ]

ovlpScore.df <- cal_ovlpScore(ident2ref, ident2query)
ovlpScore.df <- as.data.frame(ovlpScore.df)
write.table(ovlpScore.df, paste(outF,"coembed.overlapScore.tsv", sep="."), sep="\t", col.names = T, row.names = F, quote = F)

# ovlpScore.df.sel <- subset(ovlpScore.df, ovlpScore.df$ovlpScore>=0.2)
ovlpScore.mx <- reshape2::dcast(ovlpScore.df, anno1~anno2, value.var="ovlpScore", fill = 0)
ovlpScore.mx <- as.data.frame(ovlpScore.mx)
write.table(ovlpScore.mx, paste(outF,"coembed.overlapScore.mx.tsv", sep="."), sep="\t", col.names = T, row.names = F, quote = F)

ovlpScore.plot <- ovlpScore.mx
rownames(ovlpScore.plot) <- ovlpScore.plot$anno1
ovlpScore.plot$anno1 <- NULL
ovlpScore.plot <- as.matrix(ovlpScore.plot)

pdf(paste(outF, ".coembed.ovlpScore.pdf",sep=""))
pheatmap(ovlpScore.plot, main = "overlap score: row: mouse; col: human")
dev.off()

