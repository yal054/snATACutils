#!/usr/bin/env Rscript

# use new version
suppressPackageStartupMessages(library("SnapATAC"))
suppressPackageStartupMessages(library("GenomicRanges"));
suppressPackageStartupMessages(library("foreach"))
suppressPackageStartupMessages(library("doParallel"))
library("umap")
library("leiden")
library("tictoc")

inputF <- "rs1atac.datasets.list"
#down_rate <- 0.25
fragment_num <- 1000
umi_num <- 500
mito_ratio <- 1
dup_ratio <- 1
umap_ratio <- 0
pair_ratio <- 0
bin_size <- 5000
doubletsF <- "rs1atac.pmat.doubelts.txt"
black_list <- "/projects/ps-renlab/yangli/genome/mm10/mm10.blacklist.bed.gz"
pc_num <- 50
cpus <- 5
outF <- "rs1atac"
tsse2depthF <- "rs1atac.tsse2depth.txt"
seed.use <- 2019
dims <- 20
geneF <- "gene.markers.txt"
resolution <- 0.5

inputf <- read.table(inputF,sep="\t",header=F)
file.list <- as.character(inputf[,2])
sample.list <- as.character(inputf[,1])

# function for downsample datasets
downSampleSnap <- function(obj, down_rate, seed.use) {
  UseMethod("downSampleSnap", obj);
}

downSampleSnap.default <- function(obj, down_rate=0.3, seed.use=10){
    nBarcode = length(obj@barcode);
    set.seed(seed.use);
    selected_barcode = sample(obj@barcode, ceiling(nBarcode * down_rate));
    selected_idx = which(obj@barcode %in% selected_barcode);
    obj@barcode = obj@barcode[selected_idx];
    obj@file = obj@file[selected_idx];
    obj@sample = obj@sample[selected_idx];
    obj@metaData = obj@metaData[selected_idx,];
    return(obj);
}

filterDoublets <- function(obj, doublets) {
  UseMethod("filterDoublets", obj);
}

filterDoublets.default <- function(obj, doublets=doubletsF){
    doublets = read.table(doubletsF, header=F, sep="\t")
    selected_idx <- which(!paste(obj@sample, obj@barcode,sep=".") %in% doublets$V1)
    obj@barcode = obj@barcode[selected_idx];
    obj@file = obj@file[selected_idx];
    obj@sample = obj@sample[selected_idx];
    obj@metaData = obj@metaData[selected_idx,];
    return(obj);
}


tic("createSnap")
x.sp.ls = lapply(seq(file.list), function(i){
    x.sp = createSnap(file=file.list[i], sample=sample.list[i], do.par = TRUE, num.cores=cpus);
    x.sp
  })
names(x.sp.ls) = sample.list;
x.sp = Reduce(snapRbind, x.sp.ls);
x.sp@metaData["sample"] = x.sp@sample;
toc()

rm(x.sp.ls)

outfname = paste(outF, ".raw.RData",sep="")
save(x.sp, file=outfname)

outmetaf <- paste(outF, ".raw.meta.txt", sep="")
outmetamx <- x.sp@metaData
write.table(outmetamx, outmetaf, row.names=F, col.names=T, sep="\t", quote=F)


if(!is.null(doubletsF)){
message("Filter out potential doublets: ", doubletsF);
x.sp <- filterDoublets(x.sp, doubletsF)
}

tsse2depth <- read.table(tsse2depthF, sep="\t", header=T)
tsse2depth$log10UMI <- log(tsse2depth$depth+1, 10);
tsse2depth <- tsse2depth[which(paste(tsse2depth$sample, tsse2depth$barcode, sep=".") %in% paste(x.sp@sample, x.sp@barcode, sep=".")), ]


library(viridisLite);
library(ggplot2);
p = ggplot(
    tsse2depth[sample(1:nrow(tsse2depth), 20000),],
    aes(x= log10UMI, y= tsse)) +
    geom_point(size=0.1, col="grey") +
    theme_bw() +
    ggtitle("RS1 atac-seq 201908") +
    ylim(0, 50) + xlim(0, 6) +
    geom_hline(yintercept=10) + geom_vline(xintercept=3) +
    labs(x = "log10(UMI)", y="tss enrichment")
pdf(paste(outF, ".tsse2depth.raw.pdf",sep=""))
plot(p)
dev.off()

tsse2depth.sel = tsse2depth[which(tsse2depth$depth >= 3 & tsse2depth$tsse >= 10),c(1,2,5,6)];
tsse2depth.sel.idx <- which(paste(x.sp@sample, x.sp@barcode,sep=".") %in% paste(tsse2depth.sel$sample,tsse2depth.sel$barcode, sep="."))
x.sp = x.sp[tsse2depth.sel.idx,];

pdf(paste(outF, ".filter.sta.pdf",sep=""))
plotBarcode(x.sp,
              pdf.file.name=NULL,
              pdf.width=7,
              pdf.height=7,
              col="grey",
              border="grey",
              breaks=50
              );
dev.off()

outfname = paste(outF, ".filter.RData",sep="")
save(x.sp, file=outfname)

outmetaf <- paste(outF, ".filter.meta.txt", sep="")
outmetamx <- cbind(x.sp@sample, x.sp@metaData)
colnames(outmetamx)[which(colnames(outmetamx)=="x.sp@sample")] <- "sample"
suppressPackageStartupMessages(library("plyr"))
outmetamx <- join(outmetamx, tsse2depth[,c(1,2,5,6,7)], by=c("sample", "barcode"))
write.table(outmetamx, outmetaf, row.names=F, col.names=T, sep="\t", quote=F)


## 0. sampling
sampleSize <- 10000

tic("sampling landmark")
x.sp@metaData$logUMI <- log10(x.sp@metaData$UQ + 1)

row.covs.dens <- density(
    x = x.sp@metaData[,"logUMI"],
    bw = 'nrd', adjust = 1
  );

sampling_prob <- 1 / (approx(x = row.covs.dens$x, y = row.covs.dens$y, xout = x.sp@metaData[,"logUMI"])$y + .Machine$double.eps);
set.seed(2019);
idx.landmark.ds <- sort(sample(x = seq(nrow(x.sp)), size = sampleSize, prob = sampling_prob));
x.landmark.sp = x.sp[idx.landmark.ds,];
x.query.sp = x.sp[-idx.landmark.ds,];
toc()


## 1. identify usable features
tic("addBmatToSnap")
x.landmark.sp = addBmatToSnap(x.landmark.sp, bin.size=5000);
toc()

tic("makeBinary")
x.landmark.sp = makeBinary(x.landmark.sp, mat="bmat");
toc()

tic("filterBins")
black_list = read.table(black_list);
black_list.gr = GRanges(
    black_list[,1], 
    IRanges(black_list[,2], black_list[,3])
  );
idy = queryHits(
    findOverlaps(x.sp@feature, black_list.gr)
  );
if(length(idy) > 0){
    x.landmark.sp = x.landmark.sp[,-idy, mat="bmat"];
  };
x.sp

chr.exclude = seqlevels(x.landmark.sp@feature)[grep("random|chrM", seqlevels(x.landmark.sp@feature))];
idy = grep(paste(chr.exclude, collapse="|"), x.landmark.sp@feature);
if(length(idy) > 0){
    x.landmark.sp = x.landmark.sp[,-idy, mat="bmat"]
  };
x.sp

bin.cov = log10(Matrix::colSums(x.landmark.sp@bmat)+1);
pdf(paste(outF, ".BinCoverage.pdf",sep=""))
hist(
    bin.cov[bin.cov > 0], 
    xlab="log10(bin cov)", 
    main="log10(Bin Cov)", 
    col="lightblue", 
    xlim=c(0, 5)
  );
dev.off()

bin.cutoff = quantile(bin.cov[bin.cov > 0], 0.95);
idy = which(bin.cov <= bin.cutoff & bin.cov > 0);
x.landmark.sp = x.landmark.sp[, idy, mat="bmat"];
x.landmark.sp

idx = which(Matrix::rowSums(x.landmark.sp@bmat) > 1000);
x.landmark.sp = x.landmark.sp[idx,];
x.landmark.sp
toc()

#summarySnap(x.landmark.sp)
x.landmark.sp
outfname = paste(outF, ".landmark.RData",sep="")
save(x.landmark.sp, file=outfname)

#summarySnap(x.query.sp)
x.query.sp
outfname = paste(outF, ".query.RData",sep="")
save(x.query.sp, file=outfname)


# 2. embedding
tic("runDiffusionMaps")
x.landmark.sp = runDiffusionMaps(
    obj= x.landmark.sp,
    input.mat="bmat",
    num.eigs=50
  );
x.landmark.sp@metaData$landmark = 1;
toc()


# 3.extension
tic("extension")
num.files = length(sample.list);

registerDoParallel(cpus)
#x.query.ls <- foreach (i=1:num.files, .combine=snapRbind, .inorder=TRUE) %dopar% {
x.query.ls <- foreach (i=1:num.files, .inorder=TRUE) %dopar% {
print(sample.list[i])
x.query.sub <- x.query.sp[which(x.query.sp@sample == as.character(sample.list[i])), ]
x.query.sub = addBmatToSnap(x.query.sub, do.par=FALSE, num.cores=cpus, bin.size=5000);
x.query.sub = makeBinary(x.query.sub);

idy = unique(queryHits(findOverlaps(x.query.sub@feature, x.landmark.sp@feature)));
x.query.sub = x.query.sub[,idy, mat="bmat"];
x.query.sub = runDiffusionMapsExtension(
    obj1=x.landmark.sp,
    obj2=x.query.sub,
    input.mat="bmat"
  );
x.query.sub@metaData$landmark = 0;
x.query.sub = rmBmatFromSnap(x.query.sub)
return(x.query.sub)
}
closeAllConnections()

x.query.sp = Reduce(snapRbind, x.query.ls);
x.landmark.sp = rmBmatFromSnap(x.landmark.sp)
x.sp = snapRbind(x.landmark.sp, x.query.sp);
x.sp = x.sp[order(x.sp@metaData[,"sample"])];
toc()

rm(x.landmark.sp)
rm(x.query.sp)

outfname = paste(outF, ".pre.RData",sep="")
save(x.sp, file=outfname)

outmeta = paste(outF, ".pre.meta.txt",sep="")
write.table(x.sp@metaData, file=outmeta, sep="\t", quote=F, col.names=T, row.names=F)


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
    labs.subtitle=NULL,
    );

plotDimReductPW(
    obj=x.sp,
    eigs.dims=1:pc_num
    );
dev.off()


# KNN
tic("runKNN")
x.sp = runKNN(
    obj=x.sp,
    eigs.dims=1:dims,
    k=15,
    );
toc()


## Clustering

#library("leiden")
#tic("runCluster_leiden")
#x.sp = runCluster(
#    obj=x.sp,
#    tmp.folder=tempdir(),
#    louvain.lib="leiden",
#    resolution = resolution,
#    seed.use=10
#    );
#toc()

## R-igraph
tic("runCluster_R-igraph")
x.sp = runCluster(obj=x.sp, 
    tmp.folder=tempdir(), 
    louvain.lib="R-igraph",
    seed.use=10
    );
toc()

## Visulization

tic("runViz_umap")
x.sp = runViz(
    obj=x.sp, 
    tmp.folder=tempdir(),
    dims=2,
    eigs.dims=1:dims, 
    method="umap",
    seed.use=10
  );
toc()

tic("runViz_tsne")
x.sp = runViz(
    obj=x.sp,
    tmp.folder=tempdir(),
    dims=2,
    eigs.dims=1:dims,
    method="Rtsne",
    seed.use=10
  );
toc()


pdf(paste(outF, ".cluster.pdf",sep=""))
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
    method="tsne",
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
    method="tsne",
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

plotFeatureSingle(
    obj=x.sp,
    feature.value=x.sp@metaData[,"logUMI"],
    method="umap", 
    main="Read Depth",
    point.size=0.2, 
    point.shape=19, 
    down.sample=10000,
    quantiles=c(0.01, 0.99)
  );

plotViz(
    obj= x.sp,
    method="umap", 
    main="Landmark",
    point.size=0.2, 
    point.shape=19, 
    point.color=x.sp@metaData[,"landmark"], 
    text.add=FALSE,
    text.size=1.5,
    text.color="black",
    down.sample=10000,
    legend.add=TRUE
  );
dev.off()

## Heretical clustering of the clusters
# calculate the ensemble signals for each cluster
#ensemble.ls = lapply(split(seq(length(x.sp@cluster)), x.sp@cluster), function(x){
#    Matrix::colMeans(x.sp@bmat[x,])
#    })

# cluster using 1-cor as distance
#hc = hclust(as.dist(1 - cor(t(do.call(rbind, ensemble.ls)))), method="ward.D2");
#par(mfrow=c(1,1))
#pdf(paste(outF, ".hc.pdf",sep=""))
#plot(hc, hang=-1, xlab="");
#dev.off()

## Create chromatin lanscape and identify cis-elements for each cluster seperately.

outmetaf <- paste(outF, ".cluster.meta.txt", sep="")
outmetamx <- cbind(x.sp@sample, x.sp@metaData, x.sp@cluster, x.sp@umap)
write.table(outmetamx, outmetaf, row.names=F, col.names=T, sep="\t", quote=F)

outfname = paste(outF, ".cluster.RData",sep="")
save(x.sp, file=outfname)

graph.knn <- x.sp@graph@mat

outfname = paste(outF, ".knn.mmtx", sep="")
writeMM(graph.knn,file=outfname)


