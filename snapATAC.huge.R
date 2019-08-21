#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input", required=TRUE, help="input txt with the information of path of snap file and name")
parser$add_argument("--doublets", default="/projects/ps-renlab/yangli/projects/CEMBA/01.merged_dat/doublet.txt", help="input txt with doublets")
parser$add_argument("--fragment_num", default = 1000, help="cutoff of fragment.num [default %(default)s]")
parser$add_argument("--umi_num", default = 500, help="cutoff of UMI.num [default %(default)s]")
parser$add_argument("--mito_ratio", default = 1, help="cutoff of mito.ratio [default %(default)s]")
parser$add_argument("--dup_ratio", default = 1, help="cutoff of dup.ratio [default %(default)s]")
parser$add_argument("--umap_ratio", default = 0, help="cutoff of umap.ratio [default %(default)s]")
parser$add_argument("--pair_ratio", default = 0, help="cutoff of pair.ratio [default %(default)s]")
parser$add_argument("--bin_size", default = 5000, help="binSize to use [default %(default)s]")
parser$add_argument("--black_list", required=TRUE, help="black list file")
parser$add_argument("--pc_num", default = 50, help="num of PCA components [default %(default)s]")
parser$add_argument("--cpu", default = 4, help="# of cpus [default %(default)s]")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

suppressPackageStartupMessages(library("SnapATAC"))
suppressPackageStartupMessages(library("tictoc"))
suppressPackageStartupMessages(library("parallel"))
suppressPackageStartupMessages(library("GenomicRanges"))
suppressPackageStartupMessages(library("data.table"))

inputF = args$input
fragment_num = as.numeric(args$fragment_num)
umi_num = as.numeric(args$umi_num)
mito_ratio = as.numeric(args$mito_ratio)
dup_ratio = as.numeric(args$dup_ratio)
umap_ratio = as.numeric(args$umap_ratio)
pair_ratio = as.numeric(args$pair_ratio)
bin_size = as.numeric(args$bin_size)
black_list = args$black_list
doubletsF = args$doublets
pc_num = as.numeric(args$pc_num)
cpus = as.numeric(args$cpu)
outF = args$output


runJaccard2 <- function(
    obj1,
    obj2
){
    calJaccard <- function(X_i, X_j){
        A = Matrix::tcrossprod(X_i, X_j);
        bi = Matrix::rowSums(X_i);
        bj = Matrix::rowSums(X_j);
        jmat = as.matrix(A / (replicate(ncol(A), bi) + t(replicate(nrow(A), bj)) - A));
        rm(A);
        rm(bi);
        rm(bj);
        gc();    
        return(jmat);                
    }
    
    jmat = calJaccard(obj1@bmat, obj2@bmat);

    # remove large objects
    obj1@jmat@jmat = jmat;
    obj1@jmat@p1 = Matrix::rowMeans(obj1@bmat);
    obj1@jmat@p2 = Matrix::rowMeans(obj2@bmat);
    obj1@jmat@norm = FALSE;
    obj1@jmat@method = character();
    rm(jmat);
    gc();
    return(obj1);
}


# 1) collect the coverage from all samples
# 2) randomly select cells based on the coverage
# 3) upload the 5000 cells as anchor cells
# 4) identify usable features
# 5) for loop for each sample

# 1) collect the coverage from all samples
# load file list
inputf <- read.table(inputF,sep="\t",header=F)
fnames <- as.character(inputf[,2])
pnames <- as.character(inputf[,1])
num.files = nrow(inputf)

# load doublets
doublets = read.table(doubletsF, header=T, sep="\t")

tic("createSnap")
x.sp = createSnap(
    file=fnames,
    sample=pnames,
    num.cores=cpus
    );
toc()

tic("filterCells")
x.sp = filterCells(x.sp,
                   subset.names=c("fragment.num",
                                  "UMI",
                                  "mito.ratio",
                                  "dup.ratio",
                                  "umap.ratio",
                                  "pair.ratio"),
                   low.thresholds=c(fragment_num, umi_num, 0, 0, umap_ratio, pair_ratio),
                   high.thresholds=c(Inf, Inf, mito_ratio, dup_ratio, 1, 1)
                  );
toc()


singletons <- subset(doublets, doublets$predicted_doublets==FALSE)
barcodes.sel <- which(paste(x.sp@sample,x.sp@barcode,sep=".") %in% paste(singletons$sample, singletons$barcode, sep="."))

x.sp <- x.sp[barcodes.sel,]


# 2) randomly select cells based on the coverage
row.covs = log(x.sp@metaData$UQ+1,10);        
row.covs.dens <- density(x = row.covs, bw = 'nrd', adjust = 1)
sampling_prob <- 1 / (approx(x = row.covs.dens$x, y = row.covs.dens$y, xout = row.covs)$y + .Machine$double.eps)
set.seed(2019);
idx <- sort(sample(x = seq(row.covs), size = 5000, prob = sampling_prob));

# 3) upload the 5000 cells as anchor cells
x.anchor.sp = x.sp[idx,];
x.anchor.sp = addBmatToSnap(x.anchor.sp, do.par=FALSE, num.cores=1, bin.size=5000);
x.anchor.sp = makeBinary(x.anchor.sp);

# 4) identify usable features
bin.cov = log(Matrix::colSums(x.anchor.sp@bmat)+1, 10)
bin.cutoff = quantile(bin.cov[bin.cov > 0], 0.95)
idy1 = which(bin.cov >= bin.cutoff | bin.cov == 0)

black_list = read.table(black_list);
black_list.gr = GRanges(
    black_list[,1], 
    IRanges(black_list[,2], black_list[,3])
);

idy2 = queryHits(findOverlaps(x.anchor.sp@feature, black_list.gr));
idy = unique(c(idy1, idy2));
features.gr = x.anchor.sp@feature[-idy];

idy = unique(queryHits(findOverlaps(x.anchor.sp@feature, features.gr)));
x.anchor.sp = x.anchor.sp[,idy, mat="bmat"];

outfname = paste(outF, ".anchor.RData",sep="")
save(x.anchor.sp, file=outfname)


# 5)  for loop for each sample --> jaccard matrix
nmat.ls = list();
for(fname in fnames){
    print(fname);
    x.sub.sp = createSnap(fname, sample=as.character(fname));
    x.sub.sp = x.sub.sp[which(x.sub.sp@barcode %in% x.sp@barcode[x.sp@file == fname]),];
    x.sub.sp = addBmatToSnap(x.sub.sp, bin.size=5000);
    x.sub.sp = makeBinary(x.sub.sp);
    idy = unique(queryHits(findOverlaps(x.sub.sp@feature, features.gr)));
    x.sub.sp = x.sub.sp[,idy,mat="bmat"];
    x.sub.sp = runJaccard2(x.sub.sp, x.anchor.sp);
    x.sub.sp = runNormJaccard(
        obj=x.sub.sp,
        tmp.folder=tempdir(),
        method="normOVE",
        row.center=TRUE,
        row.scale=TRUE,
        low.threshold=-5,
        high.threshold=5,
        do.par=TRUE,
        ncell.chunk=1000,
        num.cores=10,
        seed.use=10
    );

    nmat.out <- x.sub.sp@jmat@jmat
    
   
    nmat.ls[[fname]] = x.sub.sp@jmat@jmat;
    rm(x.sub.sp);
    gc();
}




# make parallel
message("get nmat in parallel")
get_nmat <- function(fname, x.anchor.sp, x.sp, nmat.ls){
    print(fname);
    x.sub.sp = createSnap(fname, sample=as.character(fname));
    x.sub.sp = x.sub.sp[which(x.sub.sp@barcode %in% x.sp@barcode[x.sp@file == fname]),];
    x.sub.sp = addBmatToSnap(x.sub.sp, bin.size=5000);
    x.sub.sp = makeBinary(x.sub.sp);
    idy = unique(queryHits(findOverlaps(x.sub.sp@feature, features.gr)));
    x.sub.sp = x.sub.sp[,idy,mat="bmat"];
    x.sub.sp = runJaccard2(x.sub.sp, x.anchor.sp);
    x.sub.sp = runNormJaccard(
        obj=x.sub.sp,
        tmp.folder=tempdir(),
        method="normOVE",
        row.center=TRUE,
        row.scale=TRUE,
        low.threshold=-5,
        high.threshold=5,
        do.par=TRUE,
        ncell.chunk=1000,
        num.cores=10,
        seed.use=10
    );
    nmat.out <- x.sub.sp@jmat@jmat
    pname <- unlist(strsplit(tail(unlist(strsplit(fname, "/")),1),".snap"))
    nmat.name <- paste(outF, pname, "nmat", sep=".")
    fwrite(nmat.out, nmat.name, sep="\t", col.names=F, row.names=F, quote=F)
    nmat.ls[[fname]] = x.sub.sp@jmat@jmat;
    rm(x.sub.sp);
    gc();
}

mclapply(fnames, get_nmat, x.anchor.sp = x.anchor.sp, x.sp = x.sp, nmat.ls = nmat.ls, mc.cores = cpus)




nmat = nmat.ls[[1]];
if(length(fnames) > 1){
    for(i in seq(2, length(fnames))){
        nmat = rbind(nmat, nmat.ls[[i]]);
    }    
}
rm(nmat.ls);
gc();

x.sp@jmat@jmat = nmat;

## PCA / svd
tic("runDimReduct")
x.sp = runDimReduct(
    x.sp,
    pc.num=pc_num,
    input.mat="jmat",
    method="svd",
    center=TRUE,
    scale=FALSE,
    seed.use=10
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
    labs.subtitle=NULL,
    );

plotDimReductPW(
    obj=x.sp,
    pca.dims=1:pc_num,
    );
dev.off()

outfname = paste(outF, ".pre.RData",sep="")
save(x.sp, file=outfname)

outmeta = paste(outF, ".pre.meta.txt",sep="")
write.table(x.sp@metaData, file=outmeta, sep="\t", quote=F, col.names=T, row.names=F)


x.sp = runKNN(
    obj=x.sp,
    pca.dims=2:30,
    weight.by.sd=TRUE,
    k=50
);
    
x.sp = runCluster(
    obj=x.sp,
    tmp.folder=tempdir(),
    louvain.lib="R-igraph",
    seed.use=10,
    resolution=1
);
    
x.sp = runViz(
    obj=x.sp, 
    tmp.folder=tempdir(),
    dims=2,
    pca.dims=2:30, 
    weight.by.sd=TRUE,
    method="Rtsne",
    fast_tsne_path=NULL,
    Y.init=NULL,
    seed.use=10,
    num.cores=5
);

x.sp = runViz(
    obj=x.sp, 
    tmp.folder=tempdir(),
    dims=2,
    pca.dims=2:30, 
    weight.by.sd=TRUE,
    method="umap",
    fast_tsne_path=NULL,
    Y.init=NULL,
    seed.use=10,
    num.cores=5
);

pdf(paste(outF, ".cluster.pdf",sep=""))
plotViz(
    obj=x.sp, 
    method="tsne", 
    point.size=0.5, 
    point.shape=19, 
    point.alpha=0.8, 
    point.color="cluster", 
    text.add=FALSE,
    text.size=1.5,
    text.color="black",
    text.halo.add=TRUE,
    text.halo.color="white",
    text.halo.width=0.2,
    down.sample=10000,
    pdf.file.name=NULL,
    pdf.width=7, 
    pdf.height=7,
    legend.add=TRUE,
    main="SnapATAC"
);

plotViz(
    obj=x.sp, 
    method="umap", 
    point.size=0.5, 
    point.shape=19, 
    point.alpha=0.8, 
    point.color="cluster", 
    text.add=FALSE,
    text.size=1.5,
    text.color="black",
    text.halo.add=TRUE,
    text.halo.color="white",
    text.halo.width=0.2,
    down.sample=10000,
    pdf.file.name=NULL,
    pdf.width=7, 
    pdf.height=7,
    legend.add=TRUE,
    main="SnapATAC"
);

dev.off()

outmetaf <- paste(outF, ".cluster.meta.txt", sep="")
outmetamx <- cbind(x.sp@sample, x.sp@metaData, x.sp@cluster, x.sp@tsne, x.sp@umap)
write.table(outmetamx, outmetaf, row.names=F, col.names=T, sep="\t", quote=F)

outfname = paste(outF, ".cluster.RData",sep="")
save(x.sp, file=outfname)

