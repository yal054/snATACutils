#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input", required=TRUE, help="load preprocessed RData")
parser$add_argument("-m", "--meta", required=TRUE, help="load metatable")
parser$add_argument("-n", "--num", default = 50, help="number of cells within one cluster [default %(default)s]")
parser$add_argument("--bcv", default = 0.1, help="number of bcv [default %(default)s]")
parser$add_argument("--pval", default = 0.05, help="p-value [default %(default)s]")
parser$add_argument("--fdr", default = 0.05, help="fdr [default %(default)s]")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()
suppressPackageStartupMessages(library("SnapATAC"))
suppressPackageStartupMessages(library("data.table"))

RDataF = args$input
metaF = args$meta
N = args$num
bcv = as.numeric(args$bcv)
pval = as.numeric(args$pval)
fdr = as.numeric(args$fdr)
outF = args$output


# export marker genes
# load snap
load(RDataF)

clsN <- which(table(x.sp@cluster)>N)
#clsN <- unique(x.sp@cluster)

diffgenes.df <- data.frame(matrix(ncol = 7, nrow = 0))
x <- c("cluster", "gene", "PValue", "logFC", "logCPM", "label", "avg")
colnames(diffgenes.df) <- x

for(cls in clsN){
idx <- which(x.sp@cluster == cls)
avgExpr <- Matrix::colMeans(x.sp@gmat[idx,])

res = findDAR(
    obj=x.sp,
    mat="gmat",
    cluster.pos=cls,
    cluster.neg=NULL,
    bcv=bcv,
    fdr=fdr,
    pvalue=pval,
    test.method="exactTest",
    seed.use=10
    );

outgenes <- data.frame(cluster=cls, gene=colnames(x.sp@gmat), res[,c(3,1,2,4)], avg=avgExpr)
#outgenesf = paste(outF, ".", cls, ".diffgenes.sta.txt", sep="")
#write.table(outgenes, outgenesf, sep="\t",col.names=T, row.names=F,quote=F)

idy <- which(res$PValue<=pval & res$label == 1)

if(length(idy) > 0){
diffgenes <- outgenes[idy,]
#diffgenesf = paste(outF, ".", cls, ".diffgenes.txt", sep="")
#write.table(diffgenes, diffgenesf, sep="\t",col.names=T, row.names=F,quote=F)
diffgenes.df <- rbind(diffgenes.df, diffgenes)
}
}

diffgeneFname <- paste(outF, "diffgenes.tsv", sep=".")
fwrite(diffgenes.df, diffgeneFname, sep="\t", col.names = T, row.names = F, quote = F)


# normalize 
x.sp = scaleCountMatrix(
    x.sp,
    cov=x.sp@metaData$UQ,
    mat="gmat",
    method = "RPM"
    );

# export gmat to mtx format
gmat <- x.sp@gmat
rownames(gmat) <- paste(x.sp@sample, x.sp@barcode, sep=".")
gmat <- t(gmat)

# for every genes, scale to range 0-1
#tic("L2norm")
#gmat.rownorm <- Matrix::Diagonal(x = 1 / sqrt(Matrix::rowSums(gmat^2))) %*% gmat
#toc()
#tic("01norm")
#suppressPackageStartupMessages(library("clusterSim"))
#gmat.01norm <- data.Normalization(gmat, type="n4", normalization="row")
#toc()

suppressPackageStartupMessages(library("qlcMatrix"))
suppressPackageStartupMessages(library("tictoc"))
#tic("01norm")
#rowMin.vec <- rowMin(gmat, ignore.zero=F)
#gmat@x <- as.numeric(gmat@x - rowMin.vec[gmat@i+1])
#gmat <- Matrix::Diagonal(x = 1 / (rowMax(gmat, ignore.zero=F) - rowMin(gmat, ignore.zero=F)) ) %*% gmat
#toc()

#tic("01norm")
#rowMin.vec <- rowMin(gmat, ignore.zero=F)
#gmat.subtract <- sweep(gmat, 1, rowMin.vec, "-")
#gmat.01norm <- Matrix::Diagonal(x = 1 / (rowMax(gmat.subtract, ignore.zero=F) - rowMin(gmat.subtract, ignore.zero=F)) ) %*% gmat.subtract
#toc()

tic("01norm")
gmat <- Matrix::Diagonal(x = 1 / (rowMax(gmat, ignore.zero=F) - rowMin(gmat, ignore.zero=F)) ) %*% gmat
toc()


#rowMin.ls <- replicate(dim(foo)[2], rowMin(foo, ignore.zero=F))
#rowMin.ls <- lapply(rowMin.ls, as, "sparseMatrix")
#rowMin.sp <- do.call(cBind, rowMin.ls)
#foo - rowMin.sp

#gmat <- t(apply(gmat, MARGIN = 1, FUN = function(X) (X - min(X))/Matrix::diff(range(X))))
#gmat <- Matrix(gmat, sparse = TRUE) 

geneFname <- paste(outF, "gene.tsv", sep=".")
fwrite(as.data.frame(rownames(gmat)), geneFname, sep="\t", col.names = F, row.names = F, quote = F)
cellFname <- paste(outF, "cell.tsv", sep=".")
fwrite(as.data.frame(colnames(gmat)), cellFname, sep="\t", col.names = F, row.names = F, quote = F)
mtxFname <- paste(outF, "gmat.mtx", sep=".")
writeMM(gmat, file=mtxFname)

# export coords
barcodes.ls <- paste(x.sp@sample, x.sp@barcode, sep=".")
tsne.coords <- data.frame(cellId=barcodes.ls, x=x.sp@tsne[,1], y=x.sp@tsne[,2])
umap.coords <- data.frame(cellId=barcodes.ls, x=x.sp@umap[,1], y=x.sp@umap[,2])

tsneFname <- paste(outF, "tsne.coords.tsv", sep=".")
fwrite(tsne.coords, tsneFname, sep="\t", col.names = T, row.names = F, quote = F)
umapFname <- paste(outF, "umap.coords.tsv", sep=".")
fwrite(umap.coords, umapFname, sep="\t", col.names = T, row.names = F, quote = F)

# export metatable
sumTable <- read.table(metaF, sep="\t", header = T)
colnames(sumTable)[which(colnames(sumTable)=="datasets")] <- "sample"
metaTable <- data.frame(sample=x.sp@sample, x.sp@metaData, cluster=x.sp@cluster)
metaTable <- merge(metaTable, sumTable, by = "sample")

.libPaths("/home/yangli1/R/x86_64-pc-linux-gnu-library/3.4") 
suppressPackageStartupMessages(library("tidyr"))
metaTable <- unite(metaTable, Cell, c("sample","barcode"), sep=".")

metaFname <- paste(outF, "meta.tsv", sep=".")
fwrite(metaTable, metaFname, sep="\t", col.names = T, row.names = F, quote = F)

