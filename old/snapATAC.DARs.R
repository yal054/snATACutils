#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input", required=TRUE, help="load preprocessed RData")
parser$add_argument("-d", "--pc_dim", default = 25, help="num of PCA dims used for clustering [default %(default)s]")
parser$add_argument("--bcv", default = 0.4, help="number of bcv [default %(default)s]")
parser$add_argument("--pval", default = 0.05, help="p-value [default %(default)s]")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

suppressPackageStartupMessages(library("SnapATAC"))
library("RColorBrewer")
library("umap")

RDataF = args$input
dims = as.numeric(args$pc_dim)
bcv = as.numeric(args$bcv)
pval = as.numeric(args$pval)
outF = args$output

load(RDataF)

clsN <- max(as.numeric(x.sp@cluster))
for(i in 1:clsN){

cls_barcodes = x.sp@barcode[which(x.sp@cluster == i)];

idx <- which(rowSums(x.sp@pmat)!=0)
x.sp@barcode <- x.sp@barcode[idx]
x.sp@bmat <- x.sp@bmat[idx,]
x.sp@gmat <- x.sp@gmat[idx,]
x.sp@pmat <- x.sp@pmat[idx,]
x.sp@smat <- x.sp@smat[idx,]
x.sp@jmat <- x.sp@jmat[idx,]
x.sp@nmat <- x.sp@nmat[idx,]
x.sp@tsne <- x.sp@tsne[idx,]
x.sp@umap <- x.sp@umap[idx,]
x.sp@cluster <- x.sp@cluster[idx]

#idy = findDAR(
#    object=x.sp, 
#    mat="pmat",
#    barcodes.sel=cls_barcodes,
#    bcv=bcv, 
#    rand.seed=10, 
#    pca_dims=1:dims,
#    fdr=1e-1,
#    method="exactTest",
#    k=30
#    )

cmat = x.sp@pmat;
pca_dims = seq(ncol(x.sp@smat));
set.seed(10);

message("Identifying accessible regions using postive sample");
idx.pos = which(x.sp@barcode %in% cls_barcodes);
idx.neg = setdiff(seq(nrow(x.sp)), idx.pos);
nn.num = min(length(idx.pos), length(idx.neg));
idx.neg = sample(idx.neg, nn.num);

cmat.pos = cmat[idx.pos,];
cmat.neg = cmat[idx.neg,];

x = data.frame(Matrix::colSums(cmat.neg), Matrix::colSums(cmat.pos))
group <- factor(c(1,2));
design <- model.matrix(~group);
y <- DGEList(counts=x, group=group);
tb.pos <- exactTest(y, dispersion=bcv^2)$table;

message("Identifying accessible regions using negative sample")
neg.idx.pos = sample(seq(nrow(x.sp)), length(idx.pos))
neg.idx.neg = setdiff(seq(nrow(x.sp)), neg.idx.pos);
nn.num = min(length(neg.idx.pos), length(neg.idx.neg));
neg.idx.neg = sample(neg.idx.neg, nn.num);

neg.cmat.pos = cmat[neg.idx.pos,];
neg.cmat.neg = cmat[neg.idx.neg,];
x = data.frame(Matrix::colSums(neg.cmat.pos), Matrix::colSums(neg.cmat.neg))
group <- factor(c(1,2));
design <- model.matrix(~group);
y <- DGEList(counts=x, group=group);

tb.neg <- exactTest(y, dispersion=bcv^2)$table;

message("calculating p-value and FDR table")
fdr_table = data.frame()
for(p_i in seq(0.001, 0.05, by=0.001)){
    fdr_table = rbind(fdr_table,
        data.frame(
        PValue=p_i,
        fdr=length(which(tb.neg$PValue < p_i & tb.neg$logFC > 0)) / length(which(tb.pos$PValue < p_i & tb.pos$logFC > 0))
        )
    )
}

#if(length(which(fdr_table[,2] <= fdr)) > 0){
#        fdr.idx = max(which(fdr_table[,2] <= fdr));
#        idy <- which(tb.pos$PValue < fdr_table[fdr.idx,1] & tb.pos$logFC > 0)
#} 

idy <- which(tb.pos$logFC>=0 & tb.pos$PValue<=pval)

y = rowSums(x.sp@pmat[,idy]) / rowSums(x.sp@pmat);
y = y/max(y);

pdf(paste(outF, "cls.", i, ".seltpeaks.pdf", sep=""))
plot(x.sp@tsne, 
    col=alpha("blue", y), 
    yaxt='n', 
    xaxt="n",
    xlab="", 
    ylab="",
    cex=0.2
    );
dev.off()

peaks <- as.data.frame(x.sp@peak)[,c(1:3)];
outpeaks <- cbind(peaks, tb.pos)
seltpeaks <- outpeaks[idy,]

outpeaksf = paste(outF, "cls.", i, ".diffpeaks.txt", sep="")
write.table(outpeaks, outpeaksf, sep="\t",col.names=T, row.names=F,quote=F)

seltpeaksf = paste(outF, "cls.", i, ".seltpeaks.txt", sep="")
write.table(seltpeaks, seltpeaksf, sep="\t",col.names=T, row.names=F,quote=F)

}

