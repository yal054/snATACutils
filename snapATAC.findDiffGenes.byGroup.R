#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input", required=TRUE, help="load gmat RData")
parser$add_argument("-n", "--num", default = 300, help="number of cells within one cluster [default %(default)s]")
parser$add_argument("--clustLabel", required=TRUE, help="cluster labels used for subset, 2 column: label, cluster #")
parser$add_argument("--bcv", default = 0.1, help="number of bcv [default %(default)s]")
parser$add_argument("--pval", default = 0.01, help="p-value [default %(default)s]")
parser$add_argument("--fdr", default = 0.05, help="fdr [default %(default)s]")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

suppressPackageStartupMessages(library("SnapATAC"))
library("tictoc")

RDataF = args$input
N = as.numeric(args$num)
clustLabelsF = args$clustLabel
bcv = as.numeric(args$bcv)
pval = as.numeric(args$pval)
fdr = as.numeric(args$fdr)
outF = args$output

load(RDataF)
clustLabels.df <- read.table(clustLabelsF)

# map cluster to group id
library(plyr);
current.cluster.ids <- clustLabels.df[,2]
new.cluster.ids <- as.character(clustLabels.df[,1])

x.sp@cluster <- plyr::mapvalues(
    x = x.sp@cluster, 
    from = current.cluster.ids, 
    to = new.cluster.ids
    );

clsN <- names(which(table(x.sp@cluster)>N))

for(cls in clsN){
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


res$FDR <- p.adjust(res$PValue)
outgenes <- as.data.frame(cbind(symbol=rownames(res), res, cls))
outgenesf = paste(outF, ".", cls, ".diffgenes.sta.txt", sep="")
write.table(outgenes, outgenesf, sep="\t",col.names=F, row.names=F,quote=F)


idy <- which(res$PValue<=pval & res$FDR <= fdr)

if(length(idy) > 0){

y = SnapATAC::rowSums(x.sp[,idy,mat="gmat"], mat="gmat") / SnapATAC::rowSums(x.sp, mat="gmat") * 1000000;
# normalize to zscore
y_z = (y - mean(y)) / sd(y);

outZscore <- as.data.frame(cbind(x.sp@sample, x.sp@barcode, x.sp@cluster, y, y_z))
outZscoref = paste(outF, ".", cls, ".diffgenes.zscore.txt", sep="")
write.table(outZscore, outZscoref, sep="\t",col.names=T, row.names=F,quote=F)

outpdf = paste(outF, ".", cls, ".diffgenes.boxplot.pdf", sep="")
pdf(outpdf)
boxPlotFeature(
    obj = x.sp,
    feature = y_z,
    outline = FALSE,
    ylab = "zscore of RPM",
    main = "DARs Enrichment",
    add.point = TRUE,
    point.size = 0.2,
    point.shape = 19,
    point.alpha = 0.5,
    pdf.file.name=NULL,
    pdf.height=7,
    pdf.width=7
    );
dev.off()

genes <- colnames(x.sp@gmat)
diffgenes <- genes[idy]
diffgenesf = paste(outF, ".", cls, ".diffgenes.txt", sep="")
write.table(diffgenes, diffgenesf, sep="\t",col.names=T, row.names=F,quote=F)
}
}


