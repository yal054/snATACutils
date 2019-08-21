#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input", required=TRUE, help="load preprocessed RData")
parser$add_argument("--path_to_homer", default = "/home/yangli1/apps/homer/bin/findMotifsGenome.pl", help="path to homer [default %(default)s]")

parser$add_argument("--cpu", default = 4, help="number of cpus [default %(default)s]")
parser$add_argument("--pval", default = 0.01, help="p-value [default %(default)s]")
parser$add_argument("--fdr", default = 0.05, help="fdr [default %(default)s]")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

suppressPackageStartupMessages(library("SnapATAC"))
library("RColorBrewer")
library("edgeR")

RDataF = args$input
path_to_homer = args$path_to_homer
cpus = args@cpu


metaF = args$meta
bcv = as.numeric(args$bcv)
pval = as.numeric(args$pval)
fdr = as.numeric(args$fdr)
outF = args$output

load(RDataF)

#clsN <- unique(x.sp@cluster)
clsN <- which(table(x.sp@cluster)>300)

for(cls in clsN){

motifs = runHomer(
    x.sp[,idy_sst,"pmat"], 
    mat = "pmat",
    path.to.homer = "/projects/ps-renlab/r3fang/public_html/softwares/homer/bin/findMotifsGenome.pl",
    result.dir = "./homer/Sst",
    num.cores=5,
    genome = 'mm10',
    motif.length = 10,
    scan.size = 300,
    optimize.count = 2,
    background = 'automatic',
    local.background = FALSE,
    only.known = FALSE,
    only.denovo = FALSE,
    fdr.num = 5,
    cache = 100,
    overwrite = TRUE,
    keep.minimal = FALSE
    );













idy = findDAR(
    obj=x.sp,
    mat="pmat",
    cluster.pos=cls,
    cluster.neg=NULL,
    bcv=bcv,
    fdr=fdr,
    pvalue=pval,
    test.method="exactTest",
    seed.use=10
    );

if(length(idy) > 0){

y = SnapATAC::rowSums(x.sp[,idy,mat="pmat"], mat="pmat") / SnapATAC::rowSums(x.sp, mat="pmat") * 1000000;
# normalize to zscore
y_z = (y - mean(y)) / sd(y);

outZscore <- as.data.frame(cbind(x.sp@sample, x.sp@barcode, x.sp@cluster, y, y_z))
outZscoref = paste(outF, ".", cls, "_diffpeaks.zscore.txt", sep="")
write.table(outZscore, outZscoref, sep="\t",col.names=T, row.names=F,quote=F)

outpdf = paste(outF, ".", cls, "_diffpeaks.boxplot.pdf", sep="")
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

peaks <- as.data.frame(x.sp@peak);
diffpeaks <- peaks[idy,]
diffpeaksf = paste(outF, ".", cls, "_diffpeaks.txt", sep="")
write.table(diffpeaks, diffpeaksf, sep="\t",col.names=T, row.names=F,quote=F)
}
}


