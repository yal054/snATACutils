#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input", required=TRUE, help="load preprocessed RData")
parser$add_argument("-n", "--num", default = 300, help="number of cells within one cluster [default %(default)s]")
parser$add_argument("--bcv", default = 0.1, help="number of bcv [default %(default)s]")
parser$add_argument("--pval", default = 0.01, help="p-value [default %(default)s]")
parser$add_argument("--fdr", default = 0.05, help="fdr [default %(default)s]")
parser$add_argument("--path_to_homer", default = "/home/yangli1/apps/homer/bin/findMotifsGenome.pl", help="path to homer [default %(default)s]")
parser$add_argument("--cpu", default = 4, help="number of cpus [default %(default)s]")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

suppressPackageStartupMessages(library("SnapATAC"))

RDataF = args$input
N = args$num
bcv = as.numeric(args$bcv)
pval = as.numeric(args$pval)
fdr = as.numeric(args$fdr)
path_to_homer = args$path_to_homer
cpus = args$cpu
outF = args$output

load(RDataF)

#clsN <- unique(x.sp@cluster)
clsN <- which(table(x.sp@cluster)>N)

for(cls in clsN){
print(cls)
res = findDAR(
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

outpeaks <- data.frame(cbind(x.sp@peak, res, cls))
outpeaksf = paste(outF, ".", cls, ".diffpeaks.sta.txt", sep="")
write.table(outpeaks, outpeaksf, sep="\t",col.names=F, row.names=F,quote=F)

idy <- which(res$PValue<=pval & res$label == 1)

if(length(idy) > 0){

y = SnapATAC::rowSums(x.sp[,idy,mat="pmat"], mat="pmat") / SnapATAC::rowSums(x.sp, mat="pmat") * 1000000;
# normalize to zscore
y_z = (y - mean(y)) / sd(y);

outZscore <- as.data.frame(cbind(x.sp@sample, x.sp@barcode, x.sp@cluster, y, y_z, cls))
outZscoref = paste(outF, ".", cls, ".diffpeaks.zscore.txt", sep="")
write.table(outZscore, outZscoref, sep="\t",col.names=T, row.names=F,quote=F)

outpdf = paste(outF, ".", cls, ".diffpeaks.boxplot.pdf", sep="")
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
diffpeaksf = paste(outF, ".", cls, ".diffpeaks.txt", sep="")
write.table(diffpeaks, diffpeaksf, sep="\t",col.names=T, row.names=F,quote=F)

motif_dir <- paste(outF, ".", cls, ".motif", sep="") 
#dir.create(motif_dir)

# find motif
motifs = runHomer(
    x.sp[,idy,"pmat"],
    mat = "pmat",
    path.to.homer = path_to_homer,
    result.dir = motif_dir,
    num.cores=cpus,
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

motifsf = paste(outF, ".", cls, ".motifs.txt", sep="")
write.table(motifs, motifsf, sep="\t",col.names=T, row.names=F,quote=F)
}
}


