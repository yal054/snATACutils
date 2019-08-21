#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input", required=TRUE, help="input txt with the information of path of snap file and name")
parser$add_argument("--RData", required=TRUE, help="input RData with bmat & gmat")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

inputF = args$input
RDataF = args$RData
outF = args$output

suppressPackageStartupMessages(library("SnapATAC"))

inputf <- read.table(inputF,sep="\t",header=F) 
fnames <- as.character(inputf[,2])
pnames <- as.character(inputf[,1])

x.sp <- get(load(RDataF))

x.sp.ls <- lapply(as.list(1:length(fnames)), function(i){
    prefix_i = pnames[i]
    x1.sp = createSnap(fnames[i], metaData=TRUE);
    x1.sp = addPmatToSnap(x1.sp, fnames[i]);
    x1.sp@barcode = paste(pnames[i], x1.sp@barcode, sep=".");
    return(x1.sp)
})

names(x.sp.ls) <- pnames;

# combine all the slices as one snap object
x.sp.pmat = Reduce(snapRbind, x.sp.ls);
x.sp.pmat = makeBinary(x.sp.pmat, mat="pmat");

matchpmat <- x.sp.pmat@pmat[x.sp.pmat@barcode %in% x.sp@barcode,]
x.sp@pmat <- matchpmat
x.sp@peak <- x.sp.pmat@peak

outfname = paste(outF, ".full.RData",sep="")
save(x.sp, file=outfname)







