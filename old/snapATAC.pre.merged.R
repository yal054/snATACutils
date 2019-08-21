#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input", required=TRUE, help="input txt with the information of path of snap file and name")
parser$add_argument("--fragment_num", default = 1000, help="cutoff of fragment.num [default %(default)s]")
parser$add_argument("--mito_ratio", default = 0.1, help="cutoff of mito.ratio [default %(default)s]")
parser$add_argument("--dup_ratio", default = 0.5, help="cutoff of dup.ratio [default %(default)s]")
parser$add_argument("--umap_ratio", default = 0.8, help="cutoff of umap.ratio [default %(default)s]")
parser$add_argument("--pair_ratio", default = 0.9, help="cutoff of pair.ratio [default %(default)s]")
parser$add_argument("--bin_size", default = 5000, help="binSize to use [default %(default)s]")
parser$add_argument("--black_list", required=TRUE, help="black list file")
parser$add_argument("--pc_num", default = 50, help="num of PCA components [default %(default)s]")
parser$add_argument("--cpu", default = 5, help="# of cpus [default %(default)s]")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

suppressPackageStartupMessages(library("SnapATAC"))
library("RColorBrewer")

inputF = args$input
fragment_num = as.numeric(args$fragment_num)
mito_ratio = as.numeric(args$mito_ratio)
dup_ratio = as.numeric(args$dup_ratio)
umap_ratio = as.numeric(args$umap_ratio)
pair_ratio = as.numeric(args$pair_ratio)
bin_size = as.numeric(args$bin_size)
black_list = args$black_list
pc_num = as.numeric(args$pc_num)
cpus = as.numeric(args$cpu)
outF = args$output

suppressPackageStartupMessages(library("SnapATAC"))

inputf <- read.table(inputF,sep="\t",header=F) 
fnames <- as.character(inputf[,2])
pnames <- as.character(inputf[,1])

#x.sp.ls <- lapply(as.list(1:length(fnames)), function(i){
#    prefix_i = pnames[i]
#    x1.sp = createSnap(fnames[i], metaData=TRUE);
#    x1.sp = filterCells(x1.sp, 
#                       subset.names=c("fragment.num",
#                                      "mito.ratio",
#                                      "dup.ratio",
#                                      "umap.ratio",
#                                      "pair.ratio"),
#                       low.thresholds=c(fragment_num, 0, 0, umap_ratio, pair_ratio),
#                       high.thresholds=c(Inf, mito_ratio, dup_ratio, 1, 1)
#                      );
#    x1.sp = addBmat(x1.sp, fnames[i], binSize=bin_size);
#    x1.sp = addGmat(x1.sp, fnames[i]);
#    x1.sp@barcode = paste(pnames[i], x1.sp@barcode, sep=".");
#    return(x1.sp)
#})

x.sp.ls <- lapply(as.list(1:length(fnames)), function(i){
    prefix_i = pnames[i]
    x1.sp = createSnap(fnames[i], metaData=TRUE);
    x1.sp = addBmat(x1.sp, fnames[i], binSize=bin_size);
    x1.sp = addGmat(x1.sp, fnames[i]);
    x1.sp@barcode = paste(pnames[i], x1.sp@barcode, sep=".");
    return(x1.sp)
})


names(x.sp.ls) <- pnames;

rBind <- function(obj1, ...){
  UseMethod("rBind", obj1);
}

rBind.default <- function(obj1, obj2){
    # only the following slots can be combined
    # barcode, feature, metaData, cmat, bmat
    # among these slots, barcode, feature, cmat are enforced, the others are optional
    # the rest slots must be set to be empty
    
    if(!is.snap(obj1)){stop(paste("Error @rBind: obj1 is not a snap object!", sep=""))};
    if(!is.snap(obj2)){stop(paste("Error @rBind: obj2 is not a snap object!", sep=""))};

    # barcode from obj1 and obj2
    barcode1 = obj1@barcode;
    barcode2 = obj2@barcode;    
    
    # check barcode name, if there exists duplicate barcode raise error and exist
    if(length(unique(c(barcode1, barcode2))) < length(barcode1) + length(barcode2)){
        stop("Error: @rBind: identifcal barcodes found in obj1 and obj2!")
    }
    barcode = c(barcode1, barcode2);
    
    # check meta data
    if(nrow(obj1@metaData) > 0 && nrow(obj2@metaData) > 0){
        metaData = rbind(obj1@metaData, obj2@metaData);        
    }else{
        metaData = data.frame();
    }
    
    # check feature
    feature1 = obj1@feature;
    feature2 = obj2@feature;
    if((length(feature1) == 0) != (length(feature2) == 0)){
        stop("different feature found in obj1 and obj2!")
    }else{
        if(length(feature1) > 0){
            if(FALSE %in% (feature1$name == feature2$name)){
                stop("Error: @rBind: different feature found in obj1 and obj2!")
            }
            feature = feature1;                    
        }else{
            feature = feature1;                                
        }
    }
    
    # check peak
    peak1 = obj1@peak;
    peak2 = obj2@peak;
    if((length(peak1) == 0) != (length(peak2) == 0)){
        stop("different peak found in obj1 and obj2!")
    }else{
        if(length(peak1) > 0){
            if(FALSE %in% (peak1$name == peak2$name)){
                stop("Error: @rBind: different feature found in obj1 and obj2!")
            }
            peak = peak1;                    
        }else{
            peak = peak1;                                
        }
    }
    
    # check bmat    
    bmat1 = obj1@bmat;
    bmat2 = obj2@bmat;
    if((length(bmat1) == 0) != (length(bmat2) == 0)){
        stop("bmat has different dimentions in obj1 and obj2!")
    }else{
        bmat = Matrix::rBind(bmat1, bmat2);
    }

    # check gmat    
    gmat1 = obj1@gmat;
    gmat2 = obj2@gmat;
    if((length(gmat1) == 0) != (length(gmat2) == 0)){
        stop("gmat has different dimentions in obj1 and obj2!")
    }else{
        gmat = Matrix::rBind(gmat1, gmat2);
    }

    # check pmat    
    pmat1 = obj1@pmat;
    pmat2 = obj2@pmat;
    if((length(pmat1) == 0) != (length(pmat2) == 0)){
        stop("pmat has different dimentions in obj1 and obj2!")
    }else{
        pmat = Matrix::rBind(pmat1, pmat2);
    }

    res = newSnap();
    res@barcode = barcode;
    res@metaData = metaData;
    res@bmat = bmat;
    res@pmat = pmat;
    res@feature = feature;
    res@peak = peak;
    res@gmat = gmat;
    return(res)
}

# combine all the slices as one snap object
x.sp = Reduce(rBind, x.sp.ls);

#x.sp = filterCells(x.sp,
#                       subset.names=c("fragment.num",
#                                      "mito.ratio",
#                                      "dup.ratio",
#                                      "umap.ratio",
#                                      "pair.ratio"),
#                       low.thresholds=c(fragment_num, 0, 0, umap_ratio, pair_ratio),
#                       high.thresholds=c(Inf, mito_ratio, dup_ratio, 1, 1)
#                      );

x.sp = filterCells(x.sp,
                       subset.names="fragment.num",
                       low.thresholds=1000,
                       high.thresholds=Inf
                      );

pdf(paste(outF, ".pre.sta.pdf",sep=""))
plotBarcode(x.sp);
dev.off()

## Matrix Binarization
x.sp = makeBinary(x.sp, mat="bmat");
checkBinSize(x.sp);

## Feature Selection (filter bins according to black list and chrom)
black_list = read.table(black_list);
black_list.gr = GRanges(
                          black_list[,1],
                          IRanges(black_list[,2], black_list[,3])
                         );
idy1 = queryHits(findOverlaps(x.sp@feature, black_list.gr));
idy2 = grep("chrM|random", x.sp@feature);
idy = unique(c(idy1, idy2));
x.sp = x.sp[,-idy, mat="bmat"];

plotBinCoverage <- function(obj,...) {
  UseMethod("plotBinCoverage", obj);
}

#' @export
plotBinCoverage.default <- function(obj, ...){    
    # check the input
    if(!(class(obj)=="snap")){stop(paste("Error @plotBarcode: obj is not a snap object!", sep=""))};
    cov = Matrix::colSums(obj@bmat);
    idy = seq(ncol(colSums(obj@bmat)));    
    cov = cov[which(cov > 0)];    
    cov = log10(cov + 1);
    cov = (cov - mean(cov)) / sd(cov);
    hist(cov, ...);
}

pdf(paste(outF, ".pre.BinCoverage.pdf",sep=""))
plotBinCoverage(x.sp);
dev.off()

x.sp = filterBins(x.sp, low.threshold=-2, high.threshold=2, mat="bmat");

## Jaccard Matrix & Normlaization

idx <- which(rowSums(x.sp@bmat)!=0)
x.sp@barcode <- x.sp@barcode[idx]
x.sp@bmat <- x.sp@bmat[idx,]
x.sp@gmat <- x.sp@gmat[idx,]

x.sp = calJaccard(
    x.sp,
    mat = "bmat",
    ncell.chunk=1000,
    max.var=2000,
    seed.use=10,
    norm.method="normOVE",
    row.center=TRUE,
    row.scale=TRUE,
    low.threshold=-5,
    high.threshold=5,
    keep.jmat=TRUE,
    do.par=TRUE,
    num.cores=cpus)

## PCA
x.sp = runPCA(
    x.sp,
    pc.num=pc_num,
    input.mat="nmat",
    method="svd",
    weight.by.sd=TRUE,
    center=TRUE,
    scale=FALSE,
    seed.use=10
    )

pdf(paste(outF, ".pre.plotPCA.pdf",sep=""))
plotPCA(x.sp, method="elbow");
plotPCA(x.sp, method="pairwise");
dev.off()

outfname = paste(outF, ".pre.RData",sep="")
save(x.sp, file=outfname)

