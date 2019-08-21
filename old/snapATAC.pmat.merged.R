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

suppressPackageStartupMessages(library("SnapATAC"))
library("RColorBrewer")

inputF = args$input
RDataF = args$RData
outF = args$output

inputf <- read.table(inputF,sep="\t",header=F) 
fnames <- as.character(inputf[,2])
pnames <- as.character(inputf[,1])

load(RDataF)

x.sp.ls <- lapply(as.list(1:length(fnames)), function(i){
    prefix_i = pnames[i]
    x1.sp = createSnap(fnames[i], metaData=TRUE);
    x1.sp = addPmat(x1.sp, fnames[i]);
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
x.sp.pmat = Reduce(rBind, x.sp.ls);
x.sp.pmat = makeBinary(x.sp.pmat, mat="pmat");

matchpmat <- x.sp.pmat@pmat[x.sp.pmat@barcode %in% x.sp@barcode,]
x.sp@pmat <- matchpmat
x.sp@peak <- x.sp.pmat@peak

outfname = paste(outF, ".full.RData",sep="")
save(x.sp, file=outfname)

