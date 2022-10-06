#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input", required=TRUE, help="load clustered RData")
parser$add_argument("-p", "--pset", required=TRUE, help="column name for used peak set: subclass | celltype | MajorType | SubType")
parser$add_argument("--path_to_snap", default="/projects/ps-renlab/yangli/projects/HBA/03.peakcall/snap/", help="path to snap files")
parser$add_argument("--path_to_peak", default="/projects/ps-renlab/yangli/projects/CEMBA/01.joint_dat/rs1cemba/tracks/cCREs/", help="path to peak set for each cluster")
parser$add_argument("-a", "--attr", required=TRUE, help="attribute's column name used for regression: region | subclass | celltype | MajorRegion | SubRegion | DissectionRegion | SubType")
parser$add_argument("--min_num", default=100, help="min # of cell in attr [default %(default)s]")
parser$add_argument("--dn_num", required=FALSE, default=NULL, help="downsample # of cell in attr [default %(default)s]")
parser$add_argument("-m", "--model", default = "negbinomial", help="model used for fitting: negbinomial | binomial | quasipoisson [default %(default)s]")
parser$add_argument("-t", "--thread", default = 10, help="# of cpus [default %(default)s]")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

suppressPackageStartupMessages(library("SnapATAC"))
suppressPackageStartupMessages(library("GenomicRanges"));
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library("monocle3"))
suppressPackageStartupMessages(library("foreach"))
suppressPackageStartupMessages(library("doParallel"))
#suppressPackageStartupMessages(library("ggplot2"))
library("tictoc")

RDataF = args$input
pset = args$pset
path_to_peak = args$path_to_peak
path_to_snap = args$path_to_snap
attr = args$attr
min_num = as.numeric(args$min_num)
dn_num = as.numeric(args$dn_num)
model = args$model
cpus = as.numeric(args$thread)
outF = args$output

#---------
# init
x.sp <- get(load(RDataF))
attr.ncell <- table(x.sp@metaData[, attr])
sel.ncell <- names(attr.ncell[attr.ncell >= min_num])
idx <- which(x.sp@metaData[, attr] %in% sel.ncell)
if(length(idx) < nrow(x.sp)){
x.sp <- x.sp[idx, ]
}

x.sp@metaData$replicate <- x.sp@metaData$donor
 
#-----------------
# downsample
if( !is.null(dn_num)){

attr.list <- unique(x.sp@metaData[,attr])

idx.list <- c()
for(a in attr.list){
   idx <- which(x.sp@metaData[,attr] == a)
   if(length(idx) > dn_num){
       set.seed(2022)
       idx.sel <- sample(idx, dn_num)
   }else{
       idx.sel <- idx
   }
   idx.list <- c(idx.list, idx.sel)
}

idx.list <- sort(idx.list)
x.sp <- x.sp[idx.list, ]
}


#-----------
# load matrix
x.sp@file <- paste0(path_to_snap, x.sp@sample, ".snap")

sample.sel <- which(table(x.sp@sample)<=10)
if(length(sample.sel)>0){
idx <- which(x.sp@sample %in% names(sample.sel))
x.sp <- x.sp[-idx, ]
}

if( dim(x.sp@gmat)[1] == 0 ){
x.sp <- addPmatToSnap(x.sp)
}

save(x.sp, file=paste(outF, "snapObj.RData", sep="."))

#----------------
# filter used peaks
type.lst <- unique(x.sp@metaData[,pset])

if(length(type.lst)==1){
inf = list.files(path = path_to_peak, pattern = paste("*", "\\.", type.lst, ".bed",sep=""), full.names=T )
peak.lst <- fread(inf, sep="\t", header=F);
}else{
peak.lst <- list()
for(i in type.lst){
  inf = list.files(path = path_to_peak, pattern = paste("*", "\\.", i, ".bed",sep=""), full.names=T )
  inp = fread(inf, sep="\t", header=F);
  peak.lst[[i]] <- inp
}
peak.lst <- do.call(rbind, peak.lst)
peak.lst <- unique(peak.lst)
}

peak.gr = GRanges(
    peak.lst[[1]],
    IRanges(peak.lst[[2]], peak.lst[[3]])
  );
peak.gr$peak <- peak.lst[[4]]

#-----------------------
# construct cds obj
cell_meta <- x.sp@metaData
rownames(cell_meta) <- cell_meta$cellID

peak_anno <- data.frame(id=as.character(x.sp@peak))
rownames(peak_anno) <- peak_anno$id
peak_anno$gene_short_name <- paste("peak", 1:nrow(peak_anno), sep="_")

peak_mat <- t(x.sp@pmat)
colnames(peak_mat) <- cell_meta$cellID
rownames(peak_mat) <- peak_anno$peak
if(model == "binomial"){
peak_mat@x[peak_mat@x > 1] <- 1 # binarize peak mat
}

idx <- which(peak_anno$gene_short_name %in% peak.gr$peak)
idy = queryHits(
    findOverlaps(x.sp@peak, peak.gr)
  );
if(identical(idx, idy)){
  peak_mat <- peak_mat[idx, ]
  peak_anno <- peak_anno[idx, ]
}

cds <- new_cell_data_set(peak_mat,
                         cell_metadata = cell_meta,
                         gene_metadata = peak_anno)

cds <- detect_genes(cds)

save(cds, file=paste(outF, "cdsObj.RData", sep="."))
rm(x.sp)
rm(peak_mat)
rm(cell_meta)
rm(peak_anno)

#-----------------------
# diff peak
# u: total proportion of cells that are accessible at the ith site
# a: the membership of the jth cell in the cluster
# b: log10(total number of sites observed as accessible for the jth cell)
# e: the intercept

acc <- Matrix::colSums(exprs(cds))
cds[["log10acc"]] <- log10(acc)

# split to chunks
makeChunk <- function(x,n){
    split(x, cut(seq_along(x), n, labels = FALSE))
    }

chunks.lst <- makeChunk(1:dim(cds)[1], 100)
num.files = length(chunks.lst);

registerDoParallel(cpus)
lr_models.ls <- foreach (i=1:num.files, .inorder=TRUE) %dopar% {
print(i)
idx <- chunks.lst[[i]]
cds.used <- cds[idx,]
# full model
tic("full")
formula_str <- paste("~", attr, " + log10acc + replicate", sep="")
full_models <- fit_models(cds.used,
                                model_formula_str = formula_str,
                                expression_family=model, cores = 1)
toc()

# reduced model
tic("reduce")
reduce_models <- fit_models(cds.used,
                          model_formula_str = "~log10acc + replicate",
                          expression_family=model, cores = 1)
toc()
# likelihood ratio test
tic("compare")
lr_models <- compare_models(full_models, reduce_models)
toc()
rm(full_models)
rm(reduce_models)
return(lr_models)
}
closeAllConnections()

lr_models <- do.call(rbind, lr_models.ls)
lr_models.df <- as.data.frame(lr_models)
lr_models.df$q_value <- p.adjust(lr_models.df$p_value, method = "BH")

fwrite(lr_models.df, file=paste(outF, "LRtest.out.tsv", sep="."), quote =F, col.names = T, row.names = F, sep="\t")

