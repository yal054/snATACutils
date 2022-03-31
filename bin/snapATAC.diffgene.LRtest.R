#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input", required=TRUE, help="load clustered RData")
parser$add_argument("-g", "--gset", default="/projects/ps-renlab/yangli/genome/mm10/gencode.vM16.gene.bed", help="gene annotation in bed format")
parser$add_argument("-a", "--attr", required=TRUE, help="attribute's column name used for regression: MajorRegion | SubType | MajorType")
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
gset = args$gset
attr = args$attr
model = args$model
cpus = as.numeric(args$thread)
outF = args$output

#---------
# init
x.sp <- get(load(RDataF))
 
#-----------------------
# construct cds obj
cell_meta <- x.sp@metaData
rownames(cell_meta) <- cell_meta$cellID

ganno <- fread(gset, sep="\t", header=F)
ganno <- ganno[,c(4,7,8)]
colnames(ganno) <- c("gencode", "id", "type")
ganno <- ganno[!duplicated(ganno$id), ]

gene_anno <- data.frame(id=colnames(x.sp@gmat))
gene_anno$gene_short_name <- gene_anno$id
gene_anno <- join(gene_anno, ganno, by="id")
rownames(gene_anno) <- gene_anno$id
gene_anno <- as.data.frame(gene_anno)

gene_mat <- Matrix::t(x.sp@gmat)
colnames(gene_mat) <- cell_meta$cellID
rownames(gene_mat) <- as.character(gene_anno$id)

if(model == "binomial"){
  gene_mat@x[gene_mat@x > 1] <- 1 # binarize peak mat
}

cds <- new_cell_data_set(gene_mat,
                         cell_metadata = cell_meta,
                         gene_metadata = gene_anno)

cds <- detect_genes(cds)

save(cds, file=paste(outF, "cdsObj.RData", sep="."))
rm(x.sp)
rm(gene_mat)
rm(cell_meta)
rm(gene_anno)

#-----------------------
# diff gene
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

lr_models.df$id <- rowData(cds)$gencode
lr_models.df$type <- rowData(cds)$type

fwrite(lr_models.df, file=paste(outF, "LRtest.out.tsv", sep="."), quote =F, col.names = T, row.names = F, sep="\t")

