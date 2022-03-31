#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input", required=TRUE, help="load loom")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

suppressPackageStartupMessages(library("loomR"))
suppressPackageStartupMessages(library("Matrix"))
library("tictoc")

loomF = args$input
outF = args$output

# load loom
lfile <- connect(filename = loomF, mode = "r+")

# parse
tic("parse2mtx")
raw.matrix <- t(lfile[['matrix']][, ])
raw.matrix <- Matrix(data = raw.matrix, sparse = TRUE)

geneIDs <- lfile[['row_attrs/Gene']][]
cellIDs <- lfile[['col_attrs/CellID']][]

a <- lfile[['col_graphs/KNN/a']][] + 1
b <- lfile[['col_graphs/KNN/b']][] + 1
w <- lfile[['col_graphs/KNN/w']][]
knn <- sparseMatrix(i=a,j=b,x=w)

a <- lfile[['col_graphs/MKNN/a']][] + 1
b <- lfile[['col_graphs/MKNN/b']][] + 1
w <- lfile[['col_graphs/MKNN/w']][]
mknn <- sparseMatrix(i=a,j=b,x=w)

cellIDs.dup.idx <- which(duplicated(cellIDs))
if(length(cellIDs.dup.idx)>=1){
  cellIDs.dedup.idx <- which(!duplicated(cellIDs))
  attribute.layer <- "col"
  attribute.layer <- paste0(attribute.layer, "_attrs")
  attributes <- lfile[[attribute.layer]]$names
  attr.paths <- paste0(attribute.layer, "/", attributes)
  dims <- sapply(X = attr.paths, FUN = function(x) {
        return(length(x = lfile[[x]]$dims))
  })
  data.lst <- lapply(X = attr.paths[dims == 1], FUN = function(x) {
        return(data.frame(lfile[[x]][], stringsAsFactors = FALSE))
  })
  combined.df <- Reduce(f = cbind, x = data.lst)
  colnames(x = combined.df) <- attributes[dims == 1]
  combined.df <- combined.df[cellIDs.dedup.idx, ]
  rownames(x = combined.df) <- combined.df$CellID 
  rows.to.remove <- unname(obj = which(x = apply(X = combined.df, MARGIN = 1, FUN = function(x) {return(all(is.na(x = x)))})))
  if (length(x = rows.to.remove) > 1) {      combined.df <- combined.df[-rows.to.remove, ]}
  col.attrs <- combined.df

  # dedup matrix
  raw.matrix <- raw.matrix[, cellIDs.dedup.idx]
  knn <- knn[cellIDs.dedup.idx, cellIDs.dedup.idx]
  mknn <- mknn[cellIDs.dedup.idx, cellIDs.dedup.idx]
}else{
  col.attrs <- lfile$get.attribute.df(MARGIN = 2, col.names = "CellID")
}

row.attrs <- lfile$get.attribute.df(MARGIN = 1, row.names = "Accession")
gene.attrs <- row.attrs
cell.attrs <- col.attrs

# remove duplicated gene
idx <- which(!duplicated(row.attrs$Gene))
gene.attrs <- gene.attrs[idx, ]
rownames(gene.attrs) <- gene.attrs$Gene

raw.matrix <- raw.matrix[idx, ]
rownames(x = raw.matrix) <- rownames(gene.attrs)
colnames(x = raw.matrix) <- rownames(cell.attrs)

loom.ls <- list(raw.matrix=raw.matrix, gene.attrs=gene.attrs, cell.attrs=cell.attrs, knn=knn, mknn=mknn)
toc()

outFname <- paste(outF, "loom2mtx.rds", sep=".")
saveRDS(loom.ls, file = outFname)




