#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--RData", required=TRUE, help="clustered RData with pmat")
parser$add_argument("-n", "--num", default = 1000, help="num of cell used for downsample")
parser$add_argument("-p", "--peak", default = "/projects/ps-renlab/yangli/projects/CEMBA/01.joint_dat/rs1cemba/rs1cemba.CREs.final.bed", help="used peaks in bed format")
parser$add_argument("-g", "--gs", default="/projects/ps-renlab/yangli/genome/mm10/mm10.chrom.sizes.lite", help="genome size")
parser$add_argument("-c", "--cpus", default=4, help="# of cpus")
parser$add_argument("-s", "--seed", required=TRUE, help="seed used for shuffle")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

RDataF = args$RData
peakF = args$peak
gs = args$gs
cellN = as.numeric(args$num)
cpus = as.numeric(args$cpus)
seed = as.numeric(args$seed)
outPrefix = args$output

suppressPackageStartupMessages(library("SnapATAC"))
suppressPackageStartupMessages(library("GenomicRanges"));
suppressPackageStartupMessages(library("tictoc"))
suppressPackageStartupMessages(library("cicero"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("foreach"))
suppressPackageStartupMessages(library("doParallel"))

# load data
tic("load data")
x.sp <- get(load(RDataF))
toc()

# funcs
genBsp <- function(a, b, n, seed=2020){
  set.seed(seed)
  i <- sample(1:a, n, replace = T)
  set.seed(seed)
  j <- sample(1:b, n, replace = T)
  x <- rep(1,n)
  mx <- sparseMatrix(i=i, j=j, x=x, dims = c(a,b))
  return(mx)
}

sparse2triples <- function(m) {
 SM = summary(m)
 D1 = m@Dimnames[[1]][SM[,1]]
 D2 = m@Dimnames[[2]][SM[,2]]
 data.frame(peak=D2, cell=D1, x=m@x)
}


# filter pmat
tic("filter CREs")
peak_list = fread(peakF);
peak_list.gr = GRanges(
    peak_list[[1]],
    IRanges(peak_list[[2]], peak_list[[3]])
  );
idy = queryHits(
    findOverlaps(x.sp@peak, peak_list.gr)
  );
if(length(idy) > 0){
    x.sp = x.sp[,idy, mat="pmat"];
  };
toc()


# get pmat
pmat <- x.sp@pmat
cellid <- paste(x.sp@sample,x.sp@barcode,sep=".")
peakid <- x.sp@peak@elementMetadata$name
peakid <- gsub(":", "_", peakid)
peakid <- gsub("-", "_", peakid)
peaks <- x.sp@peak$name

colnames(pmat) <- peakid
rownames(pmat) <- cellid
rm(x.sp)


# downsample to 1000 cells
if(nrow(pmat) > cellN){
  set.seed(2020)
  idx.dn <- sample(1:nrow(pmat), cellN)
  pmat <- pmat[idx.dn, ]
#  umap_coords <- umap_coords[idx.dn, ]
}

tic("generate shuffle/permutation set")
# generate shuf pmat
# 1. permutation
m <- dim(pmat)[1]
n <- dim(pmat)[2]
k <- length(pmat@x)
pmat.perm <- genBsp(m,n,k,seed) # permutation

colnames(pmat.perm) <- colnames(pmat)
rownames(pmat.perm) <- rownames(pmat)

# 2. row shuffle (nuclei)
pmat.rowShuf <- pmat
set.seed(2020)
pmat.rowShuf <- pmat.rowShuf[sample(nrow(pmat.rowShuf)),]
colnames(pmat.rowShuf) <- colnames(pmat)
rownames(pmat.rowShuf) <- rownames(pmat)

# 3. col shuffle (peaks)
pmat.colShuf <- pmat
set.seed(2020)
pmat.colShuf <- pmat.colShuf[,sample(ncol(pmat.colShuf)), with=F]
colnames(pmat.colShuf) <- colnames(pmat)
rownames(pmat.colShuf) <- rownames(pmat)

# 4. both
pmat.shuf <- pmat
set.seed(2020)
pmat.shuf <- pmat.shuf[sample(nrow(pmat.shuf)),]
set.seed(2020)
pmat.shuf <- pmat.shuf[,sample(ncol(pmat.shuf)), with=F]
colnames(pmat.shuf) <- colnames(pmat)
rownames(pmat.shuf) <- rownames(pmat)
rm(pmat)
toc()


tic("sparse2triples")
pmat.perm.tri <- sparse2triples(pmat.perm)
rm(pmat.perm)
pmat.rowShuf.tri <- sparse2triples(pmat.rowShuf)
rm(pmat.rowShuf)
pmat.colShuf.tri <- sparse2triples(pmat.colShuf)
rm(pmat.colShuf)
pmat.shuf.tri <- sparse2triples(pmat.shuf)
rm(pmat.shuf)
toc()

# 1. run cicero
tic("run cicero")
# make parallel
pmat.tri.ls <- list(pmat.rowShuf.tri, pmat.colShuf.tri, pmat.shuf.tri, pmat.perm.tri)

registerDoParallel(cpus)
res <- foreach(i=1:4, .inorder=TRUE) %dopar% {
pmat.tri <- pmat.tri.ls[[i]]
input_cds <- make_atac_cds(pmat.tri, binarize = TRUE)
set.seed(2020)
input_cds <- detect_genes(input_cds)
input_cds <- estimate_size_factors(input_cds)
input_cds <- preprocess_cds(input_cds, method = "LSI")
input_cds <- reduce_dimension(input_cds, reduction_method = 'UMAP', preprocess_method = "LSI")

umap_coords <- reducedDims(input_cds)$UMAP
cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = umap_coords)
#cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = umap_coords, k=10) # for Gnu
sample_genome <- read.table(gs,header=F)
conns <- run_cicero(cicero_cds, sample_genome)
return(conns)
}
closeAllConnections()
toc()

outPrefix.ls <- paste(outPrefix, c("rowShuf", "colShuf", "shuf", "perm"), "conns.tsv", sep=".")
for( i in 1:length(res)){
fwrite(res[[i]], outPrefix.ls[[i]], quote = F, sep="\t", col.names = T, row.names = F)
}

