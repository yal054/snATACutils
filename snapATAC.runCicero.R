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
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

RDataF = args$RData
peakF = args$peak
gs = args$gs
cellN = as.numeric(args$num)
outPrefix = args$output

suppressPackageStartupMessages(library("SnapATAC"))
suppressPackageStartupMessages(library("GenomicRanges"));
suppressPackageStartupMessages(library("tictoc"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("cicero"))

# load data
tic("load data")
x.sp <- get(load(RDataF))
toc()

# funcs
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
#umap_coords <- x.sp@umap
#rownames(umap_coords) <- cellid

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

# remove empty col
#col.idx <- which(Matrix::colSums(pmat)>0)
#pmat <- pmat[, col.idx]


tic("sparse2triples")
pmat.tri <- sparse2triples(pmat)
toc()

# 1. run cicero
tic("run cicero")
input_cds <- make_atac_cds(pmat.tri, binarize = TRUE)
set.seed(2020)
input_cds <- detect_genes(input_cds)
input_cds <- estimate_size_factors(input_cds)
input_cds <- preprocess_cds(input_cds, method = "LSI")
input_cds <- reduce_dimension(input_cds, reduction_method = 'UMAP', preprocess_method = "LSI")

umap_coords <- reducedDims(input_cds)$UMAP
#cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = umap_coords)
cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = umap_coords, k=10) 
sample_genome <- read.table(gs,header=F)
conns <- run_cicero(cicero_cds, sample_genome)
fwrite(conns, paste(outPrefix, "conns.tsv", sep="."), quote = F, sep="\t", col.names = T, row.names = F)

# 2. Finding cis-Co-accessibility Networks (CCANS)
CCAN_assigns <- generate_ccans(conns)
fwrite(CCAN_assigns, paste(outPrefix, "ccan.tsv", sep="."), quote = F, sep="\t", col.names = T, row.names = F)
toc()

# 3. gene acticity
#gene_annotation_sample <- read.table(anno,header=T,sep="\t")

# Make a subset of the gene annotation column containing just the coordinates and the
# gene name
#gene_annotation_sub <- gene_annotation_sample[,c(1:3, 8)]

# Rename the gene symbol column to "gene"
#names(gene_annotation_sub)[4] <- "gene"

#input_cds <- annotate_cds_by_site(input_cds, gene_annotation_sub)

# generate unnormalized gene activity matrix
#unnorm_ga <- build_gene_activity_matrix(input_cds, conns)
# remove any rows/columns with all zeroes
#unnorm_ga <- unnorm_ga[!Matrix::rowSums(unnorm_ga) == 0, 
#                       !Matrix::colSums(unnorm_ga) == 0]

# make a list of num_genes_expressed
#num_genes <- pData(input_cds)$num_genes_expressed
#names(num_genes) <- row.names(pData(input_cds))

# normalize
#cicero_gene_activities <- normalize_gene_activities(unnorm_ga, num_genes)
#rda_name <- paste(outPrefix, "unnorm_ga.Rda", sep=".")
#save(input_cds, cicero_cds, unnorm_ga, num_genes, conns, cicero_gene_activities, file=rda_name)
#save(input_cds, cicero_cds, conns, CCAN_assigns, file=rda_name)

