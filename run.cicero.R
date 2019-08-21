#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-d", "--RData", required=TRUE, help="clustered RData with pmat")
parser$add_argument("-g", "--gs", required=TRUE, default="/projects/ps-renlab/yangli/genome/mm10/mm10.chrom.sizes", help="genome size")
parser$add_argument("-a", "--anno", required=TRUE, default="/projects/ps-renlab/yangli/genome/GRCm38.p6/cicero.gene_annotation_sample", help="genome anno")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

dat = args$RData
gs = args$gs
anno = args$anno
outPrefix = args$output

suppressPackageStartupMessages(library("SnapATAC"))
library("tictoc")
library("cicero")
library("data.table")

load(dat)

pmat <- x.sp@pmat
cellid <- paste(x.sp@sample,x.sp@barcode,sep=".")
peakid <- x.sp@peak@elementMetadata$name
peakid <- gsub(":", "_", peakid)
peakid <- gsub("-", "_", peakid)

colnames(pmat) <- peakid
rownames(pmat) <- cellid
pmat <- pmat[,grepl("chr", colnames(pmat))]
# writeMM(pmat, "./region.MOp.snap.pmat.mtx")

sparse2triples <- function(m) {
 SM = summary(m)
 D1 = m@Dimnames[[1]][SM[,1]]
 D2 = m@Dimnames[[2]][SM[,2]]
 data.frame(peak=D2, cell=D1, x=m@x)
}

tic("sparse2triples")
pmat.t <- sparse2triples(pmat)
toc()

outMtx <- paste(outPrefix, "pmat.mtx", sep=".")
fwrite(pmat.t, outMtx, col.names=F, row.names=F, quote=F, sep="\t")

input_cds <- make_atac_cds(pmat.t, binarize = TRUE)

set.seed(2019)
input_cds <- detectGenes(input_cds)
input_cds <- estimateSizeFactors(input_cds)

umap_coords <- x.sp@umap
rownames(umap_coords) <- cellid

cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = umap_coords)

sample_genome <- read.table(gs,header=F)
conns <- run_cicero(cicero_cds, sample_genome)

# 1. Finding cis-Co-accessibility Networks (CCANS)
CCAN_assigns <- generate_ccans(conns)

# 2. Cicero gene activity scores
gene_annotation_sample <- read.table(anno,header=T,sep="\t")

plot_connections(conns, "chr4", 136856906, 136909416, 
                gene_model = gene_annotation_sample, 
                coaccess_cutoff = .25, 
                connection_width = .5, 
                collapseTranscripts = "longest" )

# Make a subset of the gene annotation column containing just the coordinates and the
# gene name
gene_annotation_sub <- gene_annotation_sample[,c(1:3, 8)]

# Rename the gene symbol column to "gene"
names(gene_annotation_sub)[4] <- "gene"

input_cds <- annotate_cds_by_site(input_cds, gene_annotation_sub)

# generate unnormalized gene activity matrix
unnorm_ga <- build_gene_activity_matrix(input_cds, conns)

# make a list of num_genes_expressed
num_genes <- pData(input_cds)$num_genes_expressed
names(num_genes) <- row.names(pData(input_cds))

# normalize
cicero_gene_activities <- normalize_gene_activities(unnorm_ga, num_genes)
rda_name <- paste(outPrefix, "Rda", sep=".")
save(input_cds, cicero_cds, unnorm_ga, num_genes, conns, CCAN_assigns, cicero_gene_activities, file=rda_name)


# generate fake second unnormalized gene activity matrix
#   unnorm_ga2 <- build_gene_activity_matrix(input_cds, conns)

# if you had two datasets to normalize, you would pass both:
# num_genes should then include all cells from both sets

#   cicero_gene_activities <- normalize_gene_activities(list(unnorm_ga, unnorm_ga2),num_genes)

