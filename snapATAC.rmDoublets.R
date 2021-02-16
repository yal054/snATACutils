#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--RData", required=TRUE, help="load pmat RData")
parser$add_argument("--mat", default = "gmat", required=TRUE, help="pmat/bmat/gmat")
parser$add_argument("--path_to_python", default = NULL, help="path.to.python")
parser$add_argument("--path_to_condaenv", default = NULL, help="path.to.condaenv")
parser$add_argument("--path_to_venv", default = NULL, help="path.to.venv")
parser$add_argument("--rate", default = 0.05, help="expected.doublet.rate")
parser$add_argument("--min_counts", default = 3, help="min.counts")
parser$add_argument("--min_cells", default = 3L, help="min.cells")
parser$add_argument("--min_gene_variability_pctl", default = 85,  help="min.gene_variability_pctl")
parser$add_argument("--n_prin_comps", default = 30L, help="n.prin_comps")
parser$add_argument("--black_list", default="/projects/ps-renlab/yangli/genome/mm10/mm10.blacklist.bed.gz", help="black list file")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

suppressPackageStartupMessages(library("SnapATAC"))
suppressPackageStartupMessages(library("GenomicRanges"))
library("tictoc")
library("umap")

#.libPaths("/home/yangli1/apps/anaconda3/envs/r_env/lib/R/library")

RDataF = as.character(args$RData)
mat = args$mat
path_to_python = args$path_to_python
path_to_condaenv = args$path_to_condaenv
path_to_venv = args$path_to_venv
expected_doublet_rate = as.numeric(args$rate)
min_counts = as.numeric(args$min_counts)
min_cells = args$min_cells
min_gene_variability_pctl = as.numeric(args$min_gene_variability_pctl)
n_prin_comps = args$n_prin_comps
black_listF = args$black_list
outF = as.character(args$output)

#========================#
#------ rmDoublets ------#

# usage example:
# library("SnapATAC")
# load("/projects/ps-renlab/yangli/projects/CEMBA/01.doublets/1A/CEMBA180226_1A.snap.mat3.RData")
# res <- rmDoublets(x.sp, mat="pmat", path.to.venv = "/mnt/silencer2/home/yangli/apps/anaconda3/envs/venv_scrublet/")


#' Detect potential doublets using scrublet (https://github.com/AllonKleinLab/scrublet)
#'
#' Identify potential doublets using scrublet (https://github.com/AllonKleinLab/scrublet)
#' This function requires Python or virtual/Conda environments, and "scrublet" preinstalled.
#' 
#' @param obj A snap object.
#' @param mat Matrix to use for finding differential features c("bmat", "pmat", "gmat").
#' @param path.to.python Path to Python excutable file.
#' @param path.to.condaenv Path to Python in Conda environments
#' @param path.to.vene Path to Python in virtual environments
#' @param expected.doublet.rate float, optional (default: 0.06)
#'            The estimated doublet rate for the experiment.
#' @param min.counts float, optional (default: 3L)
#' Used for gene filtering prior to PCA. Genes expressed at fewer than 
#' `min_counts` in fewer than `min_cells` (see below) are excluded.
#' @param min.cells int, optional (default: 3L)
#' Used for gene filtering prior to PCA. Genes expressed at fewer than 
#' `min_counts` (see above) in fewer than `min_cells` are excluded.
#' @param min.gene_variability_pctl float, optional (default: 85.0)
#' Used for gene filtering prior to PCA. Keep the most highly variable genes
#' (in the top min_gene_variability_pctl percentile), as measured by 
#' the v-statistic [Klein et al., Cell 2015].
#' @param n.prin_comps int, optional (default: 30L)
#' Number of principal components used to embed the transcriptomes prior
#' to k-nearest-neighbor graph construction.
#' @return Return a list contains "scrubObj", "doublet_scores" and "predicted_doublets"
#' @examples
#' res <- rmDoublets(x.sp, mat="pmat")
#' @export
rmDoublets <- function(obj, mat, ...) {
  UseMethod("rmDoublets", obj);
}

#' @export
rmDoublets.default <- function(
    obj, 
    mat = c("pmat", "bmat", "gmat"),
    path.to.python = NULL,
    path.to.condaenv = NULL,
    path.to.venv = NULL,
    expected.doublet.rate = 0.06,
    min.counts = 3, 
    min.cells = 3L, 
    min.gene_variability_pctl = 85, 
    n.prin_comps = 30L
){
  
  cat("Epoch: checking input parameters ... \n", file = stderr())
    if(missing(obj)){
        stop("obj is missing")
    }else{
        if(!is.snap(obj)){
            stop("obj is not a snap object");
        }        
        if((x=nrow(obj))==0L){
            stop("obj is empty");
        }
    }
  
  ncell = nrow(obj);
        
    mat = match.arg(mat);
    mat.use = methods::slot(obj, mat);
    if((x=nrow(mat.use)) == 0L){
        stop("mat is empty, add matrix to snap object first")
    }
        
    if((x=nrow(mat.use)) != ncell){
        stop("mat has different length with cell number, re-add this matrix to snap object");
    }
  
  cat("Epoch: loading python environments and packages ... \n", file = stderr())
  # load library and python env
  library("reticulate")
  if(!is.null(path.to.python)){
    reticulate::use_python(path.to.python)
    cat("use the Python located in:", path.to.python, "\n")
  }else if(!is.null(path.to.condaenv)){
    reticulate::use_virtualenv(path.to.condaenv)
    cat("use the Python in Conda environment:", path.to.condaenv, "\n")
  }else if(!is.null(path.to.venv)){
    reticulate::use_virtualenv(path.to.venv)
    cat("use the Python in virtual environments:", path.to.venv, "\n")
  }else{
    cat("By default, uses the version of Python found on your PATH", "\n")
  }

  # convert paramter to python integer
  setSessionTimeLimit(cpu = Inf, elapsed = Inf)
  mat.use = r_to_py(mat.use)
    min.cells = as.integer(min.cells);
    n.prin_comps = as.integer(n.prin_comps);
  
  # load python packages
  os <- import("os", convert = FALSE)
  scipy.io <- import("scipy.io", convert = FALSE)
  np <- import("numpy", convert = FALSE)
  scr <- import("scrublet", convert = FALSE)
  
  # identify potential doublets
  cat("Epoch: identify potential doublets ... \n", file = stderr())
  scrub = scr$Scrublet(mat.use, 
                       expected_doublet_rate=expected.doublet.rate)
  
  out = scrub$scrub_doublets(min_counts = min.counts, 
                             min_cells = min.cells, 
                             min_gene_variability_pctl = min.gene_variability_pctl, 
                             n_prin_comps = n.prin_comps)
  # convert to R explicitly at the end
  scrub <- py_to_r(scrub)
  out <- py_to_r(out)
  
  # summary the results
  cat("Epoch: summary doublets detection results ... \n", file = stderr())
  message("Automatically set threshold at doublet score = ", round(scrub$threshold_,2));
  message("Detected doublet rate = ", round(scrub$detected_doublet_rate_ * 100,2), "%");
  message("Estimated detectable doublet fraction = ", round(scrub$detectable_doublet_fraction_ * 100, 2), "%");
  message("Overall doublet rate:");
  message("\tExpected   = ", round(scrub$expected_doublet_rate * 100, 2), "%");
  message("\tEstimated  = ", round(scrub$overall_doublet_rate_ * 100, 2), "%");

  # output
  outList <- list("scrub" = scrub, "doublet_scores" = out[[1]], "predicted_doublets" = out[[2]]);
  
  return(outList);
  }


#=======================================#
#------ function plotHistDoublets ------#

# usage example:
# plotHistDoublets(res$scrub)


#' Plot histogram of doublet scores
#'
#' Plot histogram of doublet scores for observed transcriptomes and simulated doublets 
#' The histogram for simulated doublets is useful for determining the correct doublet 
#' score threshold. 
#' 
#' @param obj A scrub object, first elements from the output list of "rmDoublets"
#' @param pdf.file.name pdf file name to save the plot [NULL].
#' @param pdf.width the width of the graphics region in inches [8].
#' @param pdf.height the height of the graphics region in inches [3].
#' @param ... Arguments passed to plot method.
#' @example
#' plotHistDoublets(res$scrub)
#' @export
plotHistDoublets <- function(scrubObj, pdf.file.name, pdf.width, pdf.height, ...) {
  UseMethod("plotHistDoublets", scrubObj);
}

#' @export
plotHistDoublets.default <- function(
  scrubObj, 
    rm.zeros=TRUE,
    pdf.file.name=NULL,
    pdf.width=8, 
    pdf.height=3, 
    ...
){
  if(missing(scrubObj)){
        stop("scrubObj is missing");
    }else{
        if(length(scrubObj$doublet_scores_obs_) == 0L){
            stop("doublet_scores is not found, please rerun function rmDoublets");
        }
        if(length(scrubObj$doublet_scores_sim_) == 0L){
            stop("doublet_scores is not found, please rerun function rmDoublets");
        }
    }
    
    if(!is.null(pdf.file.name)){
        if(file.exists(pdf.file.name)){
            warning("pdf.file already exists");
            file.remove(pdf.file.name);
        }else{
            if(!file.create(pdf.file.name)){
                stop("cannot create pdf.file, not a directory")                
            }
            file.remove(pdf.file.name);
        }    
        grDevices::pdf(pdf.file.name,width=pdf.width,height=pdf.height); 
    }
    
  threshold = scrubObj$threshold_
    doublet_scores_obs = scrubObj$doublet_scores_obs_
    doublet_scores_sim = scrubObj$doublet_scores_sim_
    doublet_scores_obs = doublet_scores_obs[doublet_scores_obs > 0]
    doublet_scores_sim = doublet_scores_sim[doublet_scores_sim > 0]
    
    doublet_scores_obs.hist = hist(doublet_scores_obs, breaks = 50, plot=F)
    doublet_scores_obs.hist$counts = log(doublet_scores_obs.hist$counts+1, 10)
    
    doublet_scores_sim.hist = hist(doublet_scores_sim, breaks = 50, plot=F)
    
    par(mfrow=c(1,2));
    
    plot(doublet_scores_obs.hist, xlim = c(0,1), col="grey", border = "grey", main = "Observed profiles", xlab = "Doublet score", ylab = "Prob. density (log10)")
    abline(v = threshold, col = "black")
    
    plot(doublet_scores_sim.hist, xlim = c(0,1), col="grey", border = "grey", main = "Simulated doublets", xlab = "Doublet score", ylab = "Prob. density")
    abline(v = threshold, col = "black")
    
    if(!is.null(pdf.file.name)){
        grDevices::dev.off()        
    }
    graphics::par(mfrow=c(1,1));
  
}



library("scales")

plotValue <- function(obj, values, viz.method, point.size, point.color, point.shape,  background.point, background.point.color, background.point.alpha, background.point.size, background.point.shape, low.value, high.value, down.sample, seed.use, plot.nrow, plot.ncol, pdf.file.name, pdf.height, pdf.width,...){
  UseMethod("plotValue", obj);
}

#' @export
plotValue.default <- function(
    obj, 
    values,
    viz.method=c("tsne", "umap"),
    point.size=0.5,
    point.color="red",
    point.shape=19,
    background.point=TRUE,
    background.point.color="grey",
    background.point.alpha=0.3,
    background.point.size=0.5,
    background.point.shape=19,
    low.value=0.0,
    high.value=1.0,
    down.sample=5000,
    seed.use=10,
    plot.nrow=3,
    plot.ncol=3,
    pdf.file.name=NULL, 
    pdf.height=7, 
    pdf.width=7,
    ...
){
    
    if(missing(obj)){
        stop("obj is missing")
    }else{
        if(!is(obj, "snap")){
            stop("obj is not a snap object")
        }
        ncell = nrow(obj);
        data.use = values;
        if((x=length(data.use)) == 0L){
            stop("value is empty, input a list of value")
        }
        if((x = length(data.use))!= ncell){
            stop("value has different number of rows from cell number, add value again")
        }
    }
    
    viz.method = match.arg(viz.method);
    viz.use = methods::slot(obj, viz.method);
    
    if((x = nrow(viz.use)) == 0L){
        stop("visulization matrix is empty, run runViz first")
    }
    
    
    if(down.sample < ncell){
        set.seed(seed.use);
        idx = sort(sample(seq(ncell), down.sample));
        data.use = data.use[idx];
        viz.use = viz.use[idx,];
    }
    
    if(!is.null(pdf.file.name)){
        if(file.exists(pdf.file.name)){
            warning("pdf.file already exists");
            file.remove(pdf.file.name);
        }else{
            if(!file.create(pdf.file.name)){
                stop("cannot create pdf.file, not a directory")                
            }
            file.remove(pdf.file.name);
        }    
        pdf(pdf.file.name,width=pdf.width,height=pdf.height); 
    }
    
    y = data.use;
    if(background.point){
        plot(viz.use, 
                  main="plotValue",
             col=scales::alpha(background.point.color, background.point.alpha), 
             cex=background.point.size,
             pch=background.point.shape,
                yaxt='n', 
                xaxt="n",
                xlab="", 
                ylab="",
                ...
             );    
        points(viz.use, 
                    col=alpha(point.color, y), 
                    cex=point.size,
                    pch=point.shape
                  );    
    }else{
        plot(viz.use, 
               main="plotValue",
                 col=alpha(point.color, y), 
                 cex=point.size,
                 pch=point.shape,
                    yaxt='n', 
                    xaxt="n",
                    xlab="", 
                    ylab="",
                 ...
                 );    
        }

    if(!is.null(pdf.file.name)){
        dev.off()        
    }
    par(mfrow=c(1,1));
}



#---------------------------

tic("loadSnap")
load(RDataF)
x.sp = addBmatToSnap(
    x.sp,
    bin.size=5000,
    num.cores=1
    );
x.sp = addGmatToSnap(
    x.sp,
    num.cores=1
    );
toc()


tic("makeBinary")
x.sp = makeBinary(x.sp, mat="bmat");
toc()

tic("filterBins")
black_list = read.table(black_listF);
black_list.gr = GRanges(
    black_list[,1],
    IRanges(black_list[,2], black_list[,3])
  );
idy = queryHits(
    findOverlaps(x.sp@feature, black_list.gr)
  );
if(length(idy) > 0){
    x.sp = x.sp[,-idy, mat="bmat"];
  };


chr.exclude = seqlevels(x.sp@feature)[grep("random|chrM", seqlevels(x.sp@feature))];
idy = grep(paste(chr.exclude, collapse="|"), x.sp@feature);
if(length(idy) > 0){
    x.sp = x.sp[,-idy, mat="bmat"]
  };

bin.cov = log10(Matrix::colSums(x.sp@bmat)+1);
bin.cutoff = quantile(bin.cov[bin.cov > 0], 0.95);
idy = which(bin.cov <= bin.cutoff & bin.cov > 0);
x.sp = x.sp[, idy, mat="bmat"];
x.sp
toc()

tic("runDiffusionMaps")
x.sp = runDiffusionMaps(
    obj= x.sp,
    input.mat="bmat",
    num.eigs=n_prin_comps
  );
toc()

tic("runKNN")
x.sp = runKNN(
    obj=x.sp,
    eigs.dims=1:n_prin_comps,
    k=15
    );
toc()

## R-igraph
tic("runCluster")
x.sp = runCluster(obj=x.sp,
    tmp.folder=tempdir(),
    louvain.lib="R-igraph",
    seed.use=10
    );
toc()

tic("runViz_umap")
x.sp = runViz(
    obj=x.sp,
    tmp.folder=tempdir(),
    dims=2,
    eigs.dims=1:n_prin_comps,
    method="umap",
    seed.use=10
  );
toc()


tic("rmDoublets")
res <- rmDoublets(
    x.sp,
    mat = mat,
    path.to.python = path_to_python,
    path.to.condaenv = path_to_condaenv,
    path.to.venv = path_to_venv,
    expected.doublet.rate = expected_doublet_rate,
    min.counts = min_counts,
    min.cells = min_cells,
    min.gene_variability_pctl = min_gene_variability_pctl,
    n.prin_comps = n_prin_comps
)
toc()

# doublet scores
#head(res$doublet_scores)

# predicted_doublets
#head(res$predicted_doublets)

# plot histogram of doublet scores
pdf(paste(outF, ".plotHistDoublets.pdf",sep=""), height=3, width=8)
plotHistDoublets(res$scrub)
dev.off()

set.seed(2020)
doublet_sim <- sample(res$scrub$doublet_scores_sim_, length(res$scrub$doublet_scores_obs_))
x.sp@metaData$doublet_scores <- res$doublet_scores
outTable <- as.data.frame(cbind(x.sp@sample, x.sp@barcode, res$doublet_scores, res$predicted_doublets, doublet_sim, res$scrub$threshold_))
colnames(outTable) <- c("sample","barcode", "doublet_scores", "predicted_doublets", "sim_doublet_scores", "threshold")
outTablename = paste(outF, ".rmDoublets.txt",sep="")
write.table(outTable, outTablename, col.names=T, row.names=F, sep="\t", quote=F)


tic("saveSnap")
outfname = paste(outF, ".qc.cluster.RData",sep="")
#outfname = paste(outF, ".rmDoublets.RData", sep="")
save(x.sp, file=outfname)
toc()

