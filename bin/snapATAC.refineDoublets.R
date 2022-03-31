#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input", required=TRUE, help="input doublets scores")
parser$add_argument("-m", "--meta", required=TRUE, help="input cluster meta")
parser$add_argument("-t", "--threshold", default = 0.3, help="threshold of doublet score")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

input = as.character(args$input)
meta = as.character(args$meta)
threshold = args$threshold
outF = as.character(args$output)

dat <- read.table(input, sep="\t", header=T)

doublet_scores_obs = dat$doublet_scores
doublet_scores_sim = dat$sim_doublet_scores
doublet_scores_obs = doublet_scores_obs[doublet_scores_obs > 0]
doublet_scores_sim = doublet_scores_sim[doublet_scores_sim > 0]
    
doublet_scores_obs.hist = hist(doublet_scores_obs, breaks = 50, plot=F)
doublet_scores_obs.hist$counts = log(doublet_scores_obs.hist$counts+1, 10)
    
doublet_scores_sim.hist = hist(doublet_scores_sim, breaks = 50, plot=F)
    

# plot histogram of doublet scores
pdf(paste(outF, ".refineHistDoublets.pdf",sep=""), height=3, width=8)

par(mfrow=c(1,2));    
plot(doublet_scores_obs.hist, xlim = c(0,1), col="grey", border = "grey", main = "Observed profiles", xlab = "Doublet score", ylab = "Prob. density (log10)")
abline(v = threshold, col = "black")
plot(doublet_scores_sim.hist, xlim = c(0,1), col="grey", border = "grey", main = "Simulated doublets", xlab = "Doublet score", ylab = "Prob. density")
abline(v = threshold, col = "black")
dev.off()

dat$predicted_doublets <- FALSE
dat[which(dat$doublet_scores >= threshold), "predicted_doublets"] <- TRUE
dat$threshold <- threshold


meta <- read.table(meta, header=T, sep="\t")
names(meta)[names(meta) == 'x.sp.sample'] <- 'sample'
names(meta)[names(meta) == 'x.sp.cluster'] <- 'cluster'
out <- merge(meta, dat, by=c("sample", "barcode"))

outname = paste(outF, ".refineDoublets.txt",sep="")
write.table(out, outname, col.names=T, row.names=F, sep="\t", quote=F)


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
        ncell = nrow(obj);
        data.use = values;
    }

    viz.method = match.arg(viz.method);
    viz.use = obj[,grepl(viz.method, names(obj))];

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
             main="plot doublets",
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
               main="plot doublets",
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


pdf(paste(outF, ".refineDoublets.pdf",sep=""))
# project doublet values on umap
use.value <- out$doublet_scores
plotValue(
    obj=out,
    values=use.value,
    viz.method="umap",
    point.size=0.3,
    point.color="red",
    point.shape=19,
    background.point=TRUE,
    background.point.color="grey",
    background.point.alpha=0.3,
    background.point.size=0.1,
    background.point.shape=19,
    low.value=0.0,
    high.value=1.0,
    down.sample=5000,
    seed.use=10,
    plot.nrow=4,
    plot.ncol=4,
    pdf.file.name=NULL,
    pdf.height=7,
    pdf.width=7
    );

# project protential doublets on umap
use.value <- as.numeric(out$predicted_doublets)
plotValue(
    obj=out,
    values=use.value,
    viz.method="umap",
    point.size=0.3,
    point.color="red",
    point.shape=19,
    background.point=TRUE,
    background.point.color="grey",
    background.point.alpha=0.3,
    background.point.size=0.1,
    background.point.shape=19,
    low.value=0.0,
    high.value=1.0,
    down.sample=5000,
    seed.use=10,
    plot.nrow=4,
    plot.ncol=4,
    pdf.file.name=NULL,
    pdf.height=7,
    pdf.width=7
    );

dev.off()


