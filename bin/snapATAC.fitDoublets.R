#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input", required=TRUE, help="input RData")
parser$add_argument("-m", "--meta", required=TRUE, help="input cluster meta")
parser$add_argument("-c", "--cutoff", default=0.25, help="default cutoff")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

input = as.character(args$input)
meta = as.character(args$meta)
defaultCutoff = as.numeric(args$cutoff)
outF = as.character(args$output)

suppressPackageStartupMessages(library("SnapATAC"))
suppressPackageStartupMessages(library("GenomicRanges"))
suppressPackageStartupMessages(library("plyr"))
library("tictoc")
library("umap")
suppressPackageStartupMessages(library("mixtools"))

#---------------------
# func
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


#---------------------
load(input)
dat <- read.table(meta, sep="\t", header=T)

doublet_scores_obs = as.numeric(as.character(dat$doublet_scores))
doublet_scores_sim = as.numeric(as.character(dat$sim_doublet_scores))
doublet_scores_obs = doublet_scores_obs[doublet_scores_obs > 0]
doublet_scores_sim = doublet_scores_sim[doublet_scores_sim > 0]

doublet_scores_obs.hist = hist(doublet_scores_obs, breaks = 50, plot=F)
doublet_scores_obs.hist$counts = log(doublet_scores_obs.hist$counts+1, 10)
doublet_scores_sim.hist = hist(doublet_scores_sim, breaks = 50, plot=F)

mixmdl = normalmixEM(doublet_scores_sim[is.finite(doublet_scores_sim)], k=2)

dist1.lambda <- mixmdl$lambda[1]
dist2.lambda <- mixmdl$lambda[2]
dist1.mu <- mixmdl$mu[1]
dist1.sigma <- mixmdl$sigma[1]
dist2.mu <- mixmdl$mu[2]
dist2.sigma <- mixmdl$sigma[2]

index.lower <- which.min(mixmdl$mu)

find.cutoff <- function(model, proba=0.5, i=index.lower) {
    ## Cutoff such that Pr[drawn from bad component] == proba
    f <- function(x) {
        proba - (model$lambda[i]*dnorm(x, model$mu[i], model$sigma[i]) /
                     (model$lambda[1]*dnorm(x, model$mu[1], model$sigma[1]) + model$lambda[2]*dnorm(x, model$mu[2], model$sigma[2])))
        }
        return(uniroot(f=f, lower=0.2, upper=1)$root)  # Careful with division by zero if changing lower and upper
}

t <- try(find.cutoff(mixmdl, proba=0.5), silent=TRUE)
if("try-error" %in% class(t)){
cutoffs <- defaultCutoff
}else{
cutoffs <- find.cutoff(mixmdl, proba=0.5)
}

type1error <- pnorm(cutoffs, mean = dist1.mu, sd = dist1.sigma, lower.tail = FALSE, log.p = FALSE)
type2error <- pnorm(cutoffs, mean = dist2.mu, sd = dist2.sigma, lower.tail = TRUE, log.p = FALSE)

threshold <- cutoffs

staout <- c(dist1.lambda, dist1.mu, dist1.sigma, dist2.lambda, dist2.mu, dist2.sigma, type1error, type2error, threshold)
names(staout) <- c("lambda1", "mu1", "sigma1", "lambda2", "mu2", "sigma2", "type1error", "type2error", "threshold")
staout <- t(as.data.frame(staout))
staname = paste(outF, ".fitDoublets.sta",sep="")
write.table(staout, staname, col.names=T, row.names=F, sep="\t", quote=F)

# plot histogram of doublet scores
pdf(paste(outF, ".fitHistDoublets.pdf",sep=""), height=6, width=8)

par(mfrow=c(2,2));
plot(doublet_scores_obs.hist, xlim = c(0,1), col="grey", border = "grey", main = "Observed profiles", xlab = "Doublet score", ylab = "Prob. density (log10)")
abline(v = threshold, col = "black")

plot(doublet_scores_sim.hist, xlim = c(0,1), col="grey", border = "grey", main = "Simulated doublets", xlab = "Doublet score", ylab = "Prob. density")
abline(v = threshold, col = "black")

plot(mixmdl,which=2)

plot(doublet_scores_sim.hist, xlim = c(0,1), col="grey", border = "grey", main = "Simulated doublets", xlab = "Doublet score", ylab = "Prob. density")
abline(v=cutoffs, col=c("red", "blue"), lty=2)
dev.off()

dat$doublet_scores <- as.numeric(as.character(dat$doublet_scores))
dat$predicted_doublets <- FALSE
dat[which(dat$doublet_scores >= threshold), "predicted_doublets"] <- TRUE
dat$threshold <- threshold
dat <- join(x.sp@metaData, dat, by=c("sample","barcode"))

outname = paste(outF, ".fitDoublets.txt",sep="")
write.table(dat, outname, col.names=T, row.names=F, sep="\t", quote=F)

pdf(paste(outF, ".fitDoublets.pdf",sep=""))
use.value <- dat$doublet_scores
plotValue(
    obj=x.sp,
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

use.value <- as.numeric(dat$predicted_doublets)
plotValue(
    obj=x.sp,
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

plotViz(
    obj= x.sp,
    method="umap",
    main="Cluster",
    point.color=x.sp@cluster,
    point.size=0.2,
    point.shape=19,
    text.add=TRUE,
    text.size=1,
    text.color="black",
    down.sample=10000,
    legend.add=FALSE
  );
dev.off()

