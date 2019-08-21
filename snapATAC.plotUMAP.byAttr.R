#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--meta", required=TRUE, help="load meta table")
parser$add_argument("--byattr", required=TRUE, help="load attr")
parser$add_argument("--down_num", default=10000, help="cutoff of downsample rate [default %(default)s]")
#parser$add_argument("--down_rate", default=0.2, help="cutoff of downsample rate [default %(default)s]")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

library("data.table")

inF = as.character(args$meta)
#down_rate = as.numeric(args$down_rate)
down_num = as.numeric(args$down_num)
attrF = as.character(args$byattr)
outF = as.character(args$output)

## Visulization

dat <- fread(inF, sep="\t", header=T)

#downSampleDat<- function(dat, down_rate=0.2, seed.use=10){
#    nBarcode = nrow(dat);
#    set.seed(seed.use);
#    selected_idx = sample(1:nBarcode, ceiling(nBarcode * down_rate));
#    dat = dat[selected_idx,];
#    return(dat);
#}

#if(!is.null(down_rate)){
#message("Downsample barcode by rate: ", down_rate, ", number to ", ceiling(nrow(dat) * down_rate));
#dat <- downSampleDat(dat, down_rate)
#}

downSampleDat<- function(dat, down_num=10000, seed.use=10){
    nBarcode = nrow(dat);
    set.seed(seed.use);
    selected_idx = sample(1:nBarcode, down_num);
    dat = dat[selected_idx,];
    return(dat);
}

if(nrow(dat)>down_num){
message("Downsample barcode number to ", down_num);
dat <- downSampleDat(dat, down_num)
}

attr.dat <- read.table(attrF, sep="\t",header=F)
colnames(attr.dat) <- c("datasets", "attr")
attr.len <- length(unique(attr.dat$attr))
colnames(dat)[1] <- "datasets"

data.use <- as.data.frame(merge(dat, attr.dat, by = "datasets"))
library(RColorBrewer)
colPanel <- colorRampPalette(brewer.pal(8,"Set1"))(attr.len)

pdf(paste(outF, ".umap2attr.pdf", sep=""), width=7, height=7)

graphics::plot(x=data.use[,"umap-1"],
                   y=data.use[,"umap-2"],
                    cex=0.5, 
                    pch=19, 
                    col=scales::alpha(colPanel[factor(data.use$attr)], 0.8),
                    xlab="",
                    ylab="",
                    yaxt='n',
                    xaxt='n',
                    axes=FALSE,
                   );
graphics::title(ylab="Dim-2", line=0.5, cex.lab=1.2, font.lab=2)
graphics::title(xlab="Dim-1", line=0.5, cex.lab=1.2, font.lab=2)
graphics::box(lwd=2);        
      
legend("bottomright", 
          legend = levels(factor(data.use$attr)), 
          col = colPanel,
          pch = 19,
          pt.cex=1, 
          bty = "n", 
          cex = 1, 
          text.col = "black",
          horiz = FALSE
)        
dev.off()

