#!/usr/bin/env Rscript

# read results from ksklearn

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input", required=TRUE, help="input statH")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

data <- read.table(args$input,sep="\t",head=F)

staoutmx <- data.frame(row.names=c("Min","Q1","Median","Mean","Q3","Max","TopWhisker","BottomWhisker","Box1","Box2","Box3","UpWhisker","DnWhisker"))
for (i in c(5,6,7)){
x <- data[,i]
boxMx <- matrix(summary(x))
rownames(boxMx) <- c("Min","Q1","Median","Mean","Q3","Max")
iqr <- IQR(x)
q1 <- summary(x)[2]
q3 <- summary(x)[5]
TopWhisker <- min(max(x), q3 + 1.5 * iqr)
BottomWhisker <- max(min(x), q1 - 1.5 * iqr)
Box1 <- boxMx["Q1",]
Box2 <- boxMx["Median",] - boxMx["Q1",]
Box3 <- boxMx["Q3",] - boxMx["Median",]
UpWhisker <- TopWhisker - boxMx["Q3",]
DnWhisker <- boxMx["Q1",] - BottomWhisker
boxMx <- rbind(boxMx,TopWhisker,BottomWhisker,Box1,Box2,Box3,UpWhisker,DnWhisker)
colnames(boxMx) <- i
staoutmx <- cbind(staoutmx,boxMx)
}
colnames(staoutmx) <- c("contributes","sparseness","entropy")

cat(median(data$V6),"\n")

k <- max(data$V3)
n <- nrow(data)
normInfoGain = 1 - sum(data$V7) / (n * log2(k))
cat(normInfoGain)

write.table(staoutmx, file=paste(args$output,".box.sta",sep=''), sep="\t", quote=F, col.names=T, row.names=T)

