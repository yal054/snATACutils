#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input", required=TRUE, help="list of summit bed")
parser$add_argument("-g", "--genome", default = "mm10", help="used genome [default %(default)s]")
parser$add_argument("--blacklist", default = "/projects/ps-renlab/yangli/genome/mm10/mm10.blacklist.bed.gz", help="blacklist in bed gz")
parser$add_argument("--chromSize", default = "/projects/ps-renlab/yangli/genome/mm10/mm10.chrom.sizes", help="chrom size")
parser$add_argument("-d", "--dir", required=TRUE, help="output dir")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("GenomicRanges"))
suppressPackageStartupMessages(library("BSgenome"))
library(tictoc)

inF = args$input
blacklistF = args$blacklist
chromF = args$chromSize
outDir = args$dir
outF = args$output

##########################
# Functions 

# read bed to gr
read2gr <- function(bedF, label){
  df <- fread(bedF, sep="\t", header=F)
  colnames(df) <- c("chr", "start", "end", "name", "score")
  df$label <- label
  gr <- GRanges(
    df$chr,
    IRanges(df$start, df$end)
  );
  mcols(gr)$score <- df$score
  mcols(gr)$name <- df$name
  mcols(gr)$label <- df$label
  return(gr)
}


# extend summit to 500 bp
extendSummit <- function(gr, size=500){
  gr <- resize(gr, width=size, fix="center")
  return(gr)
}


# filter blacklist
filter4blacklist <- function(gr, blacklistF){
  black_list = read.table(blacklistF);
  black_list.gr = GRanges(
    black_list[,1],
    IRanges(black_list[,2], black_list[,3])
  );
  
  idx = queryHits(
    findOverlaps(gr, black_list.gr)
  );
  if(length(idx) > 0){
    gr = gr[-idx]
  }
  return(gr)
}


# filter non-chromosome
filter4chrom <- function(gr, chromF){
  chrom = read.table(chromF);
  chrom.gr = GRanges(
    chrom[,1],
    IRanges(0, chrom[,2])
  );
  idx = queryHits(
    findOverlaps(gr, chrom.gr, type="within")
  );
  if(length(idx) > 0){
    gr = gr[idx]
  }
  return(gr)
}


# filter N containing regions
filter4N <- function(gr, genome="mm10"){
  genome <- getBSgenome(genome)
  nucFreq <- BSgenome::alphabetFrequency(getSeq(genome, gr))
  mcols(gr)$GC <- round(rowSums(nucFreq[,c("G","C")]) / rowSums(nucFreq),4)
  mcols(gr)$N <- round(nucFreq[,c("N")] / rowSums(nucFreq),4)
  gr[which(mcols(gr)$N < 0.001)] #Remove N Containing Peaks
  return(gr)
}

# get non-overlapped regions
#' Retreive a non-overlapping set of regions from a Genomic Ranges object
#'
#' This function returns a GRanges object containing a non-overlapping set regions derived from a supplied Genomic Ranges object.
#'
#' @param gr A `GRanges` object.
#' @param by The name of a column in `mcols(gr)` that should be used to determine how overlapping regions should be resolved.
#' The resolution of overlapping regions also depends on `decreasing`. For example, if a column named "score" is used for `by`,
#' `decreasing = TRUE` means that the highest "score" in the overlap will be retained and `decreasing = FALSE` means that the
#' lowest "score" in the overlap will be retained.
#' @param decreasing A boolean value indicating whether the values in the column indicated via `by` should be ordered in decreasing
#' order. If `TRUE`, the higher value in `by` will be retained.
#' @param verbose A boolean value indicating whether the output should include extra reporting.
#' @export
nonOverlappingGR <- function(
	gr = NULL, 
	by = "score", 
	decreasing = TRUE, 
	verbose = FALSE
  ){
  
  stopifnot(by %in% colnames(mcols(gr)))

  #-----------
  # Cluster GRanges into islands using reduce and then select based on input
  #-----------
  .clusterGRanges <- function(gr = NULL, filter = TRUE, by = "score", decreasing = TRUE){
    gr <- sort(sortSeqlevels(gr))
    r <- GenomicRanges::reduce(gr, min.gapwidth=0L, ignore.strand=TRUE)
    o <- findOverlaps(gr,r, ignore.strand = TRUE)
    mcols(gr)$cluster <- subjectHits(o)
    gr <- gr[order(mcols(gr)[,by], decreasing = decreasing),]
    gr <- gr[!duplicated(mcols(gr)$cluster),]
    gr <- sort(sortSeqlevels(gr))
    mcols(gr)$cluster <- NULL
    return(gr)
  }
  
  if(verbose){
    message("Converging", appendLF = FALSE)
  }
  i <-  0
  grConverge <- gr
  while(length(grConverge) > 0){
    if(verbose){
      message(".", appendLF = FALSE)
    }
    i <-  i + 1
    grSelect <- .clusterGRanges(
      gr = grConverge, 
      filter = TRUE, 
      by = by, 
      decreasing = decreasing)

    grConverge <- subsetByOverlaps(
      grConverge,
      grSelect, 
      invert=TRUE, 
      ignore.strand = TRUE) #blacklist selected gr
    
    if(i == 1){ #if i=1 then set gr_all to clustered
      grAll <- grSelect
    
    }else{
      grAll <- c(grAll, grSelect)
    } 

  }
  message(sprintf("Converged after %s iterations!", i))

  if(verbose){
    message("\nSelected ", length(grAll), " from ", length(gr))
  }
  grAll <- sort(sortSeqlevels(grAll))

  return(grAll)

}


# normlize to score per million

norm2spm <- function(gr, by = "score"){
  mlogp = mcols(gr)[,by]
  normmlogp = 10^6 * mlogp / sum(mlogp)
  mcols(gr)$spm <- normmlogp
  return(gr)
}


################################
# working on peaks

# load summit
summitF <- read.table(inF, sep="\t", header=F)
label.lst <- as.character(summitF$V1)
file.lst <- as.character(summitF$V2)


tic("parse summit peak set")
peak.list = lapply(seq(file.lst), function(i){
  message("working on summit set for... ", label.lst[i])
  p.gr <- read2gr(file.lst[i], label=label.lst[i])
  p.gr <- extendSummit(p.gr, size=500)
  p.gr <- filter4blacklist(p.gr, blacklistF)
  p.gr <- filter4chrom(p.gr, chromF)
  p.gr <- filter4N(p.gr, genome="mm10")
  p.gr <- nonOverlappingGR(p.gr, by = "score", decreasing = TRUE)
  p.gr <- norm2spm(p.gr, by="score")
  p.gr
})
toc()


tic("write filtered & fixed peak set")
for(i in 1:length(label.lst)){
  message("write fixed & filtered peak set to: ", outDir, label.lst[i], ".filteredNfixed.peakSet")
  outPeak <- as.data.frame(peak.list[[i]])
  outFname <- paste(outDir, label.lst[i], ".filterNfixed.peakset", sep="")
  fwrite(outPeak, file=outFname, sep="\t", quote = F, col.names = T, row.names = F)
}
toc()



tic("merge to union peak list")
# merge
merged.gr <- do.call(c, peak.list)
merged.gr <- nonOverlappingGR(merged.gr, by = "spm", decreasing = TRUE)
toc()


tic("save union peak set")
outUnion <- as.data.frame(merged.gr)
outfname = paste(outDir, outF, ".filteredNfixed.union.peakSet",sep="")
fwrite(outUnion, file=outfname, sep="\t", quote = F, col.names = T, row.names = F)
toc()

