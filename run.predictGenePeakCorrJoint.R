#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-e", "--exp", required=TRUE, help="rna exp")
parser$add_argument("-p", "--pmat", required=TRUE, help="pmat norm")
parser$add_argument("-t", "--tss", required=TRUE, default="/projects/ps-renlab/yangli/genome/mm10/gencode.vM16.tss.pcg.symbol.bed", help="tss bed")
parser$add_argument("-m", "--method", default="pearson", help="# of ranging distance")
parser$add_argument("-d", "--distance", default=1000000, help="# of ranging distance")
parser$add_argument("-s", "--seed", default=2020, help="seed used for shuffle")
parser$add_argument("-c", "--cpu", default=4, help="# of cpu")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

expF = args$exp
pmatF = args$pmat
tssF = args$tss
seed = as.numeric(args$seed)
corrMethod = args$method
dis = as.numeric(args$distance)
cpus = as.numeric(args$cpu)
outF = args$output


suppressPackageStartupMessages(library("tictoc"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("reshape2"))
suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library("GenomicRanges"))
suppressPackageStartupMessages(library("foreach"))
suppressPackageStartupMessages(library("doParallel"))

########################
# load data
tic("load data")
exp.mat <- fread(expF, sep="\t", header=T)
cre.mat <- fread(pmatF, sep="\t", header=T)

# select & sort column
exp.col <- colnames(exp.mat)[-1]
cre.col <- colnames(cre.mat)[-1]
sel.col <- intersect(exp.col, cre.col)

exp.mat <- exp.mat[, match(c("gene",sel.col), colnames(exp.mat)), with=F]
cre.mat <- cre.mat[, match(c("peak",sel.col), colnames(cre.mat)), with=F]

peaks <- strsplit(cre.mat[["peak"]], ":|-|\\|")
peaks <- do.call(rbind, peaks)
peaks.gr <- GRanges(peaks[,1],
    IRanges(as.numeric(peaks[,2]), as.numeric(peaks[,3])), name=peaks[,4]
  );
rm(peaks)

queryTSSf <- tssF
queryTSS = read.table(queryTSSf);
queryTSS.gr = GRanges(queryTSS[,1],
    IRanges(queryTSS[,2], queryTSS[,3]), name=queryTSS[,4]
  );
idx.sel <- which(queryTSS.gr$name %in% exp.mat$gene)
queryTSS.gr <- queryTSS.gr[idx.sel]
toc()

#####################
# func
predictGenePeakCorr <- function(
    exp.mat, 
    cre.mat,
    cre.gr,
    method=c("pearson", "spearman"),
    gene.name=NULL,
    gene.loci=NULL,
    seed
){
    mtd <- match.arg(method)
    seed <- seed
    gene.mat = exp.mat
    gene.val = gene.mat[which(gene.mat$gene==gene.name), 2:ncol(gene.mat)];

    data.raw = cre.mat[,2:ncol(cre.mat)];
    peak.raw = cre.gr;
    
    if(!(identical(names(gene.val), colnames(data.raw)))){
        stop("colname don't match")
    }
    idy = queryHits(findOverlaps(peak.raw, gene.loci));
    ovlpn = length(idy)
    idyexcl = which(!1:length(peak.raw) %in% idy)
    if((x=length(idy))==0L){
        stop("no peaks are found within the gene flanking window")
    }else{
        data.use = data.raw[idy,];
        peak.use = peak.raw[idy];
        # remove peaks have low coverage
        idy = which(Matrix::rowSums(data.use) > 0);
        if((x=length(idy)) == 0L){
            stop("no peaks are found within the gene flanking window")
        }else if((x=length(idy)) == 1L){
        data.use = data.frame(data.use[idy,]);
        peak.use = peak.use[idy];
        }else{
        data.use = data.use[idy,];
        peak.use = peak.use[idy]; 
        }
    }

    # shuf set
    data.shuf <- data.use
    peak.shuf <- peak.use
    colheader <- colnames(data.shuf)
    #Shuffle row-wise:
    set.seed(seed)
    data.shuf <- data.shuf[sample(nrow(data.shuf)),]
    #Shuffle column-wise:
    set.seed(seed)
    data.shuf <- data.shuf[,sample(ncol(data.shuf)), with=F]
    colnames(data.shuf) <- colheader

    # rdm set
    set.seed(seed)
#    rdm <- sample(idyexcl, 1000)
    rdm <- sample(idyexcl, ovlpn)
    if((x=length(rdm))==0L){
        stop("no peaks are found within the gene flanking window")
    }else{
        data.rdm = data.raw[rdm,];
        peak.rdm = peak.raw[rdm];
        # remove peaks have low coverage
        rdm = which(Matrix::rowSums(data.rdm) > 0);
        if((x=length(rdm)) == 0L){
            stop("no peaks are found within the gene flanking window")
        }else if((x=length(rdm)) == 1L){
        data.rdm = data.frame(data.rdm[rdm,]);
        peak.rdm = peak.rdm[rdm];
        }else{
        data.rdm = data.rdm[rdm,];
        peak.rdm = peak.rdm[rdm];
        }
    }

    # rdm shuf set
    data.rdmshuf <- data.rdm
    peak.rdmshuf <- peak.rdm
    colheader <- colnames(data.rdmshuf)
    #Shuffle row-wise:
    set.seed(seed)
    data.rdmshuf <- data.rdmshuf[sample(nrow(data.rdmshuf)),]
    #Shuffle column-wise:
    set.seed(seed)
    data.rdmshuf <- data.rdmshuf[,sample(ncol(data.rdmshuf)), with=F]
    colnames(data.rdmshuf) <- colheader 
 
    corrFunc <- function(var1, var2, method) {
    result = cor.test(var1, var2, method = method)
    data.frame(result[c("estimate","p.value","statistic")], 
             stringsAsFactors=FALSE)
    }
    
#-----------------------------------
    cat("Real set: performing statitical test on... \n", file = stderr())
      peaks.id = seq(nrow(data.use));
        corr = lapply(peaks.id, function(t){
         corrFunc(log2(as.numeric(gene.val)+1), log2(as.numeric(data.use[t,])+1), mtd)
        })
      corr <- do.call(rbind, corr)
    
    peak.use$estimate <- corr[, "estimate"]
    peak.use$statistic <- corr[, "statistic"]
    peak.use$method <- mtd
    peak.use$Pval <- corr[, "p.value"]
    peak.use$FDR <- p.adjust(peak.use$Pval, method = "BH")
    peak.use$class <- "corr"
    peak.use$gene <- gene.name
    
    cat("Done ... \n", file = stderr())

#------------------------------
    cat("Shuffle set: performing statitical test on... \n", file = stderr())
      peaks.id = seq(nrow(data.shuf));
        corr = lapply(peaks.id, function(t){
         corrFunc(log2(as.numeric(gene.val)+1), log2(as.numeric(data.shuf[t,])+1), mtd)
        })
      corr <- do.call(rbind, corr)

    peak.shuf$estimate <- corr[, "estimate"]
    peak.shuf$statistic <- corr[, "statistic"]
    peak.shuf$method <- mtd
    peak.shuf$Pval <- corr[, "p.value"]
    peak.shuf$FDR <- p.adjust(peak.shuf$Pval, method = "BH")
    peak.shuf$class <- "shuf"
    peak.shuf$gene <- gene.name

    cat("Done ... \n", file = stderr())

#--------------------------------
    cat("Random set: performing statitical test on... \n", file = stderr())
      peaks.id = seq(nrow(data.rdm));
        corr = lapply(peaks.id, function(t){
         corrFunc(log2(as.numeric(gene.val)+1), log2(as.numeric(data.rdm[t,])+1), mtd)
        })
      corr <- do.call(rbind, corr)

    peak.rdm$estimate <- corr[, "estimate"]
    peak.rdm$statistic <- corr[, "statistic"]
    peak.rdm$method <- mtd
    peak.rdm$Pval <- corr[, "p.value"]
    peak.rdm$FDR <- p.adjust(peak.rdm$Pval, method = "BH")
    peak.rdm$class <- "random"
    peak.rdm$gene <- gene.name
    cat("Done ... \n", file = stderr())


#----------------------
    cat("Rdmshuf set: performing statitical test on... \n", file = stderr())
      peaks.id = seq(nrow(data.rdmshuf));
        corr = lapply(peaks.id, function(t){
         corrFunc(log2(as.numeric(gene.val)+1), log2(as.numeric(data.rdmshuf[t,])+1), mtd)
        })
      corr <- do.call(rbind, corr)
    peak.rdmshuf$estimate <- corr[, "estimate"]
    peak.rdmshuf$statistic <- corr[, "statistic"]
    peak.rdmshuf$method <- mtd
    peak.rdmshuf$Pval <- corr[, "p.value"]
    peak.rdmshuf$FDR <- p.adjust(peak.rdmshuf$Pval, method = "BH")
    peak.rdmshuf$class <- "rdmShuf"
    peak.rdmshuf$gene <- gene.name
    cat("Done ... \n", file = stderr())

#############################
    # OUTPUT
    peak.use.df <- as.data.frame(peak.use)
    peak.shuf.df <- as.data.frame(peak.shuf)
    peak.rdm.df <- as.data.frame(peak.rdm)
    peak.rdmshuf.df <- as.data.frame(peak.rdmshuf)
    peak.df <- rbind(peak.use.df, peak.shuf.df)
    peak.df <- rbind(peak.df, peak.rdm.df)
    peak.df <- rbind(peak.df, peak.rdmshuf.df) 
    
    return(peak.df);
}



##################
# analysis
# prepare for cluster 
tic("start analysis")
registerDoParallel(cpus)
genePeakPair.df <- foreach (g=queryTSS.gr@elementMetadata$name, .combine=rbind, .inorder=TRUE, .errorhandling = 'remove') %dopar% {
print(as.character(g))
#g <- "Meis2"
TSS.loci = queryTSS.gr[which(queryTSS.gr$name==g),]
pairs.df = predictGenePeakCorr(
    exp.mat=exp.mat, 
    cre.mat=cre.mat,
    cre.gr=peaks.gr,
    method = corrMethod,
    gene.name=g,
    gene.loci=resize(TSS.loci, width=dis, fix="center"),
    seed=seed
    );
tss <- paste(TSS.loci@seqnames, ":", TSS.loci@ranges@start, "-", TSS.loci@ranges@start+1, sep="")
pairs.df$tss <- tss
return(pairs.df)
}
closeAllConnections()
toc()

tic("output results")
outfname <- paste(outF,"genePeakCorr.tsv", sep=".")
fwrite(genePeakPair.df, outfname, sep="\t", col.names = T, row.names = F, quote = F)
toc()

