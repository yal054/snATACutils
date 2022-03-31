#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input", required=TRUE, help="load snap RData")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

metaF <- args$input
outF = args$output

suppressPackageStartupMessages(library("RColorBrewer"))
suppressPackageStartupMessages(library("viridis"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("reshape2"))
suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library("ComplexHeatmap"))
suppressPackageStartupMessages(library("circlize"))
library("wesanderson")
library("dendsort")
sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))


metaF <- fread(metaF, sep="\t", header=T, na.strings=c("","NA"), quote="")

#-----------------------------------------------------
# t1: table with 2 columns: coembed labels, raw labels
# t2: table with 2 columns: coembed labels, raw labels
cal_ovlpScore <- function(t1, t2){
  t1.table <- table(t1)
  t2.table <- table(t2)
  t1.pct <- apply(t1.table, 2, function(x){x/sum(x)})
  t2.pct <- apply(t2.table, 2, function(x){x/sum(x)})
  t1.labels <- colnames(t1.pct)
  t2.labels <- colnames(t2.pct)
  ovlpScore.df <- data.frame(anno1=as.character(), anno2=as.character(), ovlpScore=as.numeric())
  for(t1.label in t1.labels){
    for(t2.label in t2.labels){
      t1.pct.df <- data.frame(t1.pct[,t1.label])
      colnames(t1.pct.df) <- "t1"
      t1.pct.df$ident <- rownames(t1.pct.df)
      t2.pct.df <- data.frame(t2.pct[,t2.label])
      colnames(t2.pct.df) <- "t2"
      t2.pct.df$ident <- rownames(t2.pct.df)
      comp.df <- join(t1.pct.df, t2.pct.df, by="ident", type="full")
      comp.df[is.na(comp.df)] <- 0
      comp.df$ident <- NULL
      comp.df <- t(comp.df)
      ovlpScore <- sum(apply(comp.df, 2, min))
      out <- data.frame(anno1=t1.label, anno2=t2.label, ovlpScore=ovlpScore)
      ovlpScore.df <- rbind(ovlpScore.df, out)
    }
  }
  return(ovlpScore.df)
}

#####################
# sub_class level
#####################

# calculate overlap
ident2rna <- data.frame(idents=metaF$coembed.idents, rna_label=metaF$ClusterName)
ident2rna <- ident2rna[complete.cases(ident2rna), ]

ident2atac <- data.frame(idents=metaF$coembed.idents, atac_label=metaF$L2cluster)
ident2atac <- ident2atac[complete.cases(ident2atac), ]

ovlpScore.df <- cal_ovlpScore(ident2rna, ident2atac)
write.table(ovlpScore.df, paste(outF,"L2cluster.overlapScore.tsv", sep="."), sep="\t", col.names = T, row.names = F, quote = F)

# plot heatmap
mapSubclass <- ovlpScore.df
colnames(mapSubclass) <- c("ClusterName", "L2cluster", "ovlpScore")
mapSubclass <- dcast(mapSubclass, L2cluster~ClusterName, value.var = "ovlpScore", fun.aggregate = identity, fill = 0)
mapSubclass <- as.data.frame(mapSubclass)
rowN <- mapSubclass$L2cluster
rownames(mapSubclass) <- rowN
mapSubclass$L2cluster <- NULL
mapSubclass <- as.matrix(sapply(mapSubclass, as.numeric))
rownames(mapSubclass) <- rowN
write.table(mapSubclass, paste(outF,"L2cluster.overlapScore.heatmap.tsv", sep="."), sep="\t", col.names = T, row.names = T, quote = F)

row_ha <- data.frame(unique(metaF[,c("L1cluster", "L2cluster", "L1Color", "L2Color")]))
row_ha <- row_ha[complete.cases(row_ha),]
rownames(row_ha) <- row_ha$L2cluster
ha_row <- HeatmapAnnotation(L2=row_ha[["L2cluster"]], L1=row_ha[["L1cluster"]], 
                            which = "row",
                            col = list(
                              L2 = setNames(unique(row_ha[,"L2Color"]), unique(row_ha[,"L2cluster"])), 
                              L1 = setNames(unique(row_ha[,"L1Color"]), unique(row_ha[,"L1cluster"]))
                            ))

#col_ha <- data.frame(unique(metaF[,c("class_label", "subclass_label", "class_color", "subclass_color")]))
#col_ha <- col_ha[complete.cases(col_ha),]
#rownames(col_ha) <- col_ha$subclass_label
#ha_col <- HeatmapAnnotation(subclass_label=col_ha[["subclass_label"]], class_label=col_ha[["class_label"]], 
#                            which = "column",
#                            col = list(
#                              subclass_label = setNames(unique(col_ha[,"subclass_color"]), unique(col_ha[,"subclass_label"])),
#                              class_label = setNames(unique(col_ha[,"class_color"]), unique(col_ha[,"class_label"]))
#                            ))

hclust_rows <- sort_hclust(hclust(dist(mapSubclass), method = "ward.D2"))
hclust_cols <- sort_hclust(hclust(dist(t(mapSubclass)),method = "ward.D2"))
col_fun = colorRamp2(seq(0, 1, length = 50), wes_palette("Zissou1", 50, type = "continuous"))
#col_fun = wes_palette("Zissou1", 50, type = "continuous")

pdf(paste(outF, ".L2cluster.overlapScore.heatmap.pdf",sep=""), width = 20, height = 10)
Heatmap(mapSubclass, 
        col = col_fun,
        cluster_columns = hclust_cols, 
        cluster_rows = hclust_rows,
        show_row_names = T,
        row_names_gp = gpar(fontsize = 8),
        show_column_names = T,
#        top_annotation = ha_col,
        left_annotation = ha_row
        )
dev.off()


#######################
# sub_cluster level
#######################

# calculate overlap
ident2rna <- data.frame(idents=metaF$coembed.idents, rna_label=metaF$ClusterName)
ident2rna <- ident2rna[complete.cases(ident2rna), ]

ident2atac <- data.frame(idents=metaF$coembed.idents, atac_label=metaF$L3cluster)
ident2atac <- ident2atac[complete.cases(ident2atac), ]

ovlpScore.df <- cal_ovlpScore(ident2rna, ident2atac)
write.table(ovlpScore.df, paste(outF,"L3cluster.overlapScore.tsv", sep="."), sep="\t", col.names = T, row.names = F, quote = F)


# plot heatmap
mapSubcluster <- ovlpScore.df
colnames(mapSubcluster) <- c("ClusterName", "L3cluster", "ovlpScore")
mapSubcluster <- dcast(mapSubcluster, L3cluster~ClusterName, value.var = "ovlpScore", fun.aggregate = identity, fill = 0)
mapSubcluster <- as.data.frame(mapSubcluster)
rowN <- mapSubcluster$L3cluster
rownames(mapSubcluster) <- rowN
mapSubcluster$L3cluster <- NULL
mapSubcluster <- as.matrix(sapply(mapSubcluster, as.numeric))
rownames(mapSubcluster) <- rowN
write.table(mapSubcluster, paste(outF,"L3cluster.overlapScore.heatmap.tsv", sep="."), sep="\t", col.names = T, row.names = T, quote = F)

row_ha <- data.frame(unique(metaF[,c("L1cluster", "L2cluster", "L3cluster", "L1Color", "L2Color", "L3Color")]))
row_ha <- row_ha[complete.cases(row_ha),]
rownames(row_ha) <- row_ha$L3cluster
ha_row <- HeatmapAnnotation(L3=row_ha[["L3cluster"]], L2=row_ha[["L2cluster"]], L1=row_ha[["L1cluster"]],
                            which = "row",
                            col = list(
                              L3 = setNames(unique(row_ha[,"L3Color"]), unique(row_ha[,"L3cluster"])),
                              L2 = setNames(unique(row_ha[,"L2Color"]), unique(row_ha[,"L2cluster"])),
                              L1 = setNames(unique(row_ha[,"L1Color"]), unique(row_ha[,"L1cluster"]))
                            ))

#col_ha <- data.frame(unique(metaF[,c("class_label", "subclass_label", "cluster_label", "class_color", "subclass_color", "cluster_color")]))
#col_ha <- col_ha[complete.cases(col_ha),]
#rownames(col_ha) <- col_ha$subcluster_label
#ha_col <- HeatmapAnnotation(subcluster_label=col_ha[["cluster_label"]], subclass_label=col_ha[["subclass_label"]], class_label=col_ha[["class_label"]],
#                            which = "column",
#                            col = list(
#                              subcluster_label = setNames(unique(col_ha[,"cluster_color"]), unique(col_ha[,"cluster_label"])),
#                              subclass_label = setNames(unique(col_ha[,"subclass_color"]), unique(col_ha[,"subclass_label"])),
#                              class_label = setNames(unique(col_ha[,"class_color"]), unique(col_ha[,"class_label"]))
#                            ))

hclust_rows <- sort_hclust(hclust(dist(mapSubcluster), method = "ward.D2"))
hclust_cols <- sort_hclust(hclust(dist(t(mapSubcluster)),method = "ward.D2"))
col_fun = colorRamp2(seq(0, 1, length = 50), wes_palette("Zissou1", 50, type = "continuous"))
#col_fun = wes_palette("Zissou1", 50, type = "continuous")

pdf(paste(outF, ".L3cluster.overlapScore.heatmap.pdf",sep=""), width = 20, height = 10)
Heatmap(mapSubcluster,
        col = col_fun,
        cluster_columns = hclust_cols,
        cluster_rows = hclust_rows,
        show_row_names = T,
        row_names_gp = gpar(fontsize = 8),
        show_column_names = T,
#        top_annotation = ha_col,
        left_annotation = ha_row
        )
dev.off()



