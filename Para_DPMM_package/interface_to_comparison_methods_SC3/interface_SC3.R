library(SingleCellExperiment)
library(SC3)
library(scater)
library(R.matlab)
library(clues)
rm(list=ls())
clus_num=3
Data_matlab<-readMat('~/Downloads/dpmm_subclusters/Multinomial/data_matrix_2_M_Set.mat')
data_new=Data_matlab$full.gene.trun.cell.trun.comb2.perm
label_new=Data_matlab$full.gene.trun.cell.trun.label.comb2.perm
size_matrix=dim(data_new)
gene_num=size_matrix[1]
cell_num=size_matrix[2]

rownames(data_new) <- paste('gene', 1:gene_num)
colnames(data_new) <- paste('cell', 1:cell_num)

sce <- SingleCellExperiment(
  assays = list(
    counts = as.matrix(data_new),
    logcounts = log2(as.matrix(data_new) + 1)
  ), 
  colData = colnames(data_new)
)


rowData(sce)$feature_symbol <- rownames(sce)
sce <- sce[!duplicated(rowData(sce)$feature_symbol), ]
isSpike(sce, "ERCC") <- grepl("ERCC", rowData(sce)$feature_symbol)

sce <- sc3(sce, ks = clus_num:clus_num, biology = TRUE)

extract1=colData(sce)
sc3_clusters=extract1@listData$sc3_3_clusters
sc3_clusters=as.vector(sc3_clusters)
label_new=as.vector(label_new)
score_sc3=randIndex(sc3_clusters,label_new)
writeMat(con="~/Downloads/dpmm_subclusters/Multinomial/R_matrix_to_sc3mat_cd56pl_cd4plcd25pl.mat", sc3_clusters=as.matrix(sc3_clusters), label_new=as.matrix(label_new))



