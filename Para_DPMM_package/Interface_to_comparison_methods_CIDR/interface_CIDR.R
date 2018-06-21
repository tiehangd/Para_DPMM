
library(R.matlab)
library(cidr)
rm(list=ls())

clus_num=3
Data_matlab<-readMat('~/Downloads/dpmm_subclusters/Multinomial/data_matrix_1_S_Set.mat')
data_new=Data_matlab$full.gene.trun.cell.trun.comb2.perm
label_new=Data_matlab$full.gene.trun.cell.trun.label.comb2.perm
size_matrix=dim(data_new)
gene_num=size_matrix[1]
cell_num=size_matrix[2]

rownames(data_new) <- paste('gene', 1:gene_num)
colnames(data_new) <- paste('cell', 1:cell_num)




#brainTags <- read.csv("~/Downloads/dpmm_subclusters/CIDR/brainTags.csv")
#rownames(brainTags) <- brainTags[,1]
#brainTags <- brainTags[,-1]

## Read in annotation
#info <- read.csv("~/Downloads/dpmm_subclusters/CIDR/SraRunTable.txt",sep="\t")
#cellType <- info$cell_type_s[match(colnames(brainTags),info$Sample_Name_s)]
#cellType <- factor(cellType)
#types <- levels(cellType)

scBrain <- scDataConstructor(as.matrix(data_new))
scBrain <- determineDropoutCandidates(scBrain)
scBrain <- wThreshold(scBrain)
scBrain <- scDissim(scBrain)
scBrain <- scPCA(scBrain)
scBrain <- nPC(scBrain)
#nCluster(scBrain)
scBrain <- scCluster(scBrain,nCluster = clus_num)
this_cluster=scBrain@clusters
ARI_CIDR <- adjustedRandIndex(scBrain@clusters,label_new)




