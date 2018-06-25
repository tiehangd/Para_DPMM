
library(SIMLR)
library(R.matlab)
library(igraph)
library(scater)
library(clues)

#data(BuettnerFlorian)
#data(ZeiselAmit)
#data(GliomasReduced)

set.seed(11111)
clus_num=3
Data_matlab<-readMat('~/Downloads/dpmm_subclusters/Multinomial/data_matrix_1_S_Set.mat')
data_new=Data_matlab$full.gene.trun.cell.trun.comb2.perm
label_new=Data_matlab$full.gene.trun.cell.trun.label.comb2.perm

data_new=as.matrix(data_new)
label_new=as.matrix(label_new)




example = SIMLR(X = data_new, c = clus_num, cores.ratio = 0)
SIMLR_clusters=example$y$cluster
score_ari = adjustedRand(as.vector(label_new), as.vector(SIMLR_clusters), randMethod= c("Rand", "HA", "MA", "FM", "Jaccard"))
writeMat(con="~/Downloads/dpmm_subclusters/Multinomial/R_matrix_to_SIMLR_cd56pl_cd4plcd25pl.mat", SIMLR_clusters=as.matrix(SIMLR_clusters), label_new=as.matrix(label_new))















