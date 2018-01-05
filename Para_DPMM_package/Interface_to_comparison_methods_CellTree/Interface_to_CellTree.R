library(cellTree);

rm(list=ls())

library(R.matlab)



Data_matlab<-readMat('../datasets/data_matrix_1_S_Set.mat')         # replace with the datasets name;

data_new=Data_matlab$full.gene.trun.cell.trun.comb2.perm

label_new=Data_matlab$full.gene.trun.cell.trun.label.comb2.perm

start_time <- Sys.time()

lda.results = compute.lda(data_new, k.topics=3, method="maptpx")     # remember to change number of cluster;

end_time <- Sys.time()

used_time<-end_time - start_time       # program execution time;

z_lda=lda.results$omega

z_lda=apply(z_lda, 1, function(x, k) which(x >= min(sort(x, decreasing = T)[1:1]), arr.ind = T), 1)

z_lda=unname(z_lda)

z_lda=as.vector(z_lda)

label_new=as.vector(label_new)     # ground truth label;

#score_cellTree=randIndex(z_lda,label_new)

writeMat(con="R_CellTree_result.mat", z_lda=as.matrix(z_lda), label_new=as.matrix(label_new))     #  clustering result is stored in the mat file;


