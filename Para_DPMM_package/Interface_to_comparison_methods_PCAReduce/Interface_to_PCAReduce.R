rm(list=ls())

library(R.matlab)
library("pcaMethods")
library("pcaReduce")
Data_matlab<-readMat('../datasets/data_matrix_1_S_Set.mat')         # replace with the datasets name;

data_new=Data_matlab$full.gene.trun.cell.trun.comb2.perm

label_new=Data_matlab$full.gene.trun.cell.trun.label.comb2.perm

D <- log2(as.matrix(data_new) + 1)
Input <- t(D)
true_cell_cls <- label_new

start_time <- Sys.time()
Output_S <- PCAreduce(Input, nbt=2, q=4, method='S')

N <- length(Output_S)
M <- dim(Output_S[[1]])[2]
K11 <- c()


for (n in 1:N){
  cls_cell <- c()
  labels <- c()
  
  for (m in 1:M){
    cls_cell <- c(cls_cell, adjustedRandIndex(Output_S[[n]][,m], true_cell_cls))
    labels <- c(labels, length(unique(Output_S[[n]][,m])))
  }
  
  K11 <- cbind(K11, cls_cell)


}

end_time <- Sys.time()
used_time<-end_time - start_time       # program execution time;

z_pcareduce=Output_S[[1]][,3]

z_pcareduce=as.vector(z_pcareduce)       # ground truth label;

label_new=as.vector(label_new)

score_pcareduce=randIndex(z_pcareduce,label_new)        # Adjusted Random Index Score;

writeMat(con="R_PCAReduce_result.mat", z_pcareduce=as.matrix(z_pcareduce), label_new=as.matrix(label_new))        #  clustering result is stored in the mat file;    







