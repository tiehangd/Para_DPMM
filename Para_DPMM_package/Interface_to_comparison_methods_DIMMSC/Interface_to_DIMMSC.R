
library(DIMMSC)

rm(list=ls())

library(R.matlab)

Data_matlab<-readMat('../datasets/data_matrix_1_S_Set.mat')    # replace with the datasets name;

data_new=Data_matlab$full.gene.trun.cell.trun.comb2.perm

label_new=Data_matlab$full.gene.trun.cell.trun.label.comb2.perm

start_time <- Sys.time()

cluster_result <- DIMMSC(data=data_new, K=3, method_cluster_intial="kmeans", method_alpha_intial="Ronning", maxiter=200, tol=1e-4, lik.tol=1e-2)

end_time <- Sys.time()

used_time<-end_time - start_time    # program execution time;

z_r=cluster_result$mem

z_r=as.vector(z_r)

label_new=as.vector(label_new)      # ground truth label;

library(flexclust)

score_bimm=randIndex(z_r,label_new)    # Adjusted Random Index Score;

writeMat(con="R_DIMMSC_result.mat", z_r=as.matrix(z_r), label_new=as.matrix(label_new)) #  clustering result is stored in the mat file;





