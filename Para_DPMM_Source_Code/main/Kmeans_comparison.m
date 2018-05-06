
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Kmeans method on the same task for comparison

load('data_matrix_to_R_cd56pl_cd4plcd25pl.mat');   % Please change the path and file name accordingly before execution;
k=3;                                               % Change the number of clusters;
gene_trun_cell_trun_comb2_perm=sparse(full_gene_trun_cell_trun_comb2_perm);

gene_trun_cell_trun_label_comb2_perm=sparse(full_gene_trun_cell_trun_label_comb2_perm);



T_gene_trun_cell_trun_comb2_perm=transpose(gene_trun_cell_trun_comb2_perm);
[cell_num,gene_num]=size(T_gene_trun_cell_trun_comb2_perm);

for i=1:cell_num
    T_gene_trun_cell_trun_comb2_perm(i,:)=(T_gene_trun_cell_trun_comb2_perm(i,:)-min(T_gene_trun_cell_trun_comb2_perm(i,:)))/(max(T_gene_trun_cell_trun_comb2_perm(i,:))-min(T_gene_trun_cell_trun_comb2_perm(i,:)));
    
end

tStart = tic;
[kmeans_idx,kmeans_C]=kmeans(T_gene_trun_cell_trun_comb2_perm,k);
tElapsed = toc(tStart);


[AR,RI,MI,HI]=RandIndex(kmeans_idx,gene_trun_cell_trun_label_comb2_perm);

