
load('data_matrix_1_S_Set.mat');
[gene_num,cell_num]=size(full_gene_trun_cell_trun_comb2_perm);
new_matrix=zeros(gene_num,cell_num);

for i =1:cell_num
    this_col=full_gene_trun_cell_trun_comb2_perm(:,i);
    total_count=sum(this_col);
    sample_size=floor(total_count*7/8);
    dis_prob=this_col/total_count;
    x = discretesample(dis_prob, sample_size);
    new_matrix(:,i)=countmember(1:gene_num,x);
end

full_gene_trun_cell_trun_comb2_perm=new_matrix;

save('../Multinomial/S_Set_seq_dep_2100.mat','full_gene_trun_cell_trun_comb2_perm','full_gene_trun_cell_trun_label_comb2_perm');