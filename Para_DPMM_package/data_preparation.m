
CD4_plus_CD25_plus=mmread('./datasets/cd4_plus_cd25_plus_regu_t/hg19/matrix.mtx');    % change the file name based on the data file downloaded from 10X genomics website;
CD4_plus_CD45_plus_CD25_minu=mmread('./datasets/cd4_plus_cd45_plus_cd25_minu_naive_t/hg19/matrix.mtx');    % change the file name based on the data file downloaded from 10X genomics website;
CD8_plus_CD45ra_plus=mmread('./datasets/cd8_plus_cd45ra_plus_naive_cytotoxic_t/hg19/matrix.mtx');    % change the file name based on the data file downloaded from 10X genomics website;

comb2=[CD4_plus_CD25_plus,CD4_plus_CD45_plus_CD25_minu,CD8_plus_CD45ra_plus];
label1_leng=size(CD4_plus_CD25_plus,2);
label2_leng=size(CD4_plus_CD45_plus_CD25_minu,2);
label3_leng=size(CD8_plus_CD45ra_plus,2);

label1=ones(label1_leng,1);
label2=ones(label2_leng,1)*2;
label3=ones(label3_leng,1)*3;


label_comb2=[label1;label2;label3];
data_num=label1_leng+label2_leng+label3_leng;
index_perm=randperm(data_num);
comb2_perm=comb2(:,index_perm);
label_comb2_perm=label_comb2(index_perm,:);

num_genes=1000;     % set the number of top variable genes;
num_cells_1=400;     % set the number of cells (first type);
num_cells_2=400;      % set the number of cells (second type);
num_cells_3=400;       % set the number of cells (third type);

rand_cell_1=floor(rand(num_cells_1,1)*label1_leng)+1;
rand_cell_2=floor(rand(num_cells_2,1)*label2_leng)+1;
rand_cell_3=floor(rand(num_cells_3,1)*label3_leng)+1;


cell_truncate_CD4_plus_CD25_plus=CD4_plus_CD25_plus(:,rand_cell_1);
cell_truncate_CD4_plus_CD45_plus_CD25_minu=CD4_plus_CD45_plus_CD25_minu(:,rand_cell_2);
cell_truncate_CD8_plus_CD45ra_plus=CD8_plus_CD45ra_plus(:,rand_cell_3);

cell_truncate_comb2=[cell_truncate_CD4_plus_CD25_plus,cell_truncate_CD4_plus_CD45_plus_CD25_minu,cell_truncate_CD8_plus_CD45ra_plus];


cell_truncate_label_1=ones(num_cells_1,1);
cell_truncate_label_2=ones(num_cells_2,1)*2;
cell_truncate_label_3=ones(num_cells_3,1)*3;

cell_truncate_label_comb2=[cell_truncate_label_1;cell_truncate_label_2;cell_truncate_label_3];

trun_num=num_cells_1+num_cells_2+num_cells_3;
trun_index_perm=randperm(trun_num);

cell_trun_comb2_perm=cell_truncate_comb2(:,trun_index_perm);
cell_trun_label_comb2_perm=cell_truncate_label_comb2(trun_index_perm,:);

save('data_matrix_interm.mat','cell_trun_comb2_perm','cell_trun_label_comb2_perm','num_genes','num_cells_1','num_cells_2','num_cells_3');
clear
load('data_matrix_interm.mat');


gene_std=std(cell_trun_comb2_perm,0,2);

[std_sorted,Index_sorted]=sort(gene_std,'descend');


gene_index_selected=Index_sorted(1:num_genes,:);

gene_trun_cell_trun_comb2_perm=cell_trun_comb2_perm(gene_index_selected,:);
gene_trun_cell_trun_label_comb2_perm=cell_trun_label_comb2_perm;

full_gene_trun_cell_trun_comb2_perm=full(gene_trun_cell_trun_comb2_perm);
full_gene_trun_cell_trun_label_comb2_perm=full(gene_trun_cell_trun_label_comb2_perm);

save('./datasets/data_matrix_1_S_Set.mat','full_gene_trun_cell_trun_comb2_perm','full_gene_trun_cell_trun_label_comb2_perm');


% "full_gene_trun_cell_trun_comb2_perm" stores the prepared UMI matrix,
% 'full' stands for full matrix; 'gene_trun_cell_trun' means the cells and 
%  genes in the matrix are truncated to the number specified in the program;
% 'comb' means the data combines three cell types;
% 'perm' means the cells are randomly permuted;
% "full_gene_trun_cell_trun_label_comb2_perm" stores the relating ground truth label;



