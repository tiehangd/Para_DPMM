function AR = Para_DPMM()

diary('para_dpmm_out.txt');
diary on;

prompt = 'Please enter dataset path:  ';
str = input(prompt,'s');
load(str);
gene_trun_cell_trun_comb2_perm=sparse(full_gene_trun_cell_trun_comb2_perm);

gene_trun_cell_trun_label_comb2_perm=sparse(full_gene_trun_cell_trun_label_comb2_perm);

if ~isdeployed 
    
    addpath('../common');
  
end

% set the random number generation seed for reproducible data
RandStream.setGlobalStream(RandStream('mt19937ar','Seed', 1));

	

% run the sampler
initialClusters = 1;
dispOn = false;
prompt = 'Please enter number of processors:  ';
proc_num=input(prompt);
numProcessors = proc_num;
useSuperclusters = false;
approximateSampling = false;
prompt = 'Please set the value of alpha:  ';
alpha_value=input(prompt);
alpha = alpha_value;
prompt = 'Please enter computing time limit (seconds):  ';
end_time=input(prompt);
endtime = end_time;
numits = 1000;

% uncomment the algorithm you want to run

z=run_dpmnmm_subclusters(gene_trun_cell_trun_comb2_perm, initialClusters, dispOn, numProcessors, ...
    useSuperclusters, approximateSampling, alpha, endtime, numits);

z=z+1;
[AR,RI,MI,HI]=RandIndex(z,gene_trun_cell_trun_label_comb2_perm);

disp('ARI');
disp(AR);
disp('RI');
disp(RI);
disp('MI');
disp(MI);
disp('HI');
disp(HI);

label=gene_trun_cell_trun_label_comb2_perm;
save('para_dpmm_result.mat','z','label','AR','RI','MI','HI');

diary off;

end


