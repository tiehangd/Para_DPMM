

#######################################################################################

Matlab Package of Para_DPMM model for single cell transcriptomic clustering introduced in paper "Parallelized Inference for Single Cell Transcriptomic Clustering with Split Merge Sampling on DPMM Model", Author: Tiehang Duan; Xiaohui Xie;

The execution of this package is tested on HPC clusters with module gcc/6.1.0, gsl/2.3 and MATLAB/r2017b loaded; 



Usage:


1. Data Preparation:

1.1 Download raw datasets from 10X genomics website (https://support.10xgenomics.com/single-cell-gene-expression/datasets) and store the files in the datasets folder;

1.2 Follow the comments in "data_preparation.m" to modify the file names based on the downloaded data files;

1.3 Follow the comments in "data_preparation.m" to set the number of cells (randomly selected) and top variable genes;

1.4 Run "data_preparation.m";

In the datasets folder, we included the mat files that are used in the paper's experiment part.



2. Package Usage:

2.1 On HPC command line:

      1> module load gcc/6.1.0

      2> module load gsl/2.3

      3> module load MATLAB/r2017b

      4> ./Para_DPMM

2.2 Follow guidelines printed in the program (an example):

      1> Please enter dataset path:  ./datasets/data_matrix_1_S_Set.mat

      2> Please enter number of processors:  16

      3> Please set the value of alpha:  1

      4> Please enter computing time limit (seconds):  40



3. Result Files:

      The training iteration and computation time log is saved in Para_DPMM_output.txt, the clustering result is saved in Para_DPMM_result.mat, where z is the clustering result of Para_DPMM model, label is the ground truth cluster label, AR is Adjusted Random Index, RI is Random Index benchmark, MI is "Mirkin's" index and HI is "Hubert's" index.


4. Comparison Methods:

      In the paper, we performed comparison with several current state of art single cell clustering methods. Most of the methods are available in the form of R package, and in the "interface to comparison methods" folders, we provide interface programs (written in R) to use these datasets with the available R packages for comparison. Please install the related R pacakges before using the interface programs.







Please feel free to use it for academic purposes.


Note: This package utilizes multiple functions from the open sourced Dirichlet Process Mixtures package(http://people.csail.mit.edu/jchang7/code.php) written by Jason Chang. 

###############################################################################################


