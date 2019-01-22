
####################################################################

Source code and datasets of Para-DPMM model for single cell transcriptomic clustering to reproduce results in paper "[Parallel Clustering of Single Cell Transcriptomic Data with Split-Merge Sampling on Dirichlet Process Mixtures](https://arxiv.org/pdf/1812.10048.pdf)", Author: Tiehang Duan; José P. Pinto; Xiaohui Xie;


# 1. Data Preparation:

In the datasets folder, we included the mat files that are used in the paper's experiment part. You can also prepare your own data following the procedures below.

## 1.1 
Download raw datasets from 10X genomics website (https://support.10xgenomics.com/single-cell-gene-expression/datasets) and store the files in the datasets folder;

## 1.2 
Follow the comments in "data_preparation.m" to modify the file names based on the downloaded data files;

## 1.3 
Follow the comments in "data_preparation.m" to set the number of cells (randomly selected) and top variable genes;

## 1.4 
Run "data_preparation.m";




# 2. Usage (two options):

## Option 1: Run Source Code


Enter the "Para_DPMM_Source_Code" Directory

We tested the code with Matlab2015a and gcc/4.8.4 (gcc/4.6.X and gcc/4.7.X should also work), the installation of following two libraries is needed:

1) GNU Scientific Library: http://www.gnu.org/software/gsl/
2) Eigen library: http://eigen.tuxfamily.org/


Steps:

(1) Install the packages mentioned above (for the Eigen library, you only need to place the unzipped files inside the "eigen" directory);

(2) Start Matlab, enter the "main" sub directory in Matlab;

(3) In Matlab, run compile_MEX.m;

(4) Run Para_DPMM.m, follow the guidelines given in the program (an example):

     1> Please enter dataset path:  data_matrix_1_S_Set.mat
     
     2> Please enter number of processors:  16
     
     3> Please set the value of alpha:  1
     
     4> Please enter computing time limit (seconds):  20

(5) Result: The training iteration and computation time log is saved in Para_DPMM_output.txt, the clustering result is saved in Para_DPMM_result.mat, where z is the clustering result of Para_DPMM model, label is the ground truth cluster label, AR is Adjusted Random Index, RI is Random Index benchmark, MI is "Mirkin's" index and HI is "Hubert's" index.




## Option 2: Run Compiled Package

The execution of compiled package is tested on HPC clusters with module gcc/6.1.0, gsl/2.3 and MATLAB/r2017b loaded; 

Steps:

### (1) On HPC command line:

      1> module load gcc/6.1.0

      2> module load gsl/2.3

      3> module load MATLAB/r2017b

      4> ./Para_DPMM

### (2) Follow guidelines printed in the console (an example):

      1> Please enter dataset path:  ./datasets/data_matrix_1_S_Set.mat

      2> Please enter number of processors:  16

      3> Please set the value of alpha:  1

      4> Please enter computing time limit (seconds):  40



### (3) Result

The training iteration and computation time log is saved in Para_DPMM_output.txt, the clustering result is saved in Para_DPMM_result.mat, where z is the clustering result of Para_DPMM model, label is the ground truth cluster label, AR is Adjusted Random Index, RI is Random Index benchmark, MI is "Mirkin's" index and HI is "Hubert's" index.




# 3. Comparison Methods

In the paper, we performed comparison with several current widely used single cell clustering methods. Most of the methods are available in the form of R package, and in the "interface to comparison methods" folders, we provide interface programs (written in R) to use these datasets with the available R packages for comparison. Please install the related R pacakges before using the interface programs.







Please feel free to use it for academic purposes. 


Note: The Para-DPMM project depend heavily on the open source Dirichlet Process Mixtures package(http://people.csail.mit.edu/jchang7/code.php) written by Jason Chang.

###############################################################################


