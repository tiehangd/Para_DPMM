
# Seurat is a very comprehensive analysis package with many additional analysis tools,
# so the interface implementation is a little tricky compared to other methods.
# please change three places when changing datasets, two input, one output.

library(Seurat)
library(dplyr)
library(Matrix)
library(R.matlab)


data_dire=c("./seurat/cd4_plus_cd25_plus/hg19/","./seurat/cd8_plus_cd45_plus/hg19/")
pbmc.data <- Read10X(data.dir = data_dire)


Data_matlab<-readMat('../datasets/data_matrix_1_S_Set.mat')

data_new_matlab=Data_matlab$full.gene.trun.cell.trun.comb2.perm

label_new=Data_matlab$full.gene.trun.cell.trun.label.comb2.perm

data_new_matlab_sparse= as(data_new_matlab, "sparseMatrix") 

gene_num=data_new_matlab_sparse@Dim[[1]]
cell_num=data_new_matlab_sparse@Dim[[2]]




pbmc.data@i=data_new_matlab_sparse@i
pbmc.data@p=data_new_matlab_sparse@p
pbmc.data@Dim=data_new_matlab_sparse@Dim
pbmc.data@x=data_new_matlab_sparse@x
pbmc.data@Dimnames[[1]]<-pbmc.data@Dimnames[[1]][1:gene_num]
pbmc.data@Dimnames[[2]]<-as.character(label_new)



dense.size <- object.size(x = as.matrix(x = pbmc.data))
sparse.size <- object.size(x = pbmc.data)
pbmc <- CreateSeuratObject(raw.data = pbmc.data, min.cells = 3, min.genes = 200, project = "10X_PBMC")
trun_cell_num=pbmc@raw.data@Dim[[2]]

label_new<-data.matrix(pbmc@raw.data@Dimnames[[2]])
label_new<-as.numeric(label_new)



#########################################





data_dire=c("./seurat/cd4_plus_cd25_plus/hg19/","./seurat/cd8_plus_cd45_plus/hg19/")
pbmc.data <- Read10X(data.dir = data_dire)


Data_matlab<-readMat('../datasets/data_matrix_1_S_Set.mat')

data_new_matlab=Data_matlab$full.gene.trun.cell.trun.comb2.perm


data_new_matlab_sparse= as(data_new_matlab, "sparseMatrix") 

gene_num=data_new_matlab_sparse@Dim[[1]]
cell_num=data_new_matlab_sparse@Dim[[2]]



pbmc.data@i=data_new_matlab_sparse@i
pbmc.data@p=data_new_matlab_sparse@p
pbmc.data@Dim=data_new_matlab_sparse@Dim
pbmc.data@x=data_new_matlab_sparse@x
pbmc.data@Dimnames[[1]]<-pbmc.data@Dimnames[[1]][1:gene_num]
pbmc.data@Dimnames[[2]]<-pbmc.data@Dimnames[[2]][1:cell_num]

pbmc <- CreateSeuratObject(raw.data = pbmc.data, min.cells = 3, min.genes = 200, project = "10X_PBMC")



##########################################




mito.genes <- grep(pattern = "^MT-", x = rownames(x = pbmc@data), value = TRUE)
percent.mito <- Matrix::colSums(pbmc@raw.data[mito.genes, ])/Matrix::colSums(pbmc@raw.data)


pbmc <- AddMetaData(object = pbmc, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = pbmc, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)

par(mfrow = c(1, 2))
GenePlot(object = pbmc, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = pbmc, gene1 = "nUMI", gene2 = "nGene")



pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

pbmc <- ScaleData(object = pbmc, vars.to.regress = c("nUMI", "percent.mito"))

pbmc <- RunPCA(object = pbmc, pc.genes = pbmc@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)

pbmc <- ProjectPCA(object = pbmc, do.print = FALSE)

pbmc <- JackStraw(object = pbmc, num.replicate = 100, do.print = FALSE)

start_time <- Sys.time()

pbmc <- FindClusters(object = pbmc, reduction.type = "pca", dims.use = 1:10, resolution = 0.8, print.output = 0, save.SNN = TRUE)

end_time <- Sys.time()

used_time= end_time - start_time

z_seurat<-data.matrix(pbmc@ident)




writeMat(con="R_Seurat_result.mat", z_seurat=as.matrix(z_seurat), label_new=as.matrix(label_new))


z_seurat=strtoi(z_seurat)
z_seurat=as.vector(z_seurat)

label_new=as.vector(label_new)

library(flexclust)

score_seurat=randIndex(z_seurat,label_new)




pbmc <- RunTSNE(object = pbmc, dims.use = 1:10, do.fast = TRUE, check_duplicates = FALSE)

TSNEPlot(object = pbmc)








