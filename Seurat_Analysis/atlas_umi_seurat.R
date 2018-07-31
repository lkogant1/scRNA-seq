#load the dataset
#seurat_object@data - retrieve data from seurat object
#nUMI = sum of non-normalized values within a cell
atlas_umi_a <- read.table("/Volumes/SeagateEpv/Yan_Lab/AHaber_ARegev_Nature_2017/Ncbi_Data/Supplemental_Files/Zipped/Duplicate/GSE92332_atlas_UMIcounts/GSE92332_atlas_UMIcounts.txt", header=TRUE, fill=FALSE)
dim(atlas_umi_a)
#15971  7216
atlas_umi_a[1:5,1:5]
#-------------------------------------------------------------------------------------------
#create seurat object
library(SingleCellExperiment)
library(Seurat)
library(mclust)
library(dplyr)
atlas_umi_seurat_a <- CreateSeuratObject(
  raw.data = atlas_umi_a,
  min.cells = 3, 
  min.genes = 200
)
#----------------------------------------------------------------------
#expression qc
#visualize gene and umi counts
#VlnPlot(
#  object = atlas_umi_seurat_a, 
#  features.plot = c("nGene", "nUMI"), 
#  nCol = 2
#)
#
#relationship between umi counts and gene counts
#GenePlot(
#  object = atlas_umi_seurat_a, 
#  gene1 = "nUMI", 
#  gene2 = "nGene"
#)

#filter cells with outlier number of read counts
atlas_umi_seurat_b <- FilterCells(
  object = atlas_umi_seurat_a, 
  subset.names = c("nUMI"), 
  high.thresholds = c(35000)
)

#----------------------------------------------------------------------
#normalize
#normalizes the gene expression measurements for each cell by the total expression,
#multiplies this by a scale factor (10,000 by default),
#and log-transforms the result
atlas_umi_seurat_c <- NormalizeData(
  object = atlas_umi_seurat_b, 
  normalization.method = "LogNormalize", 
  scale.factor = 10000
)
#------------------------------------------------------
#HVG
atlas_umi_seurat_d <- FindVariableGenes(
  object = atlas_umi_seurat_c,
  mean.function = ExpMean, 
  dispersion.function = LogVMR, 
  x.low.cutoff = 0.0125, 
  x.high.cutoff = 3, 
  y.cutoff = 0.5
)

length(x = atlas_umi_seurat_d@var.genes)

#-------------------------------------------------------------------------------------------
#remove confounding factors
atlas_umi_seurat_e <- ScaleData(
  object = atlas_umi_seurat_d, 
  vars.to.regress = c("nUMI")
)
#-------------------------------------------------------------------------------
#PCA
atlas_umi_seurat_f <- RunPCA(
  object = atlas_umi_seurat_e, 
  pc.genes = atlas_umi_seurat_e@var.genes, 
  do.print = TRUE, 
  pcs.print = 1:5, 
  genes.print = 5
)

#PrintPCA(object =atlas_umi_seurat_f, pcs.print = 1:5, genes.print = 5, use.full = FALSE)

#VizPCA(object = atlas_umi_seurat_f, pcs.use = 1:2)

#PCAPlot(object = atlas_umi_seurat_f, dim.1 = 1, dim.2 = 2)

#PCHeatmap(
#  object = atlas_umi_seurat_f, 
#  pc.use = 13:18, 
#  cells.use = 500, 
#  do.balanced = TRUE, 
#  label.columns = FALSE,
#  use.full = FALSE
#)
#significant PCs
# Significant PCs will show a 
#strong enrichment of genes with low p-values
atlas_umi_seurat_g <- JackStraw(
  object = atlas_umi_seurat_f, 
  num.replicate = 100)
#  do.print = FALSE)

#JackStrawPlot(object = atlas_umi_seurat_g, PCs = 1:20)
#PCElbowPlot(object =atlas_umi_seurat_g)


#-------------------------------------------------------------------------------
#clustering
atlas_umi_seurat_h <- FindClusters(
  object = atlas_umi_seurat_g, 
  reduction.type = "pca", 
  dims.use = 1:13, 
  resolution = 1.0, 
  print.output = 0, 
  save.SNN = TRUE
)

#PrintFindClustersParams(object = atlas_umi_seurat_h)
#table(atlas_umi_seurat_h@ident)

#tsne plot
atlas_umi_seurat_i <- RunTSNE(
  object = atlas_umi_seurat_h,
  dims.use = 1:13,
  do.fast = TRUE
)

#saveRDS(atlas_umi_seurat_i, file = "/Volumes/SeagateEpv/Yan_Lab/AHaber_ARegev_Nature_2017/Ncbi_Data/Supplemental_Files/Zipped/Duplicate/GSE92332_atlas_UMIcounts/atlas_umi_seurat_i.rds")

TSNEPlot(object = atlas_umi_seurat_i)
#markers0 <- FindMarkers(atlas_umi_seurat_i, 0)
#markers1 <- FindMarkers(atlas_umi_seurat_i, 1)
#markers2 <- FindMarkers(atlas_umi_seurat_i, 2)
#markers3 <- FindMarkers(atlas_umi_seurat_i, 3)
#markers4 <- FindMarkers(atlas_umi_seurat_i, 4)
#markers5 <- FindMarkers(atlas_umi_seurat_i, 5)
#markers6 <- FindMarkers(atlas_umi_seurat_i, 6)
#markers7 <- FindMarkers(atlas_umi_seurat_i, 7)
#markers8 <- FindMarkers(atlas_umi_seurat_i, 8)
#markers9 <- FindMarkers(atlas_umi_seurat_i, 9)
#markers10 <- FindMarkers(atlas_umi_seurat_i, 10)
#markers11 <- FindMarkers(atlas_umi_seurat_i, 11)
#markers12 <- FindMarkers(atlas_umi_seurat_i, 12)
#markers13 <- FindMarkers(atlas_umi_seurat_i, 13)
#markers14 <- FindMarkers(atlas_umi_seurat_i, 14)


#
allmarkers <- FindAllMarkers(
  object = atlas_umi_seurat_i, 
  only.pos = TRUE, 
  min.pct = 0.25, 
  thresh.use = 0.25
)

#write.table(markers0,file="/Volumes/SeagateEpv/Yan_Lab/AHaber_ARegev_Nature_2017/Ncbi_Data/Supplemental_Files/Zipped/Duplicate/GSE92332_atlas_UMIcounts/markers0.txt",quote=F)
#write.table(markers1,file="/Volumes/SeagateEpv/Yan_Lab/AHaber_ARegev_Nature_2017/Ncbi_Data/Supplemental_Files/Zipped/Duplicate/GSE92332_atlas_UMIcounts/markers1.txt",quote=F)
#write.table(markers2,file="/Volumes/SeagateEpv/Yan_Lab/AHaber_ARegev_Nature_2017/Ncbi_Data/Supplemental_Files/Zipped/Duplicate/GSE92332_atlas_UMIcounts/markers2.txt",quote=F)
#write.table(markers3,file="/Volumes/SeagateEpv/Yan_Lab/AHaber_ARegev_Nature_2017/Ncbi_Data/Supplemental_Files/Zipped/Duplicate/GSE92332_atlas_UMIcounts/markers3.txt",quote=F)
#write.table(markers4,file="/Volumes/SeagateEpv/Yan_Lab/AHaber_ARegev_Nature_2017/Ncbi_Data/Supplemental_Files/Zipped/Duplicate/GSE92332_atlas_UMIcounts/markers4.txt",quote=F)
#write.table(markers5,file="/Volumes/SeagateEpv/Yan_Lab/AHaber_ARegev_Nature_2017/Ncbi_Data/Supplemental_Files/Zipped/Duplicate/GSE92332_atlas_UMIcounts/markers5.txt",quote=F)
#write.table(markers6,file="/Volumes/SeagateEpv/Yan_Lab/AHaber_ARegev_Nature_2017/Ncbi_Data/Supplemental_Files/Zipped/Duplicate/GSE92332_atlas_UMIcounts/markers6.txt",quote=F)
#write.table(markers7,file="/Volumes/SeagateEpv/Yan_Lab/AHaber_ARegev_Nature_2017/Ncbi_Data/Supplemental_Files/Zipped/Duplicate/GSE92332_atlas_UMIcounts/markers7.txt",quote=F)
#write.table(markers8,file="/Volumes/SeagateEpv/Yan_Lab/AHaber_ARegev_Nature_2017/Ncbi_Data/Supplemental_Files/Zipped/Duplicate/GSE92332_atlas_UMIcounts/markers8.txt",quote=F)
#write.table(markers9,file="/Volumes/SeagateEpv/Yan_Lab/AHaber_ARegev_Nature_2017/Ncbi_Data/Supplemental_Files/Zipped/Duplicate/GSE92332_atlas_UMIcounts/markers9.txt",quote=F)
#write.table(markers10,file="/Volumes/SeagateEpv/Yan_Lab/AHaber_ARegev_Nature_2017/Ncbi_Data/Supplemental_Files/Zipped/Duplicate/GSE92332_atlas_UMIcounts/markers10.txt",quote=F)
#write.table(markers11,file="/Volumes/SeagateEpv/Yan_Lab/AHaber_ARegev_Nature_2017/Ncbi_Data/Supplemental_Files/Zipped/Duplicate/GSE92332_atlas_UMIcounts/markers11.txt",quote=F)
#write.table(markers12,file="/Volumes/SeagateEpv/Yan_Lab/AHaber_ARegev_Nature_2017/Ncbi_Data/Supplemental_Files/Zipped/Duplicate/GSE92332_atlas_UMIcounts/markers12.txt",quote=F)
#write.table(markers13,file="/Volumes/SeagateEpv/Yan_Lab/AHaber_ARegev_Nature_2017/Ncbi_Data/Supplemental_Files/Zipped/Duplicate/GSE92332_atlas_UMIcounts/markers13.txt",quote=F)
#write.table(markers14,file="/Volumes/SeagateEpv/Yan_Lab/AHaber_ARegev_Nature_2017/Ncbi_Data/Supplemental_Files/Zipped/Duplicate/GSE92332_atlas_UMIcounts/markers14.txt",quote=F)

#write.table(allmarkers,file="/Volumes/SeagateEpv/Yan_Lab/AHaber_ARegev_Nature_2017/Ncbi_Data/Supplemental_Files/Zipped/Duplicate/GSE92332_atlas_UMIcounts/allmarkers.txt",quote=F)


#heatmap of top five markers
top3 <- allmarkers %>% group_by(cluster) %>% top_n(3, avg_logFC)

DoHeatmap(
  object = atlas_umi_seurat_i, 
  genes.use = top3$gene, 
  slim.col.label = TRUE, 
  remove.key = TRUE)
#-------------------------------------------------------------------------------
#?
#atlas_umi_seurat_f <- RunPCA(
#  object = atlas_umi_seurat_e, 
#  pc.genes = atlas_umi_seurat_e@var.genes, 
#  do.print = TRUE, 
#  pcs.print = 1:5, 
#  genes.print = 5
#)

#atlas_umi_seurat_z <- ProjectPCA(object = atlas_umi_seurat_f, do.print = FALSE)

