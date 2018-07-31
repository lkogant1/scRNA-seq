#load 
atlas_a <- read.table("/Volumes/SeagateEpv/Yan_Lab/AHaber_ARegev_Nature_2017/Ncbi_Data/Supplemental_Files/Zipped/Duplicate/GSE92332_atlas_UMIcounts/GSE92332_atlas_UMIcounts.txt", header=TRUE, fill=FALSE)
atlas_a <- as.matrix(atlas_a)
dim(atlas_a)
#15971  7216
atlas_a[1:5,1:5]
#check ERCC spike-ins
#atlas_a[rownames(atlas_a) %in% "Ercc5"]

#-------------------------------------------------------------------------------------------
#exclude cells with fewer than 800 detected genes
col_sums_a <- as.matrix(colSums(atlas_a))
dim(col_sums_a)
#7216    1
col_sums_a[1:5,1]
x <- as.matrix(col_sums_a[(col_sums_a[,1]<=800),])
dim(x)
x
#none
x <- as.matrix(col_sums_a[(col_sums_a[,1]>800),])
dim(x)
#7216    1
#-----------------------------------------------------------------
#normalize for differences in coverage
#divide each cell by column sum 
#https://www.r-bloggers.com/function-apply-tip-1/
atlas_b <- apply(atlas_a,2,function(x){x/sum(x)})
typeof(atlas_b)
dim(atlas_b)
#15971  7216
#validate: column sums now equal to 1
col_sums_b <- as.matrix(colSums(atlas_b))
dim(col_sums_b)
col_sums_b[1:5,1]

#multiply by 10,000 to create TPM-like values
atlas_c <-apply(atlas_b,2,function(x){x*10000})
typeof(atlas_c)
dim(atlas_c)
atlas_b[1:5,1:5]
atlas_c[1:5,1:5]

#log2transform
atlas_d <- apply(atlas_c,2,function(x){log2(x+1)})
dim(atlas_d)
atlas_d[1:5,1:5]
atlas_c[1:5,1:5]

#write data
#write.table(atlas_d,"/Volumes/SeagateEpv/Yan_Lab/AHaber_ARegev_Nature_2017/Ncbi_Data/Supplemental_Files/Zipped/Duplicate/atlas_d.txt",sep="\t",quote=FALSE)

##-------------------------------------------------------------------------------------------
#remove batch effects
#install packages
#source("https://bioconductor.org/biocLite.R")
#biocLite("sva") #part of bioconductor package
library(sva)

#prep batch file
#extract header
#head -1 GSE92332_atlas_UMIcounts.txt > atlas_batch_1
#tr "\t" "\n" < atlas_batch_1 > atlas_batch_2
#wc -l atlas_batch_2
#7216
#duplicate column
#awk 'BEGIN{OFS=","}{print $1, $1}' atlas_batch_2 > atlas_batch_3
#cut -d'_' -f1,2,3 atlas_batch_3 > atlas_batch_4
#awk -F, 'BEGIN{OFS=FS}{gsub("B","",$2); print}' atlas_batch_4 > atlas_batch_5
#cut -d',' -f2 atlas_batch_5 > atlas_batch_6
#convert column to row
#awk 'BEGIN{ORS=" "}{print}' atlas_batch_6 > atlas_batch.txt
#cp atlas_batch_6 atlas_batch.txt


batch = read.table("/Volumes/SeagateEpv/Yan_Lab/AHaber_ARegev_Nature_2017/Ncbi_Data/Supplemental_Files/Zipped/Duplicate/GSE92332_atlas_UMIcounts/GSE92332_Atlas_UMIcounts_Frm_Methods/Batch_File/atlas_batch.txt")
dim(batch)
#head(batch)
batch = as.numeric(batch)

#batch correction - #default parametric adjustment mode
atlas_d_combat <- ComBat(dat=atlas_d, batch=batch, mod=NULL, par.prior=TRUE, prior.plots=FALSE)
typeof(atlas_d_combat)
atlas_d_combat[1:5,1:5]
#write.table(atlas_d_combat,"/Volumes/SeagateEpv/Yan_Lab/AHaber_ARegev_Nature_2017/Ncbi_Data/Supplemental_Files/Zipped/Duplicate/GSE92332_atlas_UMIcounts/atlas_d_combat.txt",sep="\t",quote=FALSE)

#https://bioconductor.org/packages/release/bioc/manuals/sva/man/sva.pdf
#https://bioconductor.org/packages/release/data/experiment/vignettes/bladderbatch/inst/doc/bladderbatch.pdf
#------------------------------------------------------------------------------------
#highly variable genes - HVG
#source("https://bioconductor.org/biocLite.R")
#biocLite("scran")
library(scran)

#find row numbers for spike-in data
which(rownames(atlas_d_combat) == "Ercc1") # 6470
which(rownames(atlas_d_combat) == "Ercc2") # 6473
which(rownames(atlas_d_combat) == "Ercc3") # 15118
which(rownames(atlas_d_combat) == "Ercc4") # 13757
which(rownames(atlas_d_combat) == "Ercc5") # 139
which(rownames(atlas_d_combat) == "Ercc6") #  9327
which(rownames(atlas_d_combat) == "Ercc8") # 12507
which(rownames(atlas_d_combat) == "Ercc6l") # 2462
which(rownames(atlas_d_combat) == "Ercc6l2") # 12281

is.spike <- logical(15791)
#is.spike <- logical(ngenes)
spike_row_list <- c(6470,6473,15118,13757,139,9327,12507,2462,12281)
is.spike[spike_row_list] <- TRUE
#validate
which(is.spike %in% "TRUE")

#run - df? - dependent on 'n'?
hvg_out <- improvedCV2(atlas_d_combat, is.spike,df=) #fdr column represents adjusted pvalues
head(hvg_out)
plot(hvg_out$mean, hvg_out$cv2, log="xy")
points(hvg_out$mean, hvg_out$trend, col="red", pch=16, cex=0.5)

#fdr<0.05 - hvg
hvg_out_2 <- hvg_out[which(hvg_out$FDR<0.05),]
dim(hvg_out_2) 
#2642    6

#rownames of hvg
hvg_rownames <- noquote(rownames(as.matrix(hvg_out_2)))
hvg_rownames <- as.matrix(hvg_rownames)
head(hvg_rownames)

#subset hvg from atlas_d_comba
#%in% is a match operator
atlas_e <- subset(atlas_d_combat, rownames(atlas_d_combat) %in% hvg_rownames)
dim(atlas_e)
#2642 7216
#write.table(atlas_e,"/Volumes/SeagateEpv/Yan_Lab/AHaber_ARegev_Nature_2017/Ncbi_Data/Supplemental_Files/Zipped/Duplicate/GSE92332_atlas_UMIcounts/atlas_e.txt",sep="\t",quote=FALSE)
#atlas_e <- read.table("/Volumes/SeagateEpv/Yan_Lab/AHaber_ARegev_Nature_2017/Ncbi_Data/Supplemental_Files/Zipped/Duplicate/GSE92332_atlas_UMIcounts/atlas_e.txt",header=T,fill=FALSE)
#------------------------------------------------------------------------------------
#PCA analysis and tSNE
#center values: zero centered
#scale them: scaled to have unit variance
atlas_e_2 <- scale(atlas_e,scale=TRUE,center=TRUE) #no need we can do this in the following step
class(atlas_e_2)

#randomization approximation to PCA
#randomized principal component analysis
#PCA using randomized singular value decomposition
#install.packages("rsvd")
library(rsvd)
#rpca function: returns a list
atlas_f <- rpca(atlas_e_2,k=100,center=FALSE,scale=FALSE)
#screeplot
ggscreeplot(atlas_f)
#individual factor map: plotting principal component scores
#ggindplot(atlas_f,pcs = c(1,2))
#correlation plot: correlation between original variable and PCs
#ggcorplot(atlas_f,pcs = c(1,2))

#summary(atlas_f)
#print(atlas_f)

#typeof(atlas_f)
#"list"
#number of objects in list
#length(atlas_f)
#list names
#names(atlas_f)
# "rotation" "eigvals"  "sdev"     "var"      "center"   "scale"    "x"
#list slicing
#atlas_f["rotation"]
#atlas_f["center"]
#atlas_f["x"]
#rotations are eigenvectors: PCs
#extract "rotation" object from list
atlas_f_pcs <- as.matrix(atlas_f[["rotation"]])
class(atlas_f_pcs)
dim(atlas_f_pcs)

#----------------------------------------------------------------
#identify significant PCs: 13
#estimate a number of significant PCs from permutation test
#install.packages("jackstraw")
#ERROR: dependency ‘lfa’ is not available for package ‘jackstraw’
#source("https://bioconductor.org/biocLite.R")
#biocLite("lfa")
#install.packages("jackstraw")
library(jackstraw)
#permutationPA function: permutation parallel analysis
permutationPA(atlas_f_pcs,B = 100, threshold = 0.05)

#'Barnes hut' implementation of tSNE
#install.packages("Rtsne")
library(Rtsne)
#subset data
atlas_f_pcs_b <- atlas_f_pcs[,c(1:13)]
dim(atlas_f_pcs_b)
#[1] 7216   13

atlas_tsne <- Rtsne(atlas_f_pcs_b,
                    dims=2,
                    perplexity=20,
                    max.iter=20000)

names(atlas_tsne)
#Y: matrix containing new representations for the objects
#head(atlas_tsne["Y"])
#20,000 iterations and a perplexity setting that ranged from 10 to 30
#https://www.codeproject.com/Tips/788739/Visualization-of-High-Dimensional-Data-using-t-SNE

x = noquote(rownames(atlas_f_pcs_b))
x <- as.matrix(x)
dim(as.matrix(x))
#write.table(x,"/Volumes/SeagateEpv/Yan_Lab/AHaber_ARegev_Nature_2017/Ncbi_Data/Supplemental_Files/Zipped/Duplicate/GSE92332_atlas_UMIcounts/GSE92332_Atlas_UMIcounts_Frm_Methods/x",sep="\t",quote=FALSE)
#cut -d'_' -f3 x > x2
x2 <- read.table("/Volumes/SeagateEpv/Yan_Lab/AHaber_ARegev_Nature_2017/Ncbi_Data/Supplemental_Files/Zipped/Duplicate/GSE92332_atlas_UMIcounts/GSE92332_Atlas_UMIcounts_Frm_Methods/x2", header=TRUE, fill=FALSE)
x2 <- as.matrix(x2)
x2 <- noquote(x2)

#plot(atlas_tsne$Y, col=factor(x2), pch=20, main="tSNE")
#text(atlas_tsne$Y, labels=factor(x2))

#cluster analysis 
#install.packages("cccd")
#install.packages("igraph")
library(cccd)
library(igraph)

#k-nearest-neighbor (KNN) graph on the data 
#euclidean distance
#k = 200
#'nng' function from the R package 'cccd' 
#the k-NN graph used as input to infomap.community function from the 'igraph' R package
w <- nng(x=atlas_e,k=200)
u <- cluster_infomap(w)
#membership(u)
#communities(u)

#-------------------------------------------------------------------------------------------




