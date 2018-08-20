#load 
atlas_a <- read.table("/Volumes/SeagateEpv/Yan_Lab/AHaber_ARegev_Nature_2017/Ncbi_Data/Supplemental_Files/Zipped/Duplicate/GSE92332_atlas_UMIcounts/GSE92332_atlas_UMIcounts.txt", header=TRUE, fill=FALSE)
atlas_a <- as.matrix(atlas_a)

#-------------------------------------------------------------------------------------------
#exclude cells with fewer than 800 detected genes
col_sums_a <- as.matrix(colSums(atlas_a))
dim(col_sums_a)
#-----------------------------------------------------------------
#normalize for differences in coverage
atlas_b <- apply(atlas_a,2,function(x){x/sum(x)})

#multiply by 10,000 to create TPM-like values
atlas_c <-apply(atlas_b,2,function(x){x*10000})

#log2transform
atlas_d <- apply(atlas_c,2,function(x){log2(x+1)})

##---------------------------------------------------------------
#----------------------------------------------------------------
#----------------------------------------------------------------
library(sva)
batch = read.table("/Volumes/SeagateEpv/Yan_Lab/AHaber_ARegev_Nature_2017/Ncbi_Data/Supplemental_Files/Zipped/Duplicate/GSE92332_atlas_UMIcounts/GSE92332_Atlas_UMIcounts_Frm_Methods/Batch_File/atlas_batch.txt")
batch = as.numeric(batch)

#batch correction - #default parametric adjustment mode
atlas_d_combat <- ComBat(dat=atlas_d, batch=batch, mod=NULL, par.prior=TRUE, prior.plots=FALSE)
#------------------------------------------------------------------------------------
#----------------------------------------------------------------
#----------------------------------------------------------------
#highly variable genes - HVG
library(statmod)
#calculate estimate of variance and estimate of coefficient of variation
means <- rowMeans(atlas_a)
means <- as.matrix(means)
typeof(means)
class(means)
dim(as.matrix(means))
head(as.matrix(means))

vars <- apply(atlas_a,1,var)
vars <- as.matrix(vars)
typeof(vars)
class(means)
dim(vars)
head(vars)

cv2 <- vars/means^2
cv2 <- as.matrix(cv2)
typeof(vars)
class(cv2)
dim(cv2)
head(cv2)

#exclude genes that have low means and hence high cv2 from fit
minMeanForFit <- unname( quantile( means[ which( cv2 > 100 ) ], .95 ) )
useForFit <- (means >= minMeanForFit)
#regress cv2 on 1/means using glm gamma function with log link
fit <- glmgam.fit( cbind( a0 = 1, a1tilde = 1/means[useForFit] ),cv2[useForFit] )
#noise coefficients
a0 <- unname( fit$coefficients["a0"] )
a1 <- unname( fit$coefficients["a1tilde"])
fit$coefficients

#statistical significance of deviation: genes that siginificantly deviate from the fit
df <- ncol(atlas_a) - 1
afit <- a1/means+a0
varFitRatio <- vars/(afit*means^2)
pval <- pchisq(varFitRatio*df,df=df,lower.tail=F)
adj.pval <- p.adjust(pval,"fdr")
df1 <- data.frame(rownames(atlas_a),varFitRatio,pval,adj.pval)
colnames(df1) <- c("gene","variance_fit_ratio","pval","adj.pval")
sigVariedGenes <- df1$adj.pval<5e-2; 
table(sigVariedGenes)

#subset these genes from atlas_a
siggenes <- df1[df1$adj.pval<5e-2,]
dim(as.matrix(siggenes))
siggenes_rownames <- noquote(rownames(as.matrix(siggenes)))

#subset based on rownames
atlas_e <- subset(atlas_d_combat, rownames(atlas_d_combat) %in% siggenes_rownames)
dim(atlas_e)

atlas_a_hvg <- subset(atlas_a,rownames(atlas_a) %in% siggenes_rownames)
#------------------------------------------------------------------------------------
#----------------------------------------------------------------
#----------------------------------------------------------------
#PCA analysis and tSNE
library(rsvd)
#rpca function: returns a list
atlas_f <- rpca(t(atlas_e),k=100,center=TRUE,scale=TRUE,retx=T)
atlas_f_pcs <- atlas_f$x


#identify significant PCs: 13
#library(jackstraw)
#y = sig.pcs.perm(dat=t(atlas_umis[var.genes,]), center=T, scale=T, max.pc=100, B=1000, n.cores=20,randomized=T)
#y <- permutationPA(t(atlas_e),B = 1000,n.cores=4,center=TRUE,scale=TRUE,max.pc=100,randomized=T)
atlas_f_pcs_b <- atlas_f_pcs[,c(1:13)]

#'Barnes hut' implementation of tSNE 
#install.packages("Rtsne")
library(Rtsne)
atlas_tsne <- Rtsne(atlas_f_pcs_b,
                    check_duplicates=TRUE,	
                    pca=FALSE,
                    perplexity=20,
                    initial_dims = 13,
                    max.iter=100000,
                    verbose=TRUE,
                    whiten=FALSE)

#----------------------------------------------------------------
#----------------------------------------------------------------
#----------------------------------------------------------------
#cluster analysis 
library(cccd)
library(igraph)

dist <- as.matrix(dist(atlas_f_pcs_b))
graph <- nng(x=atlas_f_pcs_b,k=200)
clustering <- cluster_infomap(graph,modularity=TRUE)
clusters <- clustering$membership

#plot tsne
brewer16 = c(brewer.pal(9, "Set1"), brewer.pal(7, "Set2"))
brewer16[6] = "khaki2"
brewer16[8] = "lightskyblue2"

display.brewer.pal(9, "Set1")
display.brewer.pal(7,"Set2")

library(ggplot2)
x = data.frame(atlas_tsne$Y, clusters)
ggplot(x, aes(x=X1, y=X2, color=factor(clusters))) + geom_point() + scale_color_manual(values=brewer16)

#----------------------------------------------------------------
#----------------------------------------------------------------
#----------------------------------------------------------------
#BINOMIAL TEST:from regev paper
atlas_a_hvg <- subset(atlas_a,rownames(atlas_c) %in% siggenes_rownames)
atlas_c_hvg <- subset(atlas_c,rownames(atlas_c) %in% siggenes_rownames)
#test for enrichments in cluser #1
#extract cluster 1 cell names
#prep data
w <- data.frame(colnames(atlas_c_hvg),clusters)
colnames(w) <- c("colnames_atlas_c_hvg","cluster")
cluster_1_cellnames <- w[w$cluster==1,"colnames_atlas_c_hvg"]
except_cluster_1_cellnames <- w[w$cluster!=1,"colnames_atlas_c_hvg"]

effect.size <- log(2)
m = apply(atlas_c_hvg[,except_cluster_1_cellnames],1,function(x) sum(x>0))
m1 =m; m1[m==0]=1;
n = apply(atlas_c_hvg[,cluster_1_cellnames],1,function(x) sum(x>0))
pv1 = pbinom(n, length(cluster_1_cellnames), m1/length(except_cluster_1_cellnames), lower.tail = FALSE) + dbinom(n, length(cluster_1_cellnames), m1/length(except_cluster_1_cellnames))

log_fold_express = log(n*length(except_cluster_1_cellnames)/(m*length(cluster_1_cellnames))) #log proportion of expressing cells
d1 <- data.frame(log.effect=log_fold_express,pval=pv1)
d1 <- subset(d1, log.effect >= effect.size)
d1 <- d1[order(d1$pval,decreasing=FALSE),]
#d1 <- d1[order(d1$pval,decreasing=FALSE),]

n1 = n; n1[n==0]=1;
pv2 = pbinom(m, length(except_cluster_1_cellnames), n1/length(cluster_1_cellnames), lower.tail=FALSE) + dbinom(m, length(except_cluster_1_cellnames), n1/length(cluster_1_cellnames))
d2 <- data.frame(log.effect=log_fold_express,pval=pv2)
d2 <- subset(d2, log.effect <= -effect.size)
d2 <- d2[order(d2$pval,decreasing=FALSE),]

d = rbind(d1, d2);
d = d[order(d$pval, decreasing=FALSE),]

posFrac.1 = apply(atlas_a[rownames(d),cluster_1_cellnames],1,function(x) round(sum(x > 0)/length(x),2))
posFrac.2 = apply(atlas_a[rownames(d),except_cluster_1_cellnames],1,function(x) round(sum(x > 0)/length(x),2))

genes.include = posFrac.1 >= 0.1

cluster_1_result = d[genes.include,]
cluster_1_result = d[order(abs(d$log.effect), decreasing=TRUE),]
#write.table(cluster_1_result,"/Users/lk/Documents/cluster_1_result.txt")
#---------------------------------------------------------------------------------------------------------------




