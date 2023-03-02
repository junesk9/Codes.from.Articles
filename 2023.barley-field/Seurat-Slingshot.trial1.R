#############################################################
########### Seurat & Slingshot analysis of Barley field mRNA-seq
######################### Junesk9 2022.12.08

#Libraries
library(tidyverse)
library(stringr) #to modulate the row-names
library(patchwork)
library(Seurat)

library(glmnet)
library(slingshot)
library(clusterExperiment)

#library(gam) #for gam
library(e1071) #for SVM
library(earth) #for MARS
library(mgcv) # for GAM, smoothing splines
library(MuMIn) # for AIC calc

###############################################
#################################Prepare inputs
rpm_f = "../00_IPSR_KIBR_Sampledata/RPM_IPSR_KIBR_LB_17_19_BLD20201007.txt"
meta_f = "../00_IPSR_KIBR_Sampledata/samplesheet.txt"
hormon_f = "../00_IPSR_KIBR_Sampledata/Hormone_IPSR_KIBR_LB_17_19_BLD20201014.csv"

####### Filter the metadata
meta <- read.table(meta_f, header=T)
rownames(meta) <- meta$SampleID
meta <- meta %>% filter(!Place %in% c("BLD"))

rpm <- read.table(rpm_f, header=T, row.names = 1)
colnames(rpm) <- str_sub(colnames(rpm), 2) ## remove the first "X" from the colnames, using stringr()
rpm <- rpm[,meta$SampleID] ## filter only columns with metadata

hormone <- read.table(hormon_f, row.names=1, header=T, sep=",")
colnames(hormone) <- str_sub(colnames(hormone), 2) 
hormone <- t(hormone)

dim(rpm)
dim(meta)
#[1] 39734  1942
#[1] 1942   10

####### Construct Seurat object
ddseq <- CreateSeuratObject(counts=rpm, project="ddseq", min.cells=100, min.features=8000) #cut few samples by detected gene No.
# check the input integrity
plot(ddseq@meta.data$nCount_RNA) #  Total read count
plot(ddseq@meta.data$nFeature_RNA) # Number of the detected genes

#add metadata to ddseq
meta2 <- meta[rownames(ddseq@meta.data),]
ddseq@meta.data$Place <- meta2$Place
ddseq@meta.data$Accession <- meta2$Accession
ddseq@meta.data$Year <- meta2$Year
ddseq@meta.data$Case <- meta2$Case
ddseq@meta.data$Date <- meta2$Date
ddseq@meta.data$Date.2 <- meta2$Date.2
ddseq@meta.data$Date.3 <- meta2$Date.3
ddseq@meta.data$Date.4 <- meta2$Date.4

# Normalisation & define Highly Varying expressions
ddseq = NormalizeData(ddseq, scale.factor = 10^6) # log2RPM
ddseq = ScaleData(ddseq,do.scale = F, do.center = T) # Centering the normalized data to improve the performance of PCA, but not scaling, which might increase noise.
ddseq <- FindVariableFeatures(ddseq, selection.method = "mvp") #, dispersion.cutoff = c(0, Inf),mean.cutoff = c(0.1, 8)) #default as "vst"
VariableFeatures <- VariableFeatures(ddseq)
VariableFeaturePlot(ddseq)

ddseq[["RNA"]]
#Assay data with 25145 features for 1940 cells
#Top 10 variable features:
#  HORVU1Hr1G027140, HORVU5Hr1G043920, HORVU6Hr1G049120, HORVU2Hr1G010690, HORVU4Hr1G031730, HORVU3Hr1G056130, HORVU3Hr1G016750,
#HORVU0Hr1G011110, HORVU0Hr1G034290, HORVU2Hr1G099830 

###########################################################################3
################################################ Clustering (PCA)
ddseq <- RunPCA(ddseq, verbose = F)
ElbowPlot(ddseq, ndims=50)
PCs = ddseq@reductions$pca@cell.embeddings
#check EigenValues
sqrt(Stdev(ddseq, reduction="pca")[1:10])
#[1] 4.553067 4.166118 3.683503 3.401620 3.175950 3.012079 2.937055 2.725826 2.574140 2.496020
sum(sqrt(Stdev(ddseq, reduction="pca")))
#97.64946

DimPlot(ddseq, reduction = "pca", group.by=c("Case"))
d1.2 <- DimPlot(ddseq, reduction = "pca", group.by=c("Accession","Case"), dims=c(1,2))
d3.4 <- DimPlot(ddseq, reduction = "pca", group.by=c("Accession","Case"), dims=c(3,4))
d5.6 <- DimPlot(ddseq, reduction = "pca", group.by=c("Accession","Case"), dims=c(5,6))
d7.8 <- DimPlot(ddseq, reduction = "pca", group.by=c("Accession","Case"), dims=c(7,8))
d9.10 <- DimPlot(ddseq, reduction = "pca", group.by=c("Accession","Case"), dims=c(9,10))

################ Manual selection of PCA-Dims
##choose #1,3,6,7,8,9
sum(sqrt(Stdev(ddseq, reduction="pca")[c(1,3,6,7,8,9)]))
#[1] 19.48567 #~20% of variations explained 
## umap
ddseq <- RunUMAP(ddseq, dims=c(1,3,6,7,8,9))
DimPlot(ddseq, reduction = "umap", group.by = c("Accession","Case"))
FeaturePlot(ddseq, reduction = "umap", features = "Date.4")

############### GLR-aided selection of PCA-Dims
lasso.model.cv <- cv.glmnet(x = ddseq@reductions$pca@cell.embeddings, y = ddseq@meta.data$Date.4, family = "gaussian", alpha = 1)
plot(lasso.model.cv, xvar="lambda", label=TRUE)
ridge.model.cv <- cv.glmnet(x = ddseq@reductions$pca@cell.embeddings, y = ddseq@meta.data$Date.4, family = "gaussian", alpha = 0)
plot(ridge.model.cv, xvar="lambda", label=TRUE)

alpha <- seq(0.01, 0.99, 0.01)
mse.df <- NULL
for (i in 1:length(alpha)) {
  m <- cv.glmnet(x = ddseq@reductions$pca@cell.embeddings, y = ddseq@meta.data$Date.4, family = "gaussian", alpha = alpha[i])
  mse.df <- rbind(mse.df, data.frame(alpha = alpha[i],
                                     mse = min(m$cvm)))
}
best.alpha <- mse.df$alpha[mse.df$mse == min(mse.df$mse)]
#[1] 0.96
en.model.cv <- cv.glmnet(x = ddseq@reductions$pca@cell.embeddings, y = ddseq@meta.data$Date.4, family = "gaussian", alpha = best.alpha)
plot(en.model.cv, xvar="lambda", label=TRUE)
en.model <- glmnet(x = ddseq@reductions$pca@cell.embeddings, y = ddseq@meta.data$Date.4, family = "gaussian", lambda = en.model.cv$lambda.1se)
plot(abs(en.model$beta), main="betas of PC-dim, EN-model", pch=19)
abline(h=0.008,lty=2,col=2)
pcs <- which(abs(en.model$beta)>0.008)
#[1]  1  3 14 17 18 20 40
sum(sqrt(Stdev(ddseq, reduction="pca")[pcs]))

### replace the PC data of Seurat object
ddseq2 <- ddseq
PC_table <- ddseq@reductions$pca@cell.embeddings
for(i in 1:length(pcs)){
  ddseq2@reductions$pca@cell.embeddings[,i] = PC_table[,pcs[i]]
}
ddseq2@reductions$pca@cell.embeddings = ddseq2@reductions$pca@cell.embeddings[,1:length(pcs)]
### Check UMAP clusters
ddseq2 <- RunUMAP(ddseq2, dims=c(1:length(pcs)))
DimPlot(ddseq2, reduction = "umap", group.by = c("Accession","Case"))
FeaturePlot(ddseq2, reduction = "umap", features = "Date.4")

umap <- ddseq2@reductions$umap@cell.embeddings

###############################################################
################################### Clustering 
ddseq2 = FindNeighbors(ddseq2, dims=1:length(pcs))
#res=0.5
ddseq2 = FindClusters(ddseq2,resolution = 0.5)
#DimPlot(ddseq2, reduction = "pca")
#DimPlot(ddseq2, reduction = "umap")
#res=0.8
ddseq2 = FindClusters(ddseq2,resolution = 0.8)
#DimPlot(ddseq2, reduction = "pca")
#DimPlot(ddseq2, reduction = "umap")
#res=1
ddseq2 = FindClusters(ddseq2,resolution = 1)
#DimPlot(ddseq2, reduction = "pca")
#DimPlot(ddseq2, reduction = "umap")
#res=1.2
ddseq2 = FindClusters(ddseq2,resolution = 1.2)
#DimPlot(ddseq2, reduction = "pca")
#DimPlot(ddseq2, reduction = "umap")
#res=1.5
ddseq2 = FindClusters(ddseq2,resolution = 1.5)
#DimPlot(ddseq2, reduction = "pca")
#DimPlot(ddseq2, reduction = "umap")
#res=2
ddseq2 = FindClusters(ddseq2,resolution = 2)
#DimPlot(ddseq2, reduction = "pca")
#DimPlot(ddseq2, reduction = "umap")

ddseq@meta.data$clsGLM_res0.5 = ddseq2@meta.data$RNA_snn_res.0.5
ddseq@meta.data$clsGLM_res0.8 = ddseq2@meta.data$RNA_snn_res.0.8
ddseq@meta.data$clsGLM_res1 = ddseq2@meta.data$RNA_snn_res.1
ddseq@meta.data$clsGLM_res1.2 = ddseq2@meta.data$RNA_snn_res.1.2
ddseq@meta.data$clsGLM_res1.5 = ddseq2@meta.data$RNA_snn_res.1.5
ddseq@meta.data$clsGLM_res2 = ddseq2@meta.data$RNA_snn_res.2

data = ddseq@meta.data


################## Slingshot analysis
Idents(ddseq2) <- rep(1, length(Idents(ddseq2))) ##no category
sce0 <- as.SingleCellExperiment(ddseq2)
sce0 <- slingshot(sce0,clusterLabels = "ident",reducedDim = "PCA")
ddseq@meta.data$Slingshot_free = sce0$slingPseudotime_1

Idents(ddseq2) <- ddseq@meta.data$Date.4
sce1 <- as.SingleCellExperiment(ddseq2)
sce1 <- slingshot(sce1,clusterLabels = "ident",reducedDim = "PCA") ##taking minutes
ddseq@meta.data$Slingshot_Date4 = sce1$slingPseudotime_1

Idents(ddseq2) <- ddseq@meta.data$Date.3
sce2 <- as.SingleCellExperiment(ddseq2)
sce2 <- slingshot(sce2,clusterLabels = "ident",reducedDim = "PCA") ##taking minutes
ddseq@meta.data$Slingshot_Date3 = sce2$slingPseudotime_1

Idents(ddseq2) <- as.integer(str_sub(ddseq@meta.data$Date.2, 2))
sce3 <- as.SingleCellExperiment(ddseq2)
sce3 <- slingshot(sce3,clusterLabels = "ident",reducedDim = "PCA")
ddseq@meta.data$Slingshot_Date2 = sce3$slingPseudotime_1

data = ddseq@meta.data
################# Data output, visualization
g100 = ggplot(data, aes(x = Date.4, y = Slingshot_free, color = Accession))+
  geom_point(alpha=0.8, position = position_jitterdodge(0.5))+
  theme_bw() 
plot(g100)

g101 = ggplot(data, aes(x = Date.4, y = Slingshot_free, color = Accession))+
  geom_point(alpha=0.8, position = position_jitterdodge(0.1))+
  facet_wrap(vars(Accession), ncol = 4L)+
  theme_bw()
plot(g101)

g102 = ggplot(data, aes(x = Date.4, y = Slingshot_Date4, color = Accession))+
  geom_point(alpha=0.8, position = position_jitterdodge(0.1))+
  facet_wrap(vars(Accession), ncol = 4L)+
  theme_bw()
plot(g102)

g103 = ggplot(data, aes(x = Date.4, y = Slingshot_Date3, color = Accession))+
  geom_point(alpha=0.8, position = position_jitterdodge(0.1))+
  facet_wrap(vars(Accession), ncol = 4L)+
  theme_bw()
plot(g103)

g104 = ggplot(data, aes(x = Date.4, y = Slingshot_Date2, color = Accession))+
  geom_point(alpha=0.8, position = position_jitterdodge(0.1))+
  facet_wrap(vars(Accession), ncol = 4L)+
  theme_bw()
plot(g104)

png("das-pseudotime_LB_PC1713_VF2000_cat2_220719.png", res=300, width = 2800, height = 5000)
g100/g101/g104/g103/g102
dev.off()

#############################################################################
################################### Slingshot
ddseq2 = FindNeighbors(ddseq2, dims=1:length(pcs))
#res=0.5
ddseq2 = FindClusters(ddseq2,resolution = 0.5)
#DimPlot(ddseq2, reduction = "pca")
#DimPlot(ddseq2, reduction = "umap")
#res=0.8
ddseq2 = FindClusters(ddseq2,resolution = 0.8)
#DimPlot(ddseq2, reduction = "pca")
#DimPlot(ddseq2, reduction = "umap")
#res=1
ddseq2 = FindClusters(ddseq2,resolution = 1)
#DimPlot(ddseq2, reduction = "pca")
#DimPlot(ddseq2, reduction = "umap")
#res=1.2
ddseq2 = FindClusters(ddseq2,resolution = 1.2)
#DimPlot(ddseq2, reduction = "pca")
#DimPlot(ddseq2, reduction = "umap")
#res=1.5
ddseq2 = FindClusters(ddseq2,resolution = 1.5)
#DimPlot(ddseq2, reduction = "pca")
#DimPlot(ddseq2, reduction = "umap")
#res=2
ddseq2 = FindClusters(ddseq2,resolution = 2)
#DimPlot(ddseq2, reduction = "pca")
#DimPlot(ddseq2, reduction = "umap")

ddseq@meta.data$clsGLM_res0.5 = ddseq2@meta.data$RNA_snn_res.0.5
ddseq@meta.data$clsGLM_res0.8 = ddseq2@meta.data$RNA_snn_res.0.8
ddseq@meta.data$clsGLM_res1 = ddseq2@meta.data$RNA_snn_res.1
ddseq@meta.data$clsGLM_res1.2 = ddseq2@meta.data$RNA_snn_res.1.2
ddseq@meta.data$clsGLM_res1.5 = ddseq2@meta.data$RNA_snn_res.1.5
ddseq@meta.data$clsGLM_res2 = ddseq2@meta.data$RNA_snn_res.2

data = ddseq@meta.data
##################################################################
################### Cluster: matching to the average pseudo-time
Idents(ddseq2) = ddseq@meta.data$clsGLM_res0.5

cls = as.character(unique(ddseq2@active.ident))
pstime_med1 = NULL
for(i in cls){
  m = median(data$Slingshot_free[ddseq2@active.ident==i])
  pstime_med1 = rbind(pstime_med1, data.frame("cluster" = i, "pseudotime" = m))
}
pstime_med1 = pstime_med1[order(pstime_med1$pseudotime),]
Idents(ddseq2) = factor(ddseq2@active.ident, levels = pstime_med1$cluster)
ddseq@meta.data$ident = ddseq2@active.ident

FeaturePlot(ddseq2, reduction = "umap", features = "Date.4")
FeaturePlot(ddseq2, reduction = "umap", features = "Slingshot_free", max.cutoff = 70)

data = ddseq@meta.data
##################### Cluster: DEG analysis




###################################################################
################## other models to identify the impact PCs
#gam
d2h <- ddseq@meta.data$Date.4
PCs <- as.data.frame(PCs)
options(na.action = "na.fail")

gam.model <- gam(d2h~s(PC_1)+s(PC_2)+s(PC_3)+s(PC_4)+s(PC_5)+s(PC_6)+s(PC_7)+s(PC_8)+s(PC_9)+s(PC_10)
                 +s(PC_11)+s(PC_12)+s(PC_13)+s(PC_14)+s(PC_15)+s(PC_16)+s(PC_17)+s(PC_18)+s(PC_19)+s(PC_20)
                 +s(PC_21)+s(PC_22)+s(PC_23)+s(PC_25)+s(PC_26)+s(PC_27)+s(PC_28)+s(PC_29)+s(PC_30)+s(PC_24)
                 +s(PC_31)+s(PC_32)+s(PC_33)+s(PC_34)+s(PC_35)+s(PC_36)+s(PC_37)+s(PC_38)+s(PC_39)+s(PC_40)
                 +s(PC_41)+s(PC_42)+s(PC_43)+s(PC_44)+s(PC_45)+s(PC_46)+s(PC_47)+s(PC_48)+s(PC_49)+s(PC_50), data=PCs)
gam.anova <- as.data.frame(summary(gam.model)$anova)
gam.anova <- gam.anova[order(gam.anova[,3]),]
gam.PCs <- row.names(gam.anova)[-51]
aic.df <- NULL
for(i in c(1:length(gam.PCs))){
  f <- paste0("d2h~ ", paste(gam.PCs[1:i], collapse="+"))
  f = as.formula(f)
  gam.tmp = gam(f, data=PCs)
  aic = summary(gam.tmp)$aic
  aic.df <- rbind(aic.df, data.frame(nPC = i, AIC = aic))
}

######### Choose PCs
plot(aic.df[,1], aic.df[,2], main="AIC by PC numbers", ylab="AIC", xlab="No. PC")
## Plot indicates first NINE are effective
gam.PCs[1:9]
#[1] "s(PC_7)"  "s(PC_16)" "s(PC_18)" "s(PC_10)" "s(PC_15)" "s(PC_1)"  "s(PC_46)" "s(PC_28)" "s(PC_27)"
######## Dredge for AICs & AWs
f <- paste0("d2h~ ", paste(gam.PCs[1:9], collapse="+"))
f = as.formula(f)
gam9.model <- gam(f, data=PCs)
gam9.dredge <- dredge(gam9.model, rank="AIC")
#Model selection table 
#(Int)  s(PC_1)  s(PC_10)  s(PC_15)   s(PC_16) s(PC_18)  s(PC_27)   s(PC_28)   s(PC_46)   s(PC_7) df   logLik     AIC   delta weight
#512 0.6065 0.008687 0.0048740 -0.005099 -1.199e-03 0.009387 -0.007487 -4.897e-03 -3.427e-03 -0.006767 11  872.137 -1722.3    0.00      1
#384 0.6065 0.008717 0.0048270 -0.004619 -1.595e-03 0.009890 -0.007580 -5.476e-03            -0.006798 10  861.176 -1702.4   19.92      0
#508 0.6065 0.008729 0.0046700           -1.371e-03 0.009263 -0.007560 -5.099e-03 -2.588e-03 -0.007030 10  843.562 -1667.1   55.15      0
gam9.model.best <- get.models(gam9.dredge, subset = 1)[1]
sum(sqrt(Stdev(ddseq, reduction="pca")[c(7,16,18,10,15,1,46,28,27)]))
#[1] 20.56857

###UMAP with given PCs
select.PCs <- sub(")","",str_sub(gam.PCs[1:9], 6))
select.PCs <- as.numeric(select.PCs)
ddseq3 <- RunUMAP(ddseq,dims=select.PCs)
DimPlot(ddseq3, reduction = "umap", group.by = c("Accession","Case"))
FeaturePlot(ddseq3, reduction = "umap", features = "Date.4")


################################################################################
###################### Re-Seurat with dynamic genes
VariableFeaturePlot(ddseq)
VariableFeatures <- VariableFeatures(ddseq)
rpm_d <- rpm[VariableFeatures,]
ddseq2 <- CreateSeuratObject(counts=rpm_d, project="ddseq",min.cells=100, min.features=3000)



