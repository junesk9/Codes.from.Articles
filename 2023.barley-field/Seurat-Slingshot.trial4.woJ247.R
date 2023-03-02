#############################################################
########### Seurat & Slingshot analysis of Barley field mRNA-seq
######################### Junesk9 2023.02.01
# For reasons,re-try the barley field mRNA-seq analysis without "J247" data.

#Libraries
library(tidyverse)
library(stringr) #to modulate the row-names
library(dplyr)
library(patchwork)
library(Seurat)
library(RColorBrewer)

library(glmnet)
library(slingshot)
library(tradeSeq)
library(pheatmap)
library(MAST)
#library(clusterExperiment)
#library(UpSetR)

library(GSEABase) ##for GSEA
library(GOstats)
library(wordcloud)

###############################################
#################################Prepare inputs
set.seed(101)
rpm_f = "../00_IPSR_KIBR_Sampledata/RPM_IPSR_KIBR_LB_17_19_BLD20201007.txt"
meta_f = "../00_IPSR_KIBR_Sampledata/samplesheet.txt"
hormon_f = "../00_IPSR_KIBR_Sampledata/Hormone_IPSR_KIBR_LB_17_19_BLD20201014.csv"

####### Filter the metadata
meta <- read.table(meta_f, header=T)
rownames(meta) <- meta$SampleID
meta <- meta %>% filter(!Accession %in% c("J247"))

rpm <- read.table(rpm_f, header=T, row.names = 1)
colnames(rpm) <- str_sub(colnames(rpm), 2) ## remove the first "X" from the colnames, using stringr()
rpm <- rpm[,meta$SampleID] ## filter only columns with metadata

dim(meta)
dim(rpm)
#[1] 1481   10
#[1] 39734  1481

####### Construct Seurat object
ddseq <- CreateSeuratObject(counts=rpm, project="ddseq", min.cells=100, min.features=8000) #cut few samples by detected gene No.
# check the input integrity
plot(ddseq@meta.data$nFeature_RNA, main="w/o J247; whole nFeature") #  Total read count
plot(ddseq@meta.data$nCount_RNA, main="w/o J247; whole nCount") # Number of the detected genes

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

f_split <- function(s) strsplit(s,"/")[[1]][2] #extract 12 from 2018/12/05
ddseq@meta.data$month <- sapply(meta2$Date, f_split)

# Normalisation & define Highly Varying expressions
ddseq = NormalizeData(ddseq, scale.factor = 10^6) # log2RPM
ddseq = ScaleData(ddseq,do.scale = F, do.center = T) # Centering the normalized data to improve the performance of PCA, but not scaling, which might increase noise.
ddseq <- FindVariableFeatures(ddseq, selection.method = "mvp") #, dispersion.cutoff = c(0, Inf),mean.cutoff = c(0.1, 8)) #default as "vst"
VariableFeatures <- VariableFeatures(ddseq)
VariableFeaturePlot(ddseq)

ddseq[["RNA"]]
#Assay data with 24220 features for 1480 cells
#Top 10 variable features:
#  HORVU1Hr1G027140, HORVU5Hr1G043920, HORVU2Hr1G010690, HORVU3Hr1G056130, HORVU4Hr1G031730, HORVU3Hr1G016750,

###########################################################################3
################################################ Clustering (PCA)
ddseq <- RunPCA(ddseq, verbose = F)
ElbowPlot(ddseq, ndims=50)
PCs = ddseq@reductions$pca@cell.embeddings
sum(sqrt(Stdev(ddseq, reduction="pca")))
#[1] 95.91991

############### GLR-aided selection of PCA-Dims
alpha <- seq(0.01, 0.99, 0.01)
mse.df <- NULL
for (i in 1:length(alpha)) {
  m <- cv.glmnet(x = ddseq@reductions$pca@cell.embeddings, y = ddseq@meta.data$Date.4, family = "gaussian", alpha = alpha[i])
  mse.df <- rbind(mse.df, data.frame(alpha = alpha[i],
                                     mse = min(m$cvm)))
}
best.alpha <- mse.df$alpha[mse.df$mse == min(mse.df$mse)]
#[1] 0.57
en.model.cv <- cv.glmnet(x = ddseq@reductions$pca@cell.embeddings, y = ddseq@meta.data$Date.4, family = "gaussian", alpha = best.alpha)
plot(en.model.cv, xvar="lambda", label=TRUE)
en.model <- glmnet(x = ddseq@reductions$pca@cell.embeddings, y = ddseq@meta.data$Date.4, family = "gaussian", lambda = en.model.cv$lambda.1se)
plot(abs(en.model$beta), main="betas of PC-dim, EN-model", pch=19)
abline(h=0.008,lty=2,col=2)
pcs <- which(abs(en.model$beta)>0.008)
#[1]  1  4 13 15 23 38 39
sum(sqrt(Stdev(ddseq, reduction="pca")[pcs]))
#[1] 16.68288

### replace the PC data of Seurat object
ddseq2 <- ddseq
PC_table <- ddseq@reductions$pca@cell.embeddings
for(i in 1:length(pcs)){
  ddseq2@reductions$pca@cell.embeddings[,i] = PC_table[,pcs[i]]
}
ddseq2@reductions$pca@cell.embeddings = ddseq2@reductions$pca@cell.embeddings[,1:length(pcs)]
### Check UMAP clusters
ddseq2 <- RunUMAP(ddseq2, dims=pcs)
DimPlot(ddseq2, reduction = "umap", group.by = c("Accession","Case","month"))
FeaturePlot(ddseq2, reduction = "umap", features = "Date.4")

umap <- ddseq2@reductions$umap@cell.embeddings
ddseq@meta.data$umap1 <- as.data.frame(umap)$UMAP_1
ddseq@meta.data$umap2 <- as.data.frame(umap)$UMAP_2

############################# Cluster finding
#pcs <- c(1,4,13,15,23,7,8)
ddseq2 = FindNeighbors(ddseq2, dims=pcs)

ddseq2 = FindClusters(ddseq2,resolution = 1)
ddseq@meta.data$clsGLM_res1 = ddseq2@meta.data$RNA_snn_res.1
DimPlot(ddseq2, reduction = "umap", label=TRUE, label.box=TRUE, 
        label.col="black", label.size=6)

data = ddseq@meta.data
write.csv(data, "woJ247-metadata.csv")
clus = ddseq2@meta.data$RNA_snn_res.1
Idents(ddseq2) <- ddseq2@meta.data$RNA_snn_res.1

#BoxPlot of Date.4 by clusters
oind <- order(as.numeric(by(data$Date.4, data$clsGLM_res1, median)))
data$clsGLM_res1 <- ordered(data$clsGLM_res1, levels=levels(data$clsGLM_res1)[oind])
boxplot(Date.4 ~ clsGLM_res1, data=data, main="Frac. D2H by clusters")
#Stacking barplot of Accs composition of clusters.
p <- ggplot(data, aes(x = clsGLM_res1, fill = Accession))
p + geom_bar() + theme_bw()

############### Slingshot analysis
pal <- c(RColorBrewer::brewer.pal(9, "Set1"), RColorBrewer::brewer.pal(8, "Set2"))

##re-plot for the same color strategy to slingshot-plots
plot(umap, col = pal[clus],  cex=.5,pch = 16)
for(i in levels(clus)){ 
  text( mean(umap[clus==i,1]),
        mean(umap[clus==i,2]), labels = i,font = 2, cex=2) }

#slingshot lineages
lineages <- getLineages(data=umap, clusterLabels=clus, start.clus=6, end.clus=c(5,10))
SlingshotDataSet(lineages)
plot(umap, col=pal[clus], pch=16)
lines(SlingshotDataSet(lineages), lwd = 3, col = 'black', show.constraints = TRUE)

#defining the principal curves (smooth trajectories)
curves <- getCurves(lineages, thresh = 0.01, stretch = 2, allow.breaks = TRUE, shrink = 0.99)
SlingshotDataSet(curves)
#lineages: 4 
#Lineage1: 6  4  2  9  1  10  
#Lineage2: 6  4  2  9  1  7  
#Lineage3: 6  4  0  11  8  5  
#Lineage4: 6  4  0  3  

#curves: 4 
#Curve1: Length: 14.43	Samples: 685.5
#Curve2: Length: 15.961	Samples: 721.26
#Curve3: Length: 28.791	Samples: 602.77
#Curve4: Length: 15.172	Samples: 483.83
plot(umap, col=pal[clus], pch=16)
lines(SlingshotDataSet(curves), lwd = 3, col = 'black')

pseudotime <- slingPseudotime(curves)[rownames(data),]
#[1] 1480    4
data <- cbind(data, pseudotime)
ddseq2@meta.data <- data[rownames(ddseq2@meta.data),]
FeaturePlot(ddseq2, reduction = "umap", features = "Lineage1")
FeaturePlot(ddseq2, reduction = "umap", features = "Lineage3")

write.csv(data, "woJ247-metadata.csv")
################################################################
######### Accession view ..
g102 = ggplot(data, aes(x = Date.4, y = Lineage2, color = month))+
  geom_point(alpha=0.8, position = position_jitterdodge(0.1))+
  facet_wrap(vars(Accession), ncol = 4L)+
  theme_bw()
plot(g102)

g104 = ggplot(data, aes(x = Date.4, y = Lineage4, color = month))+
  geom_point(alpha=0.8, position = position_jitterdodge(0.1))+
  facet_wrap(vars(Accession), ncol = 4L)+
  theme_bw()
plot(g104)

png("das-pseudotime_LB_PC1713_VF2000_cat2_220719.png", res=300, width = 2800, height = 5000)
g104/g102
dev.off()





