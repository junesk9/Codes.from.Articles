#############################################################
########### Seurat & Slingshot analysis of Barley field mRNA-seq
######################### Junesk9 2024.01.05

library(Seurat)
library(glmnet)
library(slingshot)
library(tradeSeq)
library(pheatmap)
library(MAST)
library(CytoTRACE)
library(phateR) ## for PHATE Dim-reduction

##For figure manipulation
library(ggplot2)

#For GESA analysis
library(GSEABase) 
library(GOstats) #GO/KEGG
library(KEGGREST) #KEGG-API


##Accessaries
library(dplyr)
library(beeswarm)
library(psych) #for corr pairs()


setwd("/Users/junesk9/理化学研究所　セルロース生産研究チーム Dropbox/June-Sik Kim/オオムギ圃場mRNA.ChIPseq")
set.seed(101) #prevent random outputs


########### not needed for now
#ensure the using python version
#require(reticulate)
#conda_create("myenv", python_version="3.11")
#use_condaenv("myenv")
#py_env()
#python:         /Users/junesk9/Library/r-miniconda/envs/myenv/bin/python
#libpython:      /Users/junesk9/Library/r-miniconda/envs/myenv/lib/libpython3.11.dylib
#pythonhome:     /Users/junesk9/Library/r-miniconda/envs/myenv:/Users/junesk9/Library/r-miniconda/envs/myenv
#version:        3.11.7 | packaged by conda-forge | (main, Dec 23 2023, 14:38:52) [Clang 16.0.6 ]
#numpy:           [NOT FOUND]


#emulate the ggplot hue() color palette
ggcolor <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
require(scales)
show_col(ggcolor(8)) #test ggcolor
key.palette <- hcl.colors(13, "Cividis")
key.palette <- rev(hcl.colors(100, "RdYlBu"))
show_col(key.palette)
########################################
####################### Prepare datasets
########################################
so <- readRDS("seurat.barley-field.leaf.v2.rds")
#An object of class Seurat 
#25145 features across 1940 samples within 1 assay 
so$month <- factor(so$month, levels=c(12,1,2,3,4,5))
so$Accession <- factor(so$Accession, levels=c("J247", "J064", "J647", "H602"))
so@meta.data$week <- factor(as.integer(so@meta.data$week), levels=c(1:22))

#Load the basal tables
meta <- so@meta.data #[1] 1940   23
all.rpm <- as.data.frame(so@assays$RNA@data) #[1] 25145  1940
vari.rpm <- all.rpm[VariableFeatures(so), ] #[1] 2963 1940

flor <- read.table("00_IPSR_KIBR_Sampledata/FLOR_genes.txt", header=T, row.names=1, sep="\t")
flor <- rownames(flor)
flor.rpm <- all.rpm[rownames(all.rpm) %in% flor, ] #[1]  128 1940

#lineage curve (for GAM fitting)
curve <- readRDS("seurat-curves.barley-field.leaf.v2.rds")
#class: PseudotimeOrdering 
#cellnames(1940): 10001_IH6_171215 10002_IH6_171215 ... 40485_KJ6_190426 40486_KJ6_190426
#cellData names(2): reducedDim clusterLabels
#pathnames(1): Lineage1

#CytoTRACE
results <- CytoTRACE(as.matrix(so@assays$RNA@counts), ncores = 6, subsamplesize = 1000)
so$CytoTRACE <- 1 - results$CytoTRACE

#GOSlim/KEGG DB
go_f ="./00_IPSR_KIBR_Sampledata/Barley.IBSCv2.ensembl49.GOSlim.txt"
go <- read.table(go_f, header=T, sep="\t")
go$evi = "IEA"
go <- go[,c(2,3,1)]
goFrame <- GOFrame(go, organism="barley")
goAllFrame <- GOAllFrame(goFrame)
gsc <- GeneSetCollection(goAllFrame, setType = GOCollection())
all <- unique(as.vector(go[,3]))  ## a vector of all referenced list
#[1] 27595

kegg_f <- "./00_IPSR_KIBR_Sampledata/Barley.IBSCv2.ensembl49.KEGG.txt"
kegg.df <- read.table(kegg_f, header=T, row.names=NULL, sep="\t", colClasses=c(rep("factor",3)))
kegg.df <- unique(kegg.df[,c(2,1)]) #ほけん
head(kegg.df)
#  KEGG.Pathway   Gene.stable.ID
#1        05235 HORVU0Hr1G004120
#2        05235 HORVU7Hr1G117420
#3        05235 HORVU0Hr1G001600
keggframe <- KEGGFrame(kegg.df, organism="barely")
gsck <- GeneSetCollection(keggframe, setType=KEGGCollection())
kegg.all <- unique(as.vector(kegg.df[,2]))
#[1] 2556
kegg.annt <- keggList("pathway")
kegg.annt <- data.frame(KEGGid=names(kegg.annt), desc=as.vector(kegg.annt), row.names=NULL)
#     KEGGid                                         desc
#1 map01100                           Metabolic pathways
#2 map01110        Biosynthesis of secondary metabolites
#3 map01120 Microbial metabolism in diverse environments


runKEGG1 <- function(gid_ls, prefix){
  kegg.p <- GSEAKEGGHyperGParams(
    name = "Paramaters",
    geneSetCollection = gsck,
    geneIds = gid_ls, 
    universeGeneIds = kegg.all,
    pvalueCutoff = 1,
    testDirection = "over"
  )
  kegg.res <- hyperGTest(kegg.p)
  kegg.res_gid <- geneIdsByCategory(kegg.res)
  kegg.res2 <- summary(kegg.res)
  kegg.res2$gid <- sapply(kegg.res2[,1], function(x)paste(unlist(kegg.res_gid[x]), collapse="/"))
  kegg.res2$KEGGID <- sapply(kegg.res2$KEGGID, function(x) paste0("taes",x))
  kegg.res2$Term <- kegg.annt[kegg.res2$KEGGID,]$desc
  
  write.csv(kegg.res2, paste0(prefix,".KEGG.csv"), row.names=F)
} 
#############################################
########## Correlation of the metadata
############################################
colnames(so@meta.data)
# [1] "orig.ident"       "nCount_RNA"       "nFeature_RNA"     "CytoTRACE"        "Place"           
#[6] "Accession"        "Season"           "Case"             "Date"             "Date.2"          
#[11] "Date.3"           "fD2H"             "month"            "week"             "RNA_snn_res.0.45"
#[16] "seurat_clusters"  "Lineage1"         "Lineage2"         "Lineage0"         "LineageM"        
#[21] "AvgTemp"          "LightHour"        "DayHour"         

meta <- so@meta.data
pairs.panels(meta[,c(4,12,20,21,22,23)])
pairs.panels(meta[,c(12,21,22,23)])
pairs.panels(meta15[,c(12,21,22,23)])
pairs.panels(meta[,c(4,19,12,20,21,23)], main="All")
pairs.panels(meta[meta$Accession %in% "J247",c(4,19,12,20,21,23)], main="J247")
pairs.panels(meta[meta$Accession %in% "J064",c(4,19,12,20,21,23)], main="J064")
pairs.panels(meta[meta$Accession %in% "H602",c(4,19,12,20,21,23)], main="H602")
pairs.panels(meta[meta$Accession %in% "J647",c(4,19,12,20,21,23)], main="J647")
pairs.panels(meta[meta$seurat_clusters %in% c("4","5"),c(4,19,12,20,21,23)], main="Clus4.5")

#############################################
####################### Dimensional reduction
############################################
#I. all features
AllFeatures <- rownames(so@assays$RNA@data)
#[1] 25145
so <- RunPCA(so, vervose=FALSE, approx=FALSE, npcs=50, features=AllFeatures)
ElbowPlot(so, ndims=50, reduction = "pca") + ggtitle("Whole25145, pca50")
so@reductions$pca_whole <- so@reductions$pca #save the feature
#(UMAP)
so <- RunUMAP(so, dims=c(1:50))
DimPlot(so, reduction="umap", group.by="Accession") + ggtitle("gene:whole-pca:whole-Accession")
DimPlot(so, reduction="umap", group.by="Case") + ggtitle("gene:whole-pca:whole-Case")
FeaturePlot(so, reduction="umap", features="Date.4") + ggtitle("gene:whole-pca:whole-fD2H")
FeaturePlot(so, reduction="umap", features="CytoTRACE") + ggtitle("gene:whole-pca:whole-CytoTRACE")
so@reductions$umap_whole50 <- so@reductions$umap #save the feature

#II. mvp-Variable genes only
so <- FindVariableFeatures(so, selection.method = "mvp") #, dispersion.cutoff = c(0, Inf),mean.cutoff = c(0.1, 8)) #default as "vst"
VariableFeaturePlot(so) + ggtitle("mvp-variable featreus")
VariableFeatures <- VariableFeatures(so)
#[1] 2842
so <- RunPCA(so, vervose=FALSE, approx=FALSE, npcs=50, features=VariableFeatures)
ElbowPlot(so, ndims=50, reduction = "pca") + ggtitle("Dynamic2842, pca50")
so@reductions$pca_vari <- so@reductions$pca  #save the feature
#(UMAP)
so <- RunUMAP(so, dims=c(1:50))
DimPlot(so, reduction="umap", group.by="Accession") + ggtitle("gene:vari2842-pca:whole-Accession")
DimPlot(so, reduction="umap", group.by="Case") + ggtitle("gene:vari2842-pca:whole-Case")
FeaturePlot(so, reduction="umap", features="Date.4") + ggtitle("gene:vari2842-pca:whole-fD2H")
FeaturePlot(so, reduction="umap", features="CytoTRACE") + ggtitle("gene:vari2842-pca:whole-CytoTRACE")
so@reductions$umap_vari50 <- so@reductions$umap #save the feature

#III. FLOR-guided reduction
so <- RunPCA(so, vervose=FALSE, approx=FALSE, npcs=50, features=flor)
ElbowPlot(so, ndims=50, reduction = "pca") + ggtitle("Flor149, pca50")
so@reductions$pca_flor <- so@reductions$pca  #save the feature
#(UMAP)
so <- RunUMAP(so, dims=c(1:50))
DimPlot(so, reduction="umap", group.by="Accession") + ggtitle("gene:flor149-pca:whole-Accession")
DimPlot(so, reduction="umap", group.by="Case") + ggtitle("gene:flor149-pca:whole-Case")
FeaturePlot(so, reduction="umap", features="Date.4") + ggtitle("gene:flor149-pca:whole-fD2H")
FeaturePlot(so, reduction="umap", features="CytoTRACE") + ggtitle("gene:flor149-pca:whole-CytoTRACE")
so@reductions$umap_flor50 <- so@reductions$umap #save the feature

#IV. FLOR+Dynamic genes-guided reduction
vf2 <- unique(c(VariableFeatures, flor))
#2984
so <- RunPCA(so, vervose=FALSE, approx=FALSE, npcs=50, features=vf2)
ElbowPlot(so, ndims=50, reduction = "pca") + ggtitle("Vari+Flor2984, pca50")
so@reductions$pca_vari2 <- so@reductions$pca  #save the feature
#(UMAP)
so <- RunUMAP(so, dims=c(1:50))
DimPlot(so, reduction="umap", group.by="Accession") + ggtitle("gene:Vari+Flor2984-pca:whole-Accession")
DimPlot(so, reduction="umap", group.by="Case") + ggtitle("gene:Vari+Flor2984-pca:whole-Case")
FeaturePlot(so, reduction="umap", features="Date.4") + ggtitle("gene:Vari+Flor2984-pca:whole-fD2H")
FeaturePlot(so, reduction="umap", features="CytoTRACE") + ggtitle("gene:Vari+Flor2984-pca:whole-CytoTRACE")
so@reductions$umap_vari2_50 <- so@reductions$umap #save the feature

#VI. Check Eigenvalues of PCA reductions; currently seurat command Stdev() not working properly,
sum(sqrt(sapply(as.data.frame(so@reductions$pca_whole@cell.embeddings), sd)))
#[1] 127.2552
sum(sqrt(sapply(as.data.frame(so@reductions$pca_vari@cell.embeddings), sd)))
#[1] 97.64946
sum(sqrt(sapply(as.data.frame(so@reductions$pca_flor@cell.embeddings), sd)))
#[1] 40.93352
sum(sqrt(sapply(as.data.frame(so@reductions$pca_vari2@cell.embeddings), sd))) # equal to sum(sqrt(Stdev(so, reduction="pca")[1:50]))
#[1] 98.01709

#VII. GLR-aided PC selection
#(Vari+Flor2984)
lasso.model.cv <- cv.glmnet(x = so@reductions$pca@cell.embeddings, y = so@meta.data$Date.4, family = "gaussian", alpha = 1)
plot(lasso.model.cv, xvar="lambda", label=TRUE, main="gene:Vari+Flor2984-PCA50-lasso")
ridge.model.cv <- cv.glmnet(x = so@reductions$pca@cell.embeddings, y = so@meta.data$Date.4, family = "gaussian", alpha = 0)
plot(ridge.model.cv, xvar="lambda", label=TRUE, main="gene:Vari+Flor2984-PCA50-ridge")
#(EN-model)
alpha <- seq(0.01, 0.99, 0.01)
mse.df <- NULL
for (i in 1:length(alpha)) {
  m <- cv.glmnet(x = so@reductions$pca@cell.embeddings, y = so@meta.data$Date.4, family = "gaussian", alpha = alpha[i])
  mse.df <- rbind(mse.df, data.frame(alpha = alpha[i], mse = min(m$cvm)))
}
best.alpha <- mse.df$alpha[mse.df$mse == min(mse.df$mse)]
#[1] 0.62
en.model.cv <- cv.glmnet(x = so@reductions$pca@cell.embeddings, y = so@meta.data$Date.4, family = "gaussian", alpha = best.alpha)
plot(en.model.cv, xvar="lambda", label=TRUE, main="gene:Vari+Flor2984-PCA50-EN")
en.model <- glmnet(x = so@reductions$pca@cell.embeddings, y = so@meta.data$Date.4, family = "gaussian", lambda = en.model.cv$lambda.1se)
plot(abs(en.model$beta), main="betas of PC-dim, EN-model", pch=19)
abline(h=0.006,lty=2,col=2)
PCs <- which(abs(en.model$beta)>0.006)
sum(sqrt(Stdev(so, reduction="pca")[PCs]))
#[1]  1  3 14 17 18 20 40
#17.61876

#(Re-UMAP)
so <- RunUMAP(so, dims=PCs)
DimPlot(so, reduction="umap", group.by="Accession") + ggtitle("gene:Vari+Flor2984-pca:glr7-Accession")
DimPlot(so, reduction="umap", group.by="Case") + ggtitle("gene:Vari+Flor2984-pca:glr7-Case")
FeaturePlot(so, reduction="umap", features="Date.4") + ggtitle("gene:Vari+Flor2984-pca:glr7-fD2H")
FeaturePlot(so, reduction="umap", features="CytoTRACE") + ggtitle("gene:Vari+Flor2984-pca:glr7-CytoTRACE")
so@reductions$umap_vari2_glr7 <- so@reductions$umap #save the feature


#(Whole)
lasso.model.cv <- cv.glmnet(x = so@reductions$pca@cell.embeddings, y = so@meta.data$Date.4, family = "gaussian", alpha = 1)
plot(lasso.model.cv, xvar="lambda", label=TRUE, main="gene:whole25145-PCA50-lasso")
ridge.model.cv <- cv.glmnet(x = so@reductions$pca@cell.embeddings, y = so@meta.data$Date.4, family = "gaussian", alpha = 0)
plot(ridge.model.cv, xvar="lambda", label=TRUE, main="gene:whole25145-PCA50-ridge")
#(EN-model)
alpha <- seq(0.01, 0.99, 0.01)
mse.df <- NULL
for (i in 1:length(alpha)) {
  m <- cv.glmnet(x = so@reductions$pca@cell.embeddings, y = so@meta.data$Date.4, family = "gaussian", alpha = alpha[i])
  mse.df <- rbind(mse.df, data.frame(alpha = alpha[i], mse = min(m$cvm)))
}
best.alpha <- mse.df$alpha[mse.df$mse == min(mse.df$mse)]
#[1] 0.64
en.model.cv <- cv.glmnet(x = so@reductions$pca@cell.embeddings, y = so@meta.data$Date.4, family = "gaussian", alpha = best.alpha)
plot(en.model.cv, xvar="lambda", label=TRUE, main="gene:whole25145-PCA50-EN")
en.model <- glmnet(x = so@reductions$pca@cell.embeddings, y = so@meta.data$Date.4, family = "gaussian", lambda = en.model.cv$lambda.1se)
plot(abs(en.model$beta), main="betas of PC-dim, EN-model", pch=19)
PCs <- which(abs(en.model$beta)>0.0045)
sum(sqrt(Stdev(so, reduction="pca")[PCs]))
#[1] 18.26949
#[1]  2 14 17 21 31 33 45
#(Re-UMAP)
so <- RunUMAP(so, reduction="pca", dims=PCs, umap.method = "umap-learn", metric = "correlation", n.components = 2)
DimPlot(so, reduction="umap", group.by="Accession") + ggtitle("gene:whole25145-pca:glr7-Accession")
DimPlot(so, reduction="umap", group.by="Case") + ggtitle("gene:whole25145-pca:glr7-Case")
FeaturePlot(so, reduction="umap", features="Date.4") + ggtitle("gene:whole25145-pca:glr7-fD2H")
FeaturePlot(so, reduction="umap", features="CytoTRACE") + ggtitle("gene:whole25145-pca:glr7-CytoTRACE")
so@reductions$umap_whole_glr7 <- so@reductions$umap #save the feature

##### Make decision
so@reductions$pca <- so@reductions$pca_vari2
PCs <- c(1, 3, 14, 17, 18, 20, 40)
so <- RunUMAP(so, reduction="pca", dims=PCs, n.components = 3) #, umap.method = "umap-learn", metric = "correlation")
so@reductions$umap_3D <- so@reductions$umap
so <- RunUMAP(so, reduction="pca", dims=PCs, n.components = 2) #, umap.method = "umap-learn", metric = "correlation")
so@reductions$umap_2D <- so@reductions$umap

DimPlot(so, reduction="umap_2D", group.by="Accession") + ggtitle("UMAP_2D-Accession")
DimPlot(so, reduction="umap_2D", group.by="Case") + ggtitle("MAP_2D-Case")
FeaturePlot(so, reduction="umap_2D", features="Date.4") + ggtitle("MAP_2D-fD2H")
FeaturePlot(so, reduction="umap_2D", features="CytoTRACE") + ggtitle("MAP_2D-CytoTRACE")

DimPlot(so, reduction="umap_3D", group.by="Accession") + ggtitle("UMAP_3D-Accession")
DimPlot(so, reduction="umap_3D", group.by="Case") + ggtitle("MAP_3D-Case")
FeaturePlot(so, reduction="umap_3D", features="Date.4") + ggtitle("MAP_3D-fD2H")
FeaturePlot(so, reduction="umap_3D", features="CytoTRACE") + ggtitle("MAP_3D-CytoTRACE")
so@reductions$umap <- so@reductions$umap_2D

#############################################
####################### Clustering - Lineage analysis
############################################
#I. set clusters... 
so <- FindNeighbors(so, dims=PCs)
so <- FindClusters(so, resolution=0.45)
so$seurat_clusters <- so$RNA_snn_res.0.45

#II. rename clusters by median fD2H
boxplot(so$Date.4~so$seurat_clusters, main="fD2H ~ cluster") #check the initial view
oind <- order(as.numeric(by(so$Date.4, so$seurat_clusters, mean)))
so$seurat_clusters <- ordered(so$seurat_clusters, levels=levels(so$seurat_clusters)[oind])
levels(so$seurat_clusters) <- c(1:7) # re-name the clusters from "[1] "2" "5" "1" "4" "6" "3" "0"" to 1:7

boxplot(so$Date.4~so$seurat_clusters, main="fD2H ~ cluster: re-named") #re-named
boxplot(so$CytoTRACE~so$seurat_clusters, main="CytoTRACE ~ cluster: re-named")
Idents(so) <- so$seurat_clusters
DimPlot(so, reduction="umap") + ggtitle("Cluster7 by res0.45: re-named")
#(option: accession by cluster)
ggplot(so@meta.data, aes(x = seurat_clusters, fill = Accession)) + geom_bar() + theme_bw()
#(option: month by cluster)
ggplot(so@meta.data, aes(y = Date.4, x = seurat_clusters, color = month))+ geom_point(alpha=0.8)+ theme_bw()

#III. Lineage analysis
clus <- so$seurat_clusters
umap <- so@reductions$umap@cell.embeddings
colnames(umap) <- c("UMAP_1","UMAP_2")

lineages <- getLineages(data=umap, clusterLabels=clus, start.clus="1")#, end.clus=c("13","7"))
plot(umap, col=ggcolor(7)[clus],  cex=.5, pch = 16)
lines(SlingshotDataSet(lineages), lwd = 3, col = 'black', show.constraints = TRUE)

curves <- getCurves(lineages, thresh = 0.01, stretch = 2, allow.breaks = TRUE, shrink = 0.99)
SlingshotDataSet(curves)
#class: SlingshotDataSet 
#Samples Dimensions
#1940          2
#lineages: 2 
#Lineage1: 1  2  3  4  
#Lineage2: 1  5  6  7  
#curves: 2 
#Curve1: Length: 27.763	Samples: 1650.63
#Curve2: Length: 28.084	Samples: 1195.76
lines(slingCurves(curves)[[1]], lwd = 2, lty=2, col=1)
lines(slingCurves(curves)[[2]], lwd = 2, lty=2, col=3)

pstime <-  slingPseudotime(curves)[rownames(so@meta.data),]
pstime <- as.data.frame(pstime)
pstime$clus <- so$seurat_clusters
so@meta.data <- cbind(so@meta.data, pstime)
boxplot(so$Lineage1~so$seurat_clusters)
FeaturePlot(so, reduction="umap", features=c("Lineage1","Lineage0"))
ggplot(so@meta.data, aes(y = CytoTRACE, x = Lineage0, color = Accession))+ geom_point(alpha=0.8)+ theme_bw()

LinMerge1 <- so@meta.data[so@meta.data$seurat_clusters %in% c(1,2,3,4), ]$Lineage1 - 5
names(LinMerge1) <- rownames(so@meta.data[so@meta.data$seurat_clusters %in% c(1,2,3,4), ])
LinMerge2 <- so@meta.data[so@meta.data$seurat_clusters %in% c(5,6,7),]$Lineage0
names(LinMerge2) <- rownames(so@meta.data[so@meta.data$seurat_clusters %in% c(5,6,7), ])
LineageM <- c(LinMerge1, LinMerge2)
so@meta.data$LineageM <- LineageM[rownames(so@meta.data)]
boxplot(so$LineageM~so$seurat_clusters, main="Merged Psuedotime ~ clusters")
p <- ggplot(so@meta.data, aes(x = Date.4, y = LineageM, color = month))
p + geom_point(alpha=0.8)+ theme_bw() #print1
p + geom_point(alpha=0.8, position = position_jitterdodge(0.1))+ facet_wrap(vars(Accession), ncol = 4L)+theme_bw() #print2

################################################
############ Characteristics of clusters
#1.
boxplot(meta[meta$seurat_clusters %in% "7",]$Date.4 ~ meta[meta$seurat_clusters %in% "7",]$Accession, main="fD2H~accs. Cluster7", ylim=c(0, 1.2), col=0, outline=F, ylab="fD2H", xlab="clusters")
beeswarm(meta[meta$seurat_clusters %in% "7",]$Date.4 ~ meta[meta$seurat_clusters %in% "7",]$Accession, main="fD2H~accs. Cluster7", ylim=c(0, 1.2), pch=19, cex=0.5, add=TRUE, col=ggcolor(4))
#2
FeaturePlot(so, features="Date.4", split.by="Accession", ncol=2)
#3 Stacking barplot of Accs composition of clusters.
p <- ggplot(meta, aes(x = seurat_clusters, fill = Accession))
p + geom_bar() + theme_bw()
#3a Stacking barplot wiht percentage
p <- ggplot(meta, aes(x = seurat_clusters, y="perc", fill = Accession))
p +  geom_col(position = "fill") +theme_bw()

############################################
############ Cluster DEG 
#############################################
#I. Cluster marker by MAST
integrated.markers <- FindAllMarkers(object=so, only.pos=TRUE, min.pct=0.25, test.use="MAST")
write.csv(integrated.markers,"./02_cluster-DEG/cluster.markers.MAST.csv")
integrated.markers %>%
  filter(avg_log2FC >= 0) %>%
  group_by(cluster) %>%
  top_n(5, avg_log2FC) -> top5.posi
#(visualize)
DoHeatmap(object = so, features=top5.posi$gene, size=3, label=TRUE) + NoLegend() + NoAxes()
 
#II. DEG by accession
j247.markers <- FindMarkers(so, group.by="Accession", ident.1="J247", only.pos=FALSE, min.pct=0.25, test.use="MAST")
write.csv(j247.markers,"./02_cluster-DEG/j247.markers.MAST.csv")
j247.markers %>% filter(avg_log2FC >= 1) %>% top_n(6, avg_log2FC) -> top6.posi
j247.markers %>% filter(avg_log2FC <= -1 & p_val_adj <= 0.01)  %>% top_n(-6, avg_log2FC) -> top6.neg
FeaturePlot(so, reduction="umap", features=rownames(top6.posi), ncol=3, col=c("grey","red"))
FeaturePlot(so, reduction="umap", features=rownames(top6.neg), ncol=3)

#III. DEG of cluster 2&3
deg.2vs3 <- FindMarkers(so, group.by="seurat_clusters", ident.1="2", ident.2="3", only.pos=FALSE, min.pct=0.25, test.use="MAST")
deg.2 <- FindMarkers(so, group.by="seurat_clusters", ident.1="2", only.pos=TRUE, min.pct=0.25, test.use="MAST")
deg.3 <- FindMarkers(so, group.by="seurat_clusters", ident.1="3", only.pos=TRUE, min.pct=0.25, test.use="MAST")
deg.2n3 <- FindMarkers(so, group.by="seurat_clusters", ident.1=c("2","3"), only.pos=TRUE, min.pct=0.25, test.use="MAST")
deg.4 <- FindMarkers(so, group.by="seurat_clusters", ident.1="4", only.pos=TRUE, min.pct=0.25, test.use="MAST")
deg.5 <- FindMarkers(so, group.by="seurat_clusters", ident.1="5", only.pos=TRUE, min.pct=0.25, test.use="MAST")
deg.4n5 <- FindMarkers(so, group.by="seurat_clusters", ident.1=c("4","5"), only.pos=TRUE, min.pct=0.25, test.use="MAST")
deg.4vs5 <- FindMarkers(so, group.by="seurat_clusters", ident.1="4", ident.2="5", only.pos=FALSE, min.pct=0.25, test.use="MAST")
write.csv(deg.2vs3,"./02_cluster-DEG/cl2vs3.markers.MAST.csv")
write.csv(deg.2,"./02_cluster-DEG/cl2-pos.markers.MAST.csv")
write.csv(deg.3,"./02_cluster-DEG/cl3-pos.markers.MAST.csv")
write.csv(deg.2n3,"./02_cluster-DEG/cl2n3-pos.markers.MAST.csv")
write.csv(deg.4,"./02_cluster-DEG/cl4-pos.markers.MAST.csv")
write.csv(deg.5,"./02_cluster-DEG/cl5-pos.markers.MAST.csv")
write.csv(deg.4n5,"./02_cluster-DEG/cl4n5-pos.markers.MAST.csv")
write.csv(deg.4vs5,"./02_cluster-DEG/cl4vs5.markers.MAST.csv")

#IV. DEG GO/KEGG
runGO2 <- function(gid_ls, prefix){
  out_ls = list()
  for (f in c("MF","CC","BP")){
    p <- GSEAGOHyperGParams(
      name = "Paramaters",
      geneSetCollection = gsc,
      geneIds = gid_ls,
      universeGeneIds = all,
      ontology = f,
      pvalueCutoff = 1,
      conditional = FALSE,
      testDirection = "over"
    )
    res <- hyperGTest(p)
    res_gid <- geneIdsByCategory(res)
    res2 <- summary(res)
    res2$GeneId <- sapply(res2[,1], function(x)paste(unlist(res_gid[x]), collapse="/"))
    
    out_ls = c(out_ls, list(res2))
  }
  write.csv(out_ls[[1]], paste0(prefix,".GO-MP.csv"), row.names=F)
  write.csv(out_ls[[2]], paste0(prefix,".GO-CC.csv"), row.names=F)
  write.csv(out_ls[[3]], paste0(prefix,".GO-BP.csv"), row.names=F)
  
  kegg.p <- GSEAKEGGHyperGParams(
    name = "Paramaters",
    geneSetCollection = gsck,
    geneIds = gid_ls, 
    universeGeneIds = kegg.all,
    pvalueCutoff = 1,
    testDirection = "over"
  )
  kegg.res <- hyperGTest(kegg.p)
  kegg.res_gid <- geneIdsByCategory(kegg.res)
  kegg.res2 <- summary(kegg.res)
  kegg.res2$Term <- kegg.annt[kegg.res2$KEGGID,]$annt
  kegg.res2$gid <- sapply(kegg.res2[,1], function(x)paste(unlist(kegg.res_gid[x]), collapse="/"))
  kegg.res2$KEGGID <- sapply(kegg.res2$KEGGID, function(x) paste0("map",x))
  write.csv(kegg.res2, paste0(prefix,".KEGG.csv"), row.names=FALSE)
} 
deg.sig <- subset(deg.2n3, deg.2n3$p_val_adj < 0.05 & deg.2n3$avg_log2FC >= 1)
runGO2(rownames(deg.sig), "./02_cluster-DEG/cl2n3-pos825")

deg.sig2 <- deg.2vs3[rownames(deg.2vs3) %in% rownames(deg.sig), ]
deg.sig2 <- subset(deg.sig2, deg.sig2$p_val_adj < 0.05 & abs(deg.sig2$avg_log2FC) >= 1)




############################################
############ TradeSeq-DE
#############################################
#0. Prepare the dataset, settings
curve <- readRDS("seurat-curves.barley-field.leaf.v2.rds")
counts <- as.matrix(so@assays$RNA@counts)
#counts <- subset(counts, rownames(counts) %in% VariableFeatures(so))
#[1] 2963 1940
BPPARAM <- BiocParallel::bpparam() #for Intel-Mac
BPPARAM$workers <- 4

#1. replace the curve for interest
curve@assays@data$pseudotime[,1] <- meta$fD2H
curve@assays@data$pseudotime[,1] <- meta$Lineage0
curve@assays@data$pseudotime[,1] <- meta$LineageM
curve@assays@data$pseudotime[,1] <- meta$CytoTRACE

na.cell <- !is.na(curve@assays@data$pseudotime)
cw <- as.matrix(curve@assays@data$weights[na.cell,])
pstime <- as.matrix(curve@assays@data$pseudotime[na.cell,])
count <-counts[, rownames(cw)]
#count <- meta[rownames(cw), c("AvgTemp","DayHour")]
#(for conditional fitGAM)
accs.cond <- meta[colnames(count), "Accession"]
accs.cond <-factor(accs.cond, level=c("J247", "J064", "H602","J647"))
levels(accs.cond) <- c(1,2,2,2) #set J247 as 1, the others as 2

#2. choose best-K for fitGAM(): take minutes ~ hours
icMat <- evaluateK(counts=count, pseudotime=pstime, cellWeight=cw, k=3:15, nGenes=200,verbose=F, plot=T, parallel=T, BPPARAM=BPPARAM)
bestK = 6 #fD2H, LineageM
bestK = 7 #Lineage0
bestK = 8 #CytoTRACE

  
#3. fitGAM()
gam = fitGAM(counts=count, pseudotime=pstime, cellWeight=cw, 
             nknots=bestK, parallel=T, BPPARAM=BPPARAM)
gam.assoc = associationTest(gam, lineages=T, global=T, l2fc=1)
gam.assoc <- gam.assoc[order(gam.assoc$pvalue, -gam.assoc$meanLogFC), ]
dim(subset(gam.assoc, gam.assoc$pvalue < 0.001 & abs(gam.assoc$meanLogFC) >= 1))[1]
#[1] 632 

gam.j247 = fitGAM(counts=count, pseudotime=pstime, cellWeight=cw, 
                conditions=accs.cond, nknots=bestK, parallel=T, BPPARAM=BPPARAM)
deg.j247 <- conditionTest(gam.j247, lineages=T, l2fc=1)
length(deg.j247[deg.j247$pvalue < 0.001, 1])
#[1] 342
gam.assoc$j247.pval <- deg.j247[rownames(gam), ]$pvalue
write.csv(gam.assoc, "03_TradeSeq-assoc/Tradeseq-gam.fD2H-assoc.csv")


#4. Smoothering the exp-data on the GAM
preSm <- predictSmooth(gam, gene=rownames(count), nPoints=100, tidy=FALSE) #tidy=TRUE for matches of pstime~100-linear windows
write.csv(preSm, "03_TradeSeq-assoc/GAM-LM.J2vsOther.100win2.csv")
#normalisation & standardisation of the smoothed data
x1 <- 100* preSm/ apply(preSm, 1, sum) #converted to the percentage
x2 <- 100* preSm / apply(preSm, 1, function(y) sum(abs(y))) # L1 normalisation
x3 <- preSm - apply(preSm, 1, mean) # geneMean, centralization
x4 <- (preSm - apply(preSm, 1, mean) ) / apply(preSm, 1, sd)  # Z-standardization simply x4=t(scale(t(preSm)))
x <- x4
pheatmap(x, cluster_cols = FALSE, cutree_rows = bestK, show_rownames = FALSE, show_colnames = FALSE
         , main = "fD2H smoothered, z-scaled vari3k") #(first veiw)

#5. k-mean with clusGAP() &  H-clustering
#(https://zenn.dev/rchiji/books/55a7b6a1777e4d/viewer/6b2959)
require(cluster)
cg.res <- clusGap(x, kmeans, K.max=20, B=100, verbose=interactive()) ##taking minutes
#plot(cg.res)
#cg.res
nClus = 6 #fD2H w/ z-scale
clus <- kmeans(x = x, centers = nClus, iter.max = 50)
centers <- clus$centers-apply(clus$centers, 1, mean)
hclus <- hclust(as.dist(1-cor(t(centers), method = "pearson")), method = "average")
pheatmap(x, cluster_cols = FALSE, cluster_rows = hclus, show_rownames = FALSE, show_colnames = FALSE
         , main = "fD2H smoothered, z-scaled, hclust-mean") #(hclust-mean veiw)
#exchange the "kmean cluster" to "hcluster" order
nOrder <- match(clus$cluster, hclus$order) 
Kmeans_matrix <- x[order(nOrder),]
Kmeans_cluster <- sort(nOrder)
ngenes <- as.character(table(Kmeans_cluster))
#[1] "385" "343" "712" "577" "359" "587"
preSm2 <- preSm[rownames(Kmeans_matrix), ]
gam.assoc2 <- gam.assoc[rownames(Kmeans_matrix), ]
gam.assoc2$clus <- Kmeans_cluster
gam.preSm <- cbind(preSm2, gam.assoc2)
write.csv(gam.preSm, "03_TradeSeq-assoc/Tradeseq-gam.fD2H-clus-smth.csv")

#6. heatmap.2 & Linear median plots with random-choose 200.
n_randomSample = 200
ix <- sort(sample(1:length(Kmeans_cluster), n_randomSample))
Kmeans_cluster2 <- Kmeans_cluster[ix]
Kmeans_matrix2 <- Kmeans_matrix[ix,]
#(polishing the view by removing the out-liners)
Kmeans_matrix2 <- as.matrix(Kmeans_matrix2) - apply(Kmeans_matrix2, 1, mean) #centralization
cutoff <- median(unlist(Kmeans_matrix2)) + 3*sd (unlist(Kmeans_matrix2)) #remove median +3*sd
Kmeans_matrix2[Kmeans_matrix2 >cutoff] <- cutoff 
cutoff <- median(unlist(Kmeans_matrix2)) - 3*sd (unlist(Kmeans_matrix2)) #remove median +3*sd
Kmeans_matrix2[Kmeans_matrix2 < cutoff] <- cutoff
doHM3 <- function(mtx, cls){
  #mtx <- t(scale(t(mtx))) # z-standardization
  #breaks = seq(-3, 3, by=0.1)
  ngenes <- as.character(table(cls))
  hm <- heatmap.2(x = mtx,
                  Rowv =F,
                  Colv=F,
                  dendrogram ="none",
                  density.info="none",
                  trace="none",
                  scale="none",
                  keysize=.3,
                  key= F, 
                  labRow = F,
                  labCol = F,
                  RowSideColors = rainbow(length(unique(cls)))[as.factor(cls)],
                  margins = c(8, 16), #for labels
                  srtCol=45,
                  col=key.palette,
  )
  legend.text <- paste(toupper(letters)[unique(cls)],
                       " (N=", ngenes,")",
                       sep="")
  par(lend = 1) # square line ends for the color legend
  legend(x = "bottomright", 
         legend = legend.text, # category labels
         col = rainbow(length(unique(cls))),  # color key
         lty= 1, # line style
         lwd = 10 # line width
  )
  
  legend_image <- as.raster(rev(hm$colorTable$color), ncol=1)
  rasterImage(legend_image, 0.8, 0.8, 0.85, 1)
  
  return(hm)
} #wo/ z-scaling
require(gplots)
hm <- doHM3(Kmeans_matrix2, Kmeans_cluster2)
doLP3 <- function(mtx, cls){
  mtx <- as.matrix(mtx) - apply(mtx, 1, mean) #centralization
  ymax = 0
  ymin = 0
  Ncls = length(unique(cls))
  par(mfrow = c(3,3))
  
  for(i in unique(cls)){
    submtx <- mtx[cls %in% i, ]
    submean <- apply(submtx, 2, mean)
    submax <- max(submean)
    submin <- min(submean)
    if (submax > ymax){ymax = submax}
    if (submin < ymin){ymin = submin}
  }
  for(i in unique(cls)){
    submtx <- mtx[cls %in% i, ]
    submean <- apply(submtx, 2, mean)
    title = paste("Cluster ", i, " (N=", dim(submtx)[1],")", sep="")
    if (length(submean) > 100){
      submean1 <- submean[1:100] ## for j247 vs other
      submean2 <- submean[101:200]
      submean <- as.matrix(data.frame(submean1, submean2))
      matplot(submean, type="l", lwd=4, lty=c(1,3), col=rainbow(Ncls)[i],
              ylab="Mean z-scaled exp", ylim=c(ymin, ymax), main=title)
    } else {
      matplot(submean, type="l", lwd=4, lty=1, col=rainbow(Ncls)[i],
            ylab="Mean z-scaled exp", ylim=c(ymin, ymax), main=title)
    }
    abline(h=0, lty=2)
  }
}
doLP3(Kmeans_matrix, Kmeans_cluster)
dev.off()

#7. J247-DEG heatmap.2 & Linear median plots
j247.deg <- rownames(deg.j247[deg.j247$pvalue < 0.001, ])
#[1] 342
preSm <- predictSmooth(gam.j247, gene=j247.deg, nPoints=100, tidy=FALSE)
x <- t(scale(t(preSm))) # z-scaling
############################### old-keep the clusters from the whole data.
cls <- gam.preSm[rownames(preSm), "clus"]
names(cls) <- rownames(preSm)
nOrder <- match(cls, sort(unique(cls)))
Kmeans_matrix3 <- x[order(nOrder), ]
Kmeans_cluster3 <- sort(cls)
############################### New cluster identification
cg.res <- clusGap(x, kmeans, K.max=20, B=100, verbose=interactive()) ##taking minutes
nClus = 6 #
clus <- kmeans(x = x, centers = nClus, iter.max = 50)
centers <- clus$centers-apply(clus$centers, 1, mean)
hclus <- hclust(as.dist(1-cor(t(centers), method = "pearson")), method = "average")
pheatmap(x, cluster_cols = FALSE, cluster_rows = hclus, show_rownames = FALSE, show_colnames = FALSE
         , main = "j247 vs. other, z-scaled, hclust-mean") #(hclust-mean veiw)
nOrder <- match(clus$cluster, hclus$order) 
Kmeans_matrix <- x[order(nOrder),]
Kmeans_cluster <- sort(nOrder)
ngenes <- as.character(table(Kmeans_cluster))
#[1] "60" "60" "62" "33" "40" "87"
preSm2 <- preSm[rownames(Kmeans_matrix), ]
gam.assoc2 <- gam.assoc[rownames(Kmeans_matrix), ]
gam.assoc2$clus <- Kmeans_cluster
gam.preSm <- cbind(preSm2, gam.assoc2)
write.csv(gam.preSm, "03_TradeSeq-assoc/Tradeseq-gam.j247-other.csv")
Kmeans_cluster3 <- Kmeans_cluster
Kmeans_matrix3 <- Kmeans_matrix

################## (common) visualization
Kmeans_matrix3 <- as.matrix(Kmeans_matrix3) - apply(Kmeans_matrix3, 1, mean) #centralization
cutoff <- median(unlist(Kmeans_matrix3)) + 3*sd (unlist(Kmeans_matrix3)) #remove median +3*sd
Kmeans_matrix3[Kmeans_matrix3 >cutoff] <- cutoff 
cutoff <- median(unlist(Kmeans_matrix3)) - 3*sd (unlist(Kmeans_matrix3)) #remove median +3*sd
Kmeans_matrix3[Kmeans_matrix3 < cutoff] <- cutoff
hm <- doHM3(Kmeans_matrix3, Kmeans_cluster3)
doLP3(Kmeans_matrix, Kmeans_cluster)
dev.off()


################ GO analysis
### set function to GO-GSEA
RunGO <- function(gid_ls){
  out_ls = list()
  for (f in c("MF","CC","BP")){
    p <- GSEAGOHyperGParams(
      name = "Paramaters",
      geneSetCollection = gsc,
      geneIds = gid_ls,
      universeGeneIds = all,
      ontology = f,
      pvalueCutoff = 1,
      conditional = FALSE,
      testDirection = "over"
    )
    res <- hyperGTest(p)
    res_gid <- geneIdsByCategory(res)
    res2 <- summary(res)
    res2$GeneId <- sapply(res2[,1], function(x)paste(unlist(res_gid[x]), collapse="/"))
    
    out_ls = c(out_ls, list(res2))
  }
  
  return(out_ls)
} 
out1 <- RunGO(rownames(Kmeans_matrix3))
out1.MP <- out1[[1]]
out1.CC <- out1[[2]]
out1.BP <- out1[[3]]


##################################################
################## DEG visualization & characterization





















#IV. Tradeseq - to fD2H
curve@assays@data$pseudotime[,1] <- so@meta.data$fD2H
na.cell <- !is.na(curve@assays@data$pseudotime)
cw2 <- as.matrix(curve@assays@data$weights[na.cell,])
pstime2 <- as.matrix(curve@assays@data$pseudotime[na.cell,])
counts2 <- counts[,rownames(cw2)]
acc.cond <- so@meta.data[colnames(counts2), "Accession"]
#acc.cond <-factor(acc.cond, level=c("J247", "J064", "H602","J647"))
levels(acc.cond) <- c(1,2,2,2) #set J247 as 1, the others as 2
#(GAM-modeling)
icMat3 <- evaluateK(counts = counts2, pseudotime = pstime2, cellWeight=cw2, k = 3:15, nGenes = 200, 
                    verbose = FALSE, plot = TRUE, parallel=TRUE, BPPARAM=BPPARAM)
bestK=7
gam3 <- fitGAM(counts = counts2, pseudotime = pstime2, cellWeight=cw2, 
               nknots=bestK, parallel=TRUE, BPPARAM=BPPARAM, conditions=acc.cond)
gam3.assoc <- associationTest(gam3, lineages=TRUE, global=TRUE, l2fc=1)
write.csv(gam3.assoc,"./03_TradeSeq-assoc/Tradeseq-gam.fD2H.csv")
gam3.j247 <- conditionTest(gam3, pairwise=TRUE, global=TRUE, l2fc=1)
write.csv(gam3.j247,"./03_TradeSeq-assoc/Tradeseq-gam.fD2H-j247.csv")
gam3.j247s <- subset(gam3.j247, gam3.j247$pvalue < 0.001)
#[1] 334   3
gam3.deg <- subset(gam3.assoc, gam3.assoc$pvalue < 0.001 & abs(gam3.assoc$meanLogFC) >= 1)
#[1] 563   7

#(Smoothing the exp-data on the timeline)
preSm <- predictSmooth(gam3, gene=rownames(gam3.deg), nPoints=100, tidy=FALSE) #tidy=TRUE for matches of pstime~100-linear windows
#(normalization the exp data)
deg.exp2 <- 100* preSm/ apply(preSm, 1, sum) #converted to the percentage
deg.exp3 <- 100* preSm / apply(preSm, 1, function(y) sum(abs(y))) # L1 normalisation
deg.exp4 <- preSm - apply(preSm, 1, mean) # geneMean, centralization
deg.exp5 <- (preSm - apply(preSm, 1, mean) ) / apply(preSm, 1, sd)  # gene standardization
x <- deg.exp4 ### with no reason
pheatmap(t(scale(t(x))), cluster_cols = FALSE, cutree_rows = 7,show_rownames = FALSE, show_colnames = FALSE
         , main = "fD2H smoothered, centralized DEG 563, z-scaled vari3k") #Initial VIEW w/ z-standardization

#(K-mean clustering with clusGAP)
#(https://zenn.dev/rchiji/books/55a7b6a1777e4d/viewer/6b2959)
require(cluster)
cg.res <- clusGap(x, kmeans, K.max=20, B=100, verbose=interactive()) ##taking minutes
#plot(cg.res)
#cg.res
nClusters=5
clus <- kmeans(x = x, centers = nClusters, iter.max = 50)
#stepwise cluster by the mean of each cluster
centers <- clus$centers-apply(clus$centers, 1, mean)
hclus <- hclust(as.dist(1-cor(t(centers), method = "pearson")), method = "average")
#plot(hclus)
pheatmap(t(scale(t(x))), cluster_cols = FALSE, cluster_rows = hclus, show_rownames = FALSE, show_colnames = FALSE
         , main = "fD2H smoothered, Hclust-mean, z-scaled vari3k") #Hclust-mean view
#exchange the "kmean cluster" to "hcluster" order
nOrder <- match(clus$cluster, hclus$order) 
Kmeans_matrix <- x[order(nOrder),]
Kmeans_cluster <- sort(nOrder)
ngenes <- as.character(table(Kmeans_cluster))
#[1] "30"  "61"  "113" "38"  "321"

#(visulaization by random-200)
n_randomSample = 200
ix <- sort(sample(1:length(Kmeans_cluster), n_randomSample))
Kmeans_cluster2 <- Kmeans_cluster[ix]
Kmeans_matrix2 <- Kmeans_matrix[ix,]
#[1] 400 100
#(polishing the view by removing the out-of-liners 3*sd of whole genes)
Kmeans_matrix2 <- as.matrix(Kmeans_matrix2) - apply(Kmeans_matrix2, 1, mean)
cutoff <- median(unlist(Kmeans_matrix2)) + 3*sd (unlist(Kmeans_matrix2)) #remove median +3*sd
Kmeans_matrix2[Kmeans_matrix2 >cutoff] <- cutoff 
cutoff <- median(unlist(Kmeans_matrix2)) - 3*sd (unlist(Kmeans_matrix2)) #remove median +3*sd
Kmeans_matrix2[Kmeans_matrix2 < cutoff] <- cutoff
Kmeans_matrix2 <- Kmeans_matrix[rownames(Kmeans_matrix2),]
#[1] 200 100
#(heatmap)
hm <- doHM3(Kmeans_matrix2, Kmeans_cluster2)

#(save DEG matrix with cluster-info)
preSm2 <- preSm[rownames(Kmeans_matrix), ]
gam3.deg <- gam3.deg[rownames(Kmeans_matrix), ]
gam3.deg$clus5 <- Kmeans_cluster
gam3.deg2 <- cbind(gam3.deg, preSm2)
write.csv(gam3.deg2,"./03_TradeSeq-assoc/Tradeseq-gam.fD2H-DEG.csv")


#V. Tradeseq - to LM
curve@assays@data$pseudotime[,1] <- so@meta.data$LineageM
na.cell <- !is.na(curve@assays@data$pseudotime)
cw2 <- as.matrix(curve@assays@data$weights[na.cell,])
pstime2 <- as.matrix(curve@assays@data$pseudotime[na.cell,])
counts2 <- counts[,rownames(cw2)]
#(GAM-modeling)
icMat4 <- evaluateK(counts = counts2, pseudotime = pstime2, cellWeight=cw2, k = 3:15, nGenes = 200, 
                    verbose = FALSE, plot = TRUE, parallel=TRUE, BPPARAM=BPPARAM)
bestK=7
gam4 <- fitGAM(counts = counts2, pseudotime = pstime2, cellWeight=cw2, nknots=bestK
               , parallel=TRUE, BPPARAM=BPPARAM)
gam4.assoc <- associationTest(gam4, lineages=TRUE, global=TRUE, l2fc=1)
write.csv(gam4.assoc,"./03_TradeSeq-assoc/Tradeseq-gam.LineageM.csv")
gam4.deg <- subset(gam4.assoc, gam4.assoc$pvalue < 0.001 & abs(gam4.assoc$meanLogFC) >= 1)
#[1] 549   7

#(Smoothing the exp-data on the timeline)
preSm <- predictSmooth(gam4, gene=rownames(gam4.deg), nPoints=100, tidy=FALSE) #tidy=TRUE for matches of pstime~100-linear windows
#(normalization the exp data)
deg.exp2 <- 100* preSm/ apply(preSm, 1, sum) #converted to the percentage
deg.exp3 <- 100* preSm / apply(preSm, 1, function(y) sum(abs(y))) # L1 normalisation。各行の合計を100にする。
deg.exp4 <- preSm - apply(preSm, 1, mean) # geneMean, centralization。中心化。平均を引く。
deg.exp5 <- (preSm - apply(preSm, 1, mean) ) / apply(preSm, 1, sd)  # gene standardization。標準化。平均を引いて、SDで割る。
x <- deg.exp4 ### with no reason
pheatmap(t(scale(t(x))), cluster_cols = FALSE, cutree_rows = 5,show_rownames = FALSE, show_colnames = FALSE
         , main = "LineageM smoothered, centralized DEG 549, z-scaled vari3k") #Initial VIEW w/ z-standardization

#(K-mean clustering with clusGAP)
#(https://zenn.dev/rchiji/books/55a7b6a1777e4d/viewer/6b2959)
require(cluster)
cg.res <- clusGap(x, kmeans, K.max=20, B=100, verbose=interactive()) ##taking minutes
#plot(cg.res)
#cg.res
nClusters=5
clus <- kmeans(x = x, centers = nClusters, iter.max = 100)
#stepwise cluster by the mean of each cluster
centers <- clus$centers-apply(clus$centers, 1, mean)
centers.order <- names(sort(apply(clus$centers, 1, which.max)))
hclus <- hclust(as.dist(1-cor(t(centers), method = "pearson")), method = "average")
#plot(hclus)
pheatmap(t(scale(t(x))), cluster_cols = FALSE, cluster_rows = hclus,show_rownames = FALSE, show_colnames = FALSE
         , main = "LineageM smoothered, Hclust-mean, z-scaled vari3k") #Hclust-mean view

nOrder <- match(clus$cluster, hclus$order) #re-order the "kmean cluster" to "hcluster" order
nOrder2 <- match(clus$cluster, centers.order) #re-order the "kmean cluster" to the median-peak appearance
Kmeans_matrix <- x[order(nOrder),]
Kmeans_cluster <- sort(nOrder)
ngenes <- as.character(table(Kmeans_cluster))
#[1] "44"  "67"  "92"  "53"  "293"

#(visulaization by random-200)
n_randomSample = 200
ix <- sort(sample(1:length(Kmeans_cluster), n_randomSample))
Kmeans_cluster2 <- Kmeans_cluster[ix]
Kmeans_matrix2 <- Kmeans_matrix[ix,]
#[1] 400 100
#(polishing the view by removing the out-of-liners)
Kmeans_matrix2<-as.matrix(Kmeans_matrix2)-apply(Kmeans_matrix2,1,mean)
cutoff <- median(unlist(Kmeans_matrix2)) + 3*sd (unlist(Kmeans_matrix2)) #remove median +3*sd
Kmeans_matrix2[Kmeans_matrix2 > cutoff] <- cutoff 
cutoff <- median(unlist(Kmeans_matrix2)) - 3*sd (unlist(Kmeans_matrix2)) #remove median +3*sd
Kmeans_matrix2[Kmeans_matrix2 < cutoff] <- cutoff
Kmeans_matrix3 <- x[rownames(Kmeans_matrix2),]
#[1] 200 100
#(heatmap)
hm <- doHM3(Kmeans_matrix2, Kmeans_cluster2)
hmz <- doHM3z(Kmeans_matrix3, Kmeans_cluster2)
#(save DEG matrix with cluster-info)
preSm2 <- preSm[rownames(Kmeans_matrix), ]
gam4.deg <- gam4.deg[rownames(Kmeans_matrix), ]
gam4.deg$clus5 <- Kmeans_cluster
gam4.deg2 <- cbind(gam4.deg, preSm2)
write.csv(gam4.deg2,"./03_TradeSeq-assoc/Tradeseq-gam.LineageM-DEG.csv")









##############################################
###### UMAP & PHATE by accession
j247 <- subset(so, subset = Accession == "J247")
j064 <- subset(so, subset = Accession == "J064")
h602 <- subset(so, subset = Accession == "H602")
j647 <- subset(so, subset = Accession == "J647")
other <- subset(so, subset = Accession != "J247")

j247 <- RunPCA(j247, vervose=FALSE, approx=FALSE, npcs=50, features=VariableFeatures(j247))
ElbowPlot(j247, ndims=50, reduction = "pca") + ggtitle("J247, vari3k, pca50")
j247 <- RunUMAP(j247, reduction="pca", dims=c(1:50), umap.method = "umap-learn", metric = "correlation", n.components = 2)
FeaturePlot(j247, reduction="umap", features=c("fD2H")) + ggtitle("J247, vari3k, pca50")
##(not linear, selecting PCs by fD2H GLR)
alpha <- seq(0.01, 0.99, 0.01)
mse.df <- NULL
for (i in 1:length(alpha)) {
  m <- cv.glmnet(x = j247@reductions$pca@cell.embeddings, y = j247@meta.data$fD2H, family = "gaussian", alpha = alpha[i])
  mse.df <- rbind(mse.df, data.frame(alpha = alpha[i], mse = min(m$cvm)))
}
best.alpha <- mse.df$alpha[mse.df$mse == min(mse.df$mse)]
#[1] 0.17 j247
#[1] 0.48 h602
#[1] 0.4 j064
#[1] 0.2 j647
#[1] 0.86 not-j247
en.model.cv <- cv.glmnet(x = j247@reductions$pca@cell.embeddings, y = j247@meta.data$fD2H, family = "gaussian", alpha = best.alpha)
plot(en.model.cv, xvar="lambda", label=TRUE, main="gene:Vari3K-PCA50-EN (J247)")
en.model <- glmnet(x = j247@reductions$pca@cell.embeddings, y = j247@meta.data$fD2H, family = "gaussian", lambda = en.model.cv$lambda.1se)
plot(abs(en.model$beta), main="betas of PC-dim, EN-model (J247)", pch=19)
abline(h=0.005, col=2, lty=2)
PCs <- which(abs(en.model$beta)>0.005)
#[1]  1  3  6 10 11 15 17 21 25 $J247 0.005
#[1]  1  3  7 10 14 18 23 27 43 $H602 0.007
#[1]  1  3  5  8 11 13 20 21    $J064 0.006
#[1]  1  3  7 11 12 16 26 28 35 40 $J647 0.006
#[1]  1  4 13 14 15 16 18 22 23 26 38 39 45 #not-J247 0.006
sum(sqrt(Stdev(so, reduction="pca")[PCs]))
#[1] 23.79911 vs. [1] 98.01709 #J247
#[1] 22.52408 vs. [1] 98.01709 #H602
#[1] 22.52428 vs. [1] 98.01709 #J064
#[1] 24.1218 vs. [1] 98.01709 #J647
#[1] 27.70766 vs. [1] 98.01709 #not-J247
##(UMAP again)
j247 <- RunUMAP(j247, reduction="pca", dims=PCs, umap.method = "umap-learn", metric = "correlation", n.components = 2)
FeaturePlot(j247, reduction="umap", features=c("fD2H","CytoTRACE"))
DimPlot(j247, reduction="umap", group.by=c("month","seurat_clusters"))

##(PHATE trial; optional)####
j247.exp <- t(as.data.frame(j247@assays$RNA@data))
j247.phate <- phate(j247.exp)
plot(j247.phate)
j247.phate2 <- phate(j247.exp, gamma=1, t=120)
plot(j247.phate2)
j247[["phate"]] = CreateDimReducObject(embeddings = j247.phate$embedding * 1000, 
                                        key = "PHATE_", assay = DefaultAssay(j247))
FeaturePlot(j247, reduction="phate", features="fD2H")
######





















############## Identify genes on curves
plot_differential_expression <- function(feature_id) {
  #feature_id <- pseudotime_association %>% filter(pvalue < 0.05) %>% top_n(1, -waldStat) %>% pull(feature_id)
  title <- cowplot::ggdraw() + cowplot::draw_label(feature_id, fontface="bold", x=0.9, hjust=1) + ggplot2::theme(plot.margin=margin(0,0,0,7))
  genecount <- plotGeneCount(curve_sub, count_sub, gene = feature_id[1], clusters = clus_sub, models = sce_sub) + ggplot2::theme(legend.position = "none")
  smoothers <- plotSmoothers(sce_sub, as.matrix(count_sub), gene = feature_id[1])
  cowplot::plot_grid(title, cowplot::plot_grid(genecount, smoothers), ncol = 1, rel_heights = c(0.1, 1))
}

pseudotime_association <- associationTest(sce_sub, lineages=TRUE) ##test for all lineages simultaneously.
#pseudotime_association$fdr_1 <- p.adjust(pseudotime_association$pvalue_1, method = "fdr")
#pseudotime_association$fdr_2 <- p.adjust(pseudotime_association$pvalue_2, method = "fdr")
pseudotime_association <- pseudotime_association[order(pseudotime_association$pvalue), ]
pseudotime_association$feature_id <- rownames(pseudotime_association)
#feature_id <- pseudotime_association %>% filter(fdr < 0.05) %>% top_n(1, -waldStat) %>% pull(feature_id)
plot_differential_expression(feature_id)

### Genes changing with pseudo-time
startRes <- startVsEndTest(sce_sub, global=TRUE, lineages=TRUE)
oStart <- order(startRes$waldStat, decreasing = TRUE)
sigGeneStart <- names(sce_sub)[oStart[3]]
plotSmoothers(sce_sub, count_sub, gene = sigGeneStart)
plotGeneCount(curve_sub, count_sub, gene = sigGeneStart)
#if focusing on two pseudotime points
customRes <- startVsEndTest(sce, pseudotimeValues = c(0.1, 0.8))

### Genes changing at the terminal points
endRes <- diffEndTest(sce_sub, global=TRUE, pairwise=TRUE)
o <- order(endRes$waldStat, decreasing = TRUE)
sigGene <- names(sce_sub)[o[1]]
plotSmoothers(sce_sub, count_sub, sigGene)
plotGeneCount(curve_sub, count_sub, gene = sigGene)

### Genes changing differently between lineages
patternRes <- patternTest(sce_sub, global=TRUE, pairwise=TRUE, l2fc = log2(2), eigenThresh = 0.01)
patternRes$fdr <- p.adjust(patternRes$pvalue, method="fdr")
oPat <- order(patternRes$waldStat, decreasing = TRUE)
sigGene <- names(sce_sub)[oPat[1]]
plotSmoothers(sce_sub, count_sub, gene=sigGene) + ggtitle(sigGene)
plotGeneCount(curve_sub, count_sub, gene=sigGene)













##################################
###### Subset #1 & #5 for details of their transient states
so15 <- subset(so, idents = c(1,5))
soN15 <- subset(so, idents = c(2,3,4,6,7))
#An object of class Seurat 
#25145 features across 503 samples within 1 assay 
#Active assay: RNA (25145 features, 2842 variable features)
#3 layers present: counts, data, scale.data
#14 dimensional reductions calculated: pca, pca_whole, pca_vari, umap, umap_50, umap_3D, umap_2D, umap_whole50, umap_vari50, pca_flor, umap_flor50, pca_vari2, umap_vari2_50, umap_whole_glr7
DimPlot(so15, group.by=c("Case","Accession","month","seurat_clusters"))
FeaturePlot(so15, features=c("CytoTRACE","Date.4","LineageM","AvgTemp","LightHour","SunHour"), ncol=3)

#3(compare AvgTemp in cls1&5 and others)
t15 <- so15@meta.data
t15$class <- "t15"
tA <- so@meta.data
tA$class <- "tA"
tN <- soN15@meta.data
tN$class <- "tN"

tM <- rbind(tA, t15, tN)
#[1] 3880    5
tM$class <- factor(tM$class, levels=c("tA","t15","tN"))
boxplot(tM$AvgTemp ~ tM$class + tM$month, at=c(1:3, 5:7, 9:11, 13:15, 17:19, 21:23),col=ggcolor(3), 
        names=c("","12","","","1","","","2","","","3","","","4","","","5",""), xlab="month", ylab="Avg. Air. Temp.")
beeswarm(tM$AvgTemp ~ tM$class + tM$month, at=c(1:3, 5:7, 9:11, 13:15, 17:19, 21:23),col="grey40", pch=1, cex=.7, corral="wrap", add=TRUE)
legend("topleft", fill = ggcolor(3), legend = c("whole","cls1&5","others"), horiz = T)
require(vioplot)
vioplot(tM$AvgTemp~tM$class, col=ggcolor(3), main="Avg.Temp ~ clusters")

table(tM[,c(2,3)])
#     class
#month  tA t15  tN
#12 334 208 126
#1  383  56 327
#2  408 107 301
#3  455 113 342
#4  330  19 311
#5   30   0  30
        
#3(compare LightHour in cls1&5 and others of the temperature range in 1-5oC)
tM.s1 <- subset(tM, tM$AvgTemp >= 5 & tM$AvgTemp <= 10)
#[1] 1580    5
require(vioplot)
vioplot(tM.s1$DayHour~tM.s1$class, col=ggcolor(3), main="DayHour ~ clusters; in 1-5C temp. range")
vioplot(tM.s1$LightHour~tM.s1$class, col=ggcolor(3), main="LightHour ~ clusters; in 1-5C temp. range")
table(tM.s1[,5])
#tA t15  tN 
#790 327 463 
table(tM.s1$Accession)
#J247 J064 H602 J647 
# 197  198  198  197 
table(tM.s1$Season)
#y2017 y2018 
#  216   574  <<<- Two-time different
table(meta15$Season)
#y2017 y2018 
#156   347 

#4. statistical analysis
t15l <- tM.s1[tM.s1$class %in% "t15",]$LightHour
tNl <- tM.s1[tM.s1$class %in% "tN",]$LightHour
tAl <- tM.s1[tM.s1$class %in% "tA",]$LightHour
t15d <- tM.s1[tM.s1$class %in% "t15",]$DayHour
tNd <- tM.s1[tM.s1$class %in% "tN",]$DayHour
tAd <- tM.s1[tM.s1$class %in% "tA",]$DayHour
shapiro.test(t15l)
#p-value = 2.847e-13 # not-normalty
shapiro.test(tNl)
#p-value = 3.382e-14 # not-normalty
var.test(x=t15l,y=tNl)
#p-value = 0.6384 # homogeic variance for Light
#p-value = 2.161e-05 #non-homogenic variand for DayTime data
##Welch's t-test or U test
t.test(x=t15l, y=tNl, var.equal=F, paired=F) #Welch Two Sample t-test
#p-value = 5.668e-11 #enough different for t15 vs. tN
#p-value = 2.203e-05 #enough different for t15 vs. tA
#p-value = 0.0008971 # less different for tN vs. tA
wilcox.test(t15l, tNl) #Mann-Whitney U test
#p-value = 1.433e-11 #enough different for Light data
wilcox.test(tAl, tNl)
#p-value = 0.0005559
wilcox.test(tAl, t15l)
#p-value = 1.365e-05
wilcox.test(t15d, tNd) 
#p-value = 5.459e-07 #enough different to each other
wilcox.test(t15d, tAd) 
#p-value = 0.001253
wilcox.test(tNd, tAd) 
#p-value = 0.01047
##### well, analyzing metadata is over. struck. not that fascinating

#(find markers)
mk1.5 <- FindMarkers(so, ident.1 = c(1,5), min.pct = 0.25, only.pos = TRUE, test.use="MAST")
mk1 <- FindMarkers(so, ident.1 = c(1), min.pct = 0.25, only.pos = T, test.use="MAST")
mk5 <- FindMarkers(so, ident.1 = c(5), min.pct = 0.25, only.pos = T, test.use="MAST")
write.csv(mk1.5, "MAST-marker.cls1n5.csv")
write.csv(mk1, "MAST-marker.cls1.csv")
write.csv(mk5, "MAST-marker.cls5.csv")

mk1.5 <- rownames(mk1.5[mk1.5$avg_log2FC > 1 & mk1.5$p_val_adj < 0.05,])
#[1] 52  5
mk1 <- rownames(mk1[mk1$avg_log2FC > 1 & mk1$p_val_adj < 0.05,])
#[1] 37  5
mk5 <- rownames(mk5[mk5$avg_log2FC > 1 & mk5$p_val_adj < 0.05,])
#[1] 124   5
mk1.5u <- unique(c(mk1.5, mk1, mk5))
#[1] 159
plot <- DoHeatmap(object = so, features=mk1.5u, size=3, label=TRUE) + NoLegend() + NoAxes()

degAll <- FindMarkers(so15, group.by="Accession", ident.1="J247", test.use="MAST")
write.csv(degAll, "MAST-marker.cls1n5.J247.csv")
degAll <- rownames(degAll[degAll$p_val_adj < 0.5 & degAll$avg_log2FC >=2,])
DoHeatmap(object = so, features=degAll, size=3, label=TRUE) + NoLegend() + NoAxes()



#####
degAll <- FindMarkers(tbt1, group.by="orig.ident", ident.1="fugu5_8DAS", ident.2="WT_8DAS", test.use="MAST")
degAll$cluster <- "G"
degAll$gene <- rownames(degAll)

for(i in levels(Idents(tbt1))){
  print(i)
  deg_sub <- FindMarkers(tbt1, group.by="orig.ident", ident.1="fugu5_8DAS", ident.2="WT_8DAS", test.use="MAST", subset.ident=i)
  deg_sub$cluster <- i
  deg_sub$gene <- rownames(deg_sub)
  degAll <- rbind(degAll, deg_sub)
}
dim(degAll)
#[1] 1457    7
length(unique(degAll$gene))
#[1] 634
write.csv(degAll, "./3_Marker-genes/DEG-byFugu5.all.csv", quote=FALSE, row.names=FALSE)
#####



#############################################################
######## phate DR
require(phateR)
so45i <- as.data.frame(so45@assays$RNA@data)
so45i <- t(so45i)

so45.phate <- phate(so45i)
plot(so45.phate)
so45.phate2 <- phate(so45i, gamma=0, t=120, init=so45.phate)
plot(so45.phate2)
##add the new DimReduc class to the seurat object
so45[["phate"]] <- CreateDimReducObject(embeddings = so45.phate$embedding * 100, key = "PHATE_", assay = DefaultAssay(so45))


################################################################################
######################## Lineage-DEGs in pstime heatmap

given <- sig.either


breaksList = unique(c(seq(-1,0,by=0.05), seq(0,4,by=0.2)))
breaksList = seq(-1,3, by=0.1)
color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList))

preSm <- predictSmooth(gam.j2, gene=given, nPoints=100, tidy=FALSE) #col others: others L2, other L4, J247 L2, J247 L4 (x100)
preSm.L2 <- preSm[rownames(sig.L2),c(201:300, 1:100)] # J247 vs. others (x100) for L2
preSm.L4 <- preSm[rownames(sig.L4),c(301:400, 101:200)] # J247 vs. others (x100) for L4

preSm.L2 <- t(scale(t(preSm.L2)))
preSm.L4 <- t(scale(t(preSm.L4)))

phmap.L2 <- pheatmap(preSm.L2, cluster_cols = FALSE, show_rownames = FALSE, 
                     show_colnames = FALSE, breaks = breaksList, color=color)#, cutree_rows = 7)
phmap.L4 <- pheatmap(preSm.L4, cluster_cols = FALSE, show_rownames = FALSE, 
                     show_colnames = FALSE, breaks = breaksList, color=color)#, cutree_rows = 7)
glist.L2 <- phmap.L2$tree_row$labels[phmap.L2$tree_row$order]
glist.L4 <- phmap.L4$tree_row$labels[phmap.L4$tree_row$order]
# [1] "HORVU6Hr1G038750" "HORVU7Hr1G091560" "HORVU7Hr1G000040" "HORVU2Hr1G109330" "HORVU0Hr1G003020"
#[6] "HORVU2Hr1G005110" "HORVU1Hr1G013140" "HORVU3Hr1G095200" "HORVU3Hr1G092400" "HORVU5Hr1G011780"
#[11] "HORVU7Hr1G089930"

#[1] "HORVU7Hr1G089930" "HORVU6Hr1G038750" "HORVU5Hr1G006930" "HORVU5Hr1G109720" "HORVU2Hr1G081720"
#[6] "HORVU7Hr1G106970" "HORVU2Hr1G117860" "HORVU2Hr1G066090" "HORVU7Hr1G000040" "HORVU5Hr1G076380"
#[11] "HORVU2Hr1G040860" "HORVU7Hr1G034630" "HORVU4Hr1G086050" "HORVU0Hr1G003020" "HORVU2Hr1G063800"
#[16] "HORVU1Hr1G013140" "HORVU2Hr1G005110" "HORVU7Hr1G008620" "HORVU3Hr1G092400" "HORVU2Hr1G109330"
#[21] "HORVU1Hr1G063740" "HORVU2Hr1G017370" "HORVU7Hr1G091560" "HORVU0Hr1G006890" "HORVU4Hr1G087980"
#[26] "HORVU4Hr1G084080" "HORVU3Hr1G074180" "HORVU5Hr1G112970" "HORVU5Hr1G018720" "HORVU5Hr1G011780"
#[31] "HORVU3Hr1G095200"


save.image(file='230320.RData')
load('230320.RData')


####MISC pstime heatmap for all four accs
#meta.cw <- meta[rownames(cw),]
loc.j0 <- which(meta[rownames(cw),"Accession"] == "J064")
loc.j6 <- which(meta[rownames(cw),"Accession"] == "J647")
loc.h6 <- which(meta[rownames(cw),"Accession"] == "H602")
cond.all <- cond.j2
cond.all[loc.j0] <- 3 
cond.all[loc.j6] <- 4 ##J247 to "2" , J064 to "3", J647 to "4", H602 to "1"



###MISC: pHeatmap of significant genes (FDR < 0.01)
compare_sig <- subset(compare, compare$fdr < 0.01)$Gene[1:500]
yhatSmooth <- predictSmooth(sce_sub, gene = compare_sig, nPoints = 100, tidy = FALSE)
yhatScaled <- t(scale(t(yhatSmooth)))
phmap <- pheatmap(yhatScaled, cluster_cols = FALSE, show_rownames = FALSE, 
                  show_colnames = FALSE, cutree_rows = 7, main="sigPattern 500")
cls <- as.data.frame(cutree(phmap$tree_row, k = 7))
cls$Gene <- rownames(cls) #no meaning, to keep df as df
gOrder <- compare_sig[phmap$tree_row$order]
cls <- cls[gOrder,]
table(cls[,1])[unique(cls[,1])]

####MISC: Gene expression line-plot for weeknumber & pseudotime
counts <- as.matrix(so@assays$RNA@counts)
counts <- subset(counts, rownames(counts) %in% VariableFeatures(so))
counts <- data.frame(t(counts))
#[1] 1940 28420
#HORVU0Hr1G000050 HORVU0Hr1G000160 HORVU0Hr1G000430 HORVU0Hr1G000620 
#10001_IH6_171215         1.492087        0.0000000        0.7649039         
#10002_IH6_171215         2.138252        0.4547237        0.0000000        
#10003_IH6_171215         2.295178        0.7962496        0.0000000        
#10004_IH6_171215         1.307644        0.0000000        0.3471242   
meta <- so@meta.data
pst <- meta[rownames(counts),c("Accession","week","LineageM")]

gid <- "HORVU0Hr1G003020"
pst$gid <- counts[,gid]
pst <- na.omit(pst)
#[1] 1913    4

### fit-pred function (simple multinomial regression)
lm_f <- function(exp){
  fit = lm(exp$gid~exp$LineageM + I(exp$LineageM^2))
  exp$pred <- predict(fit, newdata=exp)
  exp <- exp[order(exp$pred), ]
  return(exp)
}
pst.1 <- lm_f(subset(pst, pst$Accession == "J247"))
pst.2 <- lm_f(subset(pst, pst$Accession == "J064"))
pst.3 <- lm_f(subset(pst, pst$Accession == "J647"))
pst.4 <- lm_f(subset(pst, pst$Accession == "H602"))


#(pseudo-timeline)
plot(pst.1$gid~pst.1$LineageM, pch=19, cex=0.3, xlab="Pseudotime", ylab="RPM", 
     col="skyblue", ylim=c(0,max(pst$gid)), xlim=c(0,max(pst$LineageM)), main=gid)
points(pst.2$gid~pst.2$LineageM, pch=19, cex=0.3, col="#61D04F")
points(pst.3$gid~pst.3$LineageM, pch=19, cex=0.3, col="purple")
points(pst.4$gid~pst.4$LineageM, pch=19, cex=0.3, col="#DF536B")
lines(pst.1$LineageM, pst.1$pred, col="skyblue", lwd = 3)
lines(pst.2$LineageM, pst.2$pred, col="#61D04F", lwd = 3)
lines(pst.3$LineageM, pst.3$pred, col="purple", lwd = 3)
lines(pst.4$LineageM, pst.4$pred, col="#DF536B", lwd = 3)

#(weeknumber)
mean.1 <- tapply(pst.1$gid, pst.1$week, FUN=mean)
mean.2 <- tapply(pst.2$gid, pst.2$week, FUN=mean)
mean.3 <- tapply(pst.3$gid, pst.3$week, FUN=mean)
mean.4 <- tapply(pst.4$gid, pst.4$week, FUN=mean)

n.max <- max(pst$week) #match the lengths
length(mean.1) <- n.max
length(mean.2) <- n.max
length(mean.3) <- n.max
length(mean.4) <- n.max
mean.merge <- cbind(mean.1, mean.2, mean.3, mean.4)
p0 <- matplot(mean.merge, type ="b",  
              xlab="weeks", ylab="RPM", lwd=2, pch=19, 
              col=c("skyblue","#61D04F","purple","#DF536B"), lty=1)




