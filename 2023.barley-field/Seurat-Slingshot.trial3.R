#############################################################
########### Seurat & Slingshot analysis of Barley field mRNA-seq
######################### Junesk9 2022.12.08

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
#library(KEGGREST)

#library(gam) #for gam
#library(e1071) #for SVM
#library(earth) #for MARS
#library(mgcv) # for GAM, smoothing splines
#library(MuMIn) # for AIC calc


###############################################
#################################Prepare inputs
set.seed(101)
rpm_f = "../00_IPSR_KIBR_Sampledata/RPM_IPSR_KIBR_LB_17_19_BLD20201007.txt"
#rpm_f2 = "../00_IPSR_KIBR_Sampledata/barley-field.y19-20.morex1.v3.rpm.txt"
meta_f = "../00_IPSR_KIBR_Sampledata/samplesheet.txt"
hormon_f = "../00_IPSR_KIBR_Sampledata/Hormone_IPSR_KIBR_LB_17_19_BLD20201014.csv"
tf_f <- "../00_IPSR_KIBR_Sampledata/TranscriptionFactor_list_blastp_tophit_Gene_ID_rep_AtOsBd.txt"

####### Filter the metadata
meta <- read.table(meta_f, header=T, sep="\t")
#f_split1 <- function(s) strsplit(s,"-")[[1]][1] #extract "70338" from "70338-KJ6-200403"
#rownames(meta) <- sapply(meta$SampleID, f_split1)
rownames(meta) <- meta$SampleID
meta <- meta %>% filter(!Place %in% c("BLD"))

rpm <- read.table(rpm_f, header=T, row.names = 1)
colnames(rpm) <- str_sub(colnames(rpm), 2) ## remove the first "X" from the colnames, using stringr()
rpm <- rpm[,rownames(meta)] ## filter only columns with metadata

hormone <- read.table(hormon_f, row.names=1, header=T, sep=",")
colnames(hormone) <- str_sub(colnames(hormone), 2) 
hormone <- t(hormone)
meta <- merge(meta, hormone, by=0, all.x = TRUE)
rownames(meta) <- meta$Row.names

dim(rpm)
dim(meta)
#[1] 39734  1942
#[1] 1942   18

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
ddseq@meta.data$ABA <- meta2$ABA
ddseq@meta.data$IAA <- meta2$IAA
ddseq@meta.data$iP <- meta2$iP
ddseq@meta.data$JA <- meta2$JA
ddseq@meta.data$JAIle <- meta2$JAIle
ddseq@meta.data$SA <- meta2$SA
ddseq@meta.data$tZ <- meta2$tZ

f_split <- function(s) strsplit(s,"/")[[1]][2] #extract 12 from 2018/12/05
ddseq@meta.data$month <- sapply(meta2$Date, f_split)

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
#sum(sqrt(Stdev(ddseq, reduction="pca")[c(1,3,6,7,8,9)]))
#[1] 19.48567 #~20% of variations explained 
## umap
#ddseq <- RunUMAP(ddseq, dims=c(1,3,6,7,8,9))
#DimPlot(ddseq, reduction = "umap", group.by = c("Accession","Case"))
#FeaturePlot(ddseq, reduction = "umap", features = "Date.4")

############### GLR-aided selection of PCA-Dims
lasso.model.cv <- cv.glmnet(x = ddseq@reductions$pca@cell.embeddings, y = ddseq@meta.data$Date.4, family = "gaussian", alpha = 1)
plot(lasso.model.cv, xvar="lambda", label=TRUE)
ridge.model.cv <- cv.glmnet(x = ddseq@reductions$pca@cell.embeddings, y = ddseq@meta.data$Date.4, family = "gaussian", alpha = 0)
plot(ridge.model.cv, xvar="lambda", label=TRUE)

alpha <- seq(0.01, 0.99, 0.01)
mse.df <- NULL
for (i in 1:length(alpha)) {
  m <- cv.glmnet(x = ddseq@reductions$pca@cell.embeddings, y = ddseq@meta.data$Date.4, family = "gaussian", alpha = alpha[i])
  mse.df <- rbind(mse.df, data.frame(alpha = alpha[i], mse = min(m$cvm)))
}
best.alpha <- mse.df$alpha[mse.df$mse == min(mse.df$mse)]
#[1] 0.66
en.model.cv <- cv.glmnet(x = ddseq@reductions$pca@cell.embeddings, y = ddseq@meta.data$Date.4, family = "gaussian", alpha = best.alpha)
plot(en.model.cv, xvar="lambda", label=TRUE)
en.model <- glmnet(x = ddseq@reductions$pca@cell.embeddings, y = ddseq@meta.data$Date.4, family = "gaussian", lambda = en.model.cv$lambda.1se)
plot(abs(en.model$beta), main="betas of PC-dim, EN-model", pch=19)
abline(h=0.008,lty=2,col=2)
pcs <- which(abs(en.model$beta)>0.008)
#[1]  1  3 14 17 18 20 40
sum(sqrt(Stdev(ddseq, reduction="pca")[pcs]))
#16.18

### replace the PC data of Seurat object
ddseq2 <- ddseq
PC_table <- ddseq@reductions$pca@cell.embeddings
for(i in 1:length(pcs)){
  ddseq2@reductions$pca@cell.embeddings[,i] = PC_table[,pcs[i]]
}
ddseq2@reductions$pca@cell.embeddings = ddseq2@reductions$pca@cell.embeddings[,1:length(pcs)]
### Check UMAP clusters
ddseq2 <- RunUMAP(ddseq2, dims=c(1:length(pcs)))
DimPlot(ddseq2, reduction = "umap", group.by = c("Accession","Case","month"))
FeaturePlot(ddseq2, reduction = "umap", features = c("Date.4","ABA","IAA","SA","iP","tZ","JA","JAIle"))

umap <- ddseq2@reductions$umap@cell.embeddings
ddseq@meta.data$umap1 <- as.data.frame(umap)$UMAP_1
ddseq@meta.data$umap2 <- as.data.frame(umap)$UMAP_2


###############################################################
######################## Clustering 
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

#ddseq@meta.data$clsGLM_res0.5 = ddseq2@meta.data$RNA_snn_res.0.5
#ddseq@meta.data$clsGLM_res0.8 = ddseq2@meta.data$RNA_snn_res.0.8
#ddseq@meta.data$clsGLM_res1 = ddseq2@meta.data$RNA_snn_res.1
ddseq@meta.data$clsGLM_res1.2 = ddseq2@meta.data$RNA_snn_res.1.2
#ddseq@meta.data$clsGLM_res1.5 = ddseq2@meta.data$RNA_snn_res.1.5
#ddseq@meta.data$clsGLM_res2 = ddseq2@meta.data$RNA_snn_res.2

data = ddseq@meta.data
write.csv(data, "metadata.csv")
clus = ddseq2@meta.data$RNA_snn_res.1.2
Idents(ddseq2) <- ddseq2@meta.data$RNA_snn_res.1.2

# BoxPlot of Date.4 by clusters
oind <- order(as.numeric(by(data$Date.4, data$clsGLM_res1.2, median)))
data$clsGLM_res1.2 <- ordered(data$clsGLM_res1.2, levels=levels(data$clsGLM_res1.2)[oind])
boxplot(Date.4 ~ clsGLM_res1.2, data=data, main="Date.4 by clusters")
# Stacking barplot of Accs composition of clusters.
p <- ggplot(data, aes(x = clsGLM_res1.2, fill = Accession))
p + geom_bar() + theme_bw()

############### Slingshot analysis
pal <- c(RColorBrewer::brewer.pal(9, "Set1"), RColorBrewer::brewer.pal(8, "Set2"))

#re-draw umap-plot with plot()
plot(umap, col = pal[clus],  cex=.5,pch = 16)
for(i in levels(clus)){ 
  text( mean(umap[clus==i,1]),
        mean(umap[clus==i,2]), labels = i,font = 2) }

#Plot with estimated lineages
lineages <- getLineages(data=umap, clusterLabels=clus, start.clus="12", end.clus=c("13","7"))
SlingshotDataSet(lineages)
#Samples Dimensions
#1940          2

#lineages: 5 
#Lineage1: 12  6  4  14  11  10  0  2  16  
#Lineage2: 12  6  4  14  11  10  0  2  13  
#Lineage3: 12  6  4  14  11  10  0  2  7  
#Lineage4: 12  6  5  1  8  9  3  
#Lineage5: 12  6  4  15
plot(umap, col=pal[clus], pch=16)
lines(SlingshotDataSet(lineages), lwd = 3, col = 'black', show.constraints = TRUE)

#defining the principal curves (smooth trajectories)
curves <- getCurves(lineages, thresh = 0.01, stretch = 2, allow.breaks = TRUE, shrink = 0.99)
SlingshotDataSet(curves)
#-----
#curves: 5 
#Curve1: Length: 16.333	Samples: 854.63
#Curve2: Length: 17.094	Samples: 900.61
#Curve3: Length: 18.34	Samples: 948.4
#Curve4: Length: 20.404	Samples: 1124.93
#Curve5: Length: 9.45	Samples: 453.98

plot(umap, col=pal[clus], pch=16)
lines(SlingshotDataSet(curves), lwd = 3, col = 'black')
#or
lines(slingCurves(curves)[[1]], lwd = 3) #if you want to choose a curve to draw

pseudotime <- slingPseudotime(curves)[rownames(data),]
#[1] 1940    5
data <- cbind(data, pseudotime)
ddseq2@meta.data <- data[rownames(ddseq2@meta.data),]
FeaturePlot(ddseq2, reduction = "umap", features = "Lineage2")
FeaturePlot(ddseq2, reduction = "umap", features = "Lineage4")


write.csv(as.data.frame(ddseq2@meta.data),"metadata.csv")
############### Timeline analysis with tradeSeq/GAM
counts <- as.matrix(ddseq@assays$RNA@counts[ddseq@assays$RNA@var.features, ]) ##RPM matrix of variable genes
counts_norm <- GetAssayData(object = ddseq, slot = "scale.data") ### normalized RPM matrix
counts_norm <- subset(counts_norm, rownames(counts_norm) %in% VariableFeatures)
counts_filt <- counts[rowSums(counts > 50) > ncol(counts)/100, ] #just a smaller case
dim(counts) #dim(counts_norm)
#[1] 2842 1940
dim(counts_filt)
#[1] 1015 1940

##subset curve/counts for analyzing particular pseudotime data (Lineage2 & 4)
cw <- slingCurveWeights(curves)[,c(2,4)]
cw <- cw[rowSums(cw) > 0,] ##remove lines of weight = 0
#[1] 1865    2
pstime_sub <- slingPseudotime(SlingshotDataSet(curves), na=FALSE)[,c(2,4)] #need to fill numbers of NA, to waive the error during fitGAM
pstime_sub <- pstime_sub[rownames(cw),]
count_sub <- counts[,rownames(cw)]
clus_sub <- data[rownames(cw),]$clsGLM_res1.2
curve_sub <- SlingshotDataSet(curves[rownames(cw),c(2,4)])
curve_sub@curves <- curve_sub@curves[c(2,4)]
curve_sub@lineages <- curve_sub@lineages[c(2,4)]

#identify the best K for GAM
BPPARAM <- BiocParallel::bpparam() #for Intel-Mac 
BPPARAM$workers <- 8
BPPARAM <- BiocParallel::SnowParam(workers=8, type="SOCK", progressbar=TRUE) #for M1-Mac
register(BPPARAM)
icMat_sub <- evaluateK(counts = count_sub, pseudotime = pstime_sub, cellWeight=cw, 
                       k=3:15, nGenes=200, verbose=FALSE, plot=TRUE
                       ,parallel=TRUE, BPPARAM=BPPARAM)
#choose k as 10
sce_sub <- fitGAM(counts = count_sub, pseudotime = pstime_sub, cellWeight=cw,
                  nknots=10, parallel=TRUE, BPPARAM=BPPARAM) #few min w/ 8-threads 5-min w/ M1-single
plotGeneCount(curve_sub, count_sub, clusters = clus_sub, models = sce_sub) 

# save the curve pseudotime-linage
data <- as.data.frame(ddseq2@meta.data)
tmp <- merge(data, as.data.frame(sce_sub@colData$crv), by=0, all=TRUE)
rownames(tmp) <- tmp[,1]
data <- tmp[-1]
write.csv(data, "whole.metadata.csv")

#[Optional]for analyzing the whole lineages
icMat_whole <- evaluateK(counts = counts, sds = curves, k = 3:15, nGenes = 200, 
                   verbose = FALSE, plot = TRUE, parallel=TRUE, BPPARAM=BPPARAM)
sce <- fitGAM(counts = counts, sds = curves, nknots=10, parallel=TRUE, BPPARAM=BPPARAM) #few min w/ 8-threads
plotGeneCount(curves, counts, clusters = clus, models = sce) # confirm the fitGAM output

############## Identify genes on curves
plot_differential_expression <- function(feature_id) {
  #feature_id <- pseudotime_association %>% filter(pvalue < 0.05) %>% top_n(1, -waldStat) %>% pull(feature_id)
  title <- cowplot::ggdraw() + cowplot::draw_label(feature_id, fontface="bold", x=0.9, hjust=1) + ggplot2::theme(plot.margin=margin(0,0,0,7))
  genecount <- plotGeneCount(curve_sub, count_sub, gene = feature_id[1], clusters = clus_sub, models = sce_sub) + ggplot2::theme(legend.position = "none")
  smoothers <- plotSmoothers(sce_sub, as.matrix(count_sub), gene = feature_id[1])
  cowplot::plot_grid(title, cowplot::plot_grid(genecount, smoothers), ncol = 1, rel_heights = c(0.1, 1))
}

pseudotime_association <- associationTest(sce_sub, lineages=TRUE) ##test for all lineages simultaneously.
pseudotime_association$fdr_1 <- p.adjust(pseudotime_association$pvalue_1, method = "fdr")
pseudotime_association$fdr_2 <- p.adjust(pseudotime_association$pvalue_2, method = "fdr")
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

############### List & plot genes that no DE at the end, DE during lineages
library(ggplot2)
patternRes$Gene <- rownames(patternRes)
patternRes$pattern <- patternRes$waldStat
patternRes$log2fc <- log(patternRes$fcMedian,2)
pRes <- patternRes[, c("Gene", "pattern","fdr")]
endRes$Gene <- rownames(endRes)
endRes$end <- endRes$waldStat
eRes <- endRes[, c("Gene", "end")]

compare <- merge(pRes, eRes, by = "Gene", all = FALSE)
compare$transientScore <- 
  rank(-compare$end, ties.method = "min")^2 + rank(compare$pattern, ties.method = "random")^2
compare <- compare[order(compare$transientScore, decreasing=TRUE),]
#[1] 2842    5

#plot DE-DE
ggplot(compare, aes(x = log(pattern), y = log(end))) +
  geom_point(aes(col = transientScore)) +
  labs(x = "patternTest Wald Statistic (log scale)",
       y = "diffEndTest Wald Statistic (log scale)") +
  scale_color_continuous(low = "yellow", high = "red") +
  theme_classic()

# pHeatmap of significant genes (FDR < 0.01)
compare_sig <- subset(compare, compare$fdr < 0.01)$Gene[1:500]
yhatSmooth <- predictSmooth(sce_sub, gene = compare_sig, nPoints = 100, tidy = FALSE)
yhatScaled <- t(scale(t(yhatSmooth)))
phmap <- pheatmap(yhatScaled, cluster_cols = FALSE, show_rownames = FALSE, 
                  show_colnames = FALSE, cutree_rows = 7, main="sigPattern 500")
cls <- as.data.frame(cutree(phmap$tree_row, k = 7))
cls$Gene <- rownames(cls) #not meaning, to keep df as df
gOrder <- compare_sig[phmap$tree_row$order]
cls <- cls[gOrder,]
table(cls[,1])[unique(cls[,1])]
#7   5   6   4   2   3   1 
#12  24  47  84 232  46  55 

# (Omake) pHeatmap of flor145
flor145 <- read.table("../30_Seruat.leaf/2_FLOR145-L-plots/FLOR_genes.txt",header=T, row.names=1,sep="\t")
flor145 <- rownames(flor145)
counts <- as.matrix(ddseq@assays$RNA@counts) ##RPM matrix of whole genes
counts_norm <- GetAssayData(object = ddseq, slot = "scale.data") ### normalized RPM matrix

L2_order <- order(data$Lineage2, na.last=NA) ##new pseudotime-order; remove the pseudo-pseudo time from [pstime_sub]
L4_order <- order(data$Lineage4, na.last=NA)
L2_count = count_flor[,L2_order]
L4_count = count_flor[,L4_order]
L24_count <- cbind(L4_count[,ncol(L4_count):1], L2_count)
L24_clus <- data[colnames(L24_count),]$clsGLM_res1.2
flor_h <- heatmap(log1p(L24_count), Colv = NA, ColSideColors = pal[L24_clus], scale="none", keep.dendro=TRUE)
flor_hclus <- as.hclust(flor_h$Rowv)
#plot(flor_hclus)
#cutree(flor_hclus, k=3)

FeaturePlot(ddseq2, reduction = "umap", features = "HORVU0Hr1G003020")
FeaturePlot(ddseq2, reduction = "umap", features = "HORVU2Hr1G063800")
FeaturePlot(ddseq2, reduction = "umap", features = "HORVU7Hr1G024610")
FeaturePlot(ddseq2, reduction = "umap", features = "HORVU1Hr1G011800")
FeaturePlot(ddseq2, reduction = "umap", features = "HORVU3Hr1G010240")
FeaturePlot(ddseq2, reduction = "umap", features = "HORVU4Hr1G077850")
FeaturePlot(ddseq2, reduction = "umap", features = "HORVU6Hr1G058740")
FeaturePlot(ddseq2, reduction = "umap", features = "HORVU4Hr1G039300")
FeaturePlot(ddseq2, reduction = "umap", features = "HORVU2Hr1G023180")




######################################################
############## GSEA
#### Prepare GOSlim DB
go_f ="../00_IPSR_KIBR_Sampledata/Barley.IBSCv2.ensembl49.GOSlim.txt"
go <- read.table(go_f, header=T, sep="\t")
go$evi = "IEA"
go <- go[,c(2,3,1)]
goFrame <- GOFrame(go, organism="barley")
goAllFrame <- GOAllFrame(goFrame)
gsc <- GeneSetCollection(goAllFrame, setType = GOCollection())
all <- as.vector(go[,3])  ## a vector of all referenced list

### set function to GO-GSEA
go_gsea <- function(gid_ls, prefix){
  for (f in c("MF","CC","BP")){
    outf = paste(prefix,f,sep=".")
    outf = paste(outf,"txt",sep=".")
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
    result <- hyperGTest(p)
    
    ### Check and output the result
    print(head(summary(result))) # Check the result
    write.table(summary(result), file = outf, quote=FALSE,sep="\t",row.names=FALSE)
  }
} 

#### Prepare input gene lists
traj2_sig <- rownames(subset(cls, cls[,1] %in% c(5,6,4)))
traj4_sig <- rownames(subset(cls, cls[,1] %in% c(2,3,1)))

go_gsea(traj2_sig, "GO.traj2")
go_gsea(traj4_sig, "GO.traj4")

### Visualize as a wordcloud (not complete)
traj2_bp <- summary(result)
traj2_bp <- subset(traj2_bp, traj2_bp$Pvalue < 0.05)


###############################################################################
####################### Cluster analysis
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

######### Identify genes specific for clusters
idents_ordered <- ordered(data$clsGLM_res1.2, levels=levels(data$clsGLM_res1.2)[oind])
levels(ddseq2) <- levels(idents_ordered)
##(from https://hackmd.io/@takahirosuzuki/ByoEBNVPL)
integrated.markers <- FindAllMarkers(object=ddseq2, only.pos=FALSE, min.pct=0.25, test.use="MAST") #taking time
write.csv(integrated.markers,"integrated.markers.csv")
integrated.markers %>%
  filter(avg_log2FC >= 0) %>%
  group_by(cluster) %>%
  top_n(5, avg_log2FC) -> top5.posi
plot <- DoHeatmap(object = ddseq2, features=top5.posi$gene, size=3, label=TRUE) + NoLegend() + NoAxes()
##(see https://satijalab.org/seurat/articles/essential_commands.html)
  
## for multi-cluster markers
markers5.6 <- FindMarkers(ddseq2, ident.1 = c(5,6), min.pct = 0.25, only.pos = FALSE, test.use="MAST")
markers8.9 <- FindMarkers(ddseq2, ident.1 = c(8,9), min.pct = 0.25, only.pos = FALSE, test.use="MAST")
markers14 <- FindMarkers(ddseq2, ident.1 = 14, min.pct = 0.25, only.pos = FALSE, test.use="MAST")
markers15 <- FindMarkers(ddseq2, ident.1 = 15, min.pct = 0.25, only.pos = FALSE, test.use="MAST")
markers7n13 <- FindMarkers(ddseq2, ident.1 = 7, ident.2 = 13, min.pct = 0.25, only.pos = FALSE, test.use="MAST")

f_sig <- function(df) {tmp <- subset(df, df$p_val_adj < 0.05 & df$avg_log2FC >= log2(1.5))
                       tmp <- tmp[order(tmp$avg_log2FC, decreasing=TRUE),]
                       gid <- rownames(tmp)
                       return(gid)} #extract g_id, ordered by log2FC
gid_5.6 <- f_sig(markers5.6)
#[1] 400
gid_8.9 <- f_sig(markers8.9)
#[1] 1940
gid_14 <- f_sig(markers14)
#[1] 111
gid_15 <- f_sig(markers15)
#[1] 551
gid_7n13 <- f_sig(markers7n13)
#[1] 1161

## pHeatmap for each geneset-in-interest
data <- data[order(data$Date.2), ]
J247db <- subset(data, data$Accession == "J247")
H602db <- subset(data, data$Accession == "H602")
cnt_5.6 <- rpm[gid_5.6, c(rownames(H602db), rownames(J247db))][1:200,]
scale_5.6 <- t(scale(t(cnt_5.6)))

breaksList = seq(-1,3,by=0.1)
color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList))
ph_5.6 <- pheatmap(scale_5.6, cluster_cols = FALSE, show_rownames = FALSE, color=color,
                   show_colnames = FALSE, main="scale_5.6", breaks = breaksList, cutree_rows = 9)
rOrder <- ph_5.6$tree_row$order
gOrder <- compare_sig[phmap$tree_row$order]
cls <- cls[gOrder,]
table(cls[,1])[unique(cls[,1])]

## pHeatmap with averages.
f_MeanByWeek <- function(rpm_df, data_df){
  for (i in unique(data_df$Date.2)){
    i_points <- rownames(subset(data_df, data_df$Date.2 == i))
    rpm_sub <- rpm_df[, i_points]
    rpm_sub <- na.omit(rpm_sub)
    #print(head(rpm_sub, c(5,5)))
    rpm_df$V1 <- rowMeans(rpm_sub)
    colnames(rpm_df)[length(rpm_df)] <- i 
  }
  d_len <- length(unique(data_df$Date.2)) 
  st <- length(rpm_df) - d_len + 1
  ed <- length(rpm_df)
  #print(c(st, ed))
  rpm_df2 <- rpm_df[,st:ed]
  return(rpm_df2)
} ## Mean DF by the same week number

## for cluster 5 & 6
J247rpm <- rpm[gid_5.6, rownames(J247db)][1:300,]
H602rpm <- rpm[gid_5.6, rownames(H602db)][1:300,]
J247mean <- f_MeanByWeek(J247rpm, J247db)
H602mean <- f_MeanByWeek(H602rpm, H602db)
HJmean <- cbind(as.matrix(H602mean), as.matrix(J247mean))
HJscale <- t(scale(t(HJmean)))
HJmean_ph <- pheatmap(HJscale, cluster_cols = FALSE, show_rownames = FALSE, color=color,
                   show_colnames = FALSE, main="cls5.6 H602:J247 top300", breaks = breaksList, cutree_rows = 6)
#extract the cluster-gene id info
cls <- as.data.frame(cutree(HJmean_ph$tree_row, k = 6))
cls$Gene <- rownames(cls) #not meaning, to keep df as df
rOrder <- HJmean_ph$tree_row$order
gOrder <- rownames(J247rpm)[rOrder]
cls <- cls[gOrder,]
table(cls[,1])[unique(cls[,1])]

gids_56 <- rownames(subset(cls, cls[,1] %in% c(4,1,3)))
go_gsea(gids_56, "GO.clus56")

## for cluster 8 & 9
J247rpm <- rpm[gid_8.9, rownames(J247db)][1:300,]
H602rpm <- rpm[gid_8.9, rownames(H602db)][1:300,]
J247mean <- f_MeanByWeek(J247rpm, J247db)
H602mean <- f_MeanByWeek(H602rpm, H602db)
HJmean <- cbind(as.matrix(H602mean), as.matrix(J247mean))
HJscale <- t(scale(t(HJmean)))
HJmean_ph <- pheatmap(HJscale, cluster_cols = FALSE, show_rownames = FALSE, color=color,
                      show_colnames = FALSE, main="cls8.9 H602:J247 top300", breaks = breaksList, cutree_rows = 6)
#extract the cluster-gene id info
cls <- as.data.frame(cutree(HJmean_ph$tree_row, k = 6))
cls$Gene <- rownames(cls) #not meaning, to keep df as df
rOrder <- HJmean_ph$tree_row$order
gOrder <- rownames(J247rpm)[rOrder]
cls <- cls[gOrder,]
table(cls[,1])[unique(cls[,1])]

gids_89 <- rownames(subset(cls, cls[,1] %in% c(5,1,2)))
gids_89i <- rownames(subset(cls, cls[,1] %in% c(6,3,4)))
go_gsea(gids_89, "GO.clus89")
go_gsea(gids_89i, "GO.clus89-2")

## for cluster 14
J247rpm <- rpm[gid_14, rownames(J247db)]
H602rpm <- rpm[gid_14, rownames(H602db)]
J247mean <- f_MeanByWeek(J247rpm, J247db)
H602mean <- f_MeanByWeek(H602rpm, H602db)
HJmean <- cbind(as.matrix(H602mean), as.matrix(J247mean))
HJscale <- t(scale(t(HJmean)))
HJmean_ph <- pheatmap(HJscale, cluster_cols = FALSE, show_rownames = FALSE, color=color,
                      show_colnames = FALSE, main="cls14 H602:J247 n=111", breaks = breaksList, cutree_rows = 6)
#extract the cluster-gene id info
cls <- as.data.frame(cutree(HJmean_ph$tree_row, k = 6))
cls$Gene <- rownames(cls) #not meaning, to keep df as df
rOrder <- HJmean_ph$tree_row$order
gOrder <- rownames(J247rpm)[rOrder]
cls <- cls[gOrder,]
table(cls[,1])[unique(cls[,1])]

gids_14 <- rownames(subset(cls, cls[,1] %in% c(2,4)))
go_gsea(gids_14, "GO.clus1")


#(Omake) TF_list
tf <- read.table(tf_f, header=T, sep="\t")
tf_ls <- tf$query_acc.ver
J247rpm <- rpm[tf_ls, rownames(J247db)][1:500,]
H602rpm <- rpm[tf_ls, rownames(H602db)][1:500,]
J247mean <- f_MeanByWeek(J247rpm, J247db)
H602mean <- f_MeanByWeek(H602rpm, H602db)
HJmean <- cbind(as.matrix(H602mean), as.matrix(J247mean))
HJmean <- HJmean[rowSums(HJmean)>0,] # remove rowSums>0; to evade pheatmap error
#[1] 454   41
HJscale <- t(scale(t(HJmean)))
HJmean_ph <- pheatmap(HJscale, cluster_cols = FALSE, show_rownames = FALSE, color=color,
                      show_colnames = FALSE, main="TFs n=500", breaks = breaksList, cutree_rows = 10)


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



