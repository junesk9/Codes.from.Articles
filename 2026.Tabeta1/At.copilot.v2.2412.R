############################################################
###### COPILOT scRNA-seq data analysis
###### 2024.10.29 JS KIM
##### the next steps
############################################################

library(Seurat)
library(CytoTRACE)
library(phateR)
library(ALRA)

library(pheatmap)

############################################################
####### Some Presets 
ggcolor <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
`%notin%` <- Negate(`%in%`)

setwd("/Users/junesk9/Library/CloudStorage/Box-Box/シングルセル解析(多部田)")
set.seed(101) #prevent the random output
############################################################
####### Data loading
############################################################
tbt1 <- readRDS("Tabeta-fugu5.seurat.v2.rds")
# An object of class Seurat 
# 160942 features across 8071 samples within 7 assays 
# Active assay: SCT (21377 features, 2000 variable features)
# 3 layers present: counts, data, scale.data
# 6 other assays present: RNA, spliced_RNA, unspliced_RNA, spliced_SCT, unspliced_SCT, ALRA
# 6 dimensional reductions calculated: pca, pca_whole, pca_vst2k, umap, umap_whole_all, umap_whole_pca14

############################################################
####### Data Transformation, clustering
############################################################
#(I. PCA for all features(genes))
AllFeatures <- rownames(tbt1@assays$RNA@data)
tbt1 <- FindVariableFeatures(tbt1, selection.method = "mvp")





############################################################
####### Cell type annotation - GSE213622 (Berrio 2022)
############################################################
ref <- readRDS("4_tissue.annotation/GSE213622_coaker_atpst_singlecell_sobj.rds")






############################################################
####### Cell type annotation - GSE226079
############################################################
dd <- readRDS("4_tissue.annotation/GSE226097/GSE226097_seedling_6d_230221.rds")
# An object of class Seurat 
# 24376 features across 41314 samples within 2 assays 
# Active assay: RNA (22376 features, 2000 variable features)
# 3 layers present: counts, data, scale.data
# 1 other assay present: integrated
# 2 dimensional reductions calculated: pca, umap

#(SCTransform as a pre-process)
options(future.globals.maxSize = 2048 * 1024 ^ 2)  #increase the memory pool
dd <- SCTransform(dd, vars.to.regress = c("percent.mt","percent.cp"))
saveRDS(dd, "4_tissue.annotation/GSE226097/GSE226097_seedling_6d_230221.rds")

table(dd@meta.data$CellType)
# Epidermal        Guard Meristematic    Mesophyll        Stele 
# 19233           89         2800        14952         4240 

#(random 100-cell select by celltype)
meta <- dd@meta.data
ref500 <- c()
for(ct in unique(meta$CellType)){
  cell.sub <- rownames(subset(meta, meta$CellType == ct))
  rnd <- min(100, length(cell.sub)) ## for GC as <100 cells
  cell100 <- sample(cell.sub, rnd)
  ref500 <- c(ref500, cell100)
}
meta500 <- meta[ref500, ] #[1] 489   9
sct500 <- as.data.frame(dd@assays$SCT@data[,ref500]) #[1] 22254   489
rna500 <- as.data.frame(dd@assays$RNA@data[,ref500]) #[1] 22376   489
#write.csv(meta500, "4_tissue.annotation/GSE226097/GSE226097_seedling_6d_230221.meta500.csv")
#write.csv(sct500, "4_tissue.annotation/GSE226097/GSE226097_seedling_6d_230221.sct500.csv")
#write.csv(rna500, "4_tissue.annotation/GSE226097/GSE226097_seedling_6d_230221.rna500.csv")

#(corr-test for best 5 matching)
sct.tbt <- tbt1@assays$SCT@data #[1] 21377  8071
sct.tbt <- as.data.frame(sct.tbt[rownames(sct.tbt) %in% rownames(sct500), ]) #[1] 19453  8071
sct500 <- sct500[rownames(sct.tbt), ] #[1] 19453   489 

corDF <- as.data.frame(matrix(nrow=length(colnames(sct.tbt)), ncol=4))
rownames(corDF) <- colnames(sct.tbt)
colnames(corDF) <- c("matchCells", "matchPCC","matchTypes", "bestAnnt")
for(idx in c(1:dim(sct.tbt)[2])){
  c = colnames(sct.tbt)[idx]
  tbtEXP <- as.vector(t(sct.tbt[,c]))
  cor.vec <- c()
  for(t in colnames(sct500)){
    refEXP <- as.vector(t(sct500[,t]))
    pcc <- cor.test(tbtEXP, refEXP, method="pearson")[c(3,4)]
    pcc <- pcc[2][[1]]
    cor.vec <- c(cor.vec, pcc)
  }
  names(cor.vec) <- colnames(sct500)
  cor.vec <- cor.vec[order(cor.vec, decreasing =TRUE)] #sorting
  best5 <- cor.vec[1:5]
  best.type <- meta500[names(best5),"CellType"]
  best <- names(table(best.type)[order(table(best.type), decreasing=TRUE)])[1]
  
  corDF[c,] <- c(paste0(names(best5),collapse=";"), paste0(best5,collapse=";"), paste0(best.type,collapse=";"), best)
}
corDF <- corDF[rownames(tbt1@meta.data), ] 
#write.csv(corDF, "4_tissue.annotation/GSE226097/GSE226097_seedling_6d_230221.corDF500.csv")
tbt1@meta.data$TLee2024_best_ident <- corDF$bestAnnt
 

