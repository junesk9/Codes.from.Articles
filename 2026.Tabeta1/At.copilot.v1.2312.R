############################################################
###### COPILOT scRNA-seq data analysis
###### 2023.12.14 JS KIM
############################################################
#https://github.com/Hsu-Che-Wei/COPILOT

##### Install packages
#BiocManager::install("DropletUtils")
#install.packages("remotes")
#remotes::install_github("chris-mcginnis-ucsf/DoubletFinder")
#install.packages("devtools")
#devtools::install_github('Hsu-Che-Wei/COPILOT')

#for COPILOT
library(sf) ## for Oct.2023 issue.
library(COPILOT)
library(Seurat) ### Using an legacy version v4.4 for compatibilty with "DoubletFinder_v3"
library(tidyverse)

#For Seurat analysis
library(MAST)
library(CytoTRACE)
library(glmnet)
library(slingshot)
library(tradeSeq)
library(pheatmap)
library(phateR)
library(ALRA)

#For DEG analysis
library(clusterProfiler)
library(org.At.tair.db) #https://bi.biopapyrus.jp/rnaseq/annotation/org.at.tair.db.html
library(KEGGREST)




setwd("/Users/junesk9/Library/CloudStorage/Box-Box/シングルセル解析(多部田)/1_scKB-COPILOT")
set.seed(101) #prevent the random output

ggcolor <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
###############################################
#################### COPILOT filtering 23.12.14
################################################

##load unwanted genes [mandatory option while not applying this snRNA-seq, try anyway PMID[30913408]]
pp.genes <- as.character(read.table("../0_参考になるもの/Protoplasting_DEgene_FC2_list.txt", header=F)$V1)
#[1] "AT5G58860" "AT1G55990" "AT1G26390" "AT5G23190" "AT5G04950"

#>WT
sample.stats <- data.frame(stat = c('Sample','Name','Source','Genotype','Transgene','Treatment','Age','Timepoint','Rep','Target Cells','Date','Seq Run'), 
                           value = c("WT","WT","NA","WT","NA","untreated","8das","NA","NA","8,000","NA","NA"))
copilot(sample.name = "WT_8DAS", species.name = "Arabidopsis thaliana", transcriptome.name = "TAIR10", sample.stats = sample.stats, mt.pattern = "ATMG", 
        mt.threshold = 5, cp.pattern = "ATCG", remove.doublet = TRUE, do.seurat = TRUE, do.annotation = FALSE, unwanted.genes = pp.genes, 
        dir_to_color_scheme = "../0_参考になるもの/supp_data/color_scheme_at.RData",  min.UMI.low.quality = 100, min.UMI.high.quality = 300)

#>fugu5
sample.stats <- data.frame(stat = c('Sample','Name','Source','Genotype','Transgene','Treatment','Age','Timepoint','Rep','Target Cells','Date','Seq Run'), 
                           value = c("fugu5","fugu5","NA","fugu5","NA","untreated","8das","NA","NA","8,000","NA","NA"))
copilot(sample.name = "fugu5_8DAS", species.name = "Arabidopsis thaliana", transcriptome.name = "TAIR10", sample.stats = sample.stats, mt.pattern = "ATMG", 
        mt.threshold = 5, cp.pattern = "ATCG", remove.doublet = TRUE, do.seurat = TRUE, do.annotation = FALSE, unwanted.genes = pp.genes, 
        dir_to_color_scheme = "../0_参考になるもの/supp_data/color_scheme_at.RData",  min.UMI.low.quality = 100, min.UMI.high.quality = 300)


#########################################################
######################## Seurat 
##########################################################
setwd("/Users/junesk9/Library/CloudStorage/Box-Box/シングルセル解析(多部田)/")

##1. merge into one DB
wt <- readRDS("./WT_8DAS/WT_8DAS_COPILOT.rds")
#An object of class Seurat 
#118594 features across 2700 samples within 6 assays 
#Active assay: SCT (18291 features, 18291 variable features)
#3 layers present: counts, data, scale.data
#5 other assays present: RNA, spliced_RNA, unspliced_RNA, spliced_SCT, unspliced_SCT
#2 dimensional reductions calculated: pca, umap
fg5 <- readRDS("./fugu5_8DAS/fugu5_8DAS_COPILOT.rds")
#An object of class Seurat 
#132189 features across 5371 samples within 6 assays 
#Active assay: SCT (20702 features, 20702 variable features)
#3 layers present: counts, data, scale.data
#5 other assays present: RNA, spliced_RNA, unspliced_RNA, spliced_SCT, unspliced_SCT
#2 dimensional reductions calculated: pca, umap

#(merge)
tbt1 <- merge(wt, y=fg5, add.cell.ids=c("WT","fugu5"), project="tbt_all")
Idents(tbt1) <- tbt1$orig.ident
#(re-normalize the count data by SCTransform)
tbt1 <- SCTransform(tbt1, variable.features.n = nrow(tbt1), assay = "RNA", new.assay.name = "SCT", verbose = FALSE) 
tbt1
#An object of class Seurat 
#133843 features across 8071 samples within 6 assays 
#Active assay: SCT (20790 features, 0 variable features)
#3 layers present: counts, data, scale.data
#5 other assays present: RNA, spliced_RNA, unspliced_RNA, spliced_SCT, unspliced_SCT

# check the data integrity
plot(sort(tbt1@meta.data$nCount_RNA)) #  Total read count
plot(sort(tbt1@meta.data$nFeature_RNA)) # Number of the detected genes
saveRDS(tbt1, "./Tabeta-fugu5.seurat.rds")


####2. Re-run PCA & UMAP for multiple cases
#(I. PCA for all features(genes))
AllFeatures <- rownames(tbt1@assays$RNA@data)
tbt1 <- RunPCA(tbt1, vervose=FALSE, approx=FALSE, npcs=50, features=AllFeatures)
ElbowPlot(tbt1, ndims=50)
tbt1@reductions$pca_whole <- tbt1@reductions$pca #save the feature
#(II. PCA for top-2k dynamic genes)
tbt1 <- FindVariableFeatures(tbt1, selection.method = "mvp")#, dispersion.cutoff = c(0, Inf),mean.cutoff = c(0.1, 8)) #default as "vst"
VariableFeaturePlot(tbt1, selection.method="mvp")
VariableFeatures <- VariableFeatures(tbt1)
tbt1 <- RunPCA(tbt1, vervose=FALSE, approx=FALSE, npcs=50, features=VariableFeatures)
ElbowPlot(tbt1, ndims=50)
tbt1@reductions$pca_vst2k <- tbt1@reductions$pca #save the feature

#(III. Check eigen values for each case)
#currently, the seurat command Stdev() not working properly, need to do manually
sum(sqrt(sapply(as.data.frame(tbt1@reductions$pca_whole@cell.embeddings), sd)))
#[1] 100.2589
sum(sqrt(sapply(as.data.frame(tbt1@reductions$pca_whole@cell.embeddings), sd))[1:14])
#[1] 34.96529
sum(sqrt(sapply(as.data.frame(tbt1@reductions$pca_vst2k@cell.embeddings), sd)))
#[1] 85.378
sum(sqrt(sapply(as.data.frame(tbt1@reductions$pca_vst2k@cell.embeddings), sd))[1:13])
#[1] 28.29181

#(IV. UMAP for each case)
tbt1@reductions$pca <- tbt1@reductions$pca_vst2k
tbt1 <- RunUMAP(tbt1, dims=c(1:50))
DimPlot(tbt1, reduction="umap", group.by="orig.ident") + ggtitle("vst2k-pca50") 
tbt1@reductions$umap_vsk2k_all <- tbt1@reductions$umap
tbt1 <- RunUMAP(tbt1, dims=c(1:13))
DimPlot(tbt1, reduction="umap", group.by="orig.ident") + ggtitle("vst2k-pca13") 
tbt1@reductions$umap_vsk2k_pca13 <- tbt1@reductions$umap

tbt1@reductions$pca <- tbt1@reductions$pca_whole
tbt1 <- RunUMAP(tbt1, dims=c(1:50))
DimPlot(tbt1, reduction="umap", group.by="orig.ident") + ggtitle("whole-pca50")
tbt1@reductions$umap_whole_all <- tbt1@reductions$umap
tbt1 <- RunUMAP(tbt1, dims=c(1:14))
DimPlot(tbt1, reduction="umap", group.by="orig.ident") + ggtitle("whole-pca14")
tbt1@reductions$umap_whole_pca14 <- tbt1@reductions$umap

saveRDS(tbt1, "./Tabeta-fugu5.seurat.rds")


### 3.Clustering
tbt1 = FindNeighbors(tbt1, dims=1:14) ## for the case of whole-pc14
tbt1 = FindClusters(tbt1,resolution = 0.5) #13
#DimPlot(tbt1, reduction="umap") + ggtitle("whole-pca14; res=0.5")
#DimPlot(tbt1, reduction="pca")
#tbt1$SCT_snn_res.0.5

tbt1 = FindClusters(tbt1,resolution = 1.0) #16
#DimPlot(tbt1, reduction="umap", group.by="SCT_snn_res.1") + ggtitle("whole-pca14; res=1.0")
#DimPlot(tbt1, reduction="pca")
#tbt1$SCT_snn_res.1

tbt1 = FindClusters(tbt1,resolution = 1.5) #18
#DimPlot(tbt1, reduction="umap", group.by="SCT_snn_res.1.5") + ggtitle("whole-pca14; res=1.5")
#DimPlot(tbt1, reduction="pca")
#tbt1$SCT_snn_res.1.5

tbt1 = FindClusters(tbt1,resolution = 2.0) #22
#DimPlot(tbt1, reduction="umap", group.by="SCT_snn_res.2") + ggtitle("whole-pca14; res=2.0")
#DimPlot(tbt1, reduction="pca")
#tbt1$SCT_snn_res.2

# Res=0.5 looks good.
Idents(tbt1) <- tbt1$SCT_snn_res.0.5

###4. CytoTRACE for cellular latent time estimation
results <- CytoTRACE(as.matrix(tbt1@assays$RNA@counts), ncores = 6, subsamplesize = 1000)
tbt1$CytoTRACE <- 1 - results$CytoTRACE

#Cluster ordering by the latent time
#(ordering by median; set levels of factor())
oind <- order(as.numeric(by(tbt1$CytoTRACE, tbt1$SCT_snn_res.0.5, median)))
# [1] 13 12  3  9  8  2  5  6 11  7  1 10  4
tbt1$SCT_snn_res.0.5 <- factor(tbt1$SCT_snn_res.0.5, levels=(oind-1)) # as the cluster starts from "0".
Idents(tbt1) <- tbt1$SCT_snn_res.0.5

saveRDS(tbt1, "../Tabeta-fugu5.seurat.rds")

#(visualize)
VlnPlot(tbt1, features="CytoTRACE")
boxplot(tbt1$CytoTRACE ~ tbt1$SCT_snn_res.0.5, main="CytoTRACE by clusters (res=0.5)")

###5. find-out marker genes (w/ MAST) by the given cluster order
integrated.markers <- FindAllMarkers(object=tbt1, only.pos=TRUE, min.pct=0.25, test.use="MAST")
write.csv(integrated.markers,"./3_Marker-genes/cluster.markers.MAST.csv")
integrated.markers %>%
  filter(avg_log2FC >= 0) %>%
  group_by(cluster) %>%
  top_n(5, avg_log2FC) -> top5.posi
#(visualize)
DoHeatmap(object = ddseq2, features=top5.posi$gene, size=3, label=TRUE) + NoLegend() + NoAxes()

#(ClusterCompare: KEGG, GO, and so on)
#collect DEG gids by cluster
deg_list <- split(integrated.markers$gene, integrated.markers$cluster)
names(deg_list) <- levels(Idents(tbt1)) 
tair_list <- deg_list #need to keep for KEGG analysis
for(i in names(deg_list)){
  deg_list[[i]] <- bitr(deg_list[[i]],
                        fromType="TAIR",
                        toType="ENTREZID",
                        OrgDb="org.At.tair.db",
                        drop=TRUE) #drop for remove the unconverted gid, keep them with "NA" by "FALSE"
  deg_list[[i]] <- deg_list[[i]]$ENTREZID 
} #deg list including ENTREZID
#(run compareCluster())
GOBP <- compareCluster(geneClusters=deg_list, fun="enrichGO", ont="BP", OrgDb = "org.At.tair.db")
GOBP2 <- pairwise_termsim(GOBP)
dotplot(GOBP, group=TRUE)
emapplot(GOBP2)
#GOBP3 <- setReadable(GOBP, OrgDb="org.At.tair.db", keyType = "ENTREZID")
#(prepare output: convert ENTREZID to TAIR)
goout <- as.data.frame(GOBP@compareClusterResult)
Ent2Tair <- function(vec){
  ent = strsplit(vec, "/")[[1]]
  tair = suppressMessages(bitr(ent, fromType="ENTREZID", toType="TAIR", OrgDb="org.At.tair.db"))[,2]
  tair = paste(tair, collapse="/")
  return(tair)}
goout$geneID <- sapply(goout$geneID, Ent2Tair)
goout$GeneRatio <- sapply(goout$GeneRatio, function(x)paste0(" ",x)) #prevent the automatic date-conversion in Excel.
goout$BgRatio <- sapply(goout$BgRatio, function(x)paste0(" ",x))
write.csv(goout,"./3_Marker-genes/Cluster-marekre.GOBP.csv", quote=TRUE, row.names=FALSE)


KEGG <- compareCluster(geneClusters=tair_list, fun = "enrichKEGG", organism = "ath", pvalueCutoff=0.1)
KEGG2 <- pairwise_termsim(KEGG)
dotplot(KEGG, group=TRUE)
emapplot(KEGG2)

WP <- compareCluster(geneClusters=deg_list, fun = "enrichWP", organism = "Arabidopsis thaliana", readable=TRUE)
dotplot(WP)


#########################################################
######################## AUCell to label cell identity
############ (https://bioconductor.org/books/3.13/OSCA.basic/cell-type-annotation.html#assigning-cell-labels-from-gene-sets)
##########################################################
#(need to unload the Seurat package, so  extract the EXP matrix first)
exp <- tbt1@assays$SCT@data
alra <- tbt1@assays$ALRA@data
#(unload the clushing packages)
#unloadNamespace("sctransform") ##need to unl
#unloadNamespace("Seurat")
#(load necessary libraries)
library(GSEABase)
library(AUCell)
library(scRNAseq) ##for assay()
library(dplyr) ##for distinct()
# library(scran)
# library(scater)
# sce.zeisel <- ZeiselBrainData()
# sce.tasic <- TasicBrainData()

mgenes <- read.table("../6_AUCell_PlantscRNAdb/ara_mkgs.txt", header=F,sep="\t")
mgenes <- subset(mgenes, mgenes[,3] == "Leaf")
#mg_rev <- as.data.frame(t(rev(as.data.frame(t(mgenes)))))
#[1] 28724     4
mgenes_uniq <- distinct(mgenes, V1, .keep_all = TRUE) ##retain uniq gene_ids
#[1] 8845    4
# V1                   V2        V3                      V4
# 1 AT1G34245 Arabidopsis thaliana Cotyledon Early-stage meristemoid
# 2 AT3G06120 Arabidopsis thaliana Cotyledon Early-stage meristemoid
# 3 AT5G53210 Arabidopsis thaliana Cotyledon Early-stage meristemoid
# 4 AT5G60880 Arabidopsis thaliana Cotyledon Early-stage meristemoid
# 5 AT1G80080 Arabidopsis thaliana Cotyledon              Guard cell
# 6 AT2G46720 Arabidopsis thaliana Cotyledon              Guard cell
mgr_uniq <- distinct(mg_rev, V1, .keep_all = TRUE) ##retain uniq gene_ids
mg_list <- split(mgenes_uniq[,1], mgenes_uniq[,4])
lengths(mg_list)
# Early-stage meristem Early-stage meristemoid              Guard cell       Guard mother cell  Late-stage meristemoid 
# 1                    3421                     783                    1105                     693 
# Meristem mother cell Meristemoid mother cell          Mesophyll cell           Pavement cell       Reprogrammed cell 
# 1                     162                     887                     717                     244 
# Vascular stem cell        Young guard cell 
# 824                       7 

mg1 <- read.table("../6_AUCell_PlantscRNAdb/Berrio2022_marker1.txt", header=T)

mg_list <- split(mg1[,2], mg1[,1])
marker.set <- lapply(names(mg_list), function(x) {GeneSet(mg_list[[x]], setName=x)})
marker.set <- GeneSetCollection(marker.set)
# GeneSetCollection
# names: Early-stage meristem, Early-stage meristemoid, ..., Young guard cell (12 total)
# unique identifiers: AT1G18710, AT1G34245, ..., AT5G28540 (8845 total)
# types in collection:
#   geneIdType: NullIdentifier (1 total)
# collectionType: NullCollection (1 total)

rankings <- AUCell_buildRankings(as.matrix(exp), plotStats=FALSE, verbose=FALSE)
cell.aucs <- AUCell_calcAUC(marker.set, rankings)
results <- as.data.frame(t(assay(cell.aucs)))
results$best1 = colnames(results)[max.col(results)]
results <- results[rownames(tbt1@meta.data), ]
head(results)
# Epidermis    Mesophyll Vasculature     best1
# WT_AAACGAAGTCCGGTGT 0.303765122 4.280022e-02  0.09920075 Epidermis
# WT_AAACGAAGTCGTGGAA 0.027783064 6.802571e-02  0.00000000 Mesophyll
# WT_AAACGCTGTCTCCTGT 0.014489602 4.656449e-02  0.00000000 Mesophyll
# WT_AAAGAACCAAGCTACT 0.045969825 7.332991e-05  0.00000000 Epidermis
# WT_AAAGAACGTGATCATC 0.004050564 1.047396e-01  0.00000000 Mesophyll



par(mfrow=c(3,4))
AUCell_exploreThresholds(cell.aucs, plotHist=TRUE, assign=TRUE)


#################################################
############## Optional attempt for compareClster
#(run a custum compareCluster() with "enricher", using KEGGREST())
#(extract the data & brush-up)
#https://bi.biopapyrus.jp/rnaseq/annotation/keggrest.html
#kegg2tair <- keggLink("ath", "pathway")
#keggid <- gsub("path:", "", names(kegg2tair))
#keggid.uniq <- unique(keggid)
#KEGG2TAIR <- vector("list", length(keggid.uniq))
#names(KEGG2TAIR) <- keggid.uniq
#for(i in 1:length(KEGG2TAIR)) {
#  KEGG2TAIR[[i]] <- as.character(gsub("ath:", "", kegg2tair[keggid == keggid.uniq[i]]))
#}
kegg2tair <- keggLink("ath", "pathway")
kegg <- sapply(names(kegg2tair), function(x){strsplit(x, ":")[[1]][2]})
tair <- sapply(kegg2tair, function(x){strsplit(x, ":")[[1]][2]})
KEGG2TAIR <- data.frame(kegg,tair)
#[1] 12640     2

kegg2tair <- keggList("pathway", "ath")
name <- sapply(kegg2tair, function(x){strsplit(x, "- Arabidopsis", fixed=T)[[1]][1]})
KEGG2NAME <- data.frame(name)
KEGG2NAME$kegg <- rownames(KEGG2NAME)
#[1] 143     2
KEGG2NAME <- merge(KEGG2TAIR, KEGG2NAME, by="kegg", all=T)[,c(1,3)] 
#[1] 12640     2 (need to be the same legnth/order to KEGG2TAIR)
keggrest <- compareCluster(geneClusters = tair_list,
                           fun = "enricher", 
                           TERM2GENE = KEGG2TAIR,
                           TERM2NAME = KEGG2NAME)
dotplot(keggrest, group=TRUE) ## The same result to "KEGG" 
keggrest2 <- pairwise_termsim(keggrest)
emapplot(keggrest2)

keggout <- as.data.frame(keggrest2@compareClusterResult)
keggout$GeneRatio <- sapply(keggout$GeneRatio, function(x)paste0(" ",x)) #prevent the automatic date-conversion in Excel.
keggout$BgRatio <- sapply(keggout$BgRatio, function(x)paste0(" ",x))
write.csv(keggout,"./3_Marker-genes/Cluster-marekre.KEGG.csv", quote=TRUE, row.names=FALSE)

############## The same result of "KEGG" in-front
################################################


integrated.markers %>% 
  filter(p_val_adj < 0.01) %>% 
  filter(avg_log2FC >= 1) -> sig.markers
table(sig.markers$cluster)
#0  2  3  4  5  6  7  8  9 10 11 12 
#10  1 30  5  5 21  3 10  6  1  2 52 


# Appendix. some genes expression
cls3.top5 <- tail(top5.posi$gene, 5)
#[1] "AT1G67870" "AT5G02600" "AT4G00780" "AT5G42980" "AT1G54410"
FeaturePlot(tbt1, reduction="umap", features=cls3.top5, ncol=3)


###6. DEG by mutation
#(if necessary)
setwd("/Users/junesk9/Library/CloudStorage/Box-Box/シングルセル解析(多部田)/")

#(Globally)
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

require(tidyverse)
top3.neg <- degAll %>%
  filter(avg_log2FC <= 0) %>% 
  group_by(cluster) %>%
  top_n(3, -avg_log2FC)
top3.gene <- unique(top3.neg$gene)
#  [1] "AT5G49360" "AT3G13750" "AT2G05380" "AT1G54410" "AT1G20450" "AT5G05410" "AT2G39800" "AT3G26510" "AT4G35770" "AT5G63160" "AT5G52300"
#[12] "AT1G01320" "AT5G05690" "AT5G19140" "AT2G41430" "AT1G20440" "AT2G29290" "AT3G17800" "AT2G33830" "AT5G23010" "AT5G56870" "AT3G52720"
#[23] "AT5G27690" "AT3G15353" "AT5G64240"
FeaturePlot(tbt1, features=top3.gene, ncol=4)
DotPlot(tbt1, features=top3.gene, cols=c("blue","red"), dot.scale=8, split.by="orig.ident") + RotatedAxis()
VlnPlot(tbt1, features=top3.gene[1:8], split.by="orig.ident", pt.size=0, ncol=2)



#####################################################################
##################### Re-run Mar. 27th 2024
#####################################################################
setwd("/Users/junesk9/Library/CloudStorage/Box-Box/シングルセル解析(多部田)/")
tbt1 <- readRDS("Tabeta-fugu5.seurat.rds")
# An object of class Seurat 
# 160942 features across 8071 samples within 7 assays 
# Active assay: SCT (21377 features, 2000 variable features)
# 3 layers present: counts, data, scale.data
# 6 other assays present: RNA, spliced_RNA, unspliced_RNA, spliced_SCT, unspliced_SCT, ALRA
# 6 dimensional reductions calculated: pca, pca_whole, pca_vst2k, umap, umap_whole_all, umap_whole_pca14
meta <- tbt1@meta.data
# [1] "orig.ident"                       "nCount_RNA"                       "nFeature_RNA"                    
# [4] "nCount_spliced_RNA"               "nFeature_spliced_RNA"             "nCount_unspliced_RNA"            
# [7] "nFeature_unspliced_RNA"           "percent.mt"                       "percent.cp"                      
# [10] "nCount_SCT"                       "nFeature_SCT"                     "pANN_0.25_0.15_59"               
# [13] "DF.classifications_0.25_0.15_59"  "nCount_spliced_SCT"               "nFeature_spliced_SCT"            
# [16] "nCount_unspliced_SCT"             "nFeature_unspliced_SCT"           "SCT_snn_res.0.5"                 
# [19] "seurat_clusters"                  "pANN_0.25_0.15_238"               "DF.classifications_0.25_0.15_238"
# [22] "CytoTRACE"                        "SCT_snn_res.1"                    "SCT_snn_res.1.5"                 
# [25] "SCT_snn_res.2"                   

annt <- read.csv("0_参考になるもの/TAIR10.gene_annt.csv", header=T, row.names=1)
colnames(annt) <- c("logo","desc")

##(imputation with ALRA)
counts <- tbt1@assays$RNA@counts
#[1] 26512  8071 gene x cell

counts <- t(as.matrix(counts)) #transpose to cell x gene
c_norm <- normalize_data(counts) 
k_choice <- choose_k(c_norm) #38
#(details for k)
library(ggplot2)
library(gridExtra)
df <- data.frame(x=1:100,y=k_choice$d)
g1<-ggplot(df,aes(x=x,y=y),) + geom_point(size=1)  + geom_line(size=0.5)+ geom_vline(xintercept=k_choice$k)   + theme( axis.title.x=element_blank() ) + scale_x_continuous(breaks=seq(10,100,10)) + ylab('s_i') + ggtitle('Singular values')
df <- data.frame(x=2:100,y=diff(k_choice$d))[3:99,]
g2<-ggplot(df,aes(x=x,y=y),) + geom_point(size=1)  + geom_line(size=0.5)+ geom_vline(xintercept=k_choice$k+1)   + theme(axis.title.x=element_blank() ) + scale_x_continuous(breaks=seq(10,100,10)) + ylab('s_{i} - s_{i-1}') + ggtitle('Singular value spacings')
grid.arrange(g1,g2,nrow=1)

c_alra <- alra(c_norm, k=k_choice$k)[[3]]
# Scaling all except for 2956 columns
# 0.00% of the values became negative in the scaling process and were set to zero
# The matrix went from 2.89% nonzero to 26.95% nonzero
c_alra.t <- Matrix(t(c_alra), sparse=TRUE)
c_alra.t2 <- as(t(c_alra), "dgCMatrix")
all.equal(c_alra.t, c_alra.t2) #[1] TRUE # alternative commands, both are identical
colnames(c_alra.t) <- rownames(counts) #the col-names lost during process
tbt1[["ALRA"]] <- CreateAssayObject(data=c_alra.t)

saveRDS(tbt1, "Tabeta-fugu5.seurat.rds")
DefaultAssay(tbt1) #SCT
DefaultAssay(tbt1) <- "ALRA" ##change the expression DB


############# Tissue-origin estimation
##(tissue estimation by given genes; KIM J et al. 2021)
VlnPlot(tbt1, features=c("AT3G48740","AT5G23660")) 
#clus9/10 as Phloem parenchyma
VlnPlot(tbt1, features=c("AT1G22710","AT5G06850","AT1G79430","AT5G57350"), ncol=3)
#clus3 as companion cell
VlnPlot(tbt1, features=c("AT5G41920","AT1G77990")) 
#clus12/3 as bundle sheet
VlnPlot(tbt1, features=c("AT3G51480","AT5G19530"))
#clus10?? as xylem
VlnPlot(tbt1, features=c("AT1G29910","AT2G05100","AT3G01500"), ncol=2)
#clus12 as the central mesophyll
VlnPlot(tbt1, features=c("AT2G45190"))
#clus11, 6 to spongy mesophyll
VlnPlot(tbt1, features=c("AT1G19850","AT4G32880","AT5G61480","AT1G46480"), ncol=2)
#clus9/10 as procambium
VlnPlot(tbt1, features=c("AT1G09310","AT4G21750","AT1G68530"), ncol=2)
#clus0/6 as the epidermis
VlnPlot(tbt1, features=c("AT3G24140","AT2G46070","AT3G18040","AT1G22690"), ncol=2)
#unknown to guard cell
VlnPlot(tbt1, features=c("AT1G28230","AT3G54420")) 
#unknown to hydathode

##(Correlation-based; J Kim et al. 2021); 
ref2 <- ReadMtx(mtx="0_参考になるもの/GSE161332/GSE161332_matrix.mtx.gz",
                cells = "0_参考になるもの/GSE161332/GSE161332_barcodes.tsv.gz",
                features="0_参考になるもの/GSE161332/GSE161332_features.tsv.gz")
ref2 <- CreateSeuratObject(ref2)
#no tissue annotation available; STOPPED


##(Correlation-based; Berrio et al. 2022)
ref <- readRDS("0_参考になるもの/GSE213622_coaker_atpst_singlecell_sobj.rds")
# An object of class Seurat 
# 145021 features across 11895 samples within 6 assays 
# Active assay: integrated (14941 features, 14941 variable features)
# 2 layers present: data, scale.data
# 5 other assays present: spliced, unspliced, ambiguous, RNA, SCT
# 3 dimensional reductions calculated: pca, umap, umap3d
rmeta <- ref@meta.data
tissue1 <- unique(rmeta$predicted.id)
#[1] "Mesophyll"       "Companion Cells" "Bundle Sheath"   "Guard Cell"      "Parenchyma"      "Procambium"      "Hydathode"       "Epidermis"
tissue2 <- unique(rmeta$Phase)
#[1] "G1"  "G2M" "S"  

rcnt <- ref@assays$SCT@data #[1] 19460 11895
tcnt <- tbt1@assays$SCT@data #[1] 21377  8071
#(select diverse most-1000 diverse genes)
rdf <- as.data.frame(rcnt) #[1] 11895 18201
rdf$rowVar = apply(rdf, 1, var)
best.g <- rownames(rdf[order(-rdf$rowVar),])[c(1:1000)]
rdf <- NULL #memory plunge

gcmn <- intersect(best.g, rownames(tcnt)) #[1] 993
rcnt <- rcnt[gcmn,] #[1]  993 11895
tcnt <- tcnt[gcmn, ] #[1] 993 8071

#(rcnt subset by reduce the "mesophyll" cells 10466 -> 500)
rmeta.sub <- subset(rmeta, rmeta$predicted.id %in% "Mesophyll")
rmeta.sub <- rmeta.sub[order(-rmeta.sub$prediction.score.Mesophyll), ]
rmeso500 <- rownames(rmeta.sub)[1:500]
rother <- rownames(subset(rmeta, rmeta$predicted.id != "Mesophyll"))
rmeta.s <- rmeta[rownames(rmeta) %in% c(rother, rmeso500), ] #[1] 1929   38
#Bundle Sheath Companion Cells       Epidermis      Guard Cell       Hydathode       Mesophyll      Parenchyma      Procambium 
#564             145              21             167             163             500             187             182 
rcnt.s <- rcnt[, rownames(rmeta.s)] #[1]  993 1929

pcc.df <- data.frame(row.names=colnames(tcnt))
for(t in c(tissue1, tissue2)){pcc.df[,t] <- 0}
pcc.df$Tissue <- NA
pcc.df$Phase <- NA
#pcc.df <- read.csv("4_tissue.annotation/pcc.df.csv", row.names=1, check.names=FALSE)
#pcc.df[is.null(pcc.df)] <- 0
rcnt <- rcnt.s
rmeta <- rmeta.s
for (tcol in colnames(tcnt)){
    if (NA %in% pcc.df[tcol, c("Phase","Tissue")]){
        tvec <- as.vector(tcnt[,tcol])
        for (rcol in colnames(rcnt)){
            rvec <- as.vector(rcnt[,rcol])
            r.tissue <- rmeta[rcol, "predicted.id"]
            r.tissue2 <- rmeta[rcol, "Phase"]
            pcc <- cor.test(tvec, rvec, method="pearson")[c(3,4)]
            pcc <- pcc$estimate
    
            if (pcc.df[tcol, r.tissue] < pcc){pcc.df[tcol,r.tissue] <- pcc}
            if (pcc.df[tcol, r.tissue2] < pcc){pcc.df[tcol,r.tissue2] <- pcc}
            }
  tissue.pcc <- pcc.df[tcol, tissue1]
  tissue.idx <- which(tissue.pcc == max(tissue.pcc))
  pcc.df[tcol, "Tissue"] <- tissue1[tissue.idx]
  phase.pcc <- pcc.df[tcol, tissue2]
  phase.idx <- which(phase.pcc == max(phase.pcc))
  pcc.df[tcol, "Phase"] <- tissue2[phase.idx]
  }
  
  col.idx <- which(colnames(tcnt) == tcol)
  if (col.idx %% 100 == 0){print(col.idx)} 
}
write.csv(pcc.df, "4_tissue.annotation/pcc.df.csv")
pcc.names <- colnames(pcc.df)[1:(length(pcc.df)-2)]
pcc.names <- as.vector(sapply(pcc.names, function(x){paste0(x, "_PCC")}))
colnames(pcc.df)[1:(length(pcc.df)-2)] <- pcc.names
pcc.df <- pcc.df[rownames(meta), ] #ensuring
for(cn in colnames(pcc.df)){tbt1@meta.data[, cn] <- pcc.df[, cn]}


DimPlot(tbt1, group.by=c("Tissue","est.tissue_by.markers"))


########################################################
################ scVelo
######################################################
######### Prepare scVelo input
sr <- tbt1@assays$spliced_RNA@counts
ur <- tbt1@assays$unspliced_RNA@counts
ar <- tbt1@assays$RNA@counts

sr <- as.matrix(sr/ar)
ur <- as.matrix(ur/ar)
sr[is.nan(sr)] = 0
ur[is.nan(ur)] = 0

sg <- intersect(rownames(sr), rownames(ur))
#[1] 26512
spliced <- sr[match(sg, rownames(sr)), ]
unspliced <- ur[match(sg, rownames(ur)), ]

pca <- tbt1@reductions$pca@cell.embeddings
umap <- tbt1@reductions$umap@cell.embeddings

save(spliced, unspliced, meta, sg, pca, umap, file="tbt1_scVelo_input.RData")
#find "mono-all.scVelo.ipynb" for the following python procedure.

scv.out <- read.csv("5_scVelo/tbt1.scVelo-obs.csv", header=T, row.names=1)
#[1] 8071   31
scv.out <- scv.out[rownames(meta), ] #for ensuring
tbt1$velocity_self_transition <- scv.out$velocity_self_transition
tbt1$velocity_pseudotime <- scv.out$velocity_pseudotime
tbt1$latent_time <- scv.out$latent_time
# Function for normalizing values to range 0 to 1
range01 <- function(x) {(x - min(x))/(max(x) - min(x))}
tbt1$consensus.time <- range01((tbt1$CytoTRACE+tbt1$latent_time)/2)

saveRDS(tbt1, "Tabeta-fugu5.seurat.rds")

##############################################################
####### GOI visualization
###########################################################
goi <- c("AT5G57220","AT3G09260","AT3G16400")
#before imputation
DefaultAssay(tbt1) <- "SCT"
FeaturePlot(tbt1, features=goi)
#After imputation
DefaultAssay(tbt1) <-"ALRA"
FeaturePlot(tbt1, features=goi, cols=hcl.colors(100, "Oslo")[10:100])






########################################################
################ Appendix
######################################################
VlnPlot(tbt1, features=c("AT4G04840","AT5G62210","AT5G13930","AT1G30530","AT1G18810","AT3G51240"), ncol=3) 
#clus0 as the epidermis:adaxial pavement
VlnPlot(tbt1, features=c("AT2G05520","AT1G06360","AT1G09310","AT4G38770","AT3G51600","AT3G16370"), ncol=3) 
#clus6,0 as the epidermis:abaxial pavement
VlnPlot(tbt1, features=c("AT5G59870","AT5G10400","AT5G22880","AT3G46320","AT4G27230","AT3G53650"), ncol=3) 
#unknown for Late-S phase

meta$seurat_clusters <- Idents(tbt1)
meta$est.tissue <- "unknown"
meta[meta$seurat_clusters %in% "9",]$est.tissue <- "Phloem parenchyma"
meta[meta$seurat_clusters %in% "3",]$est.tissue <- "companion cell"
meta[meta$seurat_clusters %in% c("12","2","8"),]$est.tissue <- "mesophyll"
meta[meta$seurat_clusters %in% "10",]$est.tissue <- "procambium"
meta[meta$seurat_clusters %in% c("6","0"),]$est.tissue <- "epidermis"
tbt1@meta.data$est.tissue <- meta$est.tissue
tbt1@meta.data$seurat_clusters <- Idents(tbt1)

DimPlot(tbt1, label=T)
DimPlot(tbt1, group.by="est.tissue", cols=c(ggcolor(5),"grey")) #check the estimation


##(check GOI expression pattern)
VlnPlot(tbt1, features=c("AT5G57220","AT3G09260","AT3G16400"), ncol=3)
#All three at cluster 6 & 0; meaning epidermis

##(difference between cluster 0 vs. 6)
deg6.0 <- FindMarkers(tbt1, group.by="seurat_clusters", ident.1="6", ident.2="0", test.use="MAST")
deg <- subset(deg6.0, deg6.0$p_val_adj < 0.05 & abs(deg6.0$avg_log2FC) > 1)
#[1] 16  5
deg$annt <- annt[rownames(deg), ]
write.csv(deg, "Tabeta-fugu5.cluster6vs0.DEG.csv")

FeaturePlot(tbt1, features=rownames(deg), ncol=4)


##########################################################################
############################## subsetting by request
###########################################################################
tbt1.epi <- subset(tbt1, subset = est.tissue_by.markers == "epidermis")
# An object of class Seurat 
# 160942 features across 1904 samples within 7 assays 
# Active assay: SCT (21377 features, 2000 variable features)
# 6 other assays present: RNA, spliced_RNA, unspliced_RNA, spliced_SCT, unspliced_SCT, ALRA
# 6 dimensional reductions calculated: pca, pca_whole, pca_vst2k, umap, umap_whole_all, umap_whole_pca14
saveRDS(tbt1.epi, "../Tabeta-fugu5.subset-epidermis.seurat.rds")
DimPlot(tbt1.epi, group.by="orig.ident")





sessionInfo()
#R version 4.3.1 (2023-06-16)
#Platform: x86_64-apple-darwin20 (64-bit)
#Running under: macOS Monterey 12.7.1

#Matrix products: default
#BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
#LAPACK: /Library/Frameworks/R.framework/Versions/4.3-x86_64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0

#locale:
#  [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

#time zone: Asia/Tokyo
#tzcode source: internal

#attached base packages:
#  [1] stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#[1] COPILOT_0.1.0               sf_1.0-14                   pheatmap_1.0.12             tradeSeq_1.16.0             slingshot_2.10.0           
#[6] TrajectoryUtils_1.10.0      princurve_2.1.6             glmnet_4.1-8                Matrix_1.6-4                CytoTRACE_0.3.3            
#[11] MAST_1.28.0                 SingleCellExperiment_1.24.0 SummarizedExperiment_1.32.0 Biobase_2.62.0              GenomicRanges_1.54.1       
#[16] GenomeInfoDb_1.38.1         IRanges_2.36.0              S4Vectors_0.40.2            BiocGenerics_0.48.1         MatrixGenerics_1.14.0      
#[21] matrixStats_1.2.0           lubridate_1.9.3             forcats_1.0.0               stringr_1.5.1               dplyr_1.1.4                
#[26] purrr_1.0.2                 readr_2.1.4                 tidyr_1.3.0                 tibble_3.2.1                ggplot2_3.4.4              
#[31] tidyverse_2.0.0             SeuratObject_5.0.1          Seurat_4.4.0  
