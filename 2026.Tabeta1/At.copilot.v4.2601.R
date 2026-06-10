############################################################
###### COPILOT scRNA-seq data analysis
###### 2026.01.14 JS KIM
##### a bit more for GO/KEGG
############################################################

library(Seurat)
library(pheatmap)
library(clusterProfiler)
library(enrichplot)
library(MAST)
library(tidyverse)
library(org.At.tair.db)

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
tbt1 <- readRDS("Tabeta-fugu5.seurat.v4.rds")
# An object of class Seurat 
# 160907 features across 7960 samples within 7 assays 
# Active assay: SCT (21342 features, 2000 variable features)
# 3 layers present: counts, data, scale.data
# 6 other assays present: RNA, spliced_RNA, unspliced_RNA, spliced_SCT, unspliced_SCT, ALRA
# 9 dimensional reductions calculated: pca, pca_all, pca_vst2k, umap, umap_vst2k, umap_3d, tsne, tsne_vst2k, tsne_3d

############################################################
####### Re-Run ClusterProfiler for tissues
############################################################
tbt1$gt.tissue <- paste(tbt1$orig.ident, tbt1$cons.Tissue, sep="_")
degAll <- list()
for (tis in names(table(tbt1$cons.Tissue))){
  ident1 <- paste0("fugu5_8DAS_",tis)
  ident2 <- paste0("WT_8DAS_",tis)
  condition.diffgenes <- FindMarkers(tbt1, group.by="gt.tissue", ident.1=ident1, ident.2=ident2, min.pct=0.25, logfc.threshold=0.25,test.use="MAST")
  condition.diffgenes$tissue <- tis
  condition.diffgenes$gene <- row.names(condition.diffgenes)
  rownames(condition.diffgenes) <- NULL
  condition.diffgenes <- condition.diffgenes[,c(7,6,2,1,5,3,4)] #re-order the column
  write.csv(condition.diffgenes, file=paste0("7_DEGs/tissue-DEGs.",tis,".csv"), row.names=F, quote=F)
  degAll <- c(degAll, list(condition.diffgenes))
}

degAll <- bind_rows(degAll) #[1] 1289    7
write.csv(degAll, "7_DEGs/tissue-DEGs.merged.csv", row.names=F, quote=F)
degSig <- degAll[degAll$p_val_adj < 0.05 & degAll$avg_log2FC >= 0.25, ] #fugu5-up DEGs #[1] 129   7
table(degSig$tissue)
# Bundle.Sheath Companion.Cells       Epidermis       Mesophyll      Procambium 
#         1              21              53              51               3 

#(ClusterCompare: KEGG, GO, and so on)
#collect DEG gids by cluster
deg_list <- split(degSig$gene, degSig$tissue)
names(deg_list) <- names(table(degSig$tissue)) 
tair_list <- deg_list #need to keep for KEGG analysis
for(i in names(deg_list)){
  deg_list[[i]] <- bitr(deg_list[[i]],
                        fromType="TAIR",
                        toType="ENTREZID",
                        OrgDb="org.At.tair.db",
                        drop=TRUE) #drop for remove the unconverted gid, keep them with "NA" by "FALSE"
  deg_list[[i]] <- deg_list[[i]]$ENTREZID 
} #deg list including ENTREZID

GOBP <- compareCluster(geneClusters=deg_list, fun="enrichGO", ont="BP", OrgDb = "org.At.tair.db", pvalueCutoff=0.9)
GOMF <- compareCluster(geneClusters=deg_list, fun="enrichGO", ont="MF", OrgDb = "org.At.tair.db", pvalueCutoff=0.9)
Ent2Tair <- function(vec){
  ent = strsplit(vec, "/")[[1]]
  tair = suppressMessages(bitr(ent, fromType="ENTREZID", toType="TAIR", OrgDb="org.At.tair.db"))[,2]
  tair = paste(tair, collapse="/")
  return(tair)}

goout <- as.data.frame(GOBP@compareClusterResult)
goout$geneID <- sapply(goout$geneID, Ent2Tair)
goout$GeneRatio <- sapply(goout$GeneRatio, function(x)paste0(" ",x)) #prevent the automatic date-conversion in Excel.
goout$BgRatio <- sapply(goout$BgRatio, function(x)paste0(" ",x))
write.csv(goout,"./7_DEGs/tissue-DEGs.merged.GOBP.csv", quote=TRUE, row.names=FALSE)

goout <- as.data.frame(GOMF@compareClusterResult)
goout$geneID <- sapply(goout$geneID, Ent2Tair)
goout$GeneRatio <- sapply(goout$GeneRatio, function(x)paste0(" ",x)) #prevent the automatic date-conversion in Excel.
goout$BgRatio <- sapply(goout$BgRatio, function(x)paste0(" ",x))
write.csv(goout,"./7_DEGs/tissue-DEGs.merged.GOMF.csv", quote=TRUE, row.names=FALSE)


###################################
######### DEG again to confirm
library(MAST)

DefaultAssay(tbt1) <- "SCT"
Idents(tbt1) <- "cons.Tissue" #34 in total; 6 in epidermis
tbt1@meta.data$orig.ident <- factor(tbt1@meta.data$orig.ident, levels=unique(tbt1@meta.data$orig.ident)) #set order WT -> fugu5

deg.g <- FindMarkers(tbt1, group.by="orig.ident", ident.1="fugu5_8DAS", test.use="MAST", only.pos=FALSE, logfc.threshold = 0.1) #[1] 274   5
deg.epi <- FindMarkers(tbt1, group.by="orig.ident", ident.1="fugu5_8DAS", test.use="MAST", 
                       only.pos=FALSE, subset.ident ="Epidermis", min.pct=0.1, logfc.threshold = 0.1) #[1] 413   5
deg.bs <- FindMarkers(tbt1, group.by="orig.ident", ident.1="fugu5_8DAS", test.use="MAST", 
                      only.pos=FALSE, subset.ident ="Bundle.Sheath", min.pct=0.1, logfc.threshold = 0.1) #[1] 735   5
deg.cc <- FindMarkers(tbt1, group.by="orig.ident", ident.1="fugu5_8DAS", test.use="MAST", 
                      only.pos=FALSE, subset.ident ="Companion.Cells", min.pct=0.1, logfc.threshold = 0.1) #[1] 399   5
deg.gc <- FindMarkers(tbt1, group.by="orig.ident", ident.1="fugu5_8DAS", test.use="MAST", 
                      only.pos=FALSE, subset.ident ="Guard.Cell", min.pct=0.1, logfc.threshold = 0.1) #[1] 2567    5
deg.mes <- FindMarkers(tbt1, group.by="orig.ident", ident.1="fugu5_8DAS", test.use="MAST", 
                       only.pos=FALSE, subset.ident ="Mesophyll", min.pct=0.1, logfc.threshold = 0.1) #[1] 303   5
deg.pc <- FindMarkers(tbt1, group.by="orig.ident", ident.1="fugu5_8DAS", test.use="MAST", 
                      only.pos=FALSE, subset.ident ="Procambium", min.pct=0.1, logfc.threshold = 0.1) #[1] 527   5

deg.epi$gid <- rownames(deg.epi)
deg.bs$gid <- rownames(deg.bs)
deg.cc$gid <- rownames(deg.cc)
deg.gc$gid <- rownames(deg.gc)
deg.mes$gid <- rownames(deg.mes)
deg.pc$gid <- rownames(deg.pc)

deg.epi$Tissue <- "Epidermis"
deg.bs$Tissue <- "Bundle.Sheath"
deg.cc$Tissue <- "Companion.Cells"
deg.gc$Tissue <- "Guard.Cell"
deg.mes$Tissue <- "Mesophyll"
deg.pc$Tissue <- "Procambium"

deg.all <- rbind(deg.epi, deg.bs, deg.cc, deg.gc, deg.mes, deg.pc, make.row.names = FALSE) #[1] 4944    6
write.csv(deg.all, "7_DEGs/DEG-FvsW.byTissue.csv")
write.csv(deg.g, "7_DEGs/DEG-FvsW.global.csv")

Idents(tbt1) <- "SCT_snn_res.3_subE"

deg.e1 <- FindMarkers(tbt1, group.by="orig.ident", ident.1="fugu5_8DAS", test.use="MAST", 
                      only.pos=TRUE, subset.ident ="e1", min.pct=0.1, logfc.threshold = 0.1) #[1] 352   5
deg.e2 <- FindMarkers(tbt1, group.by="orig.ident", ident.1="fugu5_8DAS", test.use="MAST", 
                      only.pos=TRUE, subset.ident ="e2", min.pct=0.1, logfc.threshold = 0.1) #[1] 328   5
deg.e3 <- FindMarkers(tbt1, group.by="orig.ident", ident.1="fugu5_8DAS", test.use="MAST", 
                      only.pos=TRUE, subset.ident ="e3", min.pct=0.1, logfc.threshold = 0.1) #[1] 328   5
deg.e4 <- FindMarkers(tbt1, group.by="orig.ident", ident.1="fugu5_8DAS", test.use="MAST", 
                      only.pos=TRUE, subset.ident ="e4", min.pct=0.1, logfc.threshold = 0.1) #[1] 685   5
deg.e5 <- FindMarkers(tbt1, group.by="orig.ident", ident.1="fugu5_8DAS", test.use="MAST", 
                      only.pos=TRUE, subset.ident ="e5", min.pct=0.1, logfc.threshold = 0.1) #[1] 326   5
deg.e6 <- FindMarkers(tbt1, group.by="orig.ident", ident.1="fugu5_8DAS", test.use="MAST", 
                      only.pos=TRUE, subset.ident ="e6", min.pct=0.1, logfc.threshold = 0.1) #[1] 326   5
deg.e7 <- FindMarkers(tbt1, group.by="orig.ident", ident.1="fugu5_8DAS", test.use="MAST", 
                      only.pos=TRUE, subset.ident ="e7", min.pct=0.1, logfc.threshold = 0.1) #[1] 523   5
deg.e8 <- FindMarkers(tbt1, group.by="orig.ident", ident.1="fugu5_8DAS", test.use="MAST", 
                      only.pos=TRUE, subset.ident ="e8", min.pct=0.1, logfc.threshold = 0.1) #[1] 168   5
subE.markers <- FindAllMarkers(object=tbt1, only.pos=TRUE, min.pct=0.1, logfc.threshold = 0.1, test.use="MAST")

deg.e1$gid <- rownames(deg.e1)
deg.e2$gid <- rownames(deg.e2)
deg.e3$gid <- rownames(deg.e3)
deg.e4$gid <- rownames(deg.e4)
deg.e5$gid <- rownames(deg.e5)
deg.e6$gid <- rownames(deg.e6)
deg.e7$gid <- rownames(deg.e7)
deg.e8$gid <- rownames(deg.e8)

deg.e1$Tissue <- "e1"
deg.e2$Tissue <- "e2"
deg.e3$Tissue <- "e3"
deg.e4$Tissue <- "e4"
deg.e5$Tissue <- "e5"
deg.e6$Tissue <- "e6"
deg.e7$Tissue <- "e7"
deg.e8$Tissue <- "e8"

deg.subE <- rbind(deg.e1, deg.e2, deg.e3, deg.e4, deg.e5, deg.e6, deg.e7,deg.e8, make.row.names = FALSE) #[1] 2382    7
write.csv(deg.subE, "7_DEGs/DEG-FvsW.subEpi_clusters.MAST.csv", quote=F, row.names=F)
write.csv(subE.markers, "7_DEGs/subEpi_clusters.markers.MAST.csv")



sessionInfo()
# R version 4.5.2 (2025-10-31)
# Platform: aarch64-apple-darwin20
# Running under: macOS Sequoia 15.7.3
# 
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
# LAPACK: /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# time zone: Asia/Tokyo
# tzcode source: internal
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] scales_1.4.0    lubridate_1.9.4 forcats_1.0.1   stringr_1.6.0   dplyr_1.1.4     purrr_1.2.0     readr_2.1.6     tidyr_1.3.2    
# [9] tibble_3.3.0    tidyverse_2.0.0 cowplot_1.2.0   ggplot2_4.0.1  
# 
# loaded via a namespace (and not attached):
#   [1] deldir_2.0-4           pbapply_1.7-4          gridExtra_2.3          rlang_1.1.6            magrittr_2.0.4        
# [6] RcppAnnoy_0.0.23       otel_0.2.0             spatstat.geom_3.6-1    matrixStats_1.5.0      ggridges_0.5.7        
# [11] compiler_4.5.2         png_0.1-8              vctrs_0.6.5            reshape2_1.4.5         pkgconfig_2.0.3       
# [16] fastmap_1.2.0          labeling_0.4.3         promises_1.5.0         tzdb_0.5.0             jsonlite_2.0.0        
# [21] goftest_1.2-3          later_1.4.4            spatstat.utils_3.2-1   irlba_2.3.5.1          parallel_4.5.2        
# [26] cluster_2.1.8.1        R6_2.6.1               ica_1.0-3              stringi_1.8.7          RColorBrewer_1.1-3    
# [31] spatstat.data_3.1-9    reticulate_1.44.1      parallelly_1.46.0      spatstat.univar_3.1-5  lmtest_0.9-40         
# [36] scattermore_1.2        Rcpp_1.1.0             tensor_1.5.1           future.apply_1.20.1    zoo_1.8-15            
# [41] sctransform_0.4.3      httpuv_1.6.16          Matrix_1.7-4           splines_4.5.2          igraph_2.2.1          
# [46] timechange_0.3.0       tidyselect_1.2.1       rstudioapi_0.17.1      abind_1.4-8            spatstat.random_3.4-3 
# [51] codetools_0.2-20       miniUI_0.1.2           spatstat.explore_3.6-0 listenv_0.10.0         lattice_0.22-7        
# [56] plyr_1.8.9             shiny_1.12.1           withr_3.0.2            S7_0.2.1               ROCR_1.0-11           
# [61] Rtsne_0.17             future_1.68.0          fastDummies_1.7.5      survival_3.8-3         polyclip_1.10-7       
# [66] fitdistrplus_1.2-4     pillar_1.11.1          Seurat_5.4.0           KernSmooth_2.23-26     plotly_4.11.0         
# [71] generics_0.1.4         RcppHNSW_0.6.0         sp_2.2-0               hms_1.1.4              globals_0.18.0        
# [76] xtable_1.8-4           glue_1.8.0             lazyeval_0.2.2         tools_4.5.2            data.table_1.18.0     
# [81] RSpectra_0.16-2        RANN_2.6.2             dotCall64_1.2          grid_4.5.2             nlme_3.1-168          
# [86] patchwork_1.3.2        cli_3.6.5              spatstat.sparse_3.1-0  spam_2.11-3            viridisLite_0.4.2     
# [91] uwot_0.2.4             gtable_0.3.6           digest_0.6.39          progressr_0.18.0       ggrepel_0.9.6         
# [96] htmlwidgets_1.6.4      SeuratObject_5.3.0     farver_2.1.2           htmltools_0.5.9        lifecycle_1.0.4       
# [101] httr_1.4.7             mime_0.13              MASS_7.3-65           