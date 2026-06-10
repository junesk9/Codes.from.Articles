############################################################
###### COPILOT scRNA-seq data analysis
###### 2025.11.20 JS KIM
##### the next steps
############################################################

library(Seurat)
library(pheatmap)
library(clusterProfiler)
library(enrichplot)

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
tbt1 <- readRDS("Tabeta-fugu5.seurat.v3a.rds")
# An object of class Seurat 
# 160907 features across 7960 samples within 7 assays 
# Active assay: SCT (21342 features, 2000 variable features)
# 3 layers present: counts, data, scale.data
# 6 other assays present: RNA, spliced_RNA, unspliced_RNA, spliced_SCT, unspliced_SCT, ALRA
# 9 dimensional reductions calculated: pca, pca_all, pca_vst2k, umap, umap_vst2k, umap_3d, tsne, tsne_vst2k, tsne_3d


############################################################
####### Re-Run ClusterProfiler and refine the visualization
############################################################
tbt1$gt.clus <- paste(tbt1$orig.ident, tbt1$seurat_clusters, sep="_")
for (cl in names(table(tbt1$seurat_clusters))){
  ident1 <- paste0("fugu5_8DAS_",cl)
  ident2 <- paste0("WT_8DAS_",cl)
  condition.diffgenes <- FindMarkers(tbt1, group.by="gt.clus", ident.1=ident1, ident.2=ident2, min.pct=0.25, logfc.threshold=0.25,test.use="MAST")
  write.csv(condition.diffgenes, file=paste0("clus-DEGs.",cl,".csv"))
}

degAll <- read.csv("7_DEGs/AllMarkers.subEpi_cluster.csv", header=T)
degsig <- degAll[degAll$p_val_adj < 0.05 & degAll$avg_log2FC >= 0.25, ] #fugu5-up DEGs
table(degsig$cluster)
# e1  e2  e3  e4  e5  e6  e7  e8 
# 172 302 211 296 173 167 158  30 



#(ClusterCompare: KEGG, GO, and so on)
#collect DEG gids by cluster
deg_list <- split(degsig$gene, degsig$cluster)
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
GOBP <- compareCluster(geneClusters=deg_list, fun="enrichGO", ont="BP", OrgDb = "org.At.tair.db", pvalueCutoff=0.9)
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
write.csv(goout,"./7_DEGs/AllMarkers.subEpi_cluster.GOBP.csv", quote=TRUE, row.names=FALSE)


KEGG <- compareCluster(geneClusters=tair_list, fun = "enrichKEGG", organism = "ath", pvalueCutoff=0.9)
KEGG2 <- pairwise_termsim(KEGG)
dotplot(KEGG, group=TRUE)
emapplot(KEGG2)

kgout <- as.data.frame(KEGG@compareClusterResult)
kgout$geneID <- kgout$geneID
kgout$GeneRatio <- sapply(kgout$GeneRatio, function(x)paste0(" ",x))
kgout$BgRatio <- sapply(kgout$BgRatio, function(x)paste0(" ",x))
write.csv(kgout,"./7_DEGs/AllMarkers.subEpi_cluster.KEGG.csv", quote=TRUE, row.names=FALSE)

############################################################
#### Improved Visualization 
library(ggplot2)
library(cowplot)
library(tidyverse)
library(scales)

setwd("/Users/junesk9/Library/CloudStorage/Box-Box/シングルセル解析(多部田)")
set.seed(101) #prevent the random output
goout <- read.csv("7_DEGs/AllMarkers.subEpi_cluster.GOBP.csv", header=T)
kgout <- read.csv("7_DEGs/AllMarkers.subEpi_cluster.KEGG.csv", header=T)


top_terms <- goout %>%
  group_by(Cluster) %>%
  slice_max(FoldEnrichment, n = 2) %>%
  ungroup() %>%
  pull(ID) %>%
  unique()

plot_data <- goout %>%
  mutate(
    log2_FE = log2(FoldEnrichment),
    neg_log_p = -log10(pvalue),
    neg_log_p_scaled = scale(neg_log_p)[,1],
    neg_log_p_capped = pmin(neg_log_p, 20),
    log2_FE_capped = pmax(log2_FE, -3),
    log2_FE_capped = pmin(log2_FE_capped, 3),
    neg_log_p_scaled_capped = pmin(neg_log_p_scaled, 4)
  ) %>%
  filter(ID %in% top_terms) %>%
  group_by(Description) %>%
  mutate(sum_log_p = sum(abs(neg_log_p))) %>%
  mutate(max_log_p = max(neg_log_p)) %>%
  ungroup() %>%
  arrange(desc(max_log_p), desc(sum_log_p)) %>%
  mutate(Description = factor(Description, levels = rev(unique(Description))))


# Create the dotplot
DOTPLOT2 <- function(plot_data, p_max=10, fc_max=6){
  p <- ggplot(plot_data, aes(x = factor(Cluster), y = Description)) +
    geom_point(aes(size =  neg_log_p_capped, 
                   color = log2_FE)) +
    
    # Color scale: blue (depleted) to red (enriched)
    scale_color_gradientn(
      colors = c("#4169E1", "#F9E79F", "#FF1100"),
      values = rescale(c(0, fc_max/2, fc_max)),  # Beige only from -0.4 to +0.4, then transitions
      limits = c(0, fc_max),
      name = "log2(FE)"
    ) +
    
    # Size scale for significance
    scale_size_continuous(
      name = "signif.",
      range = c(1, p_max),
      breaks = c(0, round(p_max/5), round(p_max/2), p_max),
      #labels = c("1", "0.001", "1e-15", "<1e-30"),
      guide = guide_legend(override.aes = list(color = "black"))
    ) +
    
    # Themes and labels
    labs(
      title = "Gene Set Enrichment Analysis",
      subtitle = "fugu-upDEGs",
      x = "Gene Clusters",
      y = "Functional terms",
      caption = "Dot size = significance level, Color = fold enrichment"
    ) +
    theme_minimal() 
  
  # Display the plot
  print(p)
}
DOTPLOT2(plot_data, fc_max=6)




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
