####Sleuth analysis of wheat RNA-seq timecourse (Mega)
####2022.01.26
####KIM JS


library(sleuth)
s2c <- read.table("wheat.mega1.txt", header=TRUE, stringsAsFactors=FALSE)
head(s2c)
#sample  osmo geno day              path    run
#1   L8B-1  None   L8   0  3_kallisto/L8B-1  L8B-1
#2   L8B-2  None   L8   0  3_kallisto/L8B-2  L8B-2
#3   L8B-3  None   L8   0  3_kallisto/L8B-3  L8B-3
#4   L8B-4  None   L8   0  3_kallisto/L8B-4  L8B-4
#5   L8B-5  None   L8   0  3_kallisto/L8B-5  L8B-5

########################################################
############## pairwise DEG
### subset inputs
so <- s2c[c(36:41,78:83),] #day 4(Water)
so <- s2c[c(30:35,72:77),] #day 3(Water)
so <- s2c[c(24:29,66:71),] #day 1(Water)
so <- s2c[c(19:23,60:65),] #day 4(DR)
so <- s2c[c(13:18,54:59),] #day 3(DR)
so <- s2c[c(7:12,48:52),] #day 1(DR)
so <- s2c[c(1:6,42:47),] #day 0

### normalize & modeling
so <- sleuth_prep(so, ~ geno)
so <- sleuth_fit(so)
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')
results_table_lrt <- sleuth_results(so, 'reduced:full', test_type = 'lrt')
table(results_table_lrt[,"qval"] < 0.05)
# [Optional] QQ-plot
#plot_qq(so, test = 'reduced:full', test_type = 'lrt', sig_level = 0.05)

models(so)
#[  full  ]
#formula:  ~geno 
#data modeled:  obs_counts 
#transform sync'ed:  TRUE 
#coefficients:
#	(Intercept)
# 	genoN
#[  reduced  ]
#formula:  ~1 
#data modeled:  obs_counts 
#transform sync'ed:  TRUE 
#coefficients:
#  (Intercept)
so <- sleuth_wt(so, "genoN")
results_table_wt <- sleuth_results(so, 'genoN')

################ Prepare output data
#Extract common DEG from two LRT & WT analyses
d5.lrt.sig_ids <- results_table_lrt$target_id[which(results_table_lrt$qval < 0.05)]
d5.wt.sig_ids <- results_table_wt$target_id[which(results_table_wt$qval < 0.05)]
shared_ids <- d5.wt.sig_ids[d5.wt.sig_ids %in% d5.lrt.sig_ids]
shared_results <- results_table_wt[results_table_wt$target_id %in% shared_ids,]
out <- shared_results[,c(1:4)]
out$log2FC <- -log(exp(out$b),2)   #<<- Add Log2FC column
out <- subset(out, select=-c(b))  #<<- Remove b column
out <- subset(out, out$log2FC >= 1 | out$log2FC <= -1) #<<- select by FC
head(out)
#target_id         pval         qval    log2FC
#1     synTaPYL1a_GFP.1 3.139324e-91 1.965908e-86 12.283900
#2 TRAESCS7D02G549900.1 1.769738e-45 5.541226e-41  3.995749
#3 TRAESCS3B02G356500.1 1.479448e-28 2.316150e-24  2.045975
#4 TRAESCS1D02G237100.1 1.301953e-28 2.316150e-24  1.815227
#5 TRAESCS3B02G277900.1 1.887018e-27 2.363377e-23  1.577205
shared_ids <- as.vector(out$target_id)  #<<-prepare new id list

#extract raw TPM data
tpm_matrix <- sleuth_to_matrix(so, 'obs_raw', 'tpm')
d <- tpm_matrix #for sure
d <- cbind(rownames(d), data.frame(d, row.names=NULL))
colnames(d)[1] <- "target_id"
tpm <- d[d$target_id %in% shared_ids,]

#merge two tables 
tmp <- merge(tpm, out, by="target_id", sort=T)
head(tmp)
#target_id     L8B.1     L8B.2     L8B.3       L8B.4     L8B.5       L8B.6        NB.1      NB.2       NB.3        NB.4       NB.5       NB.6
#1     synTaPYL1a_GFP.1 416.12942 313.22675 403.02762 277.3162398 395.83562 345.7930286  0.04818385  0.000000   0.000000  0.04744285  0.2614068  0.0000000
#2 TRAESCS1A02G015100.1  39.60877  11.20517  27.64765  10.1765429  22.53499  13.4224305 10.46664511  6.762340   5.701287  6.88073129  5.2759830  4.8097907
#3 TRAESCS1A02G033600.1  16.46916  48.48168  29.04061  55.4791777  16.76590  35.3503572 55.15638041 79.880298 106.097881 60.50789527 77.3873432 78.5085089
#pval         qval    log2FC
#1 3.139324e-91 1.965908e-86 12.283900
#2 1.926465e-04 5.639975e-03  1.081448
#3 9.396356e-05 3.381716e-03 -1.714267
write.csv(tmp, file="sleuth-DEG.pw.csv", row.names=FALSE)

##############[Optional] Volcano-plot
#preprare "whole" output (just repeat the steps)
fdr <- results_table_wt[,c(1:4)]
fdr$log2FC <- -log(exp(fdr$b),2)   #<<- Add Log2FC column
fdr <- subset(fdr, select=-c(b))   #<<- Remove b column
all <- merge(d, fdr, by="target_id", sort=T)
vc <- all[,c("target_id","qval","log2FC")]

####### Classify DEGs
v1 <- subset(vc, vc$log2FC >= 1 & vc$qval < 0.05)
v2 <- subset(vc, vc$log2FC <= -1 & vc$qval < 0.05)
deg_ids <- c(as.vector(v1$target_id), as.vector(v2$target_id))
`%notin%` <- Negate(`%in%`) #temporary designates "not-in" func
v3 <- vc[vc$target_id %notin% deg_ids,] ## non-DEG
up.deg <- length(v1[,1])
dn.deg <- length(v2[,2])

###### Transform data
vc$qval <- -log(vc$qval, 10)
vc$qval <- round(vc$qval, digits=1)
vc$log2FC <- round(vc$log2FC, digits=1)

###### Generate & count DEGs
v1 <- vc[vc$target_id %in% as.vector(v1$target_id),]
v2 <- vc[vc$target_id %in% as.vector(v2$target_id),]
v3 <- vc[vc$target_id %in% as.vector(v3$target_id),] 
####### Reduce spots
v1 <- unique(v1)
v2 <- unique(v2)
v3 <- unique(v3)

####### Drawing plot with solid limits
plot(v3$log2FC,v3$qval, pch=20, cex=0.4, col="grey70", ylim=c(0,350), xlim=c(-15,15), xlab="Log2FC", ylab="-log10FDR")
points(v1$log2FC,v1$qval, pch=20, cex=0.5, col="#ED0422")
points(v2$log2FC,v2$qval, pch=20, cex=0.5, col="#3498DB")
abline(v=0, lty=2)

####### Add numbers of DEG on the plot
text(10,340, up.deg, pos=4, cex=1.4, col="#ED0422")
text(10,320, dn.deg, pos=4, cex=1.4, col="#3498DB")




sessionInfo()
#R version 4.1.0 (2021-05-18)
#Platform: x86_64-apple-darwin17.0 (64-bit)
#Running under: macOS Big Sur 10.16

#Matrix products: default
#LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib

#locale:
#  [1] ja_JP.UTF-8/ja_JP.UTF-8/ja_JP.UTF-8/C/ja_JP.UTF-8/ja_JP.UTF-8

#attached base packages:
#  [1] stats4    parallel  stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#[1] GO.db_3.13.0           pathview_1.32.0        GOstats_2.58.0         Category_2.58.0        Matrix_1.4-0           GSEABase_1.54.0       
#[7] graph_1.70.0           annotate_1.70.0        XML_3.99-0.8           AnnotationForge_1.34.1 AnnotationDbi_1.54.1   IRanges_2.26.0        
#[13] S4Vectors_0.30.2       Biobase_2.52.0         BiocGenerics_0.38.0    KEGGREST_1.32.0        sleuth_0.30.0         

#loaded via a namespace (and not attached):
#[1] httr_1.4.2             tidyr_1.1.4            bit64_4.0.5            splines_4.1.0          assertthat_0.2.1       BiocManager_1.30.16   
#[7] RBGL_1.68.0            blob_1.2.2             GenomeInfoDbData_1.2.6 pillar_1.6.5           RSQLite_2.2.9          lattice_0.20-45       
#[13] glue_1.6.1             digest_0.6.29          XVector_0.32.0         colorspace_2.0-2       plyr_1.8.6             pkgconfig_2.0.3       
#[19] genefilter_1.74.1      zlibbioc_1.38.0        purrr_0.3.4            xtable_1.8-4           scales_1.1.1           tibble_3.1.6          
#[25] generics_0.1.1         farver_2.1.0           ggplot2_3.3.5          ellipsis_0.3.2         cachem_1.0.6           lazyeval_0.2.2        
#[31] survival_3.2-13        magrittr_2.0.1         crayon_1.4.2           KEGGgraph_1.52.0       memoise_2.0.1          fansi_1.0.2           
#[37] tools_4.1.0            data.table_1.14.2      org.Hs.eg.db_3.13.0    lifecycle_1.0.1        matrixStats_0.61.0     stringr_1.4.0         
#[43] Rhdf5lib_1.14.2        munsell_0.5.0          Biostrings_2.60.2      compiler_4.1.0         GenomeInfoDb_1.28.4    tinytex_0.36          
#[49] rlang_0.4.12           rhdf5_2.36.0           grid_4.1.0             RCurl_1.98-1.5         rhdf5filters_1.4.0     rstudioapi_0.13       
#[55] bitops_1.0-7           labeling_0.4.2         gtable_0.3.0           curl_4.3.2             DBI_1.1.2              reshape2_1.4.4        
#[61] R6_2.5.1               dplyr_1.0.7            fastmap_1.1.0          bit_4.0.4              utf8_1.2.2             Rgraphviz_2.36.0      
#[67] stringi_1.7.6          Rcpp_1.0.8             vctrs_0.3.8            png_0.1-7              tidyselect_1.1.1       xfun_0.29             
