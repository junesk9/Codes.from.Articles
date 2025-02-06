##### Hormonome + Metabolome PCA Dendrogram
##### For Plant Biotechnology
##### 10.06.24 June-Sik Kim
#####

library(dplyr)
library(RColorBrewer)
library(beeswarm)
library(pheatmap)
library(mixOmics)

set.seed(101)
setwd("~/Library/CloudStorage/Box-Box/SciData_SI/")
`%notin%` <- Negate(`%in%`) #temp. function to "not-in"

#### Metabolome + Hormonome data
dd <- read.csv("3_metabolome-data/meta-hormone.parsed.csv", header=T, row.names=1) 
dd <- subset(dd, dd$tx %in% c("h","d")) #[1]  96 476
#(ordering by desired)
dd$CellLine <- factor(dd$CellLine, levels=c("sr1","Os","Pn","Pb")) 
dd$tx <- factor(dd$tx, levels=c("h","d"))
dd <- dd[order(dd$CellLine, dd$tx), ]

meta <- dd[,c(1:3)] #[1] 96  3
ms <- dd[,c(4:445)] #[1]  96 442
hrm <- dd[,c(446:476)] #[1] 96 31

#### PCA + scatter
palette1 <- c("#FF80FF","#FF8278", "#7FB8FA", "#04CF51") #sr, os,pn,pb
#(1. ms)
cvars <-  apply(ms, 2, var)
ms2 <- ms[ ,names(cvars[cvars > 0])] #[1]  96 360

rpca <- prcomp(x=ms2, scale=T) 
PropVar <- summary(rpca)$importance[2,]
xlab = paste0("PC1: ",round(PropVar[1]*100,1), "%")
ylab = paste0("PC2: ",round(PropVar[2]*100,1), "%")
plot(x=rpca$x[,1],y=rpca$x[,2], xlab=xlab, ylab=ylab, pch=c(17,19)[meta$tx], 
     cex=1.5, col=palette1[meta$CellLine], main="PCA: 360/442 metabolites")
abline(h=0, lty=2)
abline(v=0, lty=2)
text(-10, 5, "SR", col=palette1[1])
text(-8, 5, "Os", col=palette1[2])
text(-10, 4, "Pn", col=palette1[3])
text(-8, 4, "Pb", col=palette1[4])
#save pdf as 6x6 in

#(2. hormone)
hrm[is.na(hrm)] <- 0 
cvars <-  apply(hrm, 2, var) #[1] 96 31
hrm2 <- hrm[ ,names(cvars[cvars > 0])] #[1] 96 29

rpca <- prcomp(x=hrm2, scale=T) 
PropVar <- summary(rpca)$importance[2,]
xlab = paste0("PC1: ",round(PropVar[1]*100,1), "%")
ylab = paste0("PC2: ",round(PropVar[2]*100,1), "%")
plot(x=rpca$x[,1],y=rpca$x[,2], xlab=xlab, ylab=ylab, pch=c(17,19)[meta$tx], 
     cex=1.5, col=palette1[meta$CellLine], main="PCA: 29/32 hormones")
abline(h=0, lty=2)
abline(v=0, lty=2)
text(6, 6, "SR", col=palette1[1])
text(6, 5.5, "Os", col=palette1[2])
text(6, 5, "Pn", col=palette1[3])
text(6, 4.5, "Pb", col=palette1[4])

#### heatmap


##### CK-PnPb biosynthesis genes heatmap
ck.tpm <- read.csv("6_Orthologues/CK.Orth-set.PbPn.csv", header=T, row.names=5)
ck.meta <- ck.tpm[,c(1:4)]
ck.tpm <- ck.tpm[,c(5:dim(ck.tpm)[2])]
ck.z <- t(scale(t(ck.tpm))) #[1] 61 48
ck.z2 <- na.omit(ck.z) #[1] 55 48

palette2 <-  hcl.colors(100, "Lajolla")
palette3 <- colorRampPalette(c("#0073C3", "white", "#EF4221"))(100)

pheatmap(ck.z2, cluster_cols = F, breaks=seq(-1,3, length.out=101),
         color=palette2, show_colnames = F,
         annotation_row = ck.meta[rownames(ck.z2),"SYMBOL",drop=FALSE],
         annotation_col = meta, show_rownames = FALSE)


pheatmap(ck.z2, cluster_cols = F, breaks=seq(-1,3, length.out=101),
         color=palette3, show_colnames = F,
         annotation_row = ck.meta[rownames(ck.z2),"SYMBOL",drop=FALSE],
         annotation_col = meta, show_rownames = FALSE)




sessionInfo()
# R version 4.3.1 (2023-06-16)
# Platform: x86_64-apple-darwin20 (64-bit)
# Running under: macOS Monterey 12.7.5
# 
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
# LAPACK: /Library/Frameworks/R.framework/Versions/4.3-x86_64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0
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
#   [1] mixOmics_6.26.0     lattice_0.22-6      MASS_7.3-60         BiocManager_1.30.23 beeswarm_0.4.0      RColorBrewer_1.1-3  dplyr_1.1.4        
# [8] ggplot2_3.5.1       pheatmap_1.0.12     openxlsx2_1.7       scales_1.3.0          


