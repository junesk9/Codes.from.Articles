##### LDHeatmap 
##### 2021.12.21
##### by June-Sik Kim (https://github.com/junesk9/)


library(LDheatmap)
library(snpStats)


dd <- read.table("SV274_Exome_SNPChip.1321776.bi.FILLIN.316345Filtered.hmp.txt", skip=1, row.names=1)
##since the first row was irregular

######SNP range
range = read.table(pipe("pbpaste")) ##"GAPIT.Blink.X20_Mildew.rI.Df.tValue.StdErr.csv"

snp = as.vector(range[,1])
bp = as.vector(range[,3])

######subseting
d <- subset(dd, rownames(dd) %in% snp)
ds <- d[,-c(1:10)]
dt <- t(ds)
dt[dt == "A"] <- "A/A"
dt[dt == "T"] <- "T/T"
dt[dt == "G"] <- "G/G"
dt[dt == "C"] <- "C/C"
dt[dt == "N"] <- NA
row.names(dt) <- NULL

for(i in 1:num){dt[,i] <- as.factor(dt[,i])}
df <- as.data.frame(dt)
for(i in 1:num){df[,i] <- as.integer(df[,i])}
gdat<-as(as.matrix(df),"SnpMatrix")

require(viridis)
LDheatmap(gdat , genetic.distances = bp, flip=TRUE, col=magma(20))













sessionInfo()
R version 4.1.0 (2021-05-18)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS Big Sur 10.16

Matrix products: default
LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib

locale:
  [1] ja_JP.UTF-8/ja_JP.UTF-8/ja_JP.UTF-8/C/ja_JP.UTF-8/ja_JP.UTF-8

attached base packages:
  [1] grid      stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
  [1] dplyr_1.0.7        genetics_1.3.8.1.3 mvtnorm_1.1-3      MASS_7.3-54        gtools_3.9.2       gdata_2.18.0      
[7] combinat_0.0-8     snpStats_1.42.0    Matrix_1.4-0       survival_3.2-13    LDheatmap_1.0-4   

loaded via a namespace (and not attached):
  [1] Rcpp_1.0.7          rstudioapi_0.13     magrittr_2.0.1      BiocGenerics_0.38.0 splines_4.1.0       zlibbioc_1.38.0    
[7] tidyselect_1.1.1    lattice_0.20-45     R6_2.5.1            rlang_0.4.12        fansi_0.5.0         tools_4.1.0        
[13] parallel_4.1.0      xfun_0.29           utf8_1.2.2          tinytex_0.35        DBI_1.1.1           ellipsis_0.3.2     
[19] assertthat_0.2.1    tibble_3.1.6        lifecycle_1.0.1     crayon_1.4.2        purrr_0.3.4         vctrs_0.3.8        
[25] glue_1.5.1          compiler_4.1.0      pillar_1.6.4        generics_0.1.1      pkgconfig_2.0.3    
