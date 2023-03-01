##### GAPIT application to Arabidopsis Hirai-sensei data
##### 2022.07.13 
##### June-Sik Kim


library(GAPIT3)
set.seed(101)

###########################################################################
############################## GWAS
myG <- read.delim("NB220.tair9.gt.hmp.txt",header=FALSE)
#myY <- read.delim("input_trait.txt",header=TRUE)
#myY <- myY[,c(-2,-3)] ## remove #2 #3 columns since they are not phenotype data
myY <- read.csv("NB220.trait.csv",header=T)
#myY <- myY[,c(1,505)]

for (i in range(2:length(myY))){
  myGAPIT <- GAPIT( Y=myY[,c(1,i)],G=myG,
    kinship.cluster="average",
    kinship.group="Mean",
    kinship.algorithm="VanRaden",
    PCA.total=10,
    model=c("MLMM")
  )
}

########################################################################
############ heatmap
library(gplots)

d <- read.csv("Arabi-metabolome.gwas.log10p.100k-win.csv",header=T, row.names=1)
colors = c(seq(0,3,length=101),seq(3,10,length=101),seq(10,20,length=101))
colors = unique(colors) ## break vector should composite unique numbers
length(colors)
my_palette <- colorRampPalette(c("black", "grey20","yellow", "magenta"))(n = length(colors)-1) 

heatmap.2(as.matrix(d), col=my_palette, breaks=colors, 
          density.info="none", trace="none", dendrogram="none", 
          symm=F,symkey=F,symbreaks=T, scale="none", 
          Colv=FALSE, Rowv=FALSE, key=FALSE, 
          lwid=c(0.2,4), lhei=c(0.1,4))

##simple colorbar
plot(colors, col=my_palette)

#####################################################################
######################### QQMAN
library(qqman)

f_ls <- file.path("5_GWAS-every", list.files("5_GWAS-every","GAPIT.MLMM.Chem967"))
f_ls <- file.path("3_gwas(gapit3)_methylhistidine", list.files("5_GWAS-every","GAPIT.MLMM.Chem967"))
f_ls <- f_ls[grep("Results.csv", f_ls)]
#[1] "GAPIT.MLMM.Chem967_111.GWAS.Results.csv" "GAPIT.MLMM.Chem967_112.GWAS.Results.csv"
#[3] "GAPIT.MLMM.Chem967_113.GWAS.Results.csv" "GAPIT.MLMM.Chem967_114.GWAS.Results.csv"
#[5] "GAPIT.MLMM.Chem967_115.GWAS.Results.csv" "GAPIT.MLMM.Chem967_116.GWAS.Results.csv"

for (f in f_ls){
  prefix = unlist(strsplit(f, split="\\."))
  prefix = prefix[c(2:4)]
  prefix = paste(prefix,collapse=".")

  mht_f <- paste0(prefix,".MHT.pdf")
  qq_f <- paste0(prefix, ".QQ.pdf")

  d <- read.csv(f_ls[1], header=T)
  d_qq <- d[,c(1,2,3,4,9)]
  colnames(d_qq) <- c("SNP","CHR","BP","P","FDR")
  #1 5:15426918:C:T          5 15426918 3.833088e-62          8.199933e-57
  #2   2:497458:A:C          2   497458 1.778203e-36          1.902011e-31
  #3 4:10564911:C:T          4 10564911 8.696794e-25          6.201539e-20
  #4 2:13617193:C:T          2 13617193 2.097825e-23          1.121943e-18
  #5 2:13651441:T:A          2 13651441 3.148348e-14          1.347021e-09

  ######## set significance line
  sig <- d_qq[d_qq$FDR < 0.05,] 
  th = -log10(0.05/length(d_qq$P)) ## Bonferroni-corrected 0.05 [or, so-called Bonferroni threshold] same to the green solid line of GAPIT3 output 

  ######Lightening the data
  library(dplyr)
  d_slim = d_qq
  d_slim$BP <- round(d_slim$BP, digits=-6)
  d_slim$P <- 10^-(round(-log10(d_slim$P), digits=2))
  d_slim <- d_slim %>% distinct(d_slim$BP, d_slim$P, .keep_all = TRUE)
  ## 213925 -> 7418

  ###set graph height ylim
  ylim = max(-log10(d_slim$P))
  ylim = 10* ceiling(ylim/10)


  #####################drawing plot
  pdf(mht_f, w=4, h=2, pointsize=4, useDingbats=F)
  manhattan(d_slim, suggestiveline=FALSE, genomewideline=th, col=c("grey10","grey70"), cex=1.5, highlight=sig$SNP, ylim = c(0,ylim), xaxt='n', tck=0.02, xlab="", ylab="")
  ### Suggestiveline = 1e-5   genomewideline = 5e-8 (https://jojoshin.hatenablog.com/entry/2016/01/18/142622)
  ### [tck=0.02] to put label ticks inside
  dev.off()

  pdf(qq_f, w=1.8, h=2, pointsize=4, useDingbats=F)
  qq(d_slim$P, xaxt='n', tck=0.02, xlab="", ylab="",cex=0.8)
  dev.off()
}

sessionInfo()
#R version 4.1.3 (2022-03-10)
#Platform: x86_64-apple-darwin17.0 (64-bit)
#Running under: macOS Monterey 12.5

#Matrix products: default
#LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib

#locale:
#  [1] ja_JP.UTF-8/ja_JP.UTF-8/ja_JP.UTF-8/C/ja_JP.UTF-8/ja_JP.UTF-8

#attached base packages:
#  [1] compiler  stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#  [1] biganalytics_1.1.21  biglm_0.9-2.1        DBI_1.1.3            foreach_1.5.2        bigmemory_4.6.1      scatterplot3d_0.3-41
#[7] EMMREML_3.1          Matrix_1.4-1         ape_5.6-2            genetics_1.3.8.1.3   mvtnorm_1.1-3        MASS_7.3-58.1       
#[13] gtools_3.9.3         gdata_2.18.0.1       combinat_0.0-8       LDheatmap_1.0-6      gplots_3.1.3         multtest_2.48.0     
#[19] Biobase_2.54.0       BiocGenerics_0.40.0  GAPIT3_3.1.0         cowplot_1.1.1        patchwork_1.1.1      sp_1.5-0            
#[25] SeuratObject_4.1.0   Seurat_4.1.1         dplyr_1.0.9         
