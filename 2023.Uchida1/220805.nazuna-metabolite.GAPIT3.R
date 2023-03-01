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
