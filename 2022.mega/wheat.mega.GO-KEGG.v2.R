####Sleuth analysis of wheat RNA-seq timecourse (Mega)
####2022.01.26
####KIM JS

############################################################################
############ GO & KEGG analysis
library(GSEABase)
library(GOstats)
library(KEGGREST)
library(pathview)

## Load the KEGG ref file for later
dk <- read.table("wheat.KEGG.Ensembl49.txt", header=T, row.names=NULL, sep="\t")
dk[,2] <- substr(dk[,2],1,5)
dk <- unique(dk[,c(2,1)]) #ほけん
head(dk)
#KEGG.Pathway.and.Enzyme.ID     Gene.stable.ID
#1                      00010 TraesCS6A02G386600
#2                      00010 TraesCS6B02G425700
#3                      00010 TraesCS6D02G371200
keggframe <- KEGGFrame(dk, organism="wheat")
gsck <- GeneSetCollection(keggframe, setType=KEGGCollection())

## Load a GO file to a custum GO DB
go <- read.table("wheat.GOSlim.Ensembl52.txt", header=T, sep="\t")
go$evi = "IEA"
go <- go[,c(3,4,1)]
head(go)
#  GOSlim.GOA.Accession.s. evi     Gene.stable.ID
#1              GO:0003674 IEA TraesCS3A02G154800
#2              GO:0003824 IEA TraesCS3A02G154800
#3              GO:0016787 IEA TraesCS3A02G154800
goFrame <- GOFrame(go, organism="wheat")
goAllFrame <- GOAllFrame(goFrame)
gsc <- GeneSetCollection(goAllFrame, setType = GOCollection())
all <- as.vector(go$Gene_id)  ## a vector of all referenced list

############## MOVE TO A FOLDER OF DEG RESULT (for pairwise-DEG)
##### Prepare DEG
f <- read.table("sleuth-DEG.pw.csv",header=T,sep=",")
up <- subset(f,f$log2FC > 0)
dn <- subset(f,f$log2FC < 0)

###process deg list
deg <- up
deg[,1] <- gsub(".1","",deg[,1], fixed=TRUE)
deg[,1] <- gsub("TRAESCS","TraesCS",deg[,1]) #ほけん
up <- as.vector(deg[,1])
deg <- dn
deg[,1] <- gsub(".1","",deg[,1], fixed=TRUE)
deg[,1] <- gsub("TRAESCS","TraesCS",deg[,1]) #ほけん
dn <- as.vector(deg[,1])

###### Run GOStat
for (f in c("MF","CC","BP")){
  outf = paste("GO.UP",f,sep=".")
  outf = paste(outf,"txt",sep=".")
  p <- GSEAGOHyperGParams(
    name = "Paramaters",
    geneSetCollection = gsc,
    geneIds = up,
    universeGeneIds = all,
    ontology = f,
    pvalueCutoff = 1,
    conditional = FALSE,
    testDirection = "over"
  )
  result <- hyperGTest(p)
  
  ### Check and output the result
  head(summary(result)) # Check the result
  write.table(summary(result), file = outf, quote=FALSE,sep="\t",row.names=FALSE)
}

for (f in c("MF","CC","BP")){
  outf = paste("GO.DN",f,sep=".")
  outf = paste(outf,"txt",sep=".")
  p <- GSEAGOHyperGParams(
    name = "Paramaters",
    geneSetCollection = gsc,
    geneIds = dn,
    universeGeneIds = all,
    ontology = f,
    pvalueCutoff = 1,
    conditional = FALSE,
    testDirection = "over"
  )
  result <- hyperGTest(p)
  
  ### Check and output the result
  head(summary(result)) # Check the result
  write.table(summary(result), file = outf, quote=FALSE,sep="\t",row.names=FALSE)
}

###### Run KEGG
kegg.p <- GSEAKEGGHyperGParams(
  name = "Paramaters",
  geneSetCollection = gsck,
  geneIds = up, ## DEG list
  universeGeneIds = all,
  pvalueCutoff = 1,
  testDirection = "over"
)
kegg.result <- hyperGTest(kegg.p)
head(summary(kegg.result)) # Check the result
##### Write-down the result table
write.table(summary(kegg.result), file = "KEGG.UP.txt",quote=FALSE,sep="\t",row.names=FALSE)

kegg.p <- GSEAKEGGHyperGParams(
  name = "Paramaters",
  geneSetCollection = gsck,
  geneIds = dn, ## DEG list
  universeGeneIds = all,
  pvalueCutoff = 1,
  testDirection = "over"
)
kegg.result <- hyperGTest(kegg.p)
head(summary(kegg.result)) # Check the result
##### Write-down the result table
write.table(summary(kegg.result), file = "KEGG.DN.txt",quote=FALSE,sep="\t",row.names=FALSE)


############################################################### (Optional) view KEGG pathway
gene.id <- up
path.id <- summary(kegg.result)$KEGGID[1]
pv <- pathview(gene.data = gene.id, pathway.id = path.id, species = "osa", gene.idtype = "KEGG")
## KEGG.species for bread wheat "taes" is not working now...





####################################### For clustered outputs
fc_df <- read.table("cclust.timecourse.mean.csv", header=T, row.names=1, sep=",")


#### Prepare vectors by clusters namely degC1, degC2, degC3, ...
for(i in 1:6){
  nameA <- paste("C",i,sep="")
  tmp <- row.names(subset(fc_df, fc_df$FC_k6==i))
  tmp <- gsub(".1","", tmp, fixed=TRUE)
  tmp <- gsub("TRAESCS","TraesCS",tmp) #ほけん
  tmp <- as.vector(tmp)
  assign(nameA, tmp)
  
  for (f in c("MF","CC","BP")){
    outf = paste("GO",f,nameA,sep=".")
    outf = paste(outf,"txt",sep=".")
    p <- GSEAGOHyperGParams(
      name = "Paramaters",
      geneSetCollection = gsc,
      geneIds = tmp,
      universeGeneIds = all,
      ontology = f,
      pvalueCutoff = 0.05,
      conditional = FALSE,
      testDirection = "over"
    )
    result <- hyperGTest(p)
    
    ### Check and output the result
    head(summary(result)) # Check the result
    write.table(summary(result), file = outf, quote=FALSE,sep="\t",row.names=FALSE)
  }
  
  kegg.p <- GSEAKEGGHyperGParams(
    name = "Paramaters",
    geneSetCollection = gsck,
    geneIds = tmp, ## DEG list
    universeGeneIds = all,
    pvalueCutoff = 0.05,
    testDirection = "over"
  )
  kegg.result <- hyperGTest(kegg.p)
  head(summary(kegg.result)) # Check the result
  ##### Write-down the result table
  outf2 = paste("KEGG",nameA,sep=".")
  outf2 = paste(outf2,"txt",sep=".")
  write.table(summary(kegg.result), file = outf2, quote=FALSE,sep="\t",row.names=FALSE)
}

###[Optional] findout genes in the certain categories
cat_deg <- geneIdsByCategory(kegg.result)$'00500'
cat_tot <- geneIdUniverse(kegg.result)$'00500'



############################################################################
######################## Z-standardization & heapmap
library(pheatmap)
library(viridis) ##for a color palette
library(factoextra) ##for elbow plots finding-out optimal clusters
library(cluster) ## for Clusgap() method

#### Excel-prepared -log10 P-values 
d <- read.table(pipe("pbpaste"), header=T, row.names=1)
head(d, c(5,5))
#Day.0_UP Day.1.DR_UP Day.3.DR_UP Day.4.DR_UP Day.1.WTR_UP
#10 0.005768445   0.3795943   2.1653699 2.187600380   0.00000000
#20 0.201093919   0.4457786   0.4469951 0.820593590   0.00000000
#30 0.000000000   0.1390185   0.7507784 1.548575064   0.06179596
#40 4.753535193   4.7234162   0.7800753 0.005574549   3.52809615
#51 0.000000000   0.3531578   0.6155939 1.057759594   0.00000000
d1 <- d[,c(1:7)] ##UP-deg, L8 vs. ND data
d2 <- d[,c(1,8:13)]

#Z-standardization
d1z <- t(apply(d1, 1, scale)) #need transpose
colnames(d1z) <- colnames(d1) #z-scaling losing colnames
d1z <- na.omit(d1z) #127 rows -> 106 rows

d2z <- t(apply(d2, 1, scale))
colnames(d2z) <- colnames(d2)
d2z <- na.omit(d2z) #127 rows -> 103 rows

## visualization
#different colorscheme like magma(100), viridis(100) also should be considered
palette <- colorRampPalette(c("blue", "gray88","gray88","gray88", "red"))(100) # for custom palette
pheatmap(d1z, color = inferno(100), cluster_cols = F, clustering_method = "ward.D2", 
         clustering_distance_rows = "correlation", cutree_rows = 8)
#the plot saved as PDF as a portrait (12 x 6 in) for KEGG
#the plot saved as PDF as a portrait (10 x 6 in) for GO-MF



###### Find-out Optimal clusters
##https://www.r-bloggers.com/2019/01/10-tips-for-choosing-the-optimal-number-of-clusters/
#elbow plot 
fviz_nbclust(d1z, kmeans, method = "wss", k.max = 24) + theme_minimal() + ggtitle("the Elbow Method")
#silhouette plot
fviz_nbclust(e1z, kmeans, method = "silhouette", k.max = 24) 
+ theme_minimal() + ggtitle("The Silhouette Plot")
#gap stats
gap_stat <- clusGap(d1z, FUN = kmeans, nstart = 30, K.max = 24, B = 50)
fviz_gap_stat(gap_stat) + theme_minimal() + ggtitle("fviz_gap_stat: Gap Statistic")




sessionInfo()
#R version 4.1.3 (2022-03-10)
#Platform: x86_64-apple-darwin17.0 (64-bit)
#Running under: macOS Monterey 12.6

#Matrix products: default
#LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib

#locale:
#  [1] ja_JP.UTF-8/ja_JP.UTF-8/ja_JP.UTF-8/C/ja_JP.UTF-8/ja_JP.UTF-8

#attached base packages:
#  [1] stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#[1] cluster_2.1.4     factoextra_1.0.7  clustree_0.5.0    ggraph_2.0.6      ggplot2_3.3.6    
#[6] viridis_0.6.2     viridisLite_0.4.1 pheatmap_1.0.12  

#loaded via a namespace (and not attached):
#  [1] Rcpp_1.0.9         ggpubr_0.4.0       pillar_1.8.1       compiler_4.1.3     RColorBrewer_1.1-3
#[6] tools_4.1.3        digest_0.6.29      lifecycle_1.0.1    tibble_3.1.8       gtable_0.3.0      
#[11] pkgconfig_2.0.3    rlang_1.0.4        tidygraph_1.2.2    igraph_1.3.4       cli_3.3.0         
#[16] DBI_1.1.3          rstudioapi_0.14    ggrepel_0.9.1      gridExtra_2.3      withr_2.5.0       
#[21] dplyr_1.0.9        graphlayouts_0.8.1 generics_0.1.3     vctrs_0.4.1        grid_4.1.3        
#[26] tidyselect_1.1.2   glue_1.6.2         R6_2.5.1           rstatix_0.7.0      fansi_1.0.3       
#[31] polyclip_1.10-0    carData_3.0-5      car_3.1-0          farver_2.1.1       tweenr_2.0.1      
#[36] purrr_0.3.4        tidyr_1.2.0        magrittr_2.0.3     backports_1.4.1    scales_1.2.1      
#[41] MASS_7.3-58.1      abind_1.4-5        assertthat_0.2.1   ggforce_0.3.4      colorspace_2.0-3  
#[46] ggsignif_0.6.3     labeling_0.4.2     utf8_1.2.2         munsell_0.5.0      broom_1.0.0      