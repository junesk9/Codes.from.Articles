#############################################################
########### Seurat output visualization & characterization
######################### Junesk9 2025.09.25

############################
######## Load libraries required
###########################

library(Seurat) 


#emulate the ggplot hue() color palette
ggcolor <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
require(scales)
show_col(ggcolor(8)) #test ggcolor
key.palette <- hcl.colors(100, "Cividis")
key.palette <- rev(hcl.colors(100, "RdYlBu"))
show_col(key.palette)

setwd("/Users/junesk9/理化学研究所　セルロース生産研究チーム Dropbox/June-Sik Kim/オオムギ圃場mRNA.ChIPseq")
set.seed(101)
`%notin%` <- Negate(`%in%`) #temp. function to "not-in"
###############################
##### Load the data
################################

#1. Basal Expression data
so <- readRDS("seurat.barley-field.leaf.v3.rds")

meta <- so@meta.data #[1] 1940   35
all.rpm <- as.data.frame(so@assays$RNA@data) #[1] 25145  1940
cl45deg <- read.csv("02_cluster-DEG/cl4n5-pos.markers.MAST.csv", header=T, row.names = 1)
cl45deg <- rownames(cl45deg[cl45deg$p_val_adj < 0.05 & cl45deg$avg_log2FC >= 0.5, ])  #435


# [1] "orig.ident"          "nCount_RNA"          "nFeature_RNA"        "CytoTRACE"           "Place"               "Accession"           "Season"              "Batch"              
# [9] "Date"                "Date.2"              "Date.3"              "fD2H"                "month"               "week"                "RNA_snn_res.0.45"    "seurat_clusters"    
# [17] "Lineage1"            "Lineage2"            "Lineage0"            "LineageM"            "AvgTemp"             "LightHour"           "DayHour"             "J247_UMAP1"         
# [25] "J247_UMAP2"          "J064_UMAP1"          "J064_UMAP2"          "J647_UMAP1"          "J647_UMAP2"          "H602_UMAP1"          "H602_UMAP2"          "RoomJ247wk.VF"      
# [33] "RoomJ247wk.FLOR"     "RoomJ247wk.ALL"      "RoomJ247wk.FLOR.PCC" "RoomJ247wk.VF.PCC"   "CytoTRACE.VF"        "1wk.indoor_PCC"      "2wk.indoor_PCC"      "3wk.indoor_PCC"     
# [41] "4wk.indoor_PCC"      "5wk.indoor_PCC"      "6wk.indoor_PCC"      "7wk.indoor_PCC"      "8wk.indoor_PCC"      "z_ABA"               "z_IAA"               "z_iP"               
# [49] "z_JA"                "z_JAIle"             "z_SA"                "z_tZ"                "IAA"                 "ABA"                 "iP"                  "JA"                 
# [57] "JAIle"               "SA"                  "tZ"                  "lTemp"               "hTemp"               "SAno18KIBR"          "z_SAno18KIBR"        "LNI"        

####################### 
######### tradeseq fitGAM
#########################
library(tradeSeq)
curve <- readRDS("seurat-curves.barley-field.leaf.v2.rds")
counts <- as.matrix(so@assays$RNA@counts)
meta2 <- t(meta[colnames(counts),c("fD2H","AvgTemp","DayHour","LineageM")])
meta2["AvgTemp",] <- meta2["AvgTemp",] + 1 #prevent the minus number
meta2["LineageM",] <- meta2["LineageM",] + 1
counts <- rbind(counts, meta2)
# BPPARAM <- BiocParallel::MulticoreParam(8) 
# 
curve@assays@data$pseudotime[,1] <- meta$fD2H
na.cell <- !is.na(curve@assays@data$pseudotime)
cw <- as.matrix(curve@assays@data$weights[na.cell,])
pstime <- as.matrix(curve@assays@data$pseudotime[na.cell,])
count <-counts[, rownames(cw)]
#
#(optional) conditional (parallel) 
accs.cond <- meta[colnames(count), "Accession"]
accs.cond <-factor(accs.cond, level=c("J247", "J064", "H602","J647"))
levels(accs.cond) <- c(1,2,2,2)
#
# 
# #2. choose best-K for fitGAM(): take minutes ~ hours
# icMat <- evaluateK(counts=count, pseudotime=pstime, cellWeight=cw, k=3:15, nGenes=200,verbose=F, plot=T, parallel=T, BPPARAM=BPPARAM)
# bestK = 6 #fD2H, LineageM
# bestK = 7 #Lineage0, AvgTemp, DayHour
# bestK = 8 #CytoTRACE
# 
# #3. fitGAM()
# gam = fitGAM(counts=count, pseudotime=pstime, cellWeight=cw, nknots=bestK, parallel=T, BPPARAM=BPPARAM)
# gam.assoc = associationTest(gam, lineages=T, global=T, l2fc=1)
# gam.assoc <- gam.assoc[order(gam.assoc$pvalue, -gam.assoc$meanLogFC), ]
# dim(subset(gam.assoc, gam.assoc$pvalue < 0.001 & abs(gam.assoc$meanLogFC) >= 1))[1]
# 
# #4. Smoothering the exp-data on the GAM
# preSm <- predictSmooth(gam, gene=rownames(count), nPoints=100, tidy=FALSE) #tidy=TRUE for x-axis values
# write.csv(preSm, "03_TradeSeq-assoc/GAM-AvgTemp.All.100win.csv")

#3a. fitGAM()-conditional
gam.j247 = fitGAM(counts=count, pseudotime=pstime, cellWeight=cw, conditions=accs.cond, nknots=bestK, parallel=F)
preSm <- predictSmooth(gam.j247, gene=rownames(count), nPoints=100, tidy=FALSE) 
preSm2 <- predictSmooth(gam.j247, gene=rownames(count), nPoints=100, tidy=TRUE)
preSm2 <- subset(preSm2, preSm2$gene == "HORVU0Hr1G000020")
write.csv(preSm, "03_TradeSeq-assoc/GAM-fD2H.J2vs.100win.csv")
write.csv(preSm2, "03_TradeSeq-assoc/GAM-fD2H.J2vs.100win2.csv")


gam_temp <- read.csv("03_TradeSeq-assoc/GAM-AvgTemp.All.100win.csv", header=T, row.names=1)
gam_dayh <- read.csv("03_TradeSeq-assoc/GAM-DayHour.All.100win.csv", header=T, row.names=1)
gam_cyto <- read.csv("03_TradeSeq-assoc/GAM-CytoTRACE.All.100win.csv", row.names=1, header=T)
gam_d2h <- read.csv("03_TradeSeq-assoc/GAM-fD2H.All.100win.csv", row.names=1, header=T)

gam_temp2 <- read.csv("03_TradeSeq-assoc/GAM-AvgTemp.J2vs.100win.csv", header=T, row.names=1)
gam_dayh2 <- read.csv("03_TradeSeq-assoc/GAM-dayHour.J2vs.100win.csv", header=T, row.names=1)

#5. normalization for visualization
library(cluster) #for clusGAP
library(gplots) #for heatmap.2
library(pheatmap)

gam <- gam_temp
preSm <- gam[rownames(gam) %in% ps431, ] #c(kegg.phosyn, kegg.calvin)n cl45deg kg143 ps431
#(Normalization)
# x1 <- 100* preSm/ apply(preSm, 1, sum) #converted to the percentage
# x2 <- 100* preSm / apply(preSm, 1, function(y) sum(abs(y))) # L1 normalisation
# x3 <- preSm - apply(preSm, 1, mean) # geneMean, centralization
x4 <- (preSm - apply(preSm, 1, mean) ) / apply(preSm, 1, sd)  # Z-standardization simply x4=t(scale(t(preSm)))
mtx <- x4
pheatmap(mtx, cluster_cols = FALSE, cluster_rows=TRUE, show_colnames=FALSE, show_rownames = FALSE, main="rbcS n=26",breaks = seq(-2, 2, length.out=101))

#(K-mer clustering)
cg.res <- clusGap(mtx, kmeans, K.max=20, B=100, verbose=interactive())
nClus <- 6 #6: gam_temp 8:gam_dayh 5:gam_d2h #ps431
           #8: gam_temp                      #kg143
           #5: gam_temp                      #kegg.ps + calvin
clus <- kmeans(x = mtx, centers = nClus, iter.max = 50)
centers <- clus$centers-apply(clus$centers, 1, mean)
hclus <- hclust(as.dist(1-cor(t(centers), method = "pearson")), method = "average")
nOrder <- match(clus$cluster, hclus$order) 
Kmeans_matrix <- mtx[order(nOrder),]
Kmeans_cluster <- sort(nOrder)
# table(Kmeans_cluster)

Kmeans_matrix2 <- as.matrix(Kmeans_matrix) - apply(Kmeans_matrix, 1, mean) #centralization
#(polishing the view by removing the out-liners)
Kmeans_matrix2 <- as.matrix(Kmeans_matrix2)
cutoff <- median(unlist(Kmeans_matrix2)) + 3*sd (unlist(Kmeans_matrix2)) #remove median +3*sd
Kmeans_matrix2[Kmeans_matrix2 >cutoff] <- cutoff 
cutoff <- median(unlist(Kmeans_matrix2)) - 3*sd (unlist(Kmeans_matrix2)) #remove median +3*sd
Kmeans_matrix2[Kmeans_matrix2 < cutoff] <- cutoff

Kmeans_cluster2 <- factor(Kmeans_cluster, levels=unique(Kmeans_cluster))
levels(Kmeans_cluster2) <- c(1,2,6,3,4,5) #AvgTemp
#levels(Kmeans_cluster2) <- c(8,7,6,3,5,4,2,1) #DayHour
#levels(Kmeans_cluster2) <- c(5,4,3,1,2) #D2H
#levels(Kmeans_cluster2) <- c(1,2,8,4,3,7,5,6)
#levels(Kmeans_cluster2) <- c(2,5,1,4,3) #Kegg-avgtemp
rOrder <- order(as.vector(Kmeans_cluster2))
Kmeans_matrix2 <- Kmeans_matrix2[rOrder, ]
Kmeans_cluster2 <- factor(Kmeans_cluster2[rOrder], levels=c(1:nClus))

#(pheatmap.annt)
ngenes <- dim(Kmeans_matrix2)[1]
annt_row <- data.frame(row.names=rownames(Kmeans_matrix2), cl45=rep(1, times=ngenes), clus=Kmeans_cluster2, calvin=rep(1, times=ngenes))
annt_row[!rownames(annt_row) %in% cl45deg, "cl45"] <- NA
annt_row[!rownames(annt_row) %in% kegg.calvin, "calvin"] <- NA

label_row = rownames(Kmeans_matrix2)
label_row[!label_row %in% kegg.ps$Gene_ID[kegg.ps$EC %in% c("4.1.1.39")]]  <- NA #rubisco
label_row[!is.na(label_row)] <- "*"
label_row[is.na(label_row)] <- ""

# kegg.ps <- kegg.ps[kegg.ps$Gene_ID %in% rownames(Kmeans_matrix2), ] #119 6
# rownames(kegg.ps) <- kegg.ps$Gene_ID
# annt_row$stage1 <- kegg.ps[rownames(Kmeans_matrix2), "ST1"]
# annt_row$stage2 <- kegg.ps[rownames(Kmeans_matrix2), "ST2"]
# annt_row$stage3 <- kegg.ps[rownames(Kmeans_matrix2), "ST3"]
# annt_row$stage4 <- kegg.ps[rownames(Kmeans_matrix2), "ST4"]


key.palette <- rev(hcl.colors(100, "RdYlBu"))
pheatmap(Kmeans_matrix2, cluster_cols = FALSE, cluster_rows=FALSE, labels_row=label_row,
         show_colnames=FALSE, color = key.palette, breaks = seq(-3, 3, length.out=101),
         annotation_row = annt_row, main = paste0("go-photosynthesis on AvgTemp smoothered; n=",ngenes))


#######################
#### GO/KEGG
#######################
library(GSEABase) 
library(GOstats) #GO/KEGG
library(KEGGREST) #KEGG-API

#2. GO/KEGG daata
gsc <- readRDS("Barley-IBSCv2.GO-gsc.rds")
all <- unique(unlist(geneIds(gsc)))
ps431 <- geneIds(gsc[["GO:0015979"]])
lr251 <- geneIds(gsc[["GO:0009416"]])
st604 <- geneIds(gsc[["GO:0009628"]])
en786 <- geneIds(gsc[["GO:0006091"]])
lh41 <- geneIds(gsc[["GO:0009765"]]) #41
psI <- unique(c(geneIds(gsc[["GO:0009768"]]), geneIds(gsc[["GO:0009769"]])))


gsck <- readRDS("Barley-IBSCv2.KEGG-gsc.rds")
kegg.all <- unique(unlist(geneIds(gsck)))
kegg.calvin <- geneIds(gsck[["00710"]]) #[1] 87
kegg.phosyn <- geneIds(gsck[["00195"]]) #[1] 51
kegg.gylox <- geneIds(gsck[["00630"]]) #[1] 71

kegg.annt <- keggList(database="pathway", organism="taes")
kegg.annt <- data.frame(KEGGid=names(kegg.annt), desc=as.vector(kegg.annt))
kegg.ps <- read.csv("00_IPSR_KIBR_Sampledata/Barley.IBSCv2.ensembl49.KEGG-PS-gids.csv", header=T) 
kg143 <- unique(kegg.ps[,1])
rubisco <- kegg.ps[kegg.ps$EC %in% c("4.1.1.39"), "Gene_ID"]
phorec <- kegg.ps[kegg.ps$EC %in% c("7.1.2.2","7.1.1.6"), "Gene_ID"]
ec_all <- kegg.ps[kegg.ps$EC %in% ec, "Gene_ID"] #147


runGO2 <- function(gid_ls){
  out_ls = list()
  p <- GSEAGOHyperGParams(
      name = "Paramaters",
      geneSetCollection = gsc,
      geneIds = gid_ls,
      universeGeneIds = all,
      ontology = "BP",
      pvalueCutoff = 1,
      conditional = FALSE,
      testDirection = "over"
    )
    res <- hyperGTest(p)
    res_gid <- geneIdsByCategory(res)
    res2 <- summary(res)
    res2$GeneId <- sapply(res2[,1], function(x)paste(unlist(res_gid[x]), collapse="/"))
    out_ls = c(out_ls, list(res2))
  
  kegg.p <- GSEAKEGGHyperGParams(
    name = "Paramaters",
    geneSetCollection = gsck,
    geneIds = gid_ls, 
    universeGeneIds = kegg.all,
    pvalueCutoff = 1,
    testDirection = "over"
  )
  
  kegg.res <- hyperGTest(kegg.p)
  kegg.res_gid <- geneIdsByCategory(kegg.res)
  kegg.res2 <- summary(kegg.res)
  kegg.res2$gid <- sapply(kegg.res2[,1], function(x)paste(unlist(kegg.res_gid[x]), collapse="/"))
  kegg.res2$KEGGID <- sapply(kegg.res2$KEGGID, function(x) paste0("taes",x))
  kegg.res2$Term <- kegg.annt[kegg.res2$KEGGID,]$desc
  out_ls = c(out_ls, list(kegg.res2))
  
  return(out_ls)
} 


#######################
#### Figure 2
#######################
library(dunn.test)
library(ggpubr)
library(overlapping) 

############# U-test Fig 2C
shapiro.test(meta[meta$seurat_clusters %in% "5", "DayHour"])
# W = 0.87843, p-value = 1.51e-12 #non-parametric
shapiro.test(meta[meta$seurat_clusters %in% "4", "DayHour"])
# W = 0.76627, p-value < 2.2e-16 #non-parametric
var.test(x=meta[meta$seurat_clusters %in% "5", "DayHour"], y=meta[meta$seurat_clusters %in% "4", "DayHour"])
# F = 3.4121, num df = 227, denom df = 274, p-value < 2.2e-16  #homovariance
overlap(list(meta[meta$seurat_clusters %in% "5", "DayHour"], meta[meta$seurat_clusters %in% "4", "DayHour"]))
# [1] 0.4802114 #half-overlap
w <- wilcox.test(x=meta[meta$seurat_clusters %in% "5", "DayHour"], y=meta[meta$seurat_clusters %in% "4", "DayHour"], correct=FALSE) 
# W = 42696, p-value = 2.444e-12 #correct=TRUE
# W = 42696, p-value = 2.438e-12 #correct=FALSE

#(a bit more for effect size Z/√N) 
z <- qnorm(1-(w$p.value/2))
#[1] 7.006496
r <- z/sqrt(length(meta[meta$seurat_clusters %in% c(4,5), "DayHour"])/2)
#[1] 0.4418063

############# Multiple-test Fig 2C
kw <- kruskal.test(meta$DayHour~meta$seurat_clusters)
#Kruskal-Wallis chi-squared = 1191.4, df = 6, p-value < 2.2e-16
#(effect value eta2)
eta2 <- kw$statistic/(length(meta[,1])-1)
#                  0.6144469

#Bonferonni
bf <- pairwise.wilcox.test(meta$DayHour, meta$seurat_clusters, p.adj="bonferroni", exact=F)
# 5.1e-11
hm <- pairwise.wilcox.test(meta$DayHour, meta$seurat_clusters, p.adj="holm", exact=F) 
# 1.7e-11


############# U-test Fig 2D
shapiro.test(meta[meta$seurat_clusters %in% "5", "AvgTemp"])
# W = 0.85742, p-value = 1.022e-13 #non-parametric
shapiro.test(meta[meta$seurat_clusters %in% "4", "AvgTemp"])
# W = 0.92254, p-value = 9.264e-11 #non-parametric
var.test(x=meta[meta$seurat_clusters %in% "5", "AvgTemp"], y=meta[meta$seurat_clusters %in% "4", "AvgTemp"])
# F = 0.661, num df = 227, denom df = 274, p-value = 0.001271  #homovariance
overlap(list(meta[meta$seurat_clusters %in% "4", "AvgTemp"], meta[meta$seurat_clusters %in% "5", "AvgTemp"]))
# [1] 0.7280675 # largely overlap
w <- wilcox.test(x=meta[meta$seurat_clusters %in% "5", "AvgTemp"], y=meta[meta$seurat_clusters %in% "4", "AvgTemp"], correct=TRUE)
# W = 35924, p-value = 0.004774
z <- qnorm(1-(w$p.value/2))
# [1] 2.821868
r <- z/sqrt(length(meta[meta$seurat_clusters %in% c(4,5), "AvgTemp"])/2)
# [1] 0.1779376

############# Multiple-test Fig 2D
kw <- kruskal.test(meta$AvgTemp~meta$seurat_clusters)
#Kruskal-Wallis chi-squared = 1144.1, df = 6, p-value < 2.2e-16
#(effect value eta2)
eta2 <- kw$statistic/(length(meta[,1])-1)
#                  0.5900369 

#Bonferonni
bf <- pairwise.wilcox.test(meta$AvgTemp, meta$seurat_clusters, p.adj="bonferroni", exact=F)
# 0.10 
hm <- pairwise.wilcox.test(meta$AvgTemp, meta$seurat_clusters, p.adj="holm", exact=F) 
# 0.019
bh <- dunn.test(meta$AvgTemp, meta$seurat_clusters, method = "BH")
#0.0480

############# Multiple-testFig 2E
kw <- kruskal.test(meta$CytoTRACE~meta$seurat_clusters)
# Kruskal-Wallis chi-squared = 1331.9, df = 6, p-value < 2.2e-16
eta2 <- kw$statistic/(length(meta[,1])-1)
# 0.6869142 

#Bonferonni 4vs5
bf <- pairwise.wilcox.test(meta$CytoTRACE, meta$seurat_clusters, p.adj="bonferroni", exact=F)
# 1.3e-11
hm <- pairwise.wilcox.test(meta$CytoTRACE, meta$seurat_clusters, p.adj="holm", exact=F) 
# 3.0e-12
overlap(list(meta[meta$seurat_clusters %in% "4", "CytoTRACE"], meta[meta$seurat_clusters %in% "5", "CytoTRACE"]))
# [1] 0.6663155

w <- wilcox.test(meta[meta$seurat_clusters %in% "4", "CytoTRACE"], meta[meta$seurat_clusters %in% "5", "CytoTRACE"], correct=FALSE)
#W = 19666, p-value = 6.021e-13
z <- qnorm(1-(w$p.value/2))
#[1] 7.199987


#######################
#### Figure 3
#######################
cl23deg <- read.csv("02_cluster-DEG/cl2vs3.markers.MAST.csv", header=T, row.names = 1)
cl2v3 <- rownames(cl23deg[cl23deg$p_val_adj < 0.05 & cl23deg$avg_log2FC >= 0.5, ])  #4722
cl3v2 <- rownames(cl23deg[cl23deg$p_val_adj < 0.05 & cl23deg$avg_log2FC <= -0.5, ]) #156
go2v3 <- runGO2(cl2v3)
go3v2 <- runGO2(cl3v2)
write.csv(go2v3[[1]], "02_cluster-DEG/cl2vs3.markers.MAST.GO-BP.csv", row.names=F)
write.csv(go2v3[[2]], "02_cluster-DEG/cl2vs3.markers.MAST.KEGG.csv", row.names=F)
write.csv(go3v2[[2]], "02_cluster-DEG/cl3vs2.markers.MAST.KEGG.csv", row.names=F)
write.csv(go3v2[[1]], "02_cluster-DEG/cl3vs2.markers.MAST.GO-BP.csv", row.names=F)

cl67deg <- read.csv("02_cluster-DEG/cl6vs7.markers.MAST.csv", header=T, row.names = 1)
cl6v7 <- rownames(cl67deg[cl67deg$p_val_adj < 0.05 & cl67deg$avg_log2FC >= 0.5, ])  #[1] 2032
cl7v6 <- rownames(cl67deg[cl67deg$p_val_adj < 0.05 & cl67deg$avg_log2FC <= -0.5, ]) #[1] 6787
go6v7 <- runGO2(cl6v7)
go7v6 <- runGO2(cl7v6)
write.csv(go6v7[[1]], "02_cluster-DEG/cl6vs7.markers.MAST.GO-BP.csv", row.names=F)
write.csv(go6v7[[2]], "02_cluster-DEG/cl6vs7.markers.MAST.KEGG.csv", row.names=F)
write.csv(go7v6[[2]], "02_cluster-DEG/cl7vs6.markers.MAST.KEGG.csv", row.names=F)
write.csv(go7v6[[1]], "02_cluster-DEG/cl7vs6.markers.MAST.GO-BP.csv", row.names=F)



gam_temp2 <- read.csv("03_TradeSeq-assoc/GAM-AvgTemp.J2vs.100win.csv", header=T, row.names=1)
gam_dayh2 <- read.csv("03_TradeSeq-assoc/GAM-dayHour.J2vs.100win.csv", header=T, row.names=1)
gam_pst2 <- read.csv("03_TradeSeq-assoc/GAM-PStime.J2vs.100win.csv", header=T, row.names=1)
gam_d2h2 <- read.csv("03_TradeSeq-assoc/GAM-fD2H.J2vsOther.100win.csv", header=T, row.names=1)

gam <- gam_temp2
preSm <- gam[rownames(gam) %in% psI, ] #c(kegg.phosyn, kegg.calvin)n cl45deg kg143 ps431 c(kegg.calvin, kegg.gylox), lr251
#(Normalization)
x4 <- (preSm - apply(preSm, 1, mean) ) / apply(preSm, 1, sd)  # Z-standardization simply x4=t(scale(t(preSm)))
mtx <- x4
pheatmap(mtx, cluster_cols = FALSE, cluster_rows=TRUE, show_colnames=FALSE, show_rownames = FALSE, main="PSI.II 30 J2vs AvgTemp", breaks = seq(-2.5, 2, length.out=101))

#(K-mer clustering)
cg.res <- clusGap(mtx, kmeans, K.max=20, B=100, verbose=interactive())
nClus <- 9
#9 for PStime+LR251
#8 for temp+LR251
#6 for pstime+ps431
#3 for d2h+LR251
#4 for DayHour+LR251 pstime+kegg
clus <- kmeans(x = mtx, centers = nClus, iter.max = 50)
centers <- clus$centers-apply(clus$centers, 1, mean)
hclus <- hclust(as.dist(1-cor(t(centers), method = "pearson")), method = "average")
nOrder <- match(clus$cluster, hclus$order) 
Kmeans_matrix <- mtx[order(nOrder),]
Kmeans_cluster <- sort(nOrder)
# table(Kmeans_cluster)

Kmeans_matrix2 <- as.matrix(Kmeans_matrix) - apply(Kmeans_matrix, 1, mean) #centralization
#(polishing the view by removing the out-liners)
Kmeans_matrix2 <- as.matrix(Kmeans_matrix2)
cutoff <- median(unlist(Kmeans_matrix2)) + 3*sd (unlist(Kmeans_matrix2)) #remove median +3*sd
Kmeans_matrix2[Kmeans_matrix2 >cutoff] <- cutoff 
cutoff <- median(unlist(Kmeans_matrix2)) - 3*sd (unlist(Kmeans_matrix2)) #remove median +3*sd
Kmeans_matrix2[Kmeans_matrix2 < cutoff] <- cutoff

Kmeans_cluster2 <- factor(Kmeans_cluster, levels=unique(Kmeans_cluster))
levels(Kmeans_cluster2) <- c(1,3,2,4,5,7,6,8,9) #AvgTemp

rOrder <- order(as.vector(Kmeans_cluster2))
Kmeans_matrix2 <- Kmeans_matrix2[rOrder, ]
Kmeans_cluster2 <- factor(Kmeans_cluster2[rOrder], levels=c(1:nClus))

#(pheatmap.annt)
ngenes <- dim(Kmeans_matrix2)[1]
annt_row <- data.frame(row.names=rownames(Kmeans_matrix2), deg3v2=rep(1, times=ngenes), deg6v7=rep(1, times=ngenes), clus=Kmeans_cluster2)
annt_row[!rownames(annt_row) %in% cl3v2, "deg3v2"] <- NA
annt_row[!rownames(annt_row) %in% cl6v7, "deg6v7"] <- NA

label_row = rownames(Kmeans_matrix2)
label_row[!label_row %in% psI]  <- NA #40
label_row[!is.na(label_row)] <- "*"
label_row[is.na(label_row)] <- ""

key.palette <- rev(hcl.colors(100, "RdYlBu"))
pheatmap(Kmeans_matrix2, cluster_cols = FALSE, cluster_rows=FALSE, labels_row=label_row,
         show_colnames=FALSE, color = key.palette, breaks = seq(-4, 3, length.out=101),
         annotation_row = annt_row, main = paste0("LR251 on PStime smoothered; n=",ngenes))

#(single column scheme)
smth_dh1 <- t(gam_pst2)[,"DayHour"][1:100] 
smth_dh2 <- t(gam_pst2)[,"DayHour"][101:200] 
smth_tm1 <- t(gam_pst2)[,"AvgTemp"][1:100] + 1 #since +1 before smoothening to prevent the no-minus-input error.
smth_tm2 <- t(gam_pst2)[,"AvgTemp"][101:200]
mat <- matrix(c(smth_dh1, smth_dh2), ncol = 2)
colors <- hcl.colors(100, "Plasma")
breaks = unique(seq(10, 11, length.out=51), seq(11, 13, length.out=51))
pheatmap(mat, color = hcl.colors(100, "Blue-Yellow"),cluster_cols = FALSE,cluster_rows = FALSE, breaks=seq(9.5, 12, length.out=101), main="DayHour on pseudotime")

mat <- matrix(c(smth_tm1, smth_tm2), ncol = 2)
pheatmap(mat, color = hcl.colors(100, "Plasma"), cluster_cols = FALSE,cluster_rows = FALSE,  main="AvgTemp on pseudotime")


#######################
#### Figure 3F
#######################

meta2 <- meta[,c("fD2H","AvgTemp","DayHour","LineageM", "Accession")] #[1]   1940 5
ps.rpm <- t(all.rpm[rownames(all.rpm) %in% psI, rownames(meta2)]) #[1] 1940   28
ps.z <- scale(ps.rpm)

j.data <- cbind(meta2[meta2$Accession %in% "J247", ], ps.z[meta2$Accession %in% "J247",]) #[1] 460  33
o.data <- cbind(meta2[!meta2$Accession %in% "J247", ], ps.z[!meta2$Accession %in% "J247",]) #[1] 1480   33

plot(ps.z[rownames(j.data),1] ~ j.data$AvgTemp)

########## Slops by simple regression
library(tidyverse)
library(broom)
library(dplyr)
get_gene_slopes <- function(df, env_var) {
  # Get all gene column names
  gene_cols <- df %>% dplyr::select(starts_with("HORVU")) %>% names()
  
  # For each gene, fit linear model and extract slope
  slopes <- map_dbl(gene_cols, function(gene) {
    # Fit: gene_expression ~ environmental_variable
    # Expression = Intercept + Slope × AvgTemp + Error
    model <- lm(df[[gene]] ~ df[[env_var]])
    
    # Extract slope coefficient
    slope <- coef(model)[2]  # Second coefficient is the slope
    
    return(slope)
  })
  
  # Return tibble with gene names and slopes
  tibble(
    gene = gene_cols,
    slope = slopes
  )
}
slopes_j_temp <- get_gene_slopes(j.data, "AvgTemp")
slopes_o_temp <- get_gene_slopes(o.data, "AvgTemp")
slopes_j_hour <- get_gene_slopes(j.data, "DayHour")
slopes_o_hour <- get_gene_slopes(o.data, "DayHour")

paired_temp <- slopes_j_temp %>%
  dplyr::rename(slope_J247 = slope) %>%
  inner_join(slopes_o_temp %>% dplyr::rename(slope_Other = slope), by = "gene") %>%
  mutate(slope_diff = slope_Other - slope_J247)

paired_hour <- slopes_j_hour %>%
  dplyr::rename(slope_J247 = slope) %>%
  inner_join(slopes_o_hour %>% dplyr::rename(slope_Other = slope), by = "gene") %>%
  mutate(slope_diff = slope_Other - slope_J247)

# Wilcoxon test
wilcox_temp <- wilcox.test(paired_temp$slope_Other, paired_temp$slope_J247, paired = TRUE)
wilcox_hour <- wilcox.test(paired_hour$slope_Other, paired_hour$slope_J247, paired = TRUE)

##ggplot
combined_diffs <- bind_rows(
  paired_temp %>% mutate(Variable = "Temperature"),
  paired_hour %>% mutate(Variable = "Day hour")
) %>%
  mutate(Variable = factor(Variable, levels = c("Temperature", "Day hour")))

panel_c <- ggplot(combined_diffs, aes(x = Variable, y = slope_diff, fill = Variable)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
  geom_violin(alpha = 0.6, width = 0.7) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 2) +
  scale_fill_manual(values = c("Temperature" = "#3182bd", "Day hour" = "#e6550d")) +
  annotate("text", x = 1, y = max(combined_diffs$slope_diff) * 0.95, 
           label = "NS", size = 5, fontface = "italic") +
  annotate("text", x = 2, y = max(combined_diffs$slope_diff) * 0.95, 
           label = "*", size = 7) +
  labs(title = "Slope differences",
       x = NULL,
       y = expression(Delta * " Slope (Accession 2 - Accession 1)")) +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 11),
        axis.title = element_text(size = 10),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())

panel_c

######### Slops of individual genes
plot_comparison_temp <- bind_rows(
  df1 %>%
    dplyr::select(AvgTemp, starts_with("HORVU")) %>%
    pivot_longer(cols = starts_with("HORVU"), names_to = "gene", values_to = "expression") %>%
    mutate(Accession = "Accession 1"),
  df2 %>%
    dplyr::select(AvgTemp, starts_with("HORVU")) %>%
    pivot_longer(cols = starts_with("HORVU"), names_to = "gene", values_to = "expression") %>%
    mutate(Accession = "Accession 2")
) %>%
  mutate(gene_short = str_remove(gene, "HORVU")) %>%
  ggplot(aes(x = AvgTemp, y = expression, color = Accession)) +
  geom_point(alpha = 0.2, size = 0.15) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1.2) +
  scale_color_manual(values = c("Accession 1" = "#00bfc4", "Accession 2" = "grey50")) +
  facet_wrap(~gene_short, ncol = 7, scales = "free_y") +
  labs(title = "Temperature Response Comparison",
       x = "Air Temperature (°C)",
       y = "Expression (z-score)") +
  theme_bw() +
  theme(strip.text = element_text(size = 6, family = "mono"),
        axis.text = element_text(size = 6),
        legend.position = "bottom",
        panel.grid.minor = element_blank())


plot_comparison_dayl <- bind_rows(
  df1 %>%
    dplyr::select(DayHour, starts_with("HORVU")) %>%
    pivot_longer(cols = starts_with("HORVU"), names_to = "gene", values_to = "expression") %>%
    mutate(Accession = "Accession 1"),
  df2 %>%
    dplyr::select(DayHour, starts_with("HORVU")) %>%
    pivot_longer(cols = starts_with("HORVU"), names_to = "gene", values_to = "expression") %>%
    mutate(Accession = "Accession 2")
) %>%
  mutate(gene_short = str_remove(gene, "HORVU")) %>%
  ggplot(aes(x = DayHour, y = expression, color = Accession)) +
  geom_point(alpha = 0.2, size = 0.15) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1.2) +
  scale_color_manual(values = c("Accession 1" = "#00bfc4", "Accession 2" = "grey50")) +
  facet_wrap(~gene_short, ncol = 7, scales = "free_y") +
  labs(title = "Daylength Response Comparison",
       x = "Daylength (h)",
       y = "Expression (z-score)") +
  theme_bw() +
  theme(strip.text = element_text(size = 6, family = "mono"),
        axis.text = element_text(size = 6),
        legend.position = "bottom",
        panel.grid.minor = element_blank())


plot_comparison_temp
plot_comparison_dayl


####### And Slop diffence plot (the dumbbell plot )
paired_slopes <- paired_temp #paired_temp paired_hour

gene_order_temp <- paired_slopes %>%
  mutate(mean_slope = (slope_J247 + slope_Other) / 2) %>%
  arrange(desc(mean_slope)) %>%
  pull(gene) %>%
  str_remove("HORVU")

# Prepare data with ordered genes
all_genes_hour <- paired_slopes %>%
  mutate(gene_short = str_remove(gene, "HORVU"),
         gene_short = factor(gene_short, levels = rev(gene_order_temp)))
# Calculate x-axis range
x_range <- range(c(paired_slopes$slope_J247, paired_slopes$slope_Other))

# Panel D: Day hour - all 28 genes
panel_d_dumbbell <- ggplot(all_genes_hour, aes(y = gene_short)) +
  geom_segment(aes(x = slope_J247, xend = slope_Other, yend = gene_short),
               color = "gray50", linewidth = 0.6) +
  geom_point(aes(x = slope_J247), color = "#00bfc4", size = 4, alpha = 1) +
  geom_point(aes(x = slope_Other), color = "gray50", size = 4, alpha = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.5) +
  coord_cartesian(xlim = x_range) +
  labs(title = "Air Temp",
       x = "Slope (z-score change per hour)",
       y = NULL) +
  annotate("text", x = x_range[1] * 0.9, y = 28, 
           label = "J247", color = "#00bfc4", size = 3, fontface = "bold") +
  annotate("text", x = x_range[2] * 0.9, y = 28, 
           label = "Other", color = "gray50", size = 3, fontface = "bold") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 11),
        axis.title = element_text(size = 10),
        axis.text.y = element_text(family = "mono", size = 7),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(color = "gray90"))

panel_d_dumbbell

######### Expression under given pseudotime range
j.early.data <- j.data %>% filter(LineageM <= 15)
o.early.data <- o.data %>% filter(LineageM <= 15)
# expr_comparison <- bind_rows(
#   j.early.data %>% 
#     dplyr::select(starts_with("HORVU")) %>%
#     pivot_longer(everything(), names_to = "gene", values_to = "expression") %>%
#     mutate(Accession = "J247"),
#   o.early.data %>% 
#     dplyr::select(starts_with("HORVU")) %>%
#     pivot_longer(everything(), names_to = "gene", values_to = "expression") %>%
#     mutate(Accession = "Others")
# )
# 
# # Statistical test and effect size
# wilcox_early <- wilcox.test(expression ~ Accession, data = expr_comparison)
# n1_early <- sum(expr_comparison$Accession == "J247")
# n2_early <- sum(expr_comparison$Accession == "Others")
# r_rb <- 1 - (2 * wilcox_early$statistic) / (n1_early * n2_early) #effect factor (# Rank-biserial correlation)
# 
# expr_stats <- expr_comparison %>%
#   group_by(Accession) %>%
#   summarize(
#     mean_expr = mean(expression),
#     median_expr = median(expression)
#   )
# 
# ggplot(expr_comparison, aes(x = Accession, y = expression, fill = Accession)) +
#   geom_violin(alpha = 0.6, width = 0.7) +
#   geom_boxplot(width = 0.2, alpha = 0.8, outlier.alpha = 0.3) +
#   geom_text(data = expr_stats, 
#             aes(x = Accession, y = median_expr, 
#                 label = paste0("Median: ", round(median_expr, 3))),
#             vjust = -1.5, size = 3.5, fontface = "bold") +
#   scale_fill_manual(values = c("J247" = "#3182bd", "Others" = "#e6550d")) +
#   labs(title = "PS genes <10 PStime",
#        subtitle = paste0("Wilcoxon p < 0.001*** | Δ median = ",
#                         round(diff(expr_stats$median_expr), 3)," | ",
#                          "r = ", round(r_rb, 3)),
#        x = NULL,
#        y = "Expression (z-score)") +
#   theme_bw() +
#   theme(legend.position = "none",
#         plot.title = element_text(face = "bold", size = 11),
#         plot.subtitle = element_text(size = 9, color = "gray30"),
#         panel.grid.major.x = element_blank())
# 
# 
# p_exp <- ggplot(expr_comparison, aes(x = Accession, y = expression, fill = Accession)) +
#   geom_violin(alpha = 0.6, width = 0.7) +
#   geom_boxplot(width = 0.2, alpha = 0.8, outlier.alpha = 0.3) +
#   scale_fill_manual(values = c("J247" = "#3182bd", "Others" = "#e6550d")) +
#   labs(title = "Photosystem PS1&PS2 <10 pstime",
#        subtitle = paste0("Wilcoxon p = ", format.pval(wilcox_early$p.value, digits = 3)),
#        x = NULL,
#        y = "Expression (z-score)") +
#   theme_bw() +
#   theme(legend.position = "none",
#         plot.title = element_text(face = "bold", size = 11),
#         panel.grid.major.x = element_blank())
# 
# p_exp

#######################
#### Paired median comparison
gene_medians_acc1 <- df1 %>%
  filter(fD2H <= 0.3) %>%
  dplyr::select(starts_with("HORVU")) %>%
  summarise(across(everything(), median)) %>%
  pivot_longer(everything(), names_to = "gene", values_to = "median_expression")

gene_medians_acc2 <- df2 %>%
  filter(fD2H <= 0.3) %>%
  dplyr::select(starts_with("HORVU")) %>%
  summarise(across(everything(), median)) %>%
  pivot_longer(everything(), names_to = "gene", values_to = "median_expression")

#block2
# Create paired dataset
gene_medians_paired <- gene_medians_acc1 %>%
  dplyr::rename(median_acc1 = median_expression) %>%
  left_join(gene_medians_acc2 %>% dplyr::rename(median_acc2 = median_expression), by = "gene") %>%
  mutate(difference = median_acc2 - median_acc1)

# Paired Wilcoxon test
wilcox_paired <- wilcox.test(gene_medians_paired$median_acc2, 
                             gene_medians_paired$median_acc1, 
                             paired = TRUE)

# Prepare for plotting
gene_medians_long <- gene_medians_paired %>%
  pivot_longer(cols = c(median_acc1, median_acc2), 
               names_to = "Accession", 
               values_to = "median_expression") %>%
  mutate(Accession = ifelse(Accession == "median_acc1", "Accession 1", "Accession 2"))

# Calculate medians for the horizontal lines
medians_summary <- gene_medians_long %>%
  group_by(Accession) %>%
  summarize(median_value = median(median_expression))

# Plot with median lines
p_gene_medians <- ggplot(gene_medians_long, aes(x = Accession, y = median_expression, fill = Accession)) +
  geom_violin(alpha = 0.6, width = 0.7) +
  geom_line(aes(group = gene), alpha = 0.3, color = "gray50") +
  geom_point(alpha = 0.6, size = 2.5) +
  # Add horizontal median lines
  stat_summary(fun = median, geom = "errorbar", 
               aes(ymax = after_stat(y), ymin = after_stat(y)),
               width = 0.4, linewidth = 1, color = "black") +
  scale_fill_manual(values = c("Accession 1" = "#00bfc4", "Accession 2" = "grey50")) +
  coord_cartesian(ylim = c(-1, 1)) +
  labs(title = "Gene-level Median Expression; fD2H <= 0.8",
       subtitle = paste0("28 genes (paired) | p = ", format.pval(wilcox_paired$p.value, digits = 3)),
       x = NULL,
       y = "Median expression per gene (z-score)") +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major.x = element_blank())

p_gene_medians






sessionInfo()
# R version 4.4.2 (2024-10-31)
# Platform: aarch64-apple-darwin20
# Running under: macOS Sequoia 15.6.1
# 
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
# LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# time zone: Asia/Tokyo
# tzcode source: internal
# 
# attached base packages:
#   [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] cluster_2.1.8.1             MAST_1.32.0                 SingleCellExperiment_1.28.1 SummarizedExperiment_1.36.0 GenomicRanges_1.58.0        GenomeInfoDb_1.42.3        
# [7] MatrixGenerics_1.18.1       matrixStats_1.5.0           Seurat_5.3.0                SeuratObject_5.1.0          sp_2.2-0                    GO.db_3.20.0               
# [13] KEGGREST_1.46.0             GOstats_2.72.0              Category_2.72.0             Matrix_1.7-3                GSEABase_1.68.0             graph_1.84.1               
# [19] annotate_1.84.0             XML_3.99-0.18               AnnotationDbi_1.68.0        IRanges_2.40.1              S4Vectors_0.44.0            Biobase_2.66.0             
# [25] BiocGenerics_0.52.0         scales_1.4.0                pheatmap_1.0.13            
# 
# loaded via a namespace (and not attached):
#   [1] RColorBrewer_1.1-3      rstudioapi_0.17.1       jsonlite_2.0.0          magrittr_2.0.3          spatstat.utils_3.1-5    farver_2.1.2            zlibbioc_1.52.0        
# [8] vctrs_0.6.5             ROCR_1.0-11             memoise_2.0.1           spatstat.explore_3.5-2  RCurl_1.98-1.17         progress_1.2.3          S4Arrays_1.6.0         
# [15] htmltools_0.5.8.1       curl_7.0.0              SparseArray_1.6.2       sctransform_0.4.2       parallelly_1.45.1       KernSmooth_2.23-26      htmlwidgets_1.6.4      
# [22] ica_1.0-3               plyr_1.8.9              plotly_4.11.0           zoo_1.8-14              cachem_1.1.0            igraph_2.1.4            mime_0.13              
# [29] lifecycle_1.0.4         pkgconfig_2.0.3         R6_2.6.1                fastmap_1.2.0           GenomeInfoDbData_1.2.13 fitdistrplus_1.2-4      future_1.67.0          
# [36] shiny_1.11.1            digest_0.6.37           colorspace_2.1-1        patchwork_1.3.2         tensor_1.5.1            RSpectra_0.16-2         irlba_2.3.5.1          
# [43] RSQLite_2.4.2           progressr_0.15.1        spatstat.sparse_3.1-0   polyclip_1.10-7         abind_1.4-8             httr_1.4.7              compiler_4.4.2         
# [50] bit64_4.6.0-1           DBI_1.2.3               fastDummies_1.7.5       MASS_7.3-65             DelayedArray_0.32.0     tools_4.4.2             lmtest_0.9-40          
# [57] httpuv_1.6.16           future.apply_1.20.0     goftest_1.2-3           glue_1.8.0              nlme_3.1-168            promises_1.3.3          grid_4.4.2             
# [64] Rtsne_0.17              reshape2_1.4.4          generics_0.1.4          spatstat.data_3.1-8     gtable_0.3.6            tidyr_1.3.1             hms_1.1.3              
# [71] data.table_1.17.8       XVector_0.46.0          spatstat.geom_3.5-0     RcppAnnoy_0.0.22        ggrepel_0.9.6           RANN_2.6.2              pillar_1.11.0          
# [78] stringr_1.5.1           spam_2.11-1             RcppHNSW_0.6.0          genefilter_1.88.0       later_1.4.3             splines_4.4.2           dplyr_1.1.4            
# [85] lattice_0.22-7          deldir_2.0-4            survival_3.8-3          bit_4.6.0               tidyselect_1.2.1        RBGL_1.82.0             Biostrings_2.74.1      
# [92] miniUI_0.1.2            pbapply_1.7-4           gridExtra_2.3           scattermore_1.2         stringi_1.8.7           UCSC.utils_1.2.0        lazyeval_0.2.2         
# [99] codetools_0.2-20        tibble_3.3.0            Rgraphviz_2.50.0        cli_3.6.5               uwot_0.2.3              xtable_1.8-4            reticulate_1.43.0      
# [106] Rcpp_1.1.0              spatstat.random_3.4-1   globals_0.18.0          png_0.1-8               spatstat.univar_3.1-4   parallel_4.4.2          ggplot2_3.5.2          
# [113] blob_1.2.4              prettyunits_1.2.0       dotCall64_1.2           AnnotationForge_1.48.0  bitops_1.0-9            listenv_0.9.1           viridisLite_0.4.2      
# [120] ggridges_0.5.7          purrr_1.1.0             crayon_1.5.3            rlang_1.1.6             cowplot_1.2.0    
# 




