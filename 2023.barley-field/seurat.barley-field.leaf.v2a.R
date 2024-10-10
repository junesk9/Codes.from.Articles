#############################################################
########### Seurat output visualization & characterization
######################### Junesk9 2024.02.15

############################
######## Load libraries required
###########################

#remotes::install_version("SeuratObject", "4.1.4", repos = c("https://satijalab.r-universe.dev", getOption("repos")))
#remotes::install_version("Seurat", "4.4.0", repos = c("https://satijalab.r-universe.dev", getOption("repos")))
library(Seurat) #ver.4.4.0
library(cluster) #for clusGAP
library(gplots) #for heatmap.2
library(glmnet) #for glm

#For GESA analysis
library(GSEABase) 
library(GOstats) #GO/KEGG
library(KEGGREST) #KEGG-API

#Accessories
library(dunn.test)
library(beeswarm)
library(pheatmap)

#for SCODE
library(MASS)
library(stringr) #to modulate the row-names
library(reshape2) # for melt()


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
################################################# 
###################################### Load data
################################################
#1. Basal Expression data
so <- readRDS("seurat.barley-field.leaf.v3.rds")

meta <- so@meta.data #[1] 1940   35
all.rpm <- as.data.frame(so@assays$RNA@data) #[1] 25145  1940
wkMean <- function(rpm, meta){
  week <- meta$week
  mean.df <- data.frame(row.names=rownames(rpm))
  for(i in levels(week)){
    j247.sub <- subset(meta, meta$week == i & meta$Accession == "J247")
    j064.sub <- subset(meta, meta$week == i & meta$Accession == "J064")
    h602.sub <- subset(meta, meta$week == i & meta$Accession == "H602")
    j647.sub <- subset(meta, meta$week == i & meta$Accession == "J647")
    all.sub <- subset(meta, meta$week == i)
    j247.rpm <- rowMeans(rpm[, rownames(j247.sub)])
    j064.rpm <- rowMeans(rpm[, rownames(j064.sub)])
    h602.rpm <- rowMeans(rpm[, rownames(h602.sub)])
    j647.rpm <- rowMeans(rpm[, rownames(j647.sub)])
    all.rpm <- rowMeans(rpm[, rownames(all.sub)])
    
    sub.df <- data.frame(all.rpm, j247.rpm, j064.rpm, h602.rpm, j647.rpm)
    i <- sprintf("%02d", as.integer(i)) #convert "1" -> "01"
    title0 <- paste0("All_", i, "wk")
    title1 <- paste0("J247_", i, "wk")
    title2 <- paste0("J064_", i, "wk")
    title3 <- paste0("H602_", i, "wk")
    title4 <- paste0("J647_", i, "wk")
    title <- c(title0, title1, title2, title3, title4)
    colnames(sub.df) <- title
    
    mean.df <- cbind(mean.df, sub.df)
  }
  mean.df <- mean.df[, sort(colnames(mean.df))]
  return(mean.df)
}
all.wk <- wkMean(all.rpm, meta) #[1] 25145    110
vari.rpm <- all.rpm[VariableFeatures(so), ] #[1] 2963 1940
vari.wk <- all.wk[VariableFeatures(so), ] #[1] 2963   110

exp <- read.table("00_IPSR_KIBR_Sampledata/RPM_IPSR_KIBR_LB_17_19_BLD20201007.txt", header=T, row.names=1)
newCol <- sapply(colnames(exp), function(x){substr(x, 2, nchar(x))}) #remove the akward X from the colnames
colnames(exp) <- newCol
in.rpm <- exp[, newCol[grep("^5", newCol)]]
#[1] 39734    24
wkMean2 <- function(in.rpm){
  mean.df <- data.frame(row.names=rownames(in.rpm))
  rpt=3
  for (i in seq(1, length(colnames(in.rpm)), by=3)){
    col.st = i
    col.ed = i + rpt - 1
    wk.no = paste0("J247_wk", col.ed/rpt)
    sub.exp <- in.rpm[,c(col.st:col.ed)]
    mean = rowMeans(sub.exp)
    mean.df[,wk.no] <- mean
  }
  return(mean.df)
}
in.wk <- wkMean2(in.rpm)

flor <- read.table("00_IPSR_KIBR_Sampledata/FLOR_genes.txt", header=T, row.names=1, sep="\t")
flor <- rownames(flor) #[1]  149
flor.rpm <- all.rpm[rownames(all.rpm) %in% flor, ] #[1]  128 1940

tf_f <- "00_IPSR_KIBR_Sampledata/TranscriptionFactor_list_blastp_tophit_Gene_ID_rep_AtOsBd.txt"
tf_f <- read.table(tf_f, header=T, sep="\t")
tf <- unique(tf_f[,2]) #2567
tf.rpm <- all.rpm[rownames(all.rpm) %in% tf, ] #[1] 1468 1940
tff <- unique(c(tf, flor)) #[1] 2650
tff.rpm <- all.rpm[rownames(all.rpm) %in% tff, ] #[1] 1543 1940


gam_l0.all <- read.csv("03_TradeSeq-assoc/GAM-L0.All.100win.csv", row.names=1, header=T)
gam_l0.j2 <- read.csv("03_TradeSeq-assoc/GAM-L0.J2vsOther.100win.csv", row.names=1, header=T)
gam_lm.all <- read.csv("03_TradeSeq-assoc/GAM-LM.All.100win2.csv", row.names=1, header=T)
gam_lm.j2 <- read.csv("03_TradeSeq-assoc/GAM-LM.J2vsOther.100win2.csv", row.names=1, header=T)
gam_lm.temp <- read.csv("03_TradeSeq-assoc/GAM-LM.All.100win2.csv", row.names=1, header=T)["AvgTemp",] - 10 
gam_fD.all <- read.csv("03_TradeSeq-assoc/GAM-fD2H.All.100win.csv", row.names=1, header=T)
gam_fD.temp <- read.csv("03_TradeSeq-assoc/GAM-fD2H.All.100win2.csv", row.names=1, header=T)["AvgTemp",] - 10  #since +10 before to avoid (-)temp value.
gam_fD.j2 <- read.csv("03_TradeSeq-assoc/GAM-fD2H.J2vsOther.100win2.csv", row.names=1, header=T)
gam_CT.all <- read.csv("03_TradeSeq-assoc/GAM-CytoTRACE.All.100win.csv", row.names=1, header=T)
gam_CT.j2 <- read.csv("03_TradeSeq-assoc/GAM-CytoTRACE.J2vsOther.100win.csv", row.names=1, header=T)


#2. GO/KEGG daata
go_f ="./00_IPSR_KIBR_Sampledata/Barley.IBSCv2.ensembl49.GOSlim.txt"
go.df <- read.table(go_f, header=T, sep="\t")
go.df$evi = "IEA"
go.df <- go.df[,c(2,3,1)]
goFrame <- GOFrame(go.df, organism="barley")
goAllFrame <- GOAllFrame(goFrame)
gsc <- GeneSetCollection(goAllFrame, setType = GOCollection())
all <- unique(as.vector(go.df[,3]))  ## a vector of all referenced list
#[1] 27595

go_f ="./00_IPSR_KIBR_Sampledata/Barley.IBSCv2.ensembl49.GO.txt"
go.df <- read.table(go_f, header=T, sep="\t")
go.df$evi = "IEA"
go.df <- go.df[,c(2,3,1)]
goFrame <- GOFrame(go.df, organism="barley")
goAllFrame <- GOAllFrame(goFrame)
gsc <- GeneSetCollection(goAllFrame, setType = GOCollection())
all <- unique(as.vector(go.df[,3]))  ## a vector of all referenced list
#[1] 43051


kegg_f <- "./00_IPSR_KIBR_Sampledata/Barley.IBSCv2.ensembl49.KEGG.csv"
kegg.df <- read.csv(kegg_f, header=T, row.names=NULL) 
kegg.df$kegg.id <- sapply(kegg.df[,2], function(x){strsplit(x, "+", fixed=TRUE)[[1]][1]})
kegg.df2 <- unique(kegg.df[,c(3,1)]) #ほけん
#[1] 6062    2
head(kegg.df2)
#  KEGG.Pathway   Gene.stable.ID
#1        05235 HORVU0Hr1G004120
#2        05235 HORVU7Hr1G117420
#3        05235 HORVU0Hr1G001600
keggframe <- KEGGFrame(kegg.df2, organism="barely")
gsck <- GeneSetCollection(keggframe, setType=KEGGCollection())
kegg.all <- unique(as.vector(kegg.df2[,2]))
#[1] 2555
kegg.annt <- keggList(database="pathway", organism="taes")
kegg.annt <- data.frame(KEGGid=names(kegg.annt), desc=as.vector(kegg.annt))
rownames(kegg.annt) <- kegg.annt[,1]
kegg.annt[,2] <- sapply(kegg.annt[,2], function(x){strsplit(x, " - ")[[1]][1]})
#                  KEGGid                       desc
#taes04146 taes04146                 Peroxisome
#taes04136 taes04136                  Autophagy
#taes04148 taes04148              Efferocytosis
#taes04814 taes04814             Motor proteins

#3.  DEG dataset
cl2 <- read.csv("02_cluster-DEG/cl2-pos.markers.MAST.csv", header=T, row.names=1)
cl3 <- read.csv("02_cluster-DEG/cl3-pos.markers.MAST.csv", header=T, row.names=1)
cl2n3 <- read.csv("02_cluster-DEG/cl2n3-pos.markers.MAST.csv", header=T, row.names=1)
cl2v3 <- read.csv("02_cluster-DEG/cl2vs3.markers.MAST.csv", header=T, row.names=1)
cl4 <- read.csv("02_cluster-DEG/cl4-pos.markers.MAST.csv", header=T, row.names=1)
cl5 <- read.csv("02_cluster-DEG/cl5-pos.markers.MAST.csv", header=T, row.names=1)
cl4n5 <- read.csv("02_cluster-DEG/cl4n5-pos.markers.MAST.csv", header=T, row.names=1)
cl4v5 <- read.csv("02_cluster-DEG/cl4vs5.markers.MAST.csv", header=T, row.names=1)
j247 <- read.csv("02_cluster-DEG/j247.markers.MAST.csv", header=T, row.names=1)

cl2 <- rownames(subset(cl2, cl2$avg_log2FC >= .5 & cl2$p_val_adj < 0.05)) #5908
cl3 <- rownames(subset(cl3, cl3$avg_log2FC >= .5 & cl3$p_val_adj < 0.05)) #744
cl2n3 <- rownames(subset(cl2n3, cl2n3$avg_log2FC >= .5 & cl2n3$p_val_adj < 0.05)) #3713
cl2v3.up <- rownames(subset(cl2v3, cl2v3$avg_log2FC >= 0.5 & cl2v3$p_val_adj < 0.05)) #[1] 4722
cl2v3.dn <- rownames(subset(cl2v3, cl2v3$avg_log2FC <= -0.5 & cl2v3$p_val_adj < 0.05)) #[1] 156
cl4 <- rownames(subset(cl4, cl4$avg_log2FC >= 1 & cl4$p_val_adj < 0.05)) #37 #338 #965
cl5 <- rownames(subset(cl5, cl5$avg_log2FC >= 1 & cl5$p_val_adj < 0.05)) #124 #535  #896
cl4n5 <- rownames(subset(cl4n5, cl4n5$avg_log2FC >= .5 & cl4n5$p_val_adj < 0.05)) #52 #435 #913
cl4v5.up <- rownames(subset(cl4v5, cl4v5$avg_log2FC >= 1 & cl4v5$p_val_adj < 0.05)) #43
cl4v5.dn <- rownames(subset(cl4v5, cl4v5$avg_log2FC <= -1 & cl4v5$p_val_adj < 0.05)) #372
j247.up <- rownames(subset(j247, j247$avg_log2FC >= 1 & j247$p_val_adj < 0.05)) #409
j247.dn <- rownames(subset(j247, j247$avg_log2FC <= -1 & j247$p_val_adj < 0.05)) #325
j247.deg <- c(j247.up, j247.dn)

#3.  tradeseq-DEG dataset
time.l0 <- read.csv("03_TradeSeq-assoc/Tradeseq-gam.L0-assoc.csv", header=T, row.names=1)
time.lm <- read.csv("03_TradeSeq-assoc/Tradeseq-gam.LM-assoc.csv", header=T, row.names=1)
time.fd <- read.csv("03_TradeSeq-assoc/Tradeseq-gam.fD2H-assoc.csv", header=T, row.names=1)
time.ct <- read.csv("03_TradeSeq-assoc/Tradeseq-gam.CytoTRACE-assoc.csv", header=T, row.names=1)

################################################# 
############################# Data subsetting
################################################
deg.l0 <- rownames(subset(time.l0, time.l0$pvalue < 0.01))#1253  /2963
deg.lm <- rownames(subset(time.lm, time.lm$pvalue < 0.01))#1248  /2963
deg.fd <- rownames(subset(time.fd, time.fd$pvalue < 0.01))#1305  /2963
deg.ct <- rownames(subset(time.ct, time.ct$pvalue < 0.01))#1009  /2963

j247.l0 <- rownames(subset(time.l0, time.l0$pvalue < 0.01 & time.l0$j247.pval < 0.01)) #158   
j247.lm <- rownames(subset(time.lm, time.lm$pvalue < 0.01 & time.lm$J247.pval < 0.01)) #220   
j247.fd <- rownames(subset(time.fd, time.fd$pvalue < 0.01 & time.fd$j247.pval < 0.01)) #166   
j247.ct <- rownames(subset(time.ct, time.ct$pvalue < 0.01 & time.ct$j247.pval < 0.01)) #110   
unique(c(j247.l0, j247.lm, j247.fd, j247.ct)) #591

j247.fd[j247.fd %in% flor] #[1] "HORVU5Hr1G095630" Vrn-H1
j247.l0[j247.l0 %in% flor] #[1] "HORVU3Hr1G095240" HvOS2
j247.lm[j247.lm %in% flor] #[1] "HORVU2Hr1G063800" "HORVU0Hr1G003020" HvFUL2, HvFUL3
j247.ct[j247.ct %in% flor] #[1] "HORVU3Hr1G114970" "HORVU2Hr1G117380" HvLUX, ABF4
j247.deg[j247.deg %in% flor] #[1] "HORVU3Hr1G010240" OsRAV9

flor.l0 <- rownames(subset(time.l0, time.l0$pvalue < 0.01 & rownames(time.l0) %in% flor)) #11
flor.lm <- rownames(subset(time.lm, time.lm$pvalue < 0.01 & rownames(time.lm) %in% flor)) #14
flor.fd <- rownames(subset(time.fd, time.fd$pvalue < 0.01 & rownames(time.fd) %in% flor)) #14
flor.ct <- rownames(subset(time.ct, time.ct$pvalue < 0.01 & rownames(time.ct) %in% flor)) #13

tff.l0 <- rownames(subset(time.l0, time.l0$pvalue < 0.01 & rownames(time.l0) %in% tff)) #90
tff.lm <- rownames(subset(time.lm, time.lm$pvalue < 0.01 & rownames(time.lm) %in% tff)) #101
tff.fd <- rownames(subset(time.fd, time.fd$pvalue < 0.01 & rownames(time.fd) %in% tff)) #102
tff.ct <- rownames(subset(time.ct, time.ct$pvalue < 0.01 & rownames(time.ct) %in% tff)) #84


#################################################
######## Sequence of temp.day length by along with fD2h
#################################################
x <- read.csv("03_TradeSeq-assoc/Tradeseq-gam.fD2H-meta2.csv", header=T)
# X lineage      time    gene     yhat
# 1 1       1 0.1184858 DayHour 10.86627
# 2 2       1 0.1283923 DayHour 10.87095
# 3 3       1 0.1382989 DayHour 10.87563
# 4 4       1 0.1482055 DayHour 10.88032
# 5 5       1 0.1581121 DayHour 10.88502
x <- as.vector(x[c(1:100), "time"])
# [1] 0.1184858 0.1283923 0.1382989 ... 1.079423 1.089330 1.099237
x.df <- data.frame(row.names=x, fD2H=x, AvgTemp=NA, DayLen=NA)
d <- meta[order(meta$fD2H), ]
for(rn in rownames(d)){
  z = 0
  for (fd in rownames(x.df)){
    md <- meta[rn,]$fD2H
    at <- meta[rn,]$AvgTemp
    dl <- meta[rn,]$DayHour
    if (as.double(md) < as.double(fd) & z < 1) {
      x.df[fd, ]$AvgTemp <- paste(x.df[fd, ]$AvgTemp, at, sep=";")
      x.df[fd, ]$DayLen <- paste(x.df[fd, ]$DayLen, dl, sep=";")
      z = z + 10
    }
  }
}
funMean <- function(x){mean(as.double(strsplit(x,";")[[1]]), na.rm=TRUE)}
funMed <- function(x){median(as.double(strsplit(x,";")[[1]]), na.rm=TRUE)}
x.df$AT.mean <- NA
x.df$AT.med <- NA
x.df$DL.mean <- NA
x.df$DL.med <- NA
for (i in rownames(x.df)){
  at <- x.df[i,]$AvgTemp
  dl <- x.df[i,]$DayLen
  if (is.na(at) == FALSE){
    x.df[i,]$AT.mean <- funMean(at)
    x.df[i,]$AT.med <- funMed(at)
  }
  if (is.na(dl) == FALSE){
    x.df[i,]$DL.mean <- funMean(dl)
    x.df[i,]$DL.med <- funMed(dl)
  }
}
write.csv(x.df[,!c(2,3)], "03_TradeSeq-assoc/Tradeseq-gam.fD2H-meta3.csv")

y = x.df$AT.med
fit = lm(y ~ c(1:100) + I(c(1:100)^2))
plot(y, pch=19, cex=0.5, xlab="100-win", main="GAM-fD2H 100-win")
lines(predict(fit), col=2, lwd=2)
plot(predict(fit))


################################################# 
############################# Kmean clustering
################################################
doKmean1 <- function(preSm, gids, title){
  x4 <- (preSm - apply(preSm, 1, mean) ) / apply(preSm, 1, sd)  # Z-standardization; simply x4=t(scale(t(preSm)))
  x <- x4
  x <- x[rownames(x) %in% gids, ] 
  
  title = paste0(title, "z-scaled, first view")
  cg.res <- clusGap(x, kmeans, K.max=20, B=100, verbose=interactive()) ##taking minutes
  require(pheatmap)
  p <- pheatmap(x, cluster_cols = FALSE, show_rownames = FALSE, show_colnames = FALSE
           , main = title) #(first veiw)
  
  print(cg.res)
  par(mfrow = c(1,1))
  plot(cg.res)
  return(x)
}
doKmean2 <- function(x, nClus, title){
  clus <- kmeans(x = x, centers = nClus, iter.max = 50)
  centers <- clus$centers-apply(clus$centers, 1, mean)
  hclus <- hclust(as.dist(1-cor(t(centers), method = "pearson")), method = "average")
  
  title = paste0(title, " ,z-scaled, hclust-mean")
  p <- pheatmap(x, cluster_cols = FALSE, cluster_rows = hclus, show_rownames = FALSE, show_colnames = FALSE
           , main = title) #(hclust-mean veiw)
  #exchange the "kmean cluster" to "hcluster" order
  nOrder <- match(clus$cluster, hclus$order) 
  Kmeans_matrix <- x[order(nOrder),]
  Kmeans_cluster <- sort(nOrder)
  ngenes <- as.character(table(Kmeans_cluster))
  print(ngenes)
  
  out_ls <- list(Kmeans_matrix, Kmeans_cluster)
  return(out_ls)
}
doHM4 <- function(mtx, cls, ngenes){
  #mtx <- t(scale(t(mtx))) # z-standardization
  #breaks = seq(-3, 3, by=0.1)
  #ngenes <- as.character(table(cls))
  hm <- heatmap.2(x = mtx,
                  Rowv =F,
                  Colv=F,
                  dendrogram ="none",
                  density.info="none",
                  trace="none",
                  scale="none",
                  keysize=.3,
                  key= F, 
                  labRow = F,
                  labCol = F,
                  RowSideColors = rainbow(length(levels(cls)))[cls],
                  margins = c(8, 16), #for labels
                  srtCol=45,
                  col=key.palette,
  )
  legend.text <- paste(toupper(letters)[as.numeric(levels(cls))],
                       " (N=", ngenes,")",
                       sep="")
  par(lend = 1) # square line ends for the color legend
  legend(x = "bottomright", 
         legend = legend.text, # category labels
         col = rainbow(length(levels(cls))),  # color key
         lty= 1, # line style
         lwd = 10 # line width
  )
  
  legend_image <- as.raster(rev(hm$colorTable$color), ncol=1)
  rasterImage(legend_image, 0.8, 0.8, 0.85, 1)
  
  return(hm)
}
doLP3 <- function(mtx, cls){
  mtx <- as.matrix(mtx) - apply(mtx, 1, mean) #centralization
  ymax = 0
  ymin = 0
  Ncls = length(levels(cls))
  par(mfrow = c(3,3))
  
  for(i in unique(cls)){
    submtx <- mtx[cls %in% i, ]
    submean <- apply(submtx, 2, mean)
    submax <- max(submean)
    submin <- min(submean)
    if (submax > ymax){ymax = submax}
    if (submin < ymin){ymin = submin}
  }
  for(i in unique(cls)){
    i = as.numeric(i)
    submtx <- mtx[cls %in% i, ]
    submean <- apply(submtx, 2, mean)
    title = paste("Cluster ", toupper(letters)[as.numeric(i)], " (N=", dim(submtx)[1],")", sep="")
    if (length(submean) > 100){
      submean1 <- submean[1:100] ## for j247 vs other
      submean2 <- submean[101:200]
      submean <- as.matrix(data.frame(submean1, submean2))
      matplot(submean, type="l", lwd=3, lty=c(1,2), col=rainbow(Ncls)[i],
              ylab="Mean z-scaled exp", ylim=c(ymin, ymax), main=title)
    } else {
      matplot(submean, type="l", lwd=4, lty=1, col=rainbow(Ncls)[i],
              ylab="Mean z-scaled exp", ylim=c(ymin, ymax), main=title)
    }
    abline(h=0, lty=2)
  }
  
  par(mfrow = c(1,1))
}
doKmean3 <- function(mtx_cls, nRnd){
  Kmeans_matrix <- mtx_cls[[1]]
  Kmeans_cluster <- as.factor(mtx_cls[[2]])
  ngenes <- as.character(table(Kmeans_cluster))
  Kmeans_matrix2 <- as.matrix(Kmeans_matrix) - apply(Kmeans_matrix, 1, mean) #centralization
  out_ls <- list(Kmeans_matrix2, Kmeans_cluster)
  
  n_randomSample = nRnd
  if(dim(Kmeans_matrix2)[1] > n_randomSample){
    ix <- sort(sample(1:length(Kmeans_cluster), n_randomSample))
    Kmeans_cluster2 <- Kmeans_cluster[ix]
    Kmeans_matrix2 <- Kmeans_matrix2[ix,]
  } else {
    Kmeans_cluster2 <- Kmeans_cluster
    Kmeans_matrix2 <- Kmeans_matrix2
  }
  
  #(polishing the view by removing the out-liners)
  Kmeans_matrix2 <- as.matrix(Kmeans_matrix2)
  cutoff <- median(unlist(Kmeans_matrix2)) + 3*sd (unlist(Kmeans_matrix2)) #remove median +3*sd
  Kmeans_matrix2[Kmeans_matrix2 >cutoff] <- cutoff 
  cutoff <- median(unlist(Kmeans_matrix2)) - 3*sd (unlist(Kmeans_matrix2)) #remove median +3*sd
  Kmeans_matrix2[Kmeans_matrix2 < cutoff] <- cutoff
  
  hm <- doHM4(Kmeans_matrix2, Kmeans_cluster2, ngenes)
  par(mfrow = c(1,1))
  doLP3(Kmeans_matrix, Kmeans_cluster)
  
  par(mfrow = c(1,1))
  return(out_ls)
}
doKmean4 <- function(mtx_cls2, gid){
  Kmeans_matrix2 <- mtx_cls2[[1]]
  Kmeans_cluster2 <- mtx_cls2[[2]]
  
  ix <- sort(match(gid, rownames(Kmeans_matrix2)))
  Kmeans_matrix3 <- Kmeans_matrix2[ix, ]
  Kmeans_cluster3 <- Kmeans_cluster2[ix]
  
  ngenes <- as.character(table(Kmeans_cluster3))
  hm <- doHM4(Kmeans_matrix3, Kmeans_cluster3, ngenes)
  
  return(rownames(Kmeans_matrix3))
} #subset & HM again
doLPacc1 <- function(deg){
  deg <- deg[deg %in% rownames(all.wk)]
  rpm.sub <- all.wk[deg, ]
  for(id in deg){
    pdf_f <- paste0(id,".meanRPM-byWeek.pdf")
    new.df <- data.frame(row.names = levels(meta$week))
    for(wk in levels(meta$week)){
      wk = as.integer(wk)
      wk.st = (wk-1)*4 + 1
      wk.ed = wk.st + 3
      new.line <- rpm.sub[id,c(wk.st:wk.ed)]
      colnames(new.line) <- levels(meta$Accession)
      new.df <- rbind(new.df, new.line)
      #new.df[wk,] <- rpm.sub[id,c(wk.st:wk.ed)]
    }
    pdf(file=pdf_f, height=6, width=5)
    matplot(new.df, type="b", lwd=3, lty=1, pch=19,col=ggcolor(4),
            ylab="Mean RPM", xlab="week",main=id)
    dev.off()
  }
  
  deg <- deg[deg %in% rownames(gam_lm.all)]
  for(id in deg){
    parseDF <- function(gam.all, gam.j2){
      All <- as.vector(t(gam.all[id, ]))
      J247 <- as.vector(t(gam.j2[id, c(1:100)]))
      other <- as.vector(t(gam.j2[id, c(101:200)]))
      tmp.df <- data.frame(All, J247, other)
      return(tmp.df)
    }
    l0.df <- parseDF(gam_l0.all, gam_l0.j2)
    lm.df <- parseDF(gam_lm.all, gam_lm.j2)
    fd.df <- parseDF(gam_fD.all, gam_fD.j2)
    ct.df <- parseDF(gam_CT.all, gam_CT.j2)
    
    pdf_f2 <- paste0(id,".meanRPM-byTimeLine.pdf")
    pdf(file=pdf_f2, height=8, width=12)
    par(mfrow = c(2,2))
    matplot(fd.df, type="l", lwd=c(3,5,5), lty=1, pch=19,col=c(1, ggcolor(4)[1], "grey"),
            ylab="GAM-fitted RPM", main="fD2H-smoothed")
    matplot(lm.df, type="l", lwd=c(3,5,5), lty=1, pch=19,col=c(1, ggcolor(4)[1], "grey"),
            ylab="GAM-fitted RPM", main="LM-smoothed")
    matplot(ct.df, type="l", lwd=c(3,5,5), lty=1, pch=19,col=c(1, ggcolor(4)[1], "grey"),
            ylab="GAM-fitted RPM", main="CytoTRACE-smoothed")
    matplot(l0.df, type="l", lwd=c(3,5,5), lty=1, pch=19,col=c(1, ggcolor(4)[1], "grey"),
            ylab="GAM-fitted RPM", main="L0-smoothed")
    dev.off()
  }
  dev.off()
}
doLPacc2 <- function(deg){
  deg <- deg[deg %in% rownames(all.wk)]
  rpm.sub <- all.wk[deg, ]
  for(id in deg){
    pdf_f <- paste0(id,".meanRPM-byWeek.pdf")
    new.df <- data.frame(row.names = levels(meta$week))
    new.df$J247 <- as.vector(t(rpm.sub[id, c(67:88)]))
    new.df$J064 <- as.vector(t(rpm.sub[id, c(45:66)]))
    new.df$J647 <- as.vector(t(rpm.sub[id, c(89:110)]))
    new.df$H602 <- as.vector(t(rpm.sub[id, c(23:44)]))
    new.df$all <- as.vector(t(rpm.sub[id, c(1:22)]))
    
    pdf(file=pdf_f, height=5.5, width=5)
    matplot(new.df[,c(3,2,1,4)], type="b", ylim=c(0, max(new.df[,c(1:4)], na.rm=TRUE)),
            lwd=3, lty=1, pch=19,col=ggcolor(4), ylab="Mean RPM", xlab="week",main=id)
    dev.off()
  }
    
  #deg <- deg[deg %in% rownames(gam_lm.all)]
  #deg <- deg[deg %in% rownames(gam_lm.j2)]
  gam.deg <- gam_lm.j2[rownames(gam_lm.j2) %in% deg, ]
  gam.deg <- scale(gam.deg)
  gam.deg <- as.matrix(gam.deg) - apply(gam.deg, 1, mean)
  
  for(id in deg){
    # parseDF <- function(gam.all, gam.j2){
    #   All <- as.vector(t(gam.all[id, ]))
    #   J247 <- as.vector(t(gam.j2[id, c(1:100)]))
    #   other <- as.vector(t(gam.j2[id, c(101:200)]))
    #   tmp.df <- data.frame(other, J247)
    #   tmp.df <- scale(tmp.df)
    #   return(tmp.df)
    # }
    # lm.df <- parseDF(gam_lm.all, gam_lm.j2)
    lm.df <- data.frame(gam.deg[id, c(101:200)], gam.deg[id, c(1:100)])
    
    pdf_f2 <- paste0(id,".meanRPM-byTimeLine.pdf")
    pdf(file=pdf_f2, height=5.5, width=5)
    matplot(lm.df[,c(1:2)], type="l", ylim=c(-2,2), lwd=c(3,3), lty=1, pch=19,
            col=c("grey20", ggcolor(4)[3]), ylab="GAM-fitted RPM", main="LM-smoothed")
    abline(h=0, lty=2, lwd=0.5, col="grey")
    dev.off()
  }
}
doLP4 <- function(goi){
  goi <- goi[goi %in% rownames(all.wk)]
  goi.wk <- all.wk[goi, ]
  goi.rpm <- as.data.frame(t(all.rpm[goi, ]))
  goi.rpm$week <- meta[rownames(goi.rpm),]$week
  
  for(gid in goi){
    pdf_f <- paste0(gid,".meanRPM-byWeek.pdf")
    pdf_f2 <- paste0(gid,".meanRPM-byWeek2.pdf")
    new.df <- data.frame(row.names = levels(meta$week))
    new.df$J247 <- as.vector(t(goi.wk[gid, c(67:88)]))
    new.df$J064 <- as.vector(t(goi.wk[gid, c(45:66)]))
    new.df$J647 <- as.vector(t(goi.wk[gid, c(89:110)]))
    new.df$H602 <- as.vector(t(goi.wk[gid, c(23:44)]))
    
    sub.rpm <- goi.rpm[, c("week",gid)]
    sub.j247 <- sub.rpm[rownames(subset(meta, meta$Accession %in% "J247")), ]
    sub.j064 <- sub.rpm[rownames(subset(meta, meta$Accession %in% "J064")), ] 
    sub.j647 <- sub.rpm[rownames(subset(meta, meta$Accession %in% "J647")), ]
    sub.h602 <- sub.rpm[rownames(subset(meta, meta$Accession %in% "H602")), ]
    
    pdf(file=pdf_f, height=5.5, width=5)
    matplot(new.df[,c(1:4)], type="b", lwd=6, lty=1, pch=19, col=ggcolor(4)[c(3,2,1,4)], xlab="weeks", ylab="RPM", main=gid)
    dev.off()
    
    
    pdf(file=pdf_f2, height=5.5, width=5)
    matplot(y=sub.rpm[, gid], x=jitter(as.numeric(sub.rpm$week), factor=0.7), type="n", pch=19, cex=.7,
            xlab="weeks", ylab="RPM", main=gid, ylim=c(0, max(sub.rpm[, gid])))
    matplot(y=sub.j064[, gid], x=jitter(as.numeric(sub.j064$week), factor=0.7), type="p", pch=20, cex=.7, 
            col=paste0(ggcolor(4)[1], "33"), add=TRUE) #30% transparency "4D" #20% transparency "33"
    matplot(y=sub.j647[, gid], x=jitter(as.numeric(sub.j647$week), factor=0.7), type="p", pch=20, cex=.7, 
            col=paste0(ggcolor(4)[2], "33"), add=TRUE) #30% transparency "4D"
    matplot(y=sub.j247[, gid], x=jitter(as.numeric(sub.j247$week), factor=0.7), type="p", pch=20, cex=.7, 
            col=paste0(ggcolor(4)[3], "33"), add=TRUE) #30% transparency "4D"
    matplot(y=sub.h602[, gid], x=jitter(as.numeric(sub.h602$week), factor=0.7), type="p", pch=20, cex=.7, 
            col=paste0(ggcolor(4)[4], "33"), add=TRUE) #30% transparency "4D"
    matplot(new.df[,c(1:4)], type="b", lwd=3, lty=1, pch=19, col=ggcolor(4)[c(3,2,1,4)], add=TRUE)
    dev.off()

  }
}
runGO3 <- function(gid_ls, prefix){
  out_ls = list()
  for (f in c("MF","CC","BP")){
    p <- GSEAGOHyperGParams(
      name = "Paramaters",
      geneSetCollection = gsc,
      geneIds = gid_ls,
      universeGeneIds = all,
      ontology = f,
      pvalueCutoff = 1,
      conditional = FALSE,
      testDirection = "over"
    )
    res <- hyperGTest(p)
    res_gid <- geneIdsByCategory(res)
    res2 <- summary(res)
    res2$GeneId <- sapply(res2[,1], function(x)paste(unlist(res_gid[x]), collapse="/"))
    out_ls = c(out_ls, list(res2))
  }
  write.csv(out_ls[[1]], paste0(prefix,".GO-MP.csv"), row.names=F)
  write.csv(out_ls[[2]], paste0(prefix,".GO-CC.csv"), row.names=F)
  write.csv(out_ls[[3]], paste0(prefix,".GO-BP.csv"), row.names=F)
  
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
  
  write.csv(kegg.res2, paste0(prefix,".KEGG.csv"), row.names=F)
  
  return(out_ls)
} 



# #(7. cl2n3)
# key.palette <- rev(hcl.colors(100, "Spectral"))
# cl2n3 <- cl2n3 #[1] 3713
# prefix = "gam_lM.cl2n3"
# x <- doKmean1(gam_lm.all, cl2n3, prefix)
# nClus=3
# y <- doKmean2(x, nClus, prefix)
# z <- doKmean3(y, 200)
# hm.cls <- as.data.frame(z[[1]])
# hm.cls$clus <- z[[2]]
# write.csv(hm.cls, "03_TradeSeq-assoc/GAM_lM.cl2n3.csv")  #[1] 3713  101
# #hm.cls <- read.csv("03_TradeSeq-assoc/GAM_lM.cl2n3.csv", header=T, row.names=1)
# 
# cl2a <- rownames(hm.cls[hm.cls$clus %in% 1, ]) #1230
# cl2b <- rownames(hm.cls[hm.cls$clus %in% 2, ]) #1000
# cl2c <- rownames(hm.cls[hm.cls$clus %in% 3, ]) #1483
# go.cl2a <- runGO3(cl2a, "03_TradeSeq-assoc/GAM_lM.cl2n3-a2")
# go.cl2b <- runGO3(cl2b, "03_TradeSeq-assoc/GAM_lM.cl2n3-b2")
# go.cl2c <- runGO3(cl2c, "03_TradeSeq-assoc/GAM_lM.cl2n3-c2")
# 
# prefix = "gam_lM.cl2n3-ac"
# x <- doKmean1(gam_lm.j2, c(cl2a, cl2c), prefix)
# nClus=6
# y <- doKmean2(x, nClus, prefix)
# z <- doKmean3(y, 200)
# hm.cls <- as.data.frame(z[[1]])
# hm.cls$clus <- z[[2]]
# write.csv(hm.cls, "03_TradeSeq-assoc/GAM_lM.cl2n3-c.csv") #[1] 2713  201
# cl2ac1 <- rownames(hm.cls[hm.cls$clus %in% 1, ]) #402
# cl2ac2 <- rownames(hm.cls[hm.cls$clus %in% 2, ]) #392
# cl2ac12 <- c(cl2ac1, cl2ac2) 
# 
# go.cl2ac1 <- runGO3(cl2ac1, "03_TradeSeq-assoc/GAM_lM.cl2n3-ac1ii")
# go.cl2ac2 <- runGO3(cl2ac2, "03_TradeSeq-assoc/GAM_lM.cl2n3-ac2ii")
# go.cl2ac12 <- runGO3(cl2ac12, "03_TradeSeq-assoc/GAM_lM.cl2n3-ac12ii")

#(6. flor)
prefix = "gam_lM.j2-flor"
x <- doKmean1(gam_lm.j2, flor, "gam_lM.j2-flor")
nClus=5
y <- doKmean2(x, nClus, prefix)
z <- doKmean3(y, 200)
hm.cls <- as.data.frame(z[[1]])
hm.cls$clus <- z[[2]]
write.csv(hm.cls, "03_TradeSeq-assoc/GAM-lM.flor156.csv")  #[1] 175 101


#(4. cls4n5)
#cl45 <- unique(c(cl4, cl5, cl4n5)) #145  #845
cl45 <- cl4n5 #435 
x <- doKmean1(gam_fD.all, cl45, "gam_fD.all-cl45")
nClus=3
y <- doKmean2(x, nClus, "gam_fD.all-cl45")
dev.off()
z <- doKmean3(y, 200)
dev.off()
hm.cls <- as.data.frame(z[[1]])
hm.cls$clus <- z[[2]]
write.csv(hm.cls, "03_TradeSeq-assoc/HM-cl.all-fD2H-GAM.clus4n5.csv")  #[1] 175 101

go <- runGO3(cl45, "./03_TradeSeq-assoc/GAM-fD.cl4n5.all") #List of (GOMF GOCC GOBP KEGG)
go <- runGO3(cl45, "./03_TradeSeq-assoc/GAM-fD.cl4n5.all2") #full-depth GO
gobp <- go[[3]]
kegg <- go[[4]]
kegg.id <- strsplit(kegg[1,length(kegg)], "/", fixed=TRUE)[[1]]
## photosynthesis was a top rated GOBP 89/431 & carbon fixation was a top rated KEGG 10/87
kegg.ps <- kegg.df[kegg.df$kegg.id %in% "00710",] #[1] 87
kegg.ps[kegg.ps$Gene.stable.ID %in% kegg.id,]  #9/28
# Gene.stable.ID   KEGG.Pathway kegg.id 
# 3047 HORVU0Hr1G030760 00710+4.1.1.39   00710
# 3055 HORVU7Hr1G065780 00710+4.1.1.39   00710
# 3057 HORVU2Hr1G084630 00710+4.1.1.39   00710
# 3058 HORVU7Hr1G091620 00710+4.1.1.39   00710
# 3061 HORVU2Hr1G010670 00710+4.1.1.39   00710
# 3062 HORVU2Hr1G010630 00710+4.1.1.39   00710
# 3064 HORVU2Hr1G084620 00710+4.1.1.39   00710
# 3069 HORVU6Hr1G049260 00710+4.1.1.39   00710
# 3070 HORVU1Hr1G035720 00710+4.1.1.39   00710
# 3105 HORVU0Hr1G004830 00710+1.2.1.12   00710
go.ps <- unique(go.df[go.df$GOSlim.GOA.Accession.s. %in% "GO:0015979",3]) #[1] 431
go.goi <- strsplit(gobp[1,length(gobp)], "/", fixed=TRUE)[[1]]
table(hm.cls[go.goi,"clus"])
# 1  2  3 
# 61  9 19 

x <- doKmean1(gam_fD.all, go.ps, "gam_fD.all-Photosynthesis")
nClus=2
y <- doKmean2(x, nClus, "gam_fD.all-Photosynthesis")
dev.off()
z <- doKmean3(y, 200)
dev.off()
hm.cls <- as.data.frame(z[[1]])
hm.cls$clus <- z[[2]]
write.csv(hm.cls, "03_TradeSeq-assoc/gam_fD.all-Photosynthesis.csv")  
table(hm.cls[go.goi,"clus"])
# 1   2   3 
# 0 62 27 

x <- doKmean1(gam_fD.all, kegg.ps[,1], "gam_fD.all-KEGG-00710")
nClus=3
y <- doKmean2(x, nClus, "gam_fD.all-KEGG-00710")
dev.off()
z <- doKmean3(y, 200)
hm.cls <- as.data.frame(z[[1]])
hm.cls$clus <- z[[2]]
table(hm.cls[kegg.id,"clus"])
# 1  2  3 
# 0  0 10 
write.csv(hm.cls, "03_TradeSeq-assoc/gam_fD.all-KEGG-00710.csv")  

kegg.ps1 <- subset(kegg.df, kegg.df$KEGG.Pathway == "00710+4.1.1.39")$Gene.stable.ID
kegg.ps2 <- subset(kegg.df, kegg.df$KEGG.Pathway == "00710+1.2.1.12")$Gene.stable.ID
kegg.ps12 <- unique(c(kegg.ps1, kegg.ps2))
x <- doKmean1(gam_fD.all, kegg.ps12, "gam_fD.all-Photosynthesis") 
nClus=2
y <- doKmean2(x, nClus, "gam_fD.all-Photosynthesis")
z <- doKmean3(y, 200)
median <- rowMedians(as.matrix(all.rpm[rownames(x), ]))
hm.cls <- as.data.frame(z[[1]])
hm.cls$clus <- z[[2]]
hm.cls$median <- median[rownames(hm.cls)]
write.csv(hm.cls, "03_TradeSeq-assoc/gam_fD.all-KEGG-00710.csv")  


# #(5. cls3v2)
#(analyze correlation to AvgTemp)
goi3 <- cl2v3.dn #156   #[cl2v3.dn %in% cl3] #[1] 41
go3 <- runGO3(goi3, "03_TradeSeq-assoc/Cl2v3-dn")
goi2 <- cl2v3.up #4722   #[cl2v3.up %in% cl2] #[1] 3373
go2 <- runGO3(goi2, "03_TradeSeq-assoc/Cl2v3-up")

goi <- goi3
corr.df <- data.frame(row.names=goi)
corr.df$pccALL <- NA
corr.df$pccJ2 <- NA
corr.df$pccOTH <- NA
for(g in goi){
  Ax <- t(all.rpm[g, rownames(meta)])
  Ay <- as.vector(meta$AvgTemp)
  Ap <- cor.test(x=Ax, y=Ay, method="pearson")[4]
  
  Jm <- subset(meta, meta$Accession %in% "J247")
  Jx <- t(all.rpm[g, rownames(Jm)])
  Jy <- as.vector(Jm$AvgTemp)
  Jp <- cor.test(x=Jx, y=Jy, method="pearson")[4]
  
  Om <- subset(meta, meta$Accession %notin% "J247")
  Ox <- t(all.rpm[g, rownames(Om)])
  Oy <- as.vector(Om$AvgTemp)
  Op <- cor.test(x=Ox, y=Oy, method="pearson")[4]
  
  corr.df[g, "pccALL"] <- Ap
  corr.df[g, "pccJ2"] <- Jp
  corr.df[g, "pccOTH"] <- Op
}

goi <- go.df[go.df$GO.term.accession %in% c("GO:0009768","GO:0009769"), ]$Gene.stable.ID #29
goi <- unique(goi[goi %in% rownames(all.rpm)]) #28
coef.df <- data.frame(row.names=goi)
coef.df$All <- NA
coef.df$J247 <- NA
coef.df$Other <- NA
coef.df$J647 <- NA
coef.df$J064 <- NA
for(g in goi){
  dp.sub10A <- subset(meta, meta$AvgTemp <= 10) # subset(meta, meta$DayHour <= 11)# meta#     #[1] 1269   37
  dp.sub10J <- dp.sub10A[dp.sub10A$Accession %in% "J247", ] #[1] 316
  dp.sub10O <- dp.sub10A[dp.sub10A$Accession %notin% "J247", ] #[1] 953
  dp.sub10d <- dp.sub10A[dp.sub10A$Accession %in% "H602", ] 
  dp.sub10b <- dp.sub10A[dp.sub10A$Accession %in% "J064", ]
  dp.sub10c <- dp.sub10A[dp.sub10A$Accession %in% "J647", ]
  
  t.sub10A <- dp.sub10A$AvgTemp
  t.sub10J <- dp.sub10J$AvgTemp
  t.sub10O <- dp.sub10O$AvgTemp
  t.sub10d <- dp.sub10d$AvgTemp
  t.sub10b <- dp.sub10b$AvgTemp
  t.sub10c <- dp.sub10c$AvgTemp
  h.sub10A <- dp.sub10A$DayHour
  h.sub10J <- dp.sub10J$DayHour
  h.sub10O <- dp.sub10O$DayHour
  h.sub10d <- dp.sub10d$DayHour
  h.sub10b <- dp.sub10b$DayHour
  h.sub10c <- dp.sub10c$DayHour
  
  x.sub10A <- scale(t(all.rpm[g, rownames(dp.sub10A)]))
  x.sub10J <- x.sub10A[rownames(dp.sub10J), ]
  x.sub10O <- x.sub10A[rownames(dp.sub10O), ]
  x.sub10d <- x.sub10A[rownames(dp.sub10d), ]
  x.sub10b <- x.sub10A[rownames(dp.sub10b), ]
  x.sub10c <- x.sub10A[rownames(dp.sub10c), ]
  
  coef.df[g, "All"] <- lm(x.sub10A ~ t.sub10A)$coefficients[2]
  coef.df[g, "J247"] <- lm(x.sub10J ~ t.sub10J)$coefficients[2]
  coef.df[g, "Other"] <- lm(x.sub10O ~ t.sub10O)$coefficients[2]
  coef.df[g, "J064"] <- lm(x.sub10b ~ t.sub10b)$coefficients[2]
  coef.df[g, "J647"] <- lm(x.sub10c ~ t.sub10c)$coefficients[2]
  coef.df[g, "H602"] <- lm(x.sub10d ~ t.sub10d)$coefficients[2]
  
  # plot(x=t.sub10A, y=c(x.sub10J, x.sub10O), type="n", main=g, xlab="DayHour", ylab="z-RPM")
  # points(x=t.sub10O, y=x.sub10O, pch=19, cex=0.5, col="grey50")
  # points(x=t.sub10J, y=x.sub10J, pch=19, cex=0.5, col=ggcolor(4)[3])
  # abline(lm(x.sub10J ~ t.sub10J), col=ggcolor(4)[3], lwd=3)
  # abline(lm(x.sub10O ~ t.sub10O), col=1, lwd=3)
}
boxplot(coef.df, notch=T, pch=19, main="sub11h-DayHour", ylab="lm-coefcient")
boxplot(coef.df[,c(2,5,4,6)], col=ggcolor(4)[c(3,2,1,4)], notch=T, pch=19
        ,main="sub10C-AirTemp", ylab="lm-coefcient",ylim=c(-1,1) )
        #,main="sub11h-DayLen", ylab="lm-coefcient",ylim=c(-1,1) )
        #,main="Whole-DayLen", ylab="lm-coefcient", ylim=c(-0.5,0.7))
dunn.test(coef.df[,c(2,5,4,6)], method="bonferroni")
boxplot(coef.df[,c(2,3)], col=ggcolor(4)[c(3,2,1,4)], notch=T, pch=19
        ,main="sub11h-DayHour", ylab="lm-coefcient", ylim=c(-1, 1))

library(exactRankTests)
wilcox.exact(x=coef.df$J247,y=coef.df$Other,paired=T)

###multi-panel graph
par("mar"=c(1,1,1,1))
par(mfrow = c(5, 6))
for(g in goi){
  dp.sub10A <- subset(meta, meta$AvgTemp <= 10)# subset(meta, meta$DayHour <= 11)#     #[1] 1269   37
  dp.sub10J <- dp.sub10A[dp.sub10A$Accession %in% "J247", ] #[1] 316
  dp.sub10O <- dp.sub10A[dp.sub10A$Accession %notin% "J247", ] #[1] 953
  dp.sub10d <- dp.sub10A[dp.sub10A$Accession %in% "H602", ] 
  dp.sub10b <- dp.sub10A[dp.sub10A$Accession %in% "J064", ]
  dp.sub10c <- dp.sub10A[dp.sub10A$Accession %in% "J647", ]
  
  t.sub10A <- dp.sub10A$AvgTemp
  t.sub10J <- dp.sub10J$AvgTemp
  t.sub10O <- dp.sub10O$AvgTemp
  t.sub10d <- dp.sub10d$AvgTemp
  t.sub10b <- dp.sub10b$AvgTemp
  t.sub10c <- dp.sub10c$AvgTemp
  h.sub10A <- dp.sub10A$DayHour
  h.sub10J <- dp.sub10J$DayHour
  h.sub10O <- dp.sub10O$DayHour
  h.sub10d <- dp.sub10d$DayHour
  h.sub10b <- dp.sub10b$DayHour
  h.sub10c <- dp.sub10c$DayHour
  
  x.sub10A <- scale(t(all.rpm[g, rownames(dp.sub10A)]))
  x.sub10J <- x.sub10A[rownames(dp.sub10J), ]
  x.sub10O <- x.sub10A[rownames(dp.sub10O), ]
  x.sub10d <- x.sub10A[rownames(dp.sub10d), ]
  x.sub10b <- x.sub10A[rownames(dp.sub10b), ]
  x.sub10c <- x.sub10A[rownames(dp.sub10c), ]
  
  coef.df[g, "All"] <- lm(x.sub10A ~ t.sub10A)$coefficients[2]
  coef.df[g, "J247"] <- lm(x.sub10J ~ t.sub10J)$coefficients[2]
  coef.df[g, "Other"] <- lm(x.sub10O ~ t.sub10O)$coefficients[2]
  coef.df[g, "J064"] <- lm(x.sub10b ~ t.sub10b)$coefficients[2]
  coef.df[g, "J647"] <- lm(x.sub10c ~ t.sub10c)$coefficients[2]
  coef.df[g, "H602"] <- lm(x.sub10d ~ t.sub10d)$coefficients[2]
  
  plot(x=h.sub10A, y=c(x.sub10J, x.sub10O),  main=g, xlab="", ylab = "",  yaxt="n", xaxt="n", pch=19, 
       cex=.5, col=c(rep("#00BFC433", length(x.sub10J)), rep("#7f7f7f33", length(x.sub10O))))
  abline(lm(x.sub10J ~ h.sub10J), col="#00BFC4", lwd=6)
  abline(lm(x.sub10O ~ h.sub10O), col="grey50", lwd=6)
  
  
  # plot(x=t.sub10A, y=c(x.sub10J, x.sub10O), type="n", main=g, xlab="DayHour", ylab="z-RPM")
  # points(x=t.sub10O, y=x.sub10O, pch=19, cex=0.5, col="grey50")
  # points(x=t.sub10J, y=x.sub10J, pch=19, cex=0.5, col=ggcolor(4)[3])
  # abline(lm(x.sub10J ~ t.sub10J), col=ggcolor(4)[3], lwd=3)
  # abline(lm(x.sub10O ~ t.sub10O), col=1, lwd=3)
}



# go.cl2 <- runGO3(cl2, "./03_TradeSeq-assoc/GAM-fD.cl2.all") #[1] 5908
# go.cl3 <- runGO3(cl3, "./03_TradeSeq-assoc/GAM-fD.cl3.all") #[1] 744
# 
# c2.kegg <- go.cl2[[4]] 
# c2.g <- strsplit(c2.kegg[,length(c2.kegg)][1],"/", fixed=TRUE)[[1]]
# c2.gk <- subset(kegg.df, kegg.df$Gene.stable.ID %in% c2.g & kegg.df$kegg.id == "00680")
# table(c2.gk$KEGG.Pathway)
# # 00680+1.1.1.37 00680+1.1.1.95 00680+1.17.1.9  00680+2.1.2.1 00680+2.6.1.52 00680+2.7.1.11 00680+2.7.1.29 00680+3.1.2.12  00680+3.1.3.3 
# # 4              1              1              3              2              4              1              1              2 
# # 00680+4.1.1.31 00680+4.1.2.13 00680+4.2.1.11 00680+4.4.1.19 00680+4.4.1.22 00680+5.4.2.11 00680+5.4.2.12  00680+6.2.1.1 
# # 2              4              2              1              1              2              1              2 
# 
# c3.k <- go.cl3[[4]] 
# c3.gk <- strsplit(c3.kegg[,length(c3.kegg)][1],"/", fixed=TRUE)[[1]]
# c3.gk <- subset(kegg.df, kegg.df$Gene.stable.ID %in% c3.g & kegg.df$kegg.id == "00680")
# table(c3.gk$KEGG.Pathway)
# # 00680+1.1.1.37 00680+2.7.1.11  00680+3.1.3.3 00680+4.1.1.31 00680+4.1.2.13 00680+4.4.1.22 00680+5.4.2.11 00680+5.4.2.12 
# #1              2              1              2              3              1              2              1 
# 
# c3.mf <- go.cl3[[1]]
# c3.m1 <- strsplit(c3.mf[,length(c3.mf)][1],"/", fixed=TRUE)[[1]] #312
# c3.m2 <- strsplit(c3.mf[,length(c3.mf)][2],"/", fixed=TRUE)[[1]] #131 (of #312 the upper)
# c3.m3 <- strsplit(c3.mf[,length(c3.mf)][3],"/", fixed=TRUE)[[1]] #58
# go.m3 <- runGO3(c3.m3, "./03_TradeSeq-assoc/GAM-lm.cl3-transporter.all") 
# go.m2 <- runGO3(c3.m2, "./03_TradeSeq-assoc/GAM-lm.cl3-hydrolysis.all") 
# go.m1 <- runGO3(c3.m1, "./03_TradeSeq-assoc/GAM-lm.cl3-catalysis.all") 
# c3.123 <- unique(c(c3.m1, c3.m2, c3.m3)) #[1] 364
# go.m123 <- runGO3(c3.123, "./03_TradeSeq-assoc/GAM-lm.cl3-cattrans.all") 
# prefix = "gam_lm.cluster3-123"
# x3 <- doKmean1(gam_lm.all, c3.123, prefix)
# nClus=2
# y3 <- doKmean2(x3, nClus, prefix)
# z3 <- doKmean3(y3, 200)
# hm.cls3 <- as.data.frame(z3[[1]])
# hm.cls3$clus <- z3[[2]]
# table(hm.cls3$clus)
# #1   2 
# #152 212 
# write.csv(hm.cls3, "03_TradeSeq-assoc/gam_lm.cluster3.2.csv")  
# c3.k9 <-strsplit(c3.k[,length(c3.k)][9],"/", fixed=TRUE)[[1]]
# 
# 
# prefix = "gam_lm.cluster3"
# x <- doKmean1(gam_lm.all, cl3, prefix)
# nClus=4
# y <- doKmean2(x, nClus, prefix)
# z <- doKmean3(y, 200)
# hm.cls <- as.data.frame(z[[1]])
# hm.cls$clus <- z[[2]]
# table(hm.cls$clus)
# #   1   2   3   4 
# #166  99 262 217 
# write.csv(hm.cls, "03_TradeSeq-assoc/gam_lm.cluster3.csv")  
# 
# prefix = "gam_lm.cluster3-j247"
# x <- doKmean1(gam_lm.j2, cl3, prefix)
# nClus=4
# y <- doKmean2(x, nClus, prefix)
# z <- doKmean3(y, 200)
# hm.cls2 <- as.data.frame(z[[1]])
# hm.cls2$clus <- z[[2]]
# table(hm.cls2$clus)
# #  1   2   3   4 
# #133 213 218 180 
# write.csv(hm.cls2, "03_TradeSeq-assoc/gam_lm.cluster3-j247.2.csv")  
# 
# 
# cl3.1 <- rownames(hm.cls2[hm.cls2$clus %in% "1", ]) #133
# cl3.2 <- rownames(hm.cls2[hm.cls2$clus %in% "2", ]) #213
# cl3.3 <- rownames(hm.cls2[hm.cls2$clus %in% "3", ]) #218
# cl3.12 <- unique(c(cl3.1, cl3.2))
# cl3.123 <- unique(c(cl3.12, cl3.3))
# go.cl3.1 <- runGO3(cl3.1, "./03_TradeSeq-assoc/GAM-lm.cl3-1")
# go.cl3.2 <- runGO3(cl3.2, "./03_TradeSeq-assoc/GAM-lm.cl3-2")
# go.cl3.12 <- runGO3(cl3.12, "./03_TradeSeq-assoc/GAM-lm.cl3-12")
# go.cl3.3 <- runGO3(cl3.3, "./03_TradeSeq-assoc/GAM-lm.cl3-3")
# go.cl3.123 <- runGO3(cl3.123, "./03_TradeSeq-assoc/GAM-lm.cl3-123")
# 
# cl3.n1 <- rownames(hm.cls[hm.cls$clus %in% c("2","3","4"), ]) #578
# go.cl3.n1 <- runGO3(cl3.n1, "./03_TradeSeq-assoc/GAM-lm.cl3-n1")
# 
# 
# 
# 
# go1 <- runGO3(cl2v3.dn, "./03_TradeSeq-assoc/GAM-fD.cl2v3.dn") #[1] 156
# #go2 <- runGO3(intersect(cl2v3.dn, cl2n3), "./03_TradeSeq-assoc/GAM-fD.cl2v3.dn2") #[1] 28
# go3 <- runGO3(intersect(cl2v3.dn, cl3), "./03_TradeSeq-assoc/GAM-fD.cl2v3.dn3") #[1] 41
# go4 <- runGO3(intersect(cl2v3.dn, j247.dn), "./03_TradeSeq-assoc/GAM-fD.cl2v3.dn4") #[1] 73
# go5 <- runGO3(intersect(cl3, j247.dn), "./03_TradeSeq-assoc/GAM-fD.cl2v3.dn5") #[1] 174
# go1.bp <- go1[[3]]
# go1.st <- strsplit(go1.bp[,length(go1.bp)][4], "/", fixed=TRUE)[[1]]
# go1.li <- strsplit(go1.bp[,length(go1.bp)][2], "/", fixed=TRUE)[[1]]
# 
# go3.bp <- go3[[3]]
# go3.g <- strsplit(go3.bp[,length(go3.bp)][1], "/", fixed=TRUE)[[1]]



# x <- doKmean1(gam_fD.j2[,c(1:100)], cl45, "gam_fD.J247-cl45")
# nClus=2
# y <- doKmean2(x, nClus, "gam_fD.J247-cl45")
# dev.off()
# z <- doKmean3(y, 200)
# dev.off()
# hm.cls <- as.data.frame(z[[1]])
# hm.cls$clus <- z[[2]]
# write.csv(hm.cls, "03_TradeSeq-assoc/HM-cl.J247-fD2H-GAM.clus4n5.csv")  #[1] 143 101
# 
# group.A <- rownames(hm.cls[hm.cls$clus %in% "1",])
# group.B <- rownames(hm.cls[hm.cls$clus %in% "2",])
# runGO3(group.A, "./03_TradeSeq-assoc/GAM0-fD.cl45-HM.groupA")
# runGO3(group.B, "./03_TradeSeq-assoc/GAM0-fD.cl45-HM.groupB")
# runGO3(c(group.B, group.A), "./03_TradeSeq-assoc/GAM0-fD.cl45-HM.groupAB")
# runKEGG1(group.A, "./03_TradeSeq-assoc/GAM0-fD.cl45-HM.groupA")
# runKEGG1(group.B, "./03_TradeSeq-assoc/GAM0-fD.cl45-HM.groupB")
# runKEGG1(c(group.B, group.A), "./03_TradeSeq-assoc/GAM0-fD.cl45-HM.groupAB")
# 
# wk.field <- all.wk[group.B, ] #[1]  58 110
# wk.room <- in.wk[group.B, ]  #[1] 58  8
# wk.join <- cbind(wk.field[,c(89:110)], wk.room)
# matplot(scale(t(wk.field[,c(89:110)])), type="l", ylim=c(-3,3), main="field J247 58-gene cl4n5")
# matplot(scale(t(wk.room)), type="l", ylim=c(-3,3), main="room J247 58-gene cl4n5")
# wk.field.m <- rowMeans(scale(t(wk.field[,c(89:110)])))
# wk.room.m <- rowMeans(scale(t(wk.room)))
# 
# B.rpm <- t(all.rpm[group.B,rownames(meta)]) #[1] 1940   58
# B.scale <- t(scale(t(B.rpm)))
# plot(B.rpm[,"HORVU6Hr1G055770"] ~ meta$AvgTemp)
# plot(B.rpm[,"HORVU6Hr1G055770"] ~ meta$DayHour)
# 
# cor.test(B.rpm[,"HORVU6Hr1G055770"], meta$AvgTemp, method="pearson")[c(3,4)] #0.3575134 
# cor.test(B.scale[,"HORVU6Hr1G055770"], meta$AvgTemp, method="pearson")[c(3,4)] #0.3604606 
# 
# est <- cor.test(exp.so[,i], exp.vf[,j], method="pearson")[c(3,4)]
# estC <- est$estimate
# estP <- est$p.value
# cor_v <- c(cor_v, estC)
# pval_v <- c(pval_v, estP)


#(2. j247-lineage M, degs)
hm.deg <- unique(c(j247.deg, j247.lm)) #809
x <- doKmean1(gam_lm.j2, hm.deg, "fD2H-j247.deg2")
nClus = 4
y <- doKmean2(x, nClus, "LM-j247.deg2")
z <- doKmean3(y, 200)
hm.cls <- as.data.frame(z[[1]])
hm.cls$clus <- z[[2]]
#hm.cls <- read.csv("03_TradeSeq-assoc/GAM-lm.J247-DEG1.csv", header=T, row.names=1)
hm.cls[rownames(hm.cls) %in% c("HORVU5Hr1G095630", "HORVU2Hr1G063800", "HORVU0Hr1G003020", "HORVU3Hr1G010240"),c(1, length(hm.cls))] 
write.csv(hm.cls, "03_TradeSeq-assoc/GAM-lm.J247-DEG2.csv") #[1] 809 201
j2deg.1 <- rownames(hm.cls[hm.cls$clus %in% "1",])
j2deg.2 <- rownames(hm.cls[hm.cls$clus %in% "2",])
j2deg.3 <- rownames(hm.cls[hm.cls$clus %in% "3",])
j2deg.4 <- rownames(hm.cls[hm.cls$clus %in% "4",])
j2.1.GO <- runGO3(j2deg.1, "03_TradeSeq-assoc/GAM-lm.J247-DEG2-D2")
j2.2.GO <- runGO3(j2deg.2, "03_TradeSeq-assoc/GAM-lm.J247-DEG2-B2")
j2.3.GO <- runGO3(j2deg.3, "03_TradeSeq-assoc/GAM-lm.J247-DEG2-A2")
j2.4.GO <- runGO3(j2deg.4, "03_TradeSeq-assoc/GAM-lm.J247-DEG2-C2")
j2.34.GO <- runGO3(c(j2deg.3,j2deg.4), "03_TradeSeq-assoc/GAM-lm.J247-DEG2-AC2")
j2.12.GO <- runGO3(c(j2deg.2,j2deg.1), "03_TradeSeq-assoc/GAM-lm.J247-DEG2-BD2")




# write.csv(hm.cls, "03_TradeSeq-assoc/HM-cl.j247-LM-GAM.clus5.csv") #[1] 944 201
# doLPacc1(rownames(hm.cls))
# cl1.go <- runGO3(rownames(hm.cls[hm.cls$clus %in% "1",]), "./03_TradeSeq-assoc/HM-cl.j247-LM-GAM2.clus1")
# cl2.go <- runGO3(rownames(hm.cls[hm.cls$clus %in% "2",]), "./03_TradeSeq-assoc/HM-cl.j247-LM-GAM2.clus2")
# cl3.go <- runGO3(rownames(hm.cls[hm.cls$clus %in% "3",]), "./03_TradeSeq-assoc/HM-cl.j247-LM-GAM2.clus3")
# cl4.go <- runGO3(rownames(hm.cls[hm.cls$clus %in% "4",]), "./03_TradeSeq-assoc/HM-cl.j247-LM-GAM2.clus4")
# cl0.go <- runGO3(rownames(hm.cls), "./03_TradeSeq-assoc/HM-cl.j247-LM-GAM2.clus0")
# cl12.go <- runGO3(rownames(hm.cls[hm.cls$clus %in% c("1","2"),]), "./03_TradeSeq-assoc/HM-cl.j247-LM-GAM2.clus12")
# cl34.go <- runGO3(rownames(hm.cls[hm.cls$clus %in% c("3","4"),]), "./03_TradeSeq-assoc/HM-cl.j247-LM-GAM2.clus34")
# runGO3(rownames(hm.cls[hm.cls$clus %in% "5",]), "./03_TradeSeq-assoc/HM-cl.j247-LM-GAM.clus5")
# 
# runGO3(cl2v3.up, "./02_cluster-DEG/cl2vs3.markers.MAST-up")
# runGO3(cl2v3.dn, "./02_cluster-DEG/cl2vs3.markers.MAST-dn")
# runGO3(cl2, "./02_cluster-DEG/cl2-pos.markers.MAST")
# runGO3(cl3, "./02_cluster-DEG/cl3-pos.markers.MAST")




#################################################################
############# UMAP by accession
#################################################################
so.j2 <- subset(so, subset = Accession == "J247")
#25145 features across 460 samples within 1 assay 
so.j0 <- subset(so, subset = Accession == "J064")
#25145 features across 468 samples within 1 assay 
so.j6 <- subset(so, subset = Accession == "J647")
#25145 features across 508 samples within 1 assay 
so.h6 <- subset(so, subset = Accession == "H602")
#25145 features across 504 samples within 1 assay 

#(PCA -> select PC -> UMAP)
so.sub <- so.j2
so.sub@reductions$pca <- NULL #nullify the data
so.sub@reductions$umap <- NULL #nullify the data
VariableFeatures(so.sub) <- NULL #nullify the data

so.sub <- FindVariableFeatures(so.sub, selection.method = "mvp")
#VariableFeaturePlot(so.sub)
VariableFeatures(so.sub) <- unique(c(VariableFeatures(so.sub), flor))
af <- rownames(so.sub@assays$RNA@data)
so.sub <- RunPCA(so.sub, vervose=FALSE, approx=FALSE, npcs=50, features=VariableFeatures(so.sub))
RunEN <- function(exp, fD2H, th){
  alpha <- seq(0.01, 0.99, 0.01)
  mse.df <- NULL
  for (i in 1:length(alpha)) {
    m <- cv.glmnet(x = exp, y = fD2H, family = "gaussian", alpha = alpha[i])
    mse.df <- rbind(mse.df, data.frame(alpha = alpha[i], mse = min(m$cvm)))
  }
  best.alpha <- mse.df$alpha[mse.df$mse == min(mse.df$mse)]
  print(best.alpha)
  
  en.model.cv <- cv.glmnet(x = exp, y = fD2H, family = "gaussian", alpha = best.alpha)
  plot(en.model.cv, xvar="lambda", label=TRUE, main="EN model fitting")
  en.model <- glmnet(x = exp, y = fD2H, family = "gaussian", lambda = en.model.cv$lambda.1se)
  plot(abs(en.model$beta), main="betas of PC-dim, EN-model", pch=19)
  abline(h=th, lty=2, col=2)
  
  PCs <- which(abs(en.model$beta)>th)
  return(PCs)
  
}
PCs <- RunEN(so.sub@reductions$pca@cell.embeddings, so.sub$fD2H, 0.007)
#[1] 0.47 | [1]  1  3  6 10 11 15 17 21 25 #J247 h=0.005
#[1] 0.4 |  [1]  1  3  5  8 11 13 20 21 #J064  h=0.005
#[1] 0.2  | [1]  1  3 11 12 16 26 35 40 #J647  h=0.005
#[1] 0.48  | [1]  1  3  7 10 14 18 23 27 43 #H602 h=0.006

so.sub <- RunUMAP(so.sub, reduction="pca", dims=PCs, umap.method = "umap-learn", metric = "correlation", n.components = 2)
FeaturePlot(so.sub, reduction="umap", features="fD2H")
DimPlot(so.sub, reduction="umap", group.by=c("month","seurat_clusters", "Case"), ncol=2)

umap.sub <- as.data.frame(so.sub@reductions$umap@cell.embeddings)
umap.sub <- umap.sub[rownames(meta),]
so@meta.data$H602_UMAP1 <- umap.sub[,1]
so@meta.data$H602_UMAP2 <- umap.sub[,2]

saveRDS(so.sub, "seurat.barley-field.leaf.v2.J247.rds")
saveRDS(so.sub, "seurat.barley-field.leaf.v2.J064.rds")
saveRDS(so.sub, "seurat.barley-field.leaf.v2.J647.rds")
saveRDS(so.sub, "seurat.barley-field.leaf.v2.H602.rds")



#################################################################
############# PCC to the Indoor GrowthChamber data
#################################################################
exp <- read.table("00_IPSR_KIBR_Sampledata/RPM_IPSR_KIBR_LB_17_19_BLD20201007.txt", header=T, row.names=1)
#                   X10001_IH6_171215 X10002_IH6_171215 X10003_IH6_171215 X10004_IH6_171215 X10005_IH6_171215
# HORVU0Hr1G000020          2.871378          8.346515          5.677581          7.880693          4.925499
# HORVU0Hr1G000030          3.445654         10.361191         10.949620         13.065359          7.258630
# HORVU0Hr1G000050          3.445654          7.483082          8.921913          2.696027          5.184736
# HORVU0Hr1G000080          4.019930          4.317163          4.055415          7.258533          2.592368
# HORVU0Hr1G000100          0.000000          0.000000          0.000000          0.000000          0.000000

newCol <- sapply(colnames(exp), function(x){substr(x, 2, nchar(x))}) #remove the akward X from the colnames
colnames(exp) <- newCol

#(case of variable features)
exp.vf <- exp[VariableFeatures(so), newCol[grep("^5", newCol)]]
exp.so <- so@assays$RNA@data[VariableFeatures(so), ]
#[1] 2963   24
#[1] 2963 1940
#(case of flowering genes)
exp.so <- so@assays$RNA@data[rownames(so@assays$RNA@data) %in% flor, ]
exp.vf <- exp[rownames(exp.so),  newCol[grep("^5", newCol)]]
#[1] 128 1940
#[1] 128  24
exp.so <- so@assays$RNA@data
exp.vf <- exp[rownames(exp.so), newCol[grep("^5", newCol)]]
#[1] 25145  1940
#[1] 25145    24


#(try. 1)
maxC_val = c()
minP_val = c()
maxC_wk = c()
minP_wk = c()
wk_no <- rep(c(1:8), each=3)
for (i in 1:dim(exp.so)[2]){
  cor_v = c()
  pval_v = c()
  for (j in 1:dim(exp.vf)[2]){
    est <- cor.test(exp.so[,i], exp.vf[,j], method="pearson")[c(3,4)]
    estC <- est$estimate
    estP <- est$p.value
    cor_v <- c(cor_v, estC)
    pval_v <- c(pval_v, estP)
  }
  maxC_idx <- order(-cor_v)[1]
  minP_idx <- order(pval_v)[1]
  maxC_val = c(maxC_val, cor_v[maxC_idx])
  minP_val = c(minP_val, pval_v[minP_idx])
  maxC_wk = c(maxC_wk, wk_no[maxC_idx])
  minP_wk = c(minP_wk, wk_no[minP_idx])
}
table(maxC_wk == minP_wk) #check the matchness of two predictions

so@meta.data$RoomJ247wk.VF <-  factor(maxC_wk, levels=c(1:8))
so@meta.data$RoomJ247wk.VF.PCC <-  as.vector(maxC_val)
so@meta.data$RoomJ247wk.FLOR <-  factor(maxC_wk, levels=c(1:8))
so@meta.data$RoomJ247wk.FLOR.PCC <-  as.vector(maxC_val)
so@meta.data$RoomJ247wk.ALL <-  factor(maxC_wk, levels=c(1:8))

metaj <- subset(meta, meta$Accession == "J247")
require(beeswarm)
beeswarm(metaj$fD2H ~metaj$RoomJ247wk.FLOR, corral="wrap", pch=19, cex=0.5, pwcol=ggcolor(7)[metaj$seurat_clusters], main="FLOR129; J247")
boxplot(metaj$fD2H ~metaj$RoomJ247wk.FLOR, add=T, col=rgb(0,0,0, max=255, alpha=0), outline=F) #transparent overlap
beeswarm(metaj$fD2H ~metaj$RoomJ247wk.VF, corral="wrap", pch=19, cex=0.5, pwcol=ggcolor(7)[metaj$seurat_clusters], main="VariFeat2963; J247")
boxplot(metaj$fD2H ~metaj$RoomJ247wk.VF, add=T, col=rgb(0,0,0, max=255, alpha=0), outline=F) #transparent overlap

#(try. 2)
exp.j2 <- exp.so[, rownames(meta[meta$Accession %in% "J247", ])] #[1] 2963  460
exp.rm <- exp.vf[rownames(exp.j2), ] #[1] 2963   24

cor.df <- data.frame(row.names=colnames(exp.j2))
for(i in 1:8){
  wk=paste0(i,"wk.indoor_PCC")
  cor.df[,wk] <- rep(0, dim(cor.df)[1])
}
#[1] 460   8
for (f in rownames(cor.df)){
  f.exp <- exp.j2[,f]
  for (r in c(1:length(exp.rm))){
    pcc <- cor.test(f.exp, exp.rm[,r], method="pearson")[4]
    wk_n = as.integer((r+2)/3) 
    if(cor.df[f, wk_n] < pcc) {cor.df[f, wk_n] = pcc} 
  }
}
meta2 <- merge(meta, cor.df, by=0, all=TRUE)
rownames(meta2) <- meta2$Row.names
meta2 <- meta2[,-1]
so@meta.data <- meta2
saveRDS(so, "seurat.barley-field.leaf.v2.rds")
#(visualization)
m.pcc <- so@meta.data[, c(14, 38:45)]
m.pcc <- na.omit(m.pcc) #[1] 460   9
m.pcc <- m.pcc[order(m.pcc$week), ]
m.wk <- m.pcc[,1, drop=FALSE]
m.wk$week <- as.numeric(as.vector(m.wk$week))
m.pcc <- as.matrix(m.pcc[,-1])

colnames(m.pcc) <- c("1wk","2wk","3wk","4wk","5wk","6wk","7wk","8wk")
key.palette <- hcl.colors(100, "Plasma")
row.annt <- m.wk
annt.color <- list(week=rev(hcl.colors(20, "Greens")))
p <- pheatmap(as.matrix(m.pcc), cluster_cols = FALSE, cluster_rows = FALSE
              , show_rownames = FALSE, show_colnames = TRUE
              , annotation_row = row.annt, color = key.palette, annotation_colors = annt.color)



#################################################################
############# PHATE trial
#################################################################
require(phateR)
all <- t(as.data.frame(so@assays$RNA@data))
vf <- all[,VariableFeatures(so)]
fl <- all[, colnames(all) %in% flor]

all.phate1 <- phate(all, gamma=1)
all.phate2 <- phate(all, gamma=1, t=120)
all.phate3 <- phate(all, gamma=1, t=60)
all.phate4 <- phate(all, gamma=1, t=120, knn=20)
all.phate5 <- phate(all, gamma=-1, t=120, knn=20)
plot(all.phate5, main = "All-phate gamma=1 t=60 knn=5")

vf.phate1 <- phate(vf, gamma=1)
vf.phate2 <- phate(vf, gamma=1, t=120)
vf.phate3 <- phate(vf, gamma=1, t=60)
vf.phate4 <- phate(vf, gamma=1, t=120, knn=20)
vf.phate5 <- phate(vf, gamma=-1, t=120, knn=20)
vf.phate6 <- phate(vf, gamma=1, t=40)
plot(vf.phate2, main = "VF-phate gamma=1 t=120 knn=5")

fl.phate1 <- phate(fl, gamma=1)
fl.phate2 <- phate(fl, gamma=1, t=120)
fl.phate3 <- phate(fl, gamma=1, t=60)
fl.phate4 <- phate(fl, gamma=1, t=120, knn=20)
fl.phate5 <- phate(fl, gamma=-1, t=120, knn=20)
fl.phate6 <- phate(fl, gamma=1, t=40)
plot(vf.phate2, main = "VF-phate gamma=1 t=120 knn=5")

so[["phate1"]] = CreateDimReducObject(embeddings = all.phate1$embedding * 1000, key = "PHATE_", assay = DefaultAssay(so))
so[["phate2"]] = CreateDimReducObject(embeddings = all.phate2$embedding * 1000, key = "PHATE_", assay = DefaultAssay(so))
so[["phate3"]] = CreateDimReducObject(embeddings = all.phate4$embedding * 1000, key = "PHATE_", assay = DefaultAssay(so))
so[["phate_v1"]] = CreateDimReducObject(embeddings = vf.phate4$embedding * 1000, key = "PHATE_", assay = DefaultAssay(so))
so[["phate_v2"]] = CreateDimReducObject(embeddings = vf.phate6$embedding * 1000, key = "PHATE_", assay = DefaultAssay(so))
so[["phate_f1"]] <- CreateDimReducObject(embeddings = fl.phate1$embedding * 1000, key = "PHATE_", assay = DefaultAssay(so))
so[["phate_f2"]] <- CreateDimReducObject(embeddings = fl.phate4$embedding * 1000, key = "PHATE_", assay = DefaultAssay(so))
so[["phate_f3"]] <- CreateDimReducObject(embeddings = fl.phate3$embedding * 1000, key = "PHATE_", assay = DefaultAssay(so))

DimPlot(so, reduction="phate_f2", group.by=c("Accession","month","Batch","seurat_clusters"))
FeaturePlot(so, reduction="phate_f2", features="fD2H", split.by="Accession")




#############################  (Post Process I: choose the best z score)
prefix2 <- strsplit(prefix, "-")[[1]][1]
subpath <- list.files(".", pattern=prefix2)
rss.files <- file.path(subpath, list.files(subpath, pattern="\\RSS.csv"))
cor.files <- file.path(subpath, list.files(subpath, pattern="\\A-cor.csv"))
rss.df <- data.frame(matrix(nrow=length(read.csv(rss.files[1],header=T)[,1])))
for (f in rss.files){
  rss <- as.vector(read.csv(f,header=T))
  header <- strsplit(f,"-")[[1]][3]
  rss.df[,header] <- rss
}
rss.df <- rss.df[,c(3,4,5,6,2)] #remove the unused col#1 & reorder
#       z2       z4       z6       z8      z10
#1 107945.6 94299.96 91753.77 91324.07 91210.72
#2 107955.7 94004.40 91546.65 91313.21 91179.87
#3 107996.2 93938.25 91506.16 91273.02 91160.86
cor.df <- data.frame(matrix(nrow=length(read.csv(cor.files[1],header=T)[,1])))
for (f in cor.files){
  cor <- as.vector(read.csv(f, header=T)$value)
  header <- strsplit(f,"-")[[1]][3]
  cor.df[,header] <- cor
}
cor.df <- cor.df[,c(3,4,5,6,2)] #remove the unused col#1 & reorder
#         z2        z4        z6        z8       z10
#1 0.9999939 0.9999629 0.9948352 0.8981907 0.8560561
#2 0.9999939 0.9999629 0.9948352 0.8981907 0.8560561
#3 0.9999623 0.9999483 0.9945295 0.8948109 0.8435068
#4 0.9999623 0.9999483 0.9945295 0.8948109 0.8435068
pdf_f <- paste0(prefix2,".ssr-cor.pdf")
pdf(file=pdf_f, width=6, height=6)
par(mfrow = c(1, 2))
boxplot(rss.df, ylab="Squared-Sum of Residue", xlab="Given-z (D)", 
        main=prefix)
boxplot(cor.df, ylab="correlation of A", xlab="Given-z (D)", 
        main="Top-10 Pearson-corr")
par(mfrow = c(1, 1))
dev.off()

#(Post Process II: output network)
subpath2 <- subpath[4] #z6 for LM, L0, fD, CT, HMDEG-fD HMDEG-LM
rss.files <- file.path(subpath2, list.files(subpath2, pattern="\\.A.csv"))
best.idx <- order(rss.df$z6)[1:20] ## 20 minimum RSS values index
for (idx in best.idx){
  header <- paste0("rpt",idx)
  edge.wide <- read.csv(edge.files[best.idx[1]], header=T)
  edge.long <- gather(edge.wide, key="g2", value=header, 2:length(edge.wide))
  if (idx == best.idx[1]){
    edge.table <- edge.long
  } else {
    edge.table[,header] <- edge.long[,3]
  }
}
edge.mean <- edge.table[,c(1,2)]
edge.mean$mean <- rowMeans(edge.table[,c(3:length(edge.table))])
edge.mean <- edge.mean[order(abs(edge.mean$mean), decreasing = TRUE),]
rownames(edge.mean) <- NULL ##reset the row numbers
colnames(edge.mean) <- c("TF","target","mean")
hist(abs(edge.mean[,3])) #(first.view )

edge.top <- edge.mean[c(1:200),] ##top-200 
edge.pos <- subset(edge.top, edge.top$mean > 0)
edge.neg <- subset(edge.top, edge.top$mean < 0)
calc.node.freq <- function(df){
  #node.vec <- c(df[,1], df[,2])
  node.vec <- as.vector(df[,1])
  node.freq <- as.data.frame(table(node.vec))
  node.freq <- node.freq[order(node.freq$Freq, decreasing=TRUE),]
  rownames(node.freq) <- node.freq$node.vec
  return(node.freq)
} 

edge.top.freq <- calc.node.freq(edge.top)
edge.pos.freq <- calc.node.freq(edge.pos)
edge.neg.freq <- calc.node.freq(edge.neg)
edge.top.freq$pos <- edge.pos.freq[rownames(edge.top.freq),]$Freq
edge.top.freq$neg <- edge.neg.freq[rownames(edge.top.freq),]$Freq
edge.top.freq[is.na(edge.top.freq)] <- 0 #convert NA to 0
require(ggplot2)
top.wide <- edge.top.freq[,c(1,3,4)]
top.long <- gather(top.wide, key="np", value="count", 2:length(top.wide))
pdf_f <- paste0(subpath2,".top-nodes.pdf")
pdf(file=pdf_f, width=4.8, height=6)
g <- ggplot(top.long, aes(x=factor(node.vec, level=rev(top.wide$node.vec)), y=count, fill=np))
g <- g + geom_bar(position="stack", stat="identity") + coord_flip() + theme_bw()
g <- g + scale_fill_manual(values = c("#000000", "#FFFFFF")) + geom_col(color = "black", linewidth=0.2)
g + ggtitle("top 200-edges") + xlab("top TFs")
dev.off()

nwk_f <- paste0(subpath2, ".nwk.csv")
edge.sub <- subset(edge.mean, abs(edge.mean[,3]) >= 0.1) #significant edges only
write.csv(edge.sub, file=nwk_f, row.names = F)

#################################################################
############# Linear regression of given genes by accession
#################################################################

goi <- read.csv("02_cluster-DEG/cl2vs3.markers.MAST-dn.GO-BP.csv", header=T)
goi <- strsplit(goi$GeneId[1], "/")[[1]]
goi.x <- as.data.frame(t(scale(all.rpm[goi,])))
goi.x <- goi.x[rownames(meta),]
goi.x <- cbind(goi.x, meta)
LMfit <- function(exp.df, time, col){
  exp <- exp.df[, as.integer(col)]
  #exp <- scale(exp)
  time <- exp.df[, time]
  
  fit = lm(exp~time + I(time^2))
  pred <- predict(fit)
  std.diff <- 1 - pred[1]
  #pred = pred + std.diff
  
  out.df <- data.frame(exp, pred, time)
  return(out.df)
}
PlotLM <- function(exp.df, time, col){
  exp.df <- exp.df[order(exp.df[,time]),]
  j247x <- subset(exp.df, exp.df$Accession == "J247")
  j064x <- subset(exp.df, exp.df$Accession == "J064")
  j647x <- subset(exp.df, exp.df$Accession == "J647")
  h602x <- subset(exp.df, exp.df$Accession == "H602")
  title = paste0(colnames(exp.df[col]), "_",time)
  
  j247m <- LMfit(j247x, time, col)
  j064m <- LMfit(j064x, time, col)
  j647m <- LMfit(j647x, time, col)
  h602m <- LMfit(h602x, time, col)
  
  pred.max <- max(c(j247m$pred, j064m$pred, j647m$pred, h602m$pred))
  pred.min <- min(c(j247m$pred, j064m$pred, j647m$pred, h602m$pred))
  
  matplot(j247m$time, j247m$pred, "l", ylim=c(pred.min,pred.max), lwd=2 , main=title,col=ggcolor(4)[1])
  lines(j064m$time, j064m$pred, lwd=2 ,col=ggcolor(4)[2])
  lines(j647m$time, j647m$pred, lwd=2 ,col=ggcolor(4)[3])
  lines(h602m$time, h602m$pred, lwd=2 ,col=ggcolor(4)[4])
  abline(h=0, lty=2, col="grey")
}
PlotLM(goi.x, "DayHour", 1)

goi <- go.df



#################################################################
############# Figure preparation
#################################################################

###(ChIP-data visualization)
df <- read.table(pipe("pbpaste"), header=T)
#Accs month Prom.2k genebody  ter.2k
# 1  H602     0 4088.13  2623.22 5084.82
# 2  H602     0 1183.89  3321.30 2971.55
# 3  H602     1 3773.71  3027.20 4812.35
# 4  H602     1  908.28  1646.08 3285.77
# 5  H602     2 3382.38  1413.46 6451.83
# 6  H602     2 6441.86  3113.12 5859.43
# 7  H602     3 3759.03  2810.80 6443.56
# 8  H602     3 2558.41  4209.74 6336.28
# 9  H602     4  975.45  3931.38 4876.95
# 10 H602     4 1407.63  3531.98 5073.29
# 11 J247     0 2983.00  1347.63 3523.97
# 12 J247     0 1922.40  3937.51 8356.18
# 13 J247     1 2500.76  2496.28 4633.98
# 14 J247     1 1503.67  2292.92 3378.22
# 15 J247     2 1015.20  1676.55 2538.03
# 16 J247     2    0.00  1655.92 2914.19
# 17 J247     3  456.59  1507.21 2049.88
# 18 J247     3 2310.52  1467.10 1540.17
# 19 J247     4  687.30  1040.49 1356.29
# 20 J247     4 2204.35   609.01 2211.73

doChIPLP1 <- function(df, title){
ylim = c(min(df$genebody), max(df$genebody))
ylim = c(0, max(df$genebody))
h602 <- df[df$Accs %in% "H602", c("month","genebody")]
j247 <- df[df$Accs %in% "J247", c("month","genebody")]
h602m <- tapply(h602$genebody, h602$month, FUN=mean)
j247m <- tapply(j247$genebody, j247$month, FUN=mean)
xax <- as.numeric(names(h602m))

h602.1 <- tapply(h602$genebody, h602$month, FUN=min)
h602.2 <- tapply(h602$genebody, h602$month, FUN=max)
j247.1 <- tapply(j247$genebody, j247$month, FUN=min)
j247.2 <- tapply(j247$genebody, j247$month, FUN=max)

pdf_f = paste0(title,".lineplot1.pdf")
pdf(file=pdf_f, height=5, width=3)
matplot(y=h602m, x=names(h602m),  type="l", ylim=c(0,2.5), lwd=3, col=ggcolor(4)[4], xlab="month", ylab="RPM/nucleotide", main=title)
matplot(y=h602$genebody, x=h602$month, type="p", pch=19, cex=1, col=ggcolor(4)[4], add=T)
arrows(xax, h602.1, xax, h602.2, code=3, lwd=1, angle=90, length=.07, col=ggcolor(4)[4])
matplot(y=j247m, x=names(j247m),type="l", pch=19, cex=1, lwd=3, col=ggcolor(4)[3], add=T)
matplot(y=j247$genebody, x=j247$month, type="p", pch=19, cex=1, col=ggcolor(4)[3], add=T)
arrows(xax, j247.1, xax, j247.2, code=3, lwd=1, angle=90, length=.07, col=ggcolor(4)[3])
dev.off()

}
doChIPLP1(df, "HvFUL3_H3K36me3_genebody")


###(GAM-fitted line plot with given gene-id)
mtx <- as.matrix(gam_lm.j2) - apply(gam_lm.j2, 1, mean)
genes <- c("HORVU7Hr1G024610","HORVU4Hr1G077450","HORVU5Hr1G095630","HORVU2Hr1G063800","HORVU0Hr1G003020","HORVU3Hr1G095240","HORVU2Hr1G019900" )
gene <- "HORVU3Hr1G027590"

ymax = max(mtx[genes,])
ymin = min(mtx[genes,])
plot(as.vector(mtx[gene, c(1:100)]), type="l", col=ggcolor(4)[3], lwd=2, ylim=c(ymin, ymax),main=gene, xlab="GAM-linM", ylab="z-RPM2")
lines(as.vector(mtx[gene, c(101:200)]), type="l")
abline(h=0, lty=2, col="grey")



key.palette <- rev(hcl.colors(100, "Lajolla"))
FeaturePlot(so, reduction="umap", features="DayHour", col=key.palette)
#save PDF as 7x6 inch
key.palette <- hcl.colors(100, "Inferno")
FeaturePlot(so, reduction="umap", features="AvgTemp", col=key.palette)

FeaturePlot(so, reduction="umap", features="fD2H")

DimPlot(so, reduction="umap", group.by="Accession")
DimPlot(so, reduction="umap", group.by="month", cols=rev(ggcolor(6))[order(levels(so@meta.data$month))])
DimPlot(so, reduction="umap", group.by="Case", cols=hcl.colors(4, "Fall"))
DimPlot(so, label=TRUE, label.size=15)

key.palette <- hcl.colors(100, "Mako")
FeaturePlot(so, reduction="umap_2D", features="LineageM", col=key.palette)
#save PDF as 7x6 inch
col7 <- c()
for (i in levels(meta$seurat_clusters)){
  sub <- meta[, c("seurat_clusters","LineageM")]
  min <- min(sub[,2], na.rm=TRUE)
  max <- max(sub[,2], na.rm=TRUE)
  sub <- sub[sub$seurat_clusters %in% i, 2]
  p <- punif(median(sub, na.rm=TRUE) , min, max) ### output percentile of the given range
  p <- as.integer(p*100)
  col7 <- c(col7, p)
}
boxplot(meta$LineageM~meta$seurat_clusters, col=key.palette[col7], pch=19, cex=0.5, notch=T)
#save PDF as 6x6 inch


key.palette <- hcl.colors(100, "Viridis")
FeaturePlot(so, reduction="umap", features="Lineage0", col=key.palette)

metaj2 <- meta[meta$Accession %in% "J247",]
metaj6 <- meta[meta$Accession %in% "J647",]
metaj0 <- meta[meta$Accession %in% "J064",]
metah6 <- meta[meta$Accession %in% "H602",]
col6 <- hcl.colors(6, "Viridis")
plot(metaj2$fD2H~metaj2$Lineage0, pch=19, col=col6[metaj2$month], xlim=c(0,28), ylim=c(0,1.2), main="J247")
plot(metaj0$fD2H~metaj0$Lineage0, pch=19, col=col6[metaj0$month], xlim=c(0,28), ylim=c(0,1.2), main="J064")
plot(metaj6$fD2H~metaj6$Lineage0, pch=19, col=col6[metaj6$month], xlim=c(0,28), ylim=c(0,1.2), main="J647")
plot(metah6$fD2H~metah6$Lineage0, pch=19, col=col6[metah6$month], xlim=c(0,28), ylim=c(0,1.2), main="H602")
#save PDF as 4.5x6 inch
plot(metaj2$fD2H~metaj2$LineageM, pch=19, col=ggcolor(4)[1], xlim=c(0,28), ylim=c(0,1.2), main="J247")
plot(metaj0$fD2H~metaj0$LineageM, pch=19, col=ggcolor(4)[2], xlim=c(0,28), ylim=c(0,1.2), main="J064")
plot(metaj6$fD2H~metaj6$LineageM, pch=19, col=ggcolor(4)[3], xlim=c(0,28), ylim=c(0,1.2), main="J647")
plot(metah6$fD2H~metah6$LineageM, pch=19, col=ggcolor(4)[4], xlim=c(0,28), ylim=c(0,1.2), main="H602")
#save PDF as 4.5x6 inch


VlnPlot(so, features="CytoTRACE", group.by="Accession")
VlnPlot(so, features="Lineage0", group.by="Accession")
VlnPlot(so, features="LineageM", group.by="Accession")
VlnPlot(so, features="fD2H", group.by="Accession")
#save PDF as 5x6 inch


col7 <- c()
for (i in levels(meta$seurat_clusters)){
  sub <- meta[meta$seurat_clusters %in% i, "fD2H"]
  m <- median(sub)
  m <- as.integer(m*100)
  col7 <- c(col7, m)
}
key.palette <- rev(hcl.colors(100, "Blues 2"))
boxplot(meta$fD2H~meta$seurat_clusters, col=key.palette[col7], pch=19, cex=0.5, notch=T)
#save PDF as 6x6 inch
col7 <- c()
for (i in levels(meta$seurat_clusters)){
  sub <- meta[, c("seurat_clusters","DayHour")]
  min <- min(sub[,2])
  max <- max(sub[,2])
  sub <- sub[sub$seurat_clusters %in% i, 2]
  p <- punif(median(sub) , min, max)   ### output percentile of the given range
  p <- as.integer(p*100)
  col7 <- c(col7, p)
}
key.palette <- rev(hcl.colors(100, "Lajolla"))
boxplot(meta$DayHour~meta$seurat_clusters, col=key.palette[col7], pch=19, cex=0.5, notch=T)
#save PDF as 6x6 inch
col7 <- c()
for (i in levels(meta$seurat_clusters)){
  sub <- meta[, c("seurat_clusters","AvgTemp")]
  min <- min(sub[,2])
  max <- max(sub[,2])
  sub <- sub[sub$seurat_clusters %in% i, 2]
  p <- punif(median(sub) , min, max) ### output percentile of the given range
  p <- as.integer(p*100)
  col7 <- c(col7, p)
}
key.palette <- hcl.colors(100, "Inferno")
boxplot(meta$AvgTemp~meta$seurat_clusters, col=key.palette[col7], pch=19, cex=0.5, notch=T)
#save PDF as 6x6 inch
col7 <- c()
for (i in levels(meta$seurat_clusters)){
  sub <- meta[, c("seurat_clusters","Lineage0")]
  min <- min(sub[,2])
  max <- max(sub[,2])
  sub <- sub[sub$seurat_clusters %in% i, 2]
  p <- punif(median(sub) , min, max) ### output percentile of the given range
  p <- as.integer(p*100)
  col7 <- c(col7, p)
}
key.palette <- hcl.colors(100, "Viridis")
boxplot(meta$Lineage0~meta$seurat_clusters, col=key.palette[col7], pch=19, cex=0.5, notch=T)
#save PDF as 5x6 inch
col7 <- c()
for (i in levels(meta$seurat_clusters)){
  sub <- meta[, c("seurat_clusters","CytoTRACE.VF")]
  min <- min(sub[,2])
  max <- max(sub[,2])
  sub <- sub[sub$seurat_clusters %in% i, 2]
  p <- punif(median(sub) , min, max) ### output percentile of the given range
  p <- as.integer(p*100)
  col7 <- c(col7, p)
}
key.palette <- hcl.colors(100, "Oslo")
boxplot(meta$CytoTRACE.VF~meta$seurat_clusters, col=key.palette[col7], pch=19, cex=0.5, notch=T)
#save PDF as 5x6 inch
FeaturePlot(so, reduction="umap_2D", features="CytoTRACE.VF", col=key.palette)
#save PDF as 6x5 inch


meta.j2 <- subset(meta, meta$Accession == "J247")
beeswarm(meta.j2$fD2H ~ meta.j2$RoomJ247wk.VF, pch=19, corral="wrap", cex=0.5, col=rev(hcl.colors(16, "Emrld"))[meta.j2$RoomJ247wk.VF])
boxplot(meta.j2$fD2H ~ meta.j2$RoomJ247wk.VF, col=rgb(1,1,1, alpha=0), add=T)
#save PDF as 6x6 inch
table(meta.j2$RoomJ247wk.VF)
#1   2   3   4   5   6   7   8 
#63   0   0   0  11 122 212  52 
beeswarm(meta.j2$RoomJ247wk.VF.PCC~meta.j2$RoomJ247wk.VF, corral="wrap", pch=19, cex=0.5, ylim=c(0,1))
boxplot(meta.j2$RoomJ247wk.VF.PCC~meta.j2$RoomJ247wk.VF, col=rgb(1,1,1, alpha=0), add=T, outline=F)



so23 <- subset(so, idents=c(2,3))
so45 <- subset(so, idents=c(4,5))
DimPlot(so23, reduction="umap", cols=ggcolor(7)[c(2,3)])

col7 <- rep("lightgrey", times=7)
col7[2:3] <- ggcolor(7)[2:3]
plot(meta$AvgTemp~meta$DayHour, pch=19, col=col7)


## 3-D UMAP
require(rgl)
plot3d(so@reductions$umap_3D@cell.embeddings, xlab="UMAP_1", ylab="UMAP_2", zlab="UMAP_3", col=ggcolor(7)[meta$seurat_clusters])
for( i in 0:359 ) {
  rgl.viewpoint( i, i/4 )
  rgl.snapshot( fmt="png", sprintf("so-umap%03d.png", i)  )
}
## generate GIF from PNGs
#(convert *.png screens.gif / OsX terminal)


## metadata correlation
library(psych)
meta.sub <- subset(meta, meta$seurat_clusters %in% c(4,5))
pairs.panels(meta.sub[,c(37,19,21,23)])
x = meta.sub$DayHour
y = meta.sub$CytoTRACE.VF
col5 <- c("grey20","grey20", "grey20", ggcolor(7)[4:5] ,"grey20","grey20")
plot(y ~ x, pch=19, cex=0.5, col=col5[meta$seurat_clusters])
fit = lm(y~x + I(x^2))
pred.df <- data.frame(x, predict(fit))
pred.df <- pred.df[order(pred.df$x),]
lines(pred.df$x, pred.df$pred, lwd=2 ,col=2) 
#save PDF as 5x6 inch

x=meta$fD2H
y=meta$RoomJ247wk.VF
plot(x=x, y=y, type="p", pch=19)
y= as.numeric(y)
fit = lm(y~x + I(x^2))
pred.df <- data.frame(x, predict(fit))
pred.df <- pred.df[order(pred.df$x),]
lines(pred.df$x, pred.df$pred, lwd=2 ,col=2) 

x <- meta23$AvgTemp
y <- meta23$DayHour
plot(x=meta$AvgTemp, y=meta$DayHour,xlab="Air Temp", ylab="Day Hour", col="grey", pch=19, cex=0.5)
points(x=x, y=y, pch=19, col=ggcolor(7)[meta23$seurat_clusters])
fit = lm(y~x + I(x^2))
pred.df <- data.frame(x, predict(fit))
pred.df <- pred.df[order(pred.df$x),]
lines(pred.df$x, pred.df$pred, lwd=2 ,col=2) 
x2 = meta$AvgTemp
y2 = meta$DayHour
fit2 = lm(y2~x2 + I(x2^2))
pred.df2 <- data.frame(x2, predict(fit2))
pred.df2 <- pred.df2[order(pred.df2$x2),]
lines(pred.df2$x2, pred.df2$predict.fit2., lwd=2 ,col="grey") 



#########################################
########## RRHO trial
#########################################
##(DEG; P-value test)
cl2s <- cl2[rownames(cl2) %in% j247.deg, ] #[1] 549   5 /length(j247.deg)
cl3s <- cl3[rownames(cl3) %in% j247.deg, ] #[1] 492   5 /length(j247.deg)
dde2 <- data.frame(row.names = rownames(cl2s), Genes=rownames(cl2s), DDE=-log(cl2s$p_val)*sign(cl2s$avg_log2FC) *-1, stringsAsFactors = FALSE)
dde3 <- data.frame(row.names = rownames(cl3s), Genes=rownames(cl3s), DDE=-log(cl3s$p_val)*sign(cl3s$avg_log2FC) *-1, stringsAsFactors = FALSE)
dde2 <- na.omit(dde2) 
dde3 <- na.omit(dde3) 
l1 <- dde2[rownames(dde2) %in% rownames(dde3), ] #1] 382   2
l2 <- dde3[rownames(l1), ] #[1] 382   2
RRHO_obj <-  RRHO2_initialize(l1, l2, labels = c("Cluster2", "Cluster3"), log10.ind=TRUE, boundary=0.05)
RRHO2_heatmap(RRHO_obj)


##(GO; odds ratio test)
go1 <- read.csv("03_TradeSeq-assoc/GAM-lm.J247-DEG1-A.GO-BP.csv", header=T, row.names=1)
go2 <- read.csv("03_TradeSeq-assoc/GAM-lm.J247-DEG1-B.GO-BP.csv", header=T, row.names=1)
dde1 <- data.frame(row.names = rownames(go1), Genes=rownames(go1), DDE=-log10(go1$Pvalue) *sign(log(go1$OddsRatio)), stringsAsFactors = FALSE)
dde2 <- data.frame(row.names = rownames(go2), Genes=rownames(go2), DDE=-log10(go2$Pvalue) *sign(log(go2$OddsRatio)), stringsAsFactors = FALSE)
dde1 <- na.omit(dde1)
dde2 <- na.omit(dde2)
l1 <- dde1[rownames(dde1) %in% rownames(dde2), ] #[1] 76  2
l2 <- dde2[rownames(l1), ]
RRHO_obj2 <-  RRHO2_initialize(l1, l2, labels = c("Cluster2", "Cluster3"), log10.ind=TRUE)#, method = "fisher")
RRHO2_heatmap(RRHO_obj2)
RRHO2_vennDiagram(RRHO_obj2, type="du")
RRHO_obj2$genelist_du
# $gene_list1_du
# [1] "GO:0009791" "GO:0006629" "GO:0005975" "GO:0009056"
# $gene_list2_du
# [1] "GO:0006950"


############### FLowering genes line plot
flor #[1] 150
doLPacc2(flor)

exp <- read.table("00_IPSR_KIBR_Sampledata/RPM_IPSR_KIBR_LB_17_19_BLD20201007.txt", header=T, row.names=1)
newCol <- sapply(colnames(exp), function(x){substr(x, 2, nchar(x))}) #remove the akward X from the colnames
colnames(exp) <- newCol
head(exp, c(5,5))

in.rpm <- exp[, newCol[grep("^5", newCol)]]
doLP5 <- function(goi){
  goi.rpm <- in.rpm[rownames(in.rpm) %in% goi, ]
  for(gid in rownames(goi.rpm)){
    pdf_f <- paste0(gid,".meanRPM-indoor.pdf")
    new.df <- as.data.frame(matrix(nrow=8, ncol=3))
    new.df[1, ] <- in.rpm[gid, c(1:3)]
    new.df[2, ] <- in.rpm[gid, c(4:6)]
    new.df[3, ] <- in.rpm[gid, c(7:9)]
    new.df[4, ] <- in.rpm[gid, c(10:12)]
    new.df[5, ] <- in.rpm[gid, c(13:15)]
    new.df[6, ] <- in.rpm[gid, c(16:18)]
    new.df[7, ] <- in.rpm[gid, c(19:21)]
    new.df[8, ] <- in.rpm[gid, c(22:24)]
    
    pdf(file=pdf_f, height=5.5, width=5)
    matplot(new.df, type="p", ylim=c(0, max(new.df, na.rm=TRUE)) ,pch=20, cex=1.2, col="#00BFC4a6", main=gid, xlab="weeks", ylab="mean RPM")
    matplot(rowMeans(new.df, na.rm=TRUE), type="b", lwd=3, lty=1, pch=19, col="#00BFC4", add=T)
    dev.off()
  }
}
doLP5(flor)


sessionInfo()
# R version 4.3.1 (2023-06-16)
# Platform: x86_64-apple-darwin20 (64-bit)
# Running under: macOS Monterey 12.7.4
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
#   [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] GO.db_3.18.0         scales_1.3.0         reshape2_1.4.4       stringr_1.5.1        MASS_7.3-60          beeswarm_0.4.0      
# [7] dunn.test_1.3.5      KEGGREST_1.42.0      GOstats_2.68.0       Category_2.68.0      GSEABase_1.64.0      graph_1.80.0        
# [13] annotate_1.80.0      XML_3.99-0.16        AnnotationDbi_1.64.1 IRanges_2.36.0       S4Vectors_0.40.2     Biobase_2.62.0      
# [19] BiocGenerics_0.48.1  glmnet_4.1-8         Matrix_1.6-4         gplots_3.1.3.1       cluster_2.1.6        Seurat_5.0.3        
# [25] SeuratObject_5.0.1   sp_2.1-3            
