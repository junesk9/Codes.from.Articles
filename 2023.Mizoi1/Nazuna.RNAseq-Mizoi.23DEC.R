############## Nazna RNA-seq analysis
############## by Prof. Mizoi, Todai
############## 2023.12.24 Junesk9


library(tximport)
library(DESeq2)
library(tidyverse)
library(pheatmap)

#for GO analysis
library(AnnotationForge)
library(org.At.tair.db)
library(GSEABase)
library(GOstats)
library(KEGGREST)
#library(pathviewr)


setwd("~/Dropbox (理化学研究所　セルロース生産研究チーム)/Collaborations/23DEC_溝井ナズナRNAseq")
set.seed(101) #for data reproduciblity

############## Load RSEM data
rsem.files <- file.path("0_rsem", list.files("0_rsem", ".isoforms.results"))
f <- function(s) unlist(strsplit(s, "\\/|\\."))[2]
sample.name <- sapply(rsem.files, f)
sample.name = as.vector(sample.name)
names(rsem.files) <- sample.name

# Load Exp data
tx.exp <- tximport(rsem.files, type="rsem", txIn=TRUE, txOut=TRUE)
tx2gene <- data.frame(TXNAME = rownames(tx.exp$counts))
f2 <- function(s) unlist(strsplit(s, "[.]"))[1]
tx2gene$GENEID = unlist(lapply(tx2gene[,1], f2))
#  TXNAME    GENEID
#1 AT1G01010.1 AT1G01010
#2 AT1G01020.1 AT1G01020
#3 AT1G01020.2 AT1G01020
#4 AT1G01020.3 AT1G01020
gn.exp <- summarizeToGene(tx.exp, tx2gene, countsFromAbundance = "scaledTPM")
tpm.df <- as.data.frame(gn.exp$abundance)
cnt.df <- as.data.frame(gn.exp$counts)
#[1] 32833    42
write.csv(tpm.df, file = "gene_tpm.star.whole.csv", row.names = TRUE, quote=FALSE)
write.csv(cnt.df, file = "gene_count.star.whole.csv", row.names = TRUE, quote=FALSE)

#prepare sample metadata
sample.table <- data.frame(row.names=sample.name, rpt=sample.name)
sample.table$rpt <- sapply(sample.name, function(x)strsplit(x, "-")[[1]][3])
sample.table$geno <- sapply(sample.name, function(x)strsplit(x, "-")[[1]][1])
sample.table$time <- sapply(sample.name, function(x)strsplit(x, "-")[[1]][2])
sample.table$geno <- factor(sample.table$geno, levels=c("wt","ppk"))
sample.table$hour <- factor(rep(c(0,2,4,5,6,5,6), each=3, times=2))
#          rpt geno time hour
#ppk-0-1    1  ppk    0    0
#ppk-0-2    2  ppk    0    0
#ppk-0-3    3  ppk    0    0
#ppk-h2-1   1  ppk   h2    2
#ppk-h2-2   2  ppk   h2    2
#ppk-h2-3   3  ppk   h2    2
write.csv(sample.table, "sample.meta.csv", quote=FALSE)

#prepare gene annt data
annt.table <- read.csv("./z_useful.info/TAIR10.gene_annt.csv", header=T, row.names=1)
colnames(annt.table) <- c("symbol","annt")
rownames(annt.table) <- toupper(rownames(annt.table)) #for ensuring

#prepare GO DB
frame <- toTable(org.At.tairGO)
goframeData <- data.frame(frame$go_id, frame$Evidence, frame$gene_id)
goFrame <- GOFrame(goframeData, organism = "Arabidospis thaliana")
goAllFrame <- GOAllFrame(goFrame)
gsc <- GeneSetCollection(goAllFrame, setType = GOCollection())
all <- Lkeys(org.At.tairGO)
#[1] 27416

############ Inter-sample relationship
#dendrogram
d1 <- as.dist(1-cor(tpm.df,  method="pearson"))
d2 <-  as.dist(t(1-cor(tpm.df,  method="pearson")))
c1 <- hclust(d1, method="complete")
c2 <- hclust(d2, method="complete")
plot(c1, main="Mizoi At-RNAseq-complete") #dendrogram
#draw heatmap
heatmap(as.matrix(d1,d2),
        Colv=as.dendrogram(c2), Rowv=as.dendrogram(c1),
        scale="none", col=cm.colors(256),
        main="Mizoi At-RNAseq-complete", margin=c(6,6))
require(pheatmap)
pheatmap(as.matrix(d1,d2), Colv=as.dendrogram(c2), Rowv=as.dendrogram(c1), main="Mizoi At-RNAseq-complete")

#PCA
#(data prep)
tpm.pca <- na.omit(tpm.df)
tpm.pca <- tpm.pca[which(apply(t(tpm.pca), 2, var) != 0), ] #remove no changing gene
#[1] 28844    42
tpm.pca2 <- tpm.pca[rowSums(tpm.pca) >= length(colnames(tpm.pca)), ] #avg TPM >= 1
#[1] 19693    42

#(color palette)
sample.idx = rep(1:14, each=3)
brl14 <- rep(hcl.colors(7, palette = "Berlin"), times=2)

#(run PCA & plot; case I)
rpca1 <- prcomp(x=t(tpm.pca), scale=T) 
PropVar <- summary(rpca1)$importance[2,]
barplot(PropVar[1:10], ylab="Proportion of Variance", main="All Dynamic genes")
write.csv(rpca1$x, "PCA-All.csv", quote=FALSE)
xlab = paste0("PC1: ",round(PropVar[1]*100,2), "%")
ylab = paste0("PC2: ",round(PropVar[2]*100,2), "%")
plot(rpca1$x[,1], rpca1$x[,2], xlab=xlab, ylab=ylab, cex=1.3, pch=rep(c(16,17), each=21),
     col=brl14[sample.idx], main="PCA All Dynamic genes")
#(labeling)
avg.pc1 <- c()
avg.pc2 <- c()
for(i in seq(from=1, to=42-1, by=3)){
  pc1 <- rpca1$x[,1]
  pc2 <- rpca1$x[,2]
  i2 = i + 2
  avg1 <- mean(pc1[i:i2])
  avg2 <- mean(pc2[i:i2])
  avg.pc1 <- c(avg.pc1, avg1)
  avg.pc2 <- c(avg.pc2, avg2)
}
label.df <- data.frame(avg.pc1, avg.pc2, sample.name[seq(from=1, to=42-1, by=3)])
text(label.df[,1], label.df[,2]-5, label.df[,3], col=brl14)

#(run PCA & plot; case II)
rpca1 <- prcomp(x=t(tpm.pca2), scale=T) 
PropVar <- summary(rpca1)$importance[2,]
barplot(PropVar[1:10], ylab="Proportion of Variance", main="Dynamic >1-TPM genes")
write.csv(rpca1$x, "PCA-TPM1.csv", quote=FALSE)
xlab = paste0("PC1: ",round(PropVar[1]*100,2), "%")
ylab = paste0("PC2: ",round(PropVar[2]*100,2), "%")
plot(rpca1$x[,1], rpca1$x[,2], xlab=xlab, ylab=ylab, cex=1.3, pch=rep(c(16,17), each=21),
     col=brl14[sample.idx], main="PCA Dynamic >1-TPM genes")
#(labeling)
avg.pc1 <- c()
avg.pc2 <- c()
for(i in seq(from=1, to=42-1, by=3)){
  pc1 <- rpca1$x[,1]
  pc2 <- rpca1$x[,2]
  i2 = i + 2
  avg1 <- mean(pc1[i:i2])
  avg2 <- mean(pc2[i:i2])
  avg.pc1 <- c(avg.pc1, avg1)
  avg.pc2 <- c(avg.pc2, avg2)
}
label.df <- data.frame(avg.pc1, avg.pc2, sample.name[seq(from=1, to=42-1, by=3)])
text(label.df[,1], label.df[,2]-5, label.df[,3], col=brl14)

############ DEG analysis
## Prep. dds object
dds <- DESeqDataSetFromTximport(gn.exp, sample.table, ~geno + hour + geno:hour)
saveRDS(dds, "dds-mizoi.rds")

## Case I : time 0
sample.sub1 <- subset(sample.table, sample.table$time == "0")
Deseq_geno <- function(sample.sub){
  dds.sub <- dds[,colnames(dds) %in% rownames(sample.sub)]
  dds.sub <- DESeqDataSet(dds.sub, ~geno)
  dds.sub <- DESeq(dds.sub)
  res.sub <- results(dds.sub)
  sig.sub <- subset(res.sub, abs(res.sub$log2FoldChange) >= 1 & res.sub$padj < 0.05)
  print(dim(sig.sub))
  
  out.sub <- as.data.frame(sig.sub[,c(1,2,6)])
  out.sub <- merge(out.sub, tpm.df[, rownames(sample.sub)], by=0)
  colnames(out.sub)[1:4] <- c("gid","MeanEXP", "log2FC", "FDR")
  out.sub <- out.sub[order(out.sub$log2FC, decreasing=TRUE), ]
  out.sub$annt <- annt.table[toupper(out.sub$gid), ]$annt
  
  out.ls <- list(res.sub, out.sub)
  return(out.ls)
}
out.ls1 <- Deseq_geno(sample.sub1)
res.sub1 <- out.ls1[[1]]
out.sub1 <- out.ls1[[2]]
#[1] 1936    6
write.csv(out.sub1, "DEG.t0-ppk-wt.csv", quote=TRUE, row.names=FALSE)
PlotVC <- function(df){
  vc <- df[,c("log2FoldChange","padj")]
  colnames(vc) <- c("Log2FC","qval")
  v1 <- subset(vc, vc$qval < 0.05 & vc$Log2FC >= 1)
  v2 <- subset(vc, vc$Log2FC <= -1 & vc$qval < 0.05)
  deg_ids <- c(rownames(v1), rownames(v2))
  `%notin%` <- Negate(`%in%`) #temporary designates "not-in" func
  v3 <- vc[rownames(vc) %notin% deg_ids, ]
  up.len <- length(rownames(v1))
  dn.len <- length(rownames(v2))
  nd.len <- length(rownames(v3))
  print(up.len)
  print(dn.len)
  print(nd.len)
  
  vc$qval <- -log(vc$qval, 10)
  vc$qval <- round(vc$qval, digits=1)
  vc$Log2FC <- round(vc$Log2FC, digits=1)
  v1 <- unique(vc[rownames(v1), ])
  v2 <- unique(vc[rownames(v2), ])
  v3 <- unique(vc[rownames(v3), ])
  qval <- vc$qval[!is.infinite(vc$qval)]
  Log2FC <- vc$Log2FC[!is.infinite(vc$Log2FC)]
  max_q <- max(abs(qval), na.rm=TRUE)
  max_fc <- max(abs(Log2FC), na.rm=TRUE)
  
  plot(v3$Log2FC,v3$qval, pch=20, cex=0.4, col="grey70", ylim=c(0,max_q), 
       xlim=c(-max_fc,max_fc), xlab="Log2FC", ylab="-log10FDR")
  points(v1$Log2FC,v1$qval, pch=20, cex=0.5, col="#ED0422")
  points(v2$Log2FC,v2$qval, pch=20, cex=0.5, col="#3498DB")
  abline(v=0, lty=2)
}
PlotVC(res.sub1)

GOBP <- function(df){
  gid = df$gid
  up200 <- toupper(gid[1:200])
  dn200 <- toupper(rev(gid)[1:200])
  
  p <- GSEAGOHyperGParams(
    name = "Paramaters",
    geneSetCollection = gsc,
    geneIds = up200,
    universeGeneIds = all,
    ontology = "BP",
    pvalueCutoff = 1,
    conditional = TRUE,
    testDirection = "over"
  )
  up_res <- hyperGTest(p)
  up_gid <- geneIdsByCategory(up_res)
  up_res2 <- summary(up_res)
  up_res2$GeneId <- sapply(up_res2[,1], function(x)paste(unlist(up_gid[x]), collapse="/"))
  
  p <- GSEAGOHyperGParams(
    name = "Paramaters",
    geneSetCollection = gsc,
    geneIds = dn200,
    universeGeneIds = all,
    ontology = "BP",
    pvalueCutoff = 1,
    conditional = TRUE,
    testDirection = "over"
  )
  dn_res <- hyperGTest(p)
  dn_gid <- geneIdsByCategory(dn_res)
  dn_res2 <- summary(dn_res)
  dn_res2$GeneId <- sapply(dn_res2[,1], function(x)paste(unlist(dn_gid[x]), collapse="/"))
  
  gobp <- list(up_res2, dn_res2)
  return(gobp)
}
go.sub1 <- GOBP(out.sub1)
write.csv(go.sub1[[1]], "DEG.t0-ppk-wt.up200.GOBP.csv", quote=TRUE, row.names=FALSE)
write.csv(go.sub1[[2]], "DEG.t0-ppk-wt.dn200.GOBP.csv", quote=TRUE, row.names=FALSE)


## Case II : time h2
sample.sub2 <- subset(sample.table, sample.table$time == "h2")
out.ls2 <- Deseq_geno(sample.sub2)
#[1] 2133    6
res.sub2 <- out.ls2[[1]]
out.sub2 <- out.ls2[[2]]
write.csv(out.sub2, "DEG.t2h-ppk-wt.csv", quote=TRUE, row.names=FALSE)
PlotVC(res.sub2)
go.sub2 <- GOBP(out.sub2)
write.csv(go.sub2[[1]], "DEG.t2h-ppk-wt.up200.GOBP.csv", quote=TRUE, row.names=FALSE)
write.csv(go.sub2[[2]], "DEG.t2h-ppk-wt.dn200.GOBP.csv", quote=TRUE, row.names=FALSE)

## Case III : time h4
sample.sub3 <- subset(sample.table, sample.table$time == "h4")
out.ls3 <- Deseq_geno(sample.sub3)
#3697    6
res.sub3 <- out.ls3[[1]]
out.sub3 <- out.ls3[[2]]
write.csv(out.sub3, "DEG.t4h-ppk-wt.csv", quote=TRUE, row.names=FALSE)
PlotVC(res.sub3)
go.sub3 <- GOBP(out.sub3)
write.csv(go.sub3[[1]], "DEG.t4h-ppk-wt.up200.GOBP.csv", quote=TRUE, row.names=FALSE)
write.csv(go.sub3[[2]], "DEG.t4h-ppk-wt.dn200.GOBP.csv", quote=TRUE, row.names=FALSE)


## Case IV : time time h5
sample.sub4 <- subset(sample.table, sample.table$time == "h5")
out.ls4 <- Deseq_geno(sample.sub4)
#[1] 2428    6
res.sub4 <- out.ls4[[1]]
out.sub4 <- out.ls4[[2]]
write.csv(out.sub4, "DEG.t5h-ppk-wt.csv", quote=TRUE, row.names=FALSE)
PlotVC(res.sub4)
go.sub4 <- GOBP(out.sub4)
write.csv(go.sub4[[1]], "DEG.t5h-ppk-wt.up200.GOBP.csv", quote=TRUE, row.names=FALSE)
write.csv(go.sub4[[2]], "DEG.t5h-ppk-wt.dn200.GOBP.csv", quote=TRUE, row.names=FALSE)


## Case V : time  h6
sample.sub5 <- subset(sample.table, sample.table$time == "h6")
out.ls5 <- Deseq_geno(sample.sub5)
#[1] 971   6
res.sub5 <- out.ls5[[1]]
out.sub5 <- out.ls5[[2]]
write.csv(out.sub5, "DEG.t6h-ppk-wt.csv", quote=TRUE, row.names=FALSE)
PlotVC(res.sub5)
go.sub5 <- GOBP(out.sub5)
write.csv(go.sub5[[1]], "DEG.t6h-ppk-wt.up200.GOBP.csv", quote=TRUE, row.names=FALSE)
write.csv(go.sub5[[2]], "DEG.t6h-ppk-wt.dn200.GOBP.csv", quote=TRUE, row.names=FALSE)


## Case VI : time h4r1
sample.sub6 <- subset(sample.table, sample.table$time == "h4r1")
out.ls6 <- Deseq_geno(sample.sub6)
#[1] 2444    6
res.sub6 <- out.ls6[[1]]
out.sub6 <- out.ls6[[2]]
write.csv(out.sub6, "DEG.t4hr1-ppk-wt.csv", quote=TRUE, row.names=FALSE)
PlotVC(res.sub6)
go.sub6 <- GOBP(out.sub6)
write.csv(go.sub6[[1]], "DEG.t4hr1-ppk-wt.up200.GOBP.csv", quote=TRUE, row.names=FALSE)
write.csv(go.sub6[[2]], "DEG.t4hr1-ppk-wt.dn200.GOBP.csv", quote=TRUE, row.names=FALSE)

## Case VII : time h4r2
sample.sub7 <- subset(sample.table, sample.table$time == "h4r2")
out.ls7 <- Deseq_geno(sample.sub7)
#[1] 2455    6
#-- note: fitType='parametric', but the dispersion trend was not well captured by the
#function: y = a/x + b, and a local regression fit was automatically substituted.
#specify fitType='local' or 'mean' to avoid this message next time.
res.sub7 <- out.ls7[[1]]
out.sub7 <- out.ls7[[2]]
write.csv(out.sub7, "DEG.t4hr2-ppk-wt.csv", quote=TRUE, row.names=FALSE)
PlotVC(res.sub7)
go.sub7 <- GOBP(out.sub7)
write.csv(go.sub7[[1]], "DEG.t4hr2-ppk-wt.up200.GOBP.csv", quote=TRUE, row.names=FALSE)
write.csv(go.sub7[[2]], "DEG.t4hr2-ppk-wt.dn200.GOBP.csv", quote=TRUE, row.names=FALSE)


## Case VIII : timecourse change wt by heat
sample.sub8 <- subset(sample.table, sample.table$geno == "wt")
sample.sub8a <- sample.sub8[-c(10:15),]
sample.sub8b <- sample.sub8[c(1:15),]

Deseq_hour <- function(sample.sub){
  dds.sub <- dds[,colnames(dds) %in% rownames(sample.sub)]
  dds.sub <- DESeqDataSet(dds.sub, ~hour)
  dds.sub <- DESeq(dds.sub)
  res.sub <- results(dds.sub)
  sig.sub <- subset(res.sub, abs(res.sub$log2FoldChange) >= 1 & res.sub$padj < 0.01)
  print(dim(sig.sub))
  
  out.sub <- as.data.frame(sig.sub[,c(1,2,6)])
  out.sub <- merge(out.sub, tpm.df[, rownames(sample.sub)], by=0)
  colnames(out.sub)[1:4] <- c("gid","MeanEXP", "log2FC", "FDR")
  out.sub <- out.sub[order(out.sub$log2FC, decreasing=TRUE), ]
  out.sub$annt <- annt.table[toupper(out.sub$gid), ]$annt
  
  out.ls <- list(res.sub, out.sub)
  return(out.ls)
}
out.ls8a <- Deseq_hour(sample.sub8a)
#[1] 8457    6
out.sub8a <- out.ls8a[[2]]
write.csv(out.sub8a, "DEG.wt-timecourse-h0h6.csv", quote=TRUE, row.names=FALSE)
PlotHM <- function(df){
  df <- df[order(df$FDR),]
  df <- df[c(1:500),-c(1:4,length(df))] #top500 FDR
  rep = 3

  hrs <- colnames(df)
  hrs <- sapply(hrs, function(x) paste(strsplit(x,"-")[[1]][1:2], collapse="_"))
  #print(hrs)
  meandf <- data.frame(row.names=rownames(df))
  for (i in 1:(length(df)/rep)){
    st = (i-1)*3 + 1
    ed = st + rep -1
    subdf <- df[,c(st:ed)]
    #coln <- paste0("mean",st)
    coln <- hrs[st]
    meandf[, coln] <- rowMeans(subdf)
  }
  
  mean.scaled <- t(scale(t(as.matrix(meandf))))
  rownames(mean.scaled) <- df$gid
  phmap <- pheatmap(mean.scaled, cluster_cols=FALSE, show_rownames=FALSE, 
           show_colnames=TRUE, cutree_rows = 8, main="top-500 by padj")
  phmap
  return(phmap)
}
hm.sub8a <- PlotHM(out.sub8a)

out.ls8b <- Deseq_hour(sample.sub8b)
#[1] 5869    6
#-- note: fitType='parametric', but the dispersion trend was not well captured by the
#function: y = a/x + b, and a local regression fit was automatically substituted.
#specify fitType='local' or 'mean' to avoid this message next time.
out.sub8b <- out.ls8b[[2]]
write.csv(out.sub8b, "DEG.wt-timecourse-h0h4r2.csv", quote=TRUE, row.names=FALSE)
hm.sub8b <- PlotHM(out.sub8b)


## Case IX : timecourse change ppk by heat
sample.sub9 <- subset(sample.table, sample.table$geno == "ppk")
sample.sub9a <- sample.sub9[-c(10:15),]
sample.sub9b <- sample.sub9[c(1:15),]

out.ls9a <- Deseq_hour(sample.sub9a)
#[1] 8759    6
out.sub9a <- out.ls9a[[2]]
write.csv(out.sub9a, "DEG.ppk-timecourse-h0h6.csv", quote=TRUE, row.names=FALSE)
hm.sub9a <- PlotHM(out.sub9a)

out.ls9b <- Deseq_hour(sample.sub9b)
#[1] 5894    6
out.sub9b <- out.ls9b[[2]]
write.csv(out.sub9b, "DEG.ppk-timecourse-h0h4r2.csv", quote=TRUE, row.names=FALSE)
hm.sub9b <- PlotHM(out.sub9b)


## Case X : combined DEG by time and mutant
sample.sub0a <- sample.table[-c(10:15,31:36),]
sample.sub0b <- sample.table[c(1:15,22:36),]

Deseq_both <- function(sample.sub){
  dds.sub <- dds[,colnames(dds) %in% rownames(sample.sub)]
  dds.sub <- DESeqDataSet(dds.sub, ~geno + hour + geno:hour)
  dds.sub <- DESeq(dds.sub)
  res.sub <- results(dds.sub)
  sig.sub <- subset(res.sub, abs(res.sub$log2FoldChange) >= 1 & res.sub$padj < 0.01)
  print(dim(sig.sub))
  
  out.sub <- as.data.frame(sig.sub[,c(1,2,6)])
  out.sub <- merge(out.sub, tpm.df[, rownames(sample.sub)], by=0)
  colnames(out.sub)[1:4] <- c("gid","MeanEXP", "log2FC", "FDR")
  out.sub <- out.sub[order(out.sub$log2FC, decreasing=TRUE), ]
  out.sub$annt <- annt.table[toupper(out.sub$gid), ]$annt
  
  out.ls <- list(res.sub, out.sub)
  return(out.ls)
}
out.ls0a <- Deseq_both(sample.sub0a)
#[1] 1518    6
out.sub0a <- out.ls0a[[2]]
write.csv(out.sub0a, "DEG.combined-h0h6.csv", quote=TRUE, row.names=FALSE)
hm.sub0a <- PlotHM(out.sub0a)

out.ls0b <- Deseq_both(sample.sub0b)
#[1] 2494    6
#-- note: fitType='parametric', but the dispersion trend was not well captured by the
#function: y = a/x + b, and a local regression fit was automatically substituted.
#specify fitType='local' or 'mean' to avoid this message next time.
out.sub0b <- out.ls0b[[2]]
write.csv(out.sub0b, "DEG.combined-h0h4r2.csv", quote=TRUE, row.names=FALSE)
hm.sub0b <- PlotHM(out.sub0b)








  
 


sessionInfo()
#R version 4.3.1 (2023-06-16)
#Platform: x86_64-apple-darwin20 (64-bit)
#Running under: macOS Monterey 12.7.2

#Matrix products: default
#BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
#LAPACK: /Library/Frameworks/R.framework/Versions/4.3-x86_64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0

#locale:
#  [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

#time zone: Asia/Tokyo
#tzcode source: internal

#attached base packages:
#  [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#[1] GO.db_3.18.0                pathviewr_1.1.7             KEGGREST_1.42.0            
#[4] GOstats_2.68.0              Category_2.68.0             Matrix_1.6-4               
#[7] GSEABase_1.64.0             graph_1.80.0                annotate_1.80.0            
#[10] XML_3.99-0.16               org.At.tair.db_3.18.0       AnnotationForge_1.44.0     
#[13] AnnotationDbi_1.64.1        pheatmap_1.0.12             lubridate_1.9.3            
#[16] forcats_1.0.0               stringr_1.5.1               dplyr_1.1.4                
#[19] purrr_1.0.2                 readr_2.1.4                 tidyr_1.3.0                
#[22] tibble_3.2.1                ggplot2_3.4.4               tidyverse_2.0.0            
#[25] DESeq2_1.42.0               SummarizedExperiment_1.32.0 Biobase_2.62.0             
#[28] MatrixGenerics_1.14.0       matrixStats_1.2.0           GenomicRanges_1.54.1       
#[31] GenomeInfoDb_1.38.2         IRanges_2.36.0              S4Vectors_0.40.2           
#[34] BiocGenerics_0.48.1         tximport_1.30.0  
