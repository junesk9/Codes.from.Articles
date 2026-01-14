####### 25.Jan20th
###### MEGA wheat RNAseq
###### June-Sik Kim


library(DESeq2)
library(tximport)
library(pheatmap)

#for GO analysis
library(GSEABase)
library(GOstats)
library(KEGGREST) #for KEGG API
library(clusterProfiler)
library(enrichplot)
library(ggplot2)

setwd("~/理化学研究所　セルロース生産研究チーム Dropbox/June-Sik Kim/Collaborations/25JAN_wheat.RNA_MEGA")
set.seed(101)
key.palette <- rev(hcl.colors(100, "RdYlBu"))
#################################################
############## Load Data
#(tximport)
subpath = "2_salmon"
salmon.files = file.path(list.files(subpath), 'quant.sf')
sample.name <- sapply(salmon.files, function(x){strsplit(x, "/")[[1]][1]})
salmon.files <- sapply(salmon.files, function(x){paste0(subpath,"/",x)})
names(salmon.files) <- sample.name

#(sample.meta)
meta <- data.frame(sample = sample.name)
meta$geno <- sapply(meta$sample, function(x){strsplit(x, "-")[[1]][1]})
meta$tx <- sapply(meta$sample, function(x){strsplit(x, "-")[[1]][2]})
meta$rpt <- sapply(meta$sample, function(x){strsplit(x, "-")[[1]][3]})
rownames(meta) <- meta$sample
meta$tx <- factor(meta$tx, levels=c("Mock","ABA"))
meta$geno <- factor(meta$geno, levels=c("C","WS"))

tx.exp <- tximport(salmon.files, type = "salmon", txOut = TRUE)
TXNAME = rownames(tx.exp$counts)
GENEID <- sapply(TXNAME, function(x)strsplit(x, "[.]")[[1]][1])
GENEID <- as.vector(GENEID)
tx2gene <- data.frame(TXNAME, GENEID)

gn.exp <- summarizeToGene(tx.exp, tx2gene, countsFromAbundance = "no")
tpm.df <- as.data.frame(gn.exp$counts)
cnt.df <- as.data.frame(gn.exp$abundance)
write.csv(tpm.df, file = "wheat-Mega.2501.salmon-SAF.TPM.csv", row.names = TRUE, quote=FALSE)
write.csv(cnt.df, file = "wheat-Mega.2501.salmon-SAF.count.csv", row.names = TRUE, quote=FALSE)

#(repeat-average)
SampleAvg <- function(df){
  meta$sid <- sapply(colnames(df), function(x){paste(strsplit(x, "-")[[1]][1:2], collapse = '-')})
  sample.ids <- unique(meta$sid)
  
  mean.df <- data.frame(row.names=rownames(df))
  for (s in sample.ids){
    meta.sub <- subset(meta, meta$sid == s)
    sample.sub <- meta.sub$sample
    df.sub <- df[, colnames(df) %in% sample.sub]
    colname <- s
    mean.df[, colname] <- rowMeans(df.sub, na.rm=TRUE)
  }
  return(mean.df)
}
tpm_mean.df <- SampleAvg(tpm.df)[,c("C-Mock","WS-Mock","C-ABA","WS-ABA")]
write.csv(tpm_mean.df, file = "wheat-Mega.2501.salmon-SAF.meanTPM.csv", row.names = TRUE, quote=FALSE)

#(GENE ANNT; https://urgi.versailles.inra.fr/download/iwgsc/IWGSC_RefSeq_Annotations/v2.1/)
info.table <- read.csv("0_genome.info/iwgsc_refseqv2.1_functional_annotation.csv", header=T)
annt.table <- subset(info.table, info.table$f.type == "Human readable description")
annt.table <- annt.table[,c(1,4)]
annt.table <- annt.table[!duplicated(annt.table$g2.identifier), ]
rownames(annt.table) <- annt.table[,1]

#(KEGG dataset; convert Entrez ID to Gene ID)
kegg <- keggLink("taes", "pathway") #[1] 45658
kegg.id <- gsub("path:", "", names(kegg))
id.uniq <- unique(kegg.id) #[1] 153
kegg.db <- vector("list", length(id.uniq))
names(kegg.db) <- id.uniq
entrez <- read.table("0_genome.info/wheat-entrez.Ensembl60.txt", header=T, sep="\t") #[1] 1253335       3
entrez <- na.omit(entrez[,c(1,2)]) #[1] 1150501       2
entrez <- entrez[!duplicated(entrez),]
for(i in 1:length(kegg.db)) {
  eids <- as.character(gsub("taes:", "", kegg[kegg.id == id.uniq[i]]))
  for(id in eids){
    if(id %in% entrez[,2]){
      gids <- subset(entrez, id == entrez[,2])[,1]
      kegg.db[[i]] <- c(kegg.db[[i]], gids)
    }
  }
  kegg.db[[i]] <- unique(kegg.db[[i]])
  print(c(i, id.uniq[i], length(eids), length(kegg.db[[i]])))
}
kegg.df <- data.frame(KEGGID = rep(names(kegg.db), sapply(kegg.db, length)), GIDv1 = unlist(kegg.db)) #[1] 25044     2
write.csv(kegg.df, "./0_genome.info/wheat-entrez.KEGGREST.csv", row.names = F)


#(GO dataset)
go <- subset(info.table, info.table$f.type == "Gene Ontology")
go$evi = "IEA"
go <- go[,c(4,6,1)]
goFrame <- GOFrame(go, organism="Triticum aestium")
goAllFrame <- GOAllFrame(goFrame)
gsc <- GeneSetCollection(goAllFrame, setType = GOCollection())
all <- unique(as.vector(go[,3])) #[1] 63251
#(Generate a Org.DB for ClusterCompare)
library(AnnotationForge)
GENEID <- unique(GENEID)
fSym <- data.frame(GID=GENEID, SYMBOL=GENEID, GENENAME=GENEID)
chrs <- sapply(GENEID, function(x){substr(x, 8, 9)})
fChr <- data.frame(GID=GENEID, CHROMOSOME=chrs)
fGO <- data.frame(GID=go[,3], GO=go[,1], EVIDENCE=go[,2])

makeOrgPackage(gene_info=fSym, chromosome=fChr, go=fGO,
               version="0.1",
               maintainer="JS KIM <so@someplace.org>",
               author="JS KIM <so@someplace.org>",
               outputDir = ".",
               tax_id="4565",
               genus="Triticum",
               species="aestivum",
               goTable="go")
## The folder names "org.Taestivum.eg.db" generated 
install.packages("./org.Taestivum.eg.db", repos=NULL, type="source")

#################################################
############## Inter-sample relationship
##cut-off by TPM>1
rmeans <- rowMeans(tpm.df, na.rm=TRUE)
rnames <- names(rmeans[rmeans > 10]) # & 100 > rmeans])
#106382 -> 49041
tpm.df2 <- tpm.df[rnames,]
#[1] 63232    16

### Pearson corr 
p <- cor(tpm.df2, method="pearson")
write.csv(p, "wheat-Mega.2501.Pearson-corr.csv", quote=FALSE)

##Dendrogram
d1 <- as.dist(1-cor(tpm.df2, use="pairwise.complete.obs", method="spearman"))
d2 <-  as.dist(t(1-cor(tpm.df2, use="pairwise.complete.obs", method="spearman")))
c1 <- hclust(d1, method="complete")
c2 <- hclust(d2, method="complete")
plot(c1, main="wheat-Mega.2501-spearman-complete") #dendrogram
#draw heatmap
heatmap(as.matrix(d1,d2),
        Colv=as.dendrogram(c2), Rowv=as.dendrogram(c1),
        scale="none", col=cm.colors(256),
        main="wheat-Mega.2501-spearman-complet", margin=c(6,6))

#PCA
rpca <- prcomp(x=t(tpm.df2), scale=T) 
PropVar <- summary(rpca)$importance[2,]
barplot(PropVar[1:10], ylab="Proportion of Variance", main="wheat-Mega.2501")
write.csv(rpca$x, "wheat-Mega.2501.PCA-whole.csv", quote=FALSE)
#PCA-scatter
#display.brewer.all()
xlab = paste0("PC1: ",round(PropVar[1]*100,1), "%")
ylab = paste0("PC2: ",round(PropVar[2]*100,1), "%")
plot(x=rpca$x[,1],y=rpca$x[,2], xlab=xlab, ylab=ylab,cex=1, pch=16, main="wheat-Mega.2501.PCA-TPM>10")
sample.ids <- sapply(colnames(tpm.df), function(x){paste(strsplit(x, "-")[[1]][1:2], collapse = '-')})
text(x=rpca$x[,1],y=rpca$x[,2],sample.ids, pos=1)

#################################################
############## DEG by DESeq
dds <- DESeqDataSetFromTximport(gn.exp, meta, ~tx + geno + tx:geno)
#(1. global DEGs)
runDESeq <- function(dds){
  dds <- DESeq(dds)
  res <- results(dds)
  sig <- subset(res, abs(res$log2FoldChange) >= 1 & res$padj < 0.05)
  print(dim(sig))
  out <- as.data.frame(sig[,c(1,2,6)])
  
  out <- merge(out, tpm_mean.df[rownames(out),], by=0)
  colnames(out)[1:4] <- c("gid","MeanEXP", "log2FC", "FDR")
  out <- out[order(out$log2FC, decreasing=TRUE), ]
  
  out_ls <- list(res, out)
  return(out_ls)
}
out1 <- runDESeq(dds) #[1] 4406    8
write.csv(out1[[2]], "wheat-Mega.2501.DEG01.geno-tx.csv", quote=TRUE, row.names=FALSE)
PlotHM500 <- function(df){
  df <- df[order(df$FDR),]
  rownames(df) <- df$gid
  meandf <- df[c(1:500),c("C-Mock","WS-Mock","C-ABA","WS-ABA")] #top500 FDR
  
  mean.scaled <- t(scale(t(as.matrix(meandf))))
  phmap <- pheatmap(mean.scaled, cluster_cols=FALSE, show_rownames=FALSE, 
                    show_colnames=TRUE, main="top-500 by FDR, AvgTPM, z-scaled")
  phmap
  return(phmap)
}
hm1 <- PlotHM500(out1[[2]])

#(2. WS vs. C; mock)
meta.sub <- subset(meta, meta$tx %in% c("Mock"))
dds2 <- dds[,colnames(dds) %in% meta.sub$sample]
dds2 <- DESeqDataSet(dds2, ~geno)
out2 <- runDESeq(dds2) #[1] 6031    6
out2[[2]]$annt <- annt.table[out2[[2]][,"gid"], 2]
write.csv(out2[[2]], "wheat-Mega.2501.DEG02.geno-mock.csv", quote=TRUE, row.names=FALSE)
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
  
  title = paste(c("up:",up.len," dn:",dn.len), collapse="")
  plot(v3$Log2FC,v3$qval, pch=20, cex=0.4, col="grey70", ylim=c(0,max_q), 
       xlim=c(-max_fc,max_fc), xlab="Log2FC", ylab="-log10FDR", main=title)
  points(v1$Log2FC,v1$qval, pch=20, cex=0.5, col="#ED0422")
  points(v2$Log2FC,v2$qval, pch=20, cex=0.5, col="#3498DB")
  abline(v=0, lty=2)
}
PlotVC(out2[[1]])
hm2 <- PlotHM500(out2[[2]])

#(3. C by ABA)
meta.sub <- subset(meta, meta$geno %in% c("C"))
dds3 <- dds[,colnames(dds) %in% meta.sub$sample]
dds3 <- DESeqDataSet(dds3, ~tx)
out3 <- runDESeq(dds3) #[1] 9630    6
out3[[2]]$annt <- annt.table[out3[[2]][,"gid"], 2]
write.csv(out3[[2]], "wheat-Mega.2501.DEG03.C-ABA.csv", quote=TRUE, row.names=FALSE)
PlotVC(out3[[1]])
hm3 <- PlotHM500(out3[[2]])

#(4. WS vs. C; ABA)
meta.sub <- subset(meta, meta$tx %in% c("ABA"))
dds4 <- dds[,colnames(dds) %in% meta.sub$sample]
dds4 <- DESeqDataSet(dds4, ~geno)
out4 <- runDESeq(dds4) #[1] 3626    6
out4[[2]]$annt <- annt.table[out4[[2]][,"gid"], 2]
write.csv(out4[[2]], "wheat-Mega.2501.DEG04.geno-aba.csv", quote=TRUE, row.names=FALSE)
PlotVC(out4[[1]])
hm4 <- PlotHM500(out4[[2]])

#(5. WS by ABA)
meta.sub <- subset(meta, meta$geno %in% c("WS"))
dds5 <- dds[,colnames(dds) %in% meta.sub$sample]
dds5 <- DESeqDataSet(dds5, ~tx)
out5 <- runDESeq(dds5) #[1] 7241    6
out5[[2]]$annt <- annt.table[out5[[2]][,"gid"], 2]
write.csv(out5[[2]], "wheat-Mega.2501.DEG05.WS-ABA.csv", quote=TRUE, row.names=FALSE)
PlotVC(out5[[1]])
hm5 <- PlotHM500(out5[[2]])

#(merge 2 and 3)
up2 <- out2[[2]][out2[[2]]$log2FC > 0, "gid"]
up3 <- out3[[2]][out3[[2]]$log2FC > 0, "gid"]
up23 <- unique(c(up2, up3)) #[1] 7598 from [1] 9052

dn2 <- out2[[2]][out2[[2]]$log2FC < 0, "gid"]
dn3 <- out3[[2]][out3[[2]]$log2FC < 0, "gid"]
dn23 <- unique(c(dn2, dn3)) #[1] 4710 from [1] 5817

up23.tpm <- tpm_mean.df[up23, ]
up23.tpm$FDR <- runif(n=length(up23), min=1, max=20) #random FDR for [PlotHM500]
up23.tpm <- up23.tpm[rowMeans(up23.tpm[,c(1:4)]) > 200, ] #[1] 1864    5; select by rowMean
up.p <- PlotHM500(up23.tpm)

dn23.tpm <- tpm_mean.df[dn23, ]
dn23.tpm$FDR <- runif(n=length(dn23), min=1, max=20) #random FDR for [PlotHM500]
#dn23.tpm <- dn23.tpm[rowMeans(dn23.tpm[,c(1:4)]) > 200, ] #[1] 883   5; select by rowMean
dn.p <- PlotHM500(dn23.tpm)

#################################################
############## GO
## Conventional GSEA analysis
runGO2 <- function(gids){
  out_ls = list()
  for (f in c("BP","CC","MF")){
    p <- GSEAGOHyperGParams(
      name = "Paramaters",
      geneSetCollection = gsc,
      geneIds = gids,
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
  
  return(out_ls)
} 

##simple merge of EnrichResearts to compareCluster dotplot
en2 <- enricher(up2, TERM2GENE=go[,c(1,3)])
en3 <- enricher(up3, TERM2GENE=go[,c(1,3)])
eml <- list(up2 = en2, up3 = en3)
emc <- merge_result(eml)
dotplot(emc, group=TRUE, showCategory = 3, label_format = 30, color="p.adjust") #Visualize by dot plot

##compareCluster, Based on a custom Org.DB generation.
library(org.Taestivum.eg.db)
deg2 <- read.csv("wheat-Mega.2501.DEG02.geno-mock.csv", header=T)
deg3 <- read.csv("wheat-Mega.2501.DEG03.C-ABA.csv", header=T)
deg4 <- read.csv("wheat-Mega.2501.DEG04.geno-aba.csv", header=T)
deg5 <- read.csv("wheat-Mega.2501.DEG05.WS-ABA.csv", header=T)

up2 <- deg2[deg2$log2FC > 1 , "gid"] #1332 #[1] 2951
up3 <- deg3[deg3$log2FC > 1 , "gid"]  #3449 #[1] 6101
up4 <- deg4[deg4$log2FC > 1 , "gid"]  #408  #[1] 1087
up5 <- deg5[deg5$log2FC > 1 , "gid"]  #1997 #[1] 4340

dn2 <- deg2[deg2$log2FC <= -1 , "gid"]  #[1] 2880
dn3 <- deg3[deg3$log2FC <= -1 , "gid"]  #[1] 2937
dn4 <- deg4[deg4$log2FC <= -1 , "gid"]  #[1] 2539
dn5 <- deg5[deg5$log2FC <= -1 , "gid"]  #[1] 2901

deg_ls = list(WSmock.Cmock = up2, WSaba.Caba = up4 , Caba.Cmock = up3, WSaba.WSmock = up5)
GOBP <- compareCluster(geneClusters=deg_ls, fun="enrichGO", ont="BP", OrgDb = "org.Taestivum.eg.db", keyType = "GID", pvalueCutoff=1)
GOMF <- compareCluster(geneClusters=deg_ls, fun="enrichGO", ont="MF", OrgDb = "org.Taestivum.eg.db", keyType = "GID", pvalueCutoff=1)
GOCC <- compareCluster(geneClusters=deg_ls, fun="enrichGO", ont="CC", OrgDb = "org.Taestivum.eg.db", keyType = "GID", pvalueCutoff=1)
dotplot(GOBP, group=TRUE, showCategory = 4, label_format = 40, color="p.adjust")
treeplot(pairwise_termsim(GOBP), showCategory = 3, geneClusterPanel = "dotplot", color="p.adjust")
write.csv(as.data.frame(GOBP), "wheat-Mega.2501.upGOBP-2.csv", row.names = FALSE)
write.csv(as.data.frame(GOCC), "wheat-Mega.2501.upGOCC-2.csv", row.names = FALSE)
write.csv(as.data.frame(GOMF), "wheat-Mega.2501.upGOMF-2.csv", row.names = FALSE)

deg_ls = list(WSmock.Cmock = dn2, WSaba.Caba = dn4 , Caba.Cmock = dn3, WSaba.WSmock = dn5)
GOBP <- compareCluster(geneClusters=deg_ls, fun="enrichGO", ont="BP", OrgDb = "org.Taestivum.eg.db", keyType = "GID", pvalueCutoff=1, qvalueCutoff=0.5)
GOMF <- compareCluster(geneClusters=deg_ls, fun="enrichGO", ont="MF", OrgDb = "org.Taestivum.eg.db", keyType = "GID", pvalueCutoff=1, qvalueCutoff=0.5)
GOCC <- compareCluster(geneClusters=deg_ls, fun="enrichGO", ont="CC", OrgDb = "org.Taestivum.eg.db", keyType = "GID", pvalueCutoff=1, qvalueCutoff=0.5)
dotplot(GOBP, group=TRUE, showCategory = 4, label_format = 40, color="p.adjust")
write.csv(as.data.frame(GOBP), "wheat-Mega.2501.dnGOBP-2.csv", row.names = FALSE)
write.csv(as.data.frame(GOCC), "wheat-Mega.2501.dnGOCC-2.csv", row.names = FALSE)
write.csv(as.data.frame(GOMF), "wheat-Mega.2501.dnGOMF-2.csv", row.names = FALSE)

#################################################
############## KEGG
kegg.df <- read.csv("0_genome.info/wheat-entrez.KEGGREST.csv", header=T)
g2g3 <- read.csv("0_genome.info/02G--_03G_GeneID.csv", header=T)

G3ENTZ <- function(g3ids){
  g2ids <- unique(g2g3[g2g3[,2] %in% g3ids, "v1.1"])
  enids <- unique(entrez[entrez[,1] %in% g2ids, 2])
  return(enids)
}
up2e <- G3ENTZ(up2) #[1] 933
up3e <- G3ENTZ(up3) #[1] 1985
up4e <- G3ENTZ(up4) #[1] 350
up5e <- G3ENTZ(up5) #[1] 1500
deg_ls = list(WSmock.Cmock = up2e, WSaba.Caba = up4e , Caba.Cmock = up3e, WSaba.WSmock = up5e)
KEGG <- compareCluster(geneClusters=deg_ls, fun = "enrichKEGG", organism = "taes", pvalueCutoff=1)
KEGG@compareClusterResult$Description <- sapply(KEGG@compareClusterResult$Description, function(x){strsplit(x, " -")[[1]][1]})
dotplot(KEGG, group=TRUE, showCategory = 4, label_format = 40, color="p.adjust")
write.csv(as.data.frame(KEGG), "wheat-Mega.2501.upKEGG-2.csv", row.names = FALSE)

dn2e <- G3ENTZ(dn2) #[1] 930
dn3e <- G3ENTZ(dn3) #[1] 922
dn4e <- G3ENTZ(dn4) #[1] 795
dn5e <- G3ENTZ(dn5) #[1] 975
deg_ls = list(WSmock.Cmock = dn2e, WSaba.Caba = dn4e , Caba.Cmock = dn3e, WSaba.WSmock = dn5e)
KEGG <- compareCluster(geneClusters=deg_ls, fun = "enrichKEGG", organism = "taes", pvalueCutoff=1)
KEGG@compareClusterResult$Description <- sapply(KEGG@compareClusterResult$Description, function(x){strsplit(x, " -")[[1]][1]})
dotplot(KEGG, group=TRUE, showCategory = 4, label_format = 40, color="p.adjust")
write.csv(as.data.frame(KEGG), "wheat-Mega.2501.dnKEGG-2.csv", row.names = FALSE)

sessionInfo()
# R version 4.5.2 (2025-10-31)
# Platform: aarch64-apple-darwin20
# Running under: macOS Sequoia 15.7.1
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
#   [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] GOstats_2.76.0              Category_2.76.0             Matrix_1.7-4                GSEABase_1.72.0             graph_1.88.1               
# [6] annotate_1.88.0             XML_3.99-0.20               AnnotationDbi_1.72.0        BiocManager_1.30.27         pheatmap_1.0.13            
# [11] tximport_1.38.2             DESeq2_1.50.2               SummarizedExperiment_1.40.0 Biobase_2.70.0              MatrixGenerics_1.22.0      
# [16] matrixStats_1.5.0           GenomicRanges_1.62.1        Seqinfo_1.0.0               IRanges_2.44.0              S4Vectors_0.48.0           
# [21] BiocGenerics_0.56.0         generics_0.1.4             
# 
# loaded via a namespace (and not attached):
#   [1] RColorBrewer_1.1-3      rstudioapi_0.17.1       jsonlite_2.0.0          tidydr_0.0.6            magrittr_2.0.4          ggtangle_0.0.9         
# [7] farver_2.1.2            fs_1.6.6                vctrs_0.6.5             memoise_2.0.1           RCurl_1.98-1.17         ggtree_4.0.3           
# [13] htmltools_0.5.9         S4Arrays_1.10.1         SparseArray_1.10.8      gridGraphics_0.5-1      htmlwidgets_1.6.4       plyr_1.8.9             
# [19] cachem_1.1.0            igraph_2.2.1            lifecycle_1.0.4         pkgconfig_2.0.3         gson_0.1.0              R6_2.6.1               
# [25] fastmap_1.2.0           digest_0.6.39           aplot_0.2.9             ggnewscale_0.5.2        patchwork_1.3.2         RSQLite_2.4.5          
# [31] polyclip_1.10-7         httr_1.4.7              abind_1.4-8             compiler_4.5.2          withr_3.0.2             bit64_4.6.0-1          
# [37] fontquiver_0.2.1        S7_0.2.1                BiocParallel_1.44.0     DBI_1.2.3               ggforce_0.5.0           R.utils_2.13.0         
# [43] MASS_7.3-65             rappdirs_0.3.3          DelayedArray_0.36.0     tools_4.5.2             ape_5.8-1               scatterpie_0.2.6       
# [49] R.oo_1.27.1             glue_1.8.0              nlme_3.1-168            GOSemSim_2.36.0         grid_4.5.2              cluster_2.1.8.1        
# [55] reshape2_1.4.5          fgsea_1.36.0            gtable_0.3.6            R.methodsS3_1.8.2       tidyr_1.3.2             data.table_1.18.0      
# [61] XVector_0.50.0          ggrepel_0.9.6           pillar_1.11.1           stringr_1.6.0           yulab.utils_0.2.3       genefilter_1.92.0      
# [67] splines_4.5.2           tweenr_2.0.3            dplyr_1.1.4             treeio_1.34.0           lattice_0.22-7          survival_3.8-3         
# [73] bit_4.6.0               tidyselect_1.2.1        RBGL_1.86.0             fontLiberation_0.1.0    GO.db_3.22.0            locfit_1.5-9.12        
# [79] Biostrings_2.78.0       fontBitstreamVera_0.1.1 stringi_1.8.7           lazyeval_0.2.2          ggfun_0.2.0             codetools_0.2-20       
# [85] gdtools_0.4.4           tibble_3.3.0            qvalue_2.42.0           Rgraphviz_2.54.0        ggplotify_0.1.3         cli_3.6.5              
# [91] xtable_1.8-4            systemfonts_1.3.1       Rcpp_1.1.0              png_0.1-8               parallel_4.5.2          ggplot2_4.0.1          
# [97] blob_1.2.4              DOSE_4.4.0              AnnotationForge_1.52.0  bitops_1.0-9            tidytree_0.4.6          ggiraph_0.9.2          
# [103] scales_1.4.0            purrr_1.2.0             crayon_1.5.3            rlang_1.1.6             cowplot_1.2.0           fastmatch_1.1-6        
# [109] KEGGREST_1.50.0  
