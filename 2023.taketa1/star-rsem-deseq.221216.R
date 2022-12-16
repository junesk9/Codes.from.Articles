####### STAR-RSEM-Tximport
####### 2022.12.16
####### barley nud mutant (Prof. Takeda, IPSR)

library(tximport)
library(DESeq2)


###### Version II
###############################################
############# Load data by Tximport 
rsem.files <- file.path("./0_RSEM", list.files("./0_RSEM/", ".isoforms.results"))
head(rsem.files)
#[1] "./0_RSEM/1_ATCACGAligned.isoforms.results" "./0_RSEM/3_TTAGGCAligned.isoforms.results" "./0_RSEM/4_ACAGTGAligned.isoforms.results"
#[4] "./0_RSEM/5_GCCAATAligned.isoforms.results" "./0_RSEM/6_CAGATCAligned.isoforms.results"

f <- function(s) strsplit(strsplit(s,"[.]")[[1]][2], "/")[[1]][3]
sample_name = sapply(rsem.files, f)
sample_name = as.vector(sample_name)
names(rsem.files) <- sample_name
head(sample_name)
#[1] "1_ATCACGAligned" "3_TTAGGCAligned" "4_ACAGTGAligned" "5_GCCAATAligned" "6_CAGATCAligned"

### Load exp data
tx.exp <- tximport(rsem.files, type="rsem", txIn=TRUE, txOut=TRUE)
# Subset tximport by gene name
tx.exp2 <- lapply(tx.exp, function(x) if(is.matrix(x)) return (x[grepl("HORVU*", rownames(x)),]) else return(x))
# Prepare gene list (by removing .*)
f3 <- function(s) unlist(strsplit(s, "[.]"))[1]
tx2gene <- data.frame(TXNAME = rownames(tx.exp2$counts))
tx2gene$GENEID = unlist(lapply(tx2gene[,1], f3))
tail(tx2gene)
#TXNAME           GENEID
#251581 HORVU7Hr1G122770.2 HORVU7Hr1G122770
#251582 HORVU7Hr1G122770.3 HORVU7Hr1G122770
#251583 HORVU7Hr1G122800.1 HORVU7Hr1G122800
#251584 HORVU7Hr1G122800.2 HORVU7Hr1G122800
#251585 HORVU7Hr1G122870.1 HORVU7Hr1G122870
#251586 HORVU7Hr1G122880.1 HORVU7Hr1G122880

tx.exp2 <- summarizeToGene(tx.exp2, tx2gene, countsFromAbundance = "scaledTPM")
dim(tx.exp2$abundance)
#[1] 39734     5
write.csv(tx.exp2, file = "gene_tpm.star.csv", row.names = TRUE)
#[optional] tx.exp2 <- tx.exp$abundance[grepl("HORVU*", rownames(tx.exp$abundance)),] ## remove unwanted exp-data

############## DESeq
sample_table <- data.frame(condition = factor(c("Sky","Sky","Sac","Sac","Sac")))
rownames(sample_table) <- sample_name
sample_table
#condition
#1_ATCACGAligned       Sky
#3_TTAGGCAligned       Sky
#4_ACAGTGAligned       Sac
#5_GCCAATAligned       Sac
#6_CAGATCAligned       Sac

dds <- DESeqDataSetFromTximport(tx.exp2, sample_table, ~condition)
dds <- DESeq(dds)
res <- results(dds)
res <- res[order(res$padj),]
head(res)
#baseMean log2FoldChange     lfcSE      stat      pvalue        padj
#<numeric>      <numeric> <numeric> <numeric>   <numeric>   <numeric>
#HORVU.MOREX.r3.1HG0000030  184.3832      0.6945699  0.165974  4.184808 2.85407e-05 0.000336893
#HORVU.MOREX.r3.1HG0000040  590.5260     -0.4506211  0.180694 -2.493841 1.26369e-02 0.050848727
#HORVU.MOREX.r3.1HG0000050  325.2794      0.0381841  0.134551  0.283790 7.76572e-01 0.883735602
#HORVU.MOREX.r3.1HG0000060   80.1740     -0.1644491  0.255371 -0.643963 5.19600e-01 0.702327193
#HORVU.MOREX.r3.1HG0000070 4036.5863     -0.0824834  0.104755 -0.787396 4.31050e-01 0.630409126
#HORVU.MOREX.r3.1HG0000080   41.0374      1.1376504  0.344507  3.302252 9.59117e-04 0.006425060

################## Generate output file
res_sub <- res[,c(2,5,6)]
colnames(res_sub) <- c("Log2FC", "Pvalue", "FDR")
de <- merge(res_sub, tx.exp$abundance, by=0, all=TRUE)
de <- de[order(de$FDR),]
#de$log2FoldChange <- -1 * de$log2FoldChange ##if needed
deg <- subset(de, de$padj < 0.01)  ##DEG threshold as fdr < 0.01
deg2 <- deg[order(deg$log2FoldChange, decreasing=TRUE),]
write.csv(de, file = "gene_deg-tpm.star.csv", row.names = FALSE)


###################################################################################
########### Volcano plot
####data transform
vc <- de[,c(1,2,4)]
vc <- na.omit(vc)
colnames(vc) <- c("target_id","log2FC","qval")

####### Classify DEGs
v1 <- subset(vc, vc$log2FC >= 0 & vc$qval < 0.01)
v2 <- subset(vc, vc$log2FC <= 0 & vc$qval < 0.01)
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

#### Plot
plot(v3$log2FC,v3$qval, pch=20, cex=0.4, col="grey70", ylim=c(0,40), xlim=c(-15,15), xlab="Log2FC", ylab="-log10FDR", main="Sky vs. Sac")
points(v1$log2FC,v1$qval, pch=20, cex=0.5, col="#3498DB")
points(v2$log2FC,v2$qval, pch=20, cex=0.5, col="#ED0422")
abline(v=0, lty=2)
####### Add numbers of DEG on the plot
text(-15,39, up.deg, pos=4, cex=1.4, col="#3498DB")
text(-15,37, dn.deg, pos=4, cex=1.4, col="#ED0422")






sessionInfo()
#R version 4.2.2 (2022-10-31)
#Platform: x86_64-apple-darwin17.0 (64-bit)
#Running under: macOS Monterey 12.6

#Matrix products: default
#LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib

#locale:
#  [1] ja_JP.UTF-8/ja_JP.UTF-8/ja_JP.UTF-8/C/ja_JP.UTF-8/ja_JP.UTF-8

#attached base packages:
#  [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#[1] DESeq2_1.38.1               SummarizedExperiment_1.28.0 MatrixGenerics_1.10.0       matrixStats_0.63.0         
#[5] GenomicRanges_1.50.1        GenomeInfoDb_1.34.4         IRanges_2.32.0              S4Vectors_0.36.1           
#[9] tximport_1.26.0             Biobase_2.58.0              BiocGenerics_0.44.0        
