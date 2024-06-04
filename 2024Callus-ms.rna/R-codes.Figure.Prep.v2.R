#########################################
########### Figure preparation R code log 
########### June-Sik Kim
########### 2024.05.18


################################################################################
################################### Read processing & Mapping on Linux Bash
################################################################################

#(Trimmomatic read QC)
# $ for i in *R1_001.fastq.gz; do fi=`echo $i| cut -d S -f 1`; echo $fi; java -jar trimmomatic-0.39.jar PE -threads 30 $i ${fi}S1_L001_R2_001.fastq.gz ../1_trimmomatic/${fi}1.paired.fq.gz ../1_trimmomatic/${fi}1.unpair.fq.gz ../1_trimmomatic/${fi}2.paired.fq.gz ../1_trimmomatic/${fi}2.unpair.fq.gz ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:36; done

#(STAR indexing)
# $ STAR --runMode genomeGenerate --genomeDir . --genomeFastaFiles Oryza_sativa.IRGSP-1.0.dna.chromosome.fa --sjdbGTFfile Oryza_sativa.IRGSP-1.0.56.chr.gtf --sjdbOverhang 99 --runThreadN 30 --genomeSAindexNbases 13
#(STAR run)
# [for rice & tobacco] $ STAR --runThreadN 15 --twopassMode Basic --genomeDir ./0_ref/01_rice-star --readFilesCommand 'gunzip -c' --readFilesIn %s %s --outSAMtype BAM SortedByCoordinate --outSAMattrRGline ID:FLOWCELLID LB:%s_lib1 PL:illumina SM:%s PU:unit3 --outFileNamePrefix ./2_star/%s --quantMode GeneCounts TranscriptomeSAM" % (fa1, fa2, prefix, prefix, prefix)
# [for bamboo] $ "STAR --runThreadN 8 --twopassMode Basic --genomeDir ./0_ref/03_tobacco-star --readFilesCommand 'gunzip -c' --readFilesIn %s %s --outSAMtype BAM SortedByCoordinate --outSAMattrRGline ID:FLOWCELLID LB:%s_lib1 PL:illumina SM:%s PU:unit3 --outFileNamePrefix ./5_star-Sr/%s --quantMode GeneCounts TranscriptomeSAM --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0.15 --outFilterMatchNmin 0.15" % (fa1, fa2, prefix, prefix, prefix)

#(RSEM indexing)
# $ rsem-prepare-reference --gtf Oryza_sativa.IRGSP-1.0.56.chr.gtf -p 30 Oryza_sativa.IRGSP-1.0.dna.chromosome.fa Oryza_sativa.IRGSP-1.0.RSEM_index
# $ rsem-prepare-reference --gtf Bamboo.Hic.gtf -p 20 Bamboo.HIC.genome.rename.fasta Bamboo.HIC.RSEM_index
# $ rsem-prepare-reference --gtf Nitab-v4.5_gene_models_Chr_Edwards2017.gtf -p 6 Nitab-v4.5_genome_Chr_Edwards2017.fasta Nitab-v4.5_genome_Chr_Edwards2017.RSEM_index
#(RSEM run)
#$ for i in *toTranscriptome.out.bam; do fi=`echo $i| cut -d l -f 1`; echo $fi;  rsem-calculate-expression --alignments --paired-end -p 10 $i ../0_ref/01_rice-star/Oryza_sativa.IRGSP-1.0.RSEM_index  ./${fi}


set.seed(101)
library(dplyr)
library(tximport)
library(beeswarm)
library(stringr)

################################################################################
#################################################### Data parsing & prep
################################################################################
#1. RNAseq (Os)
path.detail = "./01_RSEM-output/01_rsem-Os/"
rsem.files <- file.path(path.detail, list.files(path.detail, ".isoforms.results"))
f <- function(s) strsplit(strsplit(s,"/")[[1]][5], "A")[[1]][1]
sample.name <- sapply(rsem.files, f)
sample.name = as.vector(sample.name)
sample.name <- gsub("-", "_", sample.name)
names(rsem.files) <- sample.name
rsem.files <- rsem.files[1:24]
meta.os <- data.frame(row.names=names(rsem.files))
meta.os$cond <- sapply(names(rsem.files), function(x){strsplit(x, "_Os_")[[1]][1]})
meta.os$was <- sapply(names(rsem.files), function(x){strsplit(x, "_Os_|wk")[[1]][2]})
meta.os$rpt <- sapply(names(rsem.files), function(x){strsplit(x, "_Os_|wk|0d")[[1]][3]})
tx.exp <- tximport(rsem.files, type="rsem", txIn=TRUE, txOut=TRUE)
tx2gene <- read.table(rsem.files[1], header=T)[,c(1:2)]
colnames(tx2gene) <- c("TXNAME","GENEID")
tx2gene <- tx2gene %>% filter(stringr::str_starts(TXNAME, 'Os')) #Prot-coding genes only "OsXXXXX"
tx.exp <- lapply(tx.exp, function(x) if(is.matrix(x)) return (x[tx2gene$TXNAME,]) else return(x))
gn.exp <- summarizeToGene(tx.exp, tx2gene, countsFromAbundance = "scaledTPM")
tpm.os <- as.data.frame(gn.exp$abundance)

#2. RNAseq (Pb)
path.detail = "./01_RSEM-output/02_rsem-Pb/"
rsem.files <- file.path(path.detail, list.files(path.detail, ".isoforms.results"))
f <- function(s) strsplit(strsplit(s,"/")[[1]][5], "A")[[1]][1]
sample.name <- sapply(rsem.files, f)
sample.name = as.vector(sample.name)
sample.name <- gsub("-", "_", sample.name)
names(rsem.files) <- sample.name
rsem.files <- rsem.files[1:24]
meta.pb <- data.frame(row.names=names(rsem.files))
meta.pb$cond <- sapply(names(rsem.files), function(x){strsplit(x, "_Pb_")[[1]][1]})
meta.pb$was <- sapply(names(rsem.files), function(x){strsplit(x, "_Pb_|wk")[[1]][2]})
meta.pb$rpt <- sapply(names(rsem.files), function(x){strsplit(x, "_Pb_|wk|0d")[[1]][3]})
tx.exp <- tximport(rsem.files, type="rsem", txIn=TRUE, txOut=TRUE)
tx2gene <- read.table(rsem.files[1], header=T)[,c(1:2)]
colnames(tx2gene) <- c("TXNAME","GENEID")
gn.exp <- summarizeToGene(tx.exp, tx2gene, countsFromAbundance = "scaledTPM")
tpm.pb <- as.data.frame(gn.exp$abundance)

#3. RNAseq (Pn)
path.detail = "./01_RSEM-output/03_rsem-Pn/"
rsem.files <- file.path(path.detail, list.files(path.detail, ".isoforms.results"))
f <- function(s) strsplit(strsplit(s,"/")[[1]][5], "A")[[1]][1]
sample.name <- as.vector(sapply(rsem.files, f))
sample.name <- gsub("-", "_", sample.name)
sample.name <- as.vector(sapply(sample.name, function(x){strsplit(x, ".", fixed=TRUE)[[1]][1]}))
names(rsem.files) <- sample.name
rsem.files <- rsem.files[1:24]
meta.pn <- data.frame(row.names=names(rsem.files))
meta.pn$cond <- sapply(names(rsem.files), function(x){strsplit(x, "_Pn_")[[1]][1]})
meta.pn$was <- sapply(names(rsem.files), function(x){strsplit(x, "_Pn_|wk")[[1]][2]})
meta.pn$rpt <- sapply(names(rsem.files), function(x){strsplit(x, "_Pn_|wk|0d")[[1]][3]})
tx.exp <- tximport(rsem.files, type="rsem", txIn=TRUE, txOut=TRUE)
tx2gene <- read.table(rsem.files[1], header=T)[,c(1:2)]
colnames(tx2gene) <- c("TXNAME","GENEID")
gn.exp <- summarizeToGene(tx.exp, tx2gene, countsFromAbundance = "scaledTPM")
tpm.pn <- as.data.frame(gn.exp$abundance)

#4. RNAseq (SR1)
path.detail = "./01_RSEM-output/04_rsem-sr1/"
rsem.files <- file.path(path.detail, list.files(path.detail, ".isoforms.results"))
f <- function(s) strsplit(strsplit(s,"/")[[1]][5], "[.]")[[1]][1]
sample.name <- sapply(rsem.files, f)
sample.name = as.vector(sample.name)
sample.name <- gsub("-", "_", sample.name)
names(rsem.files) <- sample.name
rsem.files <- rsem.files[1:24]
meta.sr <- data.frame(row.names=names(rsem.files))
meta.sr$cond <- sapply(names(rsem.files), function(x){strsplit(x, "_sr1_")[[1]][1]})
meta.sr$was <- sapply(names(rsem.files), function(x){strsplit(x, "_sr1_|wk")[[1]][2]})
meta.sr$rpt <- sapply(names(rsem.files), function(x){strsplit(x, "_sr1_|wk|0d")[[1]][3]})
tx.exp <- tximport(rsem.files, type="rsem", txIn=TRUE, txOut=TRUE)
tx2gene <- read.table(rsem.files[1], header=T)[,c(1:2)]
colnames(tx2gene) <- c("TXNAME","GENEID")
gn.exp <- summarizeToGene(tx.exp, tx2gene, countsFromAbundance = "scaledTPM")
tpm.sr <- as.data.frame(gn.exp$abundance)

#5. Metabolome + Hormonome
ms <- read.csv("meta-hormone.parsed.csv", header=T, row.names=1) #[1] 108 476
#               CellLine tx time_wk plate_position    tZ   tZR tZRPs     cZ    cZR  cZRPs DZ   DZR DZRPs    iP    iPR  iPRPs tZ7G  tZ9G  tZOG    cZOG
#hd_Os_0d2            Os hd       0          P1-C2    NA 0.502 0.481 19.102 27.896 15.081 NA 0.981    NA 0.113  8.632 41.147   NA 3.332    NA 639.713
#hd_Os_0d3            Os hd       0          P1-D2    NA 0.811 1.037 12.882 29.439 16.252 NA 0.863 0.160 0.140 14.030 82.806   NA 3.917    NA 612.047
#hd_Os_0d4            Os hd       0          P1-B2    NA 0.696 0.559 21.147 28.049 11.909 NA 1.186 0.173 0.135  8.965 44.610   NA 3.510    NA 649.635
#sg_Os_s0.1_0d1       Os sg       0          P2-C3 0.205 0.958 0.943  8.823 31.853 21.265 NA    NA 0.036 0.142  8.574 59.112   NA 4.004 2.745 681.734
#sg_Os_s0.1_0d2       Os sg       0          P2-D3 0.191 1.116 0.949  9.050 36.338 20.059 NA    NA 0.131    NA  9.165 45.827   NA 3.584 2.353 722.472
#sg_Os_s0.1_0d3       Os sg       0          P2-E3 0.246 1.239 1.172  9.608 40.070 22.249 NA    NA    NA 0.124  9.826 53.378   NA 4.440 2.353 760.821
#               tZROG   cZROG tZRPsOG cZRPsOG  DZ9G iP7G   iP9G GA1 GA4    ABA SA      JA     IAA    IAAsp  JAIle
#hd_Os_0d2      6.007 668.643   1.879 229.213 1.832   NA 50.920  NA  NA 16.418 NA      NA 154.699 11053.73 40.731
#hd_Os_0d3      7.832 709.124   2.929 309.634 1.647   NA 47.087  NA  NA 18.287 NA      NA 187.086 12051.08 31.701
#hd_Os_0d4      7.448 532.734   3.144 315.883 2.061   NA 59.225  NA  NA 15.089 NA      NA 204.185 14018.26 37.804
#sg_Os_s0.1_0d1 5.013 543.447   1.737 190.612 1.751   NA 66.145  NA  NA 13.260 NA 225.194 150.098 23761.97 50.017
#sg_Os_s0.1_0d2 5.839 523.478   2.001 198.713 1.450   NA 33.569  NA  NA 11.648 NA 687.392      NA 24637.42 66.632
#sg_Os_s0.1_0d3 6.370 592.179   2.136 218.846 1.448   NA 39.454  NA  NA 12.562 NA 715.642 141.317 25546.05 66.581

ms <- subset(ms, ms$tx %in% c("h","d"))
ms.os <- subset(ms, ms$CellLine == "Os")[,-c(1:3)]
ms.pn <- subset(ms, ms$CellLine == "Pn")[,-c(1:3)]
ms.pb <- subset(ms, ms$CellLine == "Pb")[,-c(1:3)]
ms.sr <- subset(ms, ms$CellLine == "sr1")[,-c(1:3)]

################################################################################
#################################################### Figure drawing
################################################################################
# Figure 1C
hum <- read.table(pipe("pbpaste"), header=T) #[1] 29  6
#     C1    C2    C3    C4    C5    C6    D1    D2    D3    D4    D5    D6
# 1 100.0 100.0 100.0 100.0 100.0 100.0 100.0 100.0 100.0 100.0 100.0 100.0
# 2 100.0 100.0  98.8  98.9  99.7  99.7  98.2  98.5  98.4  98.6  98.2  98.5
# 3  99.3 100.0  98.5  98.6  99.7  99.4  95.7  96.8  96.5  96.6  96.0  96.8
# 4  99.3  99.7  98.5  98.6  99.4  99.4  94.0  95.3  94.6  94.8  94.1  95.1
# 5  99.3  99.7  98.5  98.2  99.4  99.1  91.8  93.8  93.0  93.1  91.5  93.6
# 6  99.3  99.7  98.2  98.2  99.4  98.8  89.7  92.1  91.1  91.7  89.7  92.4
hum.m <- data.frame(PS=rowMeans(hum[,c(1:6)]), SS=rowMeans(hum[,c(7:12)]))
hum.m$ttest <- NA
for(r in c(1:dim(hum)[1])){
  vc <- as.vector(t(t(hum[r,])))
  p <- pairwise.t.test(x=vc, g=rep(c("C","D"), each=6), alternative ="two.sided", p.adjust.method ="none" ,paired=FALSE)$p.value
  hum.m[r, "ttest"] <- p
}

ylim=c(min(hum, na.rm=TRUE), max(hum, na.rm=TRUE)) #[1]  46.8 100.0
palette <- rep(c("#00B7EB80", "#D6006E80"), each=6)
matplot(hum, type="p", pch=19, cex=.5, col=palette, xaxt="n",
        xlab="weeks", ylab="rel. Humidity", main="in-plate waterloss by different sealing")
axis(side=1, at=c(0,7,14,21,28), labels=c(0,1,2,3,4))
matplot(hum.m[,c(1:2)], type="l", lty=1, lwd=3, col=c("#00B7EB","#D6006E"), add=TRUE)

# Figure 1H
wc <- read.table(pipe("pbpaste"), header=T) 
wc.sr <- df[,c(1:3,13:15)] #SR1
wc.os <- df[,c(4:6,16:18)] #Os
wc.pn <- df[,c(7:9,19:21)] #Pn
wc.pb <- df[,c(10:12,22:24)] #Pb

mplot <- function(df, species){
  dff.m <- data.frame(row.names=rownames(dff), ps=rowMeans(dff[,c(1:3)]), ss=rowMeans(dff[,c(4:6)]))
  main = paste0("Callus water content, ", species)
  matplot(dff, type="p", ylim=c(90,100),pch=20, cex=1, col=rep(c("#00B7EB","#D6006E"), each=3),
          main=main, xlab="week", ylab="rel. content (%)")
  matplot(dff.m, type="l", lty=1, lwd=3,  col=c("#00B7EB","#D6006E"), add=TRUE)
}
mplot(wc.sr,"SR1")
mplot(wc.os,"Os")
mplot(wc.pn,"Pn")
mplot(wc.pb,"Pb")

ttest <- function(df){
  for(i in c(1:dim(dff)[1])){
    p <- t.test(x=dff[i,c(1:3)], y=dff[i,c(4:6)], var.equal=F, paired=F)$p.value
    outline = paste0("line",i,": ",p)
    print(outline)
  }
}
ttest(wc.sr)
ttest(wc.os)
ttest(wc.pn)
ttest(wc.pb)

# Figure 2-left
tpm2pcc <- function(tpm.df){
  df.means <- rowMeans(tpm.df, na.rm=TRUE)
  df.names <- names(df.means[df.means > 1])
  tpm.df2 <- tpm.df[df.names, ]
  
  pcc.df <- cor(tpm.df2, method="pearson")
  pcc.vec = c()
  for (idx in seq(1, dim(pcc.df)[1], by=3)){
    rpt.ls <- rownames(pcc.df)[idx:(idx+2)]
    p.rpt <- pcc.df[rpt.ls, rpt.ls]
    p.vec <- c(p.rpt[1,2], p.rpt[1,3], p.rpt[2,3])
    pcc.vec <- c(pcc.vec, p.vec)
  }
  return(pcc.vec)
}
os.pcc <- tpm2pcc(tpm.os)
sr.pcc <- tpm2pcc(tpm.sr)
pn.pcc <- tpm2pcc(tpm.pn)
pb.pcc <- tpm2pcc(tpm.pb)
pcc.df <- data.frame(sr.pcc, os.pcc, pb.pcc, pn.pcc)
beeswarm(pcc.df, ylim=c(0.75,1), col=rainbow(4), pch=19, cex=0.8,corral="random", main="RNAseq tpm>1")
boxplot(pcc.df, ylim=c(0.9,1), col="#00000000", outline=F, add=TRUE)

# Figure 2-right
ms.sub <- subset(ms, ms$tx %in% c("h","d"))[,-c(1:3)]
ms.sub <- na.omit(t(ms.sub)) # [1] 453  96
ms.sub <- ms.sub[which(apply(ms.sub, 1, var) != 0), ] #[1] 371  96
p <- cor(ms.sub, method="pearson")
p.parse = list()
for (sp in unique(ms$CellLine)){
  p.sub <- p[rownames(ms[ms$CellLine %in% sp, ]), rownames(ms[ms$CellLine %in% sp, ])]
  pcc.vec = c()
  for (idx in seq(1, dim(p.sub)[1], by=3)){
    rpt.ls <- rownames(p.sub)[idx:(idx+2)]
    p.rpt <- p.sub[rpt.ls, rpt.ls]
    p.vec <- c(p.rpt[1,2], p.rpt[1,3], p.rpt[2,3])
    pcc.vec <- c(pcc.vec, p.vec)
  }
  p.parse <- c(p.parse, list(pcc.vec))
}
pcc.df <- data.frame(sr=p.parse[[4]], os=p.parse[[1]], pn=p.parse[[3]], pb=p.parse[[2]])
palette4 <- c("#F8766D", "#7CAE00" ,"#00BFC4" ,"#C77CFF")
beeswarm(pcc.df, ylim=c(0.75,1), col=palette4, pch=19, cex=0.8, corral="random", main="metabolome+hormonome")
boxplot(pcc.df, ylim=c(0.75,1), col="#00000000", outline=F, add=TRUE)

# Figure 3A
pca <- function(tpm.df){
  df.means <- rowMeans(tpm.df, na.rm=TRUE)
  df.names <- names(df.means[df.means > 1])
  tpm.df2 <- tpm.df[df.names, ]
  rpca <- prcomp(x=t(tpm.df2), scale=T) 
  PropVar <- summary(rpca)$importance[2,]
  
  out <- c(list(rpca), list(PropVar))
  return(out)
}
pca.ms <- pca(ms.sub)
rpca <- pca.ms[[1]]
PropVar <- pca.ms[[2]]
xlab = paste0("PC1: ",round(PropVar[1]*100,1), "%")
ylab = paste0("PC2: ",round(PropVar[2]*100,1), "%")
#(up)
plot(x=rpca$x[,1],y=rpca$x[,2], xlab=xlab, ylab=ylab, pch=16, col="grey100", main="MS data all-4", cex=0.1)
points(rpca$x[rownames(ms[ms$CellLine %in% "sr1", ]), c(1:2)], pch=19, cex=1, col="#F8766D")
points(rpca$x[rownames(ms[ms$CellLine %in% "Os", ]), c(1:2)], pch=19, cex=1, col="#7CAE00")
points(rpca$x[rownames(ms[ms$CellLine %in% "Pn", ]), c(1:2)], pch=19, cex=1, col="#00BFC4")
points(rpca$x[rownames(ms[ms$CellLine %in% "Pb", ]), c(1:2)], pch=19, cex=1, col="#C77CFF")
#(middle)
plot(x=rpca$x[,1],y=rpca$x[,2], xlab=xlab, ylab=ylab, pch=16, col="grey100", main="MS data all-4", cex=0.1)
points(rpca$x[rownames(ms[ms$tx %in% "h", ]), c(1:2)], pch=19, cex=1, col="#00B7EB")
points(rpca$x[rownames(ms[ms$tx %in% "d", ]), c(1:2)], pch=19, cex=1, col="#D6006E")
#(bottom)
plot(x=rpca$x[,1],y=rpca$x[,2], xlab=xlab, ylab=ylab, pch=16, col="grey100", main="MS data all-4", cex=0.1)
points(rpca$x[rownames(ms[ms$time_wk %in% "1", ]), c(1:2)], pch=19, cex=1, col="#C6622A")
points(rpca$x[rownames(ms[ms$time_wk %in% "2", ]), c(1:2)], pch=19, cex=1, col="#9C3E40")
points(rpca$x[rownames(ms[ms$time_wk %in% "3", ]), c(1:2)], pch=19, cex=1, col="#6A1F39")
points(rpca$x[rownames(ms[ms$time_wk %in% "3", ]), c(1:2)], pch=19, cex=1, col="#1D0B14")

# Figure 3B
pp.common <- rownames(tpm.pn)[rownames(tpm.pn) %in% rownames(tpm.pb)]
tpm.pp <- cbind(tpm.pn[pp.common,], tpm.pb[pp.common, ])
pca.pp <- pca(tpm.pp)
rpca <- pca.pp[[1]]
PropVar <- pca.pp[[2]]
xlab = paste0("PC1: ",round(PropVar[1]*100,1), "%")
ylab = paste0("PC2: ",round(PropVar[2]*100,1), "%")
#(up)
plot(x=rpca$x[,1],y=rpca$x[,2], xlab=xlab, ylab=ylab, pch=16, col="grey100", main="Pn+Pb; RNAseq", cex=0.1)
points(rpca$x[rownames(ms[ms$CellLine %in% "sr1", ]), c(1:2)], pch=19, cex=1, col="#F8766D")
points(rpca$x[rownames(ms[ms$CellLine %in% "Os", ]), c(1:2)], pch=19, cex=1, col="#7CAE00")
points(rpca$x[rownames(ms[ms$CellLine %in% "Pn", ]), c(1:2)], pch=19, cex=1, col="#00BFC4")
points(rpca$x[rownames(ms[ms$CellLine %in% "Pb", ]), c(1:2)], pch=19, cex=1, col="#C77CFF")
#(middle)
plot(x=rpca$x[,1],y=rpca$x[,2], xlab=xlab, ylab=ylab, pch=16, col="grey100", main="Pn+Pb; RNAseq", cex=0.1)
points(rpca$x[rownames(ms[ms$tx %in% "h", ]), c(1:2)], pch=19, cex=1, col="#00B7EB")
points(rpca$x[rownames(ms[ms$tx %in% "d", ]), c(1:2)], pch=19, cex=1, col="#D6006E")
#(bottom)
plot(x=rpca$x[,1],y=rpca$x[,2], xlab=xlab, ylab=ylab, pch=16, col="grey100", main="Pn+Pb; RNAseq", cex=0.1)
points(rpca$x[rownames(ms[ms$time_wk %in% "1", ]), c(1:2)], pch=19, cex=1, col="#C6622A")
points(rpca$x[rownames(ms[ms$time_wk %in% "2", ]), c(1:2)], pch=19, cex=1, col="#9C3E40")
points(rpca$x[rownames(ms[ms$time_wk %in% "3", ]), c(1:2)], pch=19, cex=1, col="#6A1F39")
points(rpca$x[rownames(ms[ms$time_wk %in% "3", ]), c(1:2)], pch=19, cex=1, col="#1D0B14")






sessionInfo()
# R version 4.3.2 (2023-10-31)
# Platform: aarch64-apple-darwin20 (64-bit)
# Running under: macOS Sonoma 14.4.1
# 
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
# LAPACK: /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0
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
#   [1] stringr_1.5.1   beeswarm_0.4.0  tximport_1.28.0 dplyr_1.1.4    

