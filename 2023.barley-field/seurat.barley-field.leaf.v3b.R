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


#########J247 DEGs
jdeg1 <- read.csv("03_TradeSeq-assoc/Tradeseq-gam.LM-assoc.csv", header=T, row.names=1)
jdeg2 <- read.csv("02_cluster-DEG/j247.markers.MAST.csv", header=T, row.names=1)
ftg <- read.csv("00_IPSR_KIBR_Sampledata/FLOR_genes.txt", header=T, row.names=1, sep="\t")

jd1 <- subset(jdeg1, jdeg1$J247.pval < 0.01 & abs(jdeg1$meanLogFC) >= 0.5) #[1] 289   8
jd2 <- subset(jdeg2, abs(jdeg2$avg_log2FC) >= 0.5 & jdeg2$p_val_adj < 0.01) #[1] 2388    5 #Bonferroni-correction
ftg <- ftg[rownames(ftg)%in% rownames(jd1),]
#jdeg2[c("HORVU2Hr1G063800","HORVU0Hr1G003020","HORVU7Hr1G024610"), ]

ftg[rownames(ftg) %in% rownames(jd2), ]
# Closest.homolog    TF
# HORVU0Hr1G003020 OsMADS18(LOC_Os07g41370)/HvFUL3/HvBM3  TRUE
# HORVU2Hr1G063800     BdFUL2(Bradi1g59250)/HvFUL2/HvBM8  TRUE
# HORVU3Hr1G010240                OsRAV9(LOC_Os01g04800)  TRUE
# HORVU4Hr1G008610                                HvPhyA FALSE
ftg[rownames(ftg) %in% rownames(jd1), ]
# Closest.homolog    TF
# HORVU0Hr1G003020 OsMADS18(LOC_Os07g41370)/HvFUL3/HvBM3  TRUE
# HORVU2Hr1G063800     BdFUL2(Bradi1g59250)/HvFUL2/HvBM8  TRUE
# HORVU2Hr1G072750                  BdRCN4(Bradi5g09270) FALSE
table(rownames(jd1) %in% rownames(jd2))
# FALSE  TRUE 
# 66   223 
common <- rownames(jd1)[rownames(jd1) %in% rownames(jd2)]

###############################
##### Figure S8 metadata correlation
################################
library(psych)
pairs.panels(meta[,c(4,12,20,21,23)])

###############################
##### Figure 4 qRT-PCR Barplot with points
##### Users/junesk9/Library/CloudStorage/Box-Box/Int_100100_Bioproductivity\ Informatics\ Research\ Team/オオムギ-CDSクローニング_キム上原/06_transgenic.Barley/Barley.FULox.qRT-PCR.xlsx 
################################
d <- read.table(pipe("pbpaste"), header=T)
# Sample geno     FT1_1
# 1    VC1   vc 0.9096184
# 2    VC1   vc 1.1434025
# 3    VC1   vc 1.2775089
# 4    VC2   vc 1.2863947
# 5    VC2   vc 0.9885140
# 6    VC2   vc 1.2002487
d$Sample <- factor(d$Sample, levels=unique(d$Sample))

anova(aov((d$FT1_1~d$geno)))
# Response: d$FT1_1
# Df  Sum Sq Mean Sq F value    Pr(>F)    
# d$geno     2 1670.59  835.29  82.145 8.924e-10 ***
#   Residuals 18  183.03   10.17  
TukeyHSD(aov((d$FT1_1~d$geno)))

library(dunn.test)
library(tidyverse)

shapiro.test(d[c(1:9),"FT1_1"]) #W = 0.93465, p-value = 0.5269
var.test(x=d[c(1:9),"FT1_1"], y=d[c(10:12),"FT1_1"]) #0.06453134  heterovariance

wilcox.test(x=d[c(1:9),"FT1_1"], y=d[c(10:12),"FT1_1"]) #W = 0, p-value = 0.01623
wilcox.test(x=d[c(1:9),"FT1_1"], y=d[c(13:15),"FT1_1"]) #W = 0, p-value = 0.01623
wilcox.test(x=d[c(1:9),"FT1_1"], y=d[c(16:18),"FT1_1"]) #W = 0, p-value = 0.01623
wilcox.test(x=d[c(1:9),"FT1_1"], y=d[c(19:21),"FT1_1"]) #W = 0, p-value = 0.01623

means <- tapply(d$FT1_1, d$Sample, FUN=mean)
std <- tapply(d$FT1_1, d$Sample, FUN=sd)
stdErr <- function(x) {sd(x, na.rm = TRUE) / sqrt(length(x[!is.na(x)]))}
mse <- tapply(d$FT1_1, d$Sample, FUN=stdErr)


#Draw thr bars
x <- barplot(means, ylim=c(0,1.2*max(d$FT1_1)), col=0, ylab="rel. FT1")
##Add points **jitter() can be ignored with small numbers of spot
points(jitter(rep(x, table(d$Sample)), factor = 0.3), d$FT1_1[order(d$Sample)], pch=19, cex=1)
#Add the error bar 
arrows(x, means-std, x, means+std, code=3, lwd=1, angle=90, length=0.15)

############ Days to Awn emergence
dd <- read.table(pipe("pbpaste"), header=T)
# Sample D2E
# 1      v1 195
# 2      v3 193
# 3 FUL2-10 165
# 4 FUL2-12 156
# 5  FUL3-5 167
# 6  FUL3-8 152



###############################
##### Figure S12 Boxplot of FT genes in clusters
meta2 <- meta[,c("seurat_clusters","Accession")] #[1] 1940    2
rpm2 <- t(all.rpm)
rpm2 <- rpm2[rownames(meta2), colnames(rpm2) %in% rownames(ftg)] #[1] 1940  129
ftg.rpm <- cbind(meta2, rpm2) #[1] 1940  131

boxplot(ftg.rpm$HORVU5Hr1G095630~ftg.rpm$seurat_clusters, outline=F, col=ggcolor(7))



#################################
###### Visualize the promoter motif(hexamer) z-test output
setwd("/Users/junesk9/理化学研究所　セルロース生産研究チーム\ Dropbox/June-Sik\ Kim/オオムギ圃場.publication/z_working.datasource/promoter.motif/")

out1 <- read.table("cl45markers.1000-0nt.promoter.6mer-0mis.435n.1000t.z-test.txt", header=T)
out1 <- out1[order(out1$motif) , ] #order by alphabetically
rownames(out1) <- c(1:dim(out1)[1])
# motif GeneC GeneFC     GeneZ        GeneP      GeneFDR MotifC MotifFC     MotifZ       MotifP     MotifFDR
# 1 AAAAAA   390   1.13  5.426139 2.879305e-08 4.813727e-07   3270    1.49  9.1577606 0.000000e+00 0.000000e+00
# 2 AAAAAC   333   1.11  3.382885 3.586431e-04 2.339176e-03    672    1.02  0.3038354 3.806266e-01 1.000000e+00
# 3 AAAAAG   374   1.23  7.513690 2.875478e-14 1.757904e-12   1254    1.76 14.9505414 0.000000e+00 0.000000e+00
# 4 AAAAAT   379   1.10  4.186648 1.415520e-05 1.223201e-04   1567    1.31  6.1476098 3.932962e-10 6.342288e-09
# 5 AAAACA   290   0.95 -1.600178 9.452204e-01 1.000000e+00    600    0.90 -1.9574468 9.748525e-01 1.000000e+00
# 6 AAAACC   196   0.90 -2.144467 9.840023e-01 1.000000e+00    313    0.89 -1.6085774 9.461456e-01 1.000000e+00
o1ns <- out1[out1$MotifFDR >= 0.05, ] #[1] 3297   11
o1sig <- out1[out1$MotifFDR < 0.05, ] #[1] 799  11

crt <- c("CCGACA","CCGACT","CCGACG","CCGACC","ACCGAC","TCCGAC","GCCGAC","CCCGAC") #CCGAC
gbox <- c("CACGTG","GATAAG")
ltre <- c("CCGAAA")
abre <- c("ACGTTC","ACGTGC")
crt <- rownames(out1[out1$motif %in% crt, ])
ltre <- rownames(out1[out1$motif %in% ltre, ])
abre <- rownames(out1[out1$motif %in% abre, ])
gbox <- rownames(out1[out1$motif %in% gbox, ])

plot(o1ns$MotifZ ~ rownames(o1ns), ylim=c(-9,26), pch=19, cex=0.5, col="grey80", main="hexamer abundance in 2-kb promoter")
#plot(o1ns$MotifZ ~ rownames(o1ns), ylim=c(-7,22), pch=19, cex=0.5, col="grey80", main="hexamer abundance in 1-kb promoter")
points(o1sig$MotifZ ~ rownames(o1sig), pch=19, cex=0.5, col=1)
points(out1[crt,"MotifZ"] ~ crt, pch=19, cex=1, col=4)
points(out1[ltre,"MotifZ"] ~ ltre, pch=19, cex=1, col=5)
points(out1[abre,"MotifZ"] ~ abre, pch=19, cex=1, col=3)
points(out1[gbox,"MotifZ"] ~ gbox, pch=19, cex=1, col=2)
abline(h=0, lty=2)

