##2020.02.07
## mRNA-seq C04, eer4-3, taf12
## KALLISTO-SlEUTH
## by June-Sik Kim (https://github.com/junesk9/)

library(sleuth)
s2c <- read.table("meta(pe).txt", header=TRUE, stringsAsFactors=FALSE)

##################### Whole samples QC
so <- sleuth_prep(s2c)
#........................................
#normalizing est_counts
#19332 targets passed the filter
#normalizing tpm
#merging in metadata
#summarizing bootstraps
#...............

#########PCA plot
plot_pca(so, color_by = "genotype")
plot_pca(so, color_by = "tissue")

##############################################################
############ extract TPM & generate a dendrogram
### Call TPM matrix
tpm_matrix <- sleuth_to_matrix(so, 'obs_raw', 'tpm')
d <- tpm_matrix #for sure
d <- cbind(rownames(d), data.frame(d, row.names=NULL))
colnames(d)[1] <- "target_id"

#Writedown the output (whole sample TPM)
write.csv(d, file="sleuth-whole.tpm.csv", row.names=FALSE)

### remove the gene_id column
tpm <- d[,-1]
### hclust
hc <- hclust(dist(t(tpm)))
plot(hc)


#############################################################
########## DEG analysis - pairwise
sub <- s2c[c(9:12,17:20),]
so <- sleuth_prep(sub, ~genotype)
so <- sleuth_fit(so)

##LRT 
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')
results_table_lrt <- sleuth_results(so, 'reduced:full', test_type = 'lrt')
table(results_table_lrt[,"qval"] < 0.05)
#FALSE  TRUE 
#12688  6919 

#Wald-test
so <- sleuth_wt(so, 'genotypewt')
results_table_wt <- sleuth_results(so, 'genotypewt')

################ Prepare output data
#Extract common DEG from two LRT & WT analyses
d5.lrt.sig_ids <- results_table_lrt$target_id[which(results_table_lrt$qval < 0.05)]
d5.wt.sig_ids <- results_table_wt$target_id[which(results_table_wt$qval < 0.05)]
shared_ids <- d5.wt.sig_ids[d5.wt.sig_ids %in% d5.lrt.sig_ids]
shared_results <- results_table_wt[results_table_wt$target_id %in% shared_ids,]
out <- shared_results[,c(1:4)]
out$log2FC <- -log(exp(out$b),2)   #<<- Add Log2FC column
out <- subset(out, select=-c(b))  #<<- Remove b column
out <- subset(out, out$log2FC >= 1 | out$log2FC <= -1) #<<- select by FC
head(out)
#    target_id          pval          qval    log2FC
#1 AT3G10800.1  0.000000e+00  0.000000e+00 -9.655409
#2 AT3G08970.1 1.823565e-195 1.787732e-191  3.238570
#3 AT1G36180.1 4.593123e-185 3.001912e-181  3.609536
#4 AT3G10810.1 1.405622e-181 6.890009e-178 -8.707669
#5 AT3G15980.5 5.358101e-161 2.101126e-157  1.967499
#6 AT1G56580.1 8.002464e-144 2.615072e-140  2.451167
shared_ids <- as.vector(out$target_id)  #<<-prepare new id list
#extract raw TPM data
tpm_matrix <- sleuth_to_matrix(so, 'obs_raw', 'tpm')
d <- tpm_matrix #for sure
d <- cbind(rownames(d), data.frame(d, row.names=NULL))
colnames(d)[1] <- "target_id"
tpm <- d[d$target_id %in% shared_ids,]
#merge two tables 
tmp <- merge(tpm, out, by="target_id", sort=T)

#preprare "whole" output (just repeat the steps)
fdr <- results_table_wt[,c(1:4)]
fdr$log2FC <- -log(exp(fdr$b),2)   #<<- Add Log2FC column
fdr <- subset(fdr, select=-c(b))   #<<- Remove b column
all <- merge(d, fdr, by="target_id", sort=T)
head(all)
#target_id    C04.dn1    C04.dn2    C04.dn3    C04.dn4    C04.up1     wt.dn1
#1 AT1G01010.1  40.136905  33.527067  39.950900  37.921904   4.427740  52.650009
#2 AT1G01020.1  20.095826  20.350781  18.657741  22.280924  17.629998  18.479336
#3 AT1G01030.1   4.580222   2.603089   4.914958   2.215522   4.036948   3.896069
#4 AT1G01040.2  17.800275  14.267472  15.617738  15.185346  13.052227  15.348262
#5 AT1G01050.1 141.000619 147.011658 148.218315 128.649153  86.611783 127.108282
#6 AT1G01060.1 165.038759 144.205237 171.502660 149.630251 371.469047 159.080970
#wt.dn2     wt.dn3     wt.dn4      pval      qval     log2FC
#1  42.288752  44.982875  49.053445 0.1216064 0.5932873 -0.7880612
#2  21.990412  17.505219  18.109567 0.4310590 0.8743559  0.2046483
#3   3.957633   2.619693   2.576573 0.6045723 0.9274695  0.2713206
#4  15.908892  14.576715  16.026024 0.4939813 0.8963837  0.1072924
#5 130.030470 120.836210 148.168821 0.4608092 0.8850009  0.1133518
#6 156.418509 164.466625 167.007182 0.4664104 0.8865614  0.3497851

#write-down
write.csv(tmp, file="sleuth-DEG.csv", row.names=FALSE)
write.csv(all, file="sleuth-whole.csv", row.names=FALSE)

###########################################################
################## Volcano plot
#using "all","fdr" dataset from above
all <- merge(d, fdr, by="target_id", sort=T)
vc <- all[,c("target_id","qval","log2FC")]
####Define DEG (FC>2; FDR<0.01)
vd <- tmp[,c("target_id","qval","log2FC")]
v1 <- subset(vd, vd$log2FC > 0)
v2 <- subset(vd, vd$log2FC < 0)
deg_ids <- c(as.vector(vd$target_id))
`%notin%` <- Negate(`%in%`) #temporary designates "not-in" func
v3 <- vc[vc$target_id %notin% deg_ids,] ## non-DEG

# Counting DEGs
up.deg <- length(v1[,1])
dn.deg <- length(v2[,1])
non.deg <- length(v3[,1])

### Transform v1, v2 data
v1$qval <- -log(v1$qval, 10)
v1$qval <- round(v1$qval, digits=1)
v2$qval <- -log(v2$qval, 10)
v2$qval <- round(v2$qval, digits=1)

#### Transform & Reduce the v3 data
v3$qval <- -log(v3$qval, 10)
v3$qval <- round(v3$qval, digits=1)
v3$log2FC <- round(v3$log2FC, digits=1)
v3d <- v3[,-1]
v3 <- unique(v3d)

## Check the Scatter limits
ylimit = -log(min(vd$qval[vd$qval!=0], na.rm=TRUE), 10) #non-zero minimum
xlimit = max(max(vc$log2FC, na.rm=TRUE), -min(vc$log2FC, na.rm=TRUE), max(v1$log2FC), -min(v2$log2FC))
ylimit = ceiling(ylimit)
xlimit = ceiling(xlimit)
ylimit
xlimit

# Drawing plot with solid limits
plot(v3$log2FC,v3$qval, pch=20, cex=0.4, col="grey70", ylim=c(0,150), xlim=c(-10,10), xlab="Log2FC", ylab="-log10FDR")
points(v1$log2FC,v1$qval, pch=20, cex=0.5, col="#ED0422")
points(v2$log2FC,v2$qval, pch=20, cex=0.5, col="#3498DB")
abline(v=0, lty=2)

# Add numbers of DEG on the plot
text(-10,150, up.deg, pos=4, cex=1.4, col="#ED0422")
text(-10,130, dn.deg, pos=4, cex=1.4, col="#3498DB")
