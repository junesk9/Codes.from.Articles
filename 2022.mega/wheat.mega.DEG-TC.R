####Sleuth analysis of wheat RNA-seq timecourse (Mega)
####2022.01.26
####KIM JS

library(sleuth)
library(splines)
library(cclust)
library(gplots)

rm(list = ls(all.names = TRUE)) #願掛け
s2c <- read.table("wheat.mega1.txt", header=TRUE, stringsAsFactors=FALSE)
head(s2c)

############################################################################
####### PCA Plots
so <- sleuth_prep(s2c)
ppv <- plot_pc_variance(so, use_filtered = TRUE, pca_number = 10)
PCs <- ppv$data
PCs
#PC_count        var
#1        1 74.4367444
#2        2 13.2246400
#3        3  8.0627519
#4        4  1.0730896
#5        5  0.6444067

plot_pca(so, use_filtered=TRUE, point_alpha=1, point_size=2)

###########################################################################
####### Compare multi-conditions (Time + alpha)
smpl <- s2c[c(7:41),] #L8 time + cond
so <- sleuth_prep(smpl)
so <- sleuth_fit(so, ~day + osmo + day:osmo, 'full')
so <- sleuth_fit(so, ~day + osmo, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')
results_table_lrt <- sleuth_results(so, 'reduced:full', test_type = 'lrt')
table(results_table_lrt[,"qval"] < 0.05)

#####generate output table
#DEG data
sig_ids <- results_table_lrt$target_id[which(results_table_lrt$qval < 0.05)]
sig_table <- results_table_lrt[results_table_lrt$target_id %in% sig_ids,][1:3]
head(sig_table)
#target_id         pval         qval
#1 TRAESCS5B02G406300.1 1.287014e-14 7.511271e-10
#2 TRAESCS3A02G014800.1 4.958676e-13 1.273013e-08
#3 TRAESCS5D02G199500.1 6.543708e-13 1.273013e-08
#4 TRAESCS5B02G191700.1 1.170078e-12 1.707202e-08
#Extract TPM data
tpm_matrix <- sleuth_to_matrix(so, 'obs_raw', 'tpm')
d <- tpm_matrix #for sure
d <- cbind(rownames(d), data.frame(d, row.names=NULL))
colnames(d)[1] <- "target_id"
tpm <- d[d$target_id %in% sig_ids,]
dim(tpm)
#[1] 3732   36

#Merging & outputs
out <- merge(tpm, sig_table, by="target_id", sort=T)
dim(out)
#[1] 3732   38
write.csv(out, file="DEG-timecourse.csv", row.names=FALSE)

#####################################################################
################################### Clustering by FC
#Need to express the data to fold-change against 0-hr point
#re-call the whole data & convert TPM to log2FC
so <- sleuth_prep(s2c[c(1:41),])
so_tpm <- sleuth_to_matrix(so, 'obs_raw', 'tpm')
so_tpm <- cbind(rownames(so_tpm), data.frame(so_tpm, row.names=NULL))
colnames(so_tpm)[1] <- "target_id"
deg_tpm <- so_tpm[so_tpm$target_id %in% sig_ids,]

means <- as.vector(rowMeans(deg_tpm[,c(2:7)]))
deg_fc <- deg_tpm[,-c(1:7)]
deg_fc <- log((deg_fc+1e-12)/(means+1e-12),2)  ##1e-12 to evade divide by zero error
rownames(deg_fc) <- deg_tpm$target_id #Give rowname and gene_id
d_matrix <- as.matrix(deg_fc)
d_df <- as.data.frame(d_matrix)

##########[Optional, taking time] gplot
# Set color range
colors = c(seq(-5,-0.5,length=101),seq(-0.5,0.5,length=101),seq(0.5,5,length=101))
colors = unique(colors) ## break vector should composite unique numbers
length(colors)
my_palette <- colorRampPalette(c("green", "black", "red"))(n = 300)

heatmap.2(d_matrix, col=my_palette, breaks=colors, density.info="none", trace="none", dendrogram=c("row"), symm=F,symkey=F,symbreaks=T, scale="none", Colv=FALSE)

#########cclustering by k=2 to k=12
for(i in 2:12){
  kmer <- cclust(d_matrix, i, 100, method="kmeans")
  d_df[[paste("FC_k",i,sep="")]] <- kmer$cluster
}
#Writedown the output
write.csv(d_df, file="cclust.timecourse.csv", row.names=TRUE)





#############################################################
####### Spline-based timecourse (define too dynamic genes in the same condition)
smpl <- s2c[c(42:65),] ## N, DR [failed to be modeled]
smpl <- s2c[c(1:6,24:41),] ##L8, water [failed to be modeled]
smpl <- s2c[c(1:23),] ### L8, DR [failed to be modeled]

Time <- smpl$day
Time <- as.numeric(Time) ##願掛け
full_design <- model.matrix(formula(~ ns(Time, df = 3))) ## Adjust df as (No.-2)
tail(full_design)
#20           1         0.2831131         0.8382211
#21           1         0.2831131         0.8382211
#22           1         0.2831131         0.8382211
#23           1         0.2831131         0.8382211

####### Identify DEGs
so <- sleuth_prep(smpl, full_model = full_design)
so <- sleuth_fit(so)
so <- sleuth_fit(so, formula = ~ 1, fit_name = "reduced")
so <- sleuth_lrt(so, "reduced", "full")
results_table_lrt <- sleuth_results(so, 'reduced:full', test_type = 'lrt')
table(results_table_lrt[,"qval"] < 0.05)

##################generate output table
#DEG data
sig_ids <- results_table_lrt$target_id[which(results_table_lrt$qval < 0.05)]
sig_table <- results_table_lrt[results_table_lrt$target_id %in% sig_ids,][1:3]
head(sig_table)
#target_id         pval         qval
#1 TRAESCS3B02G353000.1 1.205915e-15 4.782266e-11
#2 TRAESCS4A02G278900.1 3.211681e-15 4.782266e-11
#3 TRAESCS1A02G325400.1 5.051114e-15 4.782266e-11
#4 TRAESCS3A02G014800.1 1.772913e-15 4.782266e-11
#5 TRAESCS1A02G188700.1 5.627236e-15 4.782266e-11
#Extract TPM data
tpm_matrix <- sleuth_to_matrix(so, 'obs_raw', 'tpm')
d <- tpm_matrix #for sure
d <- cbind(rownames(d), data.frame(d, row.names=NULL))
colnames(d)[1] <- "target_id"
tpm <- d[d$target_id %in% sig_ids,]

###Generate output file
out <- merge(tpm, sig_table, by="target_id", sort=T)
dim(out)
#[1] 292  44
write.csv(out, file="DEG.timecourse.csv", row.names=FALSE)

#####################################################################
################################### Clustering by FC
#Generate FC
deg_tpm <- tpm
means <- as.vector(rowMeans(deg_tpm[,c(2:7)]))
deg_fc <- tpm[,-c(1:7)]
deg_fc <- log((deg_fc+1e-12)/(means+1e-12),2)  ##1e-12 to evade divide by zero error
rownames(deg_fc) <- deg_tpm$target_id #Give rowname and gene_id
d_matrix <- as.matrix(deg_fc)
d_df <- as.data.frame(d_matrix)

##########[Optional, taking time] gplot
# Set color range
colors = c(seq(-5,-0.5,length=101),seq(-0.5,0.5,length=101),seq(0.5,5,length=101))
colors = unique(colors) ## break vector should composite unique numbers
length(colors)
my_palette <- colorRampPalette(c("green", "black", "red"))(n = 300)

#heatmap.2(d_matrix, col=my_palette, breaks=colors, density.info="none", trace="none", dendrogram=c("row"), symm=F,symkey=F,symbreaks=T, scale="none", Colv=FALSE)

#########cclustering by k=2 to k=12
for(i in 2:12){
  kmer <- cclust(d_matrix, i, 100, method="kmeans")
  d_df[[paste("FC_k",i,sep="")]] <- kmer$cluster
}
#Writedown the output
write.csv(d_df, file="cclust.timecourse.csv", row.names=TRUE)




