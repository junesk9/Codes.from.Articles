################### SCODE attempts
################### 2023.02.21
################### KIM JS
##### Original code from "https://github.com/hmatsu1226/SCODE"

library(MASS)
library(stringr) #to modulate the row-names
library(reshape2) # for melt()

#### input files
setwd("~/Dropbox (理化学研究所　セルロース生産研究チーム)/オオムギ圃場mRNA.ChIPseq/14_SCODE")
set.seed(NULL)# do not set seed 
rpm_f = "../12_Ito.san-Seurat/vst-dynamic-genes.2842.count.csv"
tfls_f = "../00_IPSR_KIBR_Sampledata/TranscriptionFactor_list_blastp_tophit_Gene_ID_rep_AtOsBd.txt"
meta_f = "../12_Ito.san-Seurat/whole.metadata.csv"

#### Load files & set options
rpm_m = read.csv(rpm_f, header=TRUE, row.names = 1)
colnames(rpm_m) <- str_sub(colnames(rpm_m), 2) ##remove the "X" from colnames
#[1] 2842 1940
tfls = read.table(tfls_f,sep="\t", header=T)[,2]
#[1] 2567
meta_m <- read.csv(meta_f, header=TRUE, row.names=1)
#[1] 1940   24
lin4_m <- meta_m[c("Date.4", "Lineage4")]
lin2_m <- meta_m[c("Date.4", "Lineage2")]
lin4_m <- na.omit(lin4_m)
#[1] 1195    2
lin2_m <- na.omit(lin2_m)
#[1] 1031    2

gene_n = 100
cell_n = 100
z = 4
iter_n = 100

maxB <- 2.0
minB <- -10.0
###########################################
###################### main body
####1. case of dynamic TFs only (lineage 4)
X.total <- as.matrix(subset(rpm_m, rownames(rpm_m) %in% tfls)) ## gene x cells matrix

# z-sampling function needed.
sample_Z <- function(){
  for (i in 1:z){
    for(j in 1:cell_n){
      Z[i,j] <<- exp(new_B[i]*pstime[j]) + runif(1, min=-0.001, max=0.001)
    }
  }
}

# data_prep
lin_m <- lin4_m
X <- X.total[,rownames(lin_m)]
gene_n <- dim(X)[1]
cell_n <- dim(X)[2]
W <- matrix(rep(0, gene_n * z), nrow=gene_n, ncol=z)
Z <- matrix(rep(0, cell_n * z), nrow=z, ncol=cell_n)
WZ <- matrix(nrow=gene_n, ncol=cell_n)
pstime <- lin_m[,2]/max(lin_m[,2])


# variables initiations
RSS <- Inf
new_B <- rep(0, z)
old_B <- rep(0, z)
for (i in 1:z){
  new_B[i] <- runif(1, min=minB, max=maxB)
  old_B[i] <- new_B[i]
}

# W and B optimization
for (ite in 1:iter_n){
  #sampling B
  target  <- floor(runif(1, min=1, max=z+1))
  new_B[target] <- runif(1, min=minB, max=maxB)
  
  #last calc.
  if(ite == iter_n){
    for(i in 1:z){
      new_B[i] <- old_B[i]
    }
  }
  
  #Z sampling
  sample_Z()
  
  #regression
  for(i in 1:gene_n){
    X.lm <- lm(X[i,] ~ t(Z)-1)
    for(j in 1:z){
      W[i,j] <- X.lm$coefficients[j]
    }
    WZ[i,] <- W[i,] %*% Z
  }
  
  #RSS
  tmp_RSS <- sum((X-WZ)**2)
  if (tmp_RSS < RSS ){
    RSS <- tmp_RSS
    old_B[target] <- new_B[target]
  } else {
    new_B[target] <- old_B[target]
  }
}

# inter A
B <- matrix(rep(0, z * z), nrow=z, ncol=z)
for (i in 1:z){
  B[i,i] <- new_B[i]
}
invW <- ginv(W) #require.MASS()
A <- W %*% B %*% invW



###########################################
###################### main body
####2. find-out optimized RSS & corr(A) for a given z-value
rpt_n = 100
lin_m <- lin4_m
X <- X.total[,rownames(lin_m)]
gene_n <- dim(X)[1]
cell_n <- dim(X)[2]

A_list <- list()
rss_v <- vector()
st_time <- Sys.time() ## for time-measure
for(rpt in 1:rpt_n){
  W <- matrix(rep(0, gene_n * z), nrow=gene_n, ncol=z)
  Z <- matrix(rep(0, cell_n * z), nrow=z, ncol=cell_n)
  WZ <- matrix(nrow=gene_n, ncol=cell_n)
  pstime <- lin_m[,2]/max(lin_m[,2]) #normalize 0<x<1
  RSS <- Inf
  new_B <- rep(0, z)
  old_B <- rep(0, z)
  
  for (ite in 1:iter_n){
    #sampling B
    target  <- floor(runif(1, min=1, max=z+1))
    new_B[target] <- runif(1, min=minB, max=maxB)
    
    #last calc.
    if(ite == iter_n){
      for(i in 1:z){
        new_B[i] <- old_B[i]
      }
    }
    
    #Z sampling
    sample_Z()
    
    #regression
    for(i in 1:gene_n){
      X.lm <- lm(X[i,] ~ t(Z)-1)
      for(j in 1:z){
        W[i,j] <- X.lm$coefficients[j]
      }
      WZ[i,] <- W[i,] %*% Z
    }
    
    #RSS
    tmp_RSS <- sum((X-WZ)**2)
    if (tmp_RSS < RSS ){
      RSS <- tmp_RSS
      old_B[target] <- new_B[target]
    } else {
      new_B[target] <- old_B[target]
    }
  }
  
  rss_v <- c(rss_v, RSS)
  
  B <- matrix(rep(0, z * z), nrow=z, ncol=z)
  for (i in 1:z){
    B[i,i] <- new_B[i]
  }
  invW <- ginv(W) #require.MASS()
  A <- W %*% B %*% invW
  A_end = length(A_list) + 1
  A_list[[A_end]] <- A
  if(rpt%%10 == 0) {
    ed_time <- Sys.time()
    print(rpt)
    print(ed_time - st_time)
    }
}

###3. select minimum RSS values & A matrices, then calculate cor()
rss_top <- order(rss_v)[1:10]
cor_df <- data.frame(matrix(nrow=length(c(A_list[[1]]))))
for (i in rss_top){
  minA <- A_list[[i]]
  col_id = paste0("rss",i)
  cor_df[,col_id] = c(minA)
}
cor_df = cor_df[,-1]
c <- cor(cor_df, use="complete.obs", method="pearson")
c[lower.tri(c)] <- NA ## mark replicated corr values
c <- na.omit(melt(c)) ## remove the marked values
c <- subset(c, c$Var1 != c$Var2) ## remove the self-corr
c <- c[order(c$value, decreasing=TRUE),] ## sorting

#output files
prefix = paste0("z",paste0(z,paste0("-rpt",paste0(rpt_n,paste0("-itr",iter_n)))))
rss_f = paste0(prefix, ".RSS.csv")
cor_f = paste0(prefix, ".A-cor.csv")

write.csv(as.data.frame(rss_v), rss_f, row.names = F)
write.csv(c, cor_f, row.names = F)
dir.create(prefix)
for (i in 1:rpt_n){
  str_n <- str_pad(i, 4, pad = "0") ## 0 to "0001"
  A_f = paste0(prefix,paste0(".rpt",paste0(str_n,".csv")))
  A <- A_list[[i]]
  colnames(A) <- rownames(X)
  write.csv(A, paste(prefix, A_f, sep="/"), row.names=rownames(X))
}




sessionInfo()
#R version 4.2.2 (2022-10-31)
#Platform: aarch64-apple-darwin20 (64-bit)
#Running under: macOS Monterey 12.6.3

#Matrix products: default
#LAPACK: /Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/lib/libRlapack.dylib

#locale:
#  [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

#attached base packages:#
#  [1] stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#  [1] MASS_7.3-58.2
