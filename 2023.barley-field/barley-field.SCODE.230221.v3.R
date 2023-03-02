################### SCODE attempts
################### 2023.02.21
################### KIM JS
##### Original code from "https://github.com/hmatsu1226/SCODE"
##### Run by $Rscript SCODE.R [z-score]


library(MASS)
library(stringr) #to modulate the row-names
library(reshape2) # for melt()

#### input files
set.seed(NULL)# do not set seed 
rpm_f = "rpm.csv"
tfls_f = "../TranscriptionFactor_list_blastp_tophit_Gene_ID_rep_AtOsBd.txt"
meta_f = "../whole.metadata.csv"

### Arguments parsing
args <- commandArgs(trailingOnly = T)
z = as.numeric(args[1])

#### Load files & set options
rpm_m = read.csv(rpm_f, header=TRUE, row.names = 1)
#rpm_m = read.table(rpm_f, sep="\t",header=TRUE, row.names = 1)

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

#z = 2 #now input as args[1]

iter_n = 100
rpt_n = 100

maxB <- 2.0
minB <- -10.0

lin_m <- lin4_m
###########################################
###################### main body
###0. Output files
prefix = paste0("z",paste0(z,paste0("-rpt",paste0(rpt_n,paste0("-itr",iter_n)))))
rss_f = paste0(prefix, ".RSS.csv")
cor_f = paste0(prefix, ".A-cor.csv")
dir.create(prefix)
print(prefix)

####1. Prep data
#X.total <- as.matrix(subset(rpm_m, rownames(rpm_m) %in% tfls)) ## gene x cells matrix
X.total <- as.matrix(rpm_m)

# z-sampling function needed.
sample_Z <- function(){
  for (i in 1:z){
    for(j in 1:cell_n){
      Z[i,j] <<- exp(new_B[i]*pstime[j]) + runif(1, min=-0.001, max=0.001)
    }
  }
}

# resume the run if necessary
st = 1
A.files <- file.path(prefix, list.files(prefix, ".A.csv"))
f <- function(s) as.numeric(str_sub(strsplit(s,"[.]")[[1]][2], 4))
A.st <- as.vector(sapply(A.files, f))
if (length(A.st) > 0){
  A.st = max(A.st)
  #st = A.st
}

rss_v <- vector()
rss.file = file.path(rss_f)
if (length(rss.file) > 0 & st > 1){
  rss_v <- read.csv(rss_f, header=T)[,1]
  rss_v <- rss_v[1:st]
}


####2. find-out optimized RSS & corr(A) for a given z-value
X <- X.total[,rownames(lin_m)]
gene_n <- dim(X)[1]
cell_n <- dim(X)[2]

A_list <- list()
st_time <- Sys.time() ## for time-measure
for(rpt in st:rpt_n){
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
  
  ##calc B & A
  B <- matrix(rep(0, z * z), nrow=z, ncol=z)
  for (i in 1:z){
    B[i,i] <- new_B[i]
  }
  invW <- ginv(W) #require.MASS()
  A <- W %*% B %*% invW
  #A_end = length(A_list) + 1
  #A_list[[A_end]] <- A

  ##output A, B-table file
  str_n <- str_pad(rpt, 4, pad = "0") ## 1 to "0001"
  A_f = paste0(prefix,paste0(".rpt",paste0(str_n,".A.csv")))
  B_f = paste0(prefix,paste0(".rpt",paste0(str_n,".B.csv")))
  colnames(A) <- rownames(X)
  write.csv(A, paste(prefix, A_f, sep="/"), row.names=rownames(X))
  write.csv(B, paste(prefix, B_f, sep="/"))
  write.csv(as.data.frame(rss_v), rss_f, row.names = F)
  A = c() # flush memory


  ## Time measure, output log
  if(rpt%%10 == 1){
	  ed_time <- Sys.time()
	  prog = paste0(prefix, paste0("; round-",rpt))
	  print(prog)
	  print(ed_time - st_time)
  	}
}

###3. select the 10-minimum RSS values & A matrices, then calculate cor()
# prepare and verify the rss/A datasets again here
A_list <- list()
A.files <- file.path(prefix, list.files(prefix, ".A.csv"))
A.files <- sort(A.files) # for ensure
rss_v <- read.csv(rss_f, header=T)[,1]

A.files <- A.files[1:length(rss_v)] #for matching the length
for (f in A.files){
  A <- as.matrix(read.csv(f, header=T, row.names=1))
  A_end = length(A_list) + 1
  A_list[[A_end]] <- A
}

# calculate the correlation coefficients
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
write.csv(c, cor_f, row.names = F)

#dir.create(prefix)
#for (i in 1:rpt_n){
#  str_n <- str_pad(i, 4, pad = "0") ## 0 to "0001"
#  A_f = paste0(prefix,paste0(".rpt",paste0(str_n,".csv")))
#  A <- A_list[[i]]
#  colnames(A) <- rownames(X)
#  write.csv(A, paste(prefix, A_f, sep="/"), row.names=rownames(X))
#}

