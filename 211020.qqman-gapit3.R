##### Manhattan plot by qqman()
##### From output csv from GAPIT3
##### 2021.10.20 Junesk9

library(qqman)
setwd("~/Library/CloudStorage/Box-Box/Per_080451_June-Sik Kim/16. Barley.SV274.GWAS/01_gwas.gapit3/Fig5_novel")

####Source Code editing
#trace(manhattan, edit=TRUE)

f_ls = list.files(pattern = "GWAS.Results.csv$")
for(f in f_ls){
  prefix = unlist(strsplit(f, split="\\."))
  prefix = prefix[c(2:4)]
  prefix = paste(prefix,collapse=".")
  #print(prefix)
  
  out_f = paste("../MHT.",prefix,".pdf",sep="")

  #d <- read.table("GAPIT.Blink.X08_head.time.GWAS.Results.csv",header=T,sep=",")
  d <- read.table(f,header=T,sep=",")
  d_qq <- d[,c(1,2,3,4,9)]
  colnames(d_qq) <- c("SNP","CHR","BP","P","FDR")
  d_qq$CHR <- as.numeric(gsub("H","",d_qq$CHR))
  head(d_qq)
  #SNP CHR        BP            P          FDR
  #1 SNP_248177   2 166783484 1.776253e-10 5.619088e-05
  #2 SNP_211456   2  29556053 1.994305e-07 3.154443e-02
  #3 SNP_530590   3 634188674 8.227466e-07 8.675726e-02
  #4 SNP_371706   2 759829153 9.145902e-06 2.631706e-01
  #5 SNP_436706   3  39455388 1.015668e-05 2.631706e-01
  #6 SNP_897697   5 645846970 1.065716e-05 2.631706e-01
  
  ###### Graph statistics (significant SNPs (FDR<0.05); threshold line)
  sig <- d_qq[d_qq$FDR < 0.05,] 
  th = -log10(0.05/length(d_qq$P)) ## Bonferroni-corrected 0.05 [or, so-called Bonferroni threshold] same to the green solid line of GAPIT3 output 

  ######Lightening the data
  library(dplyr)
  d_slim = d_qq
  d_slim$BP <- round(d_slim$BP, digits=-6)
  d_slim$P <- 10^-(round(-log10(d_slim$P), digits=2))
  d_slim <- d_slim %>% distinct(d_slim$BP, d_slim$P, .keep_all = TRUE)
  ####n= 316345 -> 66379
  
  ###set graph height ylim
  ylim = max(-log10(d_slim$P))
  ylim = 10* ceiling(ylim/10)

  #####################drawing plot
  pdf(out_f, w=5, h=2, pointsize=4, useDingbats=F)
  manhattan(d_slim, suggestiveline=FALSE, genomewideline=th, col=c("grey10","grey70"),main=prefix, cex=2, highlight=sig$SNP, ylim = c(0,ylim), xaxt='n', tck=0.02, xlab="", ylab="")
  ### Suggestiveline = 1e-5   genomewideline = 5e-8 (https://jojoshin.hatenablog.com/entry/2016/01/18/142622)
  ### [tck=0.02] to put label ticks inside
  dev.off()
}


