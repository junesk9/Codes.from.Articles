##### 2023.02.27
#### Post-analysis of SCODE
#### Kim JS

library(stringr) #to modulate the row-names
library(reshape2) #re-form tables

setwd("./01_dynamicTFs-162-lin4")
set.seed(101)

############################################################################
######################## Collect input data
rss.files <- list.files(".", "RSS.csv")
cor.files <- list.files(".", "A-cor.csv")

##### Parse RSS files
rss.df <- data.frame(matrix(nrow=length(read.csv(rss.files[1],header=T)[,1])))
for (z in rss.files){
  rss <- as.vector(read.csv(z,header=T))
  header <- strsplit(z,"-")[[1]][1]
  rss.df[,header] <- rss
}
rss.df <- rss.df[,-1]
rss.df <- rss.df[,c(2,3,4,5,1)]

##### Parse A-corr files
cor.df <- data.frame(matrix(nrow=length(read.csv(cor.files[1],header=T)[,1])))
for (z in cor.files){
  cor <- as.vector(read.csv(z, header=T)$value)
  header <- strsplit(z,"-")[[1]][1]
  cor.df[,header] <- cor
}
cor.df <- cor.df[,-1]
cor.df <- cor.df[,c(2,3,4,5,1)] #order columns

#################################################################
###################### Visualization
title <- strsplit(getwd(),"/")[[1]]
title <- title[length(title)]

pdf_f <- paste0(title,".ssr-cor.pdf")
pdf(file=paste0("../", pdf_f), width=8, height=6)
par(mfrow = c(1, 2))
boxplot(rss.df, ylab="Squared-Sum of Residue", xlab="Given-z (D)", 
        main=title)
boxplot(cor.df, ylab="correlation of A", xlab="Given-z (D)", 
        main="Top-10 Pearson-corr")
dev.off()
par(mfrow = c(1, 1)) #for sure

####################### Network output; 1)mean of the best-20, 2)cut-off by association value.
path = "z6-rpt200-itr100/" #when the best z as 6
edge.files <- file.path(path, list.files(path, ".A.csv"))
best.idx <- order(rss.df$z6)[1:20] ## 20 minimum RSS values index

for (idx in best.idx){
  header <- paste0("rpt",idx)
  edge.wide <- read.csv(edge.files[best.idx[1]], header=T)
  edge.long <- gather(edge.wide, key="g2", value=header, 2:length(edge.wide))
  if (idx == best.idx[1]){
    edge.table <- edge.long
  } else {
    edge.table[,header] <- edge.long[,3]
  }
}
edge.mean <- edge.table[,c(1,2)]
edge.mean$mean <- rowMeans(edge.table[,c(3:length(edge.table))])
edge.mean <- edge.mean[order(abs(edge.mean$mean), decreasing = TRUE),]
rownames(edge.mean) <- NULL ##reset the row numbers
colnames(edge.mean) <- c("TF","target","mean")
## simply visualize the assoc. distribution
hist(abs(edge.mean[,3]))

edge.top <- edge.mean[c(1:200),] ##top-200 
edge.pos <- subset(edge.top, edge.top$mean > 0)
edge.neg <- subset(edge.top, edge.top$mean < 0)
calc.node.freq <- function(df){
  #node.vec <- c(df[,1], df[,2])
  node.vec <- as.vector(df[,1])
  node.freq <- as.data.frame(table(node.vec))
  node.freq <- node.freq[order(node.freq$Freq, decreasing=TRUE),]
  rownames(node.freq) <- node.freq$node.vec
  return(node.freq)
} 

edge.top.freq <- calc.node.freq(edge.top)
edge.pos.freq <- calc.node.freq(edge.pos)
edge.neg.freq <- calc.node.freq(edge.neg)

edge.top.freq$pos <- edge.pos.freq[rownames(edge.top.freq),]$Freq
edge.top.freq$neg <- edge.neg.freq[rownames(edge.top.freq),]$Freq
edge.top.freq[is.na(edge.top.freq)] <- 0 #convert NA to 0
#       node.vec   Freq pos neg
#HORVU1Hr1G013140   35  17  18
#HORVU1Hr1G072100   20  13   7
#HORVU3Hr1G055260   19   6  13
#HORVU6Hr1G070750   19   8  11
#HORVU4Hr1G054420   16   6  10
#HORVU6Hr1G028150   12   3   9
#HORVU5Hr1G109720   10   7   3

##visualize as a stacked barplot
require(ggplot2)
top.wide <- edge.top.freq[,c(1,3,4)]
top.long <- gather(top.wide, key="np", value="count", 2:length(top.wide))

pdf_f <- paste0(title,".top-nodes.pdf")
pdf(file=paste0("../", pdf_f), width=4.8, height=6)
g <- ggplot(top.long, aes(x=factor(node.vec, level=rev(top.wide$node.vec)), y=count, fill=np))
g <- g + geom_bar(position="stack", stat="identity") + coord_flip() + theme_bw()
g <- g + scale_fill_manual(values = c("#000000", "#FFFFFF")) + geom_col(color = "black", linewidth=0.2)
g + ggtitle("top 200-edges") + xlab("top TFs")
dev.off()

##output files
edge_f <- paste0(title,".edges.A01.csv")
colnames(edge.mean) <- c("TF","target","mean.Assoc")
edge.sub <- subset(edge.mean, abs(edge.mean[,3]) >= 0.1)
write.csv(edge.sub, file=paste0("../",edge_f), row.names = F)
                                      

00-edges") + xlab("top TFs")
dev.off()

##output files
edge_f <- paste0(title,".edges.A01.csv")
colnames(edge.mean) <- c("TF","target","mean.Assoc")
edge.sub <- subset(edge.mean, abs(edge.mean[,3]) >= 0.1)
write.csv(edge.sub, file=paste0("../",edge_f), row.names = F)
                                      

