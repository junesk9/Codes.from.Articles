##2020.09.29
##Both boxplot & barplot of the root measure data


########################################################
############### Boxplot
library("beeswarm")
d <- read.table(pipe("pbpaste"), header=T)
head(d)
#        WT     C04
#1 66.76425 4.11525
#2 67.83975 4.63725
#3 56.16450 4.34400
#4 63.06600 3.95775
#5 60.94650 3.71925

boxplot(d, ylim=c(0,80), main="root length (n=29)", outline=FALSE)
beeswarm(d, add=TRUE, corral="gutter", pch=19, cex=.5, col="red")


#########################################################
############### Barplot with points
library(tidyverse)
d <- read.table(pipe("pbpaste"), header=T)
dt <- d %>% gather(key=key, value=value) #Spread the table; easier way to points(), anova()

#way1
#means <- apply(d,2, FUN=mean, na.rm=TRUE)
#sds <- apply(d,2, FUN=sd, na.rm=TRUE)
#way2
means <- tapply(dt$value, dt$key, FUN=mean, na.rm=TRUE)
sds <- tapply(dt$value, dt$key, FUN=sd, na.rm=TRUE)

x <- barplot(means, ylim=c(0,1.2*max(means, na.rm=TRUE)), col=0)
##Add points **jitter() can be ignored with small numbers of spot
points(jitter(rep(x, table(dt$key)), factor = 0.3), dt$value[order(dt$key)], pch=19, cex=0.5, col=2)
arrows(x, means-sds, x, means+sds, code=3, lwd=1, angle=90, length=0.2)

##if not spread the table, use "rep(x, table(rep(colnames(d), length(d[,1]))))"


#########################################################
################# ANOVA Tukey HSD
####ANOVA w/ Tukey HSD

anova(aov(value~key, data=dt))
#Analysis of Variance Table
#
#Response: value
#Df Sum Sq Mean Sq F value    Pr(>F)    
#key        1 2735.0 2735.05  18.389 0.0001548 ***
#  Residuals 32 4759.4  148.73                      
#---

p <- TukeyHSD(aov(value~key, data=dt))
p
#Tukey multiple comparisons of means
#95% family-wise confidence level#
#
#Fit: aov(formula = value ~ key, data = dt)
#
#$key
#diff      lwr      upr     p adj
#WT-bz 17.96907 9.433761 26.50438 0.0001548

####################################################
############### Scatter DEG-DEG

d <- read.table(pipe("pbpaste"), header=T) ## 200-bottom DEG of c04
dd <- subset(d, d[,3]<0.05 & d[,2] <= -1)
d2 <- read.table(pipe("pbpaste"), header=T)
dd2 <- subset(d2, d2[,3]<0.05 & d2[,2] >= 1) ## 200-top DEG of c04

dim(dd)
#[1] 166   3
dim(dd2)
#[1] 46  3
head(d)
#log2FC nbr6.log2FC     nbr6.fdr
#1 -2.411437   -2.071384 4.298630e-05
#2 -2.417207   -2.480188 1.938350e-06
#3 -2.417239   -3.481140 5.107350e-05
#4 -2.423656   -2.594325 9.854770e-08
#5 -2.426503   -1.555867 1.143598e-02

###### non-DEG subset
library(dplyr)
du <- setdiff(d, dd)
du2 <- setdiff(d2, dd2)

###### Scatter plot
plot(du[,1],du[,2], pch=20, cex=0.5, ylim=c(-10,2), xlim=c(-10,0), col="grey80")
points(dd[,1],dd[,2], pch=20, cex=0.5, col="#3498DB")
lines(x=c(-10,-0),y=c(-10,0),col=2)
abline(h=0)
text(-10,-10, length(dd[,1]), pos=4, cex=1.4, col="#3498DB")

plot(du2[,1],du2[,2], pch=20, cex=0.5, ylim=c(-2,10), xlim=c(0,10), col="grey80")
points(dd2[,1],dd2[,2], pch=20, cex=0.5, col="#ED0422")
lines(x=c(0,10),y=c(0,10),col=2)
abline(h=0)
text(2,9, length(dd2[,1]), pos=4, cex=1.4, col="#ED0422")


