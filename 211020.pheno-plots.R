############ Various plot for Phenotype data
########### SV274 GWAS
########### 2021.10.20 Juensk9

################# input data
d <- read.table(pipe("pbpaste"), header=T)
d <- d[,1]
d <- d[!is.na(d)]

####################histogram with normal curve
h <- hist(d, breaks=10, density=10, freq=FALSE)
xfit <- seq(min(d), max(d), length=40)
yfit <- dnorm(xfit, mean = mean(d), sd = sd(d)) 
#yfit <- yfit * diff(h$mids[1:2]) * length(d) 

lines(xfit, yfit, col = "black", lwd = 2)


###################histogram version I
h <- hist(d)
h$density = h$counts/sum(h$counts)*100
plot(h,freq=FALSE, col="grey80")


##################violin plot
#install.packages("vioplot", dependencies = TRUE)

library(vioplot)
d <- d[!is.na(d)]
vioplot(d, side="right", col="grey80")
