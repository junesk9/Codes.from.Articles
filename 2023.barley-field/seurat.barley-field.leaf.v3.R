#############################################################
########### Seurat output visualization & characterization
######################### Junesk9 2024.02.15

############################
######## Load libraries required
###########################

#remotes::install_version("SeuratObject", "4.1.4", repos = c("https://satijalab.r-universe.dev", getOption("repos")))
#remotes::install_version("Seurat", "4.4.0", repos = c("https://satijalab.r-universe.dev", getOption("repos")))
library(Seurat) #ver.4.4.0
library(cluster) #for clusGAP
library(gplots) #for heatmap.2
library(glmnet) #for glm

#For GESA analysis
library(GSEABase) 
library(GOstats) #GO/KEGG
library(KEGGREST) #KEGG-API

#Accessories
library(dunn.test)
library(beeswarm)
`%notin%` <- Negate(`%in%`) #temp. function to "not-in"


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
################################################# 
###################################### Load data
################################################