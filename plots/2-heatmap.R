options(java.parameters = '-Xmx100g')
library(Seurat)
library(pheatmap)
library(xlsx)
library(dplyr)


setwd('~/research/2020single/finalplot/plot2-heatmap/')
sheetnames <- names(getSheets(loadWorkbook('markers.xlsx')))
gsvall<-readRDS('../../hypoxia_diff/gsvaall.RDS')
gsvidx <- data.frame(BC=1,CRC=2,LC=3,OV=4,PDAC=5,SCC=6)

for (i in 1:17){
sheetname <- sheetnames[i]
can <- strsplit(sheetname,'-',fixed = T)[[1]][1]
cell <- strsplit(sheetname,'-',fixed = T)[[1]][2]
info <- read.xlsx('markers.xlsx',i)

genes <- c()
for (j in colnames(info)){
  genes<-c(genes,c(na.omit(c(info[,j]))))
}

if(cell=='T'){
  pref = '../../new/fine/t/'
}else if(cell=='F'){
  pref = '../../new/fine/fib/'
}else if(cell=='M'){
  pref = '../../new/fine/myeloid/'
}

datapath <- paste0(pref,can,'-final.RDS')
cat(datapath,'\n')

data <- readRDS(datapath)
cat(data@commands$FindClusters$resolution,'\n')

gsv <- gsvall[[gsvidx[can][[1]]]]
data<-AddMetaData(data,gsv)

g<-data@meta.data[,c('seurat_clusters','Hypoxia')]
rst<-as.data.frame(summarise(group_by(g,seurat_clusters),mean_hypoxia_score=mean(Hypoxia)))

hdata<-AverageExpression(data,assays = "SCT",slot = 'data')

plotdata <- as.data.frame(hdata$SCT)

write.xlsx(x = plotdata,file = 'plotdata-data-all.xlsx',sheetName = sheetname,append = T)
write.xlsx(x = rst,file = 'plotdata-data-all.xlsx',sheetName = paste0(sheetname,'-hypoxia'),row.names = F,append = T)
}


