library(pheatmap)
library(ggplot2)
library(tidyr)
library(patchwork)

pathways <- read.csv('pathways.txt',sep='\t')$pathways
pathways <- as.vector(sapply(pathways,function(x){strsplit(x,'HALLMARK_')[[1]][2]})) %>% 
  gsub(pattern = '_',replacement = ' ')
new<-read.csv('pathways.txt',sep = '\t')$new
files <- dir('data/')

bk <- c(seq(-1.5,0,length.out = 50),seq(0.01,1.5,length.out = 50))

#bc
df <- t(read.csv('data/hallmarks-BC-cancer-ave-z.txt',sep='\t',row.names = 1,check.names = F))
df <- df[pathways,]
rownames(df)<-new
p <- pheatmap(df,cluster_rows = F,cellwidth = 12,cellheight = 12,treeheight_col = 10,
              show_rownames = F,breaks = bk,legend = F,border_color = NA,angle_col = 0)
a <- ggplotify::as.ggplot(p)

#crc
df <- t(read.csv('data/hallmarks-CRC-cancer-ave-z.txt',sep='\t',row.names = 1,check.names = F))
df <- df[pathways,]
rownames(df)<-new
p <- pheatmap(df,cluster_rows = F,cellwidth = 12,cellheight = 12,treeheight_col = 10,
              show_rownames = F,breaks = bk,legend = F,border_color = NA,angle_col = 0)
b <- ggplotify::as.ggplot(p)

#lc
df <- t(read.csv('data/hallmarks-LC-cancer-ave-z.txt',sep='\t',row.names = 1,check.names = F))
df <- df[pathways,]
rownames(df)<-new
p <- pheatmap(df,cluster_rows = F,cellwidth = 12,cellheight = 12,treeheight_col = 10,
              show_rownames = F,breaks = bk,legend = F,border_color = NA,angle_col = 0)
c <- ggplotify::as.ggplot(p)

#OV
df <- t(read.csv('data/hallmarks-OvC-cancer-ave-z.txt',sep='\t',row.names = 1,check.names = F))
df <- df[pathways,]
rownames(df)<-new
p <- pheatmap(df,cluster_rows = F,cellwidth = 12,cellheight = 12,treeheight_col = 10,
              show_rownames = F,breaks = bk,legend = F,border_color = NA,angle_col = 0)
d <- ggplotify::as.ggplot(p)

#pdac
df <- t(read.csv('data/hallmarks-PDAC-cancer-ave-z.txt',sep='\t',row.names = 1,check.names = F))
df <- df[pathways,]
rownames(df)<-new
p <- pheatmap(df,cluster_rows = F,cellwidth = 12,cellheight = 12,treeheight_col = 10,
              show_rownames = F,breaks = bk,legend = F,border_color = NA,angle_col = 0)
e <- ggplotify::as.ggplot(p)

#scc
df <- t(read.csv('data/hallmarks-SCC-cancer-ave-z.txt',sep='\t',row.names = 1,check.names = F))
df <- df[pathways,]
rownames(df)<-new
p <- pheatmap(df,cluster_rows = F,cellwidth = 12,cellheight = 12,treeheight_col = 10,
              show_rownames = T,breaks = bk,legend = T,border_color = NA,angle_col = 0)
f <- ggplotify::as.ggplot(p)

a+b+c+d+e+f+plot_layout(nrow = 1,widths = c(0.1,0.2,0.1,0.3,0.2,0.8))
