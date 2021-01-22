library(Seurat)
library(ggplot2)
library(patchwork)

coln <- colnames(read.csv('plot1.txt',sep='\t',check.names = F,row.names = 1))
info <- read.csv('plot1.txt',sep='\t',check.names = T,row.names = 1)
dataall <- readRDS('../all_umap_tsne_clusterok.RDS')
#bc <- readRDS('../../new/major/BC-new.RDS')
#crc <- readRDS('../../new/major/CRC-new.RDS')
#lc <- readRDS('../../new/major/LC-new.RDS')
#ov <- readRDS('../../new/major/OV-new.RDS')
#pdac <- readRDS('../../new/major/PDAC-new.RDS')
#scc <- readRDS('../../new/major/SCC-new.RDS')


#scc <- RunTSNE(scc,dims=1:30,verbose=T,tsne.method='FIt-SNE',nthreads=24,seed.use = 2020)

rst<-list()
idx=1

for (data in dataall){

orderidx <- rownames(data@meta.data)
for (item in colnames(info)){
  tmp <- info[,item]
  tmp <- tmp[tmp!='']
  tmp <- tmp[tmp %in% rownames(data)]
  print(tmp)
  tmpmean <- colMeans(as.data.frame(data@assays$SCT@data[tmp,]))[orderidx]
  data[[item]]<-tmpmean
}

a<-FeaturePlot(data,reduction = 'tsne',features = colnames(info),ncol = 8,cols = c('grey','red'))*NoLegend()

for (i in 1:7){
a[[i]]<-a[[i]]+theme(axis.title = element_blank(),axis.text = element_blank(),
           axis.ticks = element_blank(),plot.margin = unit(c(0,0,0,0),'cm'),
           plot.title = element_blank(),plot.title.position = 'panel')
}

b<-DimPlot(data,reduction = 'tsne')*theme(axis.title = element_blank(),axis.text = element_blank()
        ,axis.ticks = element_blank(),legend.position = 'none',
        plot.title = element_text(size=10,hjust = 0.5,face = 'bold'),plot.margin = unit(c(0,0,0,0),'cm'))
#b<-LabelClusters(b,id='ident')
p<-(a+b)+plot_layout(ncol = 8)

rst[[idx]]<-p
idx = idx+1
}

rst[[1]]/rst[[2]]/rst[[3]]/rst[[4]]/rst[[5]]/rst[[6]]

