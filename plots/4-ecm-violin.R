library(ggplot2)
library(ggpubr)
library(yyplot)
library(ggsci)
library(patchwork)
library(Seurat)
options(max.print = 50)

setwd('~/research/2020single/finalplot/plot4-scores-F/')
met<-readRDS('~/research/2020single/hypoxia_diff/metadata_rownames1.RDS')
data<-readRDS('../all_umap_tsne_clusterok.RDS')
gsvidx = data.frame(BC=1,CRC=2,LC=3,OV=4,PDAC=5,SCC=6)

strt<-function(x){
  if(is.na(x)){
    return(NA)
  }else
    if (startsWith(x,'FS'))
      return('F')
  else
    return('others')
}
#bc crc hallmark
i=1
d<-data[[i]]
d <- AddMetaData(d,met[[gsvidx[[i]]]])
m <- read.csv('../../gsva/majorsub/results/done/f/hallmarks-BC-fib-SCT-data-gsva.txt',sep='\t',row.names = 1)
d <- AddMetaData(d,m)
cat(unique(d@meta.data$CellFromTumor),'\n')
d<-subset(d,CellFromTumor!="Normal")
cat(dim(d),'\n')
d[['mtag']] <- as.vector(sapply(as.vector(d@meta.data$cell_type),FUN=strt))
md <- subset(d,mtag=='F')
hypm <- md@meta.data
p1<-ggviolin(hypm,x='cell_type',y='HALLMARK_INTERFERON_GAMMA_RESPONSE',fill='cell_type',
            title = 'Interferon gamma response')+NoLegend()+scale_fill_simpsons()+
  theme(axis.title.y = element_blank(),axis.title.x = element_blank(),axis.ticks = element_blank(),
        axis.text.x = element_blank())
g1<-set_font(p1,fontface = 'bold',size=6)
plot(g1)

p2<-ggviolin(hypm,x='cell_type',y='HALLMARK_INFLAMMATORY_RESPONSE',fill='cell_type',
             title = 'Inflammatory response')+NoLegend()+scale_fill_simpsons()+
  theme(axis.title.y = element_blank(),axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45,hjust = 1))
g2<-set_font(p2,fontface = 'bold',size=6)
plot(g2)

bc<-as_ggplot(g1)/as_ggplot(g2)+plot_layout(heights = c(1,1.15))

i=2
d<-data[[i]]
d <- AddMetaData(d,met[[gsvidx[[i]]]])
m <- read.csv('../../gsva/majorsub/results/done/f/hallmarks-CRC-fib-SCT-data-gsva.txt',sep='\t',row.names = 1)
d <- AddMetaData(d,m)
cat(unique(d@meta.data$CellFromTumor),'\n')
d<-subset(d,CellFromTumor!="Normal")
cat(dim(d),'\n')
d[['mtag']] <- as.vector(sapply(as.vector(d@meta.data$cell_type),FUN=strt))
md <- subset(d,mtag=='F')
hypm <- md@meta.data
hypm<-hypm[sort(hypm$cell_type,index.return=T)$ix,]
p3<-ggviolin(hypm,x='cell_type',y='HALLMARK_INTERFERON_GAMMA_RESPONSE',fill='cell_type',
             title = 'Interferon gamma response')+NoLegend()+scale_fill_simpsons()+
  theme(axis.title.y = element_blank(),axis.title.x = element_blank(),axis.ticks = element_blank(),
        axis.text.x = element_blank())
g3<-set_font(p3,fontface = 'bold',size=6)
plot(g3)

p4<-ggviolin(hypm,x='cell_type',y='HALLMARK_INFLAMMATORY_RESPONSE',fill='cell_type',
             title = 'Inflammatory response')+NoLegend()+scale_fill_simpsons()+
  theme(axis.title.y = element_blank(),axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45,hjust = 1))+
  scale_y_continuous(limits = c(-0.5,0.5),breaks = seq(-0.5,0.5,0.3))
g4<-set_font(p4,fontface = 'bold',size=6)
plot(g4)
crc<-as_ggplot(g3)/as_ggplot(g4)+plot_layout(heights = c(1,1.15))


#ECM remodeling score
genes <- read.csv('ecm.txt')$ECM.remodeling
#BC
i=1
can = names(gsvidx)[[i]]
d <- data[[i]]
m <- met[[i]]
d <- AddMetaData(d,m)
cat(unique(d@meta.data$CellFromTumor),'\n')
d<-subset(d,CellFromTumor!="Normal")
cat(dim(d),'\n')
d[['mtag']] <- as.vector(sapply(as.vector(d@meta.data$cell_type),FUN=strt))
md <- subset(d,mtag=='F')
mtx<-as.data.frame(md@assays$SCT@data)
ecmscores <- colMeans(mtx[genes,],na.rm = T)-colMeans(mtx)
md[['ECM']]<-ecmscores
hypm <- md@meta.data[,c('ECM','cell_type')]
hypm<-hypm[sort(hypm$cell_type,index.return=T)$ix,]
outname<-paste0(can,'-F-ECM.pdf')
p<-ggviolin(hypm,x='cell_type',y='ECM',fill='cell_type')+NoLegend()+scale_fill_simpsons()+
  theme(axis.title = element_blank(),axis.text.x = element_text(angle = 45,hjust = 1))+
  scale_y_continuous(n.breaks = 3)
g<-set_font(p,fontface = 'bold',size=6)
bcf<-as_ggplot(g)
#CRC
i=2
can = names(gsvidx)[[i]]
d <- data[[i]]
m <- met[[i]]
d <- AddMetaData(d,m)
cat(unique(d@meta.data$CellFromTumor),'\n')
d<-subset(d,CellFromTumor!="Normal")
cat(dim(d),'\n')
d[['mtag']] <- as.vector(sapply(as.vector(d@meta.data$cell_type),FUN=strt))
md <- subset(d,mtag=='F')
mtx<-as.data.frame(md@assays$SCT@data)
ecmscores <- colMeans(mtx[genes,],na.rm = T)-colMeans(mtx)
md[['ECM']]<-ecmscores
hypm <- md@meta.data[,c('ECM','cell_type')]
hypm<-hypm[sort(hypm$cell_type,index.return=T)$ix,]
outname<-paste0(can,'-F-ECM.pdf')
p<-ggviolin(hypm,x='cell_type',y='ECM',fill='cell_type')+NoLegend()+scale_fill_simpsons()+
  theme(axis.title = element_blank(),axis.text.x = element_text(angle = 45,hjust = 1))+
  scale_y_continuous(n.breaks = 4)
g<-set_font(p,fontface = 'bold',size=6)
crcf<-as_ggplot(g)
#LC
i=3
can = names(gsvidx)[[i]]
d <- data[[i]]
m <- met[[i]]
d <- AddMetaData(d,m)
cat(unique(d@meta.data$CellFromTumor),'\n')
d<-subset(d,CellFromTumor!="Normal")
cat(dim(d),'\n')
d[['mtag']] <- as.vector(sapply(as.vector(d@meta.data$cell_type),FUN=strt))
md <- subset(d,mtag=='F')
mtx<-as.data.frame(md@assays$SCT@data)
ecmscores <- colMeans(mtx[genes,],na.rm = T)-colMeans(mtx)
md[['ECM']]<-ecmscores
hypm <- md@meta.data[,c('ECM','cell_type')]
hypm<-hypm[sort(hypm$cell_type,index.return=T)$ix,]
outname<-paste0(can,'-F-ECM.pdf')
p<-ggviolin(hypm,x='cell_type',y='ECM',fill='cell_type')+NoLegend()+scale_fill_simpsons()+
  theme(axis.title = element_blank(),axis.text.x = element_text(angle = 45,hjust = 1))+
  scale_y_continuous(n.breaks = 4)
g<-set_font(p,fontface = 'bold',size=6)
lcf<-as_ggplot(g)
#OV
i=4
can = names(gsvidx)[[i]]
d <- data[[i]]
m <- met[[i]]
d <- AddMetaData(d,m)
cat(unique(d@meta.data$CellFromTumor),'\n')
d<-subset(d,CellFromTumor!="Normal")
cat(dim(d),'\n')
d[['mtag']] <- as.vector(sapply(as.vector(d@meta.data$cell_type),FUN=strt))
md <- subset(d,mtag=='F')
mtx<-as.data.frame(md@assays$SCT@data)
ecmscores <- colMeans(mtx[genes,],na.rm = T)-colMeans(mtx)
md[['ECM']]<-ecmscores
hypm <- md@meta.data[,c('ECM','cell_type')]
hypm<-hypm[sort(hypm$cell_type,index.return=T)$ix,]
outname<-paste0(can,'-F-ECM.pdf')
p<-ggviolin(hypm,x='cell_type',y='ECM',fill='cell_type')+NoLegend()+scale_fill_simpsons()+
  theme(axis.title = element_blank(),axis.text.x = element_text(angle = 45,hjust = 1))+
  scale_y_continuous(n.breaks = 4)
g<-set_font(p,fontface = 'bold',size=6)
ovf<-as_ggplot(g)
#PDAC
i=5
can = names(gsvidx)[[i]]
d <- data[[i]]
m <- met[[i]]
d <- AddMetaData(d,m)
cat(unique(d@meta.data$CellFromTumor),'\n')
d<-subset(d,CellFromTumor!="Normal")
cat(dim(d),'\n')
d[['mtag']] <- as.vector(sapply(as.vector(d@meta.data$cell_type),FUN=strt))
md <- subset(d,mtag=='F')
mtx<-as.data.frame(md@assays$SCT@data)
ecmscores <- colMeans(mtx[genes,],na.rm = T)-colMeans(mtx)
md[['ECM']]<-ecmscores
hypm <- md@meta.data[,c('ECM','cell_type')]
hypm<-hypm[sort(hypm$cell_type,index.return=T)$ix,]
outname<-paste0(can,'-F-ECM.pdf')
p<-ggviolin(hypm,x='cell_type',y='ECM',fill='cell_type')+NoLegend()+scale_fill_simpsons()+
  theme(axis.title = element_blank(),axis.text.x = element_text(angle = 45,hjust = 1))+
  scale_y_continuous(n.breaks = 4)
g<-set_font(p,fontface = 'bold',size=6)
pdacf<-as_ggplot(g)






