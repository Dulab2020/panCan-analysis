options(max.print = 50)
setwd('~/research/2020single/finalplot/plot10-tsne-cancer/')
library(Seurat)
library(ggplot2)
library(patchwork)


l<-readRDS('~/research/2020single/hypoxia_diff/gsvaall.RDS')
met<-readRDS('~/research/2020single/hypoxia_diff/metadata_rownames1.RDS')
data<-readRDS('../all_umap_tsne_clusterok-cancer.RDS')
gsvidx <- data.frame(BC=1,CRC=2,LC=3,OV=4,PDAC=5,SCC=6)
scores <- read.csv('genesets.txt',sep='\t')
ecm <- scores$all.ECM[1:67]
cce <- scores$all.CellCycle


for (can in names(gsvidx)){
  cat(can,'\n')
  i<- gsvidx[[can]]
  d <- data[[i]]
  d <- AddMetaData(d,l[[i]])
  d <- AddMetaData(d,met[[i]])
  mtx <- as.data.frame(d$SCT@data)
  ecmscore <- colMeans(mtx[ecm,],na.rm = T)-colMeans(mtx)
  ccescore <- colMeans(mtx[cce,],na.rm = T)-colMeans(mtx)
  d[['ECM']]<-scale(ecmscore)
  d[['CellCycle']]<-scale(ccescore)
  #d[['ECM']]<-ecmscore
  #d[['CellCycle']]<-ccescore
  if(can!='BC'){
  a1<-FeaturePlot(d,reduction = 'tsne',features = c('ECM'))+scale_color_gradientn(limits=c(-1,5),colours = c('grey','red'))
  a2<-FeaturePlot(d,reduction = 'tsne',features = c('CellCycle'))+scale_color_gradientn(limits=c(-1,5),colours = c('grey','red'))

  a1<-a1+NoLegend()+labs(title = can)+theme(axis.title = element_blank(),axis.text = element_blank(),
                                               axis.ticks = element_blank(),plot.margin = unit(c(0,0,0,0),'cm'),
                                               plot.title = element_text(size=10),plot.title.position = 'panel')
  a2<-a2+NoLegend()+labs(title = element_blank())+theme(title = element_blank(),axis.title.x = element_blank(),axis.text = element_blank(),
                                            axis.ticks = element_blank(),plot.margin = unit(c(0,0,0,0),'cm'),
                                            plot.title = element_text(size=10),plot.title.position = 'panel')
  
  
  b<-DimPlot(d,reduction = 'tsne')*theme(title = element_blank(),axis.title.x = element_blank(),axis.text = element_blank()
                                                                         ,axis.ticks = element_blank(),legend.position = 'none',
                                                                         plot.title = element_text(size=10,hjust = 0.5,face = 'bold'),plot.margin = unit(c(0,0,0,0),'cm'))
  
  b<-LabelClusters(b,id='ident')
  p<-a1+a2+b
  script <- paste0(tolower(can),'<-p')
  eval(parse(text = script))
  }else{
    a1<-FeaturePlot(d,reduction = 'tsne',features = c('ECM'))+scale_color_gradientn(limits=c(-1,5),colours = c('grey','red'))
    a2<-FeaturePlot(d,reduction = 'tsne',features = c('CellCycle'))+scale_color_gradientn(limits=c(-1,5),colours = c('grey','red'))
    
    a1<-a1+NoLegend()+ylab(label = 'EMT score')+labs(title = can)+theme(axis.title.x = element_blank(),axis.text = element_blank(),
                                                                        axis.ticks = element_blank(),plot.margin = unit(c(0,0,0,0),'cm'),
                                                                        plot.title = element_text(size=10),plot.title.position = 'panel')
    a2<-a2+NoLegend()+ylab(label = 'Cell Cycle score')+labs(title = element_blank())+theme(axis.title.x = element_blank(),axis.text = element_blank(),
                                                                                           axis.ticks = element_blank(),plot.margin = unit(c(0,0,0,0),'cm'),
                                                                                           plot.title = element_text(size=10),plot.title.position = 'panel')
    
    
    b<-DimPlot(d,reduction = 'tsne')*theme(axis.title.x = element_blank(),axis.text = element_blank()
                                           ,axis.ticks = element_blank(),legend.position = 'none',
                                           plot.title = element_text(size=10,hjust = 0.5,face = 'bold'),plot.margin = unit(c(0,0,0,0),'cm'))
    
    b<-LabelClusters(b,id='ident')
    b<-b+ylab("Seurat clusters")
    p<-a1+a2+b
    script <- paste0(tolower(can),'<-p')
    eval(parse(text = script))
  }
}
bc|crc|lc|ov|pdac|scc
