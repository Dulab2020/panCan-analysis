options(max.print = 50)
library(Seurat)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(ggsci)
library(yyplot)
setwd('~/research/2020single/finalplot/plot5-Tcell-scores/')
l<-readRDS('~/research/2020single/hypoxia_diff/gsvaall.RDS')
met<-readRDS('~/research/2020single/hypoxia_diff/metadata_rownames1.RDS')
data<-readRDS('../all_umap_tsne_clusterok.RDS')
gsvidx = data.frame(BC=1,CRC=2,LC=3,OV=4,PDAC=5,SCC=6)

strt<-function(x){
  if(is.na(x)){
    return(NA)
  }else
    if (startsWith(x,'CD8'))
      return('CD8')
  else
    return('others')
}

genes<-read.csv('CD8score.txt',sep='\t')

for (i in 1:6){
  can <- names(gsvidx)[i]
  cat(can,'\n')
  d <- data[[i]]
  m <- met[[i]]
  d <- AddMetaData(d,l[[i]])
  d <- AddMetaData(d,m)
  d<-subset(d,CellFromTumor!="Normal")
  d[['mtag']] <- as.vector(sapply(as.vector(d@meta.data$cell_type),FUN=strt))
  md <- subset(d,mtag=='CD8')
  mtx<-as.data.frame(md@assays$SCT@data)
  cyto <- colMeans(mtx[genes$Cytotoxicity[1:12],],na.rm = T)-colMeans(mtx)
  exh <- colMeans(mtx[genes$Exhaustion,],na.rm = T)-colMeans(mtx)
  md[['Cytotoxicity']]<-cyto
  md[['Exhaustion']]<-exh
  hypm <- md@meta.data[,c('Cytotoxicity','cell_type','Exhaustion','Hypoxia')]
  hypm<-hypm[sort(hypm$cell_type,index.return=T)$ix,]
  p<-ggviolin(hypm,x='cell_type',y='Cytotoxicity',fill='cell_type')+NoLegend()+scale_fill_simpsons()+
    theme(axis.title = element_blank(),axis.text.x = element_text(angle = 45,hjust = 1))+
    scale_y_continuous(n.breaks = 4)
  g<-set_font(p,fontface = 'bold',size=6)
  script <- paste0(can,'c<-as_ggplot(g)')
  eval(parse(text = script))
  p<-ggviolin(hypm,x='cell_type',y='Exhaustion',fill='cell_type')+NoLegend()+scale_fill_simpsons()+
    theme(axis.title = element_blank(),axis.text.x = element_text(angle = 45,hjust = 1))+
    scale_y_continuous(n.breaks = 4)
  g<-set_font(p,fontface = 'bold',size=6)
  script <- paste0(can,'e<-as_ggplot(g)')
  eval(parse(text = script))
  write.csv(hypm,paste0(can,'-data.csv'))
}



