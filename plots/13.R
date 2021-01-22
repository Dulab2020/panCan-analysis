library(ggplot2)
library(ggpubr)
library(yyplot)
library(ggsci)
library(Seurat)
library(patchwork)
options(max.print = 50)

setwd('~/research/2020single/finalplot/plot13-IFNG-gsva-M-violin/')

l<-readRDS('~/research/2020single/hypoxia_diff/gsvaall.RDS')
met<-readRDS('~/research/2020single/hypoxia_diff/metadata_rownames1.RDS')
data<-readRDS('../all_umap_tsne_clusterok.RDS')
gsvidx = data.frame(BC=1,CRC=2,LC=3,OV=4,PDAC=5,SCC=6)

strt<-function(x){
  if(is.na(x)){
    return(NA)
  }else
    if (startsWith(x,'M'))
      return('M')
  else
    return('others')
}

for (i in names(gsvidx)){
  can = i
  cat(can,'\n')
  hallmpath <- paste0('~/research/2020single/gsva/majorsub/results/done/m/hallmarks-',can,'-myeloid-SCT-data-gsva.txt')
  hallm <- read.csv(hallmpath,sep='\t')
  d <- data[[gsvidx[[i]]]]
  d <- AddMetaData(d,l[[gsvidx[[i]]]])
  d <- AddMetaData(d,hallm)
  m <- met[[gsvidx[[i]]]]
  d <- AddMetaData(d,m)
  #cat(unique(d@meta.data$CellFromTumor),'\n')
  d<-subset(d,CellFromTumor!="Normal")
  cat(dim(d),'\n')
  #high <- quantile(x = d@meta.data$Hypoxia,probs = 0.5)
  #low <- quantile(x = d@meta.data$Hypoxia,probs = 0.5)
  #d[['Hypoxia_tag']]<-sapply(d@meta.data$Hypoxia,FUN = hl,low=low,high=high)
  d[['mtag']] <- as.vector(sapply(as.vector(d@meta.data$cell_type),FUN=strt))
  md <- subset(d,mtag=='M')
  hypm <- md@meta.data[,c('HALLMARK_INTERFERON_GAMMA_RESPONSE','cell_type')]
  hypm$st<-sapply(hypm$cell_type,FUN = function(x){return(strsplit(x,'-',fixed = T)[[1]][2])})
  hypm<-hypm[sort(hypm$st,index.return=T)$ix,]
  outname<-paste0(can,'-M-hypoxia.pdf')
  p<-ggviolin(hypm,x='cell_type',y='HALLMARK_INTERFERON_GAMMA_RESPONSE',fill='cell_type')+NoLegend()+scale_fill_simpsons()+
    theme(axis.text.x = element_text(angle = 45,hjust = 1),axis.title = element_blank())
  g<-set_font(p,fontface = 'bold',size=6)
  pdf(outname,width = 3,height = 3)
  plot(g)
  dev.off()
  script <- paste0(can,"<-as_ggplot(g)")
  eval(parse(text = script))
}

p<-ggviolin(hypm,x='cell_type',y='',fill='cell_type')+NoLegend()+scale_fill_simpsons()+
  theme(axis.text.x = element_text(angle = 45,hjust = 1),axis.title = element_blank())
g<-set_font(p,fontface = 'bold',size=6)+scale_x_discrete(labels=c('M-S1','M-S3','M-S4','Mon-S7'))
SCC<-as_ggplot(g)

BC|CRC|LC|OV|PDAC|SCC











