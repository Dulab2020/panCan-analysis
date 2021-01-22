library(ggplot2)
library(ggpubr)
library(ggsci)
library(reshape2)
pathw <- read.csv('cor-gsva.txt',sep='\t')$gene2
newp <- c('WNT beta catenin signaling','E2F targets','Angiogenesis','Hypoxia','Group')
newp2 <- c('Angiogenesis','E2F targets','Hypoxia','WNT beta catenin signaling','Group')
#LC
lc2ecm <- read.csv('zscore-LC-2.csvdata.csv',row.names=1)
lc3cce <- read.csv('zscore-LC-3.csvdata.csv',row.names=1)
lc2ecm$group<-sapply(lc2ecm$group,function(x){paste0('ECM ',x)})
lc3cce$group<-sapply(lc3cce$group,function(x){paste0('Cell cycle ',x)})
lc3cce<-lc3cce[sort(lc3cce$group,decreasing = T,index.return=T)$ix,]
lc2ecm<-lc2ecm[sort(lc2ecm$group,decreasing = T,index.return=T)$ix,]

lc <- rbind(lc3cce,lc2ecm)[,c(pathw,'group')]
colnames(lc)<-newp
lc<-lc[,newp2]
m<-melt(lc,id.vars = c('Group'),variable.name = 'Pathways',value.name = 'Pathway score value')


ggboxplot(m,x='Group',y='Pathway score value',facet.by = 'Pathways',color = 'Group',
          nrow=1,size = 1.5,width = 0.6,bxp.errorbar = T,notch = T)+scale_color_simpsons()+
  stat_compare_means(comparisons = list(c(1,2),c(3,4)),label = 'p.signif')+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),
        legend.position = 'bottom',axis.title.x = element_blank(),text = element_text(size=12,face = 'bold'))

#others
files <- dir('.',pattern = 'csv$')
files <- files[-2:-3]

for (file in files){
  base <- strsplit(strsplit(file,'zscore-')[[1]][2],'.',fixed = T)[[1]][1]
  df <- read.csv(file,row.names = 1)
  df <- df[df$group!='mid',c(pathw,'group')]
  colnames(df) <- newp
  df <- df[sort(df$Group,decreasing = T,index.return=T)$ix,newp2]
  m <- melt(df,id.vars = c('Group'),variable.name = 'Pathways',value.name = 'Pathway score')
  pdf(paste0(base,'.pdf'),width = 8,height = 4)
  p<-ggboxplot(m,x='Group',y='Pathway score',facet.by = 'Pathways',color = 'Group',
            nrow=1,size = 1.5,width = 0.6,bxp.errorbar = T,notch = T)+scale_color_simpsons()+
    stat_compare_means(comparisons = list(c(1,2)),label = 'p.signif')+
    theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),
          legend.position = 'bottom',axis.title.x = element_blank(),text = element_text(size=12,face = 'bold'))
  print(p)
  dev.off()
}




