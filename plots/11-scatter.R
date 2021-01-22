setwd('~/research/2020single/finalplot/plot11-L-R-supp/supp3')
options(max.print = 50)
library(ggplot2)
library(data.table)
library(ggpubr)
library(patchwork)
##data prep
tumors = read.csv('../all_tumor_sample_sheet.csv',check.names = F)
samp <- tumors[tumors$Project=='TCGA-LUAD','Sample ID']
df <- fread('../alllog2tpmp1-genename.csv',sep=',')
mtx <- as.matrix(df[,-1])
rownames(mtx)<-df$Gene_name
rm(df)
gc()
mtx <- mtx[,samp]
mtx <- as.data.frame(mtx)
can='LUAD'
######SPP1
spp <- t(mtx['SPP1',,drop=F])
rownames(spp) <- sapply(rownames(spp),FUN=function(x){return(substring(x,1,12))})

#1.spp1-glycolysis
gly <- read.csv('../Pan-cancer.txt',sep='\t',row.names = 1)
gly <- gly[gly$cancer_type=='LUAD','glycolysis_score',drop=F]
dt1 <- merge(spp,gly,by=0)
test <- cor.test(dt1$glycolysis_score,dt1$SPP1,method = 'spearman',exact = F)
sr <- sprintf('%0.3f',test$estimate )
sp <- sprintf('%0.3e',test$p.value)
anno <- paste0("Spearman R=",sr,'\npVal=',sp)
p1<-ggscatter(dt1,x='glycolysis_score',y='SPP1')+labs(title = can)+
  stat_smooth(method = lm,color='red',fill='lightpink')+
  theme(text = element_text(size=12,face = 'bold'),plot.title = element_text(hjust = 0.5),
        axis.ticks = element_blank())+ylim(limit=c(2.5,14.5))+
  annotate(geom = 'text',label=anno,fontface='bold',hjust=0,x=-0.8,y=14)
p1
#2.spp1-SLC2A1 LDHA TGFBI
dt2 <- as.data.frame(t(mtx[c('SPP1','LDHA','SLC2A1','TGFBI'),]))

test <- cor.test(dt2$LDHA,dt2$SPP1,method = 'spearman',exact = F)
sr <- sprintf('%0.3f',test$estimate )
sp <- sprintf('%0.3e',test$p.value)
anno <- paste0("Spearman R=",sr,'\npVal=',sp)
p2_1<-ggscatter(dt2,x='LDHA',y='SPP1')+labs(title = can)+
  stat_smooth(method = lm,color='red',fill='lightpink')+
  theme(text = element_text(size=12,face = 'bold'),plot.title = element_text(hjust = 0.5),
        axis.ticks = element_blank())+ylim(limit=c(2.5,14.5))+
  annotate(geom = 'text',label=anno,fontface='bold',hjust=0,x=5,y=14)
p2_1

test <- cor.test(dt2$SLC2A1,dt2$SPP1,method = 'spearman',exact = F)
sr <- sprintf('%0.3f',test$estimate )
sp <- sprintf('%0.3e',test$p.value)
anno <- paste0("Spearman R=",sr,'\npVal=',sp)
p2_2<-ggscatter(dt2,x='SLC2A1',y='SPP1')+labs(title = can)+
  stat_smooth(method = lm,color='red',fill='lightpink')+
  theme(text = element_text(size=12,face = 'bold'),plot.title = element_text(hjust = 0.5),
        axis.ticks = element_blank())+ylim(limit=c(2.5,14.5))+
  annotate(geom = 'text',label=anno,fontface='bold',hjust=0,x=2,y=14)
p2_2

test <- cor.test(dt2$TGFBI,dt2$SPP1,method = 'spearman',exact = F)
sr <- sprintf('%0.3f',test$estimate )
sp <- sprintf('%0.3e',test$p.value)
anno <- paste0("Spearman R=",sr,'\npVal=',sp)
p2_3<-ggscatter(dt2,x='TGFBI',y='SPP1')+labs(title = can)+
  stat_smooth(method = lm,color='red',fill='lightpink')+
  theme(text = element_text(size=12,face = 'bold'),plot.title = element_text(hjust = 0.5),
        axis.ticks = element_blank())+ylim(limit=c(2.5,14.5))+
  annotate(geom = 'text',label=anno,fontface='bold',hjust=0,x=2,y=14)
p2_3

p2 <- p2_1+p2_2+p2_3
p2

#3.spp1-prolife
pro <- read.csv('../all-exh-pro.gsva.csv',row.names = 1)[samp,]
rownames(pro) <- sapply(rownames(pro),FUN=function(x){return(substring(x,1,12))})
dt3 <- merge(pro,spp,0)
test <- cor.test(dt3$Proliferation,dt3$SPP1,method = 'spearman',exact = F)
sr <- sprintf('%0.3f',test$estimate )
sp <- sprintf('%0.3e',test$p.value)
anno <- paste0("Spearman R=",sr,'\npVal=',sp)
p3<-ggscatter(dt3,x='Proliferation',y='SPP1')+labs(title = can)+
  stat_smooth(method = lm,color='red',fill='lightpink')+
  theme(text = element_text(size=12,face = 'bold'),plot.title = element_text(hjust = 0.5),
        axis.ticks = element_blank())+ylim(limit=c(2.5,14.5))+
  annotate(geom = 'text',label=anno,fontface='bold',hjust=0,x=-1,y=14)
p3

#4.spp1-Hallmark EMT
emt <- as.data.frame(t(read.csv('../hallmark-gsva.txt',row.names = 1,sep='\t',check.names = F)))
emt <- emt[sapply(samp,FUN=function(x){return(substring(x,1,12))}),'EMT markers',drop=F]
dt4 <- merge(emt,spp,0)
test <- cor.test(dt4$`EMT markers`,dt4$SPP1,method = 'spearman',exact = F)
sr <- sprintf('%0.3f',test$estimate )
sp <- sprintf('%0.3e',test$p.value)
anno <- paste0("Spearman R=",sr,'\npVal=',sp)
p4<-ggscatter(dt4,x='EMT markers',y='SPP1',xlab = 'EMT')+labs(title = can)+
  stat_smooth(method = lm,color='red',fill='lightpink')+
  theme(text = element_text(size=12,face = 'bold'),plot.title = element_text(hjust = 0.5),
        axis.ticks = element_blank())+ylim(limit=c(2.5,14.5))+
  annotate(geom = 'text',label=anno,fontface='bold',hjust=0,x=-1,y=14)
p4

######TNFSF12
tnf <- as.data.frame(t(mtx['TNFSF12',,drop=F]))
rownames(tnf) <- sapply(rownames(tnf),FUN=function(x){return(substring(x,1,12))})
#1.tnfsf12-CCL5 TIMP1 TGFBI VEGFA
dt5 <- as.data.frame(t(mtx[c('TNFSF12','CCL5','TIMP1','TGFBI','VEGFA'),]))

test <- cor.test(dt5$CCL5,dt5$TNFSF12,method = 'spearman',exact = F)
sr <- sprintf('%0.3f',test$estimate )
sp <- sprintf('%0.3e',test$p.value)
anno <- paste0("Spearman R=",sr,'\npVal=',sp)
p5_1<-ggscatter(dt5,x='CCL5',y='TNFSF12')+labs(title = can)+
  stat_smooth(method = lm,color='red',fill='lightpink')+
  theme(text = element_text(size=12,face = 'bold'),plot.title = element_text(hjust = 0.5),
        axis.ticks = element_blank())+ylim(limit=c(1,8))+
  annotate(geom = 'text',label=anno,fontface='bold',hjust=0,x=2,y=7.5)
p5_1

test <- cor.test(dt5$TIMP1,dt5$TNFSF12,method = 'spearman',exact = F)
sr <- sprintf('%0.3f',test$estimate )
sp <- sprintf('%0.3e',test$p.value)
anno <- paste0("Spearman R=",sr,'\npVal=',sp)
p5_2<-ggscatter(dt5,x='TIMP1',y='TNFSF12')+labs(title = can)+
  stat_smooth(method = lm,color='red',fill='lightpink')+
  theme(text = element_text(size=12,face = 'bold'),plot.title = element_text(hjust = 0.5),
        axis.ticks = element_blank())+ylim(limit=c(1,8))+
  annotate(geom = 'text',label=anno,fontface='bold',hjust=0,x=6,y=7.5)
p5_2

test <- cor.test(dt5$TGFBI,dt5$TNFSF12,method = 'spearman',exact = F)
sr <- sprintf('%0.3f',test$estimate )
sp <- sprintf('%0.3e',test$p.value)
anno <- paste0("Spearman R=",sr,'\npVal=',sp)
p5_3<-ggscatter(dt5,x='TGFBI',y='TNFSF12')+labs(title = can)+
  stat_smooth(method = lm,color='red',fill='lightpink')+
  theme(text = element_text(size=12,face = 'bold'),plot.title = element_text(hjust = 0.5),
        axis.ticks = element_blank())+ylim(limit=c(1,8))+
  annotate(geom = 'text',label=anno,fontface='bold',hjust=0,x=2.8,y=7.5)
p5_3

test <- cor.test(dt5$VEGFA,dt5$TNFSF12,method = 'spearman',exact = F)
sr <- sprintf('%0.3f',test$estimate )
sp <- sprintf('%0.3e',test$p.value)
anno <- paste0("Spearman R=",sr,'\npVal=',sp)
p5_4<-ggscatter(dt5,x='VEGFA',y='TNFSF12')+labs(title = can)+
  stat_smooth(method = lm,color='red',fill='lightpink')+
  theme(text = element_text(size=12,face = 'bold'),plot.title = element_text(hjust = 0.5),
        axis.ticks = element_blank())+ylim(limit=c(1,8))+
  annotate(geom = 'text',label=anno,fontface='bold',hjust=0,x=2.8,y=7.5)
p5_4

p5 <- p5_1+p5_3+p5_2+p5_4+plot_layout(nrow = 1)
p5

#2.tnfsf12-prolife
pro <- read.csv('../all-exh-pro.gsva.csv',row.names = 1)[samp,]
rownames(pro) <- sapply(rownames(pro),FUN=function(x){return(substring(x,1,12))})
dt6 <- merge(pro,tnf,0)
test <- cor.test(dt6$Proliferation,dt6$TNFSF12,method = 'spearman',exact = F)
sr <- sprintf('%0.3f',test$estimate )
sp <- sprintf('%0.3e',test$p.value)
anno <- paste0("Spearman R=",sr,'\npVal=',sp)
p6<-ggscatter(dt6,x='Proliferation',y='TNFSF12')+labs(title = can)+
  stat_smooth(method = lm,color='red',fill='lightpink')+
  theme(text = element_text(size=12,face = 'bold'),plot.title = element_text(hjust = 0.5),
        axis.ticks = element_blank())+ylim(limit=c(1,8))+
  annotate(geom = 'text',label=anno,fontface='bold',hjust=0,x=-1,y=7.5)
p6

#3.tnfsf12-Hallmark EMT
emt <- as.data.frame(t(read.csv('../hallmark-gsva.txt',row.names = 1,sep='\t',check.names = F)))
emt <- emt[sapply(samp,FUN=function(x){return(substring(x,1,12))}),'EMT markers',drop=F]
dt7 <- merge(emt,tnf,0)
test <- cor.test(dt7$`EMT markers`,dt7$TNFSF12,method = 'spearman',exact = F)
sr <- sprintf('%0.3f',test$estimate )
sp <- sprintf('%0.3e',test$p.value)
anno <- paste0("Spearman R=",sr,'\npVal=',sp)
p7<-ggscatter(dt7,x='EMT markers',y='TNFSF12',xlab = 'EMT')+labs(title = can)+
  stat_smooth(method = lm,color='red',fill='lightpink')+
  theme(text = element_text(size=12,face = 'bold'),plot.title = element_text(hjust = 0.5),
        axis.ticks = element_blank())+ylim(limit=c(1,8))+
  annotate(geom = 'text',label=anno,fontface='bold',hjust=0,x=-1,y=7.5)
p7
