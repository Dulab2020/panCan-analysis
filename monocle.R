library(monocle)
library(patchwork)
rds <- dir(pattern='*M-data.rds')

for (f in rds){
    cat(f)
    cat('\n')
    outname <- paste0(strsplit(f,'.',fixed=T)[[1]][1],'.pdf')
    df <- readRDS(f)
    pd<-readRDS(paste0(strsplit(f,'-',fixed=T)[[1]][1],'-pd.rds'))
    pd=new("AnnotatedDataFrame",pd)
    fd=new("AnnotatedDataFrame",df$fd)
    HSMM <- newCellDataSet(df$exp,
                phenoData = pd,
                featureData = fd,
                lowerDetectionLimit = 0,
                expressionFamily = negbinomial.size())
    HSMM <- estimateSizeFactors(HSMM)
    HSMM <- estimateDispersions(HSMM)
    HSMM_myo <- reduceDimension(HSMM, max_components = 2,method = 'DDRTree')
    HSMM_myo <- orderCells(HSMM_myo)
    script <- paste0(strsplit(f,'-')[[1]][1],'_HSMM_myo<-HSMM_myo')
    eval(parse(text=script))
    #pdf(outname)
    a<-plot_cell_trajectory(HSMM_myo, color_by = "Pseudotime",cell_size=0.5,show_branch_points=F)+
        scale_color_viridis_c(direction = -1,option = 'B',alpha = 0.8)+
        theme(legend.position='none',axis.title=element_blank(),
        axis.text=element_blank(),axis.ticks=element_blank())

    b<-plot_cell_trajectory(HSMM_myo, color_by = "cluster",cell_size=0.5,show_branch_points=F)+
    theme(axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),
    legend.position=c(0,1),legend.justification=c(0,1),legend.title=element_blank(),
    legend.background=element_blank(),legend.direction='vertical',legend.text=element_text(size=12))+
    guides(color=guide_legend(override.aes=list(size=2)))
    print(a+b)
    script <- paste0(strsplit(f,'-')[[1]][1],'<-a+b')
    eval(parse(text=script))
    #print(plot_cell_trajectory(HSMM_myo, color_by = "cluster")+facet_wrap(~cluster,nrow=2))
    #print(plot_cell_trajectory(HSMM_myo, color_by = "SPP1")+scale_color_viridis_c(direction = -1,option = 'B',alpha = 0.8))
    #print(plot_cell_trajectory(HSMM_myo, color_by = "Hypoxia")+scale_color_viridis_c(direction = -1,option = 'B',alpha = 0.8))
    #dev.off()
}

for (i in c('BC','CRC','LC','OV','PDAC','SCC')){
    script <- paste0('HSMM_myo <- ',i,'_HSMM_myo')
    HSMM_myo <
    a<-plot_cell_trajectory(HSMM_myo, color_by = "Pseudotime",cell_size=0.5,show_branch_points=F)+
        scale_color_viridis_c(direction = -1,option = 'B',alpha = 0.8)+
        theme(legend.position='none',axis.title=element_blank(),
        axis.text=element_blank(),axis.ticks=element_blank())

    b<-plot_cell_trajectory(HSMM_myo, color_by = "cluster",cell_size=0.5,show_branch_points=F)+
    theme(axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),
    legend.position=c(0,0),legend.justification=c(0,0),legend.title=element_blank(),
    legend.background=element_blank(),legend.direction='horizontal',legend.text=element_text(size=12))+
    guides(color=guide_legend(override.aes=list(size=2)))
    print(a+b)
    script <- paste0(strsplit(f,'-')[[1]][1],'<-a+b')
    eval(parse(text=script))
}






