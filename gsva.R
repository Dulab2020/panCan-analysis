setwd('~/research/2020single/finalplot/plot11-L-R-supp')
library(GSVA)
library(BiocParallel)
library(GSEABase)
library(data.table)

df <- fread('alllog2tpmp1-genename.csv',sep=',')
mat <- as.matrix(df[,-1])
rownames(mat)<-df$Gene_name

#saveRDS(mat,'alllog2tpmp1-mat.RDS')

gmtf <- getGmt('exh-pro.gmt',sep='\t',geneIdType=SymbolIdentifier())

gsvaOut <- gsva(mat,gmtf,verbose=T,parallel.sz=1,min.sz=1,max.sz=Inf,
                BPPARAM=MulticoreParam(workers=20,progressbar=T))


gsvaOut1 <- as.data.frame(t(gsvaOut))

write.csv(gsvaOut1,'all-exh-pro.gsva.csv')
