library(gplots) # dowload package
Genesmap=read.table("D:\\QQ.Li\\BIT\\Heatmap\\ISG-tumor-ZSQ.csv",header=T,sep=",")
names(Genesmap)

a <- Genesmap[2:15]

b <- as.matrix(a)

heatmap.2(b,Colv=F,Rowv=F,scale="row",symbreaks = T,col=bluered(255),colsep=1:20,rowsep=1:20,trace="none",labRow=Genesmap$NAME,cexRow=.9,cexCol=.9,dendrogram="none",sepwidth=c(0.02,0.02))

