library(gplots) # dowload package
Genesmap=read.table("C:\\Users\\liqin\\Desktop\\20171208\\old\\MT_vs_PB\\final_results_of_MT_vs_PB.csv",header=T,sep=",")
names(Genesmap)
a <- Genesmap[intersect(which(Genesmap$log2FoldChange > log2(1.5)) ,which(Genesmap$padj < 0.05)),0:8]
write.csv(a,file="C:\\Users\\liqin\\Desktop\\20171208\\old\\MT_vs_PB\\upregulated.csv")

a <- Genesmap[intersect(which(Genesmap$log2FoldChange < -log2(1.5)) ,which(Genesmap$padj < 0.05)),0:8]
write.csv(a,file="C:\\Users\\liqin\\Desktop\\20171208\\old\\MT_vs_PB\\downregulated.csv")
