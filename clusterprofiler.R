library(clusterProfiler)
setwd("C:\\Users\\liqin\\Desktop\\ZSQ\\clusterprofiler")
a=read.table("up.txt",header = FALSE)

require(DOSE)
require(clusterProfiler)
gene=as.character(a[,1])

ego <- enrichGO(gene=gene,'org.Mm.eg.db',ont="BP",pvalueCutoff=0.01,readable=TRUE)
write.csv(summary(ego),"GO-enrich_up.csv",row.names =F)

barplot(ego, showCategory=15,title="EnrichmentGO_up")