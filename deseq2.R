args = commandArgs(T)
setwd("D://QQ.Li/DEG/20170703") #Need Change!!!!!!!!!!!!!!!!!
library(AnnotationDbi)
library(org.Mm.eg.db) #Need Change!!!!!!!!!!!!!!!!!
columns(org.Mm.eg.db)
Rawreads <- read.csv(paste("D://QQ.Li/DEG/20170703","/gene_count_matrix.csv",sep=""), row.names=1, stringsAsFactors = F)
colnames(Rawreads) <- sub("X", "", colnames(Rawreads))
Rawreads <- data.matrix(Rawreads)
pdf("Plotting.pdf", width=13, height=10)
boxplot(log2(Rawreads+1))





# Nomalization
counts <- Rawreads
  library(DESeq2)
  SampleInfo <- read.csv(paste("D://QQ.Li/DEG/20170703","/sample_information.csv",sep=""), stringsAsFactors=T, row.names=1)
  sample = data.frame(row.names=colnames(counts), condition=SampleInfo$group) #make a data.frame
  print(sample)
  # DESeqDataSet object constructors
  dds <- DESeqDataSetFromMatrix(countData = counts, 
                                colData = sample, 
                                design = ~ condition)
  # dds$condition <- relevel(dds$condition, ref="")
  # Differential expression analysis based on the Negative Binomial
  dds <- DESeq(dds)
  # holds the count data as a matrix of non-negative integer count values
  normcounts <- counts(dds, normalized=T)
  normcounts <- ceiling(normcounts) #向上取整数
  boxplot(log2(normcounts+1))
  
  
  
  
  
# Mapping Ensemble ID into EntrezID, GeneSymbol and geneName
  symbol = mapIds(org.Mm.eg.db,
                  keys=rownames(normcounts),
                  column="SYMBOL",
                  keytype="ENSEMBL",
                  multiVals="first")
  entrez = mapIds(org.Mm.eg.db,
                  keys=rownames(normcounts),
                  column="ENTREZID",
                  keytype="ENSEMBL",
                  multiVals="first")
  genename = mapIds(org.Mm.eg.db,
                    keys=rownames(normcounts),
                    column="GENENAME",
                    keytype="ENSEMBL",
                    multiVals="first")
  mapping=cbind(as.character(rownames(normcounts)),
                 as.character(symbol),
                 as.character(entrez),
                 as.character(genename))
  colnames(mapping) <- c("Ensemble", "GeneSymbol", "Entrez", "GeneName")
  rownames(mapping) <- mapping[, 1]
  mapping <- mapping[!is.na(mapping[, 3]), ] #去除未知基因名的行
  normcounts <- normcounts[rownames(mapping), ]
  write.csv(normcounts, file=paste("D://QQ.Li/DEG/20170703", "/normcounts.csv", sep=""))
  
  
  
  
  
# Filter out nonexpressed genes
  idx <- rep(0, length(normcounts[, 1]));
  for(i in 1:length(normcounts[, 1])){
    n=length(which(normcounts[i, ]>1)) #which 返回大于1的索引
    p=n/(dim(normcounts)[2]) #dim(normcounts)[2]，矩阵的第二维度：样本的总个数
    if(p >= 0.1) idx[i] <- 1 #阈值根据样品的多少来取舍
  }
  validcounts <- normcounts[idx==1, ]
  #write.csv(validcounts)
  boxplot(log2(validcounts+1), las=2)
  
  
  
  
  
# Calculate the scaled Pearson distance
  pearDist <- function(expreset){
  # expreset: expression matrix
    norm_expres <- as.matrix(expreset)
    for(i in dim(norm_expres)[1]){
        norm=mean(as.numeric(norm_expres[i, ]))
        for(j in 1:dim(norm_expres)[2]){
            norm_expres[i, j] <- norm_expres[i, j] - norm
        } 
    }
    dists <- as.dist(1-cor(norm_expres,method="pearson"))#as.dist()的作用是将普通矩阵转化为距离结构
    return(dists)                                        #cor函数，计算相关矩阵，method："spearman","pearson" and "kendall"
  }                                                      #两个样本越相近，相关系数越高，距离越小
# Cluster of all samples
    expression=log2(validcounts+1)
    distance <- pearDist(as.matrix(expression))
    hclust <- hclust(distance, method="average")
    #pdf("D://research/Platform/RNA-seq pipeline 1/data_1/hclust.pdf", width=13, height=10)
    plot(hclust, 
         hang=-1, 
         main=paste("clustering for",dim(expression)[2],"samples with ", dim(expression)[1], " genes"), 
         xlab="Samples", 
         lwd=2, 
         font=2, 
         cex=0.6)
# PCA
    pc<-prcomp(t(expression))
    names(pc)
    rownames(pc$x)
    colors=c(rep("red", 3), rep("yellow", 3)) #Need Change!!!!!!!!!!!!!!!!!
    plot(pc$x[, 1],pc$x[, 2], pch=19,col=colors,cex=2)
    legend("top",legend=c("DC", "SD"),cex=1.5,col=c("red", "yellow"),pch=19)
    dev.off()
    #library(rgl)
    #plot3d(pc$x[,1], pc$x[,2], pc$x[,3], col=colors, size=7)
  
  
#Differential expression
    F_calculate <- function(all_samples, sample_info, experiment, control, mapping){
        selected_samples <- sample_info[sample_info[, 1] %in% c(experiment, control), ] #从SampleInfo中提取实验组和对照组
        file_name <- rownames(sample_info)[sample_info[, 1] %in% c(experiment, control)] #从SampleInfo中提取实验组和对照组的文件名
        expreset <- all_samples[, file_name] #从validcounts中提取实验组和对照组

        design <- rep(0, dim(expreset)[2])
        design[which(as.character(selected_samples)==experiment)] <- 1    
        
        colData=data.frame(row.names=colnames(expreset), condition=factor(design))
        dds <- DESeqDataSetFromMatrix(countData = expreset, 
                                      colData = colData, 
                                      design = ~ condition)
        
        dds <- DESeq(dds) 
        result <- data.frame(results(dds))
        
        result=cbind(mapping[rownames(result), c(1, 2)], as.matrix(result))
        write.csv(result, file=paste("D://QQ.Li/DEG/20170703", "/final_results_of_", experiment,"_vs_",control, ".csv", sep=""))
        cat("Accompanished")
    }
    F_calculate(validcounts, SampleInfo, "SD", "DC", mapping) #Need Change!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    
    
