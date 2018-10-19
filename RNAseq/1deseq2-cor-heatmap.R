

########################################################
#####################1.deseq2-normalize#################
########################################################
# Load the library.
suppressMessages(library(DESeq2))
suppressMessages(library(dplyr))

countData = read.table(count_table, header=TRUE, sep=",", row.names=1)
#countData = countData[,-1] #delete first column(entrezID)
colData = read.csv(coldata_file, header=T,row.names = 1)
countData <- countData[,rownames(colData)]
colData$condition = as.factor(colData$condition)
#colnames(countData) <- colData[,1]
if (all(rownames(colData) %in% colnames(countData)) == F) {
  print ("rownames(colData) != colnames(countData)")
  print ("Please check.")
  stop(call.=FALSE)
}

# Create DESEq2 dataset.
dds = DESeqDataSetFromMatrix(countData=countData, colData=colData, design = ~condition) #~在R里面用于构建公式对象，~左边为因变量，右边为自变量。

#pre filter
dds = dds[ rowSums(counts(dds)) > 1 ,] 
#Set the reference to be compared
dds$condition <- relevel(dds$condition,base_group)

# Run deseq
dds <- DESeq(dds) #标准化
# Get normalized counts and write this to a file
nc = counts(dds,normalized=TRUE)
# Turn it into a dataframe to have proper column names.
dt = data.frame("id"=rownames(nc),nc)
colnames(dt) <- c("id",colnames(nc))
# Save the normalize data matrix.
write.table(dt, file=paste(path1,"norm-count-matrix.txt",sep="/"), sep="\t",  row.name=FALSE, col.names=TRUE,quote=FALSE)

save(dds, file=paste(path1,"1Data_Input.RData",sep="/"))

########################################################
#####################2.correlation######################
########################################################
df <- read.table(file = paste(path1,"norm-count-matrix.txt",sep="/"),sep = '\t',header = T,stringsAsFactors=F,row.names = 1)
colnames(df) <- rownames(colData)
cor.df <- cor(df[,unlist(lapply(df, is.numeric))])

anno = data.frame(Condition = colData$condition)
anno$Condition = as.factor(anno$Condition)
rownames(anno) = rownames(colData)


library(pheatmap)
pdf(file=paste(path3,"heatmap_cor.pdf",sep="/"),onefile=FALSE)
pheatmap(cor.df,legend = TRUE,display_numbers = TRUE, number_format = "%.3f", number_color="black",cluster_cols = FALSE,cluster_rows = FALSE,annotation_col=anno,annotation_row=anno)
dev.off()

png(file=paste(path3,"heatmap_cor.png",sep="/"))
pheatmap(cor.df,legend = TRUE,display_numbers = TRUE, number_format = "%.3f", number_color="black",cluster_cols = FALSE,cluster_rows = FALSE,annotation_col=anno,annotation_row=anno)
dev.off()

########################################################
#####################3.heatmap##########################
########################################################
rld <- rlog(dds, blind = FALSE)
#######前20样本差异较大的，做聚类热图。不同样本中的表达情况
library("genefilter") 
library("pheatmap")
topVarGenes <- head(order(rowVars(assay(rld)), decreasing = TRUE), 50)
mat  <- assay(rld)[ topVarGenes, ] 
mat  <- mat - rowMeans(mat) 
anno = data.frame(Condition = colData$condition)
anno$Condition = as.factor(anno$Condition)
rownames(anno) = rownames(colData)

#long gene name
# pdf(file=paste(path1,"heatmap_top50.pdf",sep="/"),onefile=FALSE)
# pheatmap(mat,cellwidth=15,cellheight=5,fontsize=5,annotation_col = anno,cluster_cols = FALSE)
# dev.off()
# png(file=paste(path1,"heatmap_top50.png",sep="/"))
# pheatmap(mat,cellwidth=15,cellheight=5,fontsize=5,annotation_col = anno,cluster_cols = FALSE)
# dev.off()

#short gene name
pdf(file=paste(path1,"heatmap_top50.pdf",sep="/"),onefile=FALSE)
pheatmap(mat,annotation_col = anno,cluster_cols = FALSE)
dev.off()
png(file=paste(path1,"heatmap_top50.png",sep="/"))
pheatmap(mat,annotation_col = anno,cluster_cols = FALSE)
dev.off()


