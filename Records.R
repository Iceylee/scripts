Rscript /data1/script/deseq2+GO+KEGG/Rpipe/3GO-KEGG.R AH57973 hsa SYMBOL ENTREZID /data2/ClientData/2018_08/Lishuang/

##1 2.1
setwd("/Users/Icey/work/2018/plot/0823热图-谢华平/")
coldata_file="colData.csv"
count_table="CountMatrix4DESeq.csv"
path1="./"
path3="./"
path2="./"

Rscript /data1/script/deseq2+GO+KEGG/Rpipe/run.R CountMatrix4DESeq_SvsF.csv colData_SvsF.csv LZY_Con /data2/ClientData/2018_08/XieHuaPing/ TRUE

##3GO-KEGG.R
dbname <- 'AH57973'
kegg_org <- "hsa"
GO_KEY <- "SYMBOL" 
KEGG_NEED_KEY <- "ENTREZID" 
output_path <- "/data2/ClientData/2018_08/XuLeLe/"


##pathview.R
#跑前删除pathway文件夹。不然有可能报错:
# Error in readPNG(paste(kegg.dir, "/", pathway.name, ".png", sep = "")) :
#   libpng error: Read Error

#模式生物：sig files 中有log2FoldChange列
#非模式生物：第四列为log2FoldChange的数值。手动加列名

Rscript /data1/script/deseq2+GO+KEGG/Rpipe/pathview.R ./nonModel/3.DiffExprGene/ ./nonModel/4.GO_KEGG_Enrichment/ ENSEMBL FALSE  

Rscript /data1/script/deseq2+GO+KEGG/Rpipe/pathview.R ./Model/3.DiffExprGene/ ./Model/4.GO_KEGG_Enrichment/0.05/ ENSEMBL TRUE 


#4.GetGene2KO.R
##得到物种的gene symbol对应KO编号。（大K）
Rscript /data1/script/deseq2+GO+KEGG/Rpipe/4GetGene2KO.R AH57973 hsa ./
Rscript /data1/script/deseq2+GO+KEGG/Rpipe/4GetGene2KO.R AH59553 vda ./


##length_hist.R
Rscript /data1/script/deseq2+GO+KEGG/Rpipe/length_hist.R 1.Trinity/CD-HIT/Trinity_CD-HIT_0.9.seqLengths 3000

#5.gfold R 并入
Rscript /data1/script/deseq2+GO+KEGG/Rpipe/run.R R_input/CountMatrix4DESeq.csv R_input/colData.csv SH1 /data2/ClientData/2018_08/XuLeLe/Gfold_test/ FALSE

#5DAG.R
#XuLeLe 人 模式生物
Rscript /data1/script/deseq2+GO+KEGG/Rpipe/5DAG.R ./4.GO_KEGG_Enrichment/AllGene_GO_Annotation.txt ./3.DiffExprGene/ ./4.GO_KEGG_Enrichment TRUE

#run.R
#baseGroup 参数修改 ：直接写比较组别。多组逗号隔开
Rscript /data1/script/deseq2+GO+KEGG/Rpipe/run.R CountMatrix4DESeq.csv colData.csv  NC_vs_V ./ TRUE

#PlotEnrichHeatmapVenn.sh (chipseq 画富集热图和维恩图)
bash /data1/script/deseq2+GO+KEGG/Rpipe/PlotEnrichHeatmapVenn.sh ./ A1 A2
#在call peak路径（第一个参数）下面新建一个目录CommDiffPeaks，所有文件输出到该文件夹
