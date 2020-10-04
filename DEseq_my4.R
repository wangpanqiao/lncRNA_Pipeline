rm(list=ls())
library(DESeq2)
library(data.table)
library(dplyr)
# library(heatmaply)
library(pheatmap)
require(ggsci)
require(cowplot)
library(ggrepel)

STR=c("C","F1","H","S0","S13")
# combn(STR,2)
str_comb <- function(vector){
	n <- length(vector) 
	num=0  #保留所有组合的个数
	col=1 #用作循环叠加
	num=choose(n,2)+1 #计算组合个数
	comb_matrix <- matrix(,nrow = num-1,ncol = 2) #用矩阵保存组合结果
	res=combn(vector,2)
	m=ncol(res)
	for (L in 1:m){
		# name=res[,L]
		comb_matrix[col,] <- res[,L]#字符组合函数
		col=col+1
		if(col==num)break#当组合数量达到最终个数后，跳出循环
	}
	return(comb_matrix)
}
mycomb=str_comb(STR)

exprSet<-fread("lncRNA.rsem.count.txt", header=TRUE, sep="\t")#制作表达矩阵
des.count <- exprSet[exprSet$Type!="protein_coding",]
rna.id <- des.count[,1]
rna.type <- des.count[,2]
des.count <- as.data.frame(des.count[, round(.SD), .SDcols = -c(1:2), with = TRUE])
rownames(des.count) <- rna.id[[1]]
mycount=des.count
GROUP=read.csv("sample_file.txt",header=T, sep="\t") ##做好分组因子即可


plotvolcano=function(){
	dataset <- read.table(file = "DE.tsv", 
						header = TRUE,row.names=1, sep = "\t")
	# 设置p_value和logFC的阈值
	# cut_off_pvalue=sort(dataset$pvalue)[200]  ###倒序,decreasing=TRUE
	cut_off_pvalue = 0.0000001  #统计显著性
	cut_off_logFC = 2           #差异倍数值
	
	dataset$gene=rownames(dataset)
	
	# 根据阈值参数，上调基因设置为‘up’，下调基因设置为‘Down’，无差异设置为‘Stable’，并保存到change列中
	dataset$change = ifelse(dataset$pvalue < cut_off_pvalue & abs(dataset$log2FoldChange) >= cut_off_logFC, 
							ifelse(dataset$log2FoldChange> cut_off_logFC ,'Up','Down'),
							'Stable')
	p <- ggplot(
	# 数据、映射、颜色
	dataset, aes(x = log2FoldChange, y = -log10(pvalue), colour=change)) +
	geom_point(alpha=0.4, size=3.5) +
	scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757"))+
	# 辅助线
	geom_vline(xintercept=c(-2,2),lty=4,col="black",lwd=0.8) +
	geom_hline(yintercept = -log10(cut_off_pvalue),lty=4,col="black",lwd=0.8) +
	# 坐标轴
	labs(x="log2(fold change)",
		y="-log10 (p-value)")+
	theme_bw()+
	# 图例
	theme(plot.title = element_text(hjust = 0.5), 
			legend.position="right", 
			legend.title = element_blank())
	
	# 将需要标记的基因放置在label列(logFC >= 5)

	lable_value=sort(dataset$pvalue)[10] ###小到大,decreasing=TRUE大到小
	# print(lable_value)
	# lable_value=0.0000000001
	dataset$label <- ifelse(dataset$pvalue < lable_value & abs(dataset$log2FoldChange) >= 5,
						as.character(dataset$gene), "")
	
	p=p + geom_label_repel(data = dataset, aes(x = dataset$log2FoldChange, 
										y = -log10(dataset$pvalue), 
										label = label),
					size = 3, box.padding = unit(0.5, "lines"),
					point.padding = unit(0.8, "lines"), 
					segment.color = "black", 
					show.legend = FALSE)
	save_plot("volcanol.png", p,
			base_aspect_ratio = 1.3 # make room for figure legend
	)
}

runDE=function(count_vector,group_vector){
	dds <- suppressWarnings(DESeq2::DESeqDataSetFromMatrix(countData = count_vector,
								colData = group_vector,
								design = ~ group))
	keep <- rowSums(DESeq2::counts(dds)) >= 10
	dds <- dds[keep,]
	dds <- DESeq2::DESeq(dds)
	
	res <-  results(dds)# 提取你想要的差异分析结果，我们这里是C组对P组进行比较
	resOrdered_MY <- res[order(res$padj),]
	write.table(resOrdered_MY, file="DE.tsv", sep="\t",quote=F)
	
	resLFC <- DESeq2::lfcShrink(dds, coef=2)
	resOrdered <- resLFC[order(resLFC$pvalue),]
	vsd <- DESeq2::vst(dds, blind=FALSE)
	ntd <- DESeq2::normTransform(dds)
	
	
	### MA-plot
	# pdf('MARplot.pdf')
	png('MARplot.png',width=600*3,height=480*3,res=72*3)
	DESeq2::plotMA(resLFC, ylim=c(-2,2))
	dev.off()
	
	
	### Heatmap
	
	select <- order(rowMeans(DESeq2::counts(dds,normalized=TRUE)),
					decreasing=TRUE)[1:20]
	# df <- as.data.frame(SummarizedExperiment::colData(dds)[,"condition"])
	
	pheat1=pheatmap::pheatmap(SummarizedExperiment::assay(ntd)[select,], scale = 'row', xlab = "Sample", margins = c(60,100,40,20), row_text_angle = 45, column_text_angle = 60)
	# heatmaply::heatmaply(SummarizedExperiment::assay(ntd)[select,], scale = 'row', xlab = "Sample", margins = c(60,100,40,20), row_text_angle = 45, column_text_angle = 60) %>% layout(margin = list(l = 100, b = 50, r = 0))
	
	save_plot('Heatmap.png', pheat1, base_height = 8.5,base_width = 11, base_aspect_ratio = 1.3)  ###base_width = 11,
	
	### Correlation heatmap
	
	pheat2=pheatmap::pheatmap(cor(SummarizedExperiment::assay(ntd)[select,]), margins = c(40, 40),
						k_col = 2, k_row = 2,
						limits = c(-1,1))
	# heatmaply::heatmaply(cor(SummarizedExperiment::assay(ntd)[select,]), margins = c(40, 40),
			# k_col = 2, k_row = 2,
			# limits = c(-1,1)) %>% layout(margin = list(l = 50, b = 50))
	
	save_plot('Correlation_heatmap.png', pheat2, base_width = 11, base_height = 8.5,base_aspect_ratio = 1.3)  ###base_height = 8.5, 
	
	### Principal Component Analysis
	
	pcaData <- DESeq2::plotPCA(vsd, intgroup="group", returnData=TRUE)
	percentVar <- round(100 * attr(pcaData, "percentVar"))
	p <- ggplot(pcaData, aes(PC1, PC2, color=group)) +
	geom_point(size=3) +
	xlab(paste0("PC1: ",percentVar[1],"% variance")) +
	ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
	coord_fixed()
	save_plot('pca.tiff', p, base_height = 8.5, base_width = 11, dpi = 300, compression = 'lzw')
	save_plot('pca.pdf', p, base_height = 8.5, base_width = 11, dpi = 300)
	# ggplotly(p)
	
	### Description
	
	### Differentially expressed lncRNAs table
	
	resSig <- BiocGenerics::subset(resOrdered, padj < 0.1)
	# fwrite(as.data.table(resSig), 'DE.csv', row.names = TRUE)
	id=rownames(resSig)
	resSig1=cbind(id,resSig)
	fwrite(as.data.table(resSig1), 'DE.csv', row.names = FALSE,quote=FALSE,sep="\t",col.names=TRUE)
	# DT::datatable(head(as.data.frame(resSig), n = 80L)) %>% DT::formatRound(c('baseMean', 'log2FoldChange', 'lfcSE', "stat", 'pvalue', 'padj'), digits = 2)
	#添加火山图
	plotvolcano()
}


# setwd("/home/ysq/wpq/mydocuments/newspecies_chystivus/ISOseq/lncRNA/deg_reporter5/manual_my")
mydir=getwd()
# 提取sub
for (i in 1:dim(mycomb)[1]){
	grepstr=paste(c(as.character(paste(c(mycomb[i,1],mycomb[i,2]),collapse="_vs_")),"DE"),collapse="__")
	# print(grepstr)
	subgroup=GROUP[GROUP$group==mycomb[i,1] | GROUP$group==mycomb[i,2],]
	# print(subgroup[[1]])
	subcount=mycount[,c(subgroup[[1]])]
	# print(head(subcount))
	coldata=data.frame(row.names=subgroup$Sample, group=subgroup$group)
	# print(coldata)
	mainDir <- mydir
	subDir <- grepstr
	if (dir.exists(file.path(mainDir, subDir))){
		setwd(file.path(mainDir, subDir))
	} else {
		dir.create(file.path(mainDir, subDir))
		setwd(file.path(mainDir, subDir))
	}
	write.table(subcount,paste(c(paste(c(mycomb[i,1],mycomb[i,2]),collapse="_vs_"),".count"),collapse=""),quote=F,row.names=T,col.names=T,sep="\t")
	write.table(coldata,paste(c(paste(c(mycomb[i,1],mycomb[i,2]),collapse="_vs_"),".group"),collapse=""),quote=F,row.names=T,col.names=T,sep="\t")
	runDE(subcount,coldata)
}

