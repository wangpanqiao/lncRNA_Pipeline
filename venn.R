rm(list=ls())
library(VennDiagram)
data=read.table('novel.longRNA.txt',header=T,row.names=1,sep="\t")
head(data)
CPAT1=data[data$CPAT=="lncRNA",]
PLEK1=data[data$PLEK=="lncRNA",]
CNCI1=data[data$CNCI=="lncRNA",]
CPC21=data[data$CPC2=="lncRNA",]
newcol=c("cornflowerblue","green","yellow","darkorchid1")
venn.diagram(list(CPAT=row.names(CPAT1), PLEK=row.names(PLEK1),
                 CNCI=row.names(CNCI1), CPC2=row.names(CPC21)),col = "transparent",  ##去掉外线
             # fill=c("red","green","blue","yellow"), cat.col=c("red","green","blue","yellow"),
             fill=newcol, ### cat.col=newcol,
             alpha=c(0.5,0.5,0.5,0.5),cat.fontfamily="serif",,  ###字体  2黑体  
             imagetype = "png", category.names = c("CPAT","CPC2","PLEK","CNCI"),
             height = 600, width = 600, resolution = 100,
             cex=1.6, cat.cex=1.8, filename="VennDiagram.png")  ##cat.fontface=4斜体, 2加粗,cat.fontface= "bold"



hha=Reduce(intersect,list(CPAT=row.names(CPAT1), PLEK=row.names(PLEK1),
                 CNCI=row.names(CNCI1), CPC2=row.names(CPC21)))
# write.table(hha,"intersect.lnc.txt")# 带行名和双引号

write.table(hha,"intersect.lnc.txt",row.names=F,quote=F)


# venn.diagram(list(WDSP=wdspWD40,Pfam=pfamWD40,SMART=smartWD40),
    # resolution = 300, imagetype = "tiff", alpha=c(0.5,0.5,0.5),
    # fill=c("red","yellow","blue"), cat.fontface=4,fontfamily=3,
    # main="Intersection of WD40 genes identified by different methods",
    # main.cex = 2, main.fontface = 2, main.fontfamily = 3,
    # filename = "VennDiagram.tif")


# > venn.plot <- venn.diagram(
  # x = list(
    # X124 = data$X124,
    # C88 = data$C88,
    # H6 = data$H6,
    # SI = data$SI
  # ),
  # filename = "1-venn.tiff",
  # col = "black",
  # lty = "solid", #边框线型改为"dotted"虚线
  # lwd = 2, # 边框线的宽度
  # fill = c("cornflowerblue", "green", "yellow", "darkorchid1"),#每个圈的颜色
  # alpha = 0.50, #颜色的相对深度
  # label.col = c("black"),#填充的字体颜色
  # cex = 1.5,#填充字体大小
  # fontfamily = "serif",#填充字体
  # fontface = "bold",
  # cat.col = c("darkblue", "darkgreen", "orange", "darkorchid4"),#每个圈标注文字颜色，darkblue
  # cat.cex = 1.5,#标注文字大小
  # cat.fontface = "bold",
  # cat.fontfamily = "serif"
# )
























