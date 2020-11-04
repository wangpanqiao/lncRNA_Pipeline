rm(list = ls())
# library(data.table)
options(stringsAsFactors = F)

args=commandArgs(T)

# fileroot=strsplit(basename(args[1]),split=".")[[1]][1]

# print(fileroot)
a = read.table(args[1], header = T,sep="\t")

a=a[which(a$Type=="novel"),]
a$Type[which(a$Type=="novel")]="lncRNA"
# head(a)
# a=a[,-c(S0_1,S0_2,S0_3)]

a=subset(a,select=-c(S0_1,S0_2,S0_3))

# dat=subset(des.count, ave(rowSums(des.count[-1:-2]) > 0, Category, FUN = any))
a1=subset(a, ave(rowSums(a[-1:-2]) > 0.3, ID, FUN = any))

# dim(a)
dat = a1[,-c(1,2)]
library(stringr)
ac = data.frame(group=str_split(colnames(dat),'_',simplify = T)[,1])
rownames(ac) = colnames(dat)
M=cor(log(dat+1))

p=pheatmap::pheatmap(M,
                   annotation_col = ac,display_numbers = matrix(ifelse(M > 0.92, round(M,2), ""), nrow(M)),
                   number_color = "blue", number_size= 2
                   # breaks = seq(0, 100, length.out = 100)
                   ) 



mRNA=a[which(a$Type != "novel"),]

# dat=subset(des.count, ave(rowSums(des.count[-1:-2]) > 0, Category, FUN = any))
mRNA1=subset(mRNA, ave(rowSums(mRNA[-1:-2]) > 0.3, ID, FUN = any))

datmRNA = mRNA1[,-c(1,2)]
library(stringr)
ac1 = data.frame(group=str_split(colnames(datmRNA),'_',simplify = T)[,1])
rownames(ac1) = colnames(datmRNA)
M1=cor(log(datmRNA+1))

mRNAp=pheatmap::pheatmap(M1,
                   annotation_col = ac,display_numbers = matrix(ifelse(M1 > 0.92, round(M1,2), ""), nrow(M1)),
                   number_color = "blue", number_size= 2
				   # breaks = seq(0, 100, length.out = 100)
                   ) 

# 这个包使用的是grid图形系统而非ggplot2，所以解决方法也是不同的。通过自定义函数来生成，也可一次绘制多个对象的图形。

# save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
   # stopifnot(!missing(x))
   # stopifnot(!missing(filename))
   # pdf(filename, width=width, height=height)
   # grid::grid.newpage()
   # grid::grid.draw(x$gtable)
   # dev.off()
# }

save_pheatmap_png <- function(x, filename, width=700, height=700) {
   stopifnot(!missing(x))
   stopifnot(!missing(filename))
   png(filename, width=width, height=height)
   print(x)
   dev.off()
}

save_pheatmap_tif <- function(x, filename, width=7, height=7) {
   stopifnot(!missing(x))
   stopifnot(!missing(filename))
   # filename=paste(x,'tif',sep='.'))
   tiff(filename, compression="lzw",units="in",res=300,pointsize=8, width=width, height=height)
   print(x)
   dev.off()
}


save_pheatmap_png(p, "lnc.png")
save_pheatmap_tif(p, "lnc.tif")
save_pheatmap_png(mRNAp, "mRNA.png")
save_pheatmap_tif(mRNAp, "mRNA.tif")