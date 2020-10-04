rm(list=ls())
sink("lnc_length.txt")
library(ggplot2)
library(cowplot)
library(data.table)
require(ggsci)
theme = 'npg'

toUpperFirstLetter <- function(x) {
  paste0(toupper(substr(x, 1, 1)), substring(x, 2))
}
lncRNA <- fread("basic_charac.txt", sep="\t")
setnames(lncRNA, c('Gene', 'Transcript', 'Type', 'Potential', 'Length', 'exon_num', 'Class'))

# lncRNA[, Type := toUpperFirstLetter(sub("_", " ", Type))]

lncRNA[, Type := sub("novel", "lncRNA", Type)]
lncRNA[, Type := sub("protein_coding", "mRNA", Type)]

mRNA=lncRNA[Type=="mRNA"]
summary(mRNA$Length)
fivenum(mRNA$Length)
quantile(mRNA$Length)
length(mRNA$Length)
a <- hist(mRNA$Length, breaks = 20)
a


lncRNA1=lncRNA[Type=="lncRNA"]
summary(lncRNA1$Length)
fivenum(lncRNA1$Length)
quantile(lncRNA1$Length)

length(lncRNA1$Length)
b <- hist(lncRNA1$Length, breaks = 20, freq = F)
b




# lncRNA.gtf <- unique(lncRNA[, Length := sum(Length), by = .(Gene, Type)], by = 'Gene')   ##base

# lncRNA.gtf <- unique(lncRNA[, exon_num := sum(exon_num), by = .(Gene, Type)], by = 'Gene')   ##base

lncRNA.gtf=lncRNA
# mytext=lncRNA[carrier == "AA", lapply(.SD, mean), 
          # by = .(origin, dest, month),
          # .SDcols = c("arr_delay", "dep_delay")]

mytext=lncRNA[,lapply(.SD, mean), by = .(Type),.SDcols = c("Length")]

# print(mytext)
headtext="average length: "
lable1=paste(mytext[1,1],headtext,round(mytext[1,2],2),"bp",sep=" ")
# print(lable1)
lable2=paste(mytext[2,1],headtext,round(mytext[2,2],2),"bp",sep=" ")
# print(lable2)
# mylabel=as.character(as.expression(mytext))
mylable=paste(lable1,lable2,sep="\n\n")
print(mylable)


max.lncrna.len = 10000
# min.expressed.sample = 50
# lncRNA.gtf[Length > params$max.lncrna.len, Length := params$max.lncrna.len]
# p <- ggplot() + geom_density(data = lncRNA.gtf, aes(x = Length, colour = Type), size = 1.5) +
  # xlab('Transcript length') + ylab('Density') +
  # scale_x_continuous(breaks = seq.int(200, params$max.lncrna.len, length.out = 10),
                     # labels = c(seq.int(200, params$max.lncrna.len, length.out = 9),
                                # paste0(params$max.lncrna.len, '+')), 
                     # expand = c(0.01, 0)) +
  # scale_y_continuous(expand = c(0.01, 0)) +
  # get(paste0('scale_color_',params$theme))()
lncRNA.gtf[Length > max.lncrna.len, Length := max.lncrna.len]
p <- ggplot() + geom_density(data = lncRNA.gtf, aes(x = Length, y=..density.., colour = Type), size = 1.5,adjust=1.2) + ###adjust平滑度,size线粗细
  xlab('Transcript length') + ylab('Density') +
  scale_x_continuous(breaks = seq.int(200, max.lncrna.len, length.out = 10),
                     labels = c(seq.int(200, max.lncrna.len, length.out = 9),
                                paste(max.lncrna.len, '+ bp')), 
                     expand = c(0.01, 0)) +
  scale_y_continuous(expand = c(0.01, 0)) +
  # geom_text(x=20, y=26, aes(label=mylable), parse=T) +
  annotate("text",x=5000,y=0.00075,label=mylable,size=5)+
  # annotate("text",x=9,y=0.25,label=lable1)+  ###输入公式,parse=T,novel,短
  # annotate("text",x=12,y=0.2,label=lable2)+
  get(paste0('scale_color_',theme))() +
  theme(legend.position = c(0.85,0.85)) +
  theme(legend.title = element_text(colour="black", size=20)) +
  theme(legend.text = element_text(colour="black", size = 16)) +
  theme(axis.title.x=element_text(size=25),axis.title.y=element_text(size=25))
save_plot('length_type.png', p, base_width = 11, base_height = 8.5,base_aspect_ratio = 1.3)
save_plot('length_type.tiff', p, base_height = 8.5, base_width = 11, dpi = 300, compression = 'lzw')
# save_plot('lncRNA_exon_num_distribution_with_type.pdf', p, base_height = 8.5, base_width = 11, dpi = 300)