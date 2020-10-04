rm(list=ls())
sink("exon_hist.txt")
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

# lncRNA[, Type := toUpperFirstLetter(sub("novel", "lncRNA", Type))]
# lncRNA[, Type := toUpperFirstLetter(sub("protein_coding", "mRNA", Type))]

lncRNA[, Type := sub("novel", "lncRNA", Type)]
lncRNA[, Type := sub("protein_coding", "mRNA", Type)]

mRNA=lncRNA[Type=="mRNA"]
head(mRNA)
print(summary(mRNA$exon_num))
print(fivenum(mRNA$exon_num))
print(quantile(mRNA$exon_num))

#freq <- table(data)
freq1 <- table(mRNA$exon_num)
table2 <- prop.table(freq1)
print(freq1)
print(table2)

lncRNA1=lncRNA[Type=="lncRNA"]
print(summary(lncRNA1$exon_num))
print(fivenum(lncRNA1$exon_num))
print(quantile(lncRNA1$exon_num))
freq2 <- table(lncRNA1$exon_num)

print(freq2)
table3 <- prop.table(freq2)

print(table3)

# lncRNA.gtf <- unique(lncRNA[, Length := sum(Length), by = .(Gene, Type)], by = 'Gene')   ##base

# lncRNA.gtf <- unique(lncRNA[, exon_num := sum(exon_num), by = .(Gene, Type)], by = 'Gene')   ##base

lncRNA.gtf=lncRNA
# mytext=lncRNA[carrier == "AA", lapply(.SD, mean), 
          # by = .(origin, dest, month),
          # .SDcols = c("arr_delay", "dep_delay")]

mytext=lncRNA[,lapply(.SD, mean), by = .(Type),.SDcols = c("exon_num")]

# print(mytext)
headtext="average number: "
# headtext="mean: "
lable1=paste(mytext[1,1],headtext,round(mytext[1,2],2),sep=" ")
# print(lable1)
lable2=paste(mytext[2,1],headtext,round(mytext[2,2],2),sep=" ")
# print(lable2)
# mylabel=as.character(as.expression(mytext))
mylable=paste(lable1,lable2,sep="\n\n")
print(mylable)


max.lncrna.len = 20
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
# lncRNA.gtf[exon_num > max.lncrna.len, exon_num := max.lncrna.len]
# lncRNA.gtf[exon_num > max.lncrna.len, exon_num := max.lncrna.len]
lncRNA[exon_num > max.lncrna.len, exon_num := max.lncrna.len]
# lnc=lncRNA[Type=='lncRNA',]
# mRNA=lncRNA[Type=='mRNA',]
# seq.int(1, 22, by=2 )
# seq.int(1, 20, by=2 )
p <- ggplot(data = lncRNA, aes(x = exon_num, fill = Type)) + 
  # geom_histogram(aes(y=(..count..)/sum(..count..)),position = "dodge", bins=20 ,alpha = 0.4) +  ###,binwidth=0.5,总比
  geom_histogram(aes(y=..density..),position = "identity", bins=20 ,alpha = 0.5) + ###组内百分比
  xlab('The number of exons') + ylab('Proportion') +
  scale_x_continuous(breaks = seq.int(1, max.lncrna.len+2, by=2 ), ####,length.out = 20   个数 by=1
                     labels = c(seq.int(1, max.lncrna.len, by=2 ),  ###,length.out = 19
                                paste0(max.lncrna.len, '+')), 
                     expand = c(0.01, 0)) +
  scale_y_continuous(expand = c(0.01, 0)) +
  # geom_text(x=20, y=26, aes(label=mylable), parse=T) +
  annotate("text",x=9,y=0.25,label=mylable,size=5)+
  # annotate("text",x=9,y=0.25,label=lable1)+  ###输入公式,parse=T,novel,短
  # annotate("text",x=12,y=0.2,label=lable2)+
  get(paste0('scale_color_',theme))()+
  theme(legend.position = c(0.85,0.85)) +
  theme(legend.title = element_text(colour="black", size=20)) +
  theme(legend.text = element_text(colour="black", size = 16)) +
  theme(axis.title.x=element_text(size=25),axis.title.y=element_text(size=25))
save_plot('exon_num_distribution_with_type1.png', p, base_width = 11, base_height = 8.5,base_aspect_ratio = 1.3)
save_plot('lncRNA_exon_num_distribution_with_type1.tiff', p, base_height = 8.5, base_width = 11, dpi = 300, compression = 'lzw')
# save_plot('lncRNA_exon_num_distribution_with_type1.pdf', p, base_height = 8.5, base_width = 11, dpi = 300)