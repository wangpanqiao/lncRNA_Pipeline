library(ggplot2)
library(cowplot)
library(data.table)
# library(ggpubr)

toUpperFirstLetter <- function(x) {
  paste0(toupper(substr(x, 1, 1)), substring(x, 2))
}
lncRNA <- fread("basic_charac.txt", sep="\t")
setnames(lncRNA, c('Gene', 'Transcript', 'Type', 'Potential', 'Length', 'exon_num', 'Class'))

lncRNA <- lncRNA[lncRNA$Type!="protein_coding",]

lncRNA[, Type := toUpperFirstLetter(sub("_", " ", Type))]

cls <- lncRNA[, .(count = .N), by = Class]
colors <- c('rgb(211,94,96)', 'rgb(128,133,133)', 'rgb(144,103,167)', 'rgb(171,104,87)', 'rgb(114,147,203)')

cls1 = cls[order(cls$count, decreasing = TRUE),]

write.table(as.data.frame(cls1),'classification.tsv',quote=F,sep="\t",row.names=F,col.names=T)

print(cls1)
# myLabel = as.vector(cls1$Class)   
# myLabel = paste(myLabel, "(", round(cls1$count / sum(cls1$count) * 100, 2), "%)", sep = "")  

label_value <- paste('(', round(cls1$count/sum(cls1$count) * 100, 1), '%)', sep = '')
#label <- paste(cls1$Class, label_value, sep = '')
# cls1$labelPosition<-(df$ymax + df$ymin)/2

###饼状



blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )


p1=ggplot(data=cls1, mapping=aes(x="", y=count, fill=Class)) +
  geom_bar(stat="identity",width=0.7) + coord_polar(theta = 'y', start = 0)  + 
  blank_theme + 
  labs(x = "", y = "", title = "Percentage distribution on lncRNA classification (%)")  +
  geom_text(aes(y = cls1$count/2 + c(0, cumsum(cls1$count)[-length(cls1$count)]),x=sum(cls1$count)/2400, label = label_value),size=5) ## +
  #theme(legend.position = "none")  + 
  #geom_text(aes(y = cls1$count/2 + c(0, cumsum(cls1$count)[-length(cls1$count)]) , label = label),size=5)
  #geom_text(aes(y = cls1$count/2 + c(0, cumsum(cls1$count)[-length(cls1$count)]) , x = sum(cls1$count)/3000, label = label),size=5) 




png('test.png')
p1
dev.off()

##柱形图
p2=ggplot(cls,aes(x=Class,y=count))
myFont1="Times"
# myFont1="Arial"
ex1=expression(bold(paste("lncRNA Type")))
ex2=expression(bold(paste(PM["2.5"])))
ex3=expression(bold(paste(PM["10"])))
p2=p2+geom_bar(aes(fill=Class),stat='identity',position = 'dodge',width=.5)+
  xlab(ex1)+ylab("Number")+
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5,size = 10, color = "black", face = "bold",family="myFont1"))+
  theme(axis.title.x = element_text(size = 10, color = "black", face = "bold",family="myFont1"))+
  theme(axis.text.y = element_text(size = 10, color = "black", face = "bold",family="myFont1"))+
  theme(axis.title.y = element_text(size = 10, color = "black", face = "italic",family="myFont1"))+
  theme(legend.title=element_blank())+
  theme(legend.text=element_text(face="bold", family="myFont1", colour="black",size=8))+
  theme(legend.position=c(0.1,0.9),legend.background=element_blank())



# save_plot('classification1.png', p1, base_width = 11, base_height = 8.5,base_aspect_ratio = 1.3)
save_plot('classification2.png', p2, base_width = 11, base_height = 8.5,base_aspect_ratio = 1.3)
save_plot('lncRNA_classification2.tiff', p2, base_height = 8.5, base_width = 11, dpi = 300, compression = 'lzw')
# save_plot('lncRNA_classification2.pdf', p2, base_height = 8.5, base_width = 11, dpi = 300)