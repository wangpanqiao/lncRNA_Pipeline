library(ggplot2)
library(cowplot)
library(data.table)
# library(ggpubr)
library(ggrepel)

toUpperFirstLetter <- function(x) {
  paste0(toupper(substr(x, 1, 1)), substring(x, 2))
}
lncRNA <- fread("basic_charac.txt", sep="\t")
setnames(lncRNA, c('Gene', 'Transcript', 'Type', 'Potential', 'Length', 'exon_num', 'Class'))

lncRNA <- lncRNA[lncRNA$Type!="protein_coding",]

# lncRNA[, Type := toUpperFirstLetter(sub("_", " ", Type))]

lncRNA[, Class := sub("Exonic_Sense,Intronic_Sense", "Sense", Class)]

lncRNA[, Class := sub("Intronic_Sense,Exonic_Sense", "Sense", Class)]

lncRNA[, Class := sub("Exonic_Sense", "Sense", Class)]

lncRNA[, Class := sub("Intronic_Sense", "Intronic", Class)]

# lncRNA[, Class := sub("Bidirectional", "Intergenic", Class)]

cls <- lncRNA[, .(count = .N), by = Class]
colors <- c('rgb(211,94,96)', 'rgb(128,133,133)', 'rgb(144,103,167)', 'rgb(171,104,87)', 'rgb(114,147,203)')

cls1 = cls[order(cls$count, decreasing = TRUE),]
cls1$Class=factor(cls1$Class,ordered = TRUE, levels = c("Intergenic", "Sense", "Antisense", "Bidirectional","Intronic"))

write.table(as.data.frame(cls1),'classification.tsv',quote=F,sep="\t",row.names=F,col.names=T)

print(cls1)
# myLabel = as.vector(cls1$Class)   
# myLabel = paste(myLabel, "(", round(cls1$count / sum(cls1$count) * 100, 2), "%)", sep = "")  

label_value <- paste('(', cls1$count, ', ', round(cls1$count/sum(cls1$count) * 100, 1), '%)', sep = '')
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


# p1=ggplot(data=cls1, mapping=aes(x="", y=count, fill=Class)) +
  # geom_bar(stat="identity",width=0.7) + coord_polar(theta = 'y', start = 0)  + 
  # blank_theme + 
  # labs(x = "", y = "", title = "Percentage distribution on lncRNA classification (%)")  +
  # geom_text(aes(y = cls1$count/2 + c(0, cumsum(cls1$count)[-length(cls1$count)]),x=sum(cls1$count)/2400, label = label_value),size=5) ## +
  #theme(legend.position = "none")  + 
  #geom_text(aes(y = cls1$count/2 + c(0, cumsum(cls1$count)[-length(cls1$count)]) , label = label),size=5)
  #geom_text(aes(y = cls1$count/2 + c(0, cumsum(cls1$count)[-length(cls1$count)]) , x = sum(cls1$count)/3000, label = label),size=5) 

p = ggplot(cls1, aes(x = Class, y = count, fill = Class)) + geom_bar(stat = "identity", alpha = 0.7) + coord_polar() +
    blank_theme + 
    labs(x = "", y = "")  + ##, title = "Percentage distribution on lncRNA classification (%)"
    geom_text_repel(aes(label = label_value, y = count + 50,x=Class), face = 'bold',size = 6) +
    theme(legend.title = element_text(colour="black", size=18, face="bold")) +
    theme(legend.text = element_text(colour="black", size = 15, face = "bold"))

# png('classification3.png')
# p
# dev.off()


save_plot('classification3.png', p, base_width = 11, base_height = 8.5,base_aspect_ratio = 1.3)
# save_plot('classification3.png', p, base_width = 11, base_height = 8.5,base_aspect_ratio = 1.3)
save_plot('lncRNA_classification3.tiff', p, base_height = 8.5, base_width = 11, dpi = 300, compression = 'lzw')
# save_plot('lncRNA_classification2.pdf', p2, base_height = 8.5, base_width = 11, dpi = 300)