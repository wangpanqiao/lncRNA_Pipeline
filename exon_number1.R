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

lncRNA[, Type := toUpperFirstLetter(sub("_", " ", Type))]

# lncRNA.gtf <- unique(lncRNA[, Length := sum(Length), by = .(Gene, Type)], by = 'Gene')   ##base

# lncRNA.gtf <- unique(lncRNA[, exon_num := sum(exon_num), by = .(Gene, Type)], by = 'Gene')   ##base

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
lncRNA[exon_num > max.lncrna.len, exon_num := max.lncrna.len]
# p <- ggplot() + geom_density(data = lncRNA.gtf, aes(x = exon_num, colour = Type), size = 1.5) +
  # xlab('Transcript exon_num') + ylab('Density') +
  # scale_x_continuous(breaks = seq.int(2, max.lncrna.len, length.out = 20),
                     # labels = c(seq.int(2, max.lncrna.len, length.out = 19),
                                # paste0(max.lncrna.len, '+')), 
                     # expand = c(0.01, 0)) +
  # scale_y_continuous(expand = c(0.01, 0)) +
  # get(paste0('scale_color_',theme))()
  
p <- ggplot(data = lncRNA, aes(x = exon_num, group = Type)) + geom_density(adjust=2) + ###geom_density(aes(colour = Type ,fill=Type), alpha=0.4) +
  # geom_histogram(fill="Type", colour="grey60", size=.2,, alpha=0.4) +
  xlab('Transcript exon_num') + ylab('Density') +
  scale_x_continuous(breaks = seq.int(2, max.lncrna.len, length.out = 20),
                     labels = c(seq.int(2, max.lncrna.len, length.out = 19),
                                paste0(max.lncrna.len, '+')), 
                     expand = c(0.01, 0)) +
  scale_y_continuous(expand = c(0.01, 0)) +
  get(paste0('scale_color_',theme))()
  
save_plot('exon_num_distribution_with_type.png', p, base_width = 11, base_height = 8.5,base_aspect_ratio = 1.3)
save_plot('lncRNA_exon_num_distribution_with_type.tiff', p, base_height = 8.5, base_width = 11, dpi = 300, compression = 'lzw')
save_plot('lncRNA_exon_num_distribution_with_type.pdf', p, base_height = 8.5, base_width = 11, dpi = 300)