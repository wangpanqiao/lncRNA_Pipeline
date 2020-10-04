rm(list=ls())
library(data.table)
kallisto.files <- list.files("./",pattern="*.tsv")
map.file <- read.table("./map.file",header=F,sep="\t")
names(map.file) <- c("gene","ID","Type")

trans.kallisto.count <- c()
trans.kallisto.TPM <- c()
pc.names <- c()
for(pc in kallisto.files){
  pc.exp <- read.table(pc,header=T,sep="\t",row.names=1)
  trans.kallisto.count <- cbind(trans.kallisto.count,pc.exp[,3])
  trans.kallisto.TPM <- cbind(trans.kallisto.TPM,pc.exp[,4])
  pc.names <- rownames(pc.exp)
}
#assign gene names
rownames(trans.kallisto.count) <- pc.names
rownames(trans.kallisto.TPM) <- pc.names

# paste(unlist(strsplit(kallisto.files[1],'_'))[1],unlist(strsplit(kallisto.files[1],'_'))[2],sep="_")
pc.samples <- unlist(lapply(kallisto.files,function(x){paste(unlist(strsplit(x,"_"))[1],unlist(strsplit(x,"_"))[2],sep="_")}))
pc.samples
colnames(trans.kallisto.count) <- pc.samples
colnames(trans.kallisto.TPM) <- pc.samples

# write.table(as.data.frame(trans.kallisto.count) ,file="kallisto.count.txt",sep="\t",quote=F,row.names=T,col.names=T)
# write.table(as.data.frame(trans.kallisto.TPM) ,file="kallisto.tpm.txt",sep="\t",quote=F,row.names=T,col.names=T)


row.names(map.file)<-map.file[,2]
map.file=map.file[pc.names,]

trans.kallisto.count.merge <- cbind(map.file[,c(2,3)],trans.kallisto.count)
head(trans.kallisto.count.merge)

fwrite(trans.kallisto.count.merge ,file="kallisto.count.txt",sep="\t",quote=F,row.names=F)


trans.kallisto.TPM.merge <- cbind(map.file[,c(2,3)],trans.kallisto.TPM)

fwrite(trans.kallisto.TPM.merge ,file="kallisto.tpm.txt",sep="\t",quote=F,row.names=F)

