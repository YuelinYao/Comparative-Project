type <- "CorCor"
model <- "overlap_across_tissues"
load("1000G.eur.bim.RData")
for (chr in 1:22){
output <- cbind(chr.bim[[chr]],1)
colnames(output)[5] <- "All.SNPs"
temp.c <- read.table(paste0("annotation/",type,"_",model,"/conserved.",chr)
,header=T)
colnames(temp.c) <- "Conserved.Exclusive"
temp.o <- read.table(paste0("annotation/",type,"_",model,"/overlap.",chr)
,header=T)
colnames(temp.o) <- "Overlap"
temp.d <- read.table(paste0("annotation/",type,"_",model,"/diverged.",chr)
,header=T)
colnames(temp.d) <- "Diverged.Exclusive"
output <- cbind(output,temp.c,temp.o,temp.d)
write.table(output,paste0("annotation/",type,"_",model,"/",
type,"_",model,".",chr,".annot"),quote=F,row.names=F,col.names=T)
rm(output,temp.c,temp.o,temp.d)
}

