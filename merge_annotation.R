type <- "CorCor"
model <- "toplast_10%"
tissue <- read.csv("tissue_lookup.csv",header=T)
load("1000G.eur.bim.RData")
for (chr in 1:22){
output <- cbind(chr.bim[[chr]],1)
colnames(output)[5] <- "All.SNPs"
for (i in 1:nrow(tissue)){
temp.last <- read.table(paste0("annotation/",type,"_",model,"/",
tissue$File_Prefix[i],"-last10%.",chr),header=T)
colnames(temp.last) <- paste0(tissue$File_Prefix[i],".last10%")
temp.top <- read.table(paste0("annotation/",type,"_",model,"/",
tissue$File_Prefix[i],"-top10%.",chr),header=T)
colnames(temp.top) <- paste0(tissue$File_Prefix[i],".top10%")
output <- cbind(output,temp.last,temp.top)
rm(temp.last,temp.top)
}
write.table(output,paste0("annotation/",type,"_",model,"/",
type,"_",model,".",chr,".annot"),quote=F,row.names=F,col.names=T)
rm(output)
}
