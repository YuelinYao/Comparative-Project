t <- "CorCor"
m <- "toplast_10%"
library(stringr)
lookup <- read.csv("trait_lookup.csv",header=T)
tissue <- read.csv("tissue_lookup.csv",header=T)
output <- matrix(,nrow(lookup),40)
rownames(output) <- lookup$Plot_Name
colnames(output) <- paste0(rep(tissue$Tissue_Name,2),
rep(c(" Top10%"," Bot10%"),each=20))
output2 <- output
rw.id <- paste0(rep(tissue$File_Prefix,2),
rep(c(".top10%L2_0",".last10%L2_0"),each=20))
for (i in 1:nrow(lookup)){
temp <- read.delim(paste0("results/",t,"_",m,
"/",lookup$LDSC_File_prefix[i],".results"))
output[i,] <- temp[match(rw.id,temp$Category),"Enrichment"]
output2[i,] <- temp[match(rw.id,temp$Category),"Enrichment_p"]
}
write.csv(output,paste0(t,"_",m,"_Enrichment.csv"),quote=F)
write.csv(output2,paste0(t,"_",m,"_Enrichment_p.csv"),quote=F)
