t <- "CorCor"
m <- "overlap_across_tissues"
lookup <- read.csv("trait_lookup.csv",header=T)
output <- matrix(,nrow(lookup),3)
rownames(output) <- lookup$Plot_Name
colnames(output) <- c("Conserved Exclusive","Overlap","Diverged Exclusive")
output2 <- output
rw.id <- c("Conserved.ExclusiveL2_0","OverlapL2_0","Diverged.ExclusiveL2_0")
for (i in 1:nrow(lookup)){
temp <- read.delim(paste0("results/",t,"_",m,
"/",lookup$LDSC_File_prefix[i],".results"))
output[i,] <- temp[match(rw.id,temp$Category),"Enrichment"]
output2[i,] <- temp[match(rw.id,temp$Category),"Enrichment_p"]
}
write.csv(output,paste0(t,"_",m,"_Enrichment.csv"),quote=F)
write.csv(output2,paste0(t,"_",m,"_Enrichment_p.csv"),quote=F)
