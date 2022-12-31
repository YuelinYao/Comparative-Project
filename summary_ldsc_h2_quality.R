t <- "CorCor"
m <- "toplast_10%"
library(stringr)
lookup <- read.csv("trait_lookup.csv",header=T)
output <- matrix(,nrow(lookup),11)
rownames(output) <- lookup$Plot_Name
colnames(output) <- c("No.SNPs.removed","chi2>","No.SNPs.remain","h2","h2.se",
"lambda_gc","mean_chi2","intercept","intercept.se","ratio","ratio.se")
for (i in 1:nrow(lookup)){
temp <- read.delim(paste0("results/",t,"_",m,
"/",lookup$LDSC_File_prefix[i],".log"))
a <- which(substr(temp[,1],1,7) == "Removed")
b <- which(substr(temp[,1],1,23) == "Total Observed scale h2")
c <- which(substr(temp[,1],1,10) == "Lambda GC:")
d <- which(substr(temp[,1],1,11) == "Mean Chi^2:")
e <- which(substr(temp[,1],1,10) == "Intercept:")
f <- which(substr(temp[,1],1,6) == "Ratio:")
g <- which(substr(temp[,1],1,9) == "Ratio < 0")
output[i,1:3] <- str_extract_all(temp[a,],"\\-*\\d+\\.*\\d*")[[1]][c(1,3,4)]
output[i,4:5] <-str_extract_all(temp[b,],"\\-*\\d+\\.*\\d*")[[1]][2:3]
output[i,6] <- str_extract_all(temp[c,],"\\-*\\d+\\.*\\d*")[[1]]
output[i,7] <- str_extract_all(temp[d,],"\\-*\\d+\\.*\\d*")[[1]][2]
output[i,8:9] <- str_extract_all(temp[e,],"\\-*\\d+\\.*\\d*")[[1]]
if (length(f) != 0){
output[i,10:11] <- str_extract_all(temp[f,],"\\-*\\d+\\.*\\d*")[[1]]
} else if (length(g) != 0){
output[i,10] <- "< 0"
output[i,11] <- NA
}
rm(a,b,c,d,e,f,g)
}
write.csv(output,paste0(t,"_",m,"_qualitycheck.csv"),quote=F)
