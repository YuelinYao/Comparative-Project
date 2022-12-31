rm(list = ls()) 
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(readxl)
load("overlap.Rdata")
#########
overlap
phyper(overlap$Overlap-1 , overlap$Human, 17315-overlap$Human , overlap$Cattle,lower.tail= FALSE)

overlap$P_value<-phyper(overlap$Overlap-1 , overlap$Human, 17315-overlap$Human , overlap$Cattle,lower.tail= FALSE)

Mydata_raw_FDR <- p.adjust(overlap$P_value,method = "BH")
Mydata_raw_FDR[Mydata_raw_FDR<=0.05] <- "*"
Mydata_raw_FDR[Mydata_raw_FDR<=0.01] <- "**"
Mydata_raw_FDR[Mydata_raw_FDR<=0.001] <- "***"
names(Mydata_raw_FDR)<-rownames(overlap)

Mydata_raw_FDR


############

overlap_sort<-overlap[order(overlap$Cattle),]
overlap_sort
overlap_sort_t<-t(overlap_sort)
overlap_sort_t




tiff(file = "tissue_specific_gene_number.tiff",##reqiured to change
     res = 300, width = 5500, height = 2200,compression = "lzw")
par(mar=c(13,10,7,0.3))
my<-barplot(overlap_sort_t[1:3,],beside=T,ylim = c(0,5000),col =c("blue","red","orange"),las="2",axisnames=T,cex.axis = 2,cex.lab=2,
cex.names =2 )
legend("topleft", legend = c("Cattle","Human","Overlap") , 
       col = c("blue","red","orange") , 
       bty = "n", pch=15 , pt.cex = 2, cex = 2, horiz = FALSE, inset = c(0.05, 0.05))
title(ylab="Number of tissue-specific genes", line=6, cex.lab=2.5)

text(my[2,]+0.8,overlap_sort_t[2,]+200,Mydata_raw_FDR,cex=1.8,font = 2,col = "brown")
dev.off()



