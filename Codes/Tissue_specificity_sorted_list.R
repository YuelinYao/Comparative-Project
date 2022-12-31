rm(list = ls()) 
load("human10830_cattle4866_17315gene.Rdata")
load("15696samples_17315genes.Rdata")
load("orthologous_cattle_human.Rdata")

human_file<-dir("~/Comparative/tissue_specific_limma/human",full.names = T)
human_file
humantissue<-substr(human_file,nchar("~/Comparative/tissue_specific_limma/human/myresults_")+1,nchar(human_file)-nchar("_human.txt"))
humantissue

cattle_file<-dir("~/Comparative/tissue_specific_limma/cattle",full.names = T)
cattle_file
cattletissue<-substr(cattle_file,nchar("~/Comparative/tissue_specific_limma/cattle/myresults_")+1,nchar(cattle_file)-nchar(".txt"))
cattletissue

tissue<-unique(sort(Meta_data_cattle$tissue_new))
Results1 <- array(NA,dim = c(length(tissue),7))

i=1
rownames(Results1)<-tissue
colnames(Results1)<-c("top10%","10%-20%","20%-30%","40%-50%","70%-80%","80%-90%","bottom10%")
Results1
for (i in 1:length(humantissue)){
  cattle1<- read.table(cattle_file[i],header = T,stringsAsFactors = F)
  human1<- read.table(human_file[i],header = T,stringsAsFactors = F)

  rownames(human1)<-Orthologous_human_cattle$`Cow gene stable ID`[match(rownames(human1),Orthologous_human_cattle$`Gene stable ID`)]
  
  
  top10_human_tissue<-rownames(human1)[1:1731]
  top10_cattle_tissue<-rownames(cattle1)[1:1731]
  
  
  tissue_inter<-intersect(top10_cattle_tissue,top10_human_tissue)
  Results1[i,1]<-length(tissue_inter)
  
  top10_20_human_tissue<-rownames(human1)[1732:3462]
  
  top10_20_cattle_tissue<-rownames(cattle1)[1732:3462]
  
  tissue_inter_10_20<-intersect(top10_20_cattle_tissue,top10_20_human_tissue)
  length(tissue_inter_10_20)
  Results1[i,2]<-length(tissue_inter_10_20)
  
  top20_30_human_tissue<-rownames(human1)[3463:5193]
  top20_30_cattle_tissue<-rownames(cattle1)[3463:5193]
  
  tissue_inter_20_30<-intersect(top20_30_cattle_tissue,top20_30_human_tissue)
  Results1[i,3]<-length(tissue_inter_20_30)
  
  
  
  top40_50_human_tissue<-rownames(human1)[6926:8656]
  
  top40_50_cattle_tissue<-rownames(cattle1)[6926:8656]
 
  tissue_inter_40_50<-intersect(top40_50_cattle_tissue,top40_50_human_tissue)
  length(tissue_inter_40_50)
  Results1[i,4]<-length(tissue_inter_40_50)
  
  top70_80_human_tissue<-rownames(human1)[12121:13851]
  top70_80_cattle_tissue<-rownames(cattle1)[12121:13851]
  
  tissue_inter_70_80<-intersect(top70_80_cattle_tissue,top70_80_human_tissue)
  length(tissue_inter_70_80)
  Results1[i,5]<-length(tissue_inter_70_80)
 
  
  top80_90_human_tissue<-rownames(human1)[13853:15583]
  top80_90_cattle_tissue<-rownames(cattle1)[13853:15583]
  
  
  tissue_inter_80_90<-intersect(top80_90_cattle_tissue,top80_90_human_tissue)
  length(tissue_inter_80_90)
  
  Results1[i,6]<-length(tissue_inter_80_90)
  
  
  last10_human_tissue<-rownames(human1)[15585:17315]
  last10_cattle_tissue<-rownames(cattle1)[15585:17315]
  
  
  
  
  tissue_inter_last10<-intersect(last10_cattle_tissue,last10_human_tissue)
  length(tissue_inter_last10)
  Results1[i,7]<-length(tissue_inter_last10)
  
}
Results1
Results1_percent<-Results1/1731
Results1_percent
col<-c("#00008B","#008B8B","#8B008B","#8B0000","#828282","#00B2EE","#0000FF","#FF0000","#FF00FF","#FF8C00","#FFFF00","#006400","#00FF00",
       "#CD950C","#FF6A6A","#CCCC33","pink","#FF1493","#8B7355","#00F5FF")
tiff(file = "t-value_new.tiff",##reqiured to change
     res = 300, width = 5000, height = 2000,compression = "lzw")
par(mar=c(8,8,6,1))
my<-barplot(Results1_percent,beside = T,col=col,cex.axis = 2,cex.names = 2,ylim = c(0,0.7))
title(ylab="Percentage of shared orthologs", xlab = "Windows of genes sorted by tissue specificity",line=4, cex.lab=2.5)
legend("top", legend = tissue , 
       col = col, ncol = 4,xpd=T,
       bty = "n", pch=15 , pt.cex = 1.5, cex = 1.5, horiz = FALSE, inset = c(0.05, 0.05))
dev.off()





######################3
tissue<-unique(sort(experiment$tissue_new))
Results1 <- array(NA,dim = c(length(tissue),7))
experiment$annotation<-paste(experiment$tissue_new,experiment$species,sep="-")
variance_tissue_tpm<-aggregate(expression,list(as.factor(experiment$annotation)),mad)
variance_tissue_tpm
rownames(variance_tissue_tpm)<-variance_tissue_tpm[,1]
variance_tissue_tpm<-variance_tissue_tpm[,-1]
variance_tissue_tpm[,1:10]

i=1
rownames(Results1)<-tissue
colnames(Results1)<-c("top10%","10%-20%","20%-30%","40%-50%","70%-80%","80%-90%","bottom10%")
for (i in 1:length(tissue)){
  human_name<-paste(tissue[i],"Human",sep="-")
  cattle_name<-paste(tissue[i],"Cattle",sep="-")
  human_tissue<-sort(variance_tissue_tpm[human_name,],decreasing = T)
  cattle_tissue<-sort(variance_tissue_tpm[cattle_name,],decreasing = T)
  top10_human_tissue<-colnames(human_tissue)[1:1731]
  top10_cattle_tissue<-colnames(cattle_tissue)[1:1731]
  
  
  tissue_inter<-intersect(top10_cattle_tissue,top10_human_tissue)
  Results1[i,1]<-length(tissue_inter)
  
  top10_20_human_tissue<-colnames(human_tissue)[1732:3463]
  top10_20_cattle_tissue<-colnames(cattle_tissue)[1732:3463]
  
  tissue_inter_10_20<-intersect(top10_20_cattle_tissue,top10_20_human_tissue)
  length(tissue_inter_10_20)
  Results1[i,2]<-length(tissue_inter_10_20)
  
  top20_30_human_tissue<-colnames(human_tissue)[3463:5194]
  top20_30_cattle_tissue<-colnames(cattle_tissue)[3463:5194]
  
  tissue_inter_20_30<-intersect(top20_30_cattle_tissue,top20_30_human_tissue)
  Results1[i,3]<-length(tissue_inter_20_30)
  
  
  
  top40_50_human_tissue<-colnames(human_tissue)[6926:8657]
  top40_50_cattle_tissue<-colnames(cattle_tissue)[6926:8657]
  
  tissue_inter_40_50<-intersect(top40_50_cattle_tissue,top40_50_human_tissue)
  length(tissue_inter_40_50)
  Results1[i,4]<-length(tissue_inter_40_50)
  
  top70_80_human_tissue<-colnames(human_tissue)[12121:13852]
  top70_80_cattle_tissue<-colnames(cattle_tissue)[12121:13852]
  
  tissue_inter_70_80<-intersect(top70_80_cattle_tissue,top70_80_human_tissue)
  length(tissue_inter_70_80)
  Results1[i,5]<-length(tissue_inter_70_80)
  
  
  top80_90_human_tissue<-colnames(human_tissue)[13853:15584]
  top80_90_cattle_tissue<-colnames(cattle_tissue)[13853:15584]
  
  tissue_inter_80_90<-intersect(top80_90_cattle_tissue,top80_90_human_tissue)
  length(tissue_inter_80_90)
  
  Results1[i,6]<-length(tissue_inter_80_90)
  
  
  last10_human_tissue<-colnames(human_tissue)[15584:17315]
  last10_cattle_tissue<-colnames(cattle_tissue)[15584:17315]
  
  
  
  
  tissue_inter_last10<-intersect(last10_cattle_tissue,last10_human_tissue)
  length(tissue_inter_last10)
  Results1[i,7]<-length(tissue_inter_last10)
  
}
Results1
Results1_percent<-Results1/1731
Results1_percent
col<-c("#00008B","#008B8B","#8B008B","#8B0000","#828282","#00B2EE","#0000FF","#FF0000","#FF00FF","#FF8C00","#FFFF00","#006400","#00FF00",
       "#CD950C","#FF6A6A","#CCCC33","pink","#FF1493","#8B7355","#00F5FF")
tiff(file = "mad.tiff",##reqiured to change
     res = 300, width = 5200, height = 2000,compression = "lzw")
par(mar=c(8,8,6,1))
my<-barplot(Results1_percent,beside = T,col=col,cex.axis = 1.8,cex.names = 1.8,ylim = c(0,0.7))
title(ylab="Percentage of shared orthologs", xlab = "Windows of genes sorted by expression MAD",line=4, cex.lab=2.5)
legend("top", legend = tissue , 
       col = col, ncol = 4,xpd=T,
       bty = "n", pch=15 , pt.cex = 1.5, cex = 1.5, horiz = FALSE, inset = c(0.05, 0.05))
dev.off()
