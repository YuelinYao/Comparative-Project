rm(list = ls())
####load package
library(ggplot2)
library(ggpubr)

#####load data
load("orthologous_cattle_human.Rdata")
load("human10830_cattle4866_17315gene.Rdata")
load("15696samples_17315genes.Rdata")


#################################################
tissue<-sort(unique(Meta_data_human$tissue_new))
Results <- array(NA,dim = c(length(tissue),1))
rownames(Results)<-tissue
colnames(Results)<-"expression_gene_number"
Results

#####caculate the gene number in each tissue for human
for( i in 1: length(tissue)){
  tissue_sample <-Meta_data_human[Meta_data_human$tissue_new==tissue[i],]

  human_expression_tissue<-human_expression[match(tissue_sample$SAMPID,rownames(human_expression)),]
  human_expression_tissue[1:10,1:10]
  expression_gene_mean<-apply(human_expression_tissue,2,median)
  Results[i,]<-length(expression_gene_mean[expression_gene_mean>0.1])
}


####caculate the gene number in each tissue for cattle
tissue_cattle<-sort(unique(Meta_data_cattle$tissue_new))
Results_cattle <- array(NA,dim = c(length(tissue),1))
rownames(Results_cattle)<-tissue_cattle
colnames(Results_cattle)<-"expression_gene_number"
Results_cattle

##### loop for each tissue
for( i in 1: length(tissue)){
  tissue_sample <-Meta_data_cattle[Meta_data_cattle$tissue_new==tissue_cattle[i],]
  cattle_expression_tissue<-Cattle_expression[match(tissue_sample$Sample,rownames(Cattle_expression)),]
  cattle_expression_tissue[1:10,1:10]
  expression_gene_mean<-apply(cattle_expression_tissue,2,median)
  Results_cattle[i,]<-length(expression_gene_mean[expression_gene_mean>0.1])
}


colors<-c("#00008B","#008B8B","#8B008B","#8B0000","#828282","#00B2EE","#0000FF","#FF0000","#FF00FF","#FF8C00","#FFFF00","#006400","#00FF00",
          "#CD950C","#FF6A6A","#CCCC33","pink","#FF1493","#8B7355","#00F5FF")

Results
Results<-as.data.frame(Results)
Results_cattle<-as.data.frame(Results_cattle)
table(rownames(Results_cattle)==rownames(Results))
Results<-Results/1000
Results_cattle<-Results_cattle/1000
combine<-data.frame(human=Results$expression_gene_number,cattle=Results_cattle$expression_gene_number,Tissue=rownames(Results_cattle))
combine

tiff("human-cattle-tissue_expressed_gene_number.tiff",res = 300, width = 650, height = 600,compression = "lzw")
par(mar=c(5.1,6,6,10))
DTI1 <-ggscatter(combine, x = "human", y = "cattle",label.y = 15,
                 
                 add = "reg.line", conf.int = TRUE, color =colors , size = 2, cor.coef.size = 4,
                 col = "deeppink", #label = rownames(comb),
                 
                 xlab = "No.expressed genes \n (thousands) in human", ylab = "No.expressed genes \n (thousands) in cattle",
                 
                 col.lab="red", cex.lab=0.8,cex.axis = 0.8)+stat_cor(method = "spearman", label.x = 12, label.y = 15)
DTI1+font("xlab", size =10 , color = "black")+
  font("ylab", size = 10, color = "black")+
  font("xy.text", size = 10, color = "black")
dev.off()



