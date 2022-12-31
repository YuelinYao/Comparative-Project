rm(list = ls()) 

#############
load("integratedallgenes.Rdata")
load("human10830_cattle4866_17315gene.Rdata")
load("15696samples_17315genes.Rdata")

experiment$annotation<-paste(experiment$tissue_new,experiment$species,sep=" - ")
head(experiment)
length(table(experiment$annotation))
median_tissue_tpm<-aggregate(expression_inter,list(as.factor(experiment$annotation)),median)
rownames(median_tissue_tpm)<-median_tissue_tpm[,1]
median_tissue_tpm[1:10,1:10]
corr_mat=cor(t(median_tissue_tpm[,-1]))
dim(corr_mat)
corr_mat[1:10,1:10]


median<-median_tissue_tpm[,-1]
median[1:10,1:10]

P_value<-array(NA,dim = c(40,40))
colnames(P_value)<-rownames(median)
rownames(P_value)<-rownames(median)

Results<-P_value
i=1
j=1
a$p.value
a$estimate
for (i in 1:40){
  for (j in 1:40){
    a<-cor.test(as.numeric(median[i,]),as.numeric(median[j,]))
    P_value[i,j]<-a$p.value
    Results[i,j]<-a$estimate}
}


Adipose<-corr_mat['Adipose - Human','Adipose - Cattle']
Adrenal<-corr_mat['Adrenal - Human','Adrenal - Cattle']
Blood_Immune<-corr_mat['Blood/Immune - Human','Blood/Immune - Cattle']
Brain<-corr_mat['Brain - Human','Brain - Cattle']
Heart<-corr_mat['Heart - Human','Heart - Cattle']
Kidney<-corr_mat['Kidney - Human','Kidney - Cattle']
Large_intestine<-corr_mat['Large intestine - Human','Large intestine - Cattle']
Liver<-corr_mat['Liver - Human','Liver - Cattle']
Lung<-corr_mat['Lung - Human','Lung - Cattle']
Mammary<-corr_mat['Mammary - Human','Mammary - Cattle']
Muscle<-corr_mat['Muscle - Human','Muscle - Cattle']
Ovary<-corr_mat['Ovary - Human','Ovary - Cattle']
Pituitary<-corr_mat['Pituitary - Human','Pituitary - Cattle']
Salivary_Gland<-corr_mat['Salivary Gland - Human','Salivary Gland - Cattle']
Skin<-corr_mat['Skin - Human','Skin - Cattle']
Small_intestine<-corr_mat['Small intestine - Human','Small intestine - Cattle']
Spleen<-corr_mat['Spleen - Human','Spleen - Cattle']
Stomach<-corr_mat['Stomach - Human','Stomach - Cattle']
Testis<-corr_mat['Testis - Human','Testis - Cattle']
Uterus<-corr_mat['Uterus - Human','Uterus - Cattle']

tissue_cor<-cbind(Adipose,
                  Adrenal,
                  Blood_Immune,
                  Brain,
                  Heart,
                  Kidney,
                  Large_intestine,
                  Liver,
                  Lung,
                  Mammary,
                  Muscle,
                  Ovary,
                  Pituitary,
                  Salivary_Gland,
                  Skin,
                  Small_intestine,
                  Spleen,
                  Stomach,
                  Testis,
                  Uterus)
dim(tissue_cor)
tissue_cor
colors<-c("#00008B","#008B8B","#8B008B","#8B0000","#828282","#00B2EE","#0000FF","#FF0000","#FF00FF","#FF8C00","#FFFF00","#006400","#00FF00",
          "#CD950C","#FF6A6A","#CCCC33","pink","#FF1493","#8B7355","#00F5FF")
names(colors)=colnames(tissue_cor)

tiff("correlation_all_genes_barplot.tiff",res = 300, width = 2200, height = 1500,compression = "lzw")
par(mfrow=c(1,1),mar=c(10,6,4,4),oma = c(0,0,0,0) + 0.1)
order1 <- order(tissue_cor)
order1
tissue_cor<-tissue_cor[,order(tissue_cor)]
barplot(tissue_cor,ylim = c(0,1),legend=F,las="2",cex.axis = 1.3,cex.lab=1.6,cex.names=1.5,
        axisnames=T,col = colors[order1],ylab="Correlation of gene expression")
dev.off()



#########3

overlap_t
over<-overlap_t["Overlap",]
tissue_cor_list<-as.data.frame(tissue_cor)
tissue_cor_list
tissue_cor_list<-t(tissue_cor_list)
tissue_cor_list
over_list<-as.data.frame(over)
over_list
tissue_cor_list
tissue_cor_list<-t(tissue_cor_list)
table(rownames(over_list)==rownames(tissue_cor_list))
colnames(over_list)<-"Overlap"
colnames(tissue_cor_list)<-"Tissue_correlation"
over_list
table(rownames(tissue_cor_list)==rownames(over_list))


comb<-cbind(tissue_cor_list,over_list)
rownames(comb)<-rownames(over_list)
comb
library("ggpubr")
library(magrittr)

tiff("correlation-tissue_gene_number.tiff",res = 300, width = 2000, height = 2000,compression = "lzw")

DTI1 <-ggscatter(comb, x = "tissue_cor_list", y = "over",
                 
                 add = "reg.line", conf.int = TRUE, color =colors , size = 5, 
                 
                 cor.coef = TRUE,cor.method = "spearman", col = "deeppink", #label = rownames(comb),
                 
                 xlab = "Correlation of gene expression", ylab = "Number of tissue-specific genes",
                 
                 col.lab="red", cex.lab=3,cex.axis = 1.5)
DTI1+font("xlab", size = 20, color = "black")+
  font("ylab", size = 20, color = "black")+
  font("caption", size = 20, color = "orange")+
  font("xy.text", size = 20, color = "black", face = "bold")

dev.off()

########3
overlap_t
over<-overlap_t["Percentage",]
tissue_cor_list<-as.data.frame(tissue_cor)
tissue_cor_list
tissue_cor_list<-t(tissue_cor_list)
tissue_cor_list
over_list<-as.data.frame(over)
over_list
tissue_cor_list
tissue_cor_list<-t(tissue_cor_list)
table(rownames(over_list)==rownames(tissue_cor_list))
colnames(over_list)<-"Overlap_Percentage"
colnames(tissue_cor_list)<-"Tissue_correlation"
over_list
table(rownames(tissue_cor_list)==rownames(over_list))


comb<-cbind(tissue_cor_list,over_list)
rownames(comb)<-rownames(over_list)

library("ggpubr")
library(magrittr)

tiff("correlation-Overlap_Percentage.tiff",res = 300, width = 800, height = 800,compression = "lzw")
DTI1 <-ggscatter(comb, x = "Tissue_correlation", y = "Overlap_Percentage",
                 
                 add = "reg.line", conf.int = TRUE, color =colors , size = 3, cor.coef.size = 4,
                 
                 cor.coef = TRUE,cor.method = "spearman", col = "deeppink", #label = rownames(comb),
                 
                 xlab = "Correlation of gene expression", ylab = "Percentage of tissue-specific genes",
                 
                 col.lab="red",cex.lab=2,cex.axes=2)
DTI1+font("xlab", size = 10, color = "black")+
  font("ylab", size = 10, color = "black")+
  font("caption", size = 10, color = "orange")+
  font("xy.text", size = 10, color = "black", face = "bold")

dev.off()

#########3
overlap
colors<-c("#00008B","#008B8B","#8B008B","#8B0000","#828282","#00B2EE","#0000FF","#FF0000","#FF00FF","#FF8C00","#FFFF00","#006400","#00FF00",
          "#CD950C","#FF6A6A","#CCCC33","pink","#FF1493","#8B7355","#00F5FF")
names(colors)=rownames(overlap)

cor.test(overlap$Cattle,overlap$Human,method = "spearman")
library(ggpubr)
tiff("tissue_specific_correlation_spearman.tiff",res = 300, width = 700, height = 700,compression = "lzw")
DTI1 <-ggscatter(overlap, x = "Cattle", y = "Human",
                 
                 add = "reg.line", conf.int = TRUE, color =colors , size = 3, cor.coef.size = 4,
                 
                 cor.coef = TRUE,cor.method = "spearman", col = "deeppink", #label = rownames(comb),
                 
                 xlab = " ", ylab = " ",
                 
                 col.lab="red",cex.lab=1.5,cex.axes=1.5)
DTI1+font("xlab", size = 6, color = "black")+
  font("ylab", size = 6, color = "black")+
  font("caption", size = 6, color = "orange")+
  font("xy.text", size = 6, color = "black", face = "bold")

dev.off()



