rm(list = ls())


load("Cattle_Human_comparative.Rdata")

#BiocManager::install("MergeMaid")
library(MergeMaid)

Human_expression[1:10,1:10]
table(colnames(Human_expression)==Orthologous$`Gene stable ID`)
table(colnames(cattle_expression)==Orthologous$`cattle gene stable ID`)
#########Adipose
meta_data_human$tissue_new[meta_data_human$tissue_new=="Blood/Immune"]<-"Blood_Immune"
meta_data_cattle$tissue_new[meta_data_cattle$tissue_new=="Blood/Immune"]<-"Blood_Immune"
tissue_names<-sort(unique(meta_data_human$tissue_new))
tissue_names
for (i in 1:length(tissue_names)){
  

tissue_cattle<-meta_data_cattle[meta_data_cattle$tissue_new==tissue_names[i],]
tissue_exp_cattle<-cattle_expression[match(tissue_cattle$BioSample,rownames(cattle_expression)),]
tissue_exp_cattle<-t(tissue_exp_cattle)


tissue_human<-meta_data_human[meta_data_human$tissue_new==tissue_names[i],]
tissue_exp_human<-Human_expression[match(tissue_human$SAMPID,rownames(Human_expression)),]
tissue_exp_human<-t(tissue_exp_human)
tissue_exp_human[1:10,1:10]

rownames(tissue_exp_cattle)<-rownames(tissue_exp_human)
merged <-mergeExprs(tissue_exp_human,tissue_exp_cattle)
tissue_corcor <-intCor(merged,method="pearson",exact=F)

save(tissue_corcor,file=paste0("corcor_",tissue_names[i],".Rdata"))
}


