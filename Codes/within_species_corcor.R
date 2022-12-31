rm(list = ls())
#setwd("/Users/yaoyuelin/Downloads/cattle_GTEX")
load("cattle_Human_comparative.Rdata")

#BiocManager::install("MergeMaid")
library(MergeMaid)

#BiocManager::install("MergeMaid")


Human_expression[1:10,1:10]
table(colnames(Human_expression)==Orthologous$`Gene stable ID`)
table(colnames(cattle_expression)==Orthologous$`cattle gene stable ID`)
meta_data_human$tissue_new[meta_data_human$tissue_new=="Blood/Immune"]<-"Blood_Immune"
meta_data_cattle$tissue_new[meta_data_cattle$tissue_new=="Blood/Immune"]<-"Blood_Immune"
tissue_names<-sort(unique(meta_data_human$tissue_new))
tissue_names
#########cattle

tissue_l<-c("Adipose","Blood_Immune","Brain","Colon","Heart","Kidney")

for (i in 7:16){

  
  tissue_cattle<-meta_data_cattle[meta_data_cattle$tissue_new==tissue_names[i],]
  tissue_exp_cattle<-cattle_expression[match(tissue_cattle$BioSample,rownames(cattle_expression)),]
  tissue_exp_cattle<-t(tissue_exp_cattle)
  tissue_l<-c(tissue_l,tissue_names[i])
  tissue_l<-tissue_l[!is.na( tissue_l)]
  tissue_l
  tissue_names_r<-tissue_names[!tissue_names%in%tissue_l]
  length(tissue_names_r)
  for (j in 1:length(tissue_names_r)){
    
    tissue_names_r
    tissue_cattle_select<-meta_data_cattle[meta_data_cattle$tissue_new==tissue_names_r[j],]
    tissue_exp_cattle_select<-cattle_expression[match(tissue_cattle_select$BioSample,rownames(cattle_expression)),]
    tissue_exp_cattle_select<-t(tissue_exp_cattle_select)
    
    rownames(tissue_exp_cattle)<-rownames(tissue_exp_cattle_select)
    merged <-mergeExprs(tissue_exp_cattle,tissue_exp_cattle_select)
    tissue_corcor <-intCor(merged,method="pearson",exact=F)
    
    
    save(tissue_corcor,file=paste0("Within_cattle_",tissue_names[i],"_",tissue_names_r[j],".Rdata"))
  }
}

####human

load("cattle_Human_comparative.Rdata")

#BiocManager::install("MergeMaid")
#BiocManager::install("MergeMaid")
library(MergeMaid)
Human_expression[1:10,1:10]
table(colnames(Human_expression)==Orthologous$`Gene stable ID`)
table(colnames(cattle_expression)==Orthologous$`cattle gene stable ID`)
meta_data_human$tissue_new[meta_data_human$tissue_new=="Blood/Immune"]<-"Blood_Immune"
meta_data_cattle$tissue_new[meta_data_cattle$tissue_new=="Blood/Immune"]<-"Blood_Immune"
tissue_names<-sort(unique(meta_data_human$tissue_new))
tissue_names

tissue_l<-c("Adipose","Blood_Immune","Brain","Colon","Heart","Kidney")
tissue_l
for (k in 7:16){
k=7  
  
  tissue_human<-meta_data_human[meta_data_human$tissue_new==tissue_names[k],]
  tissue_exp_human<-Human_expression[match(tissue_human$SAMPID,rownames(Human_expression)),]
  tissue_exp_human<-t(tissue_exp_human)
  tissue_exp_human[1:10,1:10]
  tissue_l<-c(tissue_l,tissue_names[k])
  tissue_l<-tissue_l[!is.na( tissue_l)]
  tissue_l
  tissue_names_r<-tissue_names[!tissue_names%in%tissue_l]

  for (z in 1:length(tissue_names_r)){
    
    tissue_names_r
    tissue_human_select<-meta_data_human[meta_data_human$tissue_new==tissue_names_r[z],]
    tissue_exp_human_select<-Human_expression[match(tissue_human_select$SAMPID,rownames(Human_expression)),]
    tissue_exp_human_select<-t(tissue_exp_human_select)
    
    rownames(tissue_exp_human_select)<-rownames(tissue_exp_human)
    merged <-mergeExprs(tissue_exp_human,tissue_exp_human_select)
    tissue_corcor <-intCor(merged,method="pearson",exact=F)
    
    
    save(tissue_corcor,file=paste0("Within_human_",tissue_names[k],"_",tissue_names_r[z],".Rdata"))
  }
}
