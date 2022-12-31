rm(list = ls()) 

library(limma)

load("human10830_cattle4866_17315gene.Rdata")
load("15696samples_17315genes.Rdata")
load("orthologous_cattle_human.Rdata")

####cattle

Meta_data_cattle$tissue_new[Meta_data_cattle$tissue_new=="Blood/Immune"]<-"Blood_Immune"
tissue_name<-unique(sort(Meta_data_cattle$tissue_new))
Meta_data_cattle_sort<-Meta_data_cattle[order(Meta_data_cattle$tissue_new),]
cattle_expression<-t(Cattle_expression)
cattle_expression_sort<-cattle_expression[,match(Meta_data_cattle_sort$Sample,colnames(cattle_expression))]
cattle_expression_sort<-log2(cattle_expression_sort+0.25)

i=1

for (i in 1:20){
  tissue_index<-which(Meta_data_cattle_sort$tissue_new==tissue_name[i])
  df<-c(rep("others",tissue_index[1]-1),rep("tissue",length(tissue_index)),rep("others",dim(Meta_data_cattle_sort)[1]-length(tissue_index)-(tissue_index[1]-1)))
  df
  design <- model.matrix(~-1+factor(df))
  colnames(design) <- c("Others","Tissue")
  contrastmatrix <- makeContrasts(Tissue-Others, levels=design)
  fit <- lmFit(cattle_expression_sort, design) ##issue these commands to fit the model and make the contrasts
  fit2 <- contrasts.fit(fit, contrastmatrix)
  fit2 <- eBayes(fit2)
  fit2
  myresults <-topTable(fit2,coef=1, adjust="fdr", number=nrow(cattle_expression_sort))
  write.table(myresults,file=paste0("/tissue_specific_limma/cattle/myresults_",tissue_name[i],".txt"),quote = F)
  }




########human
Meta_data_human$tissue_new[Meta_data_cattle$tissue_new=="Blood/Immune"]<-"Blood_Immune"
tissue_name<-unique(sort(Meta_data_human$tissue_new))
Meta_data_human_sort<-Meta_data_human[order(Meta_data_human$tissue_new),]
Human_expression<-t(human_expression)
human_expression_sort<-Human_expression[,match(Meta_data_human_sort$SAMPID,colnames(Human_expression))]
human_expression_sort<-log2(human_expression_sort+0.25)
human_expression_sort[1:10,1:10]

for (i in 1:20){
  tissue_index<-which(Meta_data_human_sort$tissue_new==tissue_name[i])
  df<-c(rep("others",tissue_index[1]-1),rep("tissue",length(tissue_index)),rep("others",dim(Meta_data_cattle_sort)[1]-length(tissue_index)-(tissue_index[1]-1)))
  df
  design <- model.matrix(~-1+factor(df))
  colnames(design) <- c("Others","Tissue")
  contrastmatrix <- makeContrasts(Tissue-Others, levels=design)
  fit <- lmFit(human_expression_sort, design) ##issue these commands to fit the model and make the contrasts
  fit2 <- contrasts.fit(fit, contrastmatrix)
  fit2 <- eBayes(fit2)
  fit2
  myresults <-topTable(fit2,coef=1, adjust="fdr", number=nrow(human_expression_sort))
  write.table(myresults,file=paste0("/tissue_specific_limma/human/myresults_",tissue_name[i],"_human.txt"),quote = F)
}




