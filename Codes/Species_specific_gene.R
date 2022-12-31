rm(list = ls())
setwd("~/Desktop/comparative_file/initial_dataset_cattle")
load("orthologous_cattle_human.Rdata")
load("integratedallgenes.Rdata")
load("human10830_cattle4866_17315gene.Rdata")
load("15696samples_17315genes.Rdata")


tissue<-unique(sorted(Meta_data_human$tissue_new))

for (i in 1:20){
  
  Tissue<-experiment[experiment$tissue_new==tissue[i],]
  Tissue_exp<-expression_inter[match(Tissue$sample,rownames(expression_inter)),]
  tissue_index<-which(Tissue$species=="Human")
  df<-c(rep("1",length(tissue_index)),rep("2",dim(Tissue)[1]-length(tissue_index)))
  table(df)
  
  design <- model.matrix(~-1+factor(df))
  colnames(design) <- c("Human","Cattle")
  table(design)
  design
  Tissue_exp<-t(Tissue_exp)
  Tissue_exp<-normalizeBetweenArrays(Tissue_exp)
  length(rownames(design))
  rownames(design)=colnames(Tissue_exp)
  design[tissue_index[1],]
  contrastmatrix <- makeContrasts(Human-Cattle, levels=design)
  #contrastmatrix <- makeContrasts(Cattle-Human, levels=design)
  fit <- lmFit(Tissue_exp, design) ##issue these commands to fit the model and make the contrasts
  fit2 <- contrasts.fit(fit, contrastmatrix)
  fit2 <- eBayes(fit2)
  fit2
  myresults <-topTable(fit2,coef=1, adjust="fdr", number=nrow(Tissue_exp))
  myresults[1:10,]
  write.table(myresults,file = paste0("/Users/yaoyuelin/Desktop/species_specific_limma/species_specific_",tissue[i],".txt"))
  
  DEG <- myresults
  
  logFC_cutoff <-log2(1.2)
  
  DEG$change = as.factor(ifelse(DEG$adj.P.Val < 0.05 & abs(DEG$logFC) > logFC_cutoff,
                                ifelse(DEG$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')
  )
  upgene<-DEG[DEG$change =='UP',]
  up_human_tissue<-rownames(upgene)
  downgene<-DEG[DEG$change =='DOWN',]
  downgene[1:10,]
  down<-Orthologous_human_cattle$`Cow gene stable ID`[match(rownames(downgene),Orthologous_human_cattle$`Gene stable ID`)]
  
  up_cattle_tissue<-down
  
  gene_sig<-rbind(upgene,downgene)
  table(gene_sig$change)
  
  
  write.table(up_human_tissue,file=paste0("~/species_specific/up_human_",tissue[i],".txt"),append = F,quote = F,sep = "/t",row.names = F,col.names = F)
  write.table(up_cattle_tissue,file=paste0("~/species_specific/up_cattle_",tissue[i],".txt"),append = F,quote = F,sep = "/t",row.names = F,col.names = F)
  
  table(gene_sig$change)
  
}




