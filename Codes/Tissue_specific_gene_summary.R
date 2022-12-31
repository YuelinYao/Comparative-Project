rm(list = ls())
load("orthologous_cattle_human.Rdata")

####for cattle
cattle_file<-dir("/Volumes/Seagate Exp/desktop/Comparative/tissue_specific_limma/cattle",full.names = T)
cattle_file

cattle_tissue<-substr(cattle_file,nchar("/Volumes/Seagate Exp/desktop/Comparative/tissue_specific_limma/cattle/myresults_")+1,nchar(cattle_file)-nchar(".txt"))
length(cattle_tissue)

logFC_cutoff<-1.5
all_gene_cattle<-data.frame(gene=NA,tissue=NA,species="Cattle")
top10_gene_cattle<-data.frame(gene=NA,tissue=NA,species="Cattle")

for (i in 1:20){
  
  myresults<-read.table(cattle_file[i],header = T,stringsAsFactors = F) 

  DEG <- myresults
  DEG$change = as.factor(ifelse(DEG$adj.P.Val < 0.05 & abs(DEG$logFC) > logFC_cutoff,
                                ifelse(DEG$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')
  )
  
  upgene<-DEG[DEG$change =='UP',]
  upgene
  upgene_id<-rownames(upgene)
  upgene_id
  tissue_g<-data.frame(gene=upgene_id,tissue=cattle_tissue[i],species="Cattle")
  tissue_g
  all_gene_cattle<-rbind(all_gene_cattle,tissue_g)
  top10<-data.frame(gene=upgene_id[1:10],tissue=cattle_tissue[i],species="Cattle")
  top10_gene_cattle<-rbind(top10_gene_cattle,top10)
  
}

top10_gene_cattle<-top10_gene_cattle[-1,]
all_gene_cattle<-all_gene_cattle[-1,]

#save(all_gene_cattle,file="up_gene_cattle.Rdata")




#####for human
human_file<-dir("/Volumes/Seagate Exp/desktop/Comparative/tissue_specific_limma/human",full.names = T)
human_file

human_tissue<-substr(human_file,nchar("/Volumes/Seagate Exp/desktop/Comparative/tissue_specific_limma/human/myresults_")+1,nchar(human_file)-nchar("_human.txt"))
length(human_tissue)
human_tissue
logFC_cutoff<-1.5
all_gene_human<-data.frame(gene=NA,tissue=NA,species="Human")
top10_gene_human<-data.frame(gene=NA,tissue=NA,species="Human")

for (i in 1:20){
  
  myresults<-read.table(human_file[i],header = T,stringsAsFactors = F) 
  
  DEG <- myresults
  DEG$change = as.factor(ifelse(DEG$adj.P.Val < 0.05 & abs(DEG$logFC) > logFC_cutoff,
                                ifelse(DEG$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')
  )
  
  upgene<-DEG[DEG$change =='UP',]
  upgene
  upgene_id<-rownames(upgene)
  upgene_id
  tissue_g<-data.frame(gene=upgene_id,tissue=human_tissue[i],species="Human")
  tissue_g
  all_gene_human<-rbind(all_gene_human,tissue_g)
  top10<-data.frame(gene=upgene_id[1:10],tissue=human_tissue[i],species="Human")
  top10_gene_human<-rbind(top10_gene_human,top10)
  
}

top10_gene_human<-top10_gene_human[-1,]
all_gene_human<-all_gene_human[-1,]

#save(all_gene_human,file="up_gene_human.Rdata")
#save(top10_gene_cattle,top10_gene_human,file="top10_gene.Rdata")



###test overlap
all_gene_cattle
all_gene_human

tissue_name<-unique(all_gene_cattle$tissue)
results<-array(data=NA,dim=c(20,4))
rownames(results)<-tissue_name
colnames(results)<-c("Human","Cattle","Overlap","p_value")
overlap_gene<-data.frame(gene=NA,tissue=NA)
i=1
for (i in 1:length(tissue_name)){
  human_gene<-all_gene_human$gene[all_gene_human$tissue==tissue_name[i]]
  cattle_gene<-all_gene_cattle$gene[all_gene_cattle$tissue==tissue_name[i]]
  
  human_gene_trans<-Orthologous_human_cattle$`Cow gene stable ID`[Orthologous_human_cattle$`Gene stable ID`%in%human_gene]
  overlap<-intersect(human_gene_trans,cattle_gene)
  overlap_g<-data.frame(gene=overlap,tissue=tissue_name[i])
  overlap_gene<-rbind(overlap_gene,overlap_g)
  results[i,c(1,2,3)]<-c(length(human_gene),length(cattle_gene),length(overlap))
  results[i,4]<-phyper(length(overlap)-1 , length(human_gene), 17315-length(human_gene) ,length(cattle_gene),lower.tail= FALSE)
  
  
}
results

results<-as.data.frame(results)
results


