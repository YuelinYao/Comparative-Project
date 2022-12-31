rm(list = ls())

load("orthologous_cattle_human.Rdata")
load("overlap_geneid.Rdata")

library(ggpubr)
tissue_list<-as.data.frame(table(overlap_gene$tissue))[,1]
tissue_list
human_file<-dir("~/tissue_specific_limma/human",full.names = T)
cattle_file<-dir("~/tissue_specific_limma/cattle",full.names = T)
human_file
humantissue<-substr(human_file,nchar("~/tissue_specific_limma/human/myresults_")+1,nchar(human_file)-nchar("_human.txt"))
humantissue
tissue_list
cattle_file
cattletissue<-substr(cattle_file,nchar("~/tissue_specific_limma/cattle/myresults_")+1,nchar(cattle_file)-nchar(".txt"))
cattletissue

logFC_cutoff=1.5
i=1
for(i in 1:length(human_file)){
  human1 <- read.table(human_file[i],header = T,stringsAsFactors = F)
  human1$change = as.factor(ifelse(human1$adj.P.Val < 0.05 & abs(human1$logFC) > logFC_cutoff,
                                   ifelse(human1$logFC > logFC_cutoff ,'UP','DOWN'),'NOT'))
  
  human_up_genes<-rownames(human1)[human1$change=="UP"]

  
  cattle1<- read.table(cattle_file[i],header = T,stringsAsFactors = F)
  cattle1$change = as.factor(ifelse(cattle1$adj.P.Val < 0.05 & abs(cattle1$logFC) > logFC_cutoff,
                                    ifelse(cattle1$logFC > logFC_cutoff ,'UP','DOWN'),'NOT'))
  
  cattle_up_genes<-rownames(cattle1)[cattle1$change=="UP"]
  cattle_up_gene_trans<-Orthologous_human_cattle$`Gene stable ID`[match(cattle_up_genes,Orthologous_human_cattle$`Cow gene stable ID`)]
  length(cattle_up_gene_trans)
  overlap_up_genes<-intersect(cattle_up_gene_trans,human_up_genes)
  overlap_up_genes
  gene<- bitr(overlap_up_genes, fromType="ENSEMBL", toType = c("SYMBOL", "ENTREZID"), OrgDb=org.Hs.eg.db, drop = T)
  ego_BP <- enrichGO(gene = gene$ENTREZID,
                     OrgDb= org.Hs.eg.db,
                     ont = "BP",
                     pAdjustMethod = "BH",
                     minGSSize = 1,
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.05,
                     readable = TRUE)
  write.xlsx(ego_BP,file=paste0("~/Tissue-specific enrichment/overlap/",tissue_list[i],".xlsx"))
}











