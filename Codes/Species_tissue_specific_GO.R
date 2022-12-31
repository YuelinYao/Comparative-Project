library(clusterProfiler)
library(org.Bt.eg.db)
library(xlsx)
library(org.Hs.eg.db)
install.packages("org.Bt.eg.db")
BiocManager::install("org.Bt.eg.db")
load("orthologous_cattle_human.Rdata")
logFC_cutoff=1.5
human_file<-dir("~Comparative/tissue_specific_limma/human",full.names = T)
human_file
humantissue<-substr(human_file,nchar("~Comparative/tissue_specific_limma/human/myresults_")+1,nchar(human_file)-nchar("_human.txt"))
humantissue

cattle_file<-dir("~Comparative/tissue_specific_limma/cattle",full.names = T)
cattle_file
cattletissue<-substr(cattle_file,nchar("~Comparative/tissue_specific_limma/cattle/myresults_")+1,nchar(cattle_file)-nchar(".txt"))
cattletissue
i=1
for( i in 18:20){
  cattle1<- read.table(cattle_file[i],header = T,stringsAsFactors = F)
  cattle1$change = as.factor(ifelse(cattle1$adj.P.Val < 0.05 & abs(cattle1$logFC) > logFC_cutoff,
                                    ifelse(cattle1$logFC > logFC_cutoff ,'UP','DOWN'),'NOT'))
  
  cattle_up_genes<-rownames(cattle1)[cattle1$change=="UP"]
  length(cattle_up_genes)
  human1<- read.table(human_file[i],header = T,stringsAsFactors = F)
  human1$change = as.factor(ifelse(human1$adj.P.Val < 0.05 & abs(human1$logFC) > logFC_cutoff,
                                   ifelse(human1$logFC > logFC_cutoff ,'UP','DOWN'),'NOT'))
  
  human_up_genes<-rownames(human1)[human1$change=="UP"]
  human_up_genes
  
  length(human_up_genes)
  
  cattle_up_gene_trans<-Orthologous_human_cattle$`Gene stable ID`[match(cattle_up_genes,Orthologous_human_cattle$`Cow gene stable ID`)]
  cattle_up_genes_specific<-setdiff(cattle_up_gene_trans,human_up_genes)
  length(cattle_up_genes_specific)
  
  cattle_up_genes_specific<-Orthologous_human_cattle$`Cow gene stable ID`[match(cattle_up_genes_specific,Orthologous_human_cattle$`Gene stable ID`)]
  
  
  gene_cattle<- bitr(cattle_up_genes_specific, fromType="ENSEMBL", toType = c("SYMBOL", "ENTREZID"), OrgDb=org.Bt.eg.db, drop = T)
  ego_BP_cattle <- enrichGO(gene = gene_cattle$ENTREZID,
                            OrgDb= org.Bt.eg.db,
                            ont = "BP",
                            pAdjustMethod = "BH",
                            minGSSize = 1,
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.05,
                            readable = TRUE)
  
  
  
  human_up_specific<-setdiff(human_up_genes,cattle_up_gene_trans)
  
  gene_human<- bitr(human_up_specific, fromType="ENSEMBL", toType = c("SYMBOL", "ENTREZID"), OrgDb=org.Hs.eg.db, drop = T)
  ego_BP_human <- enrichGO(gene = gene_human$ENTREZID,
                           OrgDb= org.Hs.eg.db,
                           ont = "BP",
                           pAdjustMethod = "BH",
                           minGSSize = 1,
                           pvalueCutoff = 0.05,
                           qvalueCutoff = 0.05,
                           readable = TRUE)
  write.xlsx(ego_BP_cattle,file=paste0("~/species_specific_tissue_specific genes/cattle-specific/cattle_",cattletissue[i],".xlsx"),showNA=TRUE)
  write.xlsx(ego_BP_human,file=paste0("~/species_specific_tissue_specific genes/human-specific/human_",humantissue[i],".xlsx"),showNA=TRUE)
}

