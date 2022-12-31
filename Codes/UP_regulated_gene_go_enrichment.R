library(clusterProfiler)
library(org.Bt.eg.db)
library(xlsx)
library(org.Hs.eg.db)
load("orthologous_cattle_human.Rdata")

human_file<-dir("~/species_specific/human",full.names = T)
human_file
humantissue<-substr(human_file,nchar("~/species_specific/human/up_human_")+1,nchar(human_file)-nchar(".txt"))
humantissue

cattle_file<-dir("~/species_specific/cattle",full.names = T)
cattle_file
cattletissue<-substr(cattle_file,nchar("~/species_specific/cattle/up_cattle_")+1,nchar(cattle_file)-nchar(".txt"))
cattletissue
i=1
 
for( i in 1:20){
  cattle1<- read.table(cattle_file[i],header = F,stringsAsFactors = F)
  dim(cattle1)
  
  human1<- read.table(human_file[i],header = F,stringsAsFactors = F)

  gene_cattle<- bitr(cattle1$V1, fromType="ENSEMBL", toType = c("SYMBOL", "ENTREZID"), OrgDb=org.Bt.eg.db, drop = T)
  ego_BP_cattle <- enrichGO(gene = gene_cattle$ENTREZID,
                            OrgDb= org.Bt.eg.db,
                            ont = "BP",
                            pAdjustMethod = "BH",
                            minGSSize = 1,
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.05,
                            readable = TRUE)
  
  
  
  gene_human<- bitr(human1$V1, fromType="ENSEMBL", toType = c("SYMBOL", "ENTREZID"), OrgDb=org.Hs.eg.db, drop = T)
  ego_BP_human <- enrichGO(gene = gene_human$ENTREZID,
                           OrgDb= org.Hs.eg.db,
                           ont = "BP",
                           pAdjustMethod = "BH",
                           minGSSize = 1,
                           pvalueCutoff = 0.05,
                           qvalueCutoff = 0.05,
                           readable = TRUE)
  write.xlsx(ego_BP_cattle,file=paste0("~/species_specific/cattle_go/cattle_",cattletissue[i],".xlsx"),showNA=TRUE)
  write.xlsx(ego_BP_human,file=paste0("~/species_specific/human_go/human_",humantissue[i],".xlsx"),showNA=TRUE)
}

