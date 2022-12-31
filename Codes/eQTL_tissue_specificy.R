rm(list = ls())

load("orthologous_cattle_human.Rdata")
human_file<-dir("GTEx_Analysis_v8_eQTL_independent",full.names = T)
cattle_file<-dir("fine_map_eQTL",full.names = T)
human_file
cattle_file
humantissue<-substr(human_file,nchar("GTEx_Analysis_v8_eQTL_independent/")+1,nchar(human_file)-nchar(".v8.independent_eqtls.txt"))
humantissue

cattletissue<-substr(cattle_file,nchar("fine_map_eQTL/")+1,nchar(cattle_file)-nchar(".finemap.eQTL.txt"))
cattletissue
cattle_egenes<-data.frame(geneID = NA,tissue=NA)
cattle_egenes
for( i in 1: length(cattle_file)){
  cattle1 <- read.table(cattle_file[i],header = T,stringsAsFactors = F)
  colnames(cattle1)<-c("Snp.id",	"Gene",	"findmapping_factor","Sample_number",	"shape1",	"shape2",	"dummy",	"ppval",	"bpval","dist",	"npval",	"slope",	"minor_allele_frequency",	"minor_allele_count",	"minor_allele_samples"   )
  
  tissue_cattle<-cattle1
  tissue_cattle<-tissue_cattle[tissue_cattle$Gene%in%Orthologous_human_cattle$`Cow gene stable ID`,]
  tissue_cattle<-tissue_cattle[tissue_cattle$npval<0.00001,]
  tissue_cattle<-data.frame(geneID=unique(tissue_cattle$Gene),tissue=i)
  cattle_egenes<-rbind(cattle_egenes,tissue_cattle)
}

table(cattle_egenes$tissue)
cattle_gene_count<-table(cattle_egenes$geneID)
cattle_gene_count<-as.data.frame(cattle_gene_count)
cattle_eqtl<-table(cattle_gene_count$Freq)
cattle_eqtl



#####
human_egenes<-data.frame(geneID = NA,tissue=NA)
human_egenes
for( i in 1: length(human_file)){
  human1 <- read.table(human_file[i],header = T,stringsAsFactors = F)
  tissue_human<-human1
  tissue_human<-tissue_human[tissue_human$gene_id%in%Orthologous_human_cattle$`Gene stable ID version`,]
  tissue_human<-tissue_human[tissue_human$pval_nominal<0.00001,]
  tissue_human<-data.frame(geneID=unique(tissue_human$gene_id),tissue=i)
  human_egenes<-rbind(human_egenes,tissue_human)
}

human_gene_count<-table(human_egenes$geneID)

human_gene_count<-as.data.frame(human_gene_count)



human_eqtl<-table(human_gene_count$Freq)
human_eqtl
as.matrix(human_eqtl)



dim(human_gene_count)

dim(cattle_gene_count)
Orthologous<-Orthologous_human_cattle[Orthologous_human_cattle$`Cow gene stable ID`%in%cattle_gene_count$Var1,]
overlap<-intersect(human_gene_count$Var1,Orthologous$`Gene stable ID version`)
Orthologous<-Orthologous[Orthologous$`Gene stable ID version`%in%overlap,]
tissue_human_overlap<-human_gene_count[match(Orthologous$`Gene stable ID version`,human_gene_count$Var1),]
tissue_cattle_overlap<-cattle_gene_count[match(Orthologous$`Cow gene stable ID`,cattle_gene_count$Var1),]

cor.test(tissue_cattle_overlap$Freq,tissue_human_overlap$Freq,method = "spearman")
as.matrix(cattle_eqtl)

tiff(file="eQTL_tissue_distribution.tiff", width = 14, height = 7, units = "cm", res = 600, pointsize = 8,compression= "lzw")

par(mfrow=c(1,2),oma = c(0, 0, 3, 0))

barplot(cattle_eqtl,ylim = c(0,4000),col="Orange")
mtext('Number of eGenes in cattle', side = 2, cex = 1.3,  line = 2)
mtext('Number of tissues in cattle', side = 1, cex = 1.3,  line = 2)
barplot(human_eqtl,ylim = c(0,600),col="PaleGreen")
mtext('Number of eGenes in human', side = 2, cex = 1.3,  line = 2)
mtext('Number of tissues in human', side = 1, cex = 1.3,  line = 2)
mtext("Spearman correlation : 0.15, P-value: 1.909e-13", side = 3, line = 0, outer = T,font = 4,cex = 1.3)
dev.off()
