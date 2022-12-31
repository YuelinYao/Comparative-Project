#############


############
Blood_cattle<-read.table("fine_map_eQTL/Cattle_Blood.finemap.eQTL.txt")
colnames(Blood_cattle)<-c("Snp.id",	"Gene",	"findmapping_factor","Sample_number",	"shape1",	"shape2",	"dummy",	"ppval",	"bpval","dist",	"npval",	"slope",	"minor_allele_frequency",	"minor_allele_count",	"minor_allele_samples"   )
length(unique(Blood_cattle$Gene))
Blood_human<-read.table("GTEx_Analysis_v8_eQTL_independent/Human_Whole_Blood.v8.independent_eqtls.txt",header = T)

length(unique(Blood_human$gene_id))
length(unique(Blood_cattle$Gene))
colnames(Blood_human)

Blood_human<-Blood_human[Blood_human$gene_id%in%Orthologous_human_cattle$`Gene stable ID version`,]
Blood_cattle<-Blood_cattle[Blood_cattle$Gene%in%Orthologous_human_cattle$`Cow gene stable ID`,]

Blood_cattle_egene<-Blood_cattle[Blood_cattle$npval<0.00001,]
Blood_human_egene<-Blood_human[Blood_human$pval_nominal<0.00001,]

Blood_cattle_egene_none<-Blood_cattle[!Blood_cattle$npval<0.00001,]
Blood_human_egene_none<-Blood_human[!Blood_human$pval_nominal<0.00001,]


length(unique(Blood_human_egene$gene_id))
length(unique(Blood_cattle_egene$Gene))

Orthologous<-Orthologous_human_cattle[Orthologous_human_cattle$`Gene stable ID version`%in%Blood_human_egene$gene_id,]


overlap<-intersect(Blood_cattle_egene$Gene,Orthologous$`Cow gene stable ID`)
Orthologous<-Orthologous[Orthologous$`Cow gene stable ID`%in%overlap,]

dim(Orthologous)

Blood_human_unique<-Blood_human_egene[order(Blood_human_egene$pval_nominal,decreasing = F),]
Blood_human_unique<-Blood_human_unique[!duplicated(Blood_human_unique$gene_id),]
Blood_human_overlap<-Blood_human_unique[match(Orthologous$`Gene stable ID version`,Blood_human_unique$gene_id),]
Blood_human_overlap[1:10,1:15]

Blood_human_specific<-Blood_human_unique[!Blood_human_unique$gene_id%in%Orthologous$`Gene stable ID version`,]
length(unique(Blood_human_specific$gene_id))


Blood_cattle_unique<-Blood_cattle_egene[order(Blood_cattle_egene$npval,decreasing = F),]
Blood_cattle_unique<-Blood_cattle_unique[!duplicated(Blood_cattle_unique$Gene),]

Blood_cattle_overlap<-Blood_cattle_unique[match(Orthologous$`Cow gene stable ID`,Blood_cattle_unique$Gene),]

Blood_cattle_specific<-Blood_cattle_unique[!Blood_cattle_unique$Gene%in%Orthologous$`Cow gene stable ID`,]






experiment$annotation<-paste(experiment$tissue_new,experiment$species,sep="-")
variance_tissue_tpm<-aggregate(expression,list(as.factor(experiment$annotation)),mad)
variance_tissue_tpm[1:10,1:10]
rownames(variance_tissue_tpm)<-variance_tissue_tpm[,1]
variance_tissue_tpm<-variance_tissue_tpm[,-1]
variance_tissue_tpm[,1:10]
mad_cattle<-variance_tissue_tpm["Blood/Immune-Cattle",]
mad_cattle[1,1:10]
mad_human<-variance_tissue_tpm["Blood/Immune-Human",]
mad_human[1,1:10]
mad_diff<-mad_cattle-mad_human
mad_diff[1,1:10]


Blood <- t(mad_diff)
Blood<-as.data.frame(Blood)
colnames(Blood)<-"mad_diff"
dim(Orthologous)

both<-Blood[rownames(Blood)%in%Orthologous$`Gene stable ID`,]



cattle<-Blood_cattle_specific$Gene
cattle<-Orthologous_human_cattle$`Gene stable ID`[Orthologous_human_cattle$`Cow gene stable ID`%in%cattle]
cattle<-Blood[rownames(Blood)%in%cattle,]


human<-Blood_human_specific$gene_id
human<-Orthologous_human_cattle$`Gene stable ID`[Orthologous_human_cattle$`Gene stable ID version`%in%human]
human<-Blood[rownames(Blood)%in%human,]


all<-c(rownames(human),rownames(cattle),rownames(both))
neither<-Orthologous_human_cattle$`Gene stable ID`[!Orthologous_human_cattle$`Gene stable ID`%in%all]
neither<-Blood[rownames(Blood)%in%neither,]


cattle<-data.frame(diff=cattle ,class="cattle")
human<-data.frame(diff=human,class="human")
both<-data.frame(diff=both,class="both")
neither<-data.frame(diff=neither,class="neither")



all<-rbind(cattle,human,both,neither)
summary(all$diff)
my_comparisons <- list( c("both", "neither"), c("cattle", "neither"),c("neither","human"),c("cattle","human") )
hist(all$diff)
library(ggpubr)
library(ggplot2)
colors<-c("blue","green","orange","purple")
tiff(file = "eQTL_MAD.tiff",##reqiured to change
     res = 300, width = 1200, height = 600,compression = "lzw")
p<-ggplot(all, aes(x = diff,colour=class))+ylab("Cumulative frequency")+xlab("Difference in MAD")
p + stat_ecdf(geom = "step")+scale_colour_manual(name="eGene classification",values=colors)+guides(colour = guide_legend(ncol = 2))+theme_classic()+xlim(-10,10)
dev.off()




cattle_neither<-wilcox.test(all$diff[all$class=="cattle"],all$diff[all$class=="neither"],alternative = "greater")
cattle_neither
human_neither<-wilcox.test(all$diff[all$class=="human"],all$diff[all$class=="neither"],alternative = "less")
human_neither
both_cattle<-wilcox.test(all$diff[all$class=="both"],all$diff[all$class=="cattle"],alternative = "less")
both_cattle
both_human<-wilcox.test(all$diff[all$class=="both"],all$diff[all$class=="human"],alternative = "greater")
both_human
both_neither<-wilcox.test(all$diff[all$class=="both"],all$diff[all$class=="neither"],alternative = "greater")
both_neither


cattle_human<-wilcox.test(all$diff[all$class=="cattle"],all$diff[all$class=="human"],alternative = "greater")
cattle_human