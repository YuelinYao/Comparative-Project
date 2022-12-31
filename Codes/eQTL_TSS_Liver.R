load("orthologous_cattle_human.Rdata")
tissue_cattle<-read.table("fine_map_eQTL/Cattle_Liver.finemap.eQTL.txt")
colnames(tissue_cattle)<-c("Snp.id",	"Gene",	"findmapping_factor","Sample_number",	"shape1",	"shape2",	"dummy",	"ppval",	"bpval","dist",	"npval",	"slope",	"minor_allele_frequency",	"minor_allele_count",	"minor_allele_samples"   )
dim(tissue_cattle)
tissue_cattle[1:10,1:10]
length(unique(tissue_cattle$Gene))
tissue_human<-read.table("GTEx_Analysis_v8_eQTL_independent/Human_Liver.v8.independent_eqtls.txt",header = T) #4_Breast_Mammary_Tissue.v8.independent_eqtls.txt
length(unique(tissue_human$gene_id))
dim(tissue_human)


tissue_human<-tissue_human[tissue_human$gene_id%in%Orthologous_human_cattle$`Gene stable ID version`,]
dim(tissue_human)
tissue_cattle<-tissue_cattle[tissue_cattle$Gene%in%Orthologous_human_cattle$`Cow gene stable ID`,]
dim(tissue_cattle)


tissue_cattle<-tissue_cattle[tissue_cattle$npval<0.00001,]
tissue_human<-tissue_human[tissue_human$pval_nominal<0.00001,]
tissue_human$tss_distance
cattle_tss<-data.frame(tss=tissue_cattle$dist/1000,species="cattle")
human_tss<-data.frame(tss=tissue_human$tss_distance/1000,species="human")
all<-rbind(human_tss,cattle_tss)

cm<-mean(all$tss[all$species=="cattle"])
hm<-mean(all$tss[all$species=="human"])

hm
cm

colors<-c("red","blue")

library(ggplot2)
names(colors)=unique(all$species)
tiff(file = "tss_distance_Liver.tiff",##reqiured to change
     
     res = 300, width = 1000, height = 500,compression = "lzw")
p<-ggplot(all, aes(x = tss))
p + geom_density(aes(color = species))+ylab("Density")+ stat_density(aes(x=tss, colour=species),
                                                     geom="line",position="identity")+scale_colour_manual(name="Species",values=colors)+geom_vline(xintercept = c(hm,cm),col=c("red","blue"),linetype="dashed",size=0.4,show.legend = T)+theme_classic()+xlab("Distance to TSS (kb)")
dev.off()


