rm(list = ls())

load("orthologous_cattle_human.Rdata")
human_file<-dir("GTEx_Analysis_v8_eQTL_independent 2",full.names = T)
cattle_file<-dir("fine_map_eQTL 2",full.names = T)


humantissue<-substr(human_file,nchar("GTEx_Analysis_v8_eQTL_independent 2/")+1,nchar(human_file)-nchar(".v8.independent_eqtls.txt"))
cattletissue<-substr(cattle_file,nchar("fine_map_eQTL 2/")+1,nchar(cattle_file)-nchar(".finemap.eQTL.txt"))


Results <- array(NA,dim = c(length(human_file),length(cattle_file)))

rownames(Results)<-humantissue
colnames(Results)<-cattletissue

Results1 <- array(NA,dim = c(length(human_file),length(cattle_file)))


rownames(Results1)<-humantissue
colnames(Results1)<-cattletissue


Overlap_mat<-array(NA,dim = c(length(human_file),length(cattle_file)))
rownames(Overlap_mat)<-humantissue
colnames(Overlap_mat)<-cattletissue
human_egenes<-vector(length =length(human_file) )
names(human_egenes)<-humantissue
cattle_egenes<-vector(length =length(cattle_file) )
names(cattle_egenes)<-cattletissue


#######p-value
for(i in 1:length(human_file)){
  human1 <- read.table(human_file[i],header = T,stringsAsFactors = F)
  for(j in 1:length(cattle_file)){
    cattle1 <- read.table(cattle_file[j],header = T,stringsAsFactors = F)
    colnames(cattle1)<-c("Snp.id",	"Gene",	"findmapping_factor","Sample_number",	"shape1",	"shape2",	"dummy",	"ppval",	"bpval","dist",	"npval",	"slope",	"minor_allele_frequency",	"minor_allele_count",	"minor_allele_samples"   )
    
    tissue_cattle<-cattle1
    tissue_human<-human1
    tissue_human<-tissue_human[tissue_human$gene_id%in%Orthologous_human_cattle$`Gene stable ID version`,]
    tissue_cattle<-tissue_cattle[tissue_cattle$Gene%in%Orthologous_human_cattle$`Cow gene stable ID`,]
    
    #
    tissue_cattle<-tissue_cattle[tissue_cattle$npval<0.00001,]
    tissue_human<-tissue_human[tissue_human$pval_nominal<0.00001,]
    
    #
    human<-length(unique(tissue_human$gene_id))
    cattle<-length(unique(tissue_cattle$Gene))
    
    
    Orthologous<-Orthologous_human_cattle[Orthologous_human_cattle$`Cow gene stable ID`%in%tissue_cattle$Gene,]
    
    overlap<-intersect(tissue_human$gene_id,Orthologous$`Gene stable ID version`)
    Overlap<-length(overlap)
    
    phyper(Overlap-1 , human, 17315-human , cattle,lower.tail= FALSE)
    
    ########calculated 
    
    P_values <-     phyper(Overlap-1 , human, 17315-human , cattle,lower.tail= FALSE)
    Overlap_mat[i,j]<-Overlap
    Results[i,j] <- -log10(P_values) 
    human_egenes[i]<-human
    cattle_egenes[j]<-cattle
    Results1[i,j] <- P_values
    
    
  }
}




Mydata_raw_FDR <- p.adjust(Results1,method = "BH")
MyFDR_log <- -log10(Mydata_raw_FDR)
MyFDR_log <- matrix(MyFDR_log,nrow = dim(Results1)[1],byrow = F)


Mydata_raw_FDR[Mydata_raw_FDR<=0.05] <- "*"
Mydata_raw_FDR[Mydata_raw_FDR>0.05] <- " "
Mydata_raw_m <- matrix(Mydata_raw_FDR,nrow = dim(Results1)[1],byrow = F)




library(pheatmap)
breaksList = seq(0, 10, by = 1)
tiff(file = "0.00001p-value.tiff",##reqiured to change
     res = 300, width = 2500, height = 2000,compression = "lzw")
pheatmap(MyFDR_log,display_numbers = Mydata_raw_m,color = colorRampPalette(c("white","firebrick3"))(10),breaks = breaksList,cellwidth = 10,cellheight = 10)
dev.off()



humanbar<-human_egenes
hname<-substr(humantissue,nchar("Human_")+1,nchar(humantissue))
names(humanbar)<-hname


col = c("#00008B","#8B008B","#8B0000","#828282","grey","#FF0000",
        "#FF00FF","#FF8C00","#FFFF00","#006400","#00FF00","#CD950C",
        "deeppink","#FF6A6A","#CCCC33","#FF1493","#8B7355","#00F5FF") # color key

order1<-order(humanbar,decreasing = F)
humanbar<-humanbar[order(humanbar,decreasing = F)]


cattlebar<-cattle_egenes
cname<-substr(cattletissue,nchar("Cattle_")+1,nchar(cattletissue))

names(cattlebar)<-cname

col2<-c("#00008B","#8B008B","#8B0000","#CCCC33","navy","#FF1493",
        "purple","#FF0000","#FF00FF","deeppink","pink","#FF8C00","purple2",
        "#FFFF00","#006400","#00FF00","pink2","#CD950C","#FF6A6A","#8B7355","#00F5FF")

order2<-order(cattlebar,decreasing = F)
cattlebar<-cattlebar[order(cattlebar,decreasing = F)]


tiff("Number_of_eQTL_genes.tiff",res = 300, width = 2000, height = 1000,compression = "lzw")
par(mfrow=c(1,2),oma = c(0, 0, 3, 0),mar=c(6,4,1,1)+0.1,xpd=TRUE)
barplot(humanbar,main="Human",las=2,col=col[order1],cex.axis = 1,cex.names = 0.6,ylim=c(0,3000),ylab = "Number of eGene")
barplot(cattlebar,main="Cattle",las=2,col=col2[order2],cex.axis = 1,cex.names = 0.6,ylim=c(0,8000),ylab = "Number of eGene")
dev.off()





