rm(list = ls()) 
library(randomcoloR)
library(ggplot2)
library(Rtsne)
library("corrplot")
library(Seurat)
library(factoextra)
library(gridExtra)
library(grid)

## sample summary table
sample_summary<-table(experiment$tissue_new,experiment$species)
sample_summary<-t(sample_summary)
sample_summary
pdf("sample_summary.pdf",width = 20,height = 10)
grid.table(sample_summary)
dev.off()

#cattle_mapping <- read.table("data_info_tissue_breed_with_mapping_rate.txt",header = T,sep="\t",stringsAsFactors = F) 
##################
#cattle_mapping$Mapping_rate[1:10]

load("human10830_cattle4866_17315gene.Rdata")
load("15696samples_17315genes.Rdata")
load("orthologous_cattle_human.Rdata")


#save(Orthologous_human_cattle,Meta_data_cattle,Cattle_expression,Meta_data_human,human_expression,file="human10830_cattle4866_17315gene.Rdata")
#save(expression,experiment,Orthologous_human_cattle,file="15696samples_17315genes.Rdata")



#check the order
table(rownames(expression)==experiment$sample)
table(rownames(Cattle_expression)==Meta_data_cattle$Sample)
table(rownames(human_expression)==Meta_data_human$SAMPID)
table(colnames(Cattle_expression)==Orthologous_human_cattle$`Cow gene stable ID` )



Cattle_expression<-t(apply(log(Cattle_expression+0.25), MARGIN = 1, FUN = scale)) 
colnames(Cattle_expression)<-Orthologous_human_cattle$`Cow gene stable ID`
Cattle_expression[1:10,1:10]

tsne_cattle <- Rtsne(Cattle_expression,dims = 2, perplexity=100, theta=0.5, verbose=TRUE, max_iter = 1000,check_duplicates = FALSE,partial_pca = T,num_threads=50)

pdf("cattle_sample.pdf",width=9,height=6)
colors<-c("#00008B","#008B8B","#8B008B","#8B0000","#828282","#00B2EE","#0000FF","#FF0000","#FF00FF","#FF8C00","#FFFF00","#006400","#00FF00",
          "#CD950C","#FF6A6A","#CCCC33","pink","#FF1493","#8B7355","#00F5FF")
names(colors)=unique(sort(Meta_data_cattle$tissue_new))
par(mar=c(5, 5, 3, 14.1), xpd=TRUE)
plot(tsne_cattle$Y, t='n', main="",xlab="t-SNE-1",ylab="t-SNE-2",cex.lab=2)
text(tsne_cattle$Y,  labels=Meta_data_cattle$Details.of.tissue , col=colors[Meta_data_cattle$tissue_new],cex=0.3) #,
legend("topright", inset=c(-0.5,0),
       legend=unique(sort(Meta_data_cattle$tissue_new)),
       col=colors[unique(sort(Meta_data_cattle$tissue_new))],
       pt.cex=4,
       cex=1,
       text.col="black",
       ncol=1,
       pch="-",title="Tissue categories")
dev.off()


library(ggplot2)
info.norm<-as.data.frame(tsne_cattle$Y)
colnames(info.norm)<-c("tsne1","tsne2")
pdf("cattle-clustering.pdf",width=9,height=6)
colors<-c("#00008B","#008B8B","#8B008B","#8B0000","#828282","#00B2EE","#0000FF","#FF0000","#FF00FF","#FF8C00","#FFFF00","#006400","#00FF00",
          "#CD950C","#FF6A6A","#CCCC33","pink","#FF1493","#8B7355","#00F5FF")
names(colors)=unique(sort(Meta_data_cattle$tissue_new))
ggplot(info.norm, aes(x = tsne1, y = tsne2, colour = Meta_data_cattle$tissue_new)) + 
  geom_point(alpha = 0.5,size=2) + theme_bw()+
  xlab("t-SNE-1")+ylab("t-SNE-2")+
  theme(panel.background=element_blank(),legend.title=element_blank(),
        panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
        axis.title.y=element_text(size=20,color="black"),
        axis.title.x=element_text(size=20,color="black"),
        axis.text.y=element_text(size=20,color="black"),
        axis.text.x=element_text(size=20,color="black",hjust=0.95,vjust=0.2))+
  scale_colour_manual(name="Tissue_categories",values=colors)#+stat_ellipse()
dev.off()

##########

human_expression<-t(apply(log(human_expression+0.25), MARGIN = 1, FUN = scale))
colnames(human_expression)<-Orthologous_human_cattle$`Gene stable ID`
tsne_human <- Rtsne(human_expression,dims = 2, perplexity=100, theta=0.5, verbose=TRUE, max_iter = 1000,check_duplicates = FALSE,partial_pca = T,num_threads=50)


pdf("human.pdf",width=9,height=6)
colors<-c("#00008B","#008B8B","#8B008B","#8B0000","#828282","#00B2EE","#0000FF","#FF0000","#FF00FF","#FF8C00","#FFFF00","#006400","#00FF00",
          "#CD950C","#FF6A6A","#CCCC33","pink","#FF1493","#8B7355","#00F5FF")
names(colors)=unique(sort(Meta_data_human$tissue_new))
par(mar=c(5, 5, 3, 14.1), xpd=TRUE)
plot(tsne_human$Y, t='n', main="",xlab="t-SNE-1",ylab="t-SNE-2")
text(tsne_human$Y,  labels=Meta_data_human$SMTSD, col=colors[Meta_data_human$tissue_new],cex=0.3) #,
legend("topright", inset=c(-0.5,0),
       legend=unique(sort(Meta_data_human$tissue_new)),
       col=colors[unique(sort(Meta_data_human$tissue_new))],
       pt.cex=3,
       cex=0.8,
       text.col="black",
       ncol=1,
       pch="-",title="Tissue categories")
dev.off()

pdf("human_dot.pdf",width=9,height=6)
info.norm<-as.data.frame(tsne_human$Y)
colnames(info.norm)<-c("tsne1","tsne2")
colors<-c("#00008B","#008B8B","#8B008B","#8B0000","#828282","#00B2EE","#0000FF","#FF0000","#FF00FF","#FF8C00","#FFFF00","#006400","#00FF00",
          "#CD950C","#FF6A6A","#CCCC33","pink","#FF1493","#8B7355","#00F5FF")
names(colors)=unique(sort(Meta_data_human$tissue_new))
ggplot(info.norm, aes(x = tsne1, y = tsne2, colour = Meta_data_human$tissue_new)) + 
  geom_point(alpha = 0.5,size=2) + theme_bw()+
  xlab("t-SNE-1")+ylab("t-SNE-2")+
  theme(panel.background=element_blank(),legend.title=element_blank(),
        panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
        axis.title.y=element_text(size=20,color="black"),
        axis.title.x=element_text(size=20,color="black"),
        axis.text.y=element_text(size=20,color="black"),
        axis.text.x=element_text(size=20,color="black",hjust=0.95,vjust=0.2))+
  scale_colour_manual(name="Tissue_categories",values=colors)#+stat_ellipse(type="norm")
dev.off()







