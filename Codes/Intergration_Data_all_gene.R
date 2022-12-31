rm(list = ls()) 
library(randomcoloR)
library(ggplot2)
library(Rtsne)
library("corrplot")
library(Seurat)
library(dplyr)
library(patchwork)
library(pheatmap)


#load("human10830_cattle4866_17315gene.Rdata")
load("15696samples_17315genes.Rdata")
#load("orthologous_cattle_human.Rdata")
load("integratedallgenes.Rdata")

############################################
###########use Seurat to remove batch effect
############################################
table(Meta_data_cattle$tissue_new)
table(experiment$species)
expression<-t(expression+0.25)
expression<-CreateSeuratObject(expression)
expression@meta.data[["species"]]<-c(c(rep("Human",10830)),c(rep("Cattle",4866)))
expression@meta.data[["tissue_new"]]<-experiment$tissue_new
expression<-SplitObject(expression,split.by = "species")

expression<-expression[c("Human","Cattle")]
for (i in 1:length(expression)) {
  expression[[i]] <- NormalizeData(expression[[i]], verbose = FALSE)
  expression[[i]] <- FindVariableFeatures(expression[[i]], selection.method = "vst", 
                                          nfeatures = 17315, verbose = FALSE)
}

expression<-FindIntegrationAnchors(object.list =expression, dims = 1:30,anchor.features = 17315)
expression<-IntegrateData(anchorset = expression, dims = 1:30)


############# extract the expression matrix and save
expression_inter<-expression@assays[["integrated"]]@data
expression_inter<-t(as.matrix(expression_inter))
dim(expression_inter)
#save(expression_inter,experiment,file="integratedallgenes.Rdata")


#####################################################
########dimension reduction tsne ####################
#####################################################
expression_inter<-t(apply(expression_inter, MARGIN = 1, FUN = scale)) ##scale 
expression_tsne <- Rtsne(expression_inter,dims = 2, perplexity=150, theta=0.5, 
                         verbose=TRUE, max_iter = 1000,check_duplicates = FALSE,partial_pca = T,num_threads=50)

##########################################################
####################visualise tsne #######################
##########################################################
tiff(file = "Cattle_Human_TSNE1_all_gene.tiff",##reqiured to change
     res = 300, width = 2000, height = 1500,compression = "lzw")
par(mar=c(2,2,2,2))
info.norm<-as.data.frame(expression_tsne$Y)
colnames(info.norm)<-c("tsne1","tsne2")
colors<-c("red","blue")
#colors<-distinctColorPalette(length(unique(sort(experiment$species))))
names(colors)=unique(sort(experiment$species))
ggplot(info.norm, aes(x = tsne1, y = tsne2, colour = experiment$species)) + 
  geom_point(alpha = 1,size=3) + theme_bw()+
  xlab("t-SNE-1")+ylab("t-SNE-2")+
  theme(panel.background=element_blank(),legend.title=element_blank(),legend.text = element_text(size=20),
        panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
        axis.title.y=element_text(size=20,color="black"),
        axis.title.x=element_text(size=20,color="black"),
        axis.text.y=element_text(size=20,color="black"),
        axis.text.x=element_text(size=20,color="black",hjust=0.95,vjust=0.2))+
  scale_colour_manual(name="Tissue_categories",values=colors)#+stat_ellipse(type="norm")
dev.off()

#####################################
####use base R to plot ##############
#####################################
colors<-c("#00008B","#008B8B","#8B008B","#8B0000","#828282","#00B2EE","#0000FF","#FF0000","#FF00FF","#FF8C00","#FFFF00","#006400","#00FF00",
          "#CD950C","#FF6A6A","#CCCC33","pink","#FF1493","#8B7355","#00F5FF")
names(colors)=unique(sort(experiment$tissue_new))
par(mar=c(5, 5, 3, 14.1), xpd=TRUE)
plot(expression_tsne$Y, t='n', main="",xlab="t-SNE-1",ylab="t-SNE-2")
text(expression_tsne$Y,  labels=experiment$tissue_details, col=colors[experiment$tissue_new],cex=0.3) #,
legend("topright", inset=c(-0.5,0),
       legend=unique(sort(experiment$tissue_new)),
       col=colors[unique(sort(experiment$tissue_new))],
       pt.cex=3,
       cex=0.8,
       text.col="black",
       ncol=1,
       pch="-",title="Tissue categories")

###########################
tiff(file = "Cattle_Human_TSNE2_all_gene.tiff",##reqiured to change
     res = 300, width = 3000, height = 1500,compression = "lzw")
par(mar=c(2,2,5,5))
info.norm<-as.data.frame(expression_tsne$Y)
colnames(info.norm)<-c("tsne1","tsne2")
colors<-c("#00008B","#008B8B","#8B008B","#8B0000","#828282","#00B2EE","#0000FF","#FF0000","#FF00FF","#FF8C00","#FFFF00","#006400","#00FF00",
          "#CD950C","#FF6A6A","#CCCC33","pink","#FF1493","#8B7355","#00F5FF")
#colors<-distinctColorPalette(length(unique(sort(experiment$tissue))))
names(colors)=unique(sort(experiment$tissue_new))
ggplot(info.norm, aes(x = tsne1, y = tsne2, colour = experiment$tissue_new)) + 
  geom_point(alpha = 0.5,size=3) + theme_bw()+
  xlab("t-SNE-1")+ylab("t-SNE-2")+
  theme(panel.background=element_blank(),legend.title=element_blank(),legend.text = element_text(size=20),
        panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
        axis.title.y=element_text(size=20,color="black"),
        axis.title.x=element_text(size=20,color="black"),
        axis.text.y=element_text(size=20,color="black"),
        axis.text.x=element_text(size=20,color="black",hjust=0.95,vjust=0.2))+
  scale_colour_manual(name="Tissue_categories",values=colors,guide = guide_legend(ncol=2))#+stat_ellipse(type="norm")
dev.off()


##############################################
##############visualize using heatmap#########
##############################################
###1. correlation heatmap for median
expression_inter[1:10,1:10]
experiment$annotation<-paste(experiment$tissue_new,experiment$species,sep=" - ")
head(experiment)
length(table(experiment$annotation))
median_tissue_tpm<-aggregate(expression_inter,list(as.factor(experiment$annotation)),median)
rownames(median_tissue_tpm)<-median_tissue_tpm[,1]
(median_tissue_tpm)[1:10,1:10]
corr_mat=cor(t(median_tissue_tpm[,-1]))
dim(corr_mat)
head(corr_mat[1:10,1:10])
pheatmap(corr_mat,file="correlation_heatmap_median.pdf",cluster_rows = T,cluster_cols = T,fontsize_row = 15,fontsize_col = 15,cellwidth = 13,cellheight = 13,clustering_distance_rows ="correlation", clustering_distance_cols = "correlation")


###2. correlation heatmap for MAD
table(experiment$annotation)
experiment$annotation<-paste(experiment$tissue_new,experiment$species,sep=" - ")
(expression_inter)[1:10,1:10]
table(experiment$sample==rownames(expression_inter))
length(table(experiment$annotation))
mad_tissue_tpm<-aggregate(expression_inter,list(as.factor(experiment$annotation)),mad)
mad_tissue_tpm[1:10,1:10]
dim(mad_tissue_tpm)
rownames(mad_tissue_tpm)<-mad_tissue_tpm[,1]
mad_tissue_tpm[1:10,1:10]
corr_mat=cor(t(mad_tissue_tpm[,-1]))
dim(corr_mat)
head(corr_mat[1:10,1:10])
pheatmap(corr_mat,file="correlation_heatmap_mad.pdf",cluster_rows = T,cluster_cols = T,fontsize_row = 15,fontsize_col = 15,cellwidth = 13,cellheight = 13,clustering_distance_rows ="correlation", clustering_distance_cols = "correlation")

########
#########mean
expression_inter[1:10,1:10]
experiment$annotation<-paste(experiment$tissue_new,experiment$species,sep=" - ")
head(experiment)
length(table(experiment$annotation))
median_tissue_tpm<-aggregate(expression_inter,list(as.factor(experiment$annotation)),mean)
rownames(median_tissue_tpm)<-median_tissue_tpm[,1]
(median_tissue_tpm)[1:10,1:10]
corr_mat=cor(t(median_tissue_tpm[,-1]))
dim(corr_mat)
head(corr_mat[1:10,1:10])
pheatmap(corr_mat,file="correlation_heatmap_mean.pdf",cluster_rows = T,cluster_cols = T,fontsize_row = 15,fontsize_col = 15,cellwidth = 13,cellheight = 13,clustering_distance_rows ="correlation", clustering_distance_cols = "correlation")




