rm(list=ls())
library(WGCNA)

library(sva)
library(foreach)
library(doParallel)
library(BBmisc)
library(ggplot2)
library(Rtsne)


#load("orthologous_cattle_human.Rdata")
#load("integratedallgenes.Rdata")
#load("human10830_cattle4866_17315gene.Rdata")
load("15696samples_17315genes.Rdata")
experiment$annotation<-paste(experiment$tissue_new,experiment$species,sep="-")
table(experiment$tissue_new)
tissue_cattle<-experiment[experiment$annotation=="Muscle-Cattle",]
tissue_exp_cattle<-expression[match(tissue_cattle$sample,rownames(expression)),]

tissue_human<-experiment[experiment$annotation=="Muscle-Human",]
tissue_exp_human<-expression[match(tissue_human$sample,rownames(expression)),]

dim(tissue_cattle)
dim(tissue_human)
dim(tissue_exp_cattle)
dim(tissue_exp_human)
tissue_exp_human[1:10,1:10]

datExpr0<-tissue_exp_cattle



colnames(datExpr0)
#############
gsg = goodSamplesGenes(datExpr0, verbose = 3) #We first check for genes and samples with too many missing values:

gsg$allOK
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(colnames(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}
gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK


dim(datExpr0)
sampleTree = hclust(dist(datExpr0), method = "average");
#############3

# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "~/AdiposeClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
abline(h = 99000, col = "red")


#dev.off()

clust = cutreeStatic(sampleTree, cutHeight = 99000, minSize = 10)
table(clust)
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]
nGenes_cattle = ncol(datExpr)
nSamples_cattle = nrow(datExpr)
nGenes_cattle
nSamples_cattle
tissue_exp_cattle<-datExpr0
dim(tissue_exp_cattle)
##################3human
datExpr0<-tissue_exp_human
tissue_exp_human[1:10,1:10]
#############
gsg = goodSamplesGenes(datExpr0, verbose = 3) #We first check for genes and samples with too many missing values:


gsg$allOK
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(colnames(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}
gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK
dim(datExpr0)
datExpr0[1:10,1:10]
sampleTree = hclust(dist(datExpr0), method = "average");
#############3

# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)

par(cex = 0.6);
par(mar = c(0,4,2,0))

plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
abline(h = 58000, col = "red")


clust = cutreeStatic(sampleTree, cutHeight = 58000, minSize = 10)
table(clust)
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]
nGenes_cattle = ncol(datExpr)
nSamples_cattle = nrow(datExpr)

nGenes_cattle
nSamples_cattle

tissue_exp_human<-datExpr0

##########
dim(tissue_exp_human)
dim(tissue_exp_cattle)
tissue_exp_cattle[1:10,1:10]
tissue_exp_human[1:10,1:10]
l<-intersect(colnames(tissue_exp_human),colnames(tissue_exp_cattle))
length(l)
tissue_exp_cattle<-tissue_exp_cattle[,colnames(tissue_exp_cattle)%in%l]
tissue_exp_human<-tissue_exp_human[,colnames(tissue_exp_human)%in%l]


dim(tissue_exp_human)
dim(tissue_exp_cattle)
table(colnames(tissue_exp_human)==colnames(tissue_exp_cattle))
save(tissue_exp_human,tissue_exp_cattle,file="muscle-wgcna.Rdata")
####reclustering
load("muscle-wgcna.Rdata")

##############network construction and module detection
allowWGCNAThreads() 
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=1))
sft = pickSoftThreshold(tissue_exp_cattle, powerVector = powers, verbose = 5)
sft$powerEstimate ## 6
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
abline(h=0.8,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
#power=19
softPower_cattle = 6

#####human
powers = c(c(1:10), seq(from = 12, to=20, by=1))
sft = pickSoftThreshold(tissue_exp_human, powerVector = powers, verbose = 5)
sft$powerEstimate #6
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
abline(h=0.8,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
#power=19

softPower_human = 6
softPower_cattle = 6






#####cattle network
adjacency_cattle = abs(cor(tissue_exp_cattle ,use = "p"))^softPower_cattle
diag(adjacency_cattle) = 0

#human
adjacency_human = abs(cor(tissue_exp_human ,use = "p"))^softPower_human
diag(adjacency_human) = 0



## Calculation of the whole network connectivity k:


ConnectivityHuman <- apply(adjacency_human,1,sum,na.rm = TRUE)
Connectivitycattle <- apply(adjacency_cattle,1,sum,na.rm = TRUE)




## Depiction of scale-free topology:

# The black curve corresponds to scale-free topology and the red curve corresponds to #truncated scale-free topology.

par(mfrow = c(2,1))

scaleFreePlot(ConnectivityHuman, main  = paste("Human ,","power = 8"),truncated = T)
scaleFreePlot(Connectivitycattle, main  = paste("Cattle ,","power = 8"),truncated = T)


## Scaling k to lie between 0 and 1:

ConnectivityHuman = ConnectivityHuman/max(ConnectivityHuman)
Connectivitycattle = Connectivitycattle/max(Connectivitycattle)





# As a preprocessing step toward module construction, we restrict the network to genes 
#with reasonably high connectivity. This does not lead to a big loss of information since 
#module genes tend to have high connectivity. Toward this end, consider the median 
#connectivity in human and cattle:


median(ConnectivityHuman)
summary(ConnectivityHuman)
minconnections = 0.01
rest1 = Connectivitycattle>minconnections | ConnectivityHuman>minconnections
table(rest1)



AdjMatcattlerest = adjacency_cattle[rest1,rest1]
AdjMatHumanrest = adjacency_human[rest1,rest1]

#AdjMatcattlerest = adjacency_cattle
#AdjMatHumanrest = adjacency_human


#Module Construction

# The topological overlap of two nodes reflects their similarity in terms of the #commonality of the nodes they connect to, see ref. 6.

## Creating distance matrices based upon the topological overlap of genes for #humans and chimpanzees:

distTOMCattle <- TOMdist(AdjMatcattlerest)
distTOMHuman <- TOMdist(AdjMatHumanrest)



# To group genes with coherent expression profiles into modules, we use average linkage

# hierarchical clustering, which uses the topological overlap measure as dissimilarity.

## Performing average linkage hierarchical clustering using these distance matrices:

hierTOMHuman <- hclust(as.dist(distTOMHuman),method = "average")
hierTOMCattle <- hclust(as.dist(distTOMCattle),method = "average")
minModuleSize=200

par(mfrow = c(1,2))
plot(hierTOMHuman,labels = F,main = "Human")
plot(hierTOMCattle,labels = F,main = "Cattle")

## Assign colors to modules based upon the height cutoff (h1) and minimum module size #(minsize1):

# Once a dendrogram is obtained from a hierarchical clustering method, we need to #choose a height cutoff h1 to arrive at a clustering. It is a judgement call where to cut the tree

# branches. The height cut-off can be found by inspection: a height cutoff value is chosen #in the dendrogram such that some of the resulting branches correspond to the discrete #diagonal blocks


minModuleSize=200
dynamicMods = cutreeDynamic(dendro = hierTOMHuman, distM = distTOMHuman,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);


dynamicColors = labels2colors(dynamicMods)
length(dynamicColors)
table(dynamicColors)
nGenes = ncol(tissue_exp_human)
nSamples = nrow(tissue_exp_human)
dynamicMods_cattle = cutreeDynamic(dendro = hierTOMCattle, distM = distTOMCattle,
                                   deepSplit = 2, pamRespectsDendro = FALSE,
                                   minClusterSize = 300);


dynamicColors_cattle = labels2colors(dynamicMods_cattle)
table(dynamicColors_cattle)


sizeGrWindow(8,6)
tiff(file = "~/muscle_net1.tiff",##reqiured to change
     res = 300, width = 1200, height = 1200,compression = "lzw")
plotDendroAndColors(hierTOMHuman, dynamicColors, "Human colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Human network, Human colors")
dev.off()
tiff(file = "~/muscle_net2.tiff",##reqiured to change
     res = 300, width = 1200, height = 1200,compression = "lzw")
plotDendroAndColors(hierTOMCattle, dynamicColors, "Human colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Cattle network, Human colors")
dev.off()

tiff(file = "~/muscle_net3.tiff",##reqiured to change
     res = 300, width = 2000, height = 1000,compression = "lzw")

layout(matrix(c(1:4),2,2))
plotDendroAndColors(hierTOMHuman, dynamicColors, "Human colors",
                    dendroLabels = FALSE, hang = 0.03,cex.colorLabels = 1,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Human network, Human colors",setLayout = F)

plotDendroAndColors(hierTOMCattle, dynamicColors, "Human colors",
                    dendroLabels = FALSE, hang = 0.03,cex.colorLabels = 1,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Cattle network, Human colors",setLayout = F)
dev.off()







tiff(file = "~/muscle_net4.tiff",##reqiured to change
     res = 300, width = 1000, height = 1000,compression = "lzw")
plotDendroAndColors(hierTOMCattle, cbind(dynamicColors, dynamicColors_cattle), c("Human colors","Cattle corlors"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Cattle network, Cattle colors")
dev.off()




colorh = as.character(dynamicColors)




tsne_cattle <- Rtsne(as.dist(distTOMCattle),dims = 2, perplexity=100, theta=0.5, verbose=TRUE, max_iter = 1000,check_duplicates = FALSE,partial_pca = T,num_threads=50)


tiff(file = "~/cattle-clustering-muscle.tiff",##reqiured to change
     res = 300, width = 600, height = 600,compression = "lzw")
info.norm<-as.data.frame(tsne_cattle$Y)
colnames(info.norm)<-c("tsne1","tsne2")
ggplot(info.norm, aes(x = tsne1, y = tsne2)) + theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))+geom_point(alpha = 0.5,size=0.03,color=colorh) 
dev.off()


tsne_human <- Rtsne(as.dist(distTOMHuman),dims = 2, perplexity=100, theta=0.5, verbose=TRUE, max_iter = 1000,check_duplicates = FALSE,partial_pca = T,num_threads=50)

library(ggplot2)

tiff(file = "~/human-clustering-muscle.tiff",##reqiured to change
     res = 300, width = 600, height = 600,compression = "lzw")
info.norm<-as.data.frame(tsne_human$Y)
colnames(info.norm)<-c("tsne1","tsne2")
ggplot(info.norm, aes(x = tsne1, y = tsne2)) + theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))+
  geom_point(alpha = 0.5,size=0.03,color=colorh) 

dev.off()


MEs0 = moduleEigengenes(tissue_exp_human[,rest1], dynamicColors)$eigengenes
MEs0
human_expression<-(tissue_exp_human[,rest1])
dim(human_expression)
probe<-colnames(human_expression)
table(dynamicColors)
length(dynamicColors)
length(probe)
library(clusterProfiler)
library(org.Hs.eg.db)



######
module="green"
inModule = (dynamicColors==module);
modProbes = probe[inModule]
length(modProbes)
write.table(modProbes,"~/green.txt",append = F,quote = F,sep = "/t",row.names = F,col.names = F)
gene<- bitr(modProbes, fromType="ENSEMBL", toType = c("SYMBOL", "ENTREZID"), OrgDb=org.Hs.eg.db, drop = T)
ego_BP <- enrichGO(gene = gene$ENTREZID,
                   OrgDb= org.Hs.eg.db,
                   ont = "BP",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05,
                   readable = TRUE)
(ego_BP[,1:6])
write.csv(ego_BP,'~/muscle-green-all-human.csv',row.names = F)


######
module="blue"
inModule = (dynamicColors==module);
modProbes = probe[inModule]
length(modProbes)
write.table(modProbes,"~/blue.txt",append = F,quote = F,sep = "/t",row.names = F,col.names = F)
gene<- bitr(modProbes, fromType="ENSEMBL", toType = c("SYMBOL", "ENTREZID"), OrgDb=org.Hs.eg.db, drop = T)
ego_BP <- enrichGO(gene = gene$ENTREZID,
                   OrgDb= org.Hs.eg.db,
                   ont = "BP",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05,
                   readable = TRUE)
(ego_BP[,1:6])
write.csv(ego_BP,'~/muscle-blue-all-human.csv',row.names = F)

######
######
module="black"
inModule = (dynamicColors==module);
modProbes = probe[inModule]
length(modProbes)
write.table(modProbes,"~/black.txt",append = F,quote = F,sep = "/t",row.names = F,col.names = F)
gene<- bitr(modProbes, fromType="ENSEMBL", toType = c("SYMBOL", "ENTREZID"), OrgDb=org.Hs.eg.db, drop = T)
ego_BP <- enrichGO(gene = gene$ENTREZID,
                   OrgDb= org.Hs.eg.db,
                   ont = "BP",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05,
                   readable = TRUE)
(ego_BP[,1:6])
write.csv(ego_BP,'~/muscle-black-all-human.csv',row.names = F)

#######

module="brown"
inModule = (dynamicColors==module);
modProbes = probe[inModule]
length(modProbes)
write.table(modProbes,"~/brown.txt",append = F,quote = F,sep = "/t",row.names = F,col.names = F)
gene<- bitr(modProbes, fromType="ENSEMBL", toType = c("SYMBOL", "ENTREZID"), OrgDb=org.Hs.eg.db, drop = T)
ego_BP <- enrichGO(gene = gene$ENTREZID,
                   OrgDb= org.Hs.eg.db,
                   ont = "BP",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05,
                   readable = TRUE)
(ego_BP[,1:6])
write.csv(ego_BP,'~/muscle-brown-all-human.csv',row.names = F)
############3
#######

module="pink"
inModule = (dynamicColors==module);
modProbes = probe[inModule]
length(modProbes)
write.table(modProbes,"~/pink.txt",append = F,quote = F,sep = "/t",row.names = F,col.names = F)
gene<- bitr(modProbes, fromType="ENSEMBL", toType = c("SYMBOL", "ENTREZID"), OrgDb=org.Hs.eg.db, drop = T)
ego_BP <- enrichGO(gene = gene$ENTREZID,
                   OrgDb= org.Hs.eg.db,
                   ont = "BP",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05,
                   readable = TRUE)
(ego_BP[,1:6])
write.csv(ego_BP,'~/muscle-pink-all-human.csv',row.names = F)

#######
#######

module="red"
inModule = (dynamicColors==module);
modProbes = probe[inModule]
length(modProbes)
write.table(modProbes,"~/red.txt",append = F,quote = F,sep = "/t",row.names = F,col.names = F)
gene<- bitr(modProbes, fromType="ENSEMBL", toType = c("SYMBOL", "ENTREZID"), OrgDb=org.Hs.eg.db, drop = T)
ego_BP <- enrichGO(gene = gene$ENTREZID,
                   OrgDb= org.Hs.eg.db,
                   ont = "BP",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05,
                   readable = TRUE)
(ego_BP[,1:6])
write.csv(ego_BP,'~/muscle-red-all-human.csv',row.names = F)

######
module="turquoise"
inModule = (dynamicColors==module);
modProbes = probe[inModule]
length(modProbes)
write.table(modProbes,"~/turquoise.txt",append = F,quote = F,sep = "/t",row.names = F,col.names = F)
gene<- bitr(modProbes, fromType="ENSEMBL", toType = c("SYMBOL", "ENTREZID"), OrgDb=org.Hs.eg.db, drop = T)
ego_BP <- enrichGO(gene = gene$ENTREZID,
                   OrgDb= org.Hs.eg.db,
                   ont = "BP",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05,
                   readable = TRUE)
(ego_BP[,1:6])
write.csv(ego_BP,'~/muscle-turquoise-all-human.csv',row.names = F)
########
module="yellow"
inModule = (dynamicColors==module);
modProbes = probe[inModule]
length(modProbes)
write.table(modProbes,"~/yellow.txt",append = F,quote = F,sep = "/t",row.names = F,col.names = F)
gene<- bitr(modProbes, fromType="ENSEMBL", toType = c("SYMBOL", "ENTREZID"), OrgDb=org.Hs.eg.db, drop = T)
ego_BP <- enrichGO(gene = gene$ENTREZID,
                   OrgDb= org.Hs.eg.db,
                   ont = "BP",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05,
                   readable = TRUE)
(ego_BP[,1:6])
write.csv(ego_BP,'~/muscle-yellow-all-human.csv',row.names = F)


############


############
color_cor<-data.frame(correlation=c(0.52,0.54,0.41,0.03,0.22,0.21,0.48,0.14),col=c("brown","blue","turquoise","green","yellow","red","pink","black"))
color_cor<-as.data.frame(color_cor)

color_cor<-color_cor[order(color_cor$correlation),]
color_cor$col

tiff("~/correlation.tiff",res = 300, width = 1600, height = 1300,compression = "lzw")

par(mfrow=c(1,1),mar=c(5,6,4,4),oma = c(0,0,0,0) + 0.1)
barplot(color_cor$correlation,ylim = c(0,0.6),legend=F,las="2",cex.axis = 1.2,cex.lab=1.6,
        axisnames=T,col = c("green",    "black",     "red ","yellow",    "turquoise","pink","brown", "blue "),ylab="Correlation of gene connectivity")
title(xlab = "Module",line = 1.2,cex.lab=1.5)
dev.off()
color_cor$correlation[order(color_cor$correlation)]
color_cor$col[order(color_cor$correlation)]
#####
human_expression<-(tissue_exp_human[,rest1])
cattle_expression<-(tissue_exp_cattle[,rest1])

dim(human_expression)
probe<-colnames(human_expression)
table(dynamicColors)

module="black"
inModule = (dynamicColors==module);
modProbes = probe[inModule]
human_color<-human_expression[,modProbes]
cattle_color<-cattle_expression[,modProbes]
adjacency_cattle = abs(cor(cattle_color ,use = "p"))
diag(adjacency_cattle) = 0
adjacency_human = abs(cor(human_color ,use = "p"))
diag(adjacency_human) = 0
Connectivity_human <- apply(adjacency_human,1,sum,na.rm = TRUE)
Connectivity_cattle <- apply(adjacency_cattle,1,sum,na.rm = TRUE)
Connectivity_human = Connectivity_human/max(Connectivity_human)
Connectivity_cattle = Connectivity_cattle/max(Connectivity_cattle)
a<-cor.test(Connectivity_human,Connectivity_cattle,method = "s",use = "p")
a$estimate
tiff(file = "~/black.tiff",##reqiured to change
     res = 300, width = 1500, height = 1500,compression = "lzw")


plot(Connectivity_human,Connectivity_cattle,main = c("rho: ",a$estimate,"p: ",a$p.value), xlab = "Human network connectivity", ylab = "Cattle network connectivity", cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5)

abline(lm(Connectivity_cattle~Connectivity_human),col = "red",lwd = 2)
dev.off()


module="blue"
inModule = (dynamicColors==module);
modProbes = probe[inModule]
human_color<-human_expression[,modProbes]
cattle_color<-cattle_expression[,modProbes]
adjacency_cattle = abs(cor(cattle_color ,use = "p"))
diag(adjacency_cattle) = 0
adjacency_human = abs(cor(human_color ,use = "p"))
diag(adjacency_human) = 0
Connectivity_human <- apply(adjacency_human,1,sum,na.rm = TRUE)
Connectivity_cattle <- apply(adjacency_cattle,1,sum,na.rm = TRUE)
Connectivity_human = Connectivity_human/max(Connectivity_human)
Connectivity_cattle = Connectivity_cattle/max(Connectivity_cattle)
a<-cor.test(Connectivity_human,Connectivity_cattle,method = "s",use = "p")
a$estimate
tiff(file = "~/blue.tiff",##reqiured to change
     res = 300, width = 1500, height = 1500,compression = "lzw")


plot(Connectivity_human,Connectivity_cattle,main = c("rho: ",a$estimate,"p: ",a$p.value), xlab = "Human network connectivity", ylab = "Cattle network connectivity", cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5)

abline(lm(Connectivity_cattle~Connectivity_human),col = "red",lwd = 2)
dev.off()

module="brown"
inModule = (dynamicColors==module);
modProbes = probe[inModule]
human_color<-human_expression[,modProbes]
cattle_color<-cattle_expression[,modProbes]
adjacency_cattle = abs(cor(cattle_color ,use = "p"))
diag(adjacency_cattle) = 0
adjacency_human = abs(cor(human_color ,use = "p"))
diag(adjacency_human) = 0
Connectivity_human <- apply(adjacency_human,1,sum,na.rm = TRUE)
Connectivity_cattle <- apply(adjacency_cattle,1,sum,na.rm = TRUE)
Connectivity_human = Connectivity_human/max(Connectivity_human)
Connectivity_cattle = Connectivity_cattle/max(Connectivity_cattle)
a<-cor.test(Connectivity_human,Connectivity_cattle,method = "s",use = "p")
a$estimate
tiff(file = "~/brown.tiff",##reqiured to change
     res = 300, width = 1500, height = 1500,compression = "lzw")


plot(Connectivity_human,Connectivity_cattle,main = c("rho: ",a$estimate,"p: ",a$p.value), xlab = "Human network connectivity", ylab = "Cattle network connectivity", cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5)

abline(lm(Connectivity_cattle~Connectivity_human),col = "red",lwd = 2)
dev.off()


module="green"
inModule = (dynamicColors==module);
modProbes = probe[inModule]
human_color<-human_expression[,modProbes]
cattle_color<-cattle_expression[,modProbes]
adjacency_cattle = abs(cor(cattle_color ,use = "p"))
diag(adjacency_cattle) = 0
adjacency_human = abs(cor(human_color ,use = "p"))
diag(adjacency_human) = 0
Connectivity_human <- apply(adjacency_human,1,sum,na.rm = TRUE)
Connectivity_cattle <- apply(adjacency_cattle,1,sum,na.rm = TRUE)
Connectivity_human = Connectivity_human/max(Connectivity_human)
Connectivity_cattle = Connectivity_cattle/max(Connectivity_cattle)
a<-cor.test(Connectivity_human,Connectivity_cattle,method = "s",use = "p")
a$estimate
tiff(file = "~/green.tiff",##reqiured to change
     res = 300, width = 1500, height = 1500,compression = "lzw")


plot(Connectivity_human,Connectivity_cattle,main = c("rho: ",a$estimate,"p: ",a$p.value), xlab = "Human network connectivity", ylab = "Cattle network connectivity", cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5)

abline(lm(Connectivity_cattle~Connectivity_human),col = "red",lwd = 2)
dev.off()



module="grey"
inModule = (dynamicColors==module);
modProbes = probe[inModule]
human_color<-human_expression[,modProbes]
cattle_color<-cattle_expression[,modProbes]
adjacency_cattle = abs(cor(cattle_color ,use = "p"))
diag(adjacency_cattle) = 0
adjacency_human = abs(cor(human_color ,use = "p"))
diag(adjacency_human) = 0
Connectivity_human <- apply(adjacency_human,1,sum,na.rm = TRUE)
Connectivity_cattle <- apply(adjacency_cattle,1,sum,na.rm = TRUE)
Connectivity_human = Connectivity_human/max(Connectivity_human)
Connectivity_cattle = Connectivity_cattle/max(Connectivity_cattle)
a<-cor.test(Connectivity_human,Connectivity_cattle,method = "s",use = "p")
a$estimate
tiff(file = "~/grey.tiff",##reqiured to change
     res = 300, width = 1500, height = 1500,compression = "lzw")


plot(Connectivity_human,Connectivity_cattle,main = c("rho: ",a$estimate,"p: ",a$p.value), xlab = "Human network connectivity", ylab = "Cattle network connectivity", cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5)

abline(lm(Connectivity_cattle~Connectivity_human),col = "red",lwd = 2)
dev.off()


module="pink"
inModule = (dynamicColors==module);
modProbes = probe[inModule]
human_color<-human_expression[,modProbes]
cattle_color<-cattle_expression[,modProbes]
adjacency_cattle = abs(cor(cattle_color ,use = "p"))
diag(adjacency_cattle) = 0
adjacency_human = abs(cor(human_color ,use = "p"))
diag(adjacency_human) = 0
Connectivity_human <- apply(adjacency_human,1,sum,na.rm = TRUE)
Connectivity_cattle <- apply(adjacency_cattle,1,sum,na.rm = TRUE)
Connectivity_human = Connectivity_human/max(Connectivity_human)
Connectivity_cattle = Connectivity_cattle/max(Connectivity_cattle)
a<-cor.test(Connectivity_human,Connectivity_cattle,method = "s",use = "p")
a$estimate
tiff(file = "~/pink.tiff",##reqiured to change
     res = 300, width = 1500, height = 1500,compression = "lzw")


plot(Connectivity_human,Connectivity_cattle,main = c("rho: ",a$estimate,"p: ",a$p.value), xlab = "Human network connectivity", ylab = "Cattle network connectivity", cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5)

abline(lm(Connectivity_cattle~Connectivity_human),col = "red",lwd = 2)
dev.off()


module="red"
inModule = (dynamicColors==module);
modProbes = probe[inModule]
human_color<-human_expression[,modProbes]
cattle_color<-cattle_expression[,modProbes]
adjacency_cattle = abs(cor(cattle_color ,use = "p"))
diag(adjacency_cattle) = 0
adjacency_human = abs(cor(human_color ,use = "p"))
diag(adjacency_human) = 0
Connectivity_human <- apply(adjacency_human,1,sum,na.rm = TRUE)
Connectivity_cattle <- apply(adjacency_cattle,1,sum,na.rm = TRUE)
Connectivity_human = Connectivity_human/max(Connectivity_human)
Connectivity_cattle = Connectivity_cattle/max(Connectivity_cattle)
a<-cor.test(Connectivity_human,Connectivity_cattle,method = "s",use = "p")
a$estimate
tiff(file = "~/red.tiff",##reqiured to change
     res = 300, width = 1500, height = 1500,compression = "lzw")


plot(Connectivity_human,Connectivity_cattle,main = c("rho: ",a$estimate,"p: ",a$p.value), xlab = "Human network connectivity", ylab = "Cattle network connectivity", cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5)

abline(lm(Connectivity_cattle~Connectivity_human),col = "red",lwd = 2)
dev.off()


module="turquoise"
inModule = (dynamicColors==module);
modProbes = probe[inModule]
human_color<-human_expression[,modProbes]
cattle_color<-cattle_expression[,modProbes]
adjacency_cattle = abs(cor(cattle_color ,use = "p"))
diag(adjacency_cattle) = 0
adjacency_human = abs(cor(human_color ,use = "p"))
diag(adjacency_human) = 0
Connectivity_human <- apply(adjacency_human,1,sum,na.rm = TRUE)
Connectivity_cattle <- apply(adjacency_cattle,1,sum,na.rm = TRUE)
Connectivity_human = Connectivity_human/max(Connectivity_human)
Connectivity_cattle = Connectivity_cattle/max(Connectivity_cattle)
a<-cor.test(Connectivity_human,Connectivity_cattle,method = "s",use = "p")
a$estimate
tiff(file = "~/turquoise.tiff",##reqiured to change
     res = 300, width = 1500, height = 1500,compression = "lzw")


plot(Connectivity_human,Connectivity_cattle,main = c("rho: ",a$estimate,"p: ",a$p.value), xlab = "Human network connectivity", ylab = "Cattle network connectivity", cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5)

abline(lm(Connectivity_cattle~Connectivity_human),col = "red",lwd = 2)
dev.off()


module="yellow"
inModule = (dynamicColors==module);
modProbes = probe[inModule]
human_color<-human_expression[,modProbes]
cattle_color<-cattle_expression[,modProbes]
adjacency_cattle = abs(cor(cattle_color ,use = "p"))
diag(adjacency_cattle) = 0
adjacency_human = abs(cor(human_color ,use = "p"))
diag(adjacency_human) = 0
Connectivity_human <- apply(adjacency_human,1,sum,na.rm = TRUE)
Connectivity_cattle <- apply(adjacency_cattle,1,sum,na.rm = TRUE)
Connectivity_human = Connectivity_human/max(Connectivity_human)
Connectivity_cattle = Connectivity_cattle/max(Connectivity_cattle)
a<-cor.test(Connectivity_human,Connectivity_cattle,method = "s",use = "p")
a$estimate
tiff(file = "~/yellow.tiff",##reqiured to change
     res = 300, width = 1500, height = 1500,compression = "lzw")


plot(Connectivity_human,Connectivity_cattle,main = c("rho: ",a$estimate,"p: ",a$p.value), xlab = "Human network connectivity", ylab = "Cattle network connectivity", cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5)

abline(lm(Connectivity_cattle~Connectivity_human),col = "red",lwd = 2)
dev.off()

table(dynamicColors)
module="pink"
inModule = (dynamicColors==module);
modProbes = probe[inModule]
human_color<-human_expression[,modProbes]
cattle_color<-cattle_expression[,modProbes]
adjacency_cattle = abs(cor(cattle_color ,use = "p"))
diag(adjacency_cattle) = 0
adjacency_human = abs(cor(human_color ,use = "p"))
diag(adjacency_human) = 0
Connectivity_human <- apply(adjacency_human,1,sum,na.rm = TRUE)
Connectivity_cattle <- apply(adjacency_cattle,1,sum,na.rm = TRUE)
Connectivity_human = Connectivity_human/max(Connectivity_human)
Connectivity_cattle = Connectivity_cattle/max(Connectivity_cattle)
a<-cor.test(Connectivity_human,Connectivity_cattle,method = "s",use = "p")
a$estimate
tiff(file = "~/pink.tiff",##reqiured to change
     res = 300, width = 1500, height = 1500,compression = "lzw")


plot(Connectivity_human,Connectivity_cattle,main = c("rho: ",a$estimate,"p: ",a$p.value), xlab = "Human network connectivity", ylab = "Cattle network connectivity", cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5)

abline(lm(Connectivity_cattle~Connectivity_human),col = "red",lwd = 2)
dev.off()


