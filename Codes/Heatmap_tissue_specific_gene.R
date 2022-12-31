#BiocManager::install("ComplexHeatmap")
rm(list = ls())
library(pheatmap)
####################################

load("orthologous_cattle_human.Rdata")
load("human10830_cattle4866_17315gene.Rdata")
load("15696samples_17315genes.Rdata")
load("overlap_geneid.Rdata")
load("top10_gene.Rdata")
load("up_gene_cattle.Rdata")
load("overlap.Rdata")
load("up_gene_human.Rdata")





######cattle in cattle
length(unique(top10_gene_cattle$gene))
cattle_gene<-top10_gene_cattle[!duplicated(top10_gene_cattle$gene),]
colnames(Cattle_expression)[1:10]

TPM_gene<-Cattle_expression[,match(cattle_gene$gene,colnames(Cattle_expression))]
dim(TPM_gene)
TPM_gene

table(colnames(TPM_gene)==cattle_gene$gene)
sub_info<-Meta_data_cattle[order(Meta_data_cattle$tissue_new),]
table(sub_info$Sample==rownames(TPM_gene))
TPM_gene<-TPM_gene[match(sub_info$Sample,rownames(TPM_gene)),]
table(sub_info$Sample==rownames(TPM_gene))
TPM_gene<-t(log2(TPM_gene+0.25))
dim(TPM_gene)
TPM_gene[1:10,1:10]




annotation_row=data.frame(cattle_gene$tissue)
rownames(annotation_row)<-rownames(TPM_gene)

annotation_row




Tissue<-sub_info$tissue_new
annotation_mix = data.frame(Tissue)
colnames(annotation_row)<-"Tissue-specific gene"
col<-list(Tissue=c("#00008B","#008B8B","#8B008B","#8B0000","#828282","#00B2EE","#0000FF","#FF0000","#FF00FF","#FF8C00","#FFFF00","#006400","#00FF00",
                   "#CD950C","#FF6A6A","#CCCC33","pink","#FF1493","#8B7355","#00F5FF"),'Tissue-specific gene' =c("#00008B","#008B8B","#8B008B","#8B0000","#828282","#00B2EE","#0000FF","#FF0000","#FF00FF","#FF8C00","#FFFF00","#006400","#00FF00",
                                                                                                         "#CD950C","#FF6A6A","#CCCC33","pink","#FF1493","#8B7355","#00F5FF"))

names(col$Tissue) <- unique(sub_info$tissue_new)
names(col$'Tissue-specific gene')<-unique(cattle_gene$tissue)
rownames(annotation_mix)=colnames(TPM_gene)
length(rownames(annotation_mix))
length(colnames(TPM_gene))

breaksList = seq(0, 9, by = 1)
colnames(annotation_row)<-"Tissue-specific gene"

col
annotation_row
pheatmap(TPM_gene,scale = "row",cluster_rows  = F,
         annotation_col=annotation_mix,annotation_colors = col,annotation_row = annotation_row,
           cluster_cols  = F, show_rownames=F,show_colnames = F,color = colorRampPalette(c("white","firebrick3"))(9),
         main = '  ',legend = T,breaks = breaksList,legend_breaks = seq(0, 8, by = 2),annotation_names_row = F,width = 10,height = 10,
         file="top10cattle_incattle.tiff")



##########cattle_in human

human_gene<-Orthologous_human_cattle$`Gene stable ID`[match(cattle_gene$gene,Orthologous_human_cattle$`Cow gene stable ID`)]

TPM_gene<-human_expression[,match(human_gene,colnames(human_expression))]
dim(TPM_gene)

table(colnames(TPM_gene)==human_gene)
sub_info<-Meta_data_human[order(Meta_data_human$tissue_new),]
table(sub_info$SAMPID==rownames(TPM_gene))
TPM_gene<-TPM_gene[match(sub_info$SAMPID,rownames(TPM_gene)),]
table(sub_info$SAMPID==rownames(TPM_gene))
TPM_gene<-t(log2(TPM_gene+0.25))
dim(TPM_gene)



annotation_row=data.frame(cattle_gene$tissue)
rownames(annotation_row)<-rownames(TPM_gene)
colnames(annotation_row)<-"Tissue-specific gene"
annotation_row

Tissue<-sub_info$tissue_new
annotation_mix = data.frame(Tissue)


col<-list(Tissue=c("#00008B","#008B8B","#8B008B","#8B0000","#828282","#00B2EE","#0000FF","#FF0000","#FF00FF","#FF8C00","#FFFF00","#006400","#00FF00",
                   "#CD950C","#FF6A6A","#CCCC33","pink","#FF1493","#8B7355","#00F5FF"),'Tissue-specific gene' =c("#00008B","#008B8B","#8B008B","#8B0000","#828282","#00B2EE","#0000FF","#FF0000","#FF00FF","#FF8C00","#FFFF00","#006400","#00FF00",
                                                                                                                 "#CD950C","#FF6A6A","#CCCC33","pink","#FF1493","#8B7355","#00F5FF"))

names(col$Tissue) <- unique(sub_info$tissue_new)
names(col$'Tissue-specific gene')<-unique(cattle_gene$tissue)
rownames(annotation_mix)=colnames(TPM_gene)

col


breaksList = seq(0, 9, by = 1)
pheatmap(TPM_gene,scale = "row",cluster_rows  = F,
         annotation_col=annotation_mix,annotation_colors = col,annotation_row = annotation_row,
         cluster_cols  = F, show_rownames=F,show_colnames = F,color = colorRampPalette(c("white","firebrick3"))(9),
         main = '  ',legend = T,breaks = breaksList,legend_breaks = seq(0, 8, by = 2),annotation_names_row = F,width = 10,height = 10,
         file="top10_cattle_in_human.tiff")



######human_in cattle
length(unique(top10_gene_human$gene))
human_gene<-unique(top10_gene_human$gene)
TPM_gene<-human_expression[,colnames(human_expression)%in%unique(human_gene)]
dim(TPM_gene)
TPM_gene<-TPM_gene[,order(factor(colnames(TPM_gene),levels=unique(human_gene)))]
table(colnames(TPM_gene)==human_gene)
sub_info<-Meta_data_human[order(Meta_data_human$tissue_new),]
table(sub_info$SAMPID==rownames(TPM_gene))
TPM_gene<-TPM_gene[match(sub_info$SAMPID,rownames(TPM_gene)),]
table(sub_info$SAMPID==rownames(TPM_gene))
TPM_gene<-t(log2(TPM_gene+0.25))
dim(TPM_gene)


Tissue<-sub_info$tissue_new
annotation_mix = data.frame(Tissue)
col<-list(Tissue=c("#00008B","#008B8B","#8B008B","#8B0000","#828282","#00B2EE","#0000FF","#FF0000","#FF00FF","#FF8C00","#FFFF00","#006400","#00FF00",
                   "#CD950C","#FF6A6A","#CCCC33","pink","#FF1493","#8B7355","#00F5FF"))
names(col$Tissue) <- unique(sub_info$tissue_new)
rownames(annotation_mix)=colnames(TPM_gene)
length(rownames(annotation_mix))
length(colnames(TPM_gene))
breaksList = seq(0, 10, by = 1)
pheatmap(TPM_gene,scale = "row",cluster_rows  = F,
         annotation_col=annotation_mix,annotation_colors = col,
         cluster_cols  = F, show_rownames=F,show_colnames = F,color = colorRampPalette(c("white","firebrick3"))(10),
         main = '  ',legend = T,breaks = breaksList,
         file="top_10human_in_human.tiff")



##########human in cattle
human_gene<-unique(top10_gene_human$gene)
cattle_gene<-Orthologous_human_cattle$`Cow gene stable ID`[match(human_gene,Orthologous_human_cattle$`Gene stable ID`)]
length(cattle_gene)
TPM_gene<-Cattle_expression[,colnames(Cattle_expression)%in%unique(cattle_gene)]
dim(TPM_gene)
TPM_gene<-TPM_gene[,order(factor(colnames(TPM_gene),levels=unique(cattle_gene)))]
table(colnames(TPM_gene)==cattle_gene)
sub_info<-Meta_data_cattle[order(Meta_data_cattle$tissue_new),]
table(sub_info$Sample==rownames(TPM_gene))
TPM_gene<-TPM_gene[match(sub_info$Sample,rownames(TPM_gene)),]
table(sub_info$Sample==rownames(TPM_gene))
TPM_gene<-t(log2(TPM_gene+0.25))
dim(TPM_gene)

Tissue<-sub_info$tissue_new
annotation_mix = data.frame(Tissue)
col<-list(Tissue=c("#00008B","#008B8B","#8B008B","#8B0000","#828282","#00B2EE","#0000FF","#FF0000","#FF00FF","#FF8C00","#FFFF00","#006400","#00FF00",
                   "#CD950C","#FF6A6A","#CCCC33","pink","#FF1493","#8B7355","#00F5FF"))
names(col$Tissue) <- unique(sub_info$tissue_new)
rownames(annotation_mix)=colnames(TPM_gene)
length(rownames(annotation_mix))
length(colnames(TPM_gene))
breaksList = seq(0, 10, by = 1)
pheatmap(TPM_gene,scale = "row",cluster_rows  = F,
         annotation_col=annotation_mix,annotation_colors = col,
         cluster_cols  = F, show_rownames=F,show_colnames = F,color = colorRampPalette(c("white","firebrick3"))(10),
         main = '  ',legend = T,breaks = breaksList,
         file="top10_human_in_cattle.tiff")







