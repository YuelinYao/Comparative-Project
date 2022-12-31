rm(list = ls())
load("human10830_cattle4866_17315gene.Rdata")
load("orthologous_cattle_human.Rdata")
library(ggplot2)

######cattle_specific
expression_data<-data.frame(expression=NA,tissue=NA,Species=NA,gene=NA)


gene_list<-c("CSN2","CSN3","CCL27")

cattle_specific<-Orthologous_human_cattle$`Cow gene stable ID`[Orthologous_human_cattle$`Cow gene name`%in%gene_list]



i=1
for (i in 1:length(cattle_specific)){
  gene<-cattle_specific[i]
  gene_names<-Orthologous_human_cattle$`Cow gene name`[Orthologous_human_cattle$`Cow gene stable ID`%in%gene]
  gene_names
  gene_expression<-Cattle_expression[,colnames(Cattle_expression)%in%gene]
  data<-data.frame(expression=gene_expression,tissue=Meta_data_cattle$tissue_new,Species="Cattle",gene=gene_names)
  
  
  gene_tran<-Orthologous_human_cattle$`Gene stable ID`[Orthologous_human_cattle$`Cow gene stable ID`%in%gene]
  gene_expression_human<-human_expression[,colnames(human_expression)%in%gene_tran]
  data_human<-data.frame(expression=gene_expression_human,tissue=Meta_data_human$tissue_new,Species="Human",gene=gene_names)
  expression_human_cattle<-rbind(data,data_human)
  
  expression_data<-rbind(expression_data,expression_human_cattle)
  
}
table(expression_data$gene)
expression_data$expression<-log2(expression_data$expression+0.25)
expression_data<-expression_data[-1,]
expression_data$gene<-factor(expression_data$gene,levels = c("CSN2","CSN3","CCL27"))
tiff(file="CSN2_genes.tiff",width = 2500,height = 1600,res = 300)
p<-ggplot(expression_data, aes(x=tissue, y=expression, fill=Species))+geom_boxplot(position=position_dodge(1),outlier.shape = NA)+scale_fill_manual(values=c("blue", "red"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size =18 ,colour = "black"),panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour="black"),axis.title =element_text(size =25,colour = "black" ),axis.text.y = element_text(size =15,colour = "black" ),legend.text = element_text(size =15,colour = "black" ) ,legend.title =element_text(size =20,colour = "black" ),panel.background = element_rect(fill='transparent', color=NA) )+ylab(expression(paste(log[2],italic(TPM+0.25))))+xlab(" ")
p+facet_wrap( ~gene,scales = "free_y",ncol=1)+theme(
  strip.background = element_rect(color = NA,fill=NA,
                                  size=1),strip.text = element_text(face="bold.italic",size =18)
)


dev.off()
