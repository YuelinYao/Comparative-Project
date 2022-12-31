library(xlsx)

GO_cattle_mammary <- read.xlsx("~/species_specific/cattle_go/cattle_Mammary.xlsx",header = T,stringsAsFactors = F,sheetIndex = 1)
head(GO_cattle_mammary,20)

GO_cattle_mammary <- GO_cattle_mammary[c(1,5),c(3,4,7)]
GO_cattle_mammary$logFDR <- -log(GO_cattle_mammary$p.adjust,10)
GO_cattle_mammary$Species<-"Cattle"



GO_human_mammary <- read.xlsx("~/species_specific/human_go/human_Mammary.xlsx",header = T,stringsAsFactors = F,sheetIndex = 1)

head(GO_human_mammary,20)

GO_human_mammary <- GO_human_mammary[c(1,2),c(3,4,7)]
GO_human_mammary$logFDR <- -log(GO_human_mammary$p.adjust,10)
GO_human_mammary$Species<-"Human"

#####3


GO_human_testis <- read.xlsx("~/species_specific/human_go/human_Testis.xlsx",header = T,stringsAsFactors = F,sheetIndex = 1)

head(GO_human_testis,20)
GO_human_testis <- GO_human_testis[c(2,4),c(3,4,7)]
GO_human_testis$logFDR <- -log(GO_human_testis$p.adjust,10)
GO_human_testis$Species<-"Human"


GO_cattle_testis <- read.xlsx("~/species_specific/cattle_go/cattle_Testis.xlsx",header = T,stringsAsFactors = F,sheetIndex = 1)
head(GO_cattle_testis,30)
GO_cattle_testis <- GO_cattle_testis[c(1,2),c(3,4,7)]
GO_cattle_testis$logFDR <- -log(GO_cattle_testis$p.adjust,10)
GO_cattle_testis$Species<-"Cattle"





########
GO_human_Heart <- read.xlsx("~/species_specific/human_go/human_Heart.xlsx",header = T,stringsAsFactors = F,sheetIndex = 1)

head(GO_human_Heart,20)
GO_human_Heart <- GO_human_Heart[c(1,2),c(3,4,7)]
GO_human_Heart$logFDR <- -log(GO_human_Heart$p.adjust,10)
GO_human_Heart$Species<-"Human"


GO_cattle_Heart <- read.xlsx("~/species_specific/cattle_go/cattle_Heart.xlsx",header = T,stringsAsFactors = F,sheetIndex = 1)
head(GO_cattle_Heart,30)
GO_cattle_Heart <- GO_cattle_Heart[c(2,4),c(3,4,7)]
GO_cattle_Heart$logFDR <- -log(GO_cattle_Heart$p.adjust,10)
GO_cattle_Heart$Species<-"Cattle"



############3
GO_human_Brain <- read.xlsx("~/species_specific/human_go/human_Brain.xlsx",header = T,stringsAsFactors = F,sheetIndex = 1)

head(GO_human_Brain,20)
GO_human_Brain <- GO_human_Brain[c(4,9),c(3,4,7)]
GO_human_Brain$logFDR <- -log(GO_human_Brain$p.adjust,10)
GO_human_Brain$Species<-"Human"


GO_cattle_Brain <- read.xlsx("~/species_specific/cattle_go/cattle_Brain.xlsx",header = T,stringsAsFactors = F,sheetIndex = 1)
head(GO_cattle_Brain,30)
GO_cattle_Brain <- GO_cattle_Brain[c(1,2),c(3,4,7)]
GO_cattle_Brain$logFDR <- -log(GO_cattle_Brain$p.adjust,10)
GO_cattle_Brain$Species<-"Cattle"



##########3
GO_human_Ovary <- read.xlsx("~/species_specific/human_go/human_Ovary.xlsx",header = T,stringsAsFactors = F,sheetIndex = 1)

head(GO_human_Ovary,20)
GO_human_Ovary <- GO_human_Ovary[c(4,5),c(3,4,7)]
GO_human_Ovary$logFDR <- -log(GO_human_Ovary$p.adjust,10)
GO_human_Ovary$Species<-"Human"


GO_cattle_Ovary <- read.xlsx("~/species_specific/cattle_go/cattle_Ovary.xlsx",header = T,stringsAsFactors = F,sheetIndex = 1)
head(GO_cattle_Ovary,30)
GO_cattle_Ovary <- GO_cattle_Ovary[c(1,10),c(3,4,7)]
GO_cattle_Ovary$logFDR <- -log(GO_cattle_Ovary$p.adjust,10)
GO_cattle_Ovary$Species<-"Cattle"






#######3
GO_rep <- rbind(GO_cattle_testis,GO_human_testis,GO_cattle_mammary,GO_human_mammary,GO_cattle_Brain,GO_human_Brain,GO_cattle_Heart,GO_human_Heart,GO_cattle_Ovary,GO_human_Ovary)
GO_rep$Category <- factor(rep(c("Testis","Mammary","Brain","Heart","Ovary"),each=4),levels = c("Testis","Mammary","Brain","Heart","Ovary"))
GO_rep$"Description" <- factor(GO_rep$"Description",levels = rev(GO_rep$"Description"))
GO_rep$Description
GO_rep

colnames(GO_rep)<-c("Description","GeneRatio","FDR","logFDR","Species","Category")

GO_rep


library(xlsx)
library(ggplot2)
GO_rep$FDR

GO_rep
GO_rep




GO_rep$number <- factor(rev(1:nrow(GO_rep)))
## shorten the names of GO terms


tiff(file = "~/species_ENRICHMENT.tiff",res = 300, width = 3300, height = 1200,compression = "lzw")

p<-ggplot(data=GO_rep, aes(x=factor(Description,level=Description), y=logFDR,fill=Species),split='Species') +scale_fill_manual(values=c("blue", "red"))+geom_bar(stat="identity")+xlab("")+ylab(expression(paste(-log[10],italic(FDR))))+theme(axis.text.x = element_text(size=10,color = "black"),axis.text.y = element_text(size=13,colour = "black"),legend.text =element_text(size=15),legend.title = element_text(size=15),strip.text =element_text(size=15),axis.title.x =element_text(size=15) )+facet_wrap( ~Category,scales="free_x",nrow = 1)

p + coord_flip()+theme(panel.grid.major=element_line(colour=NA),panel.background = element_rect(fill='transparent', color='black'))
dev.off()
