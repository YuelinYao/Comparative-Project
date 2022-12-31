rm(list = ls()) 
library(randomcoloR)
library(ggplot2)
library(Rtsne)
library("corrplot")
library(RColorBrewer)
#load("tissue_specifiy.Rdata")
#save(tissue_pc_merge20_human,tissue_pc_merge20,file="tissue_specifiy.Rdata")
#cattle_mapping <- read.table("data_info_tissue_breed_with_mapping_rate.txt",header = T,sep="\t",stringsAsFactors = F) 
##################

load("human10830_cattle4866_17315gene.Rdata")
load("15696samples_17315genes.Rdata")
load("orthologous_cattle_human.Rdata")


####category function
category <- function(x) {
  table(cut(x,c(0,0.1,1,10,100,1000,300000)))
}


#########cattle##########3
tissue_pc_mean<-aggregate(Cattle_expression,list(as.factor(Meta_data_cattle$tissue_new)),mean)
rownames(tissue_pc_mean)<-tissue_pc_mean[,1]
tissue_pc_mean<-tissue_pc_mean[,-1]

tissue_pc_mean[1:10,1:10]
tissue_pc_mean_tf <- tissue_pc_mean >=0.1
tissue_pc_mean_tf <- as.data.frame(tissue_pc_mean_tf)
tissue_pc_mean_tf[1:10,1:10]
tissue_pc_mean_count <- colSums(tissue_pc_mean_tf)
tissue_pc_mean_count <- as.data.frame(tissue_pc_mean_count)
tissue_pc_mean_count
tissue_pc_mean_table <- table(tissue_pc_mean_count)

tissue_pc_mean_table <- as.data.frame(tissue_pc_mean_table)
dim(tissue_pc_mean_table)
tissue_pc_mean_table[1:21,]

tissue_pc_mean1 <- as.data.frame(tissue_pc_mean)
tissue_pc_mean1[tissue_pc_mean1 <0.1] = 0
tissue_pc_mean2 <- apply(tissue_pc_mean1, 2, sum)

tissue_pc_mean3 <- cbind(tissue_pc_mean_count, tissue_pc_mean2)
tissue_pc_mean3$TPM <- tissue_pc_mean3$tissue_pc_mean2/tissue_pc_mean3$tissue_pc_mean_count

tissue_pc_mean3$Category[tissue_pc_mean3$TPM == "NaN"] <- "<0.1"
tissue_pc_mean3$Category[tissue_pc_mean3$TPM >=0.1 & tissue_pc_mean3$TPM <1] <- ">0.1 and <1"
tissue_pc_mean3$Category[tissue_pc_mean3$TPM >=1 & tissue_pc_mean3$TPM <10] <- ">1 and <10"
tissue_pc_mean3$Category[tissue_pc_mean3$TPM >=10 & tissue_pc_mean3$TPM <100] <-  ">10 and <100"
tissue_pc_mean3$Category[tissue_pc_mean3$TPM >=100 & tissue_pc_mean3$TPM <1000] <-  ">100 and <1000"
tissue_pc_mean3$Category[tissue_pc_mean3$TPM >=1000] <- ">1000"



tissue_pc_cate0 <- tissue_pc_mean3[tissue_pc_mean3$tissue_pc_mean_count==0,]
tissue_pc_cate0_sum <- table(tissue_pc_cate0$Category)
tissue_pc_cate0_sum <- as.data.frame(tissue_pc_cate0_sum)
colnames(tissue_pc_cate0_sum) <- c("Category","0")
tissue_pc_cate1 <- tissue_pc_mean3[tissue_pc_mean3$tissue_pc_mean_count==1,]
tissue_pc_cate1_sum <- table(tissue_pc_cate1$Category)
tissue_pc_cate1_sum <- as.data.frame(tissue_pc_cate1_sum)
colnames(tissue_pc_cate1_sum) <- c("Category","1")
tissue_pc_cate2 <- tissue_pc_mean3[tissue_pc_mean3$tissue_pc_mean_count==2,]
tissue_pc_cate2_sum <- table(tissue_pc_cate2$Category)
tissue_pc_cate2_sum <- as.data.frame(tissue_pc_cate2_sum)
colnames(tissue_pc_cate2_sum) <- c("Category","2")
tissue_pc_cate3 <- tissue_pc_mean3[tissue_pc_mean3$tissue_pc_mean_count==3,]

tissue_pc_cate3_sum <- table(tissue_pc_cate3$Category)
tissue_pc_cate3_sum <- as.data.frame(tissue_pc_cate3_sum)
colnames(tissue_pc_cate3_sum) <- c("Category","3")
tissue_pc_cate4 <- tissue_pc_mean3[tissue_pc_mean3$tissue_pc_mean_count==4,]
tissue_pc_cate4_sum <- table(tissue_pc_cate4$Category)
tissue_pc_cate4_sum <- as.data.frame(tissue_pc_cate4_sum)
colnames(tissue_pc_cate4_sum) <- c("Category","4")
tissue_pc_cate5 <- tissue_pc_mean3[tissue_pc_mean3$tissue_pc_mean_count==5,]
tissue_pc_cate5_sum <- table(tissue_pc_cate5$Category)
tissue_pc_cate5_sum <- as.data.frame(tissue_pc_cate5_sum)
colnames(tissue_pc_cate5_sum) <- c("Category","5")


tissue_pc_cate6 <- tissue_pc_mean3[tissue_pc_mean3$tissue_pc_mean_count==6,]
tissue_pc_cate6_sum <- table(tissue_pc_cate6$Category)
tissue_pc_cate6_sum <- as.data.frame(tissue_pc_cate6_sum)
colnames(tissue_pc_cate6_sum) <- c("Category","6")
tissue_pc_cate7 <- tissue_pc_mean3[tissue_pc_mean3$tissue_pc_mean_count==7,]
tissue_pc_cate7_sum <- table(tissue_pc_cate7$Category)
tissue_pc_cate7_sum <- as.data.frame(tissue_pc_cate7_sum)
colnames(tissue_pc_cate7_sum) <- c("Category","7")
tissue_pc_cate8 <- tissue_pc_mean3[tissue_pc_mean3$tissue_pc_mean_count==8,]
tissue_pc_cate8_sum <- table(tissue_pc_cate8$Category)
tissue_pc_cate8_sum <- as.data.frame(tissue_pc_cate8_sum)
colnames(tissue_pc_cate8_sum) <- c("Category","8")
tissue_pc_cate9 <- tissue_pc_mean3[tissue_pc_mean3$tissue_pc_mean_count==9,]
tissue_pc_cate9_sum <- table(tissue_pc_cate9$Category)
tissue_pc_cate9_sum <- as.data.frame(tissue_pc_cate9_sum)
colnames(tissue_pc_cate9_sum) <- c("Category","9")
tissue_pc_cate10 <- tissue_pc_mean3[tissue_pc_mean3$tissue_pc_mean_count==10,]
tissue_pc_cate10_sum <- table(tissue_pc_cate10$Category)
tissue_pc_cate10_sum <- as.data.frame(tissue_pc_cate10_sum)
colnames(tissue_pc_cate10_sum) <- c("Category","10")
tissue_pc_cate11 <- tissue_pc_mean3[tissue_pc_mean3$tissue_pc_mean_count==11,]
tissue_pc_cate11_sum <- table(tissue_pc_cate11$Category)
tissue_pc_cate11_sum <- as.data.frame(tissue_pc_cate11_sum)
colnames(tissue_pc_cate11_sum) <- c("Category","11")
tissue_pc_cate12 <- tissue_pc_mean3[tissue_pc_mean3$tissue_pc_mean_count==12,]
tissue_pc_cate12_sum <- table(tissue_pc_cate12$Category)
tissue_pc_cate12_sum <- as.data.frame(tissue_pc_cate12_sum)
colnames(tissue_pc_cate12_sum) <- c("Category","12")
tissue_pc_cate13 <- tissue_pc_mean3[tissue_pc_mean3$tissue_pc_mean_count==13,]
tissue_pc_cate13_sum <- table(tissue_pc_cate13$Category)
tissue_pc_cate13_sum <- as.data.frame(tissue_pc_cate13_sum)
colnames(tissue_pc_cate13_sum) <- c("Category","13")
tissue_pc_cate14 <- tissue_pc_mean3[tissue_pc_mean3$tissue_pc_mean_count==14,]
tissue_pc_cate14_sum <- table(tissue_pc_cate14$Category)
tissue_pc_cate14_sum <- as.data.frame(tissue_pc_cate14_sum)
colnames(tissue_pc_cate14_sum) <- c("Category","14")
tissue_pc_cate15 <- tissue_pc_mean3[tissue_pc_mean3$tissue_pc_mean_count==15,]
tissue_pc_cate15_sum <- table(tissue_pc_cate15$Category)
tissue_pc_cate15_sum <- as.data.frame(tissue_pc_cate15_sum)
colnames(tissue_pc_cate15_sum) <- c("Category","15")
tissue_pc_cate16 <- tissue_pc_mean3[tissue_pc_mean3$tissue_pc_mean_count==16,]
tissue_pc_cate16_sum <- table(tissue_pc_cate16$Category)
tissue_pc_cate16_sum <- as.data.frame(tissue_pc_cate16_sum)
colnames(tissue_pc_cate16_sum) <- c("Category","16")
tissue_pc_cate17 <- tissue_pc_mean3[tissue_pc_mean3$tissue_pc_mean_count==17,]
tissue_pc_cate17_sum <- table(tissue_pc_cate17$Category)
tissue_pc_cate17_sum <- as.data.frame(tissue_pc_cate17_sum)
colnames(tissue_pc_cate17_sum) <- c("Category","17")
tissue_pc_cate18 <- tissue_pc_mean3[tissue_pc_mean3$tissue_pc_mean_count==18,]
tissue_pc_cate18_sum <- table(tissue_pc_cate18$Category)
tissue_pc_cate18_sum <- as.data.frame(tissue_pc_cate18_sum)
colnames(tissue_pc_cate18_sum) <- c("Category","18")
tissue_pc_cate19 <- tissue_pc_mean3[tissue_pc_mean3$tissue_pc_mean_count==19,]
tissue_pc_cate19_sum <- table(tissue_pc_cate19$Category)
tissue_pc_cate19_sum <- as.data.frame(tissue_pc_cate19_sum)
colnames(tissue_pc_cate19_sum) <- c("Category","19")
tissue_pc_cate20 <- tissue_pc_mean3[tissue_pc_mean3$tissue_pc_mean_count==20,]
tissue_pc_cate20_sum <- table(tissue_pc_cate20$Category)
tissue_pc_cate20_sum <- as.data.frame(tissue_pc_cate20_sum)
colnames(tissue_pc_cate20_sum) <- c("Category","20")



tissue_pc_merge1 <- merge(tissue_pc_cate0_sum, tissue_pc_cate1_sum, all = T)
tissue_pc_merge2 <- merge(tissue_pc_merge1, tissue_pc_cate2_sum, all = T)
tissue_pc_merge3 <- merge(tissue_pc_merge2, tissue_pc_cate3_sum, all = T)
tissue_pc_merge4 <- merge(tissue_pc_merge3, tissue_pc_cate4_sum, all = T)
tissue_pc_merge5 <- merge(tissue_pc_merge4, tissue_pc_cate5_sum, all = T)
tissue_pc_merge6 <- merge(tissue_pc_merge5, tissue_pc_cate6_sum, all = T)
tissue_pc_merge7 <- merge(tissue_pc_merge6, tissue_pc_cate7_sum, all = T)
tissue_pc_merge8 <- merge(tissue_pc_merge7, tissue_pc_cate8_sum, all = T)
tissue_pc_merge9 <- merge(tissue_pc_merge8, tissue_pc_cate9_sum, all = T)
tissue_pc_merge10 <- merge(tissue_pc_merge9, tissue_pc_cate10_sum, all = T)
tissue_pc_merge11 <- merge(tissue_pc_merge10, tissue_pc_cate11_sum, all = T)
tissue_pc_merge12 <- merge(tissue_pc_merge11, tissue_pc_cate12_sum, all = T)
tissue_pc_merge13 <- merge(tissue_pc_merge12, tissue_pc_cate13_sum, all = T)
tissue_pc_merge14 <- merge(tissue_pc_merge13, tissue_pc_cate14_sum, all = T)
tissue_pc_merge15 <- merge(tissue_pc_merge14, tissue_pc_cate15_sum, all = T)
tissue_pc_merge16 <- merge(tissue_pc_merge15, tissue_pc_cate16_sum, all = T)
tissue_pc_merge17 <- merge(tissue_pc_merge16, tissue_pc_cate17_sum, all = T)
tissue_pc_merge18 <- merge(tissue_pc_merge17, tissue_pc_cate18_sum, all = T)
tissue_pc_merge19 <- merge(tissue_pc_merge18, tissue_pc_cate19_sum, all = T)
tissue_pc_merge20 <- merge(tissue_pc_merge19, tissue_pc_cate20_sum, all = T)

tissue_pc_merge20[is.na(tissue_pc_merge20)]=0
rownames(tissue_pc_merge20) <- tissue_pc_merge20$Category
tissue_pc_merge20 <- tissue_pc_merge20[,-1]
tissue_pc_merge20
tissue_pc_mean_table


#########human

tissue_pc_mean_human<-aggregate(human_expression,list(as.factor(Meta_data_human$tissue_new)),mean)
rownames(tissue_pc_mean_human)<-tissue_pc_mean_human[,1]
tissue_pc_mean_human<-tissue_pc_mean_human[,-1]

tissue_pc_mean_human[1:10,1:10]
tissue_pc_mean_tf_human <- tissue_pc_mean_human >=0.1
tissue_pc_mean_tf_human <- as.data.frame(tissue_pc_mean_tf_human)
tissue_pc_mean_count_human <- colSums(tissue_pc_mean_tf_human)
tissue_pc_mean_count_human <- as.data.frame(tissue_pc_mean_count_human)
tissue_pc_mean_count_human
tissue_pc_mean_table_human <- table(tissue_pc_mean_count_human)
tissue_pc_mean_table_human <- as.data.frame(tissue_pc_mean_table_human)


tissue_pc_mean1_human <- as.data.frame(tissue_pc_mean_human)
tissue_pc_mean1_human[tissue_pc_mean1_human <0.1] = 0
tissue_pc_mean2_human <- apply(tissue_pc_mean1_human, 2, sum)

tissue_pc_mean3_human <- cbind(tissue_pc_mean_count_human, tissue_pc_mean2_human)
tissue_pc_mean3_human$TPM <- tissue_pc_mean3_human$tissue_pc_mean2_human/tissue_pc_mean3_human$tissue_pc_mean_count_human

tissue_pc_mean3_human$Category[tissue_pc_mean3_human$TPM == "NaN"] <- "<0.1"
tissue_pc_mean3_human$Category[tissue_pc_mean3_human$TPM >=0.1 & tissue_pc_mean3_human$TPM <1] <- ">0.1 and <1"
tissue_pc_mean3_human$Category[tissue_pc_mean3_human$TPM >=1 & tissue_pc_mean3_human$TPM <10] <- ">1 and <10"
tissue_pc_mean3_human$Category[tissue_pc_mean3_human$TPM >=10 & tissue_pc_mean3_human$TPM <100] <-  ">10 and <100"
tissue_pc_mean3_human$Category[tissue_pc_mean3_human$TPM >=100 & tissue_pc_mean3_human$TPM <1000] <-  ">100 and <1000"
tissue_pc_mean3_human$Category[tissue_pc_mean3_human$TPM >=1000] <- ">1000"



tissue_pc_cate0_human <- tissue_pc_mean3_human[tissue_pc_mean3_human$tissue_pc_mean_count_human==0,]
tissue_pc_cate0_sum_human <- table(tissue_pc_cate0_human$Category)
tissue_pc_cate0_sum_human <- as.data.frame(tissue_pc_cate0_sum_human)
colnames(tissue_pc_cate0_sum_human) <- c("Category","0")
tissue_pc_cate1_human <- tissue_pc_mean3_human[tissue_pc_mean3_human$tissue_pc_mean_count_human==1,]
tissue_pc_cate1_sum_human <- table(tissue_pc_cate1_human$Category)
tissue_pc_cate1_sum_human <- as.data.frame(tissue_pc_cate1_sum_human)
colnames(tissue_pc_cate1_sum_human) <- c("Category","1")
tissue_pc_cate2_human <- tissue_pc_mean3_human[tissue_pc_mean3_human$tissue_pc_mean_count_human==2,]
tissue_pc_cate2_sum_human <- table(tissue_pc_cate2_human$Category)
tissue_pc_cate2_sum_human <- as.data.frame(tissue_pc_cate2_sum_human)
colnames(tissue_pc_cate2_sum_human) <- c("Category","2")
tissue_pc_cate3_human <- tissue_pc_mean3_human[tissue_pc_mean3_human$tissue_pc_mean_count_human==3,]
tissue_pc_cate3_sum_human <- table(tissue_pc_cate3_human$Category)
tissue_pc_cate3_sum_human <- as.data.frame(tissue_pc_cate3_sum_human)
colnames(tissue_pc_cate3_sum_human) <- c("Category","3")
tissue_pc_cate4_human <- tissue_pc_mean3_human[tissue_pc_mean3_human$tissue_pc_mean_count_human==4,]
tissue_pc_cate4_sum_human <- table(tissue_pc_cate4_human$Category)
tissue_pc_cate4_sum_human <- as.data.frame(tissue_pc_cate4_sum_human)
colnames(tissue_pc_cate4_sum_human) <- c("Category","4")
tissue_pc_cate5_human <- tissue_pc_mean3_human[tissue_pc_mean3_human$tissue_pc_mean_count_human==5,]
tissue_pc_cate5_sum_human <- table(tissue_pc_cate5_human$Category)
tissue_pc_cate5_sum_human <- as.data.frame(tissue_pc_cate5_sum_human)
colnames(tissue_pc_cate5_sum_human) <- c("Category","5")
tissue_pc_cate6_human <- tissue_pc_mean3_human[tissue_pc_mean3_human$tissue_pc_mean_count_human==6,]
tissue_pc_cate6_sum_human <- table(tissue_pc_cate6_human$Category)

tissue_pc_cate6_sum_human <- as.data.frame(tissue_pc_cate6_sum_human)
colnames(tissue_pc_cate6_sum_human) <- c("Category","6")
tissue_pc_cate7_human <- tissue_pc_mean3_human[tissue_pc_mean3_human$tissue_pc_mean_count_human==7,]
tissue_pc_cate7_sum_human <- table(tissue_pc_cate7_human$Category)
tissue_pc_cate7_sum_human <- as.data.frame(tissue_pc_cate7_sum_human)
colnames(tissue_pc_cate7_sum_human) <- c("Category","7")
tissue_pc_cate8_human <- tissue_pc_mean3_human[tissue_pc_mean3_human$tissue_pc_mean_count_human==8,]
tissue_pc_cate8_sum_human <- table(tissue_pc_cate8_human$Category)
tissue_pc_cate8_sum_human <- as.data.frame(tissue_pc_cate8_sum_human)
colnames(tissue_pc_cate8_sum_human) <- c("Category","8")
tissue_pc_cate9_human <- tissue_pc_mean3_human[tissue_pc_mean3_human$tissue_pc_mean_count_human==9,]
tissue_pc_cate9_sum_human <- table(tissue_pc_cate9_human$Category)
tissue_pc_cate9_sum_human <- as.data.frame(tissue_pc_cate9_sum_human)
colnames(tissue_pc_cate9_sum_human) <- c("Category","9")

tissue_pc_cate10_human <- tissue_pc_mean3_human[tissue_pc_mean3_human$tissue_pc_mean_count_human==10,]
tissue_pc_cate10_sum_human <- table(tissue_pc_cate10_human$Category)
tissue_pc_cate10_sum_human <- as.data.frame(tissue_pc_cate10_sum_human)
colnames(tissue_pc_cate10_sum_human) <- c("Category","10")
tissue_pc_cate11_human <- tissue_pc_mean3_human[tissue_pc_mean3_human$tissue_pc_mean_count_human==11,]
tissue_pc_cate11_sum_human <- table(tissue_pc_cate11_human$Category)
tissue_pc_cate11_sum_human <- as.data.frame(tissue_pc_cate11_sum_human)
colnames(tissue_pc_cate11_sum_human) <- c("Category","11")

tissue_pc_cate12_human <- tissue_pc_mean3_human[tissue_pc_mean3_human$tissue_pc_mean_count_human==12,]
tissue_pc_cate12_sum_human <- table(tissue_pc_cate12_human$Category)
tissue_pc_cate12_sum_human <- as.data.frame(tissue_pc_cate12_sum_human)
colnames(tissue_pc_cate12_sum_human) <- c("Category","12")
tissue_pc_cate13_human <- tissue_pc_mean3_human[tissue_pc_mean3_human$tissue_pc_mean_count_human==13,]
tissue_pc_cate13_sum_human <- table(tissue_pc_cate13_human$Category)
tissue_pc_cate13_sum_human <- as.data.frame(tissue_pc_cate13_sum_human)
colnames(tissue_pc_cate13_sum_human) <- c("Category","13")
tissue_pc_cate14_human <- tissue_pc_mean3_human[tissue_pc_mean3_human$tissue_pc_mean_count_human==14,]
tissue_pc_cate14_sum_human <- table(tissue_pc_cate14_human$Category)
tissue_pc_cate14_sum_human <- as.data.frame(tissue_pc_cate14_sum_human)
colnames(tissue_pc_cate14_sum_human) <- c("Category","14")
tissue_pc_cate15_human <- tissue_pc_mean3_human[tissue_pc_mean3_human$tissue_pc_mean_count_human==15,]
tissue_pc_cate15_sum_human <- table(tissue_pc_cate15_human$Category)
tissue_pc_cate15_sum_human <- as.data.frame(tissue_pc_cate15_sum_human)
colnames(tissue_pc_cate15_sum_human) <- c("Category","15")
tissue_pc_cate16_human <- tissue_pc_mean3_human[tissue_pc_mean3_human$tissue_pc_mean_count_human==16,]
tissue_pc_cate16_sum_human <- table(tissue_pc_cate16_human$Category)
tissue_pc_cate16_sum_human <- as.data.frame(tissue_pc_cate16_sum_human)
colnames(tissue_pc_cate16_sum_human) <- c("Category","16")
tissue_pc_cate17_human <- tissue_pc_mean3_human[tissue_pc_mean3_human$tissue_pc_mean_count_human==17,]
tissue_pc_cate17_sum_human <- table(tissue_pc_cate17_human$Category)
tissue_pc_cate17_sum_human <- as.data.frame(tissue_pc_cate17_sum_human)
colnames(tissue_pc_cate17_sum_human) <- c("Category","17")
tissue_pc_cate18_human <- tissue_pc_mean3_human[tissue_pc_mean3_human$tissue_pc_mean_count_human==18,]
tissue_pc_cate18_sum_human <- table(tissue_pc_cate18_human$Category)
tissue_pc_cate18_sum_human <- as.data.frame(tissue_pc_cate18_sum_human)
colnames(tissue_pc_cate18_sum_human) <- c("Category","18")
tissue_pc_cate19_human <- tissue_pc_mean3_human[tissue_pc_mean3_human$tissue_pc_mean_count_human==19,]
tissue_pc_cate19_sum_human <- table(tissue_pc_cate19_human$Category)
tissue_pc_cate19_sum_human <- as.data.frame(tissue_pc_cate19_sum_human)
colnames(tissue_pc_cate19_sum_human) <- c("Category","19")
tissue_pc_cate20_human <- tissue_pc_mean3_human[tissue_pc_mean3_human$tissue_pc_mean_count_human==20,]
tissue_pc_cate20_sum_human <- table(tissue_pc_cate20_human$Category)
tissue_pc_cate20_sum_human <- as.data.frame(tissue_pc_cate20_sum_human)
colnames(tissue_pc_cate20_sum_human) <- c("Category","20")



tissue_pc_merge1_human <- merge(tissue_pc_cate0_sum_human, tissue_pc_cate1_sum_human, all = T)
tissue_pc_merge2_human <- merge(tissue_pc_merge1_human, tissue_pc_cate2_sum_human, all = T)
tissue_pc_merge3_human <- merge(tissue_pc_merge2_human, tissue_pc_cate3_sum_human, all = T)
tissue_pc_merge4_human <- merge(tissue_pc_merge3_human, tissue_pc_cate4_sum_human, all = T)
tissue_pc_merge5_human <- merge(tissue_pc_merge4_human, tissue_pc_cate5_sum_human, all = T)
tissue_pc_merge6_human <- merge(tissue_pc_merge5_human, tissue_pc_cate6_sum_human, all = T)
tissue_pc_merge7_human <- merge(tissue_pc_merge6_human, tissue_pc_cate7_sum_human, all = T)
tissue_pc_merge8_human <- merge(tissue_pc_merge7_human, tissue_pc_cate8_sum_human, all = T)
tissue_pc_merge9_human <- merge(tissue_pc_merge8_human, tissue_pc_cate9_sum_human, all = T)
tissue_pc_merge10_human <- merge(tissue_pc_merge9_human, tissue_pc_cate10_sum_human, all = T)
tissue_pc_merge11_human <- merge(tissue_pc_merge10_human, tissue_pc_cate11_sum_human, all = T)
tissue_pc_merge12_human <- merge(tissue_pc_merge11_human, tissue_pc_cate12_sum_human, all = T)
tissue_pc_merge13_human <- merge(tissue_pc_merge12_human, tissue_pc_cate13_sum_human, all = T)
tissue_pc_merge14_human <- merge(tissue_pc_merge13_human, tissue_pc_cate14_sum_human, all = T)
tissue_pc_merge15_human <- merge(tissue_pc_merge14_human, tissue_pc_cate15_sum_human, all = T)
tissue_pc_merge16_human <- merge(tissue_pc_merge15_human, tissue_pc_cate16_sum_human, all = T)
tissue_pc_merge17_human <- merge(tissue_pc_merge16_human, tissue_pc_cate17_sum_human, all = T)
tissue_pc_merge18_human <- merge(tissue_pc_merge17_human, tissue_pc_cate18_sum_human, all = T)
tissue_pc_merge19_human <- merge(tissue_pc_merge18_human, tissue_pc_cate19_sum_human, all = T)
tissue_pc_merge20_human <- merge(tissue_pc_merge19_human, tissue_pc_cate20_sum_human, all = T)

tissue_pc_merge20_human[is.na(tissue_pc_merge20_human)]=0
rownames(tissue_pc_merge20_human) <- tissue_pc_merge20_human$Category
tissue_pc_merge20_human <- tissue_pc_merge20_human[,-1]

tissue_pc_merge20_human
tissue_pc_merge20

save(tissue_pc_merge20_human,tissue_pc_merge20,file="tissue_specifiy.Rdata")

load("tissue_specifiy.Rdata")
tissue_pc_merge20
tissue_pc_merge20_human
##########barplot
tiff(file="gene_pc_cattle.tiff", width = 8, height = 7, units = "cm", res = 600, pointsize = 8,compression= "lzw")
barplot(as.matrix(tissue_pc_merge20/1000),ylim = c(0,15),
        col=brewer.pal(9,"Oranges"),
        legend = rownames(tissue_pc_merge20), yaxt='n', line=-0.5,
        cex.axis = 0.7, cex.names = 0.7, las = 1, angle = 180,
        width = 0.5, space = 0.5, beside = F,
        args.legend = list(x = 6, y= 10.654, bty = 'n', cex = 0.8, y.intersp = 1,
                           x.intersp = 0.7, text.width = 1))
axis(2, col.axis="orange", cex.axis = 0.6, lwd.ticks = 1, tck=-0.01, padj=1.5, line=0.5)
mtext('Number of genes (thousands) in cattle', side = 2, cex = 0.8,  line = 2)
mtext('Number of tissues in cattle', side = 1, cex = 0.8,  line = 2)

dev.off()

tiff(file="gene_pc_human.tiff", width = 8, height = 7, units = "cm", res = 600, pointsize = 8,compression= "lzw")
barplot(as.matrix(tissue_pc_merge20/1000),ylim = c(0,15),
        col=brewer.pal(9,"YlGn"),
        legend = rownames(tissue_pc_merge20), yaxt='n', line=-0.5,
        cex.axis = 0.7, cex.names = 0.7, las = 1, angle = 180,
        width = 0.5, space = 0.5, beside = F,
        args.legend = list(x = 6, y= 10.654, bty = 'n', cex = 0.8, y.intersp = 1,
                           x.intersp = 0.7, text.width = 1))
axis(2, col.axis="darkseagreen", cex.axis = 0.6, lwd.ticks = 1, tck=-0.01, padj=1.5, line=0.5)
mtext('Number of genes (thousands) in human', side = 2, cex = 0.8,  line = 2)
mtext('Number of tissues in human', side = 1, cex = 0.8,  line = 2)
dev.off()


######################
####make some improve 
#######################
#save(tissue_pc_merge20_human,tissue_pc_merge20,file="tissue_specifiy.Rdata")
load("tissue_specifiy.Rdata")
tissue_pc_merge20
tissue_pc_merge20_human
library(RColorBrewer)
##########barplot
tiff(file="gene_pc_cattle.tiff", width = 8, height = 8, units = "cm", res = 600, pointsize = 8,compression= "lzw")

barplot(as.matrix(tissue_pc_merge20/1000),ylim = c(0,15),
        col=brewer.pal(9,"Oranges"),yaxt='n', line=-0.5,
        cex.axis = 1, cex.names = 1, las = 1, angle = 180,
        width = 0.5, space = 0.5, beside = F)


axis(2, col.axis="black", cex.axis = 1, lwd.ticks = 1, tck=-0.01, padj=1.5, line=0.5)
mtext('Number of genes (thousands) in cattle', side = 2, cex = 1.5,  line = 2)
mtext('Number of tissues in cattle', side = 1, cex = 1.2,  line = 2)

dev.off()


tiff(file="gene_pc_human2.tiff", width = 8, height = 8, units = "cm", res = 600, pointsize = 8,compression= "lzw")
barplot(as.matrix(tissue_pc_merge20_human/1000),ylim = c(0,15),
        col=brewer.pal(9,"YlGn"),yaxt='n', line=-0.5,
        cex.axis = 1, cex.names = 1, las = 1, angle = 180,
        width = 0.5, space = 0.5, beside = F,
)
axis(2, col.axis="black", cex.axis = 1, lwd.ticks = 1, tck=-0.01, padj=1.5, line=0.5)
mtext('Number of genes (thousands) in human', side = 2, cex = 1.5,  line = 2)
mtext('Number of tissues in human', side = 1, cex = 1.2,  line = 2)

dev.off()
