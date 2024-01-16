# Fig S3
library(RColorBrewer)
library(ggplot2)
library(cutpointr)
library(plyr)
library(tidyverse)
library(ggpubr)
library(cocor)
library(pROC)

# Read in composite data 
all<-read.table("Data/Fig3/all_prs.txt",sep="\t",header=T)
somatic<-read.table("Data/suppFigs/somatic_all_pts.txt",sep="\t",header=T)
df<-merge(all,somatic[,c("normal.WXS.id","TMB","zTMB","cTMB","pTMB")],by="normal.WXS.id")
ici<-df
ici$Response<-ifelse(ici$pheno==2,"R","NR")

cols<-c("#1A281F","#B8D4E3","#8A8E91","#855A5C")
cols2<-brewer.pal(3,"Dark2")

# plot each correlation of each score type with TMB
ggplot(ici, aes(x = TMB, y = Composite_PRS,color=Response)) +
  geom_jitter(width = 0, height = 0,alpha=0.7,size=2)+
  stat_cor(method = "pearson",show.legend = F,inherit.aes = F,aes(x = TMB, y = Composite_PRS),size=6)+
  scale_x_log10() + ylab(label = "Composite IC-Index") +xlab("Log10 TMB")+
  #geom_smooth(method = "lm",inherit.aes = F,aes(x = TMB, y = NumLargeSubClones)) +
  theme_light(base_size=16) + scale_color_manual(values=cols2[2:1]) #+ xlab(label = "Germline IC-Index") + ylab(label = "Somatic IC-Index")
#ggsave("/Users/tjsears/Code/GermTime/figs_july/figS3/Composite_TMB_corr.pdf",width = 8,height = 6)

ggplot(ici, aes(x = TMB, y = Somatic_PRS,color=Response)) +
  geom_jitter(width = 0, height = 0,alpha=0.7,size=2)+
  stat_cor(method = "pearson",show.legend = F,inherit.aes = F,aes(x = TMB, y = Somatic_PRS),size=6)+
  scale_x_log10() + ylab(label = "Somatic IC-Index") +xlab("Log10 TMB")+
  theme_light(base_size=16) + scale_color_manual(values=cols2[2:1]) #+ xlab(label = "Germline IC-Index") + ylab(label = "Somatic IC-Index")
#ggsave("/Users/tjsears/Code/GermTime/figs_july/figS3/Somatic_TMB_corr.pdf",width = 8,height = 6)

ggplot(ici, aes(x = TMB, y = Germline_PRS,color=Response)) +
  geom_jitter(width = 0, height = 0,alpha=0.7,size=2)+
  stat_cor(method = "pearson",show.legend = F,inherit.aes = F,aes(x = TMB, y = Germline_PRS),size=6)+
  scale_x_log10() + ylab(label = "Germline IC-Index") +xlab("Log10 TMB")+
  theme_light(base_size=16) + scale_color_manual(values=cols2[2:1]) #+ xlab(label = "Germline IC-Index") + ylab(label = "Somatic IC-Index")
#ggsave("/Users/tjsears/Code/GermTime/figs_july/figS3/Germline_TMB_corr.pdf",width = 8,height = 6)


# plot each cohort against 

# TMB
# age
# sex
library(pROC)

ici_vanallen<-ici[ici$study_cancer_x=="vanallen",]

rocobj_vanallen_tmb <- pROC::roc(ici_vanallen$pheno,ici_vanallen$TMB)
auc_tmb <- round(pROC::auc(ici_vanallen$pheno,ici_vanallen$TMB),4)

rocobj_vanallen_age <- pROC::roc(ici_vanallen$pheno,ici_vanallen$Age_y)
auc_age <- round(pROC::auc(ici_vanallen$pheno,ici_vanallen$Age_y),4)

rocobj_vanallen_gender <- pROC::roc(ici_vanallen$pheno,ici_vanallen$Gender)
auc_gender <- round(pROC::auc(ici_vanallen$pheno,ici_vanallen$Gender),4)

rocobj_vanallen_composite <- pROC::roc(ici_vanallen$pheno,ici_vanallen$Composite_PRS)
auc_composite <- round(pROC::auc(ici_vanallen$pheno,ici_vanallen$Composite_PRS),4)

ggroc(list(call_roc_name_1 = rocobj_vanallen_tmb, call_roc_name_2 = rocobj_vanallen_age,call_roc_name_3=rocobj_vanallen_gender,call_roc_name_4=rocobj_vanallen_composite),size=1.2,aes = c("color"))+
  theme_minimal(base_size=18)+
  theme(legend.position = c(0.7, 0.3)) + scale_color_discrete(labels=c(paste("TMB\n AUC=",auc_tmb,sep=""),
                                                                       paste("Age\n AUC=",auc_age,sep=""),
                                                                       paste("Gender\n AUC=",auc_gender,sep=""),
                                                                       paste("Composite IC-Index\n AUC=",auc_composite,sep="")),type = cols) +
  labs(color=NULL) + ggtitle("Vanallen (N=110)")
#ggsave("/Users/tjsears/Code/GermTime/figs_july/figS3/ICIndex_vs_TMB_vanallen.pdf",width = 6,height = 6)

ici_miao<-ici[ici$study_cancer_x=="miao",]

rocobj_miao_tmb <- pROC::roc(ici_miao$pheno,ici_miao$TMB)
auc_tmb <- round(pROC::auc(ici_miao$pheno,ici_miao$TMB),4)

rocobj_miao_age <- pROC::roc(ici_miao$pheno,ici_miao$Age_y)
auc_age <- round(pROC::auc(ici_miao$pheno,ici_miao$Age_y),4)

rocobj_miao_gender <- pROC::roc(ici_miao$pheno,ici_miao$Gender)
auc_gender <- round(pROC::auc(ici_miao$pheno,ici_miao$Gender),4)

rocobj_miao_composite <- pROC::roc(ici_miao$pheno,ici_miao$Composite_PRS)
auc_composite <- round(pROC::auc(ici_miao$pheno,ici_miao$Composite_PRS),4)

ggroc(list(call_roc_name_1 = rocobj_miao_tmb, call_roc_name_2 = rocobj_miao_age,call_roc_name_3=rocobj_miao_gender,call_roc_name_4=rocobj_miao_composite),size=1.2,aes = c("color"))+
  theme_minimal(base_size=18)+
  theme(legend.position = c(0.7, 0.3)) + scale_color_discrete(labels=c(paste("TMB\n AUC=",auc_tmb,sep=""),
                                                                       paste("Age\n AUC=",auc_age,sep=""),
                                                                       paste("Gender\n AUC=",auc_gender,sep=""),
                                                                       paste("Composite IC-Index\n AUC=",auc_composite,sep="")),type = cols) +
  labs(color=NULL) + ggtitle("Miao (N=70)")
#ggsave("/Users/tjsears/Code/GermTime/figs_july/figS3/ICIndex_vs_TMB_miao.pdf",width = 6,height = 6)

ici_rizvi<-ici[ici$study_cancer_x=="rizvi",]

rocobj_rizvi_tmb <- pROC::roc(ici_rizvi$pheno,ici_rizvi$TMB)
auc_tmb <- round(pROC::auc(ici_rizvi$pheno,ici_rizvi$TMB),4)

rocobj_rizvi_age <- pROC::roc(ici_rizvi$pheno,ici_rizvi$Age_y)
auc_age <- round(pROC::auc(ici_rizvi$pheno,ici_rizvi$Age_y),4)

rocobj_rizvi_gender <- pROC::roc(ici_rizvi$pheno,ici_rizvi$Gender)
auc_gender <- round(pROC::auc(ici_rizvi$pheno,ici_rizvi$Gender),4)

rocobj_rizvi_composite <- pROC::roc(ici_rizvi$pheno,ici_rizvi$Composite_PRS)
auc_composite <- round(pROC::auc(ici_rizvi$pheno,ici_rizvi$Composite_PRS),4)

ggroc(list(call_roc_name_1 = rocobj_rizvi_tmb, call_roc_name_2 = rocobj_rizvi_age,call_roc_name_3=rocobj_rizvi_gender,call_roc_name_4=rocobj_rizvi_composite),size=1.2,aes = c("color"))+
  theme_minimal(base_size=18)+
  theme(legend.position = c(0.7, 0.3)) + scale_color_discrete(labels=c(paste("TMB\n AUC=",auc_tmb,sep=""),
                                                                       paste("Age\n AUC=",auc_age,sep=""),
                                                                       paste("Gender\n AUC=",auc_gender,sep=""),
                                                                       paste("Composite IC-Index\n AUC=",auc_composite,sep="")),type = cols) +
  labs(color=NULL) + ggtitle("Rizvi (N=34)")
#ggsave("/Users/tjsears/Code/GermTime/figs_july/figS3/ICIndex_vs_TMB_rizvi.pdf",width = 6,height = 6)


############## RNA ###########

# PDL1
# CTLA4 
# LAG3

# and RNA based comparison
rna<-read.table("Data/suppFigs/rna_combat_checkpointsOnly.txt",sep="\t",header=T)

rna_t<-as.data.frame(t(rna[c("PDCD1","CTLA4","CD274","HAVCR2","TIGIT","LAG3"),]))
ici_rna<-merge(ici,rna_t,by.x = "tumor.RNA.id",by.y=0)

ici_miao<-ici_rna[ici_rna$study_cancer_x=="miao",]

rocobj_miao_tmb <- pROC::roc(ici_miao$pheno,ici_miao$CD274)
auc_tmb <- round(pROC::auc(ici_miao$pheno,ici_miao$CD274),4)

rocobj_miao_age <- pROC::roc(ici_miao$pheno,ici_miao$CTLA4)
auc_age <- round(pROC::auc(ici_miao$pheno,ici_miao$CTLA4),4)

rocobj_miao_gender <- pROC::roc(ici_miao$pheno,ici_miao$PDCD1)
auc_gender <- round(pROC::auc(ici_miao$pheno,ici_miao$PDCD1),4)

rocobj_miao_LAG3 <- pROC::roc(ici_miao$pheno,ici_miao$LAG3)
auc_LAG3 <- round(pROC::auc(ici_miao$pheno,ici_miao$LAG3),4)

rocobj_miao_composite <- pROC::roc(ici_miao$pheno,ici_miao$Composite_PRS)
auc_composite <- round(pROC::auc(ici_miao$pheno,ici_miao$Composite_PRS),4)

ggroc(list(call_roc_name_1 = rocobj_miao_tmb, call_roc_name_2 = rocobj_miao_age,call_roc_name_3=rocobj_miao_gender,call_roc_name_4=rocobj_miao_LAG3,call_roc_name_5=rocobj_miao_composite),size=1.2,aes = c("color"))+
  theme_minimal(base_size=18) +
  theme(legend.position = c(0.7, 0.3)) + scale_color_discrete(labels=c(paste("PD-L1\n AUC=",auc_tmb,sep=""),
                                                                       paste("CTLA4\n AUC=",auc_age,sep=""),
                                                                       paste("PDCD1\n AUC=",auc_gender,sep=""),
                                                                       paste("LAG3\n AUC=",auc_LAG3,sep=""),
                                                                       paste("Composite IC-Index\n AUC=",auc_composite,sep="")),type = c("darkblue",cols)) +
  labs(color=NULL) + ggtitle("Miao (N=33)")
#ggsave("/Users/tjsears/Code/GermTime/figs_july/figS3/ICIndex_vs_RNA_miao.pdf",width = 6,height = 6)

ici_vanallen<-ici_rna[ici_rna$study_cancer_x=="vanallen",]

rocobj_vanallen_tmb <- pROC::roc(ici_vanallen$pheno,ici_vanallen$CD274)
auc_tmb <- round(pROC::auc(ici_vanallen$pheno,ici_vanallen$CD274),4)

rocobj_vanallen_age <- pROC::roc(ici_vanallen$pheno,ici_vanallen$CTLA4)
auc_age <- round(pROC::auc(ici_vanallen$pheno,ici_vanallen$CTLA4),4)

rocobj_vanallen_gender <- pROC::roc(ici_vanallen$pheno,ici_vanallen$PDCD1)
auc_gender <- round(pROC::auc(ici_vanallen$pheno,ici_vanallen$PDCD1),4)

rocobj_vanallen_LAG3 <- pROC::roc(ici_vanallen$pheno,ici_vanallen$LAG3)
auc_LAG3 <- round(pROC::auc(ici_vanallen$pheno,ici_vanallen$LAG3),4)

rocobj_vanallen_composite <- pROC::roc(ici_vanallen$pheno,ici_vanallen$Composite_PRS)
auc_composite <- round(pROC::auc(ici_vanallen$pheno,ici_vanallen$Composite_PRS),4)

ggroc(list(call_roc_name_1 = rocobj_vanallen_tmb, call_roc_name_2 = rocobj_vanallen_age,call_roc_name_3=rocobj_vanallen_gender,call_roc_name_4=rocobj_vanallen_LAG3,call_roc_name_5=rocobj_vanallen_composite),size=1.2,aes = c("color"))+
  theme_minimal(base_size=18) +
  theme(legend.position = c(0.7, 0.3)) + scale_color_discrete(labels=c(paste("PD-L1\n AUC=",auc_tmb,sep=""),
                                                                       paste("CTLA4\n AUC=",auc_age,sep=""),
                                                                       paste("PDCD1\n AUC=",auc_gender,sep=""),
                                                                       paste("LAG3\n AUC=",auc_LAG3,sep=""),
                                                                       paste("Composite IC-Index\n AUC=",auc_composite,sep="")),type = c("darkblue",cols)) +
  labs(color=NULL) + ggtitle("Vanallen (N=38)")
#ggsave("/Users/tjsears/Code/GermTime/figs_july/figS3/ICIndex_vs_RNA_vanallen.pdf",width = 6,height = 6)


# then do correlation with TMB? is that helpful? I suspect it will be extremely uncorrelated

# randomly generated scores with a normal distribution centered on 5
# make a histogram of the test statistic?

#simulated normal dist given mean and sd of composite PRS scores


population_mean <- mean(ici$Composite_PRS)
population_sd <- sd( ici$Composite_PRS)

ici_R<-ici[ici$Response=="R",]
ici_NR<-ici[ici$Response=="NR",]

#define upper and lower bound
lower_bound <- population_mean - population_sd
upper_bound <- population_mean + population_sd

#Create a sequence of 1000 x values based on population mean and standard deviation
x <- seq(0, 10, length = 1000) * population_sd + population_mean

#create a vector of values that shows the height of the probability distribution
#for each value in x
y <- dnorm(x, population_mean, population_sd)

ggplot(data.frame(x = c(0, 10)), aes(x = x)) +
  stat_function(fun = dnorm,args=c(population_mean,population_sd)) + theme_minimal(base_size=14) +
  geom_segment(aes(x = mean(ici_R$Composite_PRS), xend = mean(ici_R$Composite_PRS), y = -Inf, yend = Inf),
               linetype = 2)  +
  geom_segment(aes(x = mean(ici_NR$Composite_PRS), xend = mean(ici_NR$Composite_PRS), y = -Inf, yend = Inf),
               linetype = 2)  
  
ggplot(ici, aes(f0)) + 
  geom_histogram(data = lowf0, fill = "red", alpha = 0.2) + 
  geom_histogram(data = mediumf0, fill = "blue", alpha = 0.2) +
  geom_histogram(data = highf0, fill = "green", alpha = 0.2)

pval<-t.test(ici_R$Composite_PRS,ici_NR$Composite_PRS)$p.value
ggplot(ici, aes(x=Composite_PRS, fill=Response)) + ylab(label="Count")+xlab(label="IC-Index")+
  geom_histogram(alpha=0.3, position="identity",bins = 20)  + theme_minimal(base_size=16) +
  #geom_segment(aes(x = median(ici_R$Composite_PRS), xend = median(ici_R$Composite_PRS), y = -Inf, yend = Inf),
  #             linetype = 2)  +
  #geom_segment(aes(x = median(ici_NR$Composite_PRS), xend = median(ici_NR$Composite_PRS), y = -Inf, yend = Inf),
  #             linetype = 2)  + 
  scale_fill_manual(values = cols2[2:1]) + ggtitle(label="All Test Pts. IC-Index",subtitle = paste("T-Test pval =",round(pval,9)))
#ggsave("/Users/tjsears/Code/GermTime/figs_july/figS3/ICIndex_Histogram.pdf",width = 7,height = 5)

### Purity vs Ploidy corr plot
#read in data w putiy and ploidy
all_plot<-read.table("Data/suppFigs/allSomaticFeatures.txt",sep="\t",header=T)
all_plot$Response<-ifelse(all_plot$pheno==2,"R","NR")

ggplot(all_plot, aes(x = purity, y = tumorPloidy,color=Response)) +
  geom_jitter(width = 0, height = 0,alpha=0.7,size=2)+
  stat_cor(method = "pearson",show.legend = F,inherit.aes = F,aes(x = purity, y = tumorPloidy),size=6)+
  ylab(label = "Ploidy") +xlab("Purity")+
  #geom_smooth(method = "lm",inherit.aes = F,aes(x = TMB, y = NumLargeSubClones)) +
  theme_light(base_size=16) + scale_color_manual(values=cols2[2:1]) #+ xlab(label = "Germline IC-Index") + ylab(label = "Somatic IC-Index")
#ggsave("/Users/tjsears/Code/GermTime/figs_july/figS4/PurityVSPloidy.pdf",width = 8,height = 6)



