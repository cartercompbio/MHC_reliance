# plot TCGA CD4/8 ratio and general immune heatmap

library(tidyverse)
library(dplyr)
library(splitstackshape)
library(Rtsne)
library(ggplot2)
library(umap)
library(readxl)
library(dplyr)
library(gridExtra)
library(survminer)
library(ggpubr)
library(ggfortify)
library(survival)

# read in TCGA data
tcga<-read.table("Data/suppFigs/tcga_dat.txt",sep="\t",header=T)

#rewrite this file to be much much smaller
table(tcga$type)

# apply filtering procedure you used before
plot_df<-tcga
table(plot_df$type)

#rerun here
plot_df_low<-plot_df
plot_df_low<-plot_df_low[plot_df_low$PFI.time>=30,]
plot_df_low<-plot_df_low[plot_df_low$PsuedoResponse=="R",]
plot_df_low<-plot_df_low[plot_df_low$ajcc_pathologic_tumor_stage%in%c("Stage IIB","Stage IIC","Stage IIIA","Stage IIIB","Stage IIIC","Stage III","Stage IVB","Stage IVC","Stage IV","Stage IVA"),]
plot_df_low<-plot_df_low[!is.na(plot_df_low$OS.time),]
plot_df_low<-plot_df_low[!is.na(plot_df_low$OS),]

quantiles <- quantile(plot_df_low$Ratio, probs = c(0,0.32,0.66,1)) 
plot_df_low$group <- cut(plot_df_low$Ratio, breaks = quantiles, labels = c("MHC-I Reliant","Balanced","MHC-II Reliant"))
plot_df_low$group<-factor(plot_df_low$group,levels = c("MHC-I Reliant","Balanced","MHC-II Reliant"))
plot_df_low$group[is.na(plot_df_low$group)]<-"MHC-I Reliant"
plot_df_low<-plot_df_low[plot_df_low$group=="MHC-I Reliant",]

plot_df_low<-plot_df_low[plot_df_low$T.cells.CD8>=0,]
plot_df_low$CD4_CD8_Ratio=(plot_df_low$T.cells.CD4.memory.resting+plot_df_low$T.cells.CD4.memory.activated
                           +0.01)/(plot_df_low$T.cells.CD8+0.01)
plot_df_low$TcellRatio<-plot_df_low$CD4_CD8_Ratio

quantiles <- quantile(plot_df_low$CD4_CD8_Ratio, probs = c(0,0.49,1))
plot_df_low$CD4_CD8_Ratio<-cut(plot_df_low$CD4_CD8_Ratio, breaks = quantiles, labels = c("CD4/CD8 lo","CD4/CD8 hi"))
plot_df_low$CD4_CD8_Ratio[is.na(plot_df_low$CD4_CD8_Ratio)]<-"CD4/CD8 lo"

km_trt_fit <- survfit(Surv(OS.time, OS) ~ CD4_CD8_Ratio, data=plot_df_low)
res<- pairwise_survdiff(Surv(OS.time, OS) ~ CD4_CD8_Ratio, data=plot_df_low,p.adjust.method = "none",rho=0)#, p.adjust.method = "none", rho = 0)
res$p.value<-round(res$p.value, digits = 4)
res

p1<-ggsurvplot(km_trt_fit,data=plot_df_low,size=2,censor.shape="|", censor.size = 8,
               ggtheme=(theme_minimal(base_size=18)),
               pval=F,xlim=c(0,2000),break.x.by=200,ncensor.plot.height=1,
               legend.labs=c("CD4/CD8 lo","CD4/CD8 hi"),palette=c("#FEAE72","#D95F02"),
               ,xlab="Time in Days",legend="right",ylim=c(0,1),ylab="Overall Survival (%)",surv.scale = c("percent"))+ggtitle(label="CD4/CD8 Ratio hi/lo\nMHC-I reliant pts.")
p1$plot+ggplot2::annotate("label",x=1200,y = 0.8,size=6, label = paste("P =",round(res$p.value,4)))

#ggsave("/Users/tjsears/Code/GermTime/figs_july/fig7/cd4_8_MHCII_low.pdf",width = 8,height = 6)

#rerun here
plot_df_mid<-plot_df
plot_df_mid<-plot_df_mid[plot_df_mid$PFI.time>=30,]
plot_df_mid<-plot_df_mid[plot_df_mid$PsuedoResponse=="R",]
plot_df_mid<-plot_df_mid[plot_df_mid$ajcc_pathologic_tumor_stage%in%c("Stage IIB","Stage IIC","Stage IIIA","Stage IIIB","Stage IIIC","Stage III","Stage IVB","Stage IVC","Stage IV","Stage IVA"),]
plot_df_mid<-plot_df_mid[!is.na(plot_df_mid$OS.time),]
plot_df_mid<-plot_df_mid[!is.na(plot_df_mid$OS),]

quantiles <- quantile(plot_df_mid$Ratio, probs = c(0,0.33,0.66,1))
plot_df_mid$group <- cut(plot_df_mid$Ratio, breaks = quantiles, labels = c("MHC-I Reliant","Balanced","MHC-II Reliant"))
plot_df_mid$group<-factor(plot_df_mid$group,levels = c("MHC-I Reliant","Balanced","MHC-II Reliant"))
plot_df_mid$group[is.na(plot_df_mid$group)]<-"Balanced"
plot_df_mid<-plot_df_mid[plot_df_mid$group=="Balanced",]

plot_df_mid<-plot_df_mid[plot_df_mid$T.cells.CD8>=0,]
plot_df_mid$CD4_CD8_Ratio=(plot_df_mid$T.cells.CD4.memory.resting+plot_df_mid$T.cells.CD4.memory.activated
                           +0.01)/(plot_df_mid$T.cells.CD8+0.01)
plot_df_mid$TcellRatio<-plot_df_mid$CD4_CD8_Ratio

quantiles <- quantile(plot_df_mid$CD4_CD8_Ratio, probs = c(0,0.49,1))
plot_df_mid$CD4_CD8_Ratio<-cut(plot_df_mid$CD4_CD8_Ratio, breaks = quantiles, labels = c("CD4/CD8 lo","CD4/CD8 hi"))
plot_df_mid$CD4_CD8_Ratio[is.na(plot_df_mid$CD4_CD8_Ratio)]<-"CD4/CD8 lo"

km_trt_fit <- survfit(Surv(OS.time, OS) ~ CD4_CD8_Ratio, data=plot_df_mid)
res<- pairwise_survdiff(Surv(OS.time, OS) ~ CD4_CD8_Ratio, data=plot_df_mid,p.adjust.method = "none",rho=0)#, p.adjust.method = "none", rho = 0)
res$p.value<-round(res$p.value, digits = 4)
res

p1<-ggsurvplot(km_trt_fit,data=plot_df_mid,size=2,censor.shape="|", censor.size = 8,
               ggtheme=(theme_minimal(base_size=18)),
               pval=F,xlim=c(0,2000),break.x.by=200,ncensor.plot.height=1,
               legend.labs=c("CD4/CD8 lo","CD4/CD8 hi"),palette=c("#FEAE72","#D95F02"),
               ,xlab="Time in Days",legend="right",ylim=c(0,1),ylab="Overall Survival (%)",surv.scale = c("percent"))+ggtitle(label="CD4/CD8 Ratio hi/lo\nBalanced pts.")
p1$plot+ggplot2::annotate("label",x=1200,y = 0.8,size=6, label = paste("P =",round(res$p.value,4)))

#ggsave("/Users/tjsears/Code/GermTime/figs_july/fig7/cd4_8_Balanced.pdf",width = 8,height = 6)

#rerun here
plot_df_high<-plot_df
plot_df_high<-plot_df_high[plot_df_high$PFI.time>=30,]
plot_df_high<-plot_df_high[plot_df_high$PsuedoResponse=="R",]
plot_df_high<-plot_df_high[plot_df_high$ajcc_pathologic_tumor_stage%in%c("Stage IIA","Stage IIB","Stage IIC","Stage IIIA","Stage IIIB","Stage IIIC","Stage III","Stage IVB","Stage IVC","Stage IV","Stage IVA"),]
plot_df_high<-plot_df_high[!is.na(plot_df_high$OS.time),]
plot_df_high<-plot_df_high[!is.na(plot_df_high$OS),]

quantiles <- quantile(plot_df_high$Ratio, probs = c(0,0.33,0.66,1))
plot_df_high$group <- cut(plot_df_high$Ratio, breaks = quantiles, labels = c("MHC-I Reliant","Balanced","MHC-II Reliant"))
plot_df_high$group<-factor(plot_df_high$group,levels = c("MHC-I Reliant","Balanced","MHC-II Reliant"))
plot_df_high$group[is.na(plot_df_high$group)]<-"MHC-II Reliant"
plot_df_high<-plot_df_high[plot_df_high$group=="MHC-II Reliant",]
plot_df_high<-plot_df_high[(plot_df_high$numLostAlleles+plot_df_high$Escape)>0,]

plot_df_high<-plot_df_high[plot_df_high$T.cells.CD8>=0,]
plot_df_high$CD4_CD8_Ratio=(plot_df_high$T.cells.CD4.memory.resting+plot_df_high$T.cells.CD4.memory.activated
                           +0.01)/(plot_df_high$T.cells.CD8+0.01)

quantiles <- quantile(plot_df_high$CD4_CD8_Ratio, probs = c(0,0.49,1))
plot_df_high$CD4_CD8_Ratio<-cut(plot_df_high$CD4_CD8_Ratio, breaks = quantiles, labels = c("CD4/CD8 lo","CD4/CD8 hi"))
plot_df_high$CD4_CD8_Ratio[is.na(plot_df_high$CD4_CD8_Ratio)]<-"CD4/CD8 hi"

km_trt_fit <- survfit(Surv(OS.time, OS) ~ CD4_CD8_Ratio, data=plot_df_high)
res<- pairwise_survdiff(Surv(OS.time, OS) ~ CD4_CD8_Ratio, data=plot_df_high,p.adjust.method = "none",rho=0)#, p.adjust.method = "none", rho = 0)
res$p.value<-round(res$p.value, digits = 4)
res

p1<-ggsurvplot(km_trt_fit,data=plot_df_high,size=2,censor.shape="|", censor.size = 8,
               ggtheme=(theme_minimal(base_size=18)),
               pval=F,xlim=c(0,2000),break.x.by=200,ncensor.plot.height=1,
               legend.labs=c("CD4/CD8 lo","CD4/CD8 hi"),palette=c("#A5E1CE","#66C2A5"),
               ,xlab="Time in Days",legend="right",ylim=c(0,1),ylab="Overall Survival (%)",surv.scale = c("percent"))+ggtitle(label="CD4/CD8 Ratio hi/lo\nMHC-II reliant pts.")
p1$plot+ggplot2::annotate("label",x=1200,y = 0.8,size=6, label = paste("P =",round(res$p.value,4)))

#ggsave("/Users/tjsears/Code/GermTime/figs_july/fig7/cd4_8_MHCII_high.pdf",width = 8,height = 6)


