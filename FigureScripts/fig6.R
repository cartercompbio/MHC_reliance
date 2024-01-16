#Fig 6

library(RColorBrewer)
library(ggplot2)
library(ggpubr)
options( warn = -1 )

cols<-brewer.pal(n = 3, name = "Dark2")[1:2]

tcga_phbr<-read.table("Data/Fig6/fig6v4_liu_rnaseq.txt",sep="\t",header=T)
quantiles <- quantile(tcga_phbr$Ratio, probs = c(0,0.33,0.66,1))
tcga_phbr$group <- cut(tcga_phbr$Ratio, breaks = quantiles, labels = c("PHBR2/PHBR1 Low","mid","PHBR2/PHBR1 High"))
tcga_phbr$group[is.na(tcga_phbr$group)]<-"mid"
tcga_phbr$Response<-ifelse(tcga_phbr$response_crist_sd==2,"R","NR")

ICI_phbr<-read.table("Data/Fig5/fig6v4_ici_rnaseq.txt",sep="\t",header=T)
quantiles <- quantile(ICI_phbr$Ratio, probs = c(0,0.33,0.66,1))
ICI_phbr$group <- cut(ICI_phbr$Ratio, breaks = quantiles, labels = c("PHBR2/PHBR1 Low","mid","PHBR2/PHBR1 High"))
ICI_phbr$group[is.na(ICI_phbr$group)]<-"mid"
ICI_phbr$Response<-ifelse(ICI_phbr$response_crist_sd==2,"R","NR")

#replot dataframes as faceted plot
tcga_phbr_nomid<-tcga_phbr[tcga_phbr$group!="mid",]
tcga_phbr_nomid$group<-factor(tcga_phbr_nomid$group,levels = c("PHBR2/PHBR1 Low","mid","PHBR2/PHBR1 High"))
tcga_phbr_nomid$study_cancer_x<-tcga_phbr_nomid$study_cancer

ICI_phbr_nomid<-ICI_phbr[ICI_phbr$group!="mid",]
ICI_phbr_nomid$group<-factor(ICI_phbr_nomid$group,levels = c("PHBR2/PHBR1 Low","PHBR2/PHBR1 High"))

plot_long<-rbind(tcga_phbr_nomid[,c("Response","Ratio","T.cells.CD8",'group',"T.cells.follicular.helper",'T.cells.regulatory..Tregs.',
                                    "T.cells.CD4.naive","T.cells.CD4.memory.resting","T.cells.CD4.memory.activated","study_cancer_x","NK.cells.activated","NK.cells.resting")],
                 ICI_phbr_nomid[,c("Response","Ratio","T.cells.CD8",'group',"T.cells.follicular.helper",'T.cells.regulatory..Tregs.',
                                   "T.cells.CD4.naive","T.cells.CD4.memory.resting","T.cells.CD4.memory.activated","study_cancer_x","NK.cells.activated","NK.cells.resting")])

plot_long$Cohort<-c(rep("Validation (50)",nrow(ICI_phbr_nomid)),rep("Discovery (37)",nrow(tcga_phbr_nomid)))
lowest_tertile<-quantile(plot_long$T.cells.CD8,probs=c(0,0.33,1))[2]
plot_long<-plot_long[plot_long$T.cells.CD8>=lowest_tertile,] #above lowest tertile of CD8 inf.
plot_long$TcellRatio<-(plot_long$T.cells.follicular.helper+plot_long$T.cells.CD4.memory.activated+0.01)/(plot_long$T.cells.CD8-0.01)

plot_long<-plot_long[!is.na(plot_long$TcellRatio),]

library(ggh4x)
ggplot(plot_long, aes(x = group, y = TcellRatio,fill=group)) + 
  xlab(label="CD4/CD8 Ratio by MHC Reliance and Cohort") + #scale_y_continuous(trans = "log2") +
  facet_nested(~Cohort+Response)+ylab(label="CIBERSORTx CD4/CD8 Ratio")+theme_minimal(base_size=18)+
  theme(strip.background=element_rect(color="grey30", fill="grey90")) +
  geom_boxplot(position = position_dodge(width = 1),width=0.98) +
  scale_fill_manual(values=cols[2:1],labels=c('MHC-I Reliant', 'MHC-II Reliant'))+ylim(c(0,5))+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),panel.spacing = unit(0.1, "lines"),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  stat_compare_means(comparisons = list(c("PHBR2/PHBR1 Low", "PHBR2/PHBR1 High")),label = "p.signif",size=7,
                     method="wilcox.test", method.args=list(var.equal = F,alternative="less"),label.y = 4.3,tip.length = 0.02,bracket.size = 0.5)
#ggsave("/Users/tjsears/Code/GermTime/figs_july/fig6_v5/CIBERSORTx_compopsite_test.pdf",width = 8,height = 6)


##### KM CURVES ###

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


# ICI 
## liu ##
plot_df<-tcga_phbr
plot_df$Primary_Type_bin<-ifelse(plot_df$Primary_Type=="skin",1,0)

cols<-brewer.pal(3,"Pastel1")
cols2<-brewer.pal(3,"Dark2")
cols<-c("#A4E1CE","#1A9E76")

quantiles <- quantile(plot_df$Ratio, probs = c(0,0.33,0.66,1)) 
plot_df$group <- cut(plot_df$Ratio, breaks = quantiles, labels = c("PHBR2/PHBR1 Low","mid","PHBR2/PHBR1 High"))
plot_df$group<-factor(plot_df$group,levels = c("PHBR2/PHBR1 Low","PHBR2/PHBR1 High","mid"))
plot_df$group[is.na(plot_df$group)]<-"PHBR2/PHBR1 Low"
#plot_df<-plot_df[plot_df$priorMAPKTx==0,]

plot_df_high<-plot_df[plot_df$group=="PHBR2/PHBR1 High",]
plot_df_mid<-plot_df[plot_df$group=="mid",]
plot_df_low<-plot_df[plot_df$group=="PHBR2/PHBR1 Low",]
plot_df_high<-plot_df_high[(plot_df_high$Escape+plot_df_high$numLostAlleles)>0,]

quantiles <- quantile(plot_df_high$LAG3, probs = c(0,0.5,1))
plot_df_high$Lag3cut<-cut(plot_df_high$LAG3, breaks = quantiles, labels = c("LAG3-","LAG3+"))
plot_df_high$Lag3cut[is.na(plot_df_high$Lag3cut)]<-"LAG3-"
plot_df_high$Lag3Bin<-ifelse(plot_df_high$Lag3cut=="LAG3+",2,1)
plot_df_high$acc<-ifelse(plot_df_high$Lag3Bin==plot_df_high$pheno,1,0)

quantiles <- quantile(plot_df_low$LAG3, probs = c(0,0.5,1))
plot_df_low$Lag3cut<-cut(plot_df_low$LAG3, breaks = quantiles, labels = c("LAG3-","LAG3+"))
plot_df_low$Lag3cut[is.na(plot_df_low$Lag3cut)]<-"LAG3-"
plot_df_low$Lag3Bin<-ifelse(plot_df_low$Lag3cut=="LAG3+",2,1)
plot_df_low$acc<-ifelse(plot_df_low$Lag3Bin==plot_df_low$pheno,1,0)

km_trt_fit <- survfit(Surv(OS.time, OS_x) ~ Lag3cut, data=plot_df_high)
res<- pairwise_survdiff(Surv(OS.time, OS_x) ~ Lag3cut, data=plot_df_high,p.adjust.method = "none",rho=0)#, p.adjust.method = "none", rho = 0)
res$p.value<-round(res$p.value, digits = 4)
res

p1<-ggsurvplot(km_trt_fit,data=plot_df_high,size=2,censor.shape="|", censor.size = 8,legend=c(0.2,0.15),
               font.x=18,font.y=18,font.legend=18,font.title=18,
               pval=F,xlim=c(0,1600),break.x.by=200,ncensor.plot.height=1,
               legend.labs=c("LAG3- (13)","LAG3+ (12)"),palette=cols,legend.title="",
               ,xlab="Time in Days",ylim=c(0,1),ylab="Overall Survival (%)",surv.scale = c("percent"))+ggtitle(label="Validation\nLAG3 +/- MHC-II Reliant")

p1$plot+ggplot2::annotate("label",x=1400,y = 0.8,size=6, label = paste("P =",round(res$p.value,4)))

#ggsave("/Users/tjsears/Code/GermTime/figs_july/fig6_v5/LiuKM_MHCII_high.pdf",width = 6,height = 6)


km_trt_fit <- survfit(Surv(OS.time, OS_x) ~ Lag3cut, data=plot_df_low)
res<- pairwise_survdiff(Surv(OS.time, OS_x) ~ Lag3cut, data=plot_df_low,p.adjust.method = "none",rho=0)#, p.adjust.method = "none", rho = 0)
res$p.value<-round(res$p.value, digits = 4)
res

p1<-ggsurvplot(km_trt_fit,data=plot_df_low,size=2,censor.shape="|", censor.size = 8,legend=c(0.2,0.15),
               font.x=18,font.y=18,font.legend=18,font.title=18,
               pval=F,xlim=c(0,1600),break.x.by=200,ncensor.plot.height=1,
               legend.labs=c("LAG3- (13)","LAG3+ (12)"),palette=cols,legend.title="",
               ,xlab="Time in Days",ylim=c(0,1),ylab="Overall Survival (%)",surv.scale = c("percent"))+ggtitle(label="Validation\nLAG3 +/- MHC-II Reliant")
p1$plot+ggplot2::annotate("label",x=1300,y = 0.8,size=6, label = paste("P =",round(res$p.value,4)))
#ggsave("/Users/tjsears/Code/GermTime/figs_july/fig6_v5/LiuKM_MHCII_low.pdf",width = 6,height = 6)



# Hazzy plot
#hazzy plot liu
plot_df<-tcga_phbr
plot_df$Primary_Type_bin<-ifelse(plot_df$Primary_Type=="skin",1,0)
quantiles <- quantile(plot_df$Ratio, probs = c(0,0.33,0.66,1))
plot_df$group <- cut(plot_df$Ratio, breaks = quantiles, labels = c("PHBR2/PHBR1 Low","mid","PHBR2/PHBR1 High"))
plot_df$group<-factor(plot_df$group,levels = c("PHBR2/PHBR1 Low","PHBR2/PHBR1 High","mid"))
plot_df$group[is.na(plot_df$group)]<-"mid"

plot_df_high<-plot_df[plot_df$group=="PHBR2/PHBR1 High",]
plot_df_mid<-plot_df[plot_df$group=="mid",]
plot_df_low<-plot_df[plot_df$group=="PHBR2/PHBR1 Low",]
plot_df_high<-plot_df_high[(plot_df_high$Escape+plot_df_high$numLostAlleles)>0,]

doMedians=T
if (doMedians){
  
  quantiles <- quantile(plot_df_high$LAG3, probs = c(0,0.5,1))
  plot_df_high$LAG3<-cut(plot_df_high$LAG3, breaks = quantiles, labels = c("LAG3-","LAG3+"))
  plot_df_high$LAG3[is.na(plot_df_high$LAG3)]<-"LAG3-"
  
  quantiles <- quantile(plot_df_mid$LAG3, probs = c(0,0.5,1))
  plot_df_mid$LAG3<-cut(plot_df_mid$LAG3, breaks = quantiles, labels = c("LAG3-","LAG3+"))
  plot_df_mid$LAG3[is.na(plot_df_mid$LAG3)]<-"LAG3-"
  
  quantiles <- quantile(plot_df_low$LAG3, probs = c(0,0.5,1))
  plot_df_low$LAG3<-cut(plot_df_low$LAG3, breaks = quantiles, labels = c("LAG3-","LAG3+"))
  plot_df_low$LAG3[is.na(plot_df_low$LAG3)]<-"LAG3-"
  
  quantiles <- quantile(plot_df_high$CTLA4, probs = c(0,0.5,1))
  plot_df_high$CTLA4<-cut(plot_df_high$CTLA4, breaks = quantiles, labels = c("CTLA4-","CTLA4+"))
  plot_df_high$CTLA4[is.na(plot_df_high$CTLA4)]<-"CTLA4-"
  
  quantiles <- quantile(plot_df_mid$CTLA4, probs = c(0,0.5,1))
  plot_df_mid$CTLA4<-cut(plot_df_mid$CTLA4, breaks = quantiles, labels = c("CTLA4-","CTLA4+"))
  plot_df_mid$CTLA4[is.na(plot_df_mid$CTLA4)]<-"CTLA4-"
  
  quantiles <- quantile(plot_df_low$CTLA4, probs = c(0,0.5,1))
  plot_df_low$CTLA4<-cut(plot_df_low$CTLA4, breaks = quantiles, labels = c("CTLA4-","CTLA4+"))
  plot_df_low$CTLA4[is.na(plot_df_low$CTLA4)]<-"CTLA4-"
  
  quantiles <- quantile(plot_df_high$CD274, probs = c(0,0.5,1))
  plot_df_high$CD274<-cut(plot_df_high$CD274, breaks = quantiles, labels = c("CD274-","CD274+"))
  plot_df_high$CD274[is.na(plot_df_high$CD274)]<-"CD274-"
  
  quantiles <- quantile(plot_df_mid$CD274, probs = c(0,0.5,1))
  plot_df_mid$CD274<-cut(plot_df_mid$CD274, breaks = quantiles, labels = c("CD274-","CD274+"))
  plot_df_mid$CD274[is.na(plot_df_mid$CD274)]<-"CD274-"
  
  quantiles <- quantile(plot_df_low$CD274, probs = c(0,0.5,1))
  plot_df_low$CD274<-cut(plot_df_low$CD274, breaks = quantiles, labels = c("CD274-","CD274+"))
  plot_df_low$CD274[is.na(plot_df_low$CD274)]<-"CD274-"
}

Liu_rna_info<-rbind(plot_df_high,plot_df_mid,plot_df_low)

all_composite <- summary(coxph(Surv(OS.time, OS_x) ~ Primary_Type_bin+Mstage..IIIC.0..M1a.1..M1b.2..M1c.3.+Gender+priorCTLA4+
                                 LAG3+CTLA4+CD274, data = plot_df_high))
all_mid <- summary(coxph(Surv(OS.time, OS_x) ~ Primary_Type_bin+Mstage..IIIC.0..M1a.1..M1b.2..M1c.3.+Gender+LAG3+priorCTLA4+
                           CTLA4+CD274, data = plot_df_mid))
all_somatic <- summary(coxph(Surv(OS.time, OS_x) ~ Primary_Type_bin+Mstage..IIIC.0..M1a.1..M1b.2..M1c.3.+Gender+priorCTLA4+
                               LAG3+CTLA4+CD274, data = plot_df_low))

big_cox<-rbind(all_composite$coefficients[(nrow(all_composite$coefficients)-2):nrow(all_composite$coefficients),],all_mid$coefficients[(nrow(all_composite$coefficients)-2):nrow(all_composite$coefficients),],all_somatic$coefficients[(nrow(all_composite$coefficients)-2):nrow(all_composite$coefficients),])
big_cox<-as.data.frame(big_cox[,c(1,3,5)])
big_cox$`lower .95`<-big_cox$coef-2*big_cox$`se(coef)`
big_cox$`upper .95`<-big_cox$coef+2*big_cox$`se(coef)`

big_cox$group<-factor(c(rep("MHC-II Reliant",3),rep("Balanced",3),rep("MHC-I Reliant",3)),levels = c("MHC-II Reliant","Balanced","MHC-I Reliant"))
big_cox$checkpoint<-c("LAG3","CTLA4","PD-L1","LAG3","CTLA4","PD-L1","LAG3","CTLA4","PD-L1")
big_cox$checkpoint_group<-factor(paste(big_cox$checkpoint,big_cox$group),levels=paste(big_cox$checkpoint,big_cox$group))

ggplot(data = big_cox, aes(x = checkpoint_group, y = `coef`, ymin = `lower .95`, ymax = `upper .95`, colour = group)) +
  geom_hline(yintercept=0, linetype="dashed",color = "black", size=1) +
  geom_errorbar(size=1.2,color="grey40",linetype=1,width=0.2) +
  geom_point(size=10,shape=15) + geom_label(label=paste("p =",round(big_cox$`Pr(>|z|)`,4)),nudge_x = 0.3) +
  coord_flip() + theme_bw(base_size = 22) + labs(y="Hazard Ratio (OS)",x="Model") + scale_color_manual(values=cols2[c(1,3,2)])+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()
        , panel.background = element_blank()
        , axis.line = element_line(colour = "black"),legend.position="none",
        ,axis.text = element_text(color="black"))  + ggtitle(label="Validation (N=77)")
#ggsave("/Users/tjsears/Code/GermTime/figs_july/fig6_v5/Liu_hazzy.pdf",width = 6.5,height = 9)








## ICI ##
plot_df<-ICI_phbr
cols2<-brewer.pal(3,"Dark2")

quantiles <- quantile(plot_df$Ratio, probs = c(0,0.33,0.66,1))
plot_df$group <- cut(plot_df$Ratio, breaks = quantiles, labels = c("PHBR2/PHBR1 Low","mid","PHBR2/PHBR1 High"))
plot_df$group<-factor(plot_df$group,levels = c("PHBR2/PHBR1 Low","PHBR2/PHBR1 High","mid"))
plot_df$group[is.na(plot_df$group)]<-"PHBR2/PHBR1 Low"
plot_df$LAG3[is.na(plot_df$LAG3)]<-0

plot_df$OS[plot_df$OS=="True"]<-1
plot_df$OS[plot_df$OS=="Y"]<-1
plot_df$OS[plot_df$OS=="Alive"]<-1
plot_df$OS[plot_df$OS=="Dead"]<-0
plot_df$OS[plot_df$OS=="False"]<-0
plot_df$OS[plot_df$OS=="Not available"]<-0
plot_df$OS[plot_df$OS=="Died not of Melanoma"]<-0
plot_df$OS<-as.numeric(plot_df$OS)
plot_df$Tissue<-ifelse(plot_df$study_cancer_x=="miao","RCC","melanoma")
plot_df[plot_df$study_cancer_x=="miao","OS"]=1-plot_df[plot_df$study_cancer_x=="miao","OS"]

quantiles <- quantile(plot_df$Ratio, probs = c(0,0.33,0.66,1))
plot_df$group <- cut(plot_df$Ratio, breaks = quantiles, labels = c("PHBR2/PHBR1 Low","mid","PHBR2/PHBR1 High"))
plot_df$group<-factor(plot_df$group,levels = c("PHBR2/PHBR1 Low","PHBR2/PHBR1 High","mid"))
plot_df$group[is.na(plot_df$group)]<-"PHBR2/PHBR1 Low"

plot_df_high<-plot_df[plot_df$group=="PHBR2/PHBR1 High",]
plot_df_mid<-plot_df[plot_df$group=="mid",]
plot_df_low<-plot_df[plot_df$group=="PHBR2/PHBR1 Low",]
plot_df_high<-plot_df_high[(plot_df_high$Escape+plot_df_high$numLostAlleles)>0,]

quantiles <- quantile(plot_df_high$LAG3, probs = c(0,0.5,1))
plot_df_high$Lag3cut<-cut(plot_df_high$LAG3, breaks = quantiles, labels = c("LAG3-","LAG3+"))
plot_df_high$Lag3cut[is.na(plot_df_high$Lag3cut)]<-"LAG3-"
plot_df_high$Lag3Bin<-ifelse(plot_df_high$Lag3cut=="LAG3+",2,1)
plot_df_high$acc<-ifelse(plot_df_high$Lag3Bin==plot_df_high$pheno,1,0)

quantiles <- quantile(plot_df_mid$LAG3, probs = c(0,0.5,1))
plot_df_mid$Lag3cut<-cut(plot_df_mid$LAG3, breaks = quantiles, labels = c("LAG3-","LAG3+"))
plot_df_mid$Lag3cut[is.na(plot_df_mid$Lag3cut)]<-"LAG3-"
plot_df_mid$Lag3Bin<-ifelse(plot_df_mid$Lag3cut=="LAG3+",2,1)
plot_df_mid$acc<-ifelse(plot_df_mid$Lag3Bin==plot_df_mid$pheno,1,0)

quantiles <- quantile(plot_df_low$LAG3, probs = c(0,0.5,1))
plot_df_low$Lag3cut<-cut(plot_df_low$LAG3, breaks = quantiles, labels = c("LAG3-","LAG3+"))
plot_df_low$Lag3cut[is.na(plot_df_low$Lag3cut)]<-"LAG3-"
plot_df_low$Lag3Bin<-ifelse(plot_df_low$Lag3cut=="LAG3+",2,1)
plot_df_low$acc<-ifelse(plot_df_low$Lag3Bin==plot_df_low$pheno,1,0)

km_trt_fit <- survfit(Surv(OS.time, OS) ~ Lag3cut, data=plot_df_high)
res<- pairwise_survdiff(Surv(OS.time, OS) ~ Lag3cut, data=plot_df_high,p.adjust.method = "none",rho=0)#, p.adjust.method = "none", rho = 0)
res$p.value<-round(res$p.value, digits = 4)
res

p1<-ggsurvplot(km_trt_fit,data=plot_df_high,size=2,censor.shape="|", censor.size = 8,legend=c(0.2,0.15),
               pval=F,xlim=c(0,1700),break.x.by=200,ncensor.plot.height=1,
               font.x=18,font.y=18,font.legend=18,font.title=18,
               legend.labs=c("LAG3- (17)","LAG3+ (16)"),palette=cols,legend.title="",
               ,xlab="Time in Days",ylim=c(0,1),ylab="Overall Survival (%)",surv.scale = c("percent"))+ggtitle(label="Discovery\nLAG3 +/- MHC-II Reliant")
p1$plot+ggplot2::annotate("label",x=1400,y = 0.8,size=6, label = paste("P =",round(res$p.value,4)))

#ggsave("/Users/tjsears/Code/GermTime/figs_july/fig6_v5/ICI_KM_MHCII_high.pdf",width = 6,height = 6)

km_trt_fit <- survfit(Surv(OS.time, OS) ~ Lag3cut, data=plot_df_mid)
res<- pairwise_survdiff(Surv(OS.time, OS) ~ Lag3cut, data=plot_df_mid,p.adjust.method = "none",rho=0)#, p.adjust.method = "none", rho = 0)
res$p.value<-round(res$p.value, digits = 4)
res

p1<-ggsurvplot(km_trt_fit,data=plot_df_mid,size=2,censor.shape="|", censor.size = 8,legend=c(0.2,0.15),
               pval=F,xlim=c(0,1700),break.x.by=200,ncensor.plot.height=1,
               font.x=18,font.y=18,font.legend=18,font.title=18,
               legend.labs=c("LAG3- (17)","LAG3+ (16)"),palette=cols,legend.title="",
               ,xlab="Time in Days",ylim=c(0,1),ylab="Overall Survival (%)",surv.scale = c("percent"))+ggtitle(label="Discovery\nLAG3 +/- Balanced Reliant")
p1$plot+ggplot2::annotate("label",x=1200,y = 0.9,size=6, label = paste("P =",round(res$p.value,4)))

#ggsave("/Users/tjsears/Code/GermTime/figs_july/fig6_v5/ICI_KM_MHCII_mid.pdf",width = 6,height = 6)


km_trt_fit <- survfit(Surv(OS.time, OS) ~ Lag3cut, data=plot_df_low)
res<- pairwise_survdiff(Surv(OS.time, OS) ~ Lag3cut, data=plot_df_low,p.adjust.method = "none",rho=0)#, p.adjust.method = "none", rho = 0)
res$p.value<-round(res$p.value, digits = 4)
res

p1<-ggsurvplot(km_trt_fit,data=plot_df_low,size=2,censor.shape="|", censor.size = 8,legend=c(0.2,0.15),
               pval=F,xlim=c(0,1700),break.x.by=200,ncensor.plot.height=1,
               font.x=18,font.y=18,font.legend=18,font.title=18,
               legend.labs=c("LAG3- (17)","LAG3+ (16)"),palette=cols,legend.title="",
               ,xlab="Time in Days",ylim=c(0,1),ylab="Overall Survival (%)",surv.scale = c("percent"))+ggtitle(label="Discovery\nLAG3 +/- MHC-I Reliant")
p1$plot+ggplot2::annotate("label",x=1200,y = 0.9,size=6, label = paste("P =",round(res$p.value,4)))

#ggsave("/Users/tjsears/Code/GermTime/figs_july/fig6_v5/ICI_KM_MHCII_low.pdf",width = 6,height = 6)


#hazzy plot ICI

plot_df<-ICI_phbr
cols<-brewer.pal(3,"Dark2")

plot_df$OS[plot_df$OS=="True"]<-1
plot_df$OS[plot_df$OS=="Y"]<-1
plot_df$OS[plot_df$OS=="Alive"]<-1
plot_df$OS[plot_df$OS=="Dead"]<-0
plot_df$OS[plot_df$OS=="False"]<-0
plot_df$OS[plot_df$OS=="Not available"]<-0
plot_df$OS[plot_df$OS=="Died not of Melanoma"]<-0
plot_df$OS<-as.numeric(plot_df$OS)
plot_df$Tissue<-ifelse(plot_df$study_cancer_x=="miao","RCC","melanoma")
plot_df[plot_df$study_cancer_x=="miao","OS"]=1-plot_df[plot_df$study_cancer_x=="miao","OS"]

plot_df$tissue<-1
plot_df$tissue[plot_df$study_cancer_x=="riaz"]<-2
plot_df$tissue[plot_df$study_cancer_x=="vanallen"]<-3
plot_df$tissue[plot_df$study_cancer_x=="hugo"]<-4

quantiles <- quantile(plot_df$Ratio, probs = c(0,0.33,0.66,1))
print(quantiles)

plot_df$group <- cut(plot_df$Ratio, breaks = quantiles, labels = c("PHBR2/PHBR1 Low","mid","PHBR2/PHBR1 High"))
plot_df$group<-factor(plot_df$group,levels = c("PHBR2/PHBR1 Low","PHBR2/PHBR1 High","mid"))
plot_df$group[is.na(plot_df$group)]<-"PHBR2/PHBR1 High"

plot_df_high<-plot_df[plot_df$group=="PHBR2/PHBR1 High",]
plot_df_mid<-plot_df[plot_df$group=="mid",]
plot_df_low<-plot_df[plot_df$group=="PHBR2/PHBR1 Low",]
plot_df_high<-plot_df_high[(plot_df_high$Escape+plot_df_high$numLostAlleles)>0,]

doMedians=T
if (doMedians){
  
  quantiles <- quantile(plot_df_high$LAG3, probs = c(0,0.5,1))
  plot_df_high$LAG3<-cut(plot_df_high$LAG3, breaks = quantiles, labels = c("LAG3-","LAG3+"))
  plot_df_high$LAG3[is.na(plot_df_high$LAG3)]<-"LAG3-"
  
  quantiles <- quantile(plot_df_mid$LAG3, probs = c(0,0.5,1))
  plot_df_mid$LAG3<-cut(plot_df_mid$LAG3, breaks = quantiles, labels = c("LAG3-","LAG3+"))
  plot_df_mid$LAG3[is.na(plot_df_mid$LAG3)]<-"LAG3-"
  
  quantiles <- quantile(plot_df_low$LAG3, probs = c(0,0.5,1))
  plot_df_low$LAG3<-cut(plot_df_low$LAG3, breaks = quantiles, labels = c("LAG3-","LAG3+"))
  plot_df_low$LAG3[is.na(plot_df_low$LAG3)]<-"LAG3-"
  
  quantiles <- quantile(plot_df_high$CTLA4, probs = c(0,0.5,1))
  plot_df_high$CTLA4<-cut(plot_df_high$CTLA4, breaks = quantiles, labels = c("CTLA4-","CTLA4+"))
  plot_df_high$CTLA4[is.na(plot_df_high$CTLA4)]<-"CTLA4-"
  
  quantiles <- quantile(plot_df_mid$CTLA4, probs = c(0,0.5,1))
  plot_df_mid$CTLA4<-cut(plot_df_mid$CTLA4, breaks = quantiles, labels = c("CTLA4-","CTLA4+"))
  plot_df_mid$CTLA4[is.na(plot_df_mid$CTLA4)]<-"CTLA4-"
  
  quantiles <- quantile(plot_df_low$CTLA4, probs = c(0,0.5,1))
  plot_df_low$CTLA4<-cut(plot_df_low$CTLA4, breaks = quantiles, labels = c("CTLA4-","CTLA4+"))
  plot_df_low$CTLA4[is.na(plot_df_low$CTLA4)]<-"CTLA4-"
  
  quantiles <- quantile(plot_df_high$CD274, probs = c(0,0.5,1))
  plot_df_high$CD274<-cut(plot_df_high$CD274, breaks = quantiles, labels = c("CD274-","CD274+"))
  plot_df_high$CD274[is.na(plot_df_high$CD274)]<-"CD274-"
  
  quantiles <- quantile(plot_df_mid$CD274, probs = c(0,0.5,1))
  plot_df_mid$CD274<-cut(plot_df_mid$CD274, breaks = quantiles, labels = c("CD274-","CD274+"))
  plot_df_mid$CD274[is.na(plot_df_mid$CD274)]<-"CD274-"
  
  quantiles <- quantile(plot_df_low$CD274, probs = c(0,0.5,1))
  plot_df_low$CD274<-cut(plot_df_low$CD274, breaks = quantiles, labels = c("CD274-","CD274+"))
  plot_df_low$CD274[is.na(plot_df_low$CD274)]<-"CD274-"
}

all_composite <- summary(coxph(Surv(OS.time, OS) ~ Age_y+Gender+tissue+LAG3+CTLA4+CD274, data = plot_df_high))
all_mid <- summary(coxph(Surv(OS.time, OS) ~ Age_y+Gender+tissue+LAG3+CTLA4+CD274, data = plot_df_mid))
all_somatic <- summary(coxph(Surv(OS.time, OS) ~ Age_y+Gender+tissue+LAG3+CTLA4+CD274, data = plot_df_low))

#cbind all data and export rna data
ICI_rna_info<-rbind(plot_df_high,plot_df_mid,plot_df_low)

big_cox<-rbind(all_composite$coefficients[(nrow(all_composite$coefficients)-2):nrow(all_composite$coefficients),],all_mid$coefficients[(nrow(all_composite$coefficients)-2):nrow(all_composite$coefficients),],all_somatic$coefficients[(nrow(all_composite$coefficients)-2):nrow(all_composite$coefficients),])
big_cox<-as.data.frame(big_cox[,c(1,3,5)])
big_cox$`lower .95`<-big_cox$coef-2*big_cox$`se(coef)`
big_cox$`upper .95`<-big_cox$coef+2*big_cox$`se(coef)`

big_cox$group<-factor(c(rep("MHC-II Reliant",3),rep("Balanced",3),rep("MHC-I Reliant",3)),levels = c("MHC-II Reliant","Balanced","MHC-I Reliant"))
big_cox$checkpoint<-c("LAG3","CTLA4","PD-L1","LAG3","CTLA4","PD-L1","LAG3","CTLA4","PD-L1")
big_cox$checkpoint_group<-factor(paste(big_cox$checkpoint,big_cox$group),levels=paste(big_cox$checkpoint,big_cox$group))

ggplot(data = big_cox, aes(x = checkpoint_group, y = `coef`, ymin = `lower .95`, ymax = `upper .95`, colour = group)) +
  geom_hline(yintercept=0, linetype="dashed",color = "black", size=1) +
  geom_errorbar(size=1.2,color="grey40",linetype=1,width=0.2) +
  geom_point(size=10,shape=15) + geom_label(label=paste("p =",round(big_cox$`Pr(>|z|)`,4)),nudge_x = 0.3) +
  coord_flip() + theme_bw(base_size = 22) + labs(y="Hazard Ratio (OS)",x="Model") + scale_color_manual(values=cols2[c(1,3,2)])+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()
        , panel.background = element_blank()
        , axis.line = element_line(colour = "black"),legend.position="none",
        ,axis.text = element_text(color="black"))  + ggtitle(label="Discovery (N=110)") 
#ggsave("/Users/tjsears/Code/GermTime/figs_july/fig6_v5/ICI_hazzy.pdf",width = 6.5,height = 9)


# MHC-II vs MHC-I km curves

tcga_phbr<-read.table("Data/Fig6/liu_composite_phbr_v4.txt",sep="\t",header=T)
ICI_phbr<-read.table("Data/Fig6/composite_phbr_v4.txt",sep="\t",header=T)

plot_df<-ICI_phbr
plot_df$OS[plot_df$OS=="True"]<-1
plot_df$OS[plot_df$OS=="Y"]<-1
plot_df$OS[plot_df$OS=="Alive"]<-0
plot_df$OS[plot_df$OS=="Dead"]<-1
plot_df$OS[plot_df$OS=="False"]<-0
plot_df$OS[plot_df$OS=="Not available"]<-0
plot_df$OS[plot_df$OS=="Died not of Melanoma"]<-0
plot_df$OS<-as.numeric(plot_df$OS)

plot_df<-plot_df[!is.na(plot_df$OS),]
plot_df<-plot_df[!is.na(plot_df$OS.time),]

quantiles <- quantile(plot_df$Ratio, probs = c(0,0.33,0.66,1))
plot_df$group <- cut(plot_df$Ratio, breaks = quantiles, labels = c("PHBR2/PHBR1 Low","mid","PHBR2/PHBR1 High"))
plot_df$group[is.na(plot_df$group)]<-"PHBR2/PHBR1 High"
table(plot_df$group)

table(plot_df$group)
plot_df<-plot_df[!plot_df$group=="mid",]

plot_df$Tissue<-ifelse(plot_df$study_cancer_x=="miao","RCC","melanoma")
plot_df[plot_df$study_cancer_x=="miao","OS"]=1-plot_df[plot_df$study_cancer_x=="miao","OS"]

plot_df<-plot_df[!(((plot_df$numLostAlleles+plot_df$Escape)==0)&plot_df$group=="PHBR2/PHBR1 High"),]

plot_df<-plot_df[plot_df$pheno==2,]
print(dim(plot_df))

table(plot_df$group)

km_trt_fit <- survfit(Surv(OS.time, OS) ~ group, data=plot_df)
res<- pairwise_survdiff(Surv(OS.time, OS) ~ group, data=plot_df,p.adjust.method = "none",rho=0)#, p.adjust.method = "none", rho = 0)
res$p.value<-round(res$p.value, digits = 4)
res

p1<-ggsurvplot(km_trt_fit,data=plot_df,size=2,censor.shape="|", censor.size = 8,
               font.x=18,font.y=18,font.legend=16,font.title=18,legend=c(0.3,0.15),
               pval=F,xlim=c(0,1700),break.x.by=200,ncensor.plot.height=1,legend.title="",
               legend.labs=c("MHC-I Reliant (45)","MHC-II Reliant (39)"),palette=rev(cols[1:2]),
               xlab="Time in Days",ylim=c(0,1),ylab="Overall Survival (%)",surv.scale = c("percent"))+ggtitle(label="Discovery Responders\nMHC-II Reliant vs MHC-I Reliant")
p1$plot+ggplot2::annotate("label",x=1200,y = 0.92,size=6, label = paste("P =",round(res$p.value,4)))

#ggsave("/Users/tjsears/Code/GermTime/figs_july/fig6_v5/MHCII_vsMHC_I.pdf",width = 6,height = 6)

plot_df$tissue<-1
plot_df$tissue[plot_df$study_cancer_x=="riaz"]<-2
plot_df$tissue[plot_df$study_cancer_x=="vanallen"]<-3
plot_df$tissue[plot_df$study_cancer_x=="hugo"]<-4

summary(coxph(Surv(OS.time, OS) ~ Age_y+Gender+tissue+group, data = plot_df))

# liu MHC-II vs I

plot_df<-tcga_phbr
plot_df$Ratio<-plot_df$PHBR2/(plot_df$PHBR1+1)

plot_df$OS<-as.numeric(plot_df$OS_x)

plot_df<-plot_df[!is.na(plot_df$OS),]
plot_df<-plot_df[!is.na(plot_df$OS.time),]

quantiles <- quantile(plot_df$Ratio, probs = c(0,0.33,0.66,1))
plot_df$group <- cut(plot_df$Ratio, breaks = quantiles, labels = c("PHBR2/PHBR1 Low","mid","PHBR2/PHBR1 High"))
plot_df$group[is.na(plot_df$group)]<-"PHBR2/PHBR1 High"
table(plot_df$group)

table(plot_df$group)
plot_df<-plot_df[!plot_df$group=="mid",]

plot_df<-plot_df[!(((plot_df$numLostAlleles+plot_df$Escape)==0)&plot_df$group=="PHBR2/PHBR1 High"),]

plot_df<-plot_df[plot_df$response_crist_sd==2,]
print(dim(plot_df))

table(plot_df$group)

km_trt_fit <- survfit(Surv(OS.time, OS) ~ group, data=plot_df)
res<- pairwise_survdiff(Surv(OS.time, OS) ~ group, data=plot_df,p.adjust.method = "none",rho=0)#, p.adjust.method = "none", rho = 0)
res$p.value<-round(res$p.value, digits = 4)
res

p1<-ggsurvplot(km_trt_fit,data=plot_df,size=2,censor.shape="|", censor.size = 8,
               font.x=18,font.y=18,font.legend=16,font.title=18,legend=c(0.3,0.15),
               pval=F,xlim=c(0,1700),break.x.by=200,ncensor.plot.height=1,legend.title="",
               legend.labs=c("MHC-I Reliant (15)","MHC-II Reliant (15)"),palette=rev(cols[1:2]),
               xlab="Time in Days",ylim=c(0,1),ylab="Overall Survival (%)",surv.scale = c("percent"))+ggtitle(label="Validation Responders\nMHC-II Reliant vs MHC-I Reliant")
p1$plot+ggplot2::annotate("label",x=1450,y = 0.92,size=6, label = paste("P =",round(res$p.value,4)))

#ggsave("/Users/tjsears/Code/GermTime/figs_july/fig6_v5/liu_MHCII_vsMHC_I.pdf",width = 6,height = 6)

plot_df$Primary_Type_bin<-ifelse(plot_df$Primary_Type=="skin",1,0)
summary(coxph(Surv(OS.time, OS_x) ~ Primary_Type_bin+Mstage..IIIC.0..M1a.1..M1b.2..M1c.3.+Gender+group, data = plot_df))


#export tables for Kairi
tcga_phbr_export<-tcga_phbr[,c("group","normal.WXS.id","tumor.WXS.id","tumor.RNA.id","Response")]
write.table(tcga_phbr_export,"/Users/tjsears/Code/GermTime/KairiMIXCR/Liu_mapping.txt",sep="\t",quote = F,row.names = F)

ICI_phbr_export<-ICI_phbr[,c("group","normal.WXS.id","tumor.WXS.id","tumor.RNA.id","Response")]
write.table(ICI_phbr_export,"/Users/tjsears/Code/GermTime/KairiMIXCR/Discovery_mapping.txt",sep="\t",quote = F,row.names = F)

#export rna cutoffs
ICI_rna_info_export<-ICI_rna_info[,c("tumor.WXS.id","LAG3","CD274","CTLA4","Response","group")]
Liu_rna_info_export<-Liu_rna_info[,c("tumor.WXS.id","LAG3","CD274","CTLA4","Response","group")]


#write.table(ICI_rna_info_export,"/Users/tjsears/Code/GermTime/FinalSuppTables/Discovery_RNA_cutoffs.txt",sep="\t",quote = F,row.names = F)
#write.table(Liu_rna_info_export,"/Users/tjsears/Code/GermTime/FinalSuppTables/Validation_RNA_cutoffs.txt",sep="\t",quote = F,row.names = F)







