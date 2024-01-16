#plot correlation, heatmaps, gene expression and infiltration estimates 
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

x<-read.table("Data/Fig4/Rplotting_RNAseq.txt",sep="\t",header=T)
x$response_crist_sd_factor<-ifelse(x$response_crist_sd==1,"NonResponder","Responder")
x$response_crist_sd_factor<-as.factor(x$response_crist_sd_factor)
x[is.na(x)]<-0

x$PFS[x$study_cancer_x=="miao"]<-1-x$PFS[x$study_cancer_x=="miao"]
x$OS[x$study_cancer_x=="miao"]<-1-x$OS[x$study_cancer_x=="miao"]

x$PRS_group<-ifelse(x$Composite_PRS>=5,"High Composite PRS","Low Composite PRS")
x = x %>%
  group_by(PRS_group) %>%
  mutate(checkpoint_target = scale(checkpoint_target))

x = x %>%
  group_by(study_cancer_x) %>%
  mutate(CD8.T.cells = scale(CD8.T.cells))
x = x %>%
  group_by(study_cancer_x) %>%
  mutate(T.cells.CD8 = scale(T.cells.CD8))
x = x %>%
  group_by(study_cancer_x) %>%
  mutate(CD8..Tem = scale(CD8..Tem))

x$Tcellmetrics<-x$T.cells.CD8

x = x %>%
  group_by(study_cancer_x) %>%
  mutate(B.lineage = scale(B.lineage))
x = x %>%
  group_by(study_cancer_x) %>%
  mutate(T.cells.follicular.helper = scale(T.cells.follicular.helper))
x = x %>%
  group_by(study_cancer_x) %>%
  mutate(Class.switched.memory.B.cells = scale(Class.switched.memory.B.cells))

# make TLS metric

TLS_sig_df<-read.table("Data/Fig4/TLS_sig.txt",sep="\t",header=T,row.names = 1)
df_<-colMeans(TLS_sig_df)

# merge x with TLS colmeans
df_<-as.data.frame(df_)
colnames(df_)<-c("TLS_sig")
x<-merge(x,df_,by.x="tumor.RNA.id",by.y=0)

#scale TLS metric of choice
x = x %>%
  group_by(study_cancer_x) %>%
  mutate(TLS_sig = scale(TLS_sig))
x$TLS_metrics<-x$TLS_sig

library(RColorBrewer)
cols<-brewer.pal(n = 10, name = "Paired")

#boxplots arrranged

upper_thresh=quantile(x$Composite_PRS,seq(0,1,0.5))[2]
lower_thresh=quantile(x$Composite_PRS,seq(0,1,0.5))[1]
upper_thresh=5
lower_thresh=5

x_boxplot<-x#[,c("study_cancer_x","response_crist_sd","Composite_PRS","checkpoint_target","checkpoint_off_target","PDCD1_x","CTLA4","CXCL13","CD8.T.cells")]
x_boxplot<-x_boxplot[x_boxplot$Composite_PRS>=upper_thresh,]

x_boxplot = x_boxplot %>%
  group_by(study_cancer_x) %>%
  mutate(CXCL13 = scale(CXCL13))
x_boxplot = x_boxplot %>%
  group_by(study_cancer_x) %>%
  mutate(CD8.T.cells = scale(CD8.T.cells))
x_boxplot = x_boxplot %>%
  group_by(study_cancer_x) %>%
  mutate(T.cells.CD8 = scale(T.cells.CD8))
x_boxplot = x_boxplot %>%
  group_by(study_cancer_x) %>%
  mutate(T.cells = scale(T.cells))
x_boxplot = x_boxplot %>%
  group_by(study_cancer_x) %>%
  mutate(T.cells.follicular.helper = scale(T.cells.follicular.helper))
x_boxplot = x_boxplot %>%
  group_by(study_cancer_x) %>%
  mutate(CD8..T.cells = scale(CD8..T.cells))
x_boxplot = x_boxplot %>%
  group_by(study_cancer_x) %>%
  mutate(CD8..Tem = scale(CD8..Tem))
x_boxplot = x_boxplot %>%
  group_by(study_cancer_x) %>%
  mutate(B.lineage = scale(B.lineage))
x_boxplot = x_boxplot %>%
  group_by(study_cancer_x) %>%
  mutate(Class.switched.memory.B.cells = scale(Class.switched.memory.B.cells))
x_boxplot_low = x_boxplot_low %>%
  group_by(study_cancer_x) %>%
  mutate(Cytotoxic.lymphocytes = scale(Cytotoxic.lymphocytes))

x_boxplot_long<-pivot_longer(x_boxplot,cols=c("checkpoint_target","T.cells.CD8","TLS_sig"))
x_boxplot_long$response<-ifelse(x_boxplot_long$response_crist_sd==1,"False Positive","True Positive")
x_boxplot_long$response<-factor(x_boxplot_long$response,levels=c("False Positive","True Positive"))

x_boxplot_long$name[x_boxplot_long$name=="CD8.T.cells"]<-"MCP  Tcells"
x_boxplot_long$name[x_boxplot_long$name=="T.cells.CD8"]<-"CIBERSORTx \n CD8 Tcells"
x_boxplot_long$name[x_boxplot_long$name=="CD8..Tem"]<-"xCell  Tem"

x_boxplot_long$name[x_boxplot_long$name=="checkpoint_target"]<-"CTLA4/PD1"

x_boxplot_long$name[x_boxplot_long$name=="TLS_sig"]<-"TLS Signature"

x_boxplot_long$system<-"Theraputic Target"
x_boxplot_long$system[grepl("CD8",x_boxplot_long$name)]<-"CD8+ T Cells"
x_boxplot_long$system[grepl("TLS Signature",x_boxplot_long$name)]<-"TLS Signature"

x_boxplot_long$system<-factor(x_boxplot_long$system,levels=c("Theraputic Target","CD8+ T Cells","TLS Signature"))

x_boxplot_low<-x
x_boxplot_low<-x_boxplot_low[x_boxplot_low$Composite_PRS<=lower_thresh,]

x_boxplot_low = x_boxplot_low %>%
  group_by(study_cancer_x) %>%
  mutate(CXCL13 = scale(CXCL13))
x_boxplot_low = x_boxplot_low %>%
  group_by(study_cancer_x) %>%
  mutate(CD8.T.cells = scale(CD8.T.cells))
x_boxplot_low = x_boxplot_low %>%
  group_by(study_cancer_x) %>%
  mutate(T.cells.CD8 = scale(T.cells.CD8))
x_boxplot_low = x_boxplot_low %>%
  group_by(study_cancer_x) %>%
  mutate(T.cells = scale(T.cells))
x_boxplot_low = x_boxplot_low %>%
  group_by(study_cancer_x) %>%
  mutate(T.cells.follicular.helper = scale(T.cells.follicular.helper))
x_boxplot_low = x_boxplot_low %>%
  group_by(study_cancer_x) %>%
  mutate(CD8..T.cells = scale(CD8..T.cells))
x_boxplot_low = x_boxplot_low %>%
  group_by(study_cancer_x) %>%
  mutate(CD8..Tem = scale(CD8..Tem))
x_boxplot_low = x_boxplot_low %>%
  group_by(study_cancer_x) %>%
  mutate(B.lineage = scale(B.lineage))
x_boxplot_low = x_boxplot_low %>%
  group_by(study_cancer_x) %>%
  mutate(Class.switched.memory.B.cells = scale(Class.switched.memory.B.cells))
x_boxplot_low = x_boxplot_low %>%
  group_by(study_cancer_x) %>%
  mutate(Cytotoxic.lymphocytes = scale(Cytotoxic.lymphocytes))

library(ggh4x)

cols=c("#B05D4A","#579F9F")

strip=strip_themed(background_x = elem_list_rect(fill = cols))

ggplot(x_boxplot_long, aes(x = response, y = value,fill=response)) + 
  xlab(label="Predicted Responders (IC-Index>=5)\nFalse Positive (n=10) / True Positive (n=19)") +
  facet_nested(~system)+ylab(label="Z-score")+theme_minimal(base_size=22)+
  theme(strip.background=element_rect(color="grey30", fill="grey90"))+
  geom_boxplot(position = position_dodge(width = -2),width=0.98)+
  scale_fill_manual(values=cols,name="Group")+ylim(-1.4,3.5)+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),panel.spacing = unit(0.1, "lines"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  stat_compare_means(comparisons = list(c("False Positive", "True Positive")),
                     method="t.test",method.args=list(alternative="less",var.equal=F),size=5.5,bracket.size = 0.55)
#ggsave("/Users/tjsears/Code/GermTime/figs_july/fig4_v4/Top_Tertile_TME_infiltrations.pdf",width = 10,height = 5)


# hazzy plots for each category

temp2<-x_boxplot
temp1<-x_boxplot_low

cols<-brewer.pal(n = 10, name = "Paired")
cols<-c("#666666","#A6761D")

#TLS metrics
all_composite <- summary(coxph(Surv(OS.time, OS) ~ Age_y+Gender+TLS_metrics, data = temp1))
all_somatic <- summary(coxph(Surv(OS.time, OS) ~ Age_y+Gender+TLS_metrics, data = temp2))

big_cox<-rbind(all_composite$coefficients[nrow(all_somatic$coefficients),],all_somatic$coefficients[nrow(all_somatic$coefficients),])
big_cox_ci<-rbind(all_composite$conf.int[nrow(all_somatic$coefficients),],all_somatic$conf.int[nrow(all_somatic$coefficients),])
big_cox<-cbind(big_cox[,5],big_cox_ci)
big_cox<-as.data.frame(big_cox)
big_cox$group<-factor(c("Predicted NonResponder","Predicted Responder"),levels = c("Predicted NonResponder","Predicted Responder"))

ggplot(data = big_cox, aes(x = group, y = `exp(coef)`, ymin = `lower .95`, ymax = `upper .95`, colour = group)) +
  geom_errorbar(size=1.2,color="grey40",linetype=1,width=0.2) + geom_hline(yintercept=1, linetype="dashed",color = "black", size=1) +
  geom_point(size=10,shape=15) + geom_label(label=paste("p =",round(big_cox$V1,4)),nudge_x = 0.25,nudge_y=0.14,color="Black",size=6) +
  coord_flip() + theme_bw(base_size = 22) + labs(y="Hazard Ratio (OS)",x="Model") + scale_color_manual(values=cols)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()
        , panel.background = element_blank()
        , axis.line = element_line(colour = "black"),  legend.position="none"
        ,axis.text = element_text(color="black"))  + ggtitle(label="TLS Signature")

#ggsave("/Users/tjsears/Code/GermTime/figs_july/fig4_v4/TLS.pdf",width = 8,height = 6)

#tcellmetrics
all_composite <- summary(coxph(Surv(OS.time, OS) ~ Age_y+Gender+Tcellmetrics, data = temp1))
all_somatic <- summary(coxph(Surv(OS.time, OS) ~ Age_y+Gender+Tcellmetrics, data = temp2))

big_cox<-rbind(all_composite$coefficients[nrow(all_somatic$coefficients),],all_somatic$coefficients[nrow(all_somatic$coefficients),])
big_cox_ci<-rbind(all_composite$conf.int[nrow(all_somatic$coefficients),],all_somatic$conf.int[nrow(all_somatic$coefficients),])
big_cox<-cbind(big_cox[,5],big_cox_ci)
big_cox<-as.data.frame(big_cox)
big_cox$group<-factor(c("Predicted NonResponder","Predicted Responder"),levels = c("Predicted NonResponder","Predicted Responder"))

ggplot(data = big_cox, aes(x = group, y = `exp(coef)`, ymin = `lower .95`, ymax = `upper .95`, colour = group)) +
  geom_errorbar(size=1.2,color="grey40",linetype=1,width=0.2) + geom_hline(yintercept=1, linetype="dashed",color = "black", size=1) +
  geom_point(size=10,shape=15) + geom_label(label=paste("p =",round(big_cox$V1,4)),nudge_x = 0.25,color="Black",size=6) +
  coord_flip() + theme_bw(base_size = 22) + labs(y="Hazard Ratio (OS)",x="Model") + scale_color_manual(values=cols)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()
        , panel.background = element_blank()
        , axis.line = element_line(colour = "black"),  legend.position="none",
        ,axis.text = element_text(color="black"))  + ggtitle(label="CD8+ T Cell Infiltrates") 
# geom_hline(yintercept=1, linetype="dashed",color = "black", size=1) #+ geom_text(aes(0,1,label = "", vjust = ,hjust=1.2),color="black",size=8)
#ggsave("/Users/tjsears/Code/GermTime/figs_july/fig4_v4/TcellMetrics.pdf",width = 8,height = 6)

#checkpoint
all_composite <- summary(coxph(Surv(OS.time, OS) ~ Age_y+Gender+checkpoint_target, data = temp1))
all_somatic <- summary(coxph(Surv(OS.time, OS) ~ Age_y+Gender+checkpoint_target, data = temp2))

big_cox<-rbind(all_composite$coefficients[nrow(all_somatic$coefficients),],all_somatic$coefficients[nrow(all_somatic$coefficients),])
big_cox_ci<-rbind(all_composite$conf.int[nrow(all_somatic$coefficients),],all_somatic$conf.int[nrow(all_somatic$coefficients),])
big_cox<-cbind(big_cox[,5],big_cox_ci)
big_cox<-as.data.frame(big_cox)
big_cox$group<-factor(c("Predicted \nNonResponder\n(n=43)  ","Predicted \nResponder\n(n=29)  "),levels = c("Predicted \nNonResponder\n(n=43)  ","Predicted \nResponder\n(n=29)  "))

ggplot(data = big_cox, aes(x = group, y = `exp(coef)`, ymin = `lower .95`, ymax = `upper .95`, colour = group)) +
  geom_errorbar(size=1.2,color="grey40",linetype=1,width=0.2) + geom_hline(yintercept=1, linetype="dashed",color = "black", size=1) +
  geom_point(size=10,shape=15) + geom_label(label=paste("p =",round(big_cox$V1,4)),nudge_x = 0.25,nudge_y = 0.12,color="Black",size=6) +
  coord_flip() + theme_bw(base_size = 22) + labs(y="Hazard Ratio (OS)",x="Model") + scale_color_manual(values=cols)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()
        , panel.background = element_blank()
        , axis.line = element_line(colour = "black"),  legend.position="none",
        ,axis.text = element_text(color="black"))  + ggtitle(label="Theraputic Target Expression") 
#geom_hline(yintercept=1, linetype="dashed",color = "black", size=1) #+ geom_text(aes(0,1,label = "", vjust = ,hjust=1.2),color="black",size=8)
#ggsave("/Users/tjsears/Code/GermTime/figs_july/fig4_v4/TargetCheckpoint.pdf",width = 8,height = 6)

#new auc plots
library(scales)

x_boxplot$TME_score<-x_boxplot$checkpoint_target+x_boxplot$TLS_metrics#+x_boxplot$Tcellmetrics
x_boxplot_low$TME_score<-x_boxplot_low$checkpoint_target+x_boxplot_low$TLS_metrics#+x_boxplot_low$Tcellmetrics
trueTME_median<-median(c(x_boxplot$TME_score,x_boxplot_low$TME_score))

#instead of auc lets try km curves
cols<-c("#A6761D","#EAC785")

x_boxplot$TME_median<-ifelse(x_boxplot$TME_score>=trueTME_median,"High TME Score","Low TME Score")

km_trt_fit <- survfit(Surv(OS.time, OS) ~ TME_median, data=x_boxplot)
res<- pairwise_survdiff(Surv(OS.time, OS) ~ TME_median, data=x_boxplot,p.adjust.method = "none",rho=0)#, p.adjust.method = "none", rho = 0)
res$p.value<-round(res$p.value, digits = 4)
res

p1<-ggsurvplot(km_trt_fit,data=x_boxplot,size=2,censor.shape="|", censor.size = 8,
               #ggtheme=(theme_minimal(base_size=18)),
               pval=F,xlim=c(0,1800),break.x.by=200,ncensor.plot.height=1,
               legend.labs=c("Hot TIME\n(n=15)","Cold TIME\n(n=14)"),palette=cols,ggtheme = theme_classic2(base_size=18),
               ,xlab="Time in Days",legend=c(0.23,0.18),ylim=c(0,1),ylab="Overall Survival (%)",surv.scale = c("percent"))+ggtitle(label="Predicted Responders")
p1$plot+ggplot2::annotate("label",x=1550,y = 0.8,size=6, label = paste("P =",round(res$p.value,4)))
#ggsave("/Users/tjsears/Code/GermTime/figs_july/fig4_v4/PredictedResponder.pdf",width = 6.5,height = 6)

cols<-c("#666666","#D6D6D6")

x_boxplot_low$TME_median<-ifelse(x_boxplot_low$TME_score>=trueTME_median,"High TME Score","Low TME Score")

km_trt_fit <- survfit(Surv(OS.time, OS) ~ TME_median, data=x_boxplot_low)
res<- pairwise_survdiff(Surv(OS.time, OS) ~ TME_median, data=x_boxplot_low,p.adjust.method = "none",rho=0)#, p.adjust.method = "none", rho = 0)
res$p.value<-round(res$p.value, digits = 4)
res

p2<-ggsurvplot(km_trt_fit,data=x_boxplot_low,size=2,censor.shape="|", censor.size = 8,
               #ggtheme=(theme_minimal(base_size=18)),
               pval=F,xlim=c(0,1800),break.x.by=200,ncensor.plot.height=1
               ,legend.labs=c("Hot TIME\n(n=21)","Cold TIME\n(n=21)"),palette=cols,ggtheme = theme_classic2(base_size=18),
               ,xlab="Time in Days",legend=c(0.23,0.18),ylim=c(0,1),ylab="Overall Survival (%)",surv.scale = c("percent"))+ggtitle(label="Predicted NonResponders")
p2$plot+ggplot2::annotate("label",x=1550,y = 0.8,size=6, label = paste("P =",round(res$p.value,3)))

#ggsave("/Users/tjsears/Code/GermTime/figs_july/fig4_v4/PredictedNonResponder.pdf",width = 6.5,height = 6)


### OS heatmap for 4 groups ###

# cristescu style heatmap showing EITHER response rates per stratification or OS hazard ratios

# take combined cohort (x), stratify by true TME median and IC-index 5

x$TME_score<-x$checkpoint_target+x$TLS_metrics
median<-median(x$TME_score)
x$TME_median<-ifelse(x$TME_score>=median,"TME High","TME Low")
x$IC_Index_median<-ifelse(x$Composite_PRS>=5,"IC-Index High","IC-Index Low")

table(x$TME_median,x$IC_Index_median)

BL<-x[(x$IC_Index_median=="IC-Index Low"&x$TME_median=="TME Low"),]
TL<-x[(x$IC_Index_median=="IC-Index High"&x$TME_median=="TME Low"),]
BR<-x[(x$IC_Index_median=="IC-Index Low"&x$TME_median=="TME High"),]
TR<-x[(x$IC_Index_median=="IC-Index High"&x$TME_median=="TME High"),]

BL<-table(BL$pheno)
TL<-table(TL$pheno)
BR<-table(BR$pheno)
TR<-table(TR$pheno)

bl<-BL[2]/(BL[1]+BL[2])
tl<-TL[2]/(TL[1]+TL[2])
br<-BR[2]/(BR[1]+BR[2])
tr<-TR[2]/(TR[1]+TR[2])

df_quad<-data.frame(c(tl,bl),c(tr,br))
print(df_quad)

df_quad<-as.data.frame(df_quad)
df_quad[1,1]<-as.character(df_quad)

col_fun = colorRamp2(c(0, 0.5, 1), c("#F2C14E", "whitesmoke", "#4D9078"))
col_fun = colorRamp2(c(0, 0.5, 1), c("#F2C14E", "whitesmoke", "#4D9078"))

library(ComplexHeatmap)
library(circlize)


fisher.test(data.frame(rbind(BL,TR)))
fisher.test(data.frame(rbind(TL,TR)))
fisher.test(data.frame(rbind(BL,BR)))
fisher.test(data.frame(rbind(TL,BR)))
fisher.test(data.frame(rbind(BR,TR)))

small_mat<-data.frame(c("7/14 (50%)","6/21 (29%)"),c("12/15 (80%)","5/21 (24%)"))

Heatmap(df_quad, name = "mat", rect_gp = gpar(col = "black", lwd = 2.25),show_column_names = F,
        show_heatmap_legend = T,col = rev(c("#F9BFC0","#D0ECF2","#F1F1F1","#CFD7E4"
                                            )),cluster_rows = T,cluster_columns=F,show_row_dend = F,show_column_dend = F,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%s", small_mat[i, j]), x, y, gp = gpar(fontsize = 20))
        })

library(ComplexHeatmap)
library(circlize)

#### SUPPLEMENTAL FIG S5 ####

if (FALSE){
  # lil high PRS vs low composite PRS  TIME boxplot?
  cols=c("#838E83","#5DA399")
  
  x_boxplot<-x#[,c("study_cancer_x","response_crist_sd","Composite_PRS","checkpoint_target","checkpoint_off_target","PDCD1_x","CTLA4","CXCL13","CD8.T.cells")]
  
  x_boxplot = x_boxplot %>%
    group_by(study_cancer_x) %>%
    mutate(CXCL13 = scale(CXCL13))
  x_boxplot = x_boxplot %>%
    group_by(study_cancer_x) %>%
    mutate(CD8.T.cells = scale(CD8.T.cells))
  x_boxplot = x_boxplot %>%
    group_by(study_cancer_x) %>%
    mutate(T.cells.CD8 = scale(T.cells.CD8))
  x_boxplot = x_boxplot %>%
    group_by(study_cancer_x) %>%
    mutate(T.cells = scale(T.cells))
  x_boxplot = x_boxplot %>%
    group_by(study_cancer_x) %>%
    mutate(T.cells.follicular.helper = scale(T.cells.follicular.helper))
  x_boxplot = x_boxplot %>%
    group_by(study_cancer_x) %>%
    mutate(CD8..T.cells = scale(CD8..T.cells))
  x_boxplot = x_boxplot %>%
    group_by(study_cancer_x) %>%
    mutate(CD8..Tem = scale(CD8..Tem))
  x_boxplot = x_boxplot %>%
    group_by(study_cancer_x) %>%
    mutate(B.lineage = scale(B.lineage))
  x_boxplot = x_boxplot %>%
    group_by(study_cancer_x) %>%
    mutate(Class.switched.memory.B.cells = scale(Class.switched.memory.B.cells))
  
  x_boxplot_long<-x_boxplot
  
  colnames(x_boxplot_long)[colnames(x_boxplot_long)=="CD8.T.cells"]<-"MCP CD8 Tcells"
  colnames(x_boxplot_long)[colnames(x_boxplot_long)=="T.cells.CD8"]<-"CIBERSORTx \n CD8 Tcells"
  colnames(x_boxplot_long)[colnames(x_boxplot_long)=="CD8..Tem"]<-"xCell CD8 Tem"
  
  colnames(x_boxplot_long)[colnames(x_boxplot_long)=="T.cells.follicular.helper"]<-"CIBERSORTx \n TFH cells"
  colnames(x_boxplot_long)[colnames(x_boxplot_long)=="B.lineage"]<-"MCP B cells"
  colnames(x_boxplot_long)[colnames(x_boxplot_long)=="Class.switched.memory.B.cells"]<-"xCell Class-switched \n memory B cells"
  
  colnames(x_boxplot_long)[colnames(x_boxplot_long)=="checkpoint_target"]<-"CTLA4/PD1"
  
  x_boxplot_temp<-x_boxplot_long
  x_boxplot_temp$TME_score<-x_boxplot_temp$`MCP CD8 Tcells`+x_boxplot_temp$`CIBERSORTx \n CD8 Tcells`+x_boxplot_temp$`xCell CD8 Tem`
  
  x_boxplot_temp$HiLowICIndex<-ifelse(x_boxplot_temp$Composite_PRS>=5,"High IC-Index","Low IC-Index")
  ggplot(x_boxplot_temp, aes(x = HiLowICIndex, y = TME_score,fill=HiLowICIndex)) + 
    ylab(label="T Cell Infiltrates")+theme_minimal(base_size=20)+
    theme(strip.background=element_rect(color="grey30", fill="grey90"))+
    geom_boxplot(position = position_dodge(width = -2),width=0.98)+
    scale_fill_manual(values=cols[2:1])+
    stat_compare_means(comparisons = list(c("High IC-Index","Low IC-Index")),
                       method="t.test",size=6,bracket.size = 0.65)+
    theme(legend.position="none")+ggtitle("")
  
  #ggsave("/Users/tjsears/Code/GermTime/figs_july/figS4/TIME_TcellInfiltrate.pdf",width = 6,height =6)
  
  x_boxplot_temp$TME_score<-x_boxplot_temp$TLS_sig
  
  x_boxplot_temp$HiLowICIndex<-ifelse(x_boxplot_temp$Composite_PRS>=5,"High IC-Index","Low IC-Index")
  ggplot(x_boxplot_temp, aes(x = HiLowICIndex, y = TME_score,fill=HiLowICIndex)) + 
    ylab(label="TLS Infiltrates")+theme_minimal(base_size=20)+
    theme(strip.background=element_rect(color="grey30", fill="grey90"))+
    geom_boxplot(position = position_dodge(width = -2),width=0.98)+
    scale_fill_manual(values=cols[2:1])+ylim(-1,5)+
    stat_compare_means(comparisons = list(c("High IC-Index","Low IC-Index")),
                       method="t.test",size=6,bracket.size = 0.65)+
    theme(legend.position="none")+ggtitle("")
  
  #ggsave("/Users/tjsears/Code/GermTime/figs_july/figS4/TIME_TLS.pdf",width = 6,height =6)
  
  
  x_boxplot_temp$TME_score<-x_boxplot_temp$`CTLA4/PD1`
  
  x_boxplot_temp$HiLowICIndex<-ifelse(x_boxplot_temp$Composite_PRS>=5,"High IC-Index","Low IC-Index")
  ggplot(x_boxplot_temp, aes(x = HiLowICIndex, y = TME_score,fill=HiLowICIndex)) + 
    ylab(label="Target Checkpoint")+theme_minimal(base_size=20)+
    theme(strip.background=element_rect(color="grey30", fill="grey90"))+
    geom_boxplot(position = position_dodge(width = -2),width=0.98)+
    scale_fill_manual(values=cols[2:1])+
    stat_compare_means(comparisons = list(c("High IC-Index","Low IC-Index")),
                       method="t.test",size=6,bracket.size = 0.65)+
    theme(legend.position="none")+ggtitle("")
  
  #ggsave("/Users/tjsears/Code/GermTime/figs_july/figS4/TIME_CPI.pdf",width = 6,height =6)
  
}

#scatterplot of TME and IC-Index correlation if needed

# orthogonality between TME and IC-Index?
cor.test(x$TME_score,x$Composite_PRS)



# Heatmaps for Somatic and Germline PRS vs TIME infiltration :)


#SOMATIC 

x$TME_score<-x$checkpoint_target+x$TLS_metrics+x$Tcellmetrics
median<-median(x$TME_score)
x$TME_median<-ifelse(x$TME_score>=median,"TME High","TME Low")
x$IC_Index_median<-ifelse(x$Somatic_PRS>=5,"IC-Index High","IC-Index Low")

table(x$TME_median,x$IC_Index_median)

BL<-x[(x$IC_Index_median=="IC-Index Low"&x$TME_median=="TME Low"),]
TL<-x[(x$IC_Index_median=="IC-Index High"&x$TME_median=="TME Low"),]
BR<-x[(x$IC_Index_median=="IC-Index Low"&x$TME_median=="TME High"),]
TR<-x[(x$IC_Index_median=="IC-Index High"&x$TME_median=="TME High"),]

BL<-table(BL$pheno)
TL<-table(TL$pheno)
BR<-table(BR$pheno)
TR<-table(TR$pheno)

bl<-BL[2]/(BL[1]+BL[2])
tl<-TL[2]/(TL[1]+TL[2])
br<-BR[2]/(BR[1]+BR[2])
tr<-TR[2]/(TR[1]+TR[2])

df_quad<-data.frame(c(tl,bl),c(tr,br))
print(df_quad)

df_quad<-as.data.frame(df_quad)
df_quad[1,1]<-as.character(df_quad)

col_fun = colorRamp2(c(0, 0.5, 1), c("#F2C14E", "whitesmoke", "#4D9078"))
col_fun = colorRamp2(c(0, 0.5, 1), c("#F2C14E", "whitesmoke", "#4D9078"))

library(ComplexHeatmap)
library(circlize)

#Heatmap(df_quad,cluster_rows = F,cluster_columns=F,show_row_dend = F,show_column_dend = F,col=col_fun)

fisher.test(data.frame(rbind(BL,TR)))
fisher.test(data.frame(rbind(TL,TR)))
fisher.test(data.frame(rbind(BL,BR)))
fisher.test(data.frame(rbind(TL,BR)))
fisher.test(data.frame(rbind(BR,TR)))

small_mat<-data.frame(c("9/15 (60%)","4/20 (20%)"),c("9/14 (64%)","8/22 (36%)"))

Heatmap(df_quad, name = "mat", rect_gp = gpar(col = "black", lwd = 1.5),show_column_names = F,
        show_heatmap_legend = T,col = rev(c("#F9BFC0","#D0ECF2","#CFD7E4",
                                            "#F1F1F1")),cluster_rows = F,cluster_columns=F,show_row_dend = F,show_column_dend = F,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%s", small_mat[i, j]), x, y, gp = gpar(fontsize = 20))
        })



# GERMLINE

x$TME_score<-x$checkpoint_target+x$TLS_metrics+x$Tcellmetrics
median<-median(x$TME_score)
x$TME_median<-ifelse(x$TME_score>=median,"TME High","TME Low")
x$IC_Index_median<-ifelse(x$Germline_PRS>=5,"IC-Index High","IC-Index Low")

table(x$TME_median,x$IC_Index_median)

BL<-x[(x$IC_Index_median=="IC-Index Low"&x$TME_median=="TME Low"),]
TL<-x[(x$IC_Index_median=="IC-Index High"&x$TME_median=="TME Low"),]
BR<-x[(x$IC_Index_median=="IC-Index Low"&x$TME_median=="TME High"),]
TR<-x[(x$IC_Index_median=="IC-Index High"&x$TME_median=="TME High"),]

BL<-table(BL$pheno)
TL<-table(TL$pheno)
BR<-table(BR$pheno)
TR<-table(TR$pheno)

bl<-BL[2]/(BL[1]+BL[2])
tl<-TL[2]/(TL[1]+TL[2])
br<-BR[2]/(BR[1]+BR[2])
tr<-TR[2]/(TR[1]+TR[2])

df_quad<-data.frame(c(tl,bl),c(tr,br))
print(df_quad)

df_quad<-as.data.frame(df_quad)
df_quad[1,1]<-as.character(df_quad)

col_fun = colorRamp2(c(0, 0.5, 1), c("#F2C14E", "whitesmoke", "#4D9078"))
col_fun = colorRamp2(c(0, 0.5, 1), c("#F2C14E", "whitesmoke", "#4D9078"))

library(ComplexHeatmap)
library(circlize)

#Heatmap(df_quad,cluster_rows = F,cluster_columns=F,show_row_dend = F,show_column_dend = F,col=col_fun)

fisher.test(data.frame(rbind(BL,TR)))
fisher.test(data.frame(rbind(TL,TR)))
fisher.test(data.frame(rbind(BL,BR)))
fisher.test(data.frame(rbind(TL,BR)))
fisher.test(data.frame(rbind(BR,TR)))

small_mat<-data.frame(c("7/19 (37%)","6/16 (38%)"),c("11/17 (65%)","6/19 (32%)"))

Heatmap(df_quad, name = "mat", rect_gp = gpar(col = "black", lwd = 1.5),show_column_names = F,
        show_heatmap_legend = T,col = rev(c("#F9BFC0","#F1F1F1","#D0ECF2","#CFD7E4"
        )),cluster_rows = F,cluster_columns=F,show_row_dend = F,show_column_dend = F,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%s", small_mat[i, j]), x, y, gp = gpar(fontsize = 20))
        })





