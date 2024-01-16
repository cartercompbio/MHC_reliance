# Model performance figure

library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(ggh4x)
library(RColorBrewer)
cols<-brewer.pal(n = 3, name = "Dark2")

ici<-read.table("Data/Fig3/all_prs.txt",sep="\t",header=T)
data_long <- ici %>%
  pivot_longer(cols=c("Germline_PRS", "Somatic_PRS", "Composite_PRS"))

data_long$name<-gsub("_PRS","",data_long$name)
data_long$name<-factor(data_long$name,levels=c("Somatic","Germline","Composite"))

data_long$Response<-ifelse(data_long$response_crist_sd==1,"NonResponder","Responder")
data_long$Response<-as.factor(data_long$Response)

my_comparisons <- list(c("NonResponder","Responder"))

data_long_Vanallen<-data_long[data_long$study_cancer_x=="vanallen",]
data_long_Miao<-data_long[data_long$study_cancer_x=="miao",]
data_long_Rizvi<-data_long[data_long$study_cancer_x=="rizvi",]

model_cols=c("#BCD3F2", "#7EA0B7","#B9CBB5")
ridiculous_strips <- strip_themed(
  # Horizontal strips
  background_x = elem_list_rect(fill = c("#BCD3F2", "#7EA0B7","#B9CBB5")), #blues
  #background_x = elem_list_rect(fill = c("#FFB997", "#F67E7D","#843B62")),
  
  text_x = elem_list_text(colour = c("black", "black","black"),
                          face = c("bold", "bold","bold")),
  by_layer_x = F,
)

ggplot(data_long_Vanallen, aes(x = Response, y = value,fill=Response)) + 
  facet_nested(~name,strip = ridiculous_strips,scales="free_x")+ylab(label="IC-Index")+theme_bw(base_size=20)+
  theme(strip.background=element_rect(color="grey30", fill="grey90"))+
  geom_boxplot(position = position_dodge(width = -2),width=0.98)+
  scale_fill_manual(values=cols[2:1])+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),panel.spacing = unit(0.1, "lines"),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  stat_compare_means(comparisons = list(c("NonResponder", "Responder")),
                     method="t.test",size=6,bracket.size = 0.65)+ylim(c(0,10.5))+
  theme(legend.position="none")+ggtitle("Vanallen (NR=81,R=29)")
#ggsave("/Users/tjsears/Code/GermTime/figs_july/fig3/vanallen_box.pdf",width = 6,height = 6)

ggplot(data_long_Miao, aes(x = Response, y = value,fill=Response)) + 
  facet_nested(~name,strip = ridiculous_strips,scales="free_x")+ylab(label="IC-Index")+theme_bw(base_size=20)+
  theme(strip.background=element_rect(color="grey30", fill="grey90"))+
  geom_boxplot(position = position_dodge(width = -2),width=0.98)+
  scale_fill_manual(values=cols[2:1])+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),panel.spacing = unit(0.1, "lines"),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  stat_compare_means(comparisons = list(c("NonResponder", "Responder")),
                     method="t.test",size=6,bracket.size = 0.65)+
  theme(legend.position="none")+ggtitle("Miao (NR=27,R=43)")+ylim(c(0,10.5))
#ggsave("/Users/tjsears/Code/GermTime/figs_july/fig3/miao_box.pdf",width = 6,height =6)

ggplot(data_long_Rizvi, aes(x = Response, y = value,fill=Response)) + 
  facet_nested(~name,strip = ridiculous_strips,scales="free_x")+ylab(label="IC-Index")+theme_bw(base_size=20)+
  theme(strip.background=element_rect(color="grey30", fill="grey90"))+
  geom_boxplot(position = position_dodge(width = -2),width=0.98)+
  scale_fill_manual(values=cols[2:1])+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),panel.spacing = unit(0.1, "lines"),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  stat_compare_means(comparisons = list(c("NonResponder", "Responder")),
                     method="t.test",size=6,bracket.size = 0.65)+
  theme(legend.position="none")+ggtitle("Rizvi (NR=12,R=22)")+ylim(c(0,10.5))
#ggsave("/Users/tjsears/Code/GermTime/figs_july/fig3/rizvi_box.pdf",width = 6,height = 6)

# barplot of cliffs D

# strip down to cohort specific data
ici_vanallen<-ici[ici$study_cancer_x=="vanallen",]
ici_vanallen_R<-ici_vanallen[ici_vanallen$pheno==2,]
ici_vanallen_NR<-ici_vanallen[ici_vanallen$pheno==1,]

ici_miao<-ici[ici$study_cancer_x=="miao",]
ici_miao_R<-ici_miao[ici_miao$pheno==2,]
ici_miao_NR<-ici_miao[ici_miao$pheno==1,]

ici_rizvi<-ici[ici$study_cancer_x=="rizvi",]
ici_rizvi_R<-ici_rizvi[ici_rizvi$pheno==2,]
ici_rizvi_NR<-ici_rizvi[ici_rizvi$pheno==1,]

# get cliff's D for each of these
library(effsize)

vg<-cliff.delta(ici_vanallen_R$Germline_PRS, ici_vanallen_NR$Germline_PRS)
vs<-cliff.delta(ici_vanallen_R$Somatic_PRS, ici_vanallen_NR$Somatic_PRS)
vc<-cliff.delta(ici_vanallen_R$Composite_PRS, ici_vanallen_NR$Composite_PRS)

mg<-cliff.delta(ici_miao_R$Germline_PRS, ici_miao_NR$Germline_PRS)
ms<-cliff.delta(ici_miao_R$Somatic_PRS, ici_miao_NR$Somatic_PRS)
mc<-cliff.delta(ici_miao_R$Composite_PRS, ici_miao_NR$Composite_PRS)

rg<-cliff.delta(ici_rizvi_R$Germline_PRS, ici_rizvi_NR$Germline_PRS)
rs<-cliff.delta(ici_rizvi_R$Somatic_PRS, ici_rizvi_NR$Somatic_PRS)
rc<-cliff.delta(ici_rizvi_R$Composite_PRS, ici_rizvi_NR$Composite_PRS)

# barplot of each delta

cliffplot<-data.frame(delta=c(vg$estimate,vs$estimate,vc$estimate,mg$estimate,ms$estimate,
                              mc$estimate,rg$estimate,rs$estimate,rc$estimate),
                      size=c(vg$magnitude,vs$magnitude,vc$magnitude,mg$magnitude,ms$magnitude,
                             mc$magnitude,rg$magnitude,rs$magnitude,rc$magnitude),
                      group=factor(c(rep("Vanallen",3),rep("Miao",3),rep("Rizvi",3)),levels=c("Vanallen","Rizvi","Miao")),
                      model=factor(rep(c("Germline","Somatic","Composite"),3),levels=c("Somatic","Germline","Composite"))
                      )

ggplot(cliffplot, aes(x = model, y = delta,fill=size)) + 
  facet_nested(~group,scales="free_x")+ylab(label="Cliff's Effect Size")+theme_bw(base_size=20)+
  theme(strip.background=element_rect(color="grey30", fill="grey90"))+
  geom_col(width=0.9,color="black")+  labs(fill = "Effect Size")+
  #scale_fill_manual(values=cols[2:1])+
  theme(panel.spacing = unit(0.35, "lines"),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+xlab(label=NULL)+
  theme(axis.text.x=element_text(angle=45,vjust=1.1,hjust = 1))+ggtitle("Cliff's Delta NR vs R")+ylim(c(0,0.75))

#ggsave("/Users/tjsears/Code/GermTime/figs_july/fig3/CliffsBarPlotV1.pdf",width = 6,height = 5)

ggplot(cliffplot, aes(x = model, y = delta,fill=model)) + 
  facet_nested(~group,scales="free_x")+ylab(label="Cliff's Effect Size")+theme_bw(base_size=20)+
  theme(strip.background=element_rect(color="grey30", fill="grey90"))+
  geom_col(width=0.9,color="black")+  labs(fill = "Model")+
  scale_fill_manual(values=model_cols)+
  xlab(label=NULL)+
  theme(panel.spacing = unit(0.35, "lines"),axis.text.x=element_blank(),axis.ticks.x=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ggtitle("Cliff's Delta NR vs R")+ylim(c(0,0.75))

#ggsave("/Users/tjsears/Code/GermTime/figs_july/fig3/CliffsBarPlotV2.pdf",width = 6,height = 5)

#vanallen
ici_vanallen<-ici[ici$study_cancer_x=="vanallen",]

rocobj_vanallen_somatic <- pROC::roc(ici_vanallen$pheno,ici_vanallen$Somatic_PRS)
auc_somatic <- round(pROC::auc(ici_vanallen$pheno,ici_vanallen$Somatic_PRS),4)

rocobj_vanallen_germline <- pROC::roc(ici_vanallen$pheno,ici_vanallen$Germline_PRS)
auc_germline <- round(pROC::auc(ici_vanallen$pheno,ici_vanallen$Germline_PRS),4)

rocobj_vanallen_composite <- pROC::roc(ici_vanallen$pheno,ici_vanallen$Composite_PRS)
auc_composite <- round(pROC::auc(ici_vanallen$pheno,ici_vanallen$Composite_PRS),4)

library(pROC)
ggroc(list(call_roc_name_1 = rocobj_vanallen_somatic, call_roc_name_2 = rocobj_vanallen_germline,rocobj_vanallen_composite),size=1.2,aes = c("color"))+
  theme_bw(base_size=20)+
  theme(legend.position = c(0.7, 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_color_discrete(labels=c(paste("Somatic IC-Index \n AUC=",auc_somatic,sep=""),
                                                                         paste("Germline IC-Index \n AUC=",auc_germline,sep=""),
                                                                          paste("Composite IC-Index \n AUC=",auc_composite,sep="")),type = model_cols) +
  labs(color=NULL) #+ ggtitle("Vanallen Classification")
#ggsave("/Users/tjsears/Code/GermTime/figs_july/fig3/AUC_vanallen.pdf",width = 6,height = 6)


#miao
ici_miao<-ici[ici$study_cancer_x=="miao",]

rocobj_miao_somatic <- pROC::roc(ici_miao$pheno,ici_miao$Somatic_PRS)
auc_somatic <- round(pROC::auc(ici_miao$pheno,ici_miao$Somatic_PRS),4)

rocobj_miao_germline <- pROC::roc(ici_miao$pheno,ici_miao$Germline_PRS)
auc_germline <- round(pROC::auc(ici_miao$pheno,ici_miao$Germline_PRS),4)

rocobj_miao_composite <- pROC::roc(ici_miao$pheno,ici_miao$Composite_PRS)
auc_composite <- round(pROC::auc(ici_miao$pheno,ici_miao$Composite_PRS),4)

library(pROC)
ggroc(list(call_roc_name_1 = rocobj_miao_somatic, call_roc_name_2 = rocobj_miao_germline,rocobj_miao_composite),size=1.2,aes = c("color"))+
  theme_bw(base_size=20)+
  theme(legend.position = c(0.7, 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_color_discrete(labels=c(paste("Somatic IC-Index \n AUC=",auc_somatic,sep=""),
                                                                       paste("Germline IC-Index \n AUC=",auc_germline,sep=""),
                                                                       paste("Composite IC-Index \n AUC=",auc_composite,sep="")),type = model_cols) +
  labs(color=NULL) #+ ggtitle("miao Classification")
#ggsave("/Users/tjsears/Code/GermTime/figs_july/fig3/AUC_miao.pdf",width = 6,height = 6)


#rizvi
ici_rizvi<-ici[ici$study_cancer_x=="rizvi",]

rocobj_rizvi_somatic <- pROC::roc(ici_rizvi$pheno,ici_rizvi$Somatic_PRS)
auc_somatic <- round(pROC::auc(ici_rizvi$pheno,ici_rizvi$Somatic_PRS),4)

rocobj_rizvi_germline <- pROC::roc(ici_rizvi$pheno,ici_rizvi$Germline_PRS)
auc_germline <- round(pROC::auc(ici_rizvi$pheno,ici_rizvi$Germline_PRS),4)

rocobj_rizvi_composite <- pROC::roc(ici_rizvi$pheno,ici_rizvi$Composite_PRS)
auc_composite <- round(pROC::auc(ici_rizvi$pheno,ici_rizvi$Composite_PRS),4)

library(pROC)
ggroc(list(call_roc_name_1 = rocobj_rizvi_somatic, call_roc_name_2 = rocobj_rizvi_germline,rocobj_rizvi_composite),size=1.2,aes = c("color"))+
  theme_bw(base_size=20)+
  theme(legend.position = c(0.7, 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_color_discrete(labels=c(paste("Somatic IC-Index \n AUC=",auc_somatic,sep=""),
                                                                       paste("Germline IC-Index \n AUC=",auc_germline,sep=""),
                                                                       paste("Composite IC-Index \n AUC=",auc_composite,sep="")),type = model_cols) +
  labs(color=NULL) #+ ggtitle("rizvi Classification")
#ggsave("/Users/tjsears/Code/GermTime/figs_july/fig3/AUC_rizvi.pdf",width = 6,height = 6)


# plot final graphs
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
library(gt)

ici_km<-ici

ici_km$quartile <- ntile(ici_km$Composite_PRS, 3)
ici_km$PFS[ici_km$study_cancer_x=="miao"]<-1-ici_km$PFS[ici_km$study_cancer_x=="miao"]

km_trt_fit <- survfit(Surv(PFS.time, PFS) ~ quartile, data=ici_km)
res<- pairwise_survdiff(Surv(PFS.time, PFS) ~ quartile, data=ici_km,p.adjust.method = "none",rho=0)#, p.adjust.method = "none", rho = 0)
res<-as.data.frame(res$p.value)

p<-ggsurvplot(km_trt_fit,data=ici_km,size=2,censor.shape="|", censor.size = 8,
              legend=c(0.8,0.8),
              ,pval=F,xlim=c(0,800),break.x.by=200,ncensor.plot.height=1,title="Composite IC-Index Tertiles"
              ,legend.labs=c("Lowest (n=72)","Middle (n=71)","Highest (n=71)"),palette=c(rev(cols[c(1,3,2)]))
              ,xlab="Time in Days",ylim=c(0,1),ylab="Progression-Free Survival (%)",surv.scale = c("percent"))
p#+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#ggsave("/Users/tjsears/Code/GermTime/figs_july/fig3/km.pdf",width = 5,height = 5)

colnames(res)<-c("Lowest","Middle")
rownames(res)<-c("Middle","Highest")
res<-round(res,6)

my_table <- gt(res,rownames_to_stub = T) %>%
  tab_style(
    style = cell_text(size = 12),
    locations = cells_body()
  )
my_table

#ggsave("/Users/tjsears/Code/GermTime/figs_july/fig3/table.pdf",width = 4,height = 4,plot=as.grob(my_table))


# hazzy plot
ici_km<-ici
ici_km$PFS[ici_km$study_cancer_x=="miao"]<-1-ici_km$PFS[ici_km$study_cancer_x=="miao"]

vanallen_composite <- summary(coxph(Surv(PFS.time, PFS) ~ Gender+Age_y+Composite_PRS, data = ici_vanallen))
vanallen_somatic <- summary(coxph(Surv(PFS.time, PFS) ~ Gender+Age_y+Somatic_PRS, data = ici_vanallen))
vanallen_germline <- summary(coxph(Surv(PFS.time, PFS) ~ Gender+Age_y+Germline_PRS, data = ici_vanallen))

miao_composite <- summary(coxph(Surv(PFS.time, PFS) ~ Gender+Age_y+Composite_PRS, data = ici_miao))
miao_somatic <- summary(coxph(Surv(PFS.time, PFS) ~ Gender+Age_y+Somatic_PRS, data = ici_miao))
miao_germline <- summary(coxph(Surv(PFS.time, PFS) ~ Gender+Age_y+Germline_PRS, data = ici_miao))

rizvi_composite <- summary(coxph(Surv(PFS.time, PFS) ~ Gender+Age_y+Composite_PRS, data = ici_rizvi))
rizvi_somatic <- summary(coxph(Surv(PFS.time, PFS) ~ Gender+Age_y+Somatic_PRS, data = ici_rizvi))
rizvi_germline <- summary(coxph(Surv(PFS.time, PFS) ~ Gender+Age_y+Germline_PRS, data = ici_rizvi))


#all samples

all_composite <- summary(coxph(Surv(PFS.time, PFS) ~ Gender+Age_y+study_cancer_x+Composite_PRS, data = ici_km))
all_somatic <- summary(coxph(Surv(PFS.time, PFS) ~ Gender+Age_y+study_cancer_x+Somatic_PRS, data = ici_km))
all_germline <- summary(coxph(Surv(PFS.time, PFS) ~ Gender+Age_y+study_cancer_x+Germline_PRS, data = ici_km))

c_cox<-coxph(Surv(PFS.time, PFS) ~ Gender+Age_y+Composite_PRS, data = ici_km,x=T)
s_cox<-coxph(Surv(PFS.time, PFS) ~ Gender+Age_y+Somatic_PRS, data = ici_km,x=T)
g_cox<-coxph(Surv(PFS.time, PFS) ~ Gender+Age_y+Germline_PRS, data = ici_km,x=T)

big_cox<-rbind(all_composite$coefficients[5,],all_somatic$coefficients[5,],all_germline$coefficients[5,])
big_cox_ci<-rbind(all_composite$conf.int[5,],all_somatic$conf.int[5,],all_germline$conf.int[5,])
big_cox<-cbind(big_cox[,5],big_cox_ci)
big_cox<-as.data.frame(big_cox)
big_cox$group<-factor(c("Composite","Somatic","Germline"),levels = c("Composite","Somatic","Germline"))

ggplot(data = big_cox, aes(x = group, y = `exp(coef)`, ymin = `lower .95`, ymax = `upper .95`, colour = group)) +
  geom_errorbar(size=1.2,color="grey40",linetype=1,width=0.2) +
  geom_point(size=12,shape=15) + geom_label(label=paste("p =",round(big_cox$V1,5)),nudge_x = 0.25,color="black") +
  coord_flip() + theme_bw(base_size = 22) + labs(y="Hazard Ratio",x="Model") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()
        , panel.background = element_blank()
        , axis.line = element_line(colour = "black"),  legend.position="none",
        ,axis.text = element_text(color="black")) + scale_color_manual(values = (model_cols[c(3,1,2)])) +
  geom_hline(yintercept=1, linetype="dashed",color = "black", size=1) #+ geom_text(aes(0,1,label = "", vjust = ,hjust=1.2),color="black",size=8)

#ggsave("/Users/tjsears/Code/GermTime/figs_july/fig3/hazzy.pdf",width = 6,height = 6)


#intra model comparison
library(nonnestcox)

plrtest(c_cox, s_cox, nested=F)

# which patients are being reclassified? are we boosting PPV or NPV?

library(cutpointr)
library(pROC)
library(DTComPair)

cutoff1=4.92654 #generated from cutpointr
cutoff2=5.49969
cutoff3=6.23346

cutoff_x<-mean(c(cutoff1,cutoff2,cutoff3))

ici_temp<-ici
ici_temp$somatic_cutoff<-ifelse(ici_temp$Somatic_PRS>=cutoff1,1,0)
ici_temp$germline_cutoff<-ifelse(ici_temp$Germline_PRS>=cutoff2,1,0)
ici_temp$composite_cutoff<-ifelse(ici_temp$Composite_PRS>=cutoff3,1,0)
ici_temp$pheno<-ici_temp$pheno-1

b1 <- tab.paired(d=pheno, y1=somatic_cutoff, y2=germline_cutoff, data=ici_temp)
a<-pv.wgs(b1) 

b2 <- tab.paired(d=pheno, y1=somatic_cutoff, y2=composite_cutoff, data=ici_temp)
b<-pv.wgs(b2) 

b3 <- tab.paired(d=pheno, y1=germline_cutoff, y2=composite_cutoff, data=ici_temp)
c<-pv.wgs(b3) 

#plot barplot of ppv and npv?

#make df of ppv and npv respectively..?

ppv_plot<-data.frame(factor(c("Somatic","Germline","Composite"),levels=c("Somatic","Germline","Composite")),c(a$ppv[1],a$ppv[2],b$ppv[2]))
colnames(ppv_plot)<-c("Model","ppv")
stat.test<-data.frame(group1=c("Germline","Somatic","Composite"),group2=c("Somatic","Composite","Germline"),p.adj=c(a$ppv[5],b$ppv[5],c$ppv[5]))
stat.test$y.position<-rev(c(1,0.9,0.8))
stat.test$Model<-ppv_plot$Model
stat.test$p.adj<-round(stat.test$p.adj,4)

library(ggh4x)
ggplot(ppv_plot, aes(x = Model,y = ppv,fill=Model)) + 
  xlab(label=NULL) + #theme_minimal(base_size=22)+
  theme(strip.background=element_rect(color="grey30", fill="grey90"))+theme_bw(base_size = 20) +
  geom_col(color="black",width=0.95) + ylab(label="Positive Predictive Value") +
  scale_fill_manual(values=(model_cols))+ylim(0,1.05)+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),panel.spacing = unit(0.1, "lines"),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  stat_pvalue_manual(stat.test, label = "p.adj",size=7,bracket.size = 1)
#ggsave("/Users/tjsears/Code/GermTime/figs_july/fig3/ppv.pdf",width = 6,height = 6)


ppv_plot<-data.frame(factor(c("Somatic","Germline","Composite"),levels=c("Somatic","Germline","Composite")),c(a$npv[1],a$npv[2],b$npv[2]))
colnames(ppv_plot)<-c("Model","ppv")
stat.test<-data.frame(group1=c("Germline","Somatic","Composite"),group2=c("Somatic","Composite","Germline"),p.adj=c(a$npv[5],b$npv[5],c$npv[5]))
stat.test$y.position<-rev(c(1,0.9,0.8))
stat.test$Model<-ppv_plot$Model
stat.test$p.adj<-round(stat.test$p.adj,4)

library(ggh4x)
ggplot(ppv_plot, aes(x = Model,y = ppv,fill=Model)) + 
  xlab(label=NULL) + #theme_minimal(base_size=22)+
  theme(strip.background=element_rect(color="grey30", fill="grey90"))+theme_bw(base_size = 20) +
  geom_col(color="black",width=0.95) + ylab(label="Negative Predictive Value") +
  scale_fill_manual(values=(model_cols))+ylim(0,1.05)+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),panel.spacing = unit(0.1, "lines"),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  stat_pvalue_manual(stat.test, label = "p.adj",size=7,bracket.size = 1)
#ggsave("/Users/tjsears/Code/GermTime/figs_july/fig3/npv.pdf",width = 6,height = 6)



# auc stuff
table(ici_temp$somatic_cutoff,ici_temp$pheno)
table(ici_temp$germline_cutoff,ici_temp$pheno)
table(ici_temp$composite_cutoff,ici_temp$pheno)

somatic<-pROC::roc(ici, pheno,Somatic_PRS)
germline<-pROC::roc(ici, pheno,Germline_PRS)
composite<-pROC::roc(ici, pheno,Composite_PRS)

print(somatic)
print(germline)
print(composite)

roc.test(germline,somatic)

cutpointr(ici_temp,Somatic_PRS,pheno)
cutpointr(ici_temp,Germline_PRS,pheno)
cutpointr(ici_temp,Composite_PRS,pheno)


# ICI plots

# read in table

# plot PRS germline -> somatic -> composite, color by response, split by response OR response?


ici<-read.table("Data/Fig3/all_prs.txt",sep="\t",header=T)

plot_df<-ici[,c("study_cancer","normal.WXS.id","pheno","Germline_PRS","Somatic_PRS","Composite_PRS")]

plot_df_vanallen<-plot_df#[plot_df$study_cancer=="vanallen",]
plot_df_vanallen<-aggregate(cbind(plot_df_vanallen$Somatic_PRS,plot_df_vanallen$Germline_PRS,plot_df_vanallen$Composite_PRS) ~ pheno, data = plot_df_vanallen, FUN = mean, na.rm = TRUE)
colnames(plot_df_vanallen)<-c("Response","Somatic_PRS","Germline_PRS","Composite_PRS")

# Load the required libraries
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(cocor)

# Convert data to long format
data_long <- plot_df_vanallen %>%
  pivot_longer(cols=c("Germline_PRS", "Somatic_PRS", "Composite_PRS"))

data_long$name<-factor(data_long$name,levels=c("Somatic_PRS","Germline_PRS","Composite_PRS"))
data_long$Response<-c("NonResponder","NonResponder","NonResponder","Responder","Responder","Responder")

# Create the plot
ggplot(data_long, aes(x = name, y = value, color = Response)) +
  geom_point() +
  geom_line() +
  labs(x = "Variable", y = "Change", title = "Change between Continuous Variables by Group") +
  scale_color_manual(values = c("red", "blue")) +   geom_line(position = position_dodge(width = 0.4), aes(group = Response)) +
  theme_minimal()

# Display the plot


library(RColorBrewer)
cols<-brewer.pal(n = 3, name = "Dark2")

# correlation dotplots for each of the conditions

plot_df$Response<-ifelse(plot_df$pheno==1,"NonResponder","Responder")
plot_df$Response<-as.factor(plot_df$Response)  

ggplot(plot_df, aes(x = Germline_PRS, y = Somatic_PRS,color=Response)) +
  geom_jitter(width = 0, height = 0,alpha=0.7,size=2) +
  stat_cor(method = "pearson",show.legend = F,inherit.aes = F,aes(x = Germline_PRS, y = Somatic_PRS),size=6) +
  geom_smooth(method = "lm",inherit.aes = F,aes(x = Germline_PRS, y = Somatic_PRS),color = cols[3]) + 
  theme_light(base_size=16) + scale_color_manual(values=cols[2:1]) + xlab(label = "Germline IC-Index") + ylab(label = "Somatic IC-Index") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#ggsave("/Users/tjsears/Code/GermTime/figs_july/fig3/SomaticGermlineCorrelation.pdf",width = 8,height = 6)

ggplot(plot_df, aes(x = Composite_PRS, y = Somatic_PRS,color=Response)) +
  geom_jitter(width = 0, height = 0,alpha=0.7,size=2) +
  stat_cor(method = "pearson",show.legend = F,inherit.aes = F,aes(x = Composite_PRS, y = Somatic_PRS),size=6) +
  geom_smooth(method = "lm",inherit.aes = F,aes(x = Composite_PRS, y = Somatic_PRS),color = cols[3]) +
  theme_light(base_size=16) + scale_color_manual(values=cols[2:1])
#ggsave("/Users/tjsears/Code/GermTime/figs_july/fig3/SomaticCompositeCorrelation.pdf",width = 8,height = 6)

ggplot(plot_df, aes(x = Composite_PRS, y = Germline_PRS,color=Response)) +
  geom_jitter(width = 0, height = 0,alpha=0.7,size=2) +
  stat_cor(method = "pearson",show.legend = F,inherit.aes = F,aes(x = Composite_PRS, y = Germline_PRS),size=6) +
  geom_smooth(method = "lm",inherit.aes = F,aes(x = Composite_PRS, y = Germline_PRS),color = cols[3]) +
  theme_light(base_size=16) + scale_color_manual(values=cols[2:1])
#ggsave("/Users/tjsears/Code/GermTime/figs_july/fig3/CompositeGermlineCorrelation.pdf",width = 8,height = 6)


# Make colored density plots of each score for the correlation plot

ici$Response<-ifelse(ici$pheno==2,"R","NR")

library(plyr)
mu <- ddply(ici, "Response", summarise, grp.mean=mean(Somatic_PRS))
head(mu)

ggplot(ici, aes(x=Somatic_PRS, fill=Response)) +
  geom_density(alpha=0.4)+geom_vline(data=mu, aes(xintercept=grp.mean, color=Response),
                                     linetype="dashed",linewidth=1)+scale_fill_manual(values = cols[2:1])+xlab(label=NULL)+
  theme_bw(base_size = 22)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position="none")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())
#ggsave("/Users/tjsears/Code/GermTime/figs_july/fig3/SomaticICIndexDensity.pdf",width = 8,height = 2)

mu <- ddply(ici, "Response", summarise, grp.mean=mean(Germline_PRS))
head(mu)

ggplot(ici, aes(x=Germline_PRS, fill=Response)) +
  geom_density(alpha=0.4)+geom_vline(data=mu, aes(xintercept=grp.mean, color=Response),
                                     linetype="dashed",linewidth=1)+scale_fill_manual(values = cols[2:1])+xlab(label=NULL)+
  theme_bw(base_size = 22)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position="none")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())
#ggsave("/Users/tjsears/Code/GermTime/figs_july/fig3/GermlineICIndexDensity.pdf",width = 8,height = 2)


# sensitivity/specificity plot across model types
# or PPV / NPV

cutoff=5
#cutoff should reflect maximum AUC for each condition?
#so, 9 different cutoffs? one for each model x cohort
# again, to maximize AUC

#use each best cutoff to generate sensitivity and specificity for each instance

#plot in a connected dotplot and compare means, highlight key differences
library(cutpointr)

plot_df_miao<-plot_df[plot_df$study_cancer=="miao",]
miao_somatic_sensitivity<-cutpointr(plot_df_miao, Somatic_PRS, pheno, 
                                    method = maximize_loess_metric, metric = youden)$sensitivity
miao_germline_sensitivity<-cutpointr(plot_df_miao, Germline_PRS, pheno, 
                                     method = maximize_loess_metric, metric = youden)$sensitivity
miao_composite_sensitivity<-cutpointr(plot_df_miao, Composite_PRS, pheno, 
                                      method = maximize_loess_metric, metric = youden)$sensitivity


plot_df_vanallen<-plot_df[plot_df$study_cancer=="vanallen",]
vanallen_somatic_sensitivity<-cutpointr(plot_df_vanallen, Somatic_PRS, pheno, 
                                        method = maximize_loess_metric, metric = youden)$sensitivity
vanallen_germline_sensitivity<-cutpointr(plot_df_vanallen, Germline_PRS, pheno, 
                                         method = maximize_loess_metric, metric = youden)$sensitivity
vanallen_composite_sensitivity<-cutpointr(plot_df_vanallen, Composite_PRS, pheno, 
                                          method = maximize_loess_metric, metric = youden)$sensitivity


plot_df_rizvi<-plot_df[plot_df$study_cancer=="rizvi",]

rizvi_somatic_sensitivity<-cutpointr(plot_df_rizvi, Somatic_PRS, pheno, 
                                     method = maximize_loess_metric, metric = youden)$sensitivity
rizvi_germline_sensitivity<-cutpointr(plot_df_rizvi, Germline_PRS, pheno, 
                                      method = maximize_loess_metric, metric = youden)$sensitivity
rizvi_composite_sensitivity<-cutpointr(plot_df_rizvi, Composite_PRS, pheno, 
                                       method = maximize_loess_metric, metric = youden)$sensitivity



plot_df_miao<-plot_df[plot_df$study_cancer=="miao",]
miao_somatic_specificity<-cutpointr(plot_df_miao, Somatic_PRS, pheno, 
                                    method = maximize_loess_metric, metric = youden)$specificity
miao_germline_specificity<-cutpointr(plot_df_miao, Germline_PRS, pheno, 
                                     method = maximize_loess_metric, metric = youden)$specificity
miao_composite_specificity<-cutpointr(plot_df_miao, Composite_PRS, pheno, 
                                      method = maximize_loess_metric, metric = youden)$specificity


plot_df_vanallen<-plot_df[plot_df$study_cancer=="vanallen",]
vanallen_somatic_specificity<-cutpointr(plot_df_vanallen, Somatic_PRS, pheno, 
                                        method = maximize_loess_metric, metric = youden)$specificity
vanallen_germline_specificity<-cutpointr(plot_df_vanallen, Germline_PRS, pheno, 
                                         method = maximize_loess_metric, metric = youden)$specificity
vanallen_composite_specificity<-cutpointr(plot_df_vanallen, Composite_PRS, pheno, 
                                          method = maximize_loess_metric, metric = youden)$specificity


plot_df_rizvi<-plot_df[plot_df$study_cancer=="rizvi",]

rizvi_somatic_specificity<-cutpointr(plot_df_rizvi, Somatic_PRS, pheno, 
                                     method = maximize_loess_metric, metric = youden)$specificity
rizvi_germline_specificity<-cutpointr(plot_df_rizvi, Germline_PRS, pheno, 
                                      method = maximize_loess_metric, metric = youden)$specificity
rizvi_composite_specificity<-cutpointr(plot_df_rizvi, Composite_PRS, pheno, 
                                       method = maximize_loess_metric, metric = youden)$specificity


combined_sens_spec_df<-rbind(c(rizvi_somatic_specificity,rizvi_germline_specificity,rizvi_composite_specificity,
                               miao_somatic_specificity,miao_germline_specificity,miao_composite_specificity,
                               vanallen_somatic_specificity,vanallen_germline_specificity,vanallen_composite_specificity),
                             c(rizvi_somatic_sensitivity,rizvi_germline_sensitivity,rizvi_composite_sensitivity,
                               miao_somatic_sensitivity,miao_germline_sensitivity,miao_composite_sensitivity,
                               vanallen_somatic_sensitivity,vanallen_germline_sensitivity,vanallen_composite_sensitivity))

rownames(combined_sens_spec_df)<-c("specificity","sensitivity")
combined_sens_spec_df<-as.data.frame(t(combined_sens_spec_df))
combined_sens_spec_df$Cohort<-c("rizvi","rizvi","rizvi","miao","miao","miao","vanallen","vanallen","vanallen")
combined_sens_spec_df$Model<-rep(c("somatic","germline","composite"),3)

combined_sens_spec_df_plot<-pivot_longer(combined_sens_spec_df,cols=c("specificity","sensitivity"))

ggplot(combined_sens_spec_df_plot, aes(x = Model, y = value, fill = name)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Models", y = "Performance", fill = "") +
  scale_fill_manual(values = c("sensitivity" = "skyblue", "specificity" = "lightgreen")) +
  ggtitle("Comparison of Sensitivity and Specificity") +
  theme_minimal()

#ggsave("/Users/tjsears/Code/GermTime/figs_july/fig3/SensSpec.pdf",width = 8,height = 6)
















