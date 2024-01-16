library(RColorBrewer)
library(ggplot2)
library(ggpubr)
# read in data

data2=read.table("Data/Fig5/tfh_HLA_pr.txt",sep="\t",header=T)[,c(3,4)]
data=read.table("Data/Fig5/tfh_HLA_sd.txt",sep="\t",header=T)[,c(3,4)]

sd=as.data.frame((data2-data)[,1])
sd$Response=rep("Stable Disease",nrow(sd))
sd$Group=c("Neither","PDCD1","TFH","Both")
colnames(sd)[1]<-"Patients"
sd$Value<-round((sd$Patients/sum(sd$Patients))*100,2)

partial=as.data.frame(data2[,2])
partial$Response=rep("Partial/Complete Responders",nrow(partial))
partial$Group=c("Neither","PDCD1","TFH","Both")
colnames(partial)[1]<-"Patients"
partial$Value<-round((partial$Patients/sum(partial$Patients))*100,2)

nonresponders=as.data.frame(data[,1])
nonresponders$Response=rep("Nonresponders",nrow(nonresponders))
nonresponders$Group=c("Neither","PDCD1","TFH","Both")
colnames(nonresponders)[1]<-"Patients"
nonresponders$Value<-round((nonresponders$Patients/sum(nonresponders$Patients))*100,2)

responder_data<-rbind(nonresponders,sd,partial)

responder_data$Group<-factor(responder_data$Group,levels=c("Neither","PDCD1","TFH","Both"))
responder_data$Response<-factor(responder_data$Response,levels=rev(c("Nonresponders","Stable Disease","Partial/Complete Responders")))

sum_vec=rep(c(sum(responder_data[responder_data$Group=="Neither",1]),sum(responder_data[responder_data$Group=="PDCD1",1]),sum(responder_data[responder_data$Group=="TFH",1]),sum(responder_data[responder_data$Group=="Both",1])),3)
responder_data$Percent<-round((responder_data$Patients/sum_vec)*100,2)

responder_test<-cbind(nonresponders,sd,partial)[,c(1,5,9)]
chisquared_value<-chisq.test(responder_test)


#add N to labels
group_sizes<-aggregate(Patients~Group,data=responder_data,FUN=sum)

#pairwise chisquared labels?
chisq.test(responder_test[c(3,4),])


#Set labels
library(ggpubr)
my_comparisons <- list( c("Neither", "Both"), c("Both", "PDCD1"), c("Both", "TFH"),c("Neither","PDCD1"),c("PDCD1","TFH"),c("TFH","Neither") )

chi.test <- function(a, b) {
  return(chisq.test(cbind(a, b)))
}


ggplot(responder_data, aes(fill=Response, y=Percent, x=Group,label=Percent)) + #scale_fill_manual(values = c())  +
  geom_bar(position="stack", stat="identity") +theme_bw(base_size = 20) +xlab(label = NULL)+ylab(label = NULL) +
  geom_signif(comparisons = my_comparisons,test = "chi.test",map_signif_level=F,y_position=c(125,119,113,99,110,104),tip_length = 0) +
  ggtitle(label = "Responder by SNP grouping",subtitle = paste("ChiSquared pval=",round(chisquared_value$p.value,6),sep="")) + geom_text(size = 5, position = position_stack(vjust = 0.5)) +
  theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),panel.grid.major.x = element_blank(),panel.grid.major.y = element_blank()) +
  scale_x_discrete(labels=c("Neither" = paste("Neither\n","(n=",group_sizes$Patients[1],")",sep="")
                            , "PDCD1" = paste("HLA Damage\n","(n=",group_sizes$Patients[2],")",sep="")
                            ,"TFH" = paste("TFH SNP\n","(n=",group_sizes$Patients[3],")",sep="")
                            ,"Both" = paste("Both\n","(n=",group_sizes$Patients[4],")",sep="")))
#ggsave("/Users/tjsears/Code/GermTime/figs_july/fig5_v3/HLA_TFH.pdf",width = 10,height = 8)



# read in tcga phbr
tcga_phbr<-read.table("Data/Fig5/composite_phbr_v4.txt",sep="\t",header=T)
Liu<-read.table("Data/Fig5/liu_composite_phbr_v4.txt",sep="\t",header=T)

#export tables for supp
plot_df<-tcga_phbr
quantiles <- quantile(plot_df$Ratio, probs = c(0,0.33,0.66,1))
plot_df$group <- cut(plot_df$Ratio, breaks = quantiles, labels = c("MHC-I Reliant","Balanced","MHC-II Reliant"))
plot_df$group<-factor(plot_df$group,levels = c("MHC-I Reliant","Balanced","MHC-II Reliant"))
plot_df$group[is.na(plot_df$group)]<-"MHC-I Reliant"
plot_df$`Defects in MHC-I antigen presentation pathway`<-ifelse(plot_df$Escape+plot_df$numLostAlleles==0,"yes","no")
plot_df$Response<-ifelse(plot_df$pheno==2,"R","NR")

tcga_phbr_export<-plot_df[,c("tumor.WXS.id","PHBR1","PHBR2","group","Defects in MHC-I antigen presentation pathway","Response")]
#write.table(tcga_phbr_export,"/Users/tjsears/Code/GermTime/FinalSuppTables/DiscoveryPHBR.txt",sep="\t",quote = F,row.names = F)

plot_df<-Liu
quantiles <- quantile(plot_df$Ratio, probs = c(0,0.33,0.66,1))
plot_df$group <- cut(plot_df$Ratio, breaks = quantiles, labels = c("MHC-I Reliant","Balanced","MHC-II Reliant"))
plot_df$group<-factor(plot_df$group,levels = c("MHC-I Reliant","Balanced","MHC-II Reliant"))
plot_df$group[is.na(plot_df$group)]<-"MHC-I Reliant"
plot_df$`Defects in MHC-I antigen presentation pathway`<-ifelse(plot_df$Escape+plot_df$numLostAlleles==0,"yes","no")
plot_df$Response<-ifelse(plot_df$pheno==2,"R","NR")

liu_export<-plot_df[,c("tumor.WXS.id","PHBR1","PHBR2","group","Defects in MHC-I antigen presentation pathway","Response")]
#write.table(liu_export,"/Users/tjsears/Code/GermTime/FinalSuppTables/ValidationPHBR.txt",sep="\t",quote = F,row.names = F)

###########################################
####### hazard plot for MHC groups ########
###########################################

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

cols2<-brewer.pal(3,"Pastel2")

plot_df<-tcga_phbr

# Feature importance
featureImportance<-read.table("Data/Fig5/featureImportanceMatrix.txt",sep="\t",header=T,row.names = 1)
colnames(featureImportance)[1]<-"TFH-cells"
rownames(featureImportance)[1]<-"TFH-cells"

cols=brewer.pal(5,"Set2")
library(reshape2)
rownames(featureImportance)<-c("TFH QTL "      ,   "ERAP2"       ,      "DCTN5"        ,     "DHFR"          ,    "Subclonal TMB",
                                "ImmunoEditing"   ,    "GPLD1"      ,       "Immune Escape"      ,      "MHC Damage" ,   "APP"       ,
                                "ITGB2"         ,    "FCGR3B"     ,       "PDCD1"       ,      "TREX1" ,            "VAMP8"            ,
                                "FCGR2B"        ,    "CTSW"       ,       "ERAP1"       ,      "CTSS" ,             "ITH",
                                "FAM167A"       ,    "LYZ"        ,       "FPR1"        ,      "VAMP3"  )
colnames(featureImportance)<-c("TFH QTL "      ,   "ERAP2"       ,      "DCTN5"        ,     "DHFR"          ,    "Subclonal TMB",
                               "ImmunoEditing"   ,    "GPLD1"      ,       "Immune Escape"      ,      "MHC Damage" ,   "APP"       ,
                               "ITGB2"         ,    "FCGR3B"     ,       "PDCD1"       ,      "TREX1" ,            "VAMP8"            ,
                               "FCGR2B"        ,    "CTSW"       ,       "ERAP1"       ,      "CTSS" ,             "ITH",
                               "FAM167A"       ,    "LYZ"        ,       "FPR1"        ,      "VAMP3"  )

feet_long<-melt(featureImportance)
feet_long$AltVar<-colnames(featureImportance)
feet_long$FinalVar<-paste(feet_long$variable,"/",feet_long$AltVar)

feet_long<-feet_long[order(feet_long$value,decreasing = T),]
feet_long<-feet_long[c(1,3,5,7,9,11,13,15,17,19,21,23,25),]
feet_long$FinalVar<-factor(feet_long$FinalVar,levels=c(feet_long$FinalVar))
feet_long$InteractionType<-c("Immune Infiltration / Antigen Presentation","Immune Infiltration / Immune Signaling","Antigen Processing",
                             "Antigen Processing","Immune Infiltration / Antigen Processing","Antigen Processing / Antigen Presentation",
                             "Immune Infiltration / Antigen Processing","Antigen Processing / Antigen Presentation","Antigen Processing / Antigen Presentation",
                             "Immune Infiltration / Antigen Presentation","Antigen Processing / Antigen Presentation","Immune Infiltration / Antigen Processing",
                             "Antigen Processing")

ggplot(feet_long, aes(x = FinalVar, y = value, fill = InteractionType)) + ylab(label="SHAP Interaction Value")+
  geom_bar(stat = "identity") + theme_minimal(base_size = 16) + xlab(label=NULL)+
  scale_fill_manual(values = cols) + theme(axis.text.x = element_text(angle = 70, vjust = 1.07, hjust=1),
                                           panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#ggsave("/Users/tjsears/Code/GermTime/figs_july/fig5_v3/SimplifiedInteraction.pdf",width = 10,height = 5)


# boxplot of hla damage and TFH snp area plot based on MHC classification
cols<-rev(c("#B4CCE3","#FBB4AE"))

plot_df<-tcga_phbr
plot_df$Response<-ifelse(plot_df$pheno==2,"R","NR")
plot_df$shapGroup<-plot_df$group
plot_df$MHCdmg<-ifelse(plot_df$numLostAlleles>0.0075,1,0)

quantiles <- quantile(plot_df$Ratio, probs = c(0,0.33,0.66,1))
plot_df$group <- cut(plot_df$Ratio, breaks = quantiles, labels = c("MHC-I Reliant","Balanced","MHC-II Reliant"))
plot_df$group<-factor(plot_df$group,levels = c("MHC-I Reliant","Balanced","MHC-II Reliant"))
plot_df$group[is.na(plot_df$group)]<-"MHC-I Reliant"

library(ggh4x)
#reformat plot_df such that each column represents a % of the total proportion of pts containing the tfh snp
plot_df_short<-plot_df
plot_df_short$aggCol<-paste(plot_df_short$group,plot_df_short$Response)

#possibly reduce tfh snp to binary
x<-table(plot_df_short$aggCol)
plot_df_short<-aggregate(plot_df_short$MHCdmg,by=list(plot_df_short$aggCol),FUN=sum)
plot_df_short$MHCdmg<-plot_df_short$x/x
plot_df_short$Response<-c("NR",'R',"NR",'R',"NR",'R')
plot_df_short$group<-factor(c("Balanced",'Balanced',"MHC-I Reliant",'MHC-I Reliant',"MHC-II Reliant",'MHC-II Reliant'),
                            levels=c("MHC-I Reliant",'Balanced','MHC-II Reliant'))

ggplot(plot_df_short, aes(x = Response,y = MHCdmg,fill=Response,group=Response)) + 
  xlab(label=NULL) + #scale_y_continuous(trans = "log2") +
  facet_nested(~group)+ylab(label="Proportion of patients\nwith MHC-I damage")+theme_minimal(base_size=14)+
  theme(strip.background=element_rect(color="grey30", fill="grey90"))+theme_bw(base_size = 20) +
  geom_col(position='stack',color="black") +
  scale_fill_manual(values=cols) + 
  theme(panel.spacing = unit(0.1, "lines"),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  stat_compare_means(comparisons = list(c("NR","R")),size=5.2,tip.length = 0.01,bracket.size = 0.7,
                     method="t.test", method.args=list(var.equal = F),label.y = c(0.24,0.24,0.24),data = plot_df) #+ ylim(0,0.5)
#ggsave("/Users/tjsears/Code/GermTime/figs_july/fig5_v3/HLAdmg_compopsiteV2.pdf",width = 8,height = 6)

# TFH AF
plot_df<-tcga_phbr
plot_df$Response<-ifelse(plot_df$pheno==2,"R","NR")
plot_df$shapGroup<-plot_df$group

quantiles <- quantile(plot_df$Ratio, probs = c(0,0.33,0.66,1))
plot_df$group <- cut(plot_df$Ratio, breaks = quantiles, labels = c("MHC-I Reliant","Balanced","MHC-II Reliant"))
plot_df$group<-factor(plot_df$group,levels = c("MHC-I Reliant","Balanced","MHC-II Reliant"))
plot_df$group[is.na(plot_df$group)]<-"MHC-I Reliant"
plot_df<-plot_df[!((plot_df$Escape+plot_df$numLostAlleles)==0&plot_df$group=="MHC-II Reliant"),]

library(ggh4x)
#reformat plot_df such that each column represents a % of the total proportion of pts containing the tfh snp
plot_df_short<-plot_df
plot_df_short$aggCol<-paste(plot_df_short$group,plot_df_short$Response)

#possibly reduce tfh snp to binary
x<-table(plot_df_short$aggCol)
plot_df_short<-aggregate(plot_df_short$TFH_snp,by=list(plot_df_short$aggCol),FUN=sum)
plot_df_short$TFH_snp<-plot_df_short$x/x
plot_df_short$Response<-c("NR",'R',"NR",'R',"NR",'R')
plot_df_short$group<-factor(c("Balanced",'Balanced',"MHC-I Reliant",'MHC-I Reliant',"MHC-II Reliant",'MHC-II Reliant'),
                            levels=c("MHC-I Reliant",'Balanced','MHC-II Reliant'))

ggplot(plot_df_short, aes(x = Response,y = TFH_snp,fill=Response,group=Response)) + 
  xlab(label=NULL) + 
  facet_nested(~group)+ylab(label="TFH qTL Allele Frequency")+theme_minimal(base_size=14)+
  theme(strip.background=element_rect(color="grey30", fill="grey90"))+theme_bw(base_size = 20) +
  geom_col(position='stack',color="black") +
  scale_fill_manual(values=cols) + 
  theme(panel.spacing = unit(0.1, "lines"),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  stat_compare_means(comparisons = list(c("NR","R")),size=5.5,tip.length = 0.01,bracket.size = 0.7,
                     method="t.test", method.args=list(var.equal = T),label.y = c(0.33,0.33,0.33),data = plot_df) #+ ylim(0,0.5)
#ggsave("/Users/tjsears/Code/GermTime/figs_july/fig5_v3/TFHsnp_compopsite.pdf",width = 8,height = 6)

# boxplot of TFH infiltration by response by group

ICI_phbr<-read.table("Data/Fig5/fig6v4_ici_rnaseq.txt",sep="\t",header=T)
plot_df<-ICI_phbr
plot_df$Response<-ifelse(plot_df$pheno==2,"R","NR")

quantiles <- quantile(plot_df$Ratio, probs = c(0,0.33,0.66,1))
plot_df$group <- cut(plot_df$Ratio, breaks = quantiles, labels = c("MHC-I Reliant","Balanced","MHC-II Reliant"))
plot_df$group<-factor(plot_df$group,levels = c("MHC-I Reliant","Balanced","MHC-II Reliant"))
plot_df$group[is.na(plot_df$group)]<-"MHC-I Reliant"
plot_df<-plot_df[!((plot_df$Escape+plot_df$numLostAlleles)==0&plot_df$group=="MHC-II Reliant"),]

library(ggh4x)
ggplot(plot_df, aes(x = Response,y = T.cells.follicular.helper,fill=Response,group=Response)) + 
  xlab(label=NULL) + #scale_y_continuous(trans = "log2") +
  facet_nested(~group)+ylab(label="TFH cell infiltration")+theme_minimal(base_size=14)+
  theme(strip.background=element_rect(color="grey30", fill="grey90"))+theme_bw(base_size = 20) +
  geom_boxplot(position = position_dodge(width = 1),width=0.98) +
  scale_fill_manual(values=(cols)) +
  theme(panel.spacing = unit(0.1, "lines"),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  stat_compare_means(comparisons = list(c("NR","R")),size=7,bracket.size = 0.7,
                     method="t.test", method.args=list(var.equal = F),label.y = c(3.5,3.5,3.5)) 
#ggsave("/Users/tjsears/Code/GermTime/figs_july/figS5/tfh_inf_compopsite.pdf",width = 8,height = 6)

#auto-generate TLS signature
ICI_phbr<-read.table("Data/Fig5/fig6v5_ici_rnaseq.txt",sep="\t",header=T)
plot_df<-ICI_phbr
plot_df$Response<-ifelse(plot_df$pheno==2,"R","NR")

quantiles <- quantile(plot_df$Ratio, probs = c(0,0.33,0.66,1))
plot_df$group <- cut(plot_df$Ratio, breaks = quantiles, labels = c("MHC-I Reliant","Balanced","MHC-II Reliant"))
plot_df$group<-factor(plot_df$group,levels = c("MHC-I Reliant","Balanced","MHC-II Reliant"))
plot_df$group[is.na(plot_df$group)]<-"MHC-I Reliant"
plot_df<-plot_df[!((plot_df$Escape+plot_df$numLostAlleles)==0&plot_df$group=="MHC-II Reliant"),]

TLS_sig_df<-read.table("Data/Fig4/TLS_sig.txt",sep="\t",header=T,row.names = 1)
df_<-colMeans(TLS_sig_df)

#merge TLS results back into df
df_<-as.data.frame(df_)
colnames(df_)<-c("TLS_sig")
plot_df2<-merge(plot_df,df_,by.x="tumor.RNA.id",by.y=0)
plot_df2$tissue<-ifelse(plot_df2$study_cancer_x=="miao",1,0)

plot_df2 = plot_df2 %>%
  group_by(study_cancer_x) %>%
  mutate(TLS_sig = scale(TLS_sig))

plot_df2$TLS_sig<-plot_df2$TLS_sig+abs(min(plot_df2$TLS_sig))+1
  
library(ggh4x)
ggplot(plot_df2, aes(x = Response,y = TLS_sig,fill=Response,group=Response)) + 
  xlab(label=NULL) + scale_y_continuous(trans = "log2") +
  facet_nested(~group)+ylab(label="Log2 TLS signature")+theme_minimal(base_size=14)+
  theme(strip.background=element_rect(color="grey30", fill="grey90"))+theme_bw(base_size = 20) +
  geom_boxplot(position = position_dodge(width = 1),width=0.98) +
  scale_fill_manual(values=cols) +
  theme(panel.spacing = unit(0.1, "lines"),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  stat_compare_means(comparisons = list(c("NR","R")),size=7,label.y = c(1.6,1.6,1.6),bracket.size = 0.7,
                     method="wilcox.test", method.args=list(var.equal = T))
#ggsave("/Users/tjsears/Code/GermTime/figs_july/fig5_v3/TLS_compopsite.pdf",width = 8,height = 6)

#between macro groups
ggplot(plot_df2, aes(x = group,y = TLS_sig,fill=group)) + 
  xlab(label=NULL) +
  ylab(label="TLS signature")+theme_minimal(base_size=14)+
  theme(strip.background=element_rect(color="grey30", fill="grey90"))+theme_bw(base_size = 20) +
  geom_boxplot(position = position_dodge(width = 1),width=0.98) +
  scale_fill_manual(values=brewer.pal(3,"Dark2")[c(2,3,1)]) +
  theme(panel.spacing = unit(0.1, "lines"))+
  stat_compare_means(comparisons = list(c("MHC-I Reliant","MHC-II Reliant"),c("Balanced","MHC-II Reliant"),c("MHC-I Reliant","Balanced")),size=7,
                     method="wilcox.test", method.args=list(var.equal = T))

# use pval from above fig to manually annotate plot :(
# done

# Scatterplot of PHBR scores with colored regions indicating reliance
tcga_phbr<-read.table("Data/Fig5/composite_phbr_v4.txt",sep="\t",header=T)
plot_df<-tcga_phbr
plot_df$Response<-ifelse(plot_df$pheno==2,"R","NR")
plot_df$PHBR1<-log2(plot_df$PHBR1)
plot_df$PHBR2<-log2(plot_df$PHBR2)

#get slope of PHBR1

#poly_df<-data.frame(x=c(0,0,11),y=c(1,10,10))
poly_df<-data.frame(x=c(0,0,12.2),y=c(2.21,10,10))

#poly_df2<-data.frame(x=c(0,1,12,12,11,0),y=c(0,0,9,10,10,1))
poly_df2<-data.frame(x=c(0,12.2,12.2,0),y=c(2.21,10,7.61,1.94))

#poly_df3<-data.frame(x=c(1,12,12),y=c(0,0,9))
poly_df3<-data.frame(x=c(0,0,12.2,12.2),y=c(1.94,1,1,7.61))

#slope1 = 0.45 b1=2.21  upper coord = x = 12.20727
#slope2 = 0.66 b2=1.9432 upper coord = y = 7.61
quantiles <- quantile(plot_df$Ratio, probs = c(0,0.33333,0.66666,1))
plot_df$group <- cut(plot_df$Ratio, breaks = quantiles, labels = c("MHC-I Reliant","Balanced","MHC-II Reliant"))

ggplot(plot_df, aes(x = PHBR1, y = PHBR2,color="grey95")) +
  geom_jitter(width = 2.5, height = 0,show.legend = F) +
  labs(x = "log2 MHC-I Neoantigens", y = "log2 MHC-II Neoantigens") +
  ggtitle("MHC Response Type") + theme(legend.position = "none") + 
  theme_bw(base_size=16) + scale_color_manual(values="grey30")  + xlim(c(0,12.3))+ ylim(1,10)+
  geom_polygon(data=poly_df,aes(x=x,y=y),inherit.aes = F,fill=c("#1A9E76"),alpha = 0.2) +
  geom_polygon(data=poly_df3,aes(x=x,y=y),inherit.aes = F,fill=c("#D95E01"),alpha = 0.2)+
  geom_polygon(data=poly_df2,aes(x=x,y=y),inherit.aes = F,fill=c("#766FB3"),alpha = 0.2)+
  #geom_segment(aes(x = 0, xend = 11, y = 0, yend = 11),colour = "grey40",linetype=2)+
  #geom_segment(aes(x = 0, xend = 12, y = 2.21, yend = 9),colour = "grey40",linetype=2)+
  annotate("label", x = 4, y = 8.5, label = "MHC-II Reliant",size=5) +
  annotate("label", x = 11, y = 8.2, label = "Balanced",size=5) +
  annotate("label", x = 9.5, y = 3, label = "MHC-I Reliant",size=5) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
#ggsave("/Users/tjsears/Code/GermTime/figs_july/fig5_v3/response_schematic.pdf",width = 5,height = 5)



