library(RColorBrewer)
library(ggplot2)
library(ggpubr)

#### SUPPLEMENTAL Fig S5 ####

#plot response rates per group

tcga_phbr<-read.table("Data/suppFigs/composite_phbr_v4.txt",sep="\t",header=T)
plot_df<-tcga_phbr
quantiles <- quantile(plot_df$Ratio, probs = c(0,0.33333,0.66666,1))
plot_df$group <- cut(plot_df$Ratio, breaks = quantiles, labels = c("MHC-I Reliant","Balanced","MHC-II Reliant"))
plot_df$group<-factor(plot_df$group,levels = c("MHC-I Reliant","Balanced","MHC-II Reliant"))
plot_df$group[is.na(plot_df$group)]<-"MHC-I Reliant"

sd<-as.data.frame(table(plot_df$group,plot_df$response_crist_sd))

library(ggplot2)

sd$Response=c(rep("NonResponders",3),rep("Responders",3))
sd$Group=sd$Var1
colnames(sd)[3]<-"Patients"
sd$Value<-c(round((sd$Patients[1]/sum(sd$Patients[c(1,4)]))*100,2),round((sd$Patients[2]/sum(sd$Patients[c(2,5)]))*100,2),round((sd$Patients[3]/sum(sd$Patients[c(3,6)]))*100,2),
            round((sd$Patients[4]/sum(sd$Patients[c(4,1)]))*100,2),round((sd$Patients[5]/sum(sd$Patients[c(5,2)]))*100,2),round((sd$Patients[6]/sum(sd$Patients[c(6,3)]))*100,2))

#add N to labels
group_sizes<-aggregate(Patients~Group,data=sd,FUN=sum)

#Set labels
library(ggpubr)
my_comparisons <- list( c("MHC-II Reliant", "MHC-I Reliant"), c("MHC-II Reliant", "Balanced"), c("MHC-I Reliant", "Balanced"))

chi.test <- function(a, b) {
  return(chisq.test(cbind(a, b)))
}

ggplot(sd, aes(fill=Response, y=Value, x=Group,label=Value)) + scale_fill_manual(values = cols<-rev(c("#B4CCE3","#FBB4AE")))  +
  geom_bar(position="stack", stat="identity",color='black') +theme_bw(base_size = 20) +xlab(label = NULL)+ylab(label = NULL) +
  geom_signif(comparisons = my_comparisons,test = "chi.test",map_signif_level=F,y_position=c(105,114,110),tip_length = 0.04) +
  ggtitle(label = "Response by Reliance Grouping") + geom_text(size = 5, position = position_stack(vjust = 0.5)) +
  theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),panel.grid.major.x = element_blank(),panel.grid.major.y = element_blank()) +
  scale_x_discrete(labels=c("MHC-II Reliant" = paste("MHC-II Reliant\n","(n=",group_sizes$Patients[3],")",sep="")
                            , "MHC-I Reliant" = paste("MHC-I Reliant\n","(n=",group_sizes$Patients[1],")",sep="")
                            ,"Balanced" = paste("Balanced\n","(n=",group_sizes$Patients[2],")",sep="")))
#ggsave("/Users/tjsears/Code/GermTime/figs_july/figS5/RelianceResponseRates_Partial.pdf",width = 8,height = 6)


#########################
# MHC Reliance dynamics #
#########################

cols<-brewer.pal(3,"Dark2")

plot_df<-tcga_phbr
quantiles <- quantile(plot_df$Ratio, probs = c(0,0.33,0.66,1))
plot_df$group <- cut(plot_df$Ratio, breaks = quantiles, labels = c("MHC-I Reliant","Balanced","MHC-II Reliant"))
plot_df$group<-factor(plot_df$group,levels = c("MHC-I Reliant","Balanced","MHC-II Reliant"))
plot_df$group[is.na(plot_df$group)]<-"MHC-I Reliant"
plot_df<-plot_df[!((plot_df$Escape+plot_df$numLostAlleles)==0&plot_df$group=="MHC-II Reliant"),]

plot_df<-plot_df[order(plot_df$Ratio),]
plot_df<-plot_df[!duplicated(plot_df$individual),]

plot_df$individual<-factor(plot_df$individual,levels=plot_df$individual)
rownames(plot_df)<-plot_df$individual

b <- ggplot(plot_df, aes(x=individual, y=Ratio, fill=group,color=group)) +
  scale_fill_discrete(name="MHC Reliance Group") + scale_color_manual(values=cols[c(2,3,1)]) +
  labs(list(title = "Waterfall plot for changes in QoL scores", x = NULL, y = "Change from baseline (%) in QoL score")) +
  theme_classic(base_size = 14) + 
  theme(axis.line.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.title.y = element_text(face="bold",angle=90),panel.background = element_blank()) +
  coord_cartesian(ylim = c(0,3)) + ggtitle(label="Neoantigen Presentability Ratios")
b <- b + geom_bar(stat="identity", width=0.7, position = position_dodge(width=0.4))
b

#ggsave("/Users/tjsears/Code/GermTime/figs_july/figS5/waterfallPHBR.pdf",width = 9,height = 3)

# check how common MHC damage is in each group

plot_df<-tcga_phbr
quantiles <- quantile(plot_df$Ratio, probs = c(0,0.33,0.66,1))
plot_df$group <- cut(plot_df$Ratio, breaks = quantiles, labels = c("MHC-I Reliant","Balanced","MHC-II Reliant"))
plot_df$group<-factor(plot_df$group,levels = c("MHC-I Reliant","Balanced","MHC-II Reliant"))
plot_df$group[is.na(plot_df$group)]<-"MHC-I Reliant"

plot_df$AntigenPresentationDamage<-plot_df$Escape+plot_df$numLostAlleles

ggplot(plot_df, aes(x = group,y = AntigenPresentationDamage,fill=group)) + 
  xlab(label=NULL) +
  ylab(label="MHC-I Presentation\nPathway Damage\n(proportion total muts)")+theme_minimal(base_size=14)+
  theme(strip.background=element_rect(color="grey30", fill="grey90"))+theme_bw(base_size =18) +
  geom_boxplot(position = position_dodge(width = 1),width=0.98) +
  scale_fill_manual(values=brewer.pal(3,"Dark2")[c(2,3,1)]) +
  theme(panel.spacing = unit(0.1, "lines"),)+ylim(0,0.11)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()
        , panel.background = element_blank()
        , axis.line = element_line(colour = "black"),legend.position="none",
        ,axis.text = element_text(color="black"))+
  stat_compare_means(comparisons = list(c("MHC-I Reliant","MHC-II Reliant"),c("Balanced","MHC-II Reliant"),c("MHC-I Reliant","Balanced")),size=7,
                     method="t.test", method.args=list(var.equal = T))
#ggsave("/Users/tjsears/Code/GermTime/figs_july/figS5/MHC_I_damage_cont.pdf",width = 5.5,height = 4)

# stacked barplot of presence of MHC-I Presentation pathway damage

plot_df$MHC_APP_bin<-ifelse(plot_df$AntigenPresentationDamage>0,"Damaged","Intact")

sd<-as.data.frame(table(plot_df$group,plot_df$MHC_APP_bin))

library(ggplot2)

sd$Response=c(rep("Damaged",3),rep("Intact",3))
sd$Group=sd$Var1
colnames(sd)[3]<-"Patients"
sd$Value<-c(round((sd$Patients[1]/sum(sd$Patients[c(1,4)]))*100,2),round((sd$Patients[2]/sum(sd$Patients[c(2,5)]))*100,2),round((sd$Patients[3]/sum(sd$Patients[c(3,6)]))*100,2),
            round((sd$Patients[4]/sum(sd$Patients[c(4,1)]))*100,2),round((sd$Patients[5]/sum(sd$Patients[c(5,2)]))*100,2),round((sd$Patients[6]/sum(sd$Patients[c(6,3)]))*100,2))

#add N to labels
group_sizes<-aggregate(Patients~Group,data=sd,FUN=sum)

#Set labels
library(ggpubr)
my_comparisons <- list( c("MHC-II Reliant", "MHC-I Reliant"), c("MHC-II Reliant", "Balanced"), c("MHC-I Reliant", "Balanced"))

chi.test <- function(a, b) {
  return(chisq.test(cbind(a, b)))
}

ggplot(sd, aes(fill=Response, y=Value, x=Group,label=Value)) + scale_fill_manual(values = cols<-(c("#B05D4A","#579F9F")))  +
  geom_bar(position="stack", stat="identity",color='black') +theme_bw(base_size = 16) +xlab(label = NULL)+ylab(label = NULL) +
  geom_signif(comparisons = my_comparisons,test = "chi.test",map_signif_level=F,y_position=c(101,112,109),tip_length = 0.04) +
  ggtitle(label = "MHC-I APP Damage by Reliance Group") + geom_text(size = 5, position = position_stack(vjust = 0.5)) +
  theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),panel.grid.major.x = element_blank(),panel.grid.major.y = element_blank()) +
  scale_x_discrete(labels=c("MHC-II Reliant" = paste("MHC-II Reliant\n","(n=",group_sizes$Patients[3],")",sep="")
                            , "MHC-I Reliant" = paste("MHC-I Reliant\n","(n=",group_sizes$Patients[1],")",sep="")
                            ,"Balanced" = paste("Balanced\n","(n=",group_sizes$Patients[2],")",sep="")))
#ggsave("/Users/tjsears/Code/GermTime/figs_july/figS5/MHC_I_damage.pdf",width = 5.2,height = 4)
