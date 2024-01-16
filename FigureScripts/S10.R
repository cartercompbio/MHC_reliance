#S10
#################################################################
#### checkpoint gene assoc with MHC groupings! supp fig 10ish ####
#################################################################
library(ggplot2)

ICI_phbr<-read.table("Data/Fig5/fig6v4_ici_rnaseq.txt",sep="\t",header=T)
quantiles <- quantile(ICI_phbr$Ratio, probs = c(0,0.33,0.66,1))
ICI_phbr$group <- cut(ICI_phbr$Ratio, breaks = quantiles, labels = c("MHC-I Reliant","Balanced","MHC-II Reliant"))
ICI_phbr$group[is.na(ICI_phbr$group)]<-"MHC-I Reliant"
ICI_phbr$Response<-ifelse(ICI_phbr$response_crist_sd==2,"R","NR")

ICI_phbr$Neoantigens<-ICI_phbr$PHBR1+ICI_phbr$PHBR2

library(tidyr)
cols<-brewer.pal(3,"Dark2")

ICI_phbr_longer<-pivot_longer(ICI_phbr,cols = c("LAG3","CD274","CTLA4"))

p<-ggplot(ICI_phbr_longer, aes(x = group,y = value,fill=group,group=group)) + ggtitle(label = "Discovery") +
  xlab(label=NULL) + 
  facet_nested(~name)+ylab(label="Checkpoint Expression")+theme_minimal(base_size=14)+
  theme(strip.background=element_rect(color="grey30", fill="grey90"))+theme_bw(base_size = 20) +
  geom_boxplot(position = position_dodge(width = 1),width=0.98) +
  scale_fill_manual(values=cols[c(2,3,1)]) + 
  theme(panel.spacing = unit(0.1, "lines"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+
  stat_compare_means(comparisons = list(c("MHC-I Reliant","MHC-II Reliant")),size=7,bracket.size = 0.7,
                     method="t.test", method.args=list(var.equal = F),label.y = c(12,12,12)) 
p

#ggsave("/Users/tjsears/Code/GermTime/figs_july/figS5/DiscoveryGermlineAssoc.pdf",width = 8,height = 6)

ICI_phbr_longer<-pivot_longer(ICI_phbr,cols = c("Neoantigens"))

p<-ggplot(ICI_phbr_longer, aes(x = group,y = value,fill=group,group=group)) + 
  xlab(label=NULL) + #scale_y_continuous(trans = "log2") 
  facet_nested(~name)+ylab(label="Neoantigens")+theme_minimal(base_size=14)+ggtitle(label = "Discovery")+
  theme(strip.background=element_rect(color="grey30", fill="grey90"))+theme_bw(base_size = 20) +
  geom_boxplot(position = position_dodge(width = 1),width=0.98) +
  scale_fill_manual(values=cols[c(2,3,1)]) + 
  theme(panel.spacing = unit(0.1, "lines"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+
  stat_compare_means(comparisons = list(c("MHC-I Reliant","MHC-II Reliant")
                                        ,c("MHC-I Reliant","Balanced")
                                        ,c("Balanced","MHC-II Reliant")),size=7,bracket.size = 0.7,
                     method="t.test", method.args=list(var.equal = F),label.y = c(1100,1300,1500)) 
p
#ggsave("/Users/tjsears/Code/GermTime/figs_july/figS5/DiscoveryNeoantigen.pdf",width = 8,height = 6)


tcga_phbr<-read.table("Data/Fig6/fig6v4_liu_rnaseq.txt",sep="\t",header=T)
quantiles <- quantile(tcga_phbr$Ratio, probs = c(0,0.33,0.66,1))
tcga_phbr$group <- cut(tcga_phbr$Ratio, breaks = quantiles, labels = c("MHC-I Reliant","Balanced","MHC-II Reliant"))
tcga_phbr$group[is.na(tcga_phbr$group)]<-"MHC-I Reliant"
tcga_phbr$Response<-ifelse(tcga_phbr$response_crist_sd==2,"R","NR")

tcga_phbr$Neoantigens<-tcga_phbr$PHBR1+tcga_phbr$PHBR2

library(tidyr)
cols<-brewer.pal(3,"Dark2")

tcga_phbr_longer<-pivot_longer(tcga_phbr,cols = c("LAG3","CD274","CTLA4"))

p<-ggplot(tcga_phbr_longer, aes(x = group,y = value,fill=group,group=group)) + 
  xlab(label=NULL) + #scale_y_continuous(trans = "log2") 
  facet_nested(~name)+ylab(label="Checkpoint Expression")+theme_minimal(base_size=14)+ggtitle(label = "Validation")+
  theme(strip.background=element_rect(color="grey30", fill="grey90"))+theme_bw(base_size = 20) +
  geom_boxplot(position = position_dodge(width = 1),width=0.98) +
  scale_fill_manual(values=cols[c(2,3,1)]) + 
  theme(panel.spacing = unit(0.1, "lines"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+
  stat_compare_means(comparisons = list(c("MHC-I Reliant","MHC-II Reliant")),size=7,bracket.size = 0.7,
                     method="t.test", method.args=list(var.equal = F),label.y = c(150,150,150)) 
p
#ggsave("/Users/tjsears/Code/GermTime/figs_july/figS5/ValidationGermlineAssoc.pdf",width = 8,height = 6)

tcga_phbr_longer<-pivot_longer(tcga_phbr,cols = c("Neoantigens"))

p<-ggplot(tcga_phbr_longer, aes(x = group,y = value,fill=group,group=group)) + 
  xlab(label=NULL) + #scale_y_continuous(trans = "log2") 
  facet_nested(~name)+ylab(label="Neoantigens")+theme_minimal(base_size=14)+ggtitle(label = "Validation")+
  theme(strip.background=element_rect(color="grey30", fill="grey90"))+theme_bw(base_size = 20) +
  geom_boxplot(position = position_dodge(width = 1),width=0.98) +
  scale_fill_manual(values=cols[c(2,3,1)]) + 
  theme(panel.spacing = unit(0.1, "lines"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+
  stat_compare_means(comparisons = list(c("MHC-I Reliant","MHC-II Reliant")
                                        ,c("MHC-I Reliant","Balanced")
                                        ,c("Balanced","MHC-II Reliant")),size=7,bracket.size = 0.7,
                     method="t.test", method.args=list(var.equal = F),label.y = c(540,650,700)) 
p
#ggsave("/Users/tjsears/Code/GermTime/figs_july/figS5/ValidationNeoantigen.pdf",width = 8,height = 6)

