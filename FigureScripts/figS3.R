# supp fig s2

library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr) 

ici<-read.table("/Users/tjsears/Code/GermTime/RPlottingScripts_final/Data/Fig3/all_prs.txt",sep="\t",header=T)
df_all<-read.table("/Users/tjsears/Code/GermTime/RPlottingScripts_final/Data/Fig2/combined_all_pts.txt",sep="\t",header=T)

# somatic heatmap
  
df<-read.table("/Users/tjsears/Code/GermTime/RPlottingScripts_final/Data/suppFigs/somatic_all_pts.txt",sep="\t",header=T)
df_all_gender<-df_all[,c("normal.WXS.id","Gender")]
df<-merge(df,df_all_gender,by= "normal.WXS.id",how="left")

df_temp<-merge(ici,df[,c("normal.WXS.id","TMB","zTMB","cTMB","pTMB")],by="normal.WXS.id")

features<-c("zTMB","ClonalTMB","FractionSubclonal","immune_DNDS","T.cell.frac","Age","numLostAlleles","Escape","MHC_Class1",
            "DNA_Repair","NumLargeSubClones","X9pLoss","PT_pathway")

#df$study_cancer[df$study_cancer=="cristescu_melanoma"]<-"cristescu"
#df$study_cancer[df$study_cancer=="cristescu_urothelial"]<-"cristescu"
#df$study_cancer[df$study_cancer=="cristescu_hnscc"]<-"cristescu"

# Initialize empty vectors to store coefficients and p-values

final_results<-as.data.frame(matrix(ncol=4))
colnames(final_results)<-c("Predictor","Coefficient","P_Value","Study")

studies<-unique(df$study_cancer)

for (study in studies){
  
  df_temp<-df[df$study_cancer==study,]
  coefficients <- numeric()
  p_values <- numeric()
  
  for (col_name in features) {
    lm_result <- glm(pheno ~ scale(df_temp[, col_name])+Age+Gender, data = df_temp)
    coefficients <- c(coefficients, coef(lm_result)[2])  # Store the coefficient
    p_values <- c(p_values, summary(lm_result)$coef[2, 4])  # Store the p-value
  }
  
  study_name<-rep(study,length(features))
  # Create a dataframe to store the results
  results_df <- data.frame(Predictor = features, 
                           Coefficient = coefficients,
                           P_Value = p_values,
                           Study=study_name)
  
  final_results<-rbind(final_results,results_df)
}

final_results<-final_results[-1,]
library(ComplexHeatmap)

final_results_wide<-pivot_wider(final_results,names_from = c("Study"),values_from = c("Coefficient","P_Value"))
final_results_wide_coef<-final_results_wide[,1:8]

final_results$Predictor<-factor(final_results$Predictor,levels=rev(features))
final_results$Study<-factor(final_results$Study,levels=c("hugo","riaz","snyder","vanallen","miao","rizvi","cristescu"))

final_results$logPval<-abs(log2(final_results$P_Value))

#remake as a heatmap with marked categories for fig 2 V2.

#I think I will remove SHAP plot too?

#annotate with pvals and groupings--probably remove nonsig stuff? or at least the checkpoints.
library(RColorBrewer)
library(ComplexHeatmap)
cols<-brewer.pal(n=3,"Dark2")
cols<-cols[c(2,3,1)]
names(cols)<-c("MHC-I Reliant","Balanced","MHC-II Reliant")

ha = HeatmapAnnotation(Group =c("MHC-I Reliant","Balanced","MHC-II Reliant"),col=list(Group=cols),
                       annotation_legend_param = list(NULL))

plot_df_heatmap<-final_results_wide
#plot_df_heatmap<-plot_df_heatmap[match(c("Sigs160.Bindea_Tfh.cells","FractionSubclonal","NumLargeSubClones",
#                                         "numLostAlleles","Escape","MHC_Class1","immune_DNDS",
#                                         "PDCD1","GPLD1","CTSS","TREX1","FPR1","ITGB2","CTSW","FCGR3B","FCGR2B","FAM167A",
#                                         "VAMP8","ERAP2","DCTN5","VAMP3","DHFR","LYZ","ERAP1"),plot_df_heatmap$Predictor),]
breaks<-factor(c("Immune\nInfiltration",rep("Immunogenicity",2),rep("Immune\nEvasion",4),rep("Immune\nSignaling",10),
                 rep("Antigen\nPresentation",7)),levels=c("Immune\nInfiltration","Immunogenicity","Immune\nEvasion","Immune\nSignaling","Antigen\nPresentation"))

plot_df_heatmap$Predictor[plot_df_heatmap$Predictor=="NumLargeSubClones"]<-"Intratumoral Heterogeneity"
plot_df_heatmap$Predictor[plot_df_heatmap$Predictor=="Escape"]<-"Immune Escape"
plot_df_heatmap$Predictor[plot_df_heatmap$Predictor=="MHC_Class1"]<-"Antigen Presentation Pathway"
#plot_df_heatmap$Predictor[plot_df_heatmap$Predictor=="Sigs160.Bindea_Tfh.cells"]<-"TFH Infiltration"
plot_df_heatmap$Predictor[plot_df_heatmap$Predictor=="numLostAlleles"]<-"Class-I HLA Damage"
plot_df_heatmap$Predictor[plot_df_heatmap$Predictor=="immune_DNDS"]<-"ImmunoEditing"
plot_df_heatmap$Predictor[plot_df_heatmap$Predictor=="FractionSubclonal"]<-"Fraction of TMB Subclonal"
plot_df_heatmap$Predictor[plot_df_heatmap$Predictor=="FractionClonal"]<-"Fraction of TMB Clonal"

plot_df_heatmap_coef<-data.frame(plot_df_heatmap[,c(2:10)],row.names=plot_df_heatmap$Predictor)
plot_df_heatmap_pval<-data.frame(plot_df_heatmap[,c(11:19)],row.names=plot_df_heatmap$Predictor)

colnames(plot_df_heatmap_coef)<-str_to_title(gsub("Coefficient_","",colnames(plot_df_heatmap_coef)))
#order by total coef? or average coef?
#lets split by effect honestly... also time to rebrand some of these variable names

# categories:
# immune escape
# antigen presentation 
# antigen processing
# immune signaling
# immune infiltration

#make annotation for each
#make mapping string? so we can use it for other plots too?

#plot_df_heatmap_coef<-plot_df_heatmap_coef[order(rowMeans(plot_df_heatmap_pval)),] 
#plot_df_heatmap_pval<-plot_df_heatmap_pval[order(rowMeans(plot_df_heatmap_pval)),]

#notes so far:
# replace Pvals with shapes (for now doing X shape... looks ok?)
# add annotation colors and a key
# remove x axis text
# rename several of the variables
# also need tissue type annotations
# maybe do green yellow assoc colors--green=good after all

ha = HeatmapAnnotation(Tissue = c("Melanoma","NSCLC","RCC","Urothelial","HNSCC","Melanoma","Melanoma","Melanoma","Melanoma"),
                       col = list(Tissue = c("Melanoma" = "brown", "NSCLC" = "#9CADCE", "RCC" = "#FFCA3A",
                                             "Urothelial"="skyblue","HNSCC"="grey15")),gp = gpar(col = "black"))

ha2 = rowAnnotation(Source = c(rep("Germline",1),rep("Somatic",6),rep("Germline",17)),
                    col = list(Source = c("Germline" = "#52B2CF", "Somatic"="#CE7DA5")))

library(circlize)
col_fun = colorRamp2(c(-0.15, 0, 0.15), c("#F2C14E", "whitesmoke", "#4D9078"))

Heatmap(plot_df_heatmap_coef,cluster_rows = T,cluster_columns=T,show_row_dend = F,show_column_dend = F,#left_annotation = ha2,
        row_title_rot = 0,top_annotation = ha,rect_gp = gpar(col = "white", lwd = 1.5),
        heatmap_legend_param = list(title = "ICI Assoc."),col=col_fun,
        cell_fun = function(j, i, x, y, width, height, fill) {
          #if(plot_df_heatmap_pval[i, j] <=0.1 & plot_df_heatmap_pval[i, j] >0.05)
            #  grid.circle(x = x, y = y, r = 0.3*height,default.units="npc", 
            #              gp = gpar(fill = "purple",size=1))
            #grid.points(x=x,y=y,pch=1)
          if(plot_df_heatmap_pval[i, j] <=0.05)
            #grid.rect(x = x, y = y, r = 0.5*height,default.units="npc", 
            #            gp = gpar(fill = "yellow",size=1))
            #grid.text("x", x, y,vjust=0.3,gp=gpar(lwd=5,fill="black"))
            grid.points(x=x,y=y,pch=4)
          #if(plot_df_heatmap_pval[i, j] <=0.01)
          #  grid.circle(x = x, y = y, r = 0.5*height,default.units="npc", 
          #              gp = gpar(fill = "orange",size=1))
        }) 

#have to save heatmap manually
#ggsave("/Users/tjsears/Code/GermTime/figs_july/figS2/SomaticHeatmap.pdf",width = 8,height = 6)


# show shap feature importance for all models?

  # they are listed in the ppt.
  

# correlation of TMB with frac subclonal and ITH

df$Response<-as.factor(ifelse(df$pheno==2,"R","NR"))

df_plot<-df
df$NumLargeSubClones<-round(df$NumLargeSubClones)
ggplot(df, aes(x = TMB, y = NumLargeSubClones,fill=Response)) +
  geom_jitter(width = 0, height = 0,alpha=0.7,size=2)+
  stat_cor(method = "pearson",show.legend = F,inherit.aes = F,aes(x = TMB, y = NumLargeSubClones),size=6)+
  #geom_smooth(method = "lm",inherit.aes = F,aes(x = TMB, y = NumLargeSubClones)) +
  theme_light(base_size=16) + scale_color_manual(values=cols[2:1]) #+ xlab(label = "Germline IC-Index") + ylab(label = "Somatic IC-Index")
#ggsave("/Users/tjsears/Code/GermTime/figs_july/figS2/ITH_TMB_corr.pdf",width = 8,height = 6)

  
df_plot<-df
df$FractionSubclonal[df$FractionSubclonal>1]<-0.9 # max it out at 0.9 for plotting

ggplot(df, aes(x = TMB, y = FractionSubclonal,fill=Response)) +
  geom_jitter(width = 0, height = 0,alpha=0.7,size=2)+
  stat_cor(method = "pearson",show.legend = F,inherit.aes = F,aes(x = TMB, y = FractionSubclonal),size=6)+
  #geom_smooth(method = "lm",inherit.aes = F,aes(x = TMB, y = NumLargeSubClones)) +
  theme_light(base_size=16) + scale_color_manual(values=cols[2:1]) #+ xlab(label = "Germline IC-Index") + ylab(label = "Somatic IC-Index")
#ggsave("/Users/tjsears/Code/GermTime/figs_july/figS2/FractionSubclonal_TMB_corr.pdf",width = 8,height = 6)


#make germline heatmap

df_all<-read.table("Data/Fig2/combined_all_pts.txt",sep="\t",header=T)

# somatic heatmap

df<-read.table("Data/suppFigs/germline_all_pts.txt",sep="\t",header=T)
df_all_gender<-df_all[,c("normal.WXS.id","Gender","Age_y")]
df<-merge(df,df_all_gender,by.y= "normal.WXS.id",by.x="FID",how="left")

features<-colnames(df)[2:24]

#df$study_cancer[df$study_cancer=="cristescu_melanoma"]<-"cristescu"
#df$study_cancer[df$study_cancer=="cristescu_urothelial"]<-"cristescu"
#df$study_cancer[df$study_cancer=="cristescu_hnscc"]<-"cristescu"

# Initialize empty vectors to store coefficients and p-values

final_results<-as.data.frame(matrix(ncol=4))
colnames(final_results)<-c("Predictor","Coefficient","P_Value","Study")

studies<-unique(df$study)

for (study in studies){
  
  df_temp<-df[df$study==study,]
  coefficients <- numeric()
  p_values <- numeric()
  
  for (col_name in features) {
    lm_result <- glm(pheno ~ scale(df_temp[, col_name])+Age_y+Gender, data = df_temp)
    coefficients <- c(coefficients, coef(lm_result)[2])  # Store the coefficient
    p_values <- c(p_values, summary(lm_result)$coef[2, 4])  # Store the p-value
  }
  
  study_name<-rep(study,length(features))
  # Create a dataframe to store the results
  results_df <- data.frame(Predictor = features, 
                           Coefficient = coefficients,
                           P_Value = p_values,
                           Study=study_name)
  
  final_results<-rbind(final_results,results_df)
}


final_results<-final_results[-1,]
library(ComplexHeatmap)

final_results_wide<-pivot_wider(final_results,names_from = c("Study"),values_from = c("Coefficient","P_Value"))
final_results_wide_coef<-final_results_wide[,1:8]

final_results$Predictor<-factor(final_results$Predictor,levels=rev(features))
final_results$Study<-factor(final_results$Study,levels=c("hugo","riaz","snyder","vanallen","miao","rizvi","cristescu"))

final_results$logPval<-abs(log2(final_results$P_Value))


#remake as a heatmap with marked categories for fig 2 V2.

#I think I will remove SHAP plot too?

#annotate with pvals and groupings--probably remove nonsig stuff? or at least the checkpoints.
library(RColorBrewer)
library(ComplexHeatmap)
cols<-brewer.pal(n=3,"Dark2")
cols<-cols[c(2,3,1)]
names(cols)<-c("MHC-I Reliant","Balanced","MHC-II Reliant")

ha = HeatmapAnnotation(Group =c("MHC-I Reliant","Balanced","MHC-II Reliant"),col=list(Group=cols),
                       annotation_legend_param = list(NULL))

plot_df_heatmap<-final_results_wide

plot_df_heatmap_coef<-data.frame(plot_df_heatmap[,c(2:10)],row.names=plot_df_heatmap$Predictor)
plot_df_heatmap_pval<-data.frame(plot_df_heatmap[,c(11:19)],row.names=plot_df_heatmap$Predictor)

colnames(plot_df_heatmap_coef)<-str_to_title(gsub("Coefficient_","",colnames(plot_df_heatmap_coef)))


ha = HeatmapAnnotation(Tissue = c("Melanoma","Melanoma","Melanoma","Melanoma","Melanoma","RCC","Urothelial","HNSCC","Melanoma"),
                       col = list(Tissue = c("Melanoma" = "brown", "NSCLC" = "#9CADCE", "RCC" = "#FFCA3A",
                                             "Urothelial"="skyblue","HNSCC"="grey15")),gp = gpar(col = "black"))

ha2 = rowAnnotation(Source = c(rep("Germline",1),rep("Somatic",6),rep("Germline",17)),
                    col = list(Source = c("Germline" = "#52B2CF", "Somatic"="#CE7DA5")))

library(circlize)
col_fun = colorRamp2(c(-0.15, 0, 0.15), c("#F2C14E", "whitesmoke", "#4D9078"))

Heatmap(plot_df_heatmap_coef,cluster_rows = T,cluster_columns=F,show_row_dend = F,show_column_dend = F,#left_annotation = ha2,
        row_title_rot = 0,top_annotation = ha,rect_gp = gpar(col = "white", lwd = 1.5),
        heatmap_legend_param = list(title = "ICI Assoc."),col=col_fun,
        cell_fun = function(j, i, x, y, width, height, fill) {
          #if(plot_df_heatmap_pval[i, j] <=0.1 & plot_df_heatmap_pval[i, j] >0.05)
          #  grid.circle(x = x, y = y, r = 0.3*height,default.units="npc", 
          #              gp = gpar(fill = "purple",size=1))
          #grid.points(x=x,y=y,pch=1)
          if(plot_df_heatmap_pval[i, j] <=0.05)
            #grid.rect(x = x, y = y, r = 0.5*height,default.units="npc", 
            #            gp = gpar(fill = "yellow",size=1))
            #grid.text("x", x, y,vjust=0.3,gp=gpar(lwd=5,fill="black"))
            grid.points(x=x,y=y,pch=4)
          #if(plot_df_heatmap_pval[i, j] <=0.01)
          #  grid.circle(x = x, y = y, r = 0.5*height,default.units="npc", 
          #              gp = gpar(fill = "orange",size=1))
        }) 













