# Fig 2. Linear outcome of each feature across all cohorts
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr) 

df<-read.table("/Users/tjsears/Code/GermTime/RPlottingScripts_final/Data/Fig2/combined_all_pts.txt",sep="\t",header=T)

#get list of columns used in the final model

features<-c("Sigs160.Bindea_Tfh.cells","ERAP2","DCTN5","DHFR","GPLD1","FractionSubclonal","immune_DNDS","numLostAlleles",
            "Escape","MHC_Class1","ITGB2","PDCD1","FCGR3B","CTSW","FPR1","VAMP8","FCGR2B","TREX1","CTSS",
            "FAM167A","ERAP1","NumLargeSubClones","LYZ","VAMP3")

df$study_cancer_x[df$study_cancer_x=="cristescu_melanoma"]<-"cristescu"
df$study_cancer_x[df$study_cancer_x=="cristescu_urothelial"]<-"cristescu"
df$study_cancer_x[df$study_cancer_x=="cristescu_hnscc"]<-"cristescu"

# Initialize empty vectors to store coefficients and p-values

final_results<-as.data.frame(matrix(ncol=5))
colnames(final_results)<-c("Predictor","Coefficient","P_Value","T_Value","Study")

studies<-unique(df$study_cancer_x)

for (study in studies){
  
  df_temp<-df[df$study_cancer_x==study,]
  coefficients <- numeric()
  p_values <- numeric()
  t_values <- numeric()
  
  for (col_name in features) {
    lm_result <- glm(pheno ~ scale(df_temp[, col_name])+Age_y+Gender, data = df_temp)
    coefficients <- c(coefficients, coef(lm_result)[2])  # Store the coefficient
    p_values <- c(p_values, summary(lm_result)$coef[2, 4])
    t_values <- c(t_values, summary(lm_result)$coef[2, 3])# Store the p-value
  }
  
  study_name<-rep(study,length(features))
  # Create a dataframe to store the results
  results_df <- data.frame(Predictor = features, 
                           Coefficient = coefficients,
                           P_Value = p_values,
                           T_Value=t_values,
                           Study=study_name)
  
  final_results<-rbind(final_results,results_df)
}

final_results<-final_results[-1,]
library(ComplexHeatmap)

final_results_wide<-pivot_wider(final_results,names_from = c("Study"),values_from = c("Coefficient","P_Value","T_Value"))
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
plot_df_heatmap<-plot_df_heatmap[match(c("Sigs160.Bindea_Tfh.cells","FractionSubclonal","NumLargeSubClones",
                                   "numLostAlleles","Escape","MHC_Class1","immune_DNDS",
                                   "PDCD1","GPLD1","CTSS","FPR1","ITGB2","CTSW","FCGR3B","FCGR2B","FAM167A",
                                   "VAMP8","ERAP2","DCTN5","VAMP3","LYZ","ERAP1","TREX1","DHFR"),plot_df_heatmap$Predictor),]
breaks<-factor(c("Immune\nInfiltration",rep("Immunogenicity",2),rep("Immune\nEvasion",4),rep("Immune\nSignaling",9),
          rep("Antigen\nProcessing &\nPresentation",6),rep("DNA\nRepair &\nReplication",2)),levels=c("Immune\nInfiltration","Immunogenicity","Immune\nEvasion","Immune\nSignaling","DNA\nRepair &\nReplication","Antigen\nProcessing &\nPresentation"))

plot_df_heatmap$Predictor[plot_df_heatmap$Predictor=="NumLargeSubClones"]<-"Intratumoral Heterogeneity"
plot_df_heatmap$Predictor[plot_df_heatmap$Predictor=="Escape"]<-"Immune Escape"
plot_df_heatmap$Predictor[plot_df_heatmap$Predictor=="MHC_Class1"]<-"Antigen Presentation Pathway"
plot_df_heatmap$Predictor[plot_df_heatmap$Predictor=="Sigs160.Bindea_Tfh.cells"]<-"TFH Infiltration"
plot_df_heatmap$Predictor[plot_df_heatmap$Predictor=="numLostAlleles"]<-"Class-I HLA Damage"
plot_df_heatmap$Predictor[plot_df_heatmap$Predictor=="immune_DNDS"]<-"ImmunoEditing"
plot_df_heatmap$Predictor[plot_df_heatmap$Predictor=="FractionSubclonal"]<-"Fraction of TMB Subclonal"

plot_df_heatmap_coef<-data.frame(plot_df_heatmap[,c(2:8)],row.names=plot_df_heatmap$Predictor)
plot_df_heatmap_pval<-data.frame(plot_df_heatmap[,c(9:15)],row.names=plot_df_heatmap$Predictor)

colnames(plot_df_heatmap_coef)<-str_to_title(gsub("Coefficient_","",colnames(plot_df_heatmap_coef)))
#order by total coef? or average coef?
#lets split by effect honestly... 

# categories:
# immune escape
# antigen presentation 
# antigen processing
# immune signaling
# immune infiltration

#make annotation for each
#make mapping string? so we can use it for other plots too?


ha = HeatmapAnnotation(Tissue = c("Melanoma","Melanoma","Melanoma","Melanoma","NSCLC","RCC","Multi-Tissue"),
                       col = list(Tissue = c("Melanoma" = "brown", "NSCLC" = "#9CADCE", "RCC" = "#FFCA3A",
                                          "Multi-Tissue"="grey70")),gp = gpar(col = "black"))

ha2 = rowAnnotation(Source = c(rep("Germline",1),rep("Somatic",6),rep("Germline",17)),
                       col = list(Source = c("Germline" = "#52B2CF", "Somatic"="#CE7DA5")))


library(circlize)
col_fun = colorRamp2(c(-0.15, 0, 0.15), c("#F2C14E", "whitesmoke", "#4D9078"))

Heatmap(plot_df_heatmap_coef,cluster_rows = T,cluster_columns=T,show_row_dend = F,show_column_dend = F,left_annotation = ha2,
        row_split=breaks,row_title_rot = 0,top_annotation = ha,rect_gp = gpar(col = "white", lwd = 1.5),
        heatmap_legend_param = list(title = "ICB Assoc."),col=col_fun,
        cell_fun = function(j, i, x, y, width, height, fill) {
          if(plot_df_heatmap_pval[i, j] <=0.1 & plot_df_heatmap_pval[i, j] >0.05)
          #  grid.circle(x = x, y = y, r = 0.3*height,default.units="npc", 
          #              gp = gpar(fill = "purple",size=1))
            grid.points(x=x,y=y,pch=1)
          if(plot_df_heatmap_pval[i, j] <=0.05)
            #grid.rect(x = x, y = y, r = 0.5*height,default.units="npc", 
            #            gp = gpar(fill = "yellow",size=1))
            #grid.text("x", x, y,vjust=0.3,gp=gpar(lwd=5,fill="black"))
            grid.points(x=x,y=y,pch=4)
          #if(plot_df_heatmap_pval[i, j] <=0.01)
          #  grid.circle(x = x, y = y, r = 0.5*height,default.units="npc", 
          #              gp = gpar(fill = "orange",size=1))
        }) 
#ggsave("/Users/tjsears/Code/GermTime/figs_july/fig2/TestHeatmap.pdf",width = 8,height = 6)







# feature importance change

# this should be cool, just do connected dot plot thing...?
pval_means<-rowMeans(plot_df_heatmap_pval)
pval_means<-pval_means[order(pval_means)]
pval_means<-data.frame(pval_means,Arank=c(24:1))

nonlinear_fi<-data.frame(fi=c("TFH Infiltration","ERAP2","DCTN5","DHFR","GPLD1","Fraction of TMB Subclonal","ImmunoEditing","Class-I HLA Damage","Immune Escape",
                              "Antigen Presentation Pathway","ITGB2","PDCD1","FCGR3B","CTSW","FPR1","VAMP8","FCGR2B",
                              "TREX1","CTSS","FAM167A","ERAP1","Intratumoral Heterogeneity","LYZ","VAMP3"),Brank=c(24:1))

linear_fi<-data.frame(fi=c(names(pval_means)),pval_means)
merged_fi<-merge(nonlinear_fi,linear_fi,by.x="fi",by.y=0)
merged_fi<-merged_fi[,c("fi","Arank","Brank")]

data_long <- tidyr::gather(merged_fi, key = "Variable", value = "Value", -fi)
data_long$x1 = c(rep(1 + 0.2, 24), rep(2 - 0.2, 24))
data_long$Value<-as.factor(data_long$Value)
data_long$VariableInt<-ifelse(data_long$Variable=="Brank",2,1)
data_long$Value<-data_long$Value

#data_long<-data_long[nrow(data_long):1,]
# Create the dot plot with lines connecting the dots
cols<-brewer.pal(12,"Set3")

ggplot(data_long, aes(x=VariableInt, y=Value, group=fi, label=fi)) +
  geom_path(aes(x=x1), arrow = arrow(length = unit(0.02,"npc")), 
            size=1, color=c(cols,cols,cols,cols)) +
  geom_label(size=3) + xlim(c(0.8,2.2))+
  theme_minimal() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank())

#ggsave("/Users/tjsears/Code/GermTime/figs_july/figS2/LtoNL_featureimportance.pdf",width = 6,height = 8)


##########################################################################
# PCA of cohort similarity based on composite associations with response #
##########################################################################


final_results_wide

pca_res <- prcomp(t(final_results_wide[,c("Coefficient_hugo"  ,    "Coefficient_snyder"   
                                        , "Coefficient_riaz"   ,   "Coefficient_vanallen"  ,"Coefficient_rizvi"  
                                        , "Coefficient_miao"   ,   "Coefficient_cristescu")]), scale. = F)
pca_res <- prcomp(t(final_results_wide[,c("P_Value_hugo"  ,    "P_Value_snyder"   
                                          , "P_Value_riaz"   ,   "P_Value_vanallen"  ,"P_Value_rizvi"  
                                          , "P_Value_miao"   ,   "P_Value_cristescu")]), scale. = F)
pca_res <- prcomp(t(final_results_wide[,c("T_Value_hugo"  ,    "T_Value_snyder"   
                                          , "T_Value_riaz"   ,   "T_Value_vanallen"  ,"T_Value_rizvi"  
                                          , "T_Value_miao"   ,   "T_Value_cristescu")]), scale. = F)

plot_df<-as.data.frame(pca_res$x)
plot_df$Cohort<-c("Hugo","Snyder","Riaz","Vanallen","Rizvi","Miao","Cristescu")
plot_df$Tissue<-c("Melanoma","Melanoma","Melanoma","Melanoma","NSCLC","RCC","PanCan")

library(paletteer)
cols<-brewer.pal("Dark2",n=7)
cols<-paletteer_d("ggsci::lanonc_lancet")

ggplot(plot_df, aes(x = PC1, y = PC3,color=Cohort)) +  scale_shape_manual(values=c(18,15,19,17)) +
  geom_jitter(width = 0, height = 0,alpha=0.9,size=6,aes(shape=Tissue)) + xlab(label = "PC1") +  ylab(label = "PC2")+
  #stat_cor(method = "pearson",show.legend = F,inherit.aes = F,aes(x = Composite_PRS, y = Germline_PRS),size=6) +
  #geom_smooth(method = "lm",inherit.aes = F,aes(x = Composite_PRS, y = Germline_PRS),color = cols[3]) +
  theme_light(base_size=20) + xlim(-3,7) + ylim(-3,3)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank() )+ scale_color_manual(values=cols) 

#ggsave("/Users/tjsears/Code/GermTime/figs_july/fig2/PCAplot.png",width = 7,height = 5)
















