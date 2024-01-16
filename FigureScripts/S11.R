#S11
######################################
#### CIBERSORT ASSOCIATIONS STUFF ####
######################################

#read in cibersort and merge with each mhc dataframe
library(ggplot2)
library(RColorBrewer)

#maybe look into using other cibersort file? the non bc one?
cols<-brewer.pal(n = 3, name = "Dark2")[1:2]

tcga_phbr<-read.table("Data/Fig6/fig6v4_liu_rnaseq.txt",sep="\t",header=T)
quantiles <- quantile(tcga_phbr$Ratio, probs = c(0,0.33,0.66,1))
tcga_phbr$group <- cut(tcga_phbr$Ratio, breaks = quantiles, labels = c("MHC-I Reliant","Balanced","MHC-II Reliant"))
tcga_phbr$group[is.na(tcga_phbr$group)]<-"MHC-I Reliant"
tcga_phbr$Response<-ifelse(tcga_phbr$response_crist_sd==2,"R","NR")

ICI_phbr<-read.table("Data/Fig5/fig6v4_ici_rnaseq.txt",sep="\t",header=T)
quantiles <- quantile(ICI_phbr$Ratio, probs = c(0,0.33,0.66,1))
ICI_phbr$group <- cut(ICI_phbr$Ratio, breaks = quantiles, labels = c("MHC-I Reliant","Balanced","MHC-II Reliant"))
ICI_phbr$group[is.na(ICI_phbr$group)]<-"MHC-I Reliant"
ICI_phbr$Response<-ifelse(ICI_phbr$response_crist_sd==2,"R","NR")

cibersort_cols<-c("B.cells.naive","B.cells.memory","Plasma.cells_y","T.cells.CD8","T.cells.CD4.naive",
                  "T.cells.CD4.memory.resting","T.cells.CD4.memory.activated","T.cells.follicular.helper",
                  "T.cells.regulatory..Tregs.","T.cells.gamma.delta","NK.cells.resting","NK.cells.activated","Monocytes_y","Macrophages.M0",
                  "Macrophages.M1_y","Macrophages.M2_y","Dendritic.cells.resting","Dendritic.cells.activated","Mast.cells.resting",
                  "Mast.cells.activated","Eosinophils_y","Neutrophils","group","Response")

cols<-brewer.pal(n=3,"Pastel2")
cols2<-brewer.pal(n=3,"Dark2")#[c(2,3,1)]

library(ComplexHeatmap)

# first, MHC-I R vs NR
plot_df<-ICI_phbr[,cibersort_cols]
plot_df<-plot_df[plot_df$group=="MHC-I Reliant",]
colnames(plot_df)<-gsub("_y","",colnames(plot_df))
colnames(plot_df)<-gsub("[.]"," ",colnames(plot_df))
colnames(plot_df)<-gsub("Tregs","",colnames(plot_df))

names(cols)<-c("MHC-I Reliant")

ha1 = HeatmapAnnotation(group=plot_df$Response,col=list(Group=cols[1:2]),
                        annotation_legend_param = list(NULL))
ha = HeatmapAnnotation(Category=plot_df$group,col=list(Category = c("MHC-I Reliant" = "#D95E01")),
                       annotation_legend_param = list(NULL))

Heatmap(t(plot_df[1:22]),top_annotation = ha,column_split = plot_df$Response,cluster_column_slices = F,cluster_rows = F,
        cluster_columns  =T,show_row_dend = F,show_column_dend = F,show_column_names = F)

# first, Balanced R vs NR
plot_df<-ICI_phbr[,cibersort_cols]
plot_df<-plot_df[plot_df$group=="Balanced",]
colnames(plot_df)<-gsub("_y","",colnames(plot_df))
colnames(plot_df)<-gsub("[.]"," ",colnames(plot_df))
colnames(plot_df)<-gsub("Tregs","",colnames(plot_df))

names(cols)<-c("Balanced")

ha1 = HeatmapAnnotation(group=plot_df$Response,col=list(Group=cols[1:2]),
                        annotation_legend_param = list(NULL))
ha = HeatmapAnnotation(Category=plot_df$group,col=list(Category = c("Balanced" = "#766FB3")),
                       annotation_legend_param = list(NULL))

Heatmap(t(plot_df[1:22]),top_annotation = ha,column_split = plot_df$Response,cluster_column_slices = F,cluster_rows = F,
        cluster_columns  =T,show_row_dend = F,show_column_dend = F,show_column_names = F)

# first, MHC-II Reliant R vs NR
plot_df<-ICI_phbr[,cibersort_cols]
plot_df<-plot_df[plot_df$group=="MHC-II Reliant",]
colnames(plot_df)<-gsub("_y","",colnames(plot_df))
colnames(plot_df)<-gsub("[.]"," ",colnames(plot_df))
colnames(plot_df)<-gsub("Tregs","",colnames(plot_df))

names(cols)<-c("MHC-II Reliant")

ha1 = HeatmapAnnotation(group=plot_df$Response,col=list(Group=cols[1:2]),
                        annotation_legend_param = list(NULL))
ha = HeatmapAnnotation(Category=plot_df$group,col=list(Category = c("MHC-II Reliant" = "#1A9E76")),
                       annotation_legend_param = list(NULL))

Heatmap(t(plot_df[1:22]),top_annotation = ha,column_split = plot_df$Response,cluster_column_slices = F,cluster_rows = F,
        cluster_columns  =T,show_row_dend = F,show_column_dend = F,show_column_names = F)

library(dplyr)

plot_df<-ICI_phbr[,cibersort_cols]
for (celltype in cibersort_cols[1:22]){
  for (group in c("MHC-I Reliant","Balanced","MHC-II Reliant")){
    
    plot_df2<-plot_df[plot_df$group==group,]
    
    x<-plot_df2[plot_df2$Response=="NR",]
    y<-plot_df2[plot_df2$Response=="R",]
    
    #print(mean(x[,celltype]))
    #print(mean(y[,celltype]))
    pval<-t.test(unlist(x[,celltype]),unlist((y[,celltype])))
    if (!is.na(pval$p.value)&pval$p.value<=0.05){
      print(c(celltype,group))
      print(pval)
    }
  }
}


###
# Repeat for validation
###

cibersort_cols<-c("B.cells.naive","B.cells.memory","Plasma.cells_y","T.cells.CD8","T.cells.CD4.naive",
                  "T.cells.CD4.memory.resting","T.cells.CD4.memory.activated","T.cells.follicular.helper",
                  "T.cells.regulatory..Tregs.","T.cells.gamma.delta","NK.cells.resting","NK.cells.activated","Monocytes_y","Macrophages.M0",
                  "Macrophages.M1_y","Macrophages.M2_y","Dendritic.cells.resting","Dendritic.cells.activated","Mast.cells.resting",
                  "Mast.cells.activated","Eosinophils_y","Neutrophils_y","group","Response")

# first, MHC-I R vs NR
plot_df<-tcga_phbr[,cibersort_cols]
plot_df<-plot_df[plot_df$group=="MHC-I Reliant",]
colnames(plot_df)<-gsub("_y","",colnames(plot_df))
colnames(plot_df)<-gsub("[.]"," ",colnames(plot_df))
colnames(plot_df)<-gsub("Tregs","",colnames(plot_df))

names(cols)<-c("MHC-I Reliant")

ha1 = HeatmapAnnotation(group=plot_df$Response,col=list(Group=cols[1:2]),
                        annotation_legend_param = list(NULL))
ha = HeatmapAnnotation(Category=plot_df$group,col=list(Category = c("MHC-I Reliant" = "#D95E01")),
                       annotation_legend_param = list(NULL))

Heatmap(t(plot_df[1:22]),top_annotation = ha,column_split = plot_df$Response,cluster_column_slices = F,cluster_rows = F,
        cluster_columns  =T,show_row_dend = F,show_column_dend = F,show_column_names = F)

# first, Balanced R vs NR
plot_df<-tcga_phbr[,cibersort_cols]
plot_df<-plot_df[plot_df$group=="Balanced",]
colnames(plot_df)<-gsub("_y","",colnames(plot_df))
colnames(plot_df)<-gsub("[.]"," ",colnames(plot_df))
colnames(plot_df)<-gsub("Tregs","",colnames(plot_df))

names(cols)<-c("Balanced")

ha1 = HeatmapAnnotation(group=plot_df$Response,col=list(Group=cols[1:2]),
                        annotation_legend_param = list(NULL))
ha = HeatmapAnnotation(Category=plot_df$group,col=list(Category = c("Balanced" = "#766FB3")),
                       annotation_legend_param = list(NULL))

Heatmap(t(plot_df[1:22]),top_annotation = ha,column_split = plot_df$Response,cluster_column_slices = F,cluster_rows = F,
        cluster_columns  =T,show_row_dend = F,show_column_dend = F,show_column_names = F)

# first, MHC-II Reliant R vs NR
plot_df<-tcga_phbr[,cibersort_cols]
plot_df<-plot_df[plot_df$group=="MHC-II Reliant",]
colnames(plot_df)<-gsub("_y","",colnames(plot_df))
colnames(plot_df)<-gsub("[.]"," ",colnames(plot_df))
colnames(plot_df)<-gsub("Tregs","",colnames(plot_df))

names(cols)<-c("MHC-II Reliant")

ha1 = HeatmapAnnotation(group=plot_df$Response,col=list(Group=cols[1:2]),
                        annotation_legend_param = list(NULL))
ha = HeatmapAnnotation(Category=plot_df$group,col=list(Category = c("MHC-II Reliant" = "#1A9E76")),
                       annotation_legend_param = list(NULL))

Heatmap(t(plot_df[1:22]),top_annotation = ha,column_split = plot_df$Response,cluster_column_slices = F,cluster_rows = F,
        cluster_columns  =T,show_row_dend = F,show_column_dend = F,show_column_names = F)

library(dplyr)

plot_df<-tcga_phbr[,cibersort_cols]
for (celltype in cibersort_cols[1:22]){
  for (group in c("MHC-I Reliant","Balanced","MHC-II Reliant")){
    
    plot_df2<-plot_df[plot_df$group==group,]
    
    x<-plot_df2[plot_df2$Response=="NR",]
    y<-plot_df2[plot_df2$Response=="R",]
    
    #print(mean(x[,celltype]))
    #print(mean(y[,celltype]))
    pval<-t.test(unlist(x[,celltype]),unlist((y[,celltype])))
    if (!is.na(pval$p.value)&pval$p.value<=0.05){
      print(c(celltype,group))
      print(pval)
    }
  }
}

