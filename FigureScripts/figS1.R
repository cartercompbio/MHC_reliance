# Fig S1


#generated from private germline data

# TCGA tfh cell correlation:
# corr,pval 0.05970531607024294, 0.015160998746864598

# ICI tfh cell correlation

# corr,pval 0.18413946742108553, 0.027699037274271275



# Read in corr charts from each plot

ICI_corr<-read.table("Data/suppFigs/ICI_corr_chart.txt",sep="\t",header=T)
ICI_corr<-ICI_corr[86:107,] #front is just repeats of TCGA oops
ICI_corr<-rbind(ICI_corr,c(22,0.18413946742108553, 0.027699037274271275,"TFH SNP/Infiltration"))

tcga_corr<-read.table("Data/suppFigs/ICI_corr_chart.txt",sep="\t",header=T)
tcga_corr<-tcga_corr[1:22,] #front is just repeats of TCGA oops
tcga_corr<-rbind(tcga_corr,c(22,0.05970531607024294, 0.015160998746864598,"TFH SNP/Infiltration"))

corr<-cbind(ICI_corr,tcga_corr)
plot_df_heatmap_coef<-as.data.frame(corr[,c(6,2)])
plot_df_heatmap_coef$corr<-as.numeric(plot_df_heatmap_coef$corr)
plot_df_heatmap_coef$corr.1<-as.numeric(plot_df_heatmap_coef$corr.1)
colnames(plot_df_heatmap_coef)<-c("TCGA","Discovery")
rownames(plot_df_heatmap_coef)<-corr$SNP_name

plot_df_heatmap_pval<-corr[,c(7,3)]
plot_df_heatmap_pval$pval<-as.numeric(plot_df_heatmap_pval$pval)
plot_df_heatmap_pval$pval.1<-as.numeric(plot_df_heatmap_pval$pval.1)

  
library(circlize)
col_fun = colorRamp2(c(-0.15, 0, 0.15), c("#F2C14E", "whitesmoke", "#4D9078"))

library(ComplexHeatmap)
Heatmap(plot_df_heatmap_coef,cluster_rows = T,cluster_columns=F,show_row_dend = F,show_column_dend = F,#left_annotation = ha2,
          row_title_rot = 0,rect_gp = gpar(col = "white", lwd = 1.5),
          heatmap_legend_param = list(title = "SNP Score Corr."),col=col_fun,
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
