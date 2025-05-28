library(openxlsx)
library(pheatmap)

dea_counts_matrix <- as.data.frame(read_xlsx('data/DEA_SNCA_vs_WT_Giosue.xlsx'))

rownames(dea_counts_matrix) <- dea_counts_matrix$...1

dea_counts_matrix$...1 <- NULL


features_of_interest <- c('Gfap','Vim','Serping1','Ggta1','Iigp1','Gbp2','Fbln5','Fkbp5','Psmb8','Srgn')

dea_counts_subset <- dea_counts_matrix[features_of_interest,]

ComplexHeatmap::pheatmap(as.matrix(dea_counts_subset),
         scale = "row",
         color = rev(RColorBrewer::brewer.pal(n =4, name = "RdBu")),
         #color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         main = "Heatmap of Gene Expression",
         cluster_rows=FALSE, cluster_cols=FALSE,
         column_names_side=c('top'), angle_col = c('45'),
         row_names_side=c('left'),
         )

