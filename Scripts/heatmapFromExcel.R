library(pheatmap)
library(readxl)
library(viridis)
library(ggplot2)
library(tidyverse)

dea_counts_matrix <- as.data.frame(read_xlsx('data/DEA_SNCA_vs_WT_Giosue.xlsx'))

rownames(dea_counts_matrix) <- dea_counts_matrix$...1

dea_counts_matrix$...1 <- NULL


features_of_interest <- c('Gfap','Vim','Serping1','Ggta1','Iigp1','Gbp2','Fbln5','Fkbp5','Psmb8','Srgn')
features_of_interest <- c('Ccl8','Ccl4','Ccl7','Ccl5','Cxcl10','Cxcl13','Cxcl16','Ccl3','Ccl2',
                          'Cx3cl1','Ccl1')


dea_counts_subset <- dea_counts_matrix[features_of_interest,]

ComplexHeatmap::pheatmap(as.matrix(dea_counts_subset),
         scale = "row",
         #color = rev(RColorBrewer::brewer.pal(n =4, name = "RdBu")),
         color = inferno(256),
         #color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         main = "Heatmap of Gene Expression",
         cluster_rows=FALSE, cluster_cols=FALSE,
         column_names_side=c('top'), angle_col = c('0'),
         row_names_side=c('left'),
         heatmap_legend_param = list(title = NULL)
         )

# Heatmap with ggplot
dea_counts_subset_scaled <- t(scale(t(dea_counts_subset)))
dea_counts_subset_modified <- data.frame(Gene=rownames(dea_counts_subset_scaled),dea_counts_subset_scaled)
dea_counts_subset_long <- pivot_longer(dea_counts_subset_modified, cols = -Gene)
dea_counts_subset_long$Gene <- factor(dea_counts_subset_long$Gene, levels = rev(features_of_interest))  # set genes as factor to order the rows

ggplot(dea_counts_subset_long, aes(x=name, y=Gene, fill=value)) +
  geom_tile(width=0.99, height=0.99) +
  scale_fill_viridis(option='rocket', discrete=F) +
  #scale_fill_distiller(palette = "PuOr") +
  theme_minimal() +
  theme(text=element_text(size=20), legend.title=element_blank(), legend.key.height = unit(2, "cm")) +
  scale_x_discrete(position = "top") +
  labs(x=NULL,y=NULL)



