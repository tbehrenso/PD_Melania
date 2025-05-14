library(viridis)

V4.2_MYELOID <- c(0,2,5,10,13,16,17,18,22,27,29,38)
V4.2_TCELL <- c(1,6,7,8,9,11,12,24,42)
V4.2_BCELL <- c(3,4,19,28,30,40)


# ----------------------------------------------------------------------------------------
# Quality Control
VlnPlot(SYN_merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", 'percent.ribo'), ncol = 2)
FeatureScatter(SYN_merged, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")


VizDimLoadings(SYN_merged, dims = 15:20, reduction = "pca")
DimPlot(SYN_merged, reduction = "pca") + NoLegend()

# ----------------------------------------------------------------------------------------
# UMAP

DimPlot(SYN_merged, reduction = 'umap', group.by = c('seurat_clusters'), label = T)

DimPlot(SYN_merged, reduction = 'umap', group.by = c('orig.ident', 'seurat_clusters'))

DimPlot(SYN_merged, reduction = 'umap', split.by = 'orig.ident')

# ----------------------------------------------------------------------------------------
# Investigate PCA
DimHeatmap(SYN_merged, dims = 1:19, cells = 500, balanced = TRUE)

ElbowPlot(SYN_merged, ndims=50)

# ----------------------------------------------------------------------------------------
# Violin Plot
VlnPlot(SYN_merged, features = "Olfml3", split.by = "orig.ident", group.by = "seurat_clusters") 
#  + geom_boxplot(width=0.1, fill="white")


# ----------------------------------------------------------------------------------------
# Various Gene Visualizations

# DotPlot
DotPlot(SYN_merged, features = c('Foxp3','Cd4','Il2ra'),idents=c(1,3,5,9,10,12,24,31)) + RotatedAxis() #Treg
DotPlot(SYN_merged, features = c('Cd44','Sell','Il7r','Cd69','Itgae','Ctla4'),idents=c(1,3,5,9,10,12,24,31)) + RotatedAxis() # Tmem

# Example split by Condition
DotPlot(SYN_merged, features = c("Cd3d", "Cd3e", "Cd3g","Cd8b1","Cd5","Cd28","Cd8a","Il7r","Sell","Tcf7","Txk","S1pr1","Trat1","Lef1","Satb1"),
        idents=c(1,6,7,8,9,11,12,24,25,28,30,42), split.by='orig.ident', cols=c('red','blue')) + 
  RotatedAxis() + ggtitle('CD8 - Naive-like')

# from gene list PowerPoint
genes_vTom <- read_xlsx('/Users/tbehr/Desktop/genes_vTom.xlsx')

DotPlot(SYN_merged, features = hierarchical_clustering_scaled(genes_vTom$M1[!is.na(genes_vTom$M1)], V4.2_MYELOID),
        idents=V4.2_MYELOID) + RotatedAxis() + ggtitle('M1')
DotPlot(SYN_merged, features = genes_vTom$M2[!is.na(genes_vTom$M2)],
        idents=V4.2_MYELOID) + RotatedAxis() + ggtitle('M2')
DotPlot(SYN_merged, features = genes_vTom$`inferferon related (M1?)`[!is.na(genes_vTom$`inferferon related (M1?)`)],
        idents=V4.2_MYELOID) + RotatedAxis() + ggtitle('Interferon Related (M1?)')
DotPlot(SYN_merged, features = genes_vTom$migration[!is.na(genes_vTom$migration)],
        idents=V4.2_MYELOID) + RotatedAxis() + ggtitle('Migration')
DotPlot(SYN_merged, features = genes_vTom$`reactive microglia`[!is.na(genes_vTom$`reactive microglia`)],
        idents=V4.2_MYELOID) + RotatedAxis() + ggtitle('Reactive Microglia')
DotPlot(SYN_merged, features = genes_vTom$`infiltrating macrophages`[!is.na(genes_vTom$`infiltrating macrophages`)],
        idents=V4.2_MYELOID) + RotatedAxis() + ggtitle('Infiltrating Microglia')

DotPlot(SYN_merged, features = genes_vTom$Th1[!is.na(genes_vTom$Th1)],
        idents=V4.2_TCELL) + RotatedAxis() + ggtitle('Th1')
DotPlot(SYN_merged, features = genes_vTom$Th2[!is.na(genes_vTom$Th2)],
        idents=V4.2_TCELL) + RotatedAxis() + ggtitle('Th2')
DotPlot(SYN_merged, features = genes_vTom$Th17[!is.na(genes_vTom$Th17)],
        idents=V4.2_TCELL) + RotatedAxis() + ggtitle('Th17')
DotPlot(SYN_merged, features = genes_vTom$Tfh[!is.na(genes_vTom$Tfh)],
        idents=V4.2_TCELL) + RotatedAxis() + ggtitle('Tfh')
DotPlot(SYN_merged, features = genes_vTom$Treg[!is.na(genes_vTom$Treg)],
        idents=V4.2_TCELL) + RotatedAxis() + ggtitle('Treg')
DotPlot(SYN_merged, features = genes_vTom$`Tr1 (these loose foxp3)`[!is.na(genes_vTom$`Tr1 (these loose foxp3)`)],
        idents=V4.2_TCELL) + RotatedAxis() + ggtitle('Tr1')
DotPlot(SYN_merged, features = genes_vTom$NKT[!is.na(genes_vTom$NKT)],
        idents=V4.2_TCELL) + RotatedAxis() + ggtitle('NKT')
DotPlot(SYN_merged, features = genes_vTom$`T cells naive`[!is.na(genes_vTom$`T cells naive`)],
        idents=V4.2_TCELL) + RotatedAxis() + ggtitle('T Cell Naive')
DotPlot(SYN_merged, features = genes_vTom$`T cells exhausted`[!is.na(genes_vTom$`T cells exhausted`)],
        idents=V4.2_TCELL) + RotatedAxis() + ggtitle('T Cell Exhausted')

DotPlot(SYN_merged, features = genes_vTom$Plasmacells[!is.na(genes_vTom$Plasmacells)],
        idents=V4.2_BCELL) + RotatedAxis() + ggtitle('Plasma Cells')
DotPlot(SYN_merged, features = genes_vTom$`B cells naive`[!is.na(genes_vTom$`B cells naive`)],
        idents=V4.2_BCELL) + RotatedAxis() + ggtitle('B Cells Naive')
DotPlot(SYN_merged, features = genes_vTom$`B cells memory`[!is.na(genes_vTom$`B cells memory`)],
        idents=V4.2_BCELL) + RotatedAxis() + ggtitle('B Cells Memory')

# Correlation Plot
FeatureScatter(SYN_merged_relabeled, feature1='Cd8a', feature2='Saa3')

# Feature on UMAP Plot
FeaturePlot(SYN_merged, features='Cd8a')
# w/ Coexpression
FeaturePlot(SYN_merged, features = c('Itgam','Itgae'), blend = TRUE)
# Interactive (not working)
InteractivePlot <- FeaturePlot(SYN_merged_relabeled, features = 'C1qa')
HoverLocator(plot = InteractivePlot, information = FetchData(SYN_merged_relabeled, vars = c('orig.ident', 'PC_1', 'nFeature_RNA')))

# Heatmap
DoHeatmap(subset(SYN_merged, downsample=100, idents=c('0','2','6','13','14','17','21','23','25')), 
          features = rownames(head(FindMarkers(SYN_merged, ident.1 = '25', ident.2=c('0','2','6','13','14','17','21','23'), only.pos=T), 30)),
          size = 3, slot='scale.data') + scale_fill_viridis(option='magma') + ggtitle('Cluster 25: Top Markers') # Top Markers for Cluster 25

DoHeatmap(subset(SYN_merged, downsample = 100, idents=c('0_SNCA','0_SNCA_PLX','2_SNCA','2_SNCA_PLX','5_SNCA','5_SNCA_PLX','10_SNCA','10_SNCA_PLX','13_SNCA','13_SNCA_PLX','16_SNCA','16_SNCA_PLX',
                                      '17_SNCA','17_SNCA_PLX','18_SNCA','18_SNCA_PLX','22_SNCA','22_SNCA_PLX','27_SNCA','27_SNCA_PLX','29_SNCA','29_SNCA_PLX','38_SNCA','38_SNCA_PLX'
                                      )), features = c('P2ry12','Tmem119','Sall1','Fcrl2','Gpr34','Hexb','Olfml3','Clec7a','Siglech','Siglece','Ptprc'),
          size = 3, slot='scale.data') + ggtitle('Microglia Markers') + scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n =4, name = "RdBu"))) # Microglia

# Heatmaps from Excel
DoHeatmap(subset(SYN_merged, downsample = 100, idents = V4.2_MYELOID), 
          features = hierarchical_clustering_scaled(genes_vTom$M1[!is.na(genes_vTom$M1)], V4.2_MYELOID),
          size = 3, slot='scale.data') + ggtitle('Macrophage Markers') + scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n =4, name = "RdBu")))



# Heatmap of AVERAGE expression
SYN_averages <- AverageExpression(SYN_merged, return.seurat = T)

Idents(SYN_averages) <- sapply(Idents(SYN_averages), function(x) substr(x, 2, nchar(as.character(x))))

#ident_subset <- c('0_SNCA','0_SNCA_PLX','2_SNCA','2_SNCA_PLX','5_SNCA','5_SNCA_PLX','10_SNCA','10_SNCA_PLX','13_SNCA','13_SNCA_PLX','16_SNCA','16_SNCA_PLX',
#                         '17_SNCA','17_SNCA_PLX','18_SNCA','18_SNCA_PLX','22_SNCA','22_SNCA_PLX','27_SNCA','27_SNCA_PLX','29_SNCA','29_SNCA_PLX','38_SNCA','38_SNCA_PLX')
ident_subset <- V4.2_MYELOID
idents_subset <- gsub('_', '-', ident_subset)
idents_subset <- paste0('g', idents_subset)

features_of_interest <- genes_vTom$`infiltrating macrophages`

DoHeatmap(subset(SYN_averages, idents = V4.2_MYELOID), 
          features = hierarchical_clustering_scaled(features_of_interest[!is.na(features_of_interest)], V4.2_MYELOID), label = T ,draw.lines = F)  +
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n =4, name = "RdBu")))


# Ridgemap
RidgePlot(SYN_merged_relabeled, features = 'C1qa', ncol = 2)
































