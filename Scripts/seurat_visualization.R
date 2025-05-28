library(viridis)

V4.2_MYELOID <- c(0,2,5,10,13,16,17,18,22,27,29,38)
V4.2_TCELL <- c(1,6,7,8,9,11,12,24,42)
V4.2_BCELL <- c(3,4,19,28,30,40)
genes_vTom <- read_xlsx('/Users/tbehr/Desktop/genes_vTom.xlsx', sheet = 'Main')

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

DimPlot(SYN_merged_relabeled, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

# ----------------------------------------------------------------------------------------
# Investigate PCA
DimHeatmap(SYN_merged, dims = 1:19, cells = 500, balanced = TRUE)

ElbowPlot(SYN_merged, ndims=50)

# ----------------------------------------------------------------------------------------
# Violin Plot
VlnPlot(SYN_merged, features = "Sema6d", split.by = "orig.ident", group.by = "seurat_clusters")
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

# from gene list Excel
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

DotPlot(SYN_merged, features = genes_vTom$`Phagocytosis and lipid metabolism`[!is.na(genes_vTom$`Phagocytosis and lipid metabolism`)],
        idents=V4.2_MYELOID) + RotatedAxis() + ggtitle('Phagocytosis Markers')

DotPlot(SYN_merged, features = c('Anpep','Spp1','Mrc1','Chil3','Mfge8','Axl','Ccl6','Ccdc80','Angptl4','Serpine1','Mgll'),
        idents=V4.2_MYELOID) + RotatedAxis() + ggtitle('Notable Cluster 29 Genes')

DotPlot(SYN_merged, features = c('Pdcd1','Ctla4','Havcr2','Lag3','Tigit','Cd244a','Cd160','Eomes','Tox','Nr4a1','Bcl2','Ifng','Tnf','Il2','Gzmb','Prf1'),
        idents=c('1_SNCA','1_SNCA_PLX','8_SNCA','8_SNCA_PLX')) + RotatedAxis() + ggtitle('Exhausted T-cells')

DotPlot(SYN_merged, features = genes_vTom$Immunosuppresive[!is.na(genes_vTom$Immunosuppresive)],
        idents=c(rbind(paste0(V4.2_MYELOID, '_SNCA'), paste0(V4.2_MYELOID, '_SNCA_PLX')))) + RotatedAxis() + ggtitle('Immunosuppresive')
DotPlot(SYN_merged, features = genes_vTom$`Inflammatory Markers`[!is.na(genes_vTom$`Inflammatory Markers`)],
        idents=c(rbind(paste0(V4.2_MYELOID, '_SNCA'), paste0(V4.2_MYELOID, '_SNCA_PLX')))) + RotatedAxis() + ggtitle('Inflammatory')
DotPlot(SYN_merged, features = genes_vTom$Phagocytosis[!is.na(genes_vTom$Phagocytosis)],
        idents=c(rbind(paste0(V4.2_MYELOID, '_SNCA'), paste0(V4.2_MYELOID, '_SNCA_PLX')))) + RotatedAxis() + ggtitle('Phagocytosis')

DotPlot(SYN_merged_relabeled, features = c('Pdcd1','Ctla4','Havcr2','Lag3','Tigit','Cd244a','Cd160','Eomes','Tox','Nr4a1','Bcl2','Ifng','Tnf','Il2','Gzmb','Prf1'),
        idents=c('CD8_SNCA','CD8_SNCA_PLX', 'Treg_SNCA','Treg_SNCA_PLX', 'CD4_SNCA','CD4_SNCA_PLX'), scale = F) + RotatedAxis() + ggtitle('Exhausted T-cells')
DotPlot(SYN_merged_relabeled, features = c('Cd3e','Cd3g','Cd3d','Cd8a','Cd8b1','Pdcd1','Ctla4','Havcr2','Lag3','Tigit','Eomes','Tox','Nr4a1','Bcl2','Tbx21','Prdm1','Nfat5'),
        idents=c('CD8_SNCA','CD8_SNCA_PLX'), scale=F) + RotatedAxis() + ggtitle('Exhausted T-cells')
DotPlot(SYN_merged_relabeled, features = c('Ccr1','Ccr5','Havcr2','Stat1','Tnf','Il2','Il4','Il5','Il13','Gata3','Ccr4','Ccr8'),
        idents=c('CD4','Th2','Treg','CD8'), scale=F) + RotatedAxis() + ggtitle('Th1 and Th2')

# Correlation Plot
FeatureScatter(SYN_merged_relabeled, feature1='Cd8a', feature2='Saa3')

# Feature on UMAP Plot
FeaturePlot(SYN_merged, features='Cd8a') + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "OrRd")))
FeaturePlot(SYN_merged, features='Anpep', cols = c('grey','red'))
FeaturePlot(SYN_merged, features='Ccl5', cols = c('#c7c7c7','red'), order = T, pt.size = 1.5)

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

DoHeatmap(subset(SYN_merged_relabeled, downsample = 200, idents=c('Microglia','Mo/M\u03D5')), features = c('Olfml3','Cx3cr1','Ccl12','F10','Hp','Gda'),
          size = 3, slot='scale.data') + ggtitle('TITLE') + scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n =4, name = "RdBu"))) # Myeloid Merged Clusters
DoHeatmap(subset(SYN_averages, downsample = 100, idents=c('Microglia','Mo/M\u03D5')), features = c('Olfml3','Cx3cr1','Ccl12','Tmem119','P2ry12','Gpr34','Fcrl2',
                                                                                                   'Cd163','Ms4a7','Gda','Ccr2','Cd1d1','F7','Hp','Fabp5','Mmp12','F10',
                                                                                                   'Mrc1'),
          size = 3, slot='scale.data', draw.line=F) + ggtitle('TITLE') + scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n =4, name = "RdBu"))) # Microglia

DoHeatmap(subset(SYN_merged_relabeled, downsample = 100, idents=c('CD4','CD8','Treg','NK','B-cells','Plasma Cells','DC')), features = c('Cd3e','Cd3g','Cd3d','Cd4','Cd8a'),
          size = 3, slot='scale.data', draw.line=F) + ggtitle('TITLE') + scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n =4, name = "RdBu"))) # T cells Merged

DoHeatmap(subset(aggregate_ifnb, downsample = 100, idents=c('Mo/M\u03D5','')), features = c('Angptl4','Mgll','Ear2','Serpine1','Ccl6','Chil3','Mrc1','Ccdc80','Plet1','Igf1','Lpl',
                                                                                                     'Tgfb1','Atp6v0d2','Cd68','Mfge8','Mertk','Anpep',
                                                                                                     'Lamp1','Atg7','Tfeb','Ctss','Gba1','Vcp','Bag3','Acp5',
                                                                                                     'Itgax','Cd86','Cxcl10','Ccl5','Ccl4',
                                                                                                     'Ccl3','Stat1','Tnf','Cd40','Nos2','Cd36'),
          size = 3, slot='scale.data', draw.line=F) + ggtitle('Cluster 29') + scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n =4, name = "RdBu"))) # T Cells - Various

DoHeatmap(subset(SYN_averages, downsample = 100, idents=c('Mo/M\u03D5','PLX-M\u03D5')), features = c('Spp1','Ccl6','Axl','Mrc1','Tgfb1','Lgals9','Pros1','Stat6',
                                                                                                     'Nos2','Tnf','Ccl3','Ccl4','Ccl5','Cxcl10','Stat1'),
          size = 3, slot='scale.data', draw.line=F) + ggtitle('TITLE') + scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n =4, name = "RdBu"))) # Immunosuppresive 

DoHeatmap(subset(SYN_merged_relabeled, downsample = 500, idents=c('CD4_SNCA','CD4_SNCA_PLX')), features = c('Il4','Il5','Il13','Gata3','Ccr4','Ccr8'),
          size = 3, slot='scale.data', draw.line=F) + ggtitle('Th2') + scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n =4, name = "RdBu"))) # Th2 markers 

DoHeatmap(subset(SYN_merged_relabeled, downsample = 500, idents=c('CD4_SNCA','CD4_SNCA_PLX')), features = c(genes_vTom$Th1[1:16]),
          size = 3, slot='scale.data', draw.line=F) + ggtitle('Th1') + scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n =4, name = "RdBu"))) # Th1 markers 


# Heatmaps from Excel
DoHeatmap(subset(SYN_merged, downsample = 100, idents = V4.2_MYELOID), 
          features = hierarchical_clustering_scaled(genes_vTom$M1[!is.na(genes_vTom$M1)], V4.2_MYELOID),
          size = 3, slot='scale.data') + ggtitle('Macrophage Markers') + scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n =4, name = "RdBu")))



# Heatmap of AVERAGE expression
SYN_averages <- AverageExpression(SYN_merged, return.seurat = T)

Idents(SYN_averages) <- sapply(Idents(SYN_averages), function(x) substr(x, 2, nchar(as.character(x))))

#ident_subset <- c('0_SNCA','0_SNCA_PLX','2_SNCA','2_SNCA_PLX','5_SNCA','5_SNCA_PLX','10_SNCA','10_SNCA_PLX','13_SNCA','13_SNCA_PLX','16_SNCA','16_SNCA_PLX',
#                         '17_SNCA','17_SNCA_PLX','18_SNCA','18_SNCA_PLX','22_SNCA','22_SNCA_PLX','27_SNCA','27_SNCA_PLX','29_SNCA','29_SNCA_PLX','38_SNCA','38_SNCA_PLX')
ident_subset <- c(rbind(paste0(V4.2_MYELOID, '_SNCA'), paste0(V4.2_MYELOID, '_SNCA_PLX')))
idents_subset <- gsub('_', '-', ident_subset)
idents_subset <- paste0('g', idents_subset)

features_of_interest <- genes_vTom$Immunosuppresive_v2
features_of_interest <- genes_vTom$`Inflammatory Markers_v2`
features_of_interest <- genes_vTom$Phagocytosis_Reduced
features_of_interest <- c('Pdcd1','Ctla4','Havcr2','Lag3','Tigit','Cd244a','Cd160','Eomes','Tox','Nr4a1','Bcl2','Ifng','Tnf','Il2','Gzmb','Prf1')
features_of_interest <- c('Cd3e','Cd3g','Cd3d','Cd8a','Cd8b1','Pdcd1','Ctla4','Havcr2','Lag3','Tigit','Eomes','Tox','Bcl2','Tbx21','Prdm1') # CD8
features_of_interest <- c('Cd3e','Cd3g','Cd3d','Cd4','Pdcd1','Ctla4','Havcr2','Lag3','Tigit','Eomes','Tox','Bcl2','Nr4a1','Tbx21','Prdm1','Nfat5') # CD4
features_of_interest <- genes_vTom$Treg[1:17]
features_of_interest <- c('Foxp3','Il2ra','Ctla4','Lag3','Icos','Cd69','Tnfrsf1b','Sell','Ccr7','Lrrc32','Lap3')

DoHeatmap(subset(SYN_averages, idents = V4.2_MYELOID),  # with hierarchical clustering
          features = hierarchical_clustering_scaled(features_of_interest[!is.na(features_of_interest)], V4.2_MYELOID), label = T ,draw.lines = F)  +
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n =4, name = "RdBu")))

DoHeatmap(subset(SYN_averages, idents = c(rbind(paste0(V4.2_MYELOID, '_SNCA'), paste0(V4.2_MYELOID, '_SNCA_PLX')))),  # for split between SNCA and PLX
          features=features_of_interest[!is.na(features_of_interest)], label = T ,draw.lines = F)  +
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n =4, name = "RdBu"))) + ggtitle('Phagocytosis')

DoHeatmap(subset(SYN_averages, idents = c('Mo/M\u03D5','PLX-M\u03D5')),  # for split between SNCA and PLX but for MERGED clusters
          features=features_of_interest[!is.na(features_of_interest)], label = T ,draw.lines = F)  +
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n =4, name = "RdBu"))) + ggtitle('Phagocytosis')

DoHeatmap(subset(aggregate_ifnb, idents = c('CD8_SNCA','CD8_SNCA-PLX')),  # for split between SNCA and PLX but for MERGED clusters, CD8
          features=features_of_interest[!is.na(features_of_interest)], label = T ,draw.lines = F)  +
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n =4, name = "RdBu"))) + ggtitle('Exhausted T-cell')

DoHeatmap(subset(aggregate_ifnb, idents = c('CD4_SNCA','CD4_SNCA-PLX')),  # for split between SNCA and PLX but for MERGED clusters, CD4
          features=features_of_interest[!is.na(features_of_interest)], label = T ,draw.lines = F)  +
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n =4, name = "RdBu"))) + ggtitle('Exhausted T-cell')

DoHeatmap(subset(aggregate_ifnb, idents = c('CD8_SNCA','CD8_SNCA-PLX','Treg_SNCA','Treg_SNCA-PLX','CD4_SNCA','CD4_SNCA-PLX')),  # for highlighting Treg
          features=features_of_interest[!is.na(features_of_interest)], label = T ,draw.lines = F)  +
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n =4, name = "RdBu"))) + ggtitle('Treg')

DoHeatmap(subset(aggregate_ifnb, idents = c('CD4_SNCA','CD4_SNCA-PLX','Th2_SNCA','Th2_SNCA-PLX','Treg_SNCA','Treg_SNCA-PLX')),  # for highlighting Treg
          features=c('Ccr1','Ccr5','Havcr2','Stat1','Tnf','Il2','Il4','Il5','Il13','Gata3','Ccr4','Ccr8'), label = T ,draw.lines = F)  +
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n =4, name = "RdBu"))) + ggtitle('Th1 and Th2')

DoHeatmap(subset(SYN_merged_relabeled, downsample=100, idents = c('CD4_SNCA','CD4_SNCA_PLX','Th2_SNCA','Th2_SNCA_PLX','Treg_SNCA','Treg_SNCA_PLX')),
          features=c('Ccr1','Ccr5','Havcr2','Stat1','Tnf','Il2','Il4','Il5','Il13','Gata3','Ccr4','Ccr8'), label = T ,draw.lines = F)  +
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n =4, name = "RdBu"))) + ggtitle('Th1 and Th2')

DoHeatmap(subset(SYN_averages, idents = c('g6-SNCA','g6-SNCA-PLX','g7-SNCA','g7-SNCA-PLX','Th2-SNCA','Th2-SNCA-PLX','g24-SNCA','g24-SNCA-PLX')),
          features=c('Ccr1','Ccr5','Havcr2','Stat1','Tnf','Il2','Il4','Il5','Il13','Gata3','Ccr4','Ccr8'), label = T ,draw.lines = F)  +
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n =4, name = "RdBu"))) + ggtitle('Th1 and Th2')

DoHeatmap(subset(SYN_averages, idents = c('6_SNCA','6_SNCA_PLX','Th2_SNCA','Th2_SNCA_PLX','7_SNCA','7_SNCA_PLX','24_SNCA','24_SNCA_PLX')),
          features=c('Ccr1','Ccr5','Havcr2','Stat1','Tnf','Il2','Il4','Il5','Il13','Gata3','Ccr4','Ccr8'), label = T ,draw.lines = F)  +
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n =4, name = "RdBu"))) + ggtitle('Th1 and Th2')

DoHeatmap(subset(aggregate_ifnb,idents = c('Mo/M\u03D5_SNCA','Mo/M\u03D5_SNCA-PLX')),
          features=c('Alcam','Hif1a','Lgals9','Irf8','Tgfb1','Tgfbr1','Cd163','Itgae','Spp1','Ccl17','Ccl22','Sema4a','Tnfsf4'), label = T ,draw.lines = F)  +
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n =4, name = "RdBu"))) + ggtitle('Macrophage')

DoHeatmap(subset(aggregate_ifnb,idents = c('Mo/M\u03D5_SNCA','Mo/M\u03D5_SNCA-PLX')),
          features=c('Alcam','Hif1a','Lgals9','Irf8','Tgfb1','Tgfbr1','Cd163','Itgae','Spp1','Ccl17','Ccl22','Sema4a','Tnfsf4'), label = T ,draw.lines = F)  +
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n =4, name = "RdBu"))) + ggtitle('Macrophage')

DoHeatmap(subset(SYN_averages,idents = c('10-SNCA','10-SNCA-PLX','16-SNCA','16-SNCA-PLX','17-SNCA','17-SNCA-PLX','22-SNCA','22-SNCA-PLX','27-SNCA','27-SNCA-PLX')),
          features=c('Alcam','Hif1a','Lgals9','Irf8','Tgfb1','Tgfbr1','Cd163','Itgae','Spp1','Ccl17','Ccl22','Sema4a','Tnfsf4'), label = T ,draw.lines = F)  +
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n =4, name = "RdBu"))) + ggtitle('Macrophage')

DoHeatmap(subset(aggregate_ifnb,idents = c('CD8_SNCA','CD8_SNCA-PLX')),
          features=c('Cd3d','Cd3e','Cd8a','Ctla4','Pdcd1','Lag3','Havcr2','Cd38','Tox','Prdm1','Ppargc1b','Cd5','Cd44','Entpd1','Irf4'), label = T ,draw.lines = F)  +
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n =4, name = "RdBu"))) + ggtitle('CD8')

# Ridgemap
RidgePlot(SYN_merged_relabeled, features = 'C1qa', ncol = 2)

# AddModuleScore

SYN_merged <- AddModuleScore(SYN_merged, features=c('Il4','Il5','Il13','Gata3','Ccr4','Ccr8'), name='Th2')
umap_th1 <- FeaturePlot(SYN_merged, features='Th21', label=T, pt.size = 1.5, order = T) + ggtitle('Th2', subtitle = paste0(c('Il4','Il5','Il13','Gata3','Ccr4','Ccr8'),collapse = ','))

SYN_merged <- AddModuleScore(SYN_merged, features=c('Ccr1','Ccr5','Ifng','Havcr2','Stat1','Tnf','Il2','Cxcr3','Tbx21','Il12rb2','Il18r1','Ifngr1','Il27ra','Fasl'), name='Th1')
FeaturePlot(SYN_merged, features='Th11', label=T, pt.size = 1.5, order = T) + 
  ggtitle('Th1',subtitle=paste0(c('Ccr1','Ccr5','Ifng','Havcr2','Stat1','Tnf','Il2','Cxcr3','Tbx21','Il12rb2','Il18r1','Ifngr1','Il27ra','Fasl'),collapse=','))
#+ scale_colour_gradientn(colours = rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))



# --------------------------------------
# Create Manual Subclusters


selected_cells <- CellSelector(umap_th1)


























