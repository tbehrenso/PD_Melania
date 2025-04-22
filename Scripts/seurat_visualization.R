

# ----------------------------------------------------------------------------------------
# Quality Control
VlnPlot(SYN_merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", 'percent.ribo'), ncol = 2)
FeatureScatter(SYN_merged, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")


VizDimLoadings(SYN_merged, dims = 1:2, reduction = "pca")
DimPlot(SYN_merged, reduction = "pca") + NoLegend()

# ----------------------------------------------------------------------------------------
# UMAP

DimPlot(SYN_merged, reduction = 'umap', group.by = c('orig.ident', 'seurat_clusters'))

DimPlot(SYN_merged, reduction = 'umap', split.by = 'orig.ident')

# ----------------------------------------------------------------------------------------
# Investigate PCA
DimHeatmap(SYN_merged, dims = 1:19, cells = 500, balanced = TRUE)

ElbowPlot(SYN_merged, ndims=50)

# ----------------------------------------------------------------------------------------
# Violin Plot
VlnPlot(SYN_merged, features = "Tmem119", split.by = "orig.ident", group.by = "seurat_clusters") 
#  + geom_boxplot(width=0.1, fill="white")



# ----------------------------------------------------------------------------------------
# Various Gene Visualizations

library(viridis)

# DotPlot
DotPlot(SYN_merged, features = c('Cd3e', 'Cd4', 'Cd8a', 'Batf','Irf4','Rora','Rorg','Stat3','Il17a','Il17f','Il21','Il22','Ccr4','Il1r1','Il6ra','Il21r','Il23r','Tgfbr2'),
        idents=c(1,3,5,9,10,12,24,31)) + RotatedAxis()
DotPlot(SYN_merged, features = c('P2ry12','Tmem119','Sparc','Clec12a','Ccr2','Fcrl2','Ifitm10','Ifitm3','Cd74','H2-Aa','H2-Eb1','Apoe','Saa3','Fth1','Lyz2',
                                 'C1qa','C1qb','C1qc','C3ar1','C5ar1','Trem2','Sall1', 'Hexb','Gpr34','Olfml3','Tgfbr1','Mertk','Pros1','Tyro3', 'Adgre1', 
                                 'Ccl6','Siglech', 'Lyve1','Cd44','Mrc1','Cd68','Itgam','Ly6g', 'Arg1','Cd163','Retnla','Pdcd1lg2','Pparg','Csf1r','Cd40','Cd80','Cd86'),
        idents=c(0,2,6,11,13,14,15,17,21,22,23,25)) + RotatedAxis() # microglia vs macrophages (from meeting)
DotPlot(SYN_merged, features = c('Gata3', 'Irf4', 'Stat5a', 'Stat6', 'Il4','Il5','Il9','Il10','Il13','Il21'),idents=c(1,3,5,9,10,12,24,31)) + RotatedAxis()
DotPlot(SYN_merged, features = c('Ahr','Batf','Stat3','Ccr4','Ccr6','Ccr10','Il22','Tnfa'),idents=c(1,3,5,9,10,12,24,31)) + RotatedAxis()
DotPlot(SYN_merged, features = c('Saa3','Apoe','Lyz2','C1qb','Il1b','Cxcl2','Cd74','H2-Eb1','Aif1','Ckb','Ramp1','Tmem176a','Capg','Ermn'), 
        cols = c("red", "blue"), split.by = 'orig.ident') + RotatedAxis() # highly diff.exp.
DotPlot(SYN_merged, features = c('Foxp3','Cd4','Il2ra'),idents=c(1,3,5,9,10,12,24,31)) + RotatedAxis() #Treg
DotPlot(SYN_merged, features = c('Cd44','Sell','Il7r','Cd69','Itgae','Ctla4'),idents=c(1,3,5,9,10,12,24,31)) + RotatedAxis() # Tmem

DotPlot(SYN_merged, features = c('Saa3','Apoe','Lyz2','Cd74','H2-Ab1')) + RotatedAxis()  # over-represented markers

DotPlot(SYN_merged, features = c('Ly6c2', 'Ly6c1', 'Ccr2', 'Cd74', 'Chil3', 'Retnlg', 'Mrc1', 'F13a1', 'Clec4e')) + RotatedAxis() # should not be in microglia
DotPlot(SYN_merged, features = c('Ly6c2', 'Gzmb', 'Il2ra', 'Nos2', 'Oas3', 'Ms4a8a', 'Arg2', 'Trem1', 'Ccr2', 'Vim', 'Ifi204', 'S100a10', 'Msrb1'), 
        idents=c(0,2,6,11,13,14,15,17,21,22,23,25)) + RotatedAxis() # higher in macrophage (Spiteri 2021) (human)
DotPlot(SYN_merged, features = c('Bmpr1a', 'Il12b','Gas6','Tnf','Cd74','Ccl12','Csf1','Ly86','Bst2','H2-aa','H2-ab1','Ifnb1','Stat1','Tlr2','Tlr3'), 
        idents=c(0,2,6,11,13,14,15,17,21,22,23,25)) + RotatedAxis() # higher in microgila (Spiteri 2021) (human)

# Correlation Plot
FeatureScatter(SYN_merged_relabeled, feature1='Cd8a', feature2='Saa3')

# Feature on UMAP Plot
FeaturePlot(SYN_merged, features='Tmem119')
# w/ Coexpression
FeaturePlot(SYN_merged, features = c('Itgam','Itgae'), blend = TRUE)
# Interactive (not working)
InteractivePlot <- FeaturePlot(SYN_merged_relabeled, features = 'C1qa')
HoverLocator(plot = InteractivePlot, information = FetchData(SYN_merged_relabeled, vars = c('orig.ident', 'PC_1', 'nFeature_RNA')))

# Heatmap
DoHeatmap(subset(SYN_merged_relabeled, downsample = 100), features = c('Cd3d', 'Cd3e','Cd3g', 'Cd4', 'Cd8a'), size = 3)

DoHeatmap(subset(SYN_merged, idents=c('11','15','27')), features = c('Fscn1','Cacnb3','Cd63','Tmem123','Cd40','Cd80','Cd83','Cd86', 'Sod2','Mt2',
                                                                     'Cxcl9','Cxcl10','Cxcl11','Ccl5','Ccl22', 'Cxcl12',
                                                                     'Cd209a','Fcer1g','Itgax','Itgam','Lamp2','Ly6c1',
                                                                     'Clec9a','Bst2','Siglech','Itgae'),
          size = 3, slot='scale.data') + scale_fill_viridis(option='magma') # Dendritic Cells
DoHeatmap(subset(SYN_merged, downsample=100, idents=c('0','2','6','13','14','17','21','23','25')), 
          features = rownames(head(FindMarkers(SYN_merged, ident.1 = '25', ident.2=c('0','2','6','13','14','17','21','23'), only.pos=T), 30)),
          size = 3, slot='scale.data') + scale_fill_viridis(option='magma') + ggtitle('Cluster 25: Top Markers') # Top Markers for Cluster 25

# Ridgemap
RidgePlot(SYN_merged_relabeled, features = 'C1qa', ncol = 2)
































