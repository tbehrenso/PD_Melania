library(dplyr)
library(Seurat)
library(patchwork)
library(hdf5r)
library(ggplot2)
library(readxl)
source('Scripts/PD_functions.R')

SCTRANSFORM <- FALSE
SOUPX <- FALSE

# Load the SNCA dataset
SNCA_out <- Read10X_h5(filename = "Results/SNCA_outs/filtered_feature_bc_matrix.h5")
SNCA_PLX_out <- Read10X_h5(filename = "Results/SNCA_PLX_outs/filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
SNCA <- CreateSeuratObject(counts = SNCA_out, project = "SNCA", min.cells = 3, min.features = 500)
SNCA_PLX <- CreateSeuratObject(counts = SNCA_PLX_out, project = "SNCA_PLX", min.cells = 3, min.features = 500)

# LOad soupX Seurat objects
if(SOUPX){
  SNCA <- SNCA_adjusted
  SNCA_PLX <- SNCA_PLX_adjusted
}

# Merge (without integration)
SYN_merged <- merge(x = SNCA, y = SNCA_PLX)

# Remove high levels of mitochondrial DNA
SYN_merged[["percent.mt"]] <- PercentageFeatureSet(SYN_merged, pattern = "^mt-")
SYN_merged[["percent.ribo"]] <- PercentageFeatureSet(SYN_merged, "^Rp[sl]")


SYN_merged <- subset(SYN_merged, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 15)

if(SCTRANSFORM){
  SYN_merged <- SCTransform(SYN_merged, vars.to.regress = 'percent.mt')
} else {
  # Analyze conditions separately with standard workflow
  SYN_merged <- NormalizeData(SYN_merged)
  SYN_merged <- FindVariableFeatures(SYN_merged)
  SYN_merged <- ScaleData(SYN_merged, features = rownames(SYN_merged))  # specifying features so it also scales non-highly variable genes (for heatmap)
}

SYN_merged <- RunPCA(SYN_merged)


# ---------------------------------------------------------------
# Integration Analysis
# ---------------------------------------------------------------

# to overwrite default maximum allowed of globals
options(future.globals.maxSize = 75000 * 1024^2)

# Integrate layers
if(SCTRANSFORM){
  SYN_merged <- IntegrateLayers(object = SYN_merged, method = RPCAIntegration, orig.reduction = "pca", new.reduction = "integrated.rpca", normalization.method = 'SCT',
#                                assay='SCT',
                                verbose = FALSE)
  
  
  
} else {
  # NOTE: RPCA integration is faster, but more conservative
  SYN_merged <- IntegrateLayers(object = SYN_merged, method = RPCAIntegration, orig.reduction = "pca", new.reduction = "integrated.rpca",
                                verbose = FALSE)
  SYN_merged[["RNA"]] <- JoinLayers(SYN_merged[['RNA']])
  
}


SYN_merged <- FindNeighbors(SYN_merged, reduction = "integrated.rpca", dims = 1:20)
SYN_merged <- FindClusters(SYN_merged, resolution = 2.2)

SYN_merged <- RunUMAP(SYN_merged, dims= 1:20, reduction = 'integrated.rpca')

if(SCTRANSFORM){
  # run to be compatible with FindMarkers
  SYN_merged <- PrepSCTFindMarkers(SYN_merged)
}



#### Finding differentially expressed features (cluster biomarkers) (part of STANDARD workflow)

# find all markers for clusters
clusterX.markers <- FindMarkers(SYN_merged, ident.1 = '8', only.pos=T)
head(clusterX.markers, n = 50)

if(F){
  SYN_markers <- FindAllMarkers(SYN_merged)
  
  SYN_markers_filtered <- SYN_markers %>% 
    group_by(cluster) %>% 
    dplyr::filter(avg_log2FC > 1)
}


# Get cells per cluster (between both conditions)
table(SYN_merged@meta.data$RNA_snn_res.1, SYN_merged@meta.data$orig.ident)




#----------------------------------------------
#  Annotate cell types

# new names (5 per row)
new.cluster.ids <- c('Microglia',
                     'CD8','Microglia','B-cells','Plasma Cells','Microglia',
                     'CD4','CD4','CD8','NA','Mo/M\u03D5',
                     'Treg','NK','Microglia','Oligodendrocytes','DC',
                     'Mo/M\u03D5','Mo/M\u03D5','Microglia','Plasma Cells','Endothelial',
                     'DC','Mo/M\u03D5','Oligodendrocytes','CD4','NA',
                     'Neutrophils','Mo/M\u03D5','Plasma Cells','Mo/M\u03D5','B-cells',
                     'DC','DC','NA','NA','Pericytes',
                     'NA','NA','Microglia','Mast Cells','B-cells',
                     'Ependymal','Erythrocytes')
new.cluster.ids <- c('Microglia',
                     'CD8','Microglia','B-cells','Plasma Cells','Microglia',
                     'CD4','CD4','CD8','NA','Mo/M\u03D5',
                     'Treg','NK','Microglia','Oligodendrocytes','DC',
                     'Mo/M\u03D5','Mo/M\u03D5','Microglia','Plasma Cells','Endothelial',
                     'DC','Mo/M\u03D5','Oligodendrocytes','CD4','NA',
                     'Neutrophils','Mo/M\u03D5','Plasma Cells','PLX-M\u03D5','B-cells',
                     'DC','DC','NA','NA','Pericytes',
                     'NA','NA','Microglia','Mast Cells','B-cells',
                     'Ependymal','Erythrocytes')
names(new.cluster.ids) <- levels(SYN_merged)
SYN_merged_relabeled <- RenameIdents(SYN_merged, new.cluster.ids)
DimPlot(SYN_merged_relabeled, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

# for plotting more nicely, with legend
DimPlot(SYN_merged_relabeled, reduction = "umap", label = TRUE, label.size = 4.5) + xlab("UMAP 1") + ylab("UMAP 2") +
  theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) + guides(colour = guide_legend(override.aes = list(size = 10)))



# --------------------------------------------
# Identify conserved cell type markers
library(metap)

nk.markers <- FindConservedMarkers(SYN_merged_relabeled, ident.1 = "NK", grouping.var = "orig.ident", verbose = FALSE)
head(nk.markers, 20)


markers.to.plot <- c('Cd3d', 'Cd3e','Cd3g', 'Cd4', 'Cd8a')
DotPlot(SYN_merged, features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8, split.by = 'orig.ident') +
  RotatedAxis()


# --------------------------------------------
# Identify differential expressed genes across conditions
library(cowplot)
theme_set(theme_cowplot())

aggregate_ifnb <- AggregateExpression(SYN_merged_relabeled, group.by = c('ident', 'orig.ident'), return.seurat = TRUE)
genes.to.label = c('Saa3','Apoe','Lyz2','Ckb','Ramp1', 'C1qb','Il1b','Cxcl2','Cd74','Inpp4b','Themis','Itpkb')

p1 <- CellScatter(aggregate_ifnb, 'CD4_SNCA', 'CD4_SNCA-PLX', highlight = genes.to.label)
p2 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)

p3 <- CellScatter(aggregate_ifnb, "CD8_SNCA", "CD4_SNCA-PLX", highlight = genes.to.label)
p4 <- LabelPoints(plot = p3, points = genes.to.label, repel = TRUE)

p2 + p4


# For Base Clusters
seurat_clusters_factor <- factor(SYN_merged$seurat_clusters)

SYN_merged$celltype.stim <- paste(Idents(SYN_merged), SYN_merged$orig.ident, sep = "_")
Idents(SYN_merged) <- factor(SYN_merged$celltype.stim, levels = c(rbind(paste0(c(0:42),'_SNCA'), paste0(c(0:42),'_SNCA_PLX'))))
#Idents(SYN_merged) <- 'celltype.stim'
condition.response <- FindMarkers(SYN_merged, ident.1 = "5_SNCA", ident.2 = "5_SNCA_PLX", verbose = FALSE)
head(condition.response, n = 30)

if(F){
  library(openxlsx)
  
  for(cluster in 0:42){
    print(cluster)
    
    condition.response <- FindMarkers(SYN_merged, ident.1 = paste0(cluster,'_SNCA'), ident.2 = paste0(cluster,'_SNCA_PLX'), verbose = FALSE)
    
    #write.xlsx(condition.response, paste0('Results/SvSPLX/SvSPLX_', cluster, '.xlsx'))
    write.csv(condition.response, file = paste0('Results/SvSPLX/v4_SvSPLX_', cluster, '.csv'), quote = F)
  }
}

# for Relabeled clusters
SYN_merged_relabeled$celltype.stim <- paste(Idents(SYN_merged_relabeled), SYN_merged_relabeled$orig.ident, sep = "_")
Idents(SYN_merged_relabeled) <- SYN_merged_relabeled$celltype.stim
condition.response <- FindMarkers(SYN_merged_relabeled, ident.1 = "CD4_SNCA", ident.2 = "CD4_SNCA_PLX", verbose = FALSE)
head(condition.response, n = 50)



