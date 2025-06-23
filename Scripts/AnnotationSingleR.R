library(SingleR)
library(celldex)
library(SingleCellExperiment)

if(SCTRANSFORM){
  SYN_merged <- JoinLayers(SYN_merged, assay = 'RNA')
}

# load reference dataset
immgen_ref <- celldex::ImmGenData()
mouserna_ref <- celldex::MouseRNAseqData()

# convert Seurat object to SCE object (using raw clusters)
SYN_sce <- as.SingleCellExperiment(SYN_merged)

# add labels to Seurat object
pred.cnts <- SingleR::SingleR(test = SYN_sce, ref = immgen_ref, labels = immgen_ref$label.main)
lbls.keep <- table(pred.cnts$labels)>10
SYN_merged$SingleR.labels <- ifelse(lbls.keep[pred.cnts$labels], pred.cnts$labels, 'Other')
DimPlot(SYN_merged, reduction='umap', group.by='SingleR.labels')

# add labels to Seurat object (fine)
pred.cnts.fine <- SingleR::SingleR(test = SYN_sce, ref = immgen_ref, labels = immgen_ref$label.fine)
lbls.keep.fine <- table(pred.cnts.fine$labels)>10
SYN_merged$SingleR.labels.fine <- ifelse(lbls.keep.fine[pred.cnts.fine$labels], pred.cnts.fine$labels, 'Other')
DimPlot(SYN_merged, reduction='umap', group.by='SingleR.labels.fine', split.by='orig.ident')

# add labels to Seurat object (fine, MouseRNAseqData)
pred.cnts.fine.mouserna <- SingleR::SingleR(test = SYN_sce, ref = mouserna_ref, labels = mouserna_ref$label.fine)
lbls.keep.fine.mouserna <- table(pred.cnts.fine.mouserna$labels)>10
SYN_merged$SingleR.labels.fine.mouserna <- ifelse(lbls.keep.fine.mouserna[pred.cnts.fine.mouserna$labels], pred.cnts.fine.mouserna$labels, 'Other')
DimPlot(SYN_merged, reduction='umap', group.by='SingleR.labels.fine.mouserna')

# T-cell labels
lbls.keep.fine.Tcells <- lbls.keep.fine
lbls.keep.fine.Tcells[!grepl("T cell", names(lbls.keep.fine.Tcells))] <- FALSE
SYN_merged$SingleR.labels.fine.Tcells <- ifelse(lbls.keep.fine.Tcells[pred.cnts.fine$labels], pred.cnts.fine$labels, 'Other')
DimPlot(SYN_merged, reduction='umap', group.by='SingleR.labels.fine.Tcells')


plotScoreHeatmap(pred.cnts)
plotScoreHeatmap(pred.cnts.fine)

plotDeltaDistribution(pred.cnts, ncol = 3)


plotMarkerHeatmap(pred.cnts, SYN_sce, label="Microglia")