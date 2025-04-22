# ------------------------------------------------------------------------
# DoubletFinder
# (run after UMAP)
# ------------------------------------------------------------------------

library(DoubletFinder)


SYN_split <- SplitObject(SYN_merged, split.by='orig.ident')

for (i in 1:length(SYN_split)){
  print(paste0("Sample ",i))
  
  # Pre-process seurat object with standard seurat workflow
  SYN_sample <- SYN_split[[i]]
  # SYN_sample <- NormalizeData(SYN_split[[i]])
  # SYN_sample <- FindVariableFeatures(SYN_sample)
  # SYN_sample <- ScaleData(SYN_sample)
  # SYN_sample <- RunPCA(SYN_sample)
  # 
  # # finish pre-processing
  # SYN_sample <- RunUMAP(SYN_sample, dims = 1:30)
  # SYN_sample <- FindNeighbors(object = SYN_sample, dims = 1:30)              
  # SYN_sample <- FindClusters(object = SYN_sample, resolution = 1.0)
  
  # pK identification (no ground-truth)
  sweep.list <- paramSweep(SYN_sample, PCs = 1:30, sct=FALSE)
  sweep.stats <- summarizeSweep(sweep.list, GT=FALSE)  # GT stands for "ground truth"
  bcmvn <- find.pK(sweep.stats)
  
  # Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
  bcmvn.max <- bcmvn[which.max(bcmvn$BCmetric),]
  optimal.pk <- bcmvn.max$pK
  optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]
  
  ## Homotypic doublet proportion estimate
  annotations <- SYN_sample@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations) 
  nExp.poi <- round(0.075 * nrow(SYN_sample@meta.data)) 
  nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop))
  
  # run DoubletFinder
  SYN_sample <- doubletFinder(seu = SYN_sample, 
                              PCs = 1:30, 
                              pK = 0.1,
                              nExp = nExp.poi.adj, sct = FALSE
  )
  metadata <- SYN_sample@meta.data
  colnames(metadata)[9] <- "doublet_finder"
  SYN_sample@meta.data <- metadata 
  
  DimPlot(SYN_sample, reduction = "umap", group.by = c('DF.classifications_0.25_0.1_675'))
  
  # subset and save
  mouse.singlets <- subset(SYN_sample, doublet_finder == "Singlet")
  SYN_split[[i]] <- mouse.singlets
  remove(mouse.singlets)
  
}