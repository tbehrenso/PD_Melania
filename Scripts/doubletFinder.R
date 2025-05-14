# ------------------------------------------------------------------------
# DoubletFinder
# (run after UMAP)
# ------------------------------------------------------------------------

library(DoubletFinder)

SOUPX <- TRUE

if(SOUPX){
  SYN_split <- list(SNCA_adjusted, SNCA_PLX_adjusted)
} else {
  SYN_split <- SplitObject(SYN_merged, split.by='orig.ident')
}


for (i in 1:length(SYN_split)){
  print(paste0("Sample ",i))
  
  SYN_sample <- SYN_split[[i]]
  
  SYN_sample[["percent.mt"]] <- PercentageFeatureSet(SYN_sample, pattern = "^mt-")
  SYN_sample[["percent.ribo"]] <- PercentageFeatureSet(SYN_sample, "^Rp[sl]")
  SYN_sample <- subset(SYN_sample, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10)
  
  if(SOUPX){
    # Pre-process seurat object with standard seurat workflow
    SYN_sample <- NormalizeData(SYN_split[[i]])
    SYN_sample <- FindVariableFeatures(SYN_sample)
    SYN_sample <- ScaleData(SYN_sample)
    SYN_sample <- RunPCA(SYN_sample)
    # 
    # # finish pre-processing
    SYN_sample <- RunUMAP(SYN_sample, dims = 1:20)
    #SYN_sample <- FindNeighbors(object = SYN_sample, dims = 1:20)              
    #SYN_sample <- FindClusters(object = SYN_sample, resolution = 1.0)
  }

  
  # pK identification (no ground-truth)
  sweep.list <- paramSweep(SYN_sample, PCs = 1:20, sct=FALSE)
  sweep.stats <- summarizeSweep(sweep.list, GT=FALSE)  # GT stands for "ground truth"
  bcmvn <- find.pK(sweep.stats)
  
  ggplot(bcmvn, aes(pK, BCmetric, group=1)) +
    geom_point() +
    geom_line()
  
  # Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
  bcmvn.max <- bcmvn[which.max(bcmvn$BCmetric),]
  optimal.pk <- bcmvn.max$pK
  optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]
  
  ## Homotypic doublet proportion estimate
  #annotations <- SYN_sample@meta.data$seurat_clusters
  #homotypic.prop <- modelHomotypic(annotations)
  homotypic.prop <- 0.3
  nExp.poi <- round(0.05 * nrow(SYN_sample@meta.data)) 
  nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop))
  
  # run DoubletFinder
  SYN_sample <- doubletFinder(seu = SYN_sample, 
                              PCs = 1:20, 
                              pK = optimal.pk,
                              nExp = nExp.poi.adj, sct = FALSE
  )
  metadata <- SYN_sample@meta.data
  colnames(metadata)[length(colnames(metadata))] <- "doublet_finder"
  SYN_sample@meta.data <- metadata
  
  DimPlot(SYN_sample, reduction = "umap", group.by = c('doublet_finder'))
  
  # subset and save
  mouse.singlets <- subset(SYN_sample, doublet_finder == "Singlet")
  SYN_split[[i]] <- mouse.singlets
  remove(mouse.singlets)
  
}

