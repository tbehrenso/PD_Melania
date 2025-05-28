
# ------------------------------------------------------------------------
# Hierarchical Clustering of Gene Set
# ------------------------------------------------------------------------

hierarchical_clustering <- function(gene_subset, idents_subset){
  
  SYN_subset <- subset(SYN_merged, idents = idents_subset)
  
  avg_exp <- AverageExpression(SYN_subset, features = gene_subset)$RNA
  
  d <- dist(avg_exp)
  hc <- hclust(d, method='complete')
  
  gene_order <- hc$labels[hc$order]
  
  return(gene_order)
}

hierarchical_clustering_scaled <- function(gene_subset, idents_subset){
  
  SYN_subset <- subset(SYN_merged, idents = idents_subset)
  
  scaled_data <- GetAssayData(SYN_merged, layer='scale.data')
  
  gene_subset <- gene_subset[gene_subset %in% rownames(SYN_merged)]
  
  scaled_genes <- scaled_data[gene_subset, ]
  
  cell_idents <- Idents(SYN_subset)
  
  # Compute average scaled expression per cluster
  avg_scaled <- sapply(levels(cell_idents), function(cluster) {
    cells_in_cluster <- WhichCells(SYN_subset, idents = cluster)
    rowMeans(scaled_genes[, cells_in_cluster, drop = FALSE])
  })
  
  
  d <- dist(avg_scaled)
  hc <- hclust(d)
  
  gene_order <- hc$labels[hc$order]
  
  return(gene_order)
}



# ------------------------------------------------------------------------
# Remove RP from markers
# ------------------------------------------------------------------------

remove_Rp_from_markers <- function(markers_df){
  markers_df <- markers_df[!grepl('Rp', rownames(markers_df)),]
  
  return(markers_df)
}

# ------------------------------------------------------------------------
# Create various Seurat objects
# ------------------------------------------------------------------------

add_cluster_from_barcodes <- function(seurat_object, barcodes, new_cluster_name){
  new_clusters <- as.character(Idents(seurat_object))
  new_clusters[colnames(seurat_object) %in% barcodes] <- new_cluster_name
  Idents(seurat_object) <- new_clusters
  
  return(seurat_object)
}






















