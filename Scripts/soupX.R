# ------------------------------------------------------------------------
# Remove ambient RNA with SoupX
# ------------------------------------------------------------------------
library(SoupX)
library(glmGamPoi)

sc_SNCA = load10X('C:/Users/tbehr/Desktop/SanRaffaele/Projects/PD_Melania/Results/SNCA_outs/')
sc_PLX = load10X('C:/Users/tbehr/Desktop/SanRaffaele/Projects/PD_Melania/Results/SNCA_PLX_outs/')

raw_SNCA <- Read10X_h5('C:/Users/tbehr/Desktop/SanRaffaele/Projects/PD_Melania/Results/SNCA_outs/raw_feature_bc_matrix.h5')
raw_PLX <- Read10X_h5('C:/Users/tbehr/Desktop/SanRaffaele/Projects/PD_Melania/Results/SNCA_PLX_outs/raw_feature_bc_matrix.h5')

filt_SNCA <- Read10X_h5('C:/Users/tbehr/Desktop/SanRaffaele/Projects/PD_Melania/Results/SNCA_outs/filtered_feature_bc_matrix.h5')
filt_PLX <- Read10X_h5('C:/Users/tbehr/Desktop/SanRaffaele/Projects/PD_Melania/Results/SNCA_PLX_outs/filtered_feature_bc_matrix.h5')

srat_SNCA <- CreateSeuratObject(counts = Read10X_h5('C:/Users/tbehr/Desktop/SanRaffaele/Projects/PD_Melania/Results/SNCA_outs/filtered_feature_bc_matrix.h5'))
srat_PLX <- CreateSeuratObject(counts = Read10X_h5('C:/Users/tbehr/Desktop/SanRaffaele/Projects/PD_Melania/Results/SNCA_PLX_outs/filtered_feature_bc_matrix.h5'))

srat <- srat_PLX

srat <- SYN_split[[2]]


srat    <- SCTransform(srat, verbose = F)
srat    <- RunPCA(srat, verbose = F)
srat    <- RunUMAP(srat, dims = 1:30, verbose = F)
srat    <- FindNeighbors(srat, dims = 1:30, verbose = F)
srat    <- FindClusters(srat, verbose = T)


meta    <- srat@meta.data
umap    <- srat@reductions$umap@cell.embeddings
sc_PLX  <- setClusters(sc_PLX, setNames(meta$seurat_clusters, rownames(meta)))
sc_PLX  <- setDR(sc_PLX, umap)

# autoestimate ambient RNA profile
sc_PLX  <- autoEstCont(sc_PLX)

# genes highest in the background
head(sc_PLX$soupProfile[order(sc_PLX$soupProfile$est, decreasing = T), ], n = 50)

out <- adjustCounts((sc_PLX))

# identify which genes have been most strongly decreased
cntSoggy = rowSums(sc_PLX$toc > 0)
cntStrained = rowSums(out > 0)
mostZeroed = tail(sort((cntSoggy - cntStrained)/cntSoggy), n = 20)
mostZeroed


plotChangeMap(sc_PLX, out, 'Trpm3')


# can create Seurat object out of the 'out' object
SNCA_PLX_adjusted <- CreateSeuratObject(out)




#### ------------------------
# Alternate Workflow
# --------------------------

whichCondition <- 'PLX'

if(whichCondition=='SNCA'){
  srat <- srat_SNCA
  soup.channel  <- SoupChannel(raw_SNCA, filt_SNCA)
} else if(whichCondition=='PLX'){
  srat <- srat_PLX
  soup.channel  <- SoupChannel(raw_PLX, filt_PLX)
}

srat <- NormalizeData(srat)
srat <- FindVariableFeatures(srat)
srat <- ScaleData(srat, features = rownames(srat))

srat <- RunPCA(srat, verbose = F)
srat <- RunUMAP(srat, dims = 1:20, verbose = F)
srat <- FindNeighbors(srat, dims = 1:20)
srat <- FindClusters(srat, verbose = T)


meta <- srat@meta.data
umap <- srat@reductions$umap@cell.embeddings
soup.channel <- setClusters(soup.channel, setNames(meta$seurat_clusters, rownames(meta)))
soup.channel <- setDR(soup.channel, umap)
head(meta)


soup.channel  <- autoEstCont(soup.channel)

head(soup.channel$soupProfile[order(soup.channel$soupProfile$est, decreasing = T), ], n = 20)

adj.matrix  <- adjustCounts(soup.channel)

# create Seurat object from output
if(whichCondition=='SNCA'){
  SNCA_adjusted <- CreateSeuratObject(adj.matrix, project = "SNCA", min.cells = 3, min.features = 500)
} else if(whichCondition=='PLX'){
  SNCA_PLX_adjusted <- CreateSeuratObject(adj.matrix, project = "SNCA_PLX", min.cells = 3, min.features = 500)
}






















