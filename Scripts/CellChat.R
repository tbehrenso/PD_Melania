library(NMF)
library(circlize)
library(ComplexHeatmap)
library(CellChat)
library(Seurat)

SYN_merged <- readRDS('SYN_merged_FINALv4.RData')
SYN_merged <- FindClusters(SYN_merged, resolution = 2.2)
SYN_merged$RNA_snn_res.1 <- NULL

# Split Seurat Object
SYN_split <- SplitObject(SYN_merged, split.by = 'orig.ident')

SYN_SNCA <- SYN_split[[1]]
SYN_PLX <- SYN_split[[2]]
SYN_PLX <- subset(SYN_PLX, idents='42', invert=T)  # remove cluster 42 (no corresponding cells in SNCA condition)
SYN_PLX <- subset(SYN_PLX, idents='Erythrocytes', invert=T)

# rename clusters (for CellChat, a cluster cannot be called just '0')
SYN_SNCA$clusters_renamed <- paste0('c',SYN_SNCA$seurat_clusters)
SYN_PLX$clusters_renamed <- paste0('c',SYN_PLX$seurat_clusters)

# # manually start from a Seurat object (not used)
# data.input <- SYN_merged[["RNA"]]$data # normalized data matrix
# labels <- Idents(SYN_merged)
# meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels


# directly create CellChat object
cellChat_SNCA <- createCellChat(object = SYN_SNCA, group.by = "clusters_renamed", assay = "RNA")
cellChat_PLX <- createCellChat(object = SYN_PLX, group.by = "clusters_renamed", assay = "RNA")

# ---------------------------------------------------------------
# For analyzing a single dataset
# ---------------------------------------------------------------

# SELECT SAMPLE
cellChat <- cellChat_PLX

# Set ligand-receptor interaction database
CellChatDB <- CellChatDB.mouse      
# exclude Non-protin signalling (can use different subset)
CellChatDB.use <- subsetDB(CellChatDB)

cellChat@DB <- CellChatDB.use


### Preprocessing
# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellChat) # This step is necessary even if using the whole database

#future::plan("multisession", workers = 4) # do parallel

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

cellchat <- smoothData(cellchat, adj = PPI.mouse)


### Inference of Cell-Cell Communication Network

cellchat <- computeCommunProb(cellchat, type = "triMean", population.size=F)  # by default, requires 10 cells per group

# Change order of groups
# cellchat <- updateClusterLabels(cellchat, new.order = paste0('c',0:41))
cellchat <- updateClusterLabels(cellchat, new.order = c('c1','c8','c12','c11','c24','c6','c7','c9',
                                                        'c4','c19','c28','c3','c40','c30','c20','c25','c21','c15','c31',
                                                        'c18','c0','c38','c2','c5','c29','c13','c10','c27','c16','c22','c17',
                                                        'c14','c23','c33','c34','c35','c41','c39','c26','c32','c36','c37'))

# keep only interactions with at least 10 cells
cellchat <- filterCommunication(cellchat, min.cells = 10)

# Subset communication (get dataframe)
df.net <- subsetCommunication(cellchat, sources.use = c('c18'), targets.use = c('c8'))
df.net <- subsetCommunication(cellchat, signaling = c("WNT", "TGFb"))

# Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)

#Calculate aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)

# Visualze aggregated network
groupSize <- as.numeric(table(cellchat@idents))

par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

# Visualize each cluster individually
mat <- cellchat@net$weight
par(mfrow = c(2,4), xpd=TRUE)
for (i in 1:8) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

# Circle plots
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = 'CXCL', layout = "circle")


# ---------------------------------------------------------------    # START BY LOADING cellchat_XXX_groupedClusters.RData (already analyzed versions)
# For combining multiple datasets (SNCA + SNCA_PLX)
# ---------------------------------------------------------------

data.dir <- './Results/CellChat'
setwd(data.dir)

cellChat_SNCA <- netAnalysis_computeCentrality(cellChat_SNCA)
cellChat_PLX<- netAnalysis_computeCentrality(cellChat_PLX)

# Merge cellchat objects of each condition
object.list <- list(SNCA = cellChat_SNCA, PLX = cellChat_PLX)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

# Compare total number of interactions
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

# Heatmap
gg1 <- netVisual_heatmap(cellchat)
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
gg1 + gg2

# Simplify Groups
# group.cellType <- c('CD8','CD8','NK','Treg','Treg','CD4','CD4','ProliferT',
#                     'Plasma','Plasma','Plasma','BCell','BCell','BCell','Endothilial','Proliferating','Dendritic','Dendritic','Dendritic',
#                     'Microglia18','Microglia','Microglia','Microglia','Microglia','Macrophage29','Mo/Mph','Mo/Mph','Mo/Mph','Mo/Mph','Mo/Mph','Mo/Mph',
#                     'Oligo/Astro','Oligo/Astro','Astrocytes','Astrocytes','Pericyte','Ependymal','NA','NA','NA','NA','NA')
# group.cellType <- factor(group.cellType, levels = unique(group.cellType))
# object.list <- lapply(object.list, function(x) {mergeInteractions(x, group.cellType)})
# cellchat <- mergeCellChat(object.list, add.names = names(object.list))
#
# weight.max <- getMaxWeight(object.list, slot.name = c("idents", "net", "net"), attribute = c("idents","count", "count.merged"))
# par(mfrow = c(1,2), xpd=TRUE)
# for (i in 1:length(object.list)) {
#   netVisual_circle(object.list[[i]]@net$count.merged, weight.scale = T, label.edge= T, edge.weight.max = weight.max[3],
#                    edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
# }

# Differential
# par(mfrow = c(1,2), xpd=TRUE)
# netVisual_diffInteraction(cellchat, weight.scale = T, measure = "count.merged", label.edge = T)
# netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight.merged", label.edge = T)


# "Create circle plots for each cell type separately
mat <- object.list$PLX@net$weight
par(mfrow = c(3,6), xpd=TRUE, mai=c(1,0.1,0.1,0.1))
for (i in 1:nrow(mat)){
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i,] <- mat[i,]
  groupSize <- as.numeric(table(cellchat@idents$PLX))
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

# Identify cell populations with significant changes in sending or receiving signals
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
patchwork::wrap_plots(plots = gg)

gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "CD8", signaling.exclude = 'MIF')
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "CD4", signaling.exclude = 'MIF')
patchwork::wrap_plots(plots = list(gg1,gg2))


## Identify signaling groups based on their functional similarity
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")

cellchat <- netEmbedding(cellchat, type = "functional")

cellchat <- netClustering(cellchat, type = "functional")

netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)

# Compute and visualize the pathway distance in the learned joint manifold
rankSimilarity(cellchat, type = "functional")


# Compare the overall information flow of each signaling pathway or ligand-receptor pair
gg1 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = F, do.stat = TRUE)
gg1 + gg2

# Compare outgoing (or incoming) signaling patterns associated with each cell population
library(ComplexHeatmap)
i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 20, height = 30)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 20, height = 30)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))




# ---------------------------------------------------------------------------------------------
##### Part III: Identify the up-gulated and down-regulated signaling ligand-receptor pairs

netVisual_bubble(cellchat, sources.use = c('Mo/M\u03D5','PLX-M\u03D5','Microglia'), targets.use = c('CD8'),  comparison = c(1, 2), angle.x = 45)

gg1 <- netVisual_bubble(cellchat, sources.use = c('Mo/M\u03D5','PLX-M\u03D5','Microglia'), targets.use = c('CD8'),  comparison = c(1, 2), max.dataset = 2,
                        title.name = "Increased signaling in PLX", angle.x = 45, remove.isolate = T)
gg2 <- netVisual_bubble(cellchat, sources.use = c('Mo/M\u03D5','PLX-M\u03D5','Microglia'), targets.use = c('CD8'),  comparison = c(1, 2), max.dataset = 1,
                        title.name = "Decreased signaling in PLX", angle.x = 45, remove.isolate = T)
gg1 + gg2

## Identify dysfunctional signaling by using differential expression analysis
# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "PLX"
# define a char name used for storing the results of differential expression analysis
features.name = paste0(pos.dataset, ".merged")

# perform differential expression analysis 
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name,
                                       only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.05,thresh.p = 0.05, group.DE.combined = TRUE) 

# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name, variable.all = TRUE)
# extract the ligand-receptor pairs with upregulated ligands in LS
net.up <- subsetCommunication(cellchat, net = net, datasets = "PLX",ligand.logFC = 0.05, receptor.logFC = NULL)
# extract the ligand-receptor pairs with upregulated ligands and upregulated receptors in NL, i.e.,downregulated in LS
net.down <- subsetCommunication(cellchat, net = net, datasets = "SNCA",ligand.logFC = -0.05, receptor.logFC = NULL)

gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)

df <- findEnrichedSignaling(object.list[[2]], features = c("Spp1"), idents = c('c29'), pattern ="outgoing")

## Visualize the identified up-regulated and down-regulated signaling ligand-receptor pairs
pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = c('CD8'), targets.use = c('Mo/M\u03D5','PLX-M\u03D5','Microglia'), comparison = c(1, 2),
                        angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = c('CD8'), targets.use = c('Mo/M\u03D5','PLX-M\u03D5','Microglia'), comparison = c(1, 2),
                        angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
gg1 + gg2

library(wordcloud)
computeEnrichmentScore(net.down, species = 'mouse', variable.both = TRUE) # enriched in SNCA
computeEnrichmentScore(net.up, species = 'mouse', variable.both = TRUE) # enriched in PLX


# ---------------------------------------------------------------------------------------------
# Part IV: Visually compare cell-cell communication using Hierarchy plot, Circle plot or Chord diagram

# pathways in both conditions
intersect(object.list$SNCA@netP$pathways, object.list$PLX@netP$pathways)

#circle plot
pathways.show <- c("CD45") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

#heatmap
pathways.show <- c("TGFb") 
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))

# Chord diagram
pathways.show <- c("SEMA6") 
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(object.list)[i]))
}

### Plot violin plot for genes related to pathway (with Seurat wrapper)
cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("NL", "LS")) # set factor level
plotGeneExpression(cellchat, signaling = "CD40", split.by = "datasets", colors.ggplot = T, type = "violin")











