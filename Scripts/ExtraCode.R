
# ------------------------------------------------------------------------
# Create heatmap from Excel file
# ------------------------------------------------------------------------
library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)

DEA_SNCA_WT <- read_xlsx('data/DEA_SNCA_vs_WT.xlsx')
colnames(DEA_SNCA_WT)[1] <- 'Gene'

GENE_LIST <- c('Il6','Ccl2','Ccl3','Ccl4','Ccl5','Ccl7','Ccl8','Cxcl1','Cxcl10','Cxcl13') # Ccl17 NOT present

DEA_subset <- DEA_SNCA_WT[DEA_SNCA_WT$Gene %in% GENE_LIST,]

DEA_subset_long <- DEA_subset %>%
  pivot_longer(
    cols = -Gene,
    names_to = "Sample",
    values_to = "Expression"
  )

DEA_subset_long <- DEA_subset_long %>%
  mutate(
    Condition = ifelse(grepl("SNCA", Sample), "SNCA", "WT")
  )

ggplot(DEA_subset_long, aes(x = Sample, y = Gene, fill = Expression)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "blue", high = "red") +
  theme_minimal() +
  labs(title = "Gene Expression Heatmap", x = "Sample", y = "Gene") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# ------------------------------------------------------------------------
# Plot Table output for Cells per Cluster
# ------------------------------------------------------------------------

table_df <- as.data.frame(table(SYN_merged@meta.data$RNA_snn_res.1, SYN_merged@meta.data$orig.ident))

ggplot(table_df, aes(x=Var1, y=Freq, fill=Var2)) +
  geom_bar(stat='identity', position=position_dodge()) +
  theme_minimal() +
  xlab('Cluster') + ylab('Count')



# ------------------------------------------------------------------------
# Export Seurat object to Loupe
# ------------------------------------------------------------------------

library(loupeR)

create_loupe_from_seurat(SYN_merged)




# ------------------------------------------------------------------------
# Remove ambient RNA with SoupX
# ------------------------------------------------------------------------
library(SoupX)
library(glmGamPoi)

sc_SNCA = load10X('C:/Users/tbehr/Desktop/SanRaffaele/Projects/PD_Melania/Results/SNCA_outs/')
sc_PLX = load10X('C:/Users/tbehr/Desktop/SanRaffaele/Projects/PD_Melania/Results/SNCA_PLX_outs/')


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













