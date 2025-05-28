
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

table_df <- as.data.frame(table(SYN_merged@meta.data$seurat_clusters, SYN_merged@meta.data$orig.ident))

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
# "Reset" Seurat object from DoubletFinder
# ------------------------------------------------------------------------

SYN_merged <- merge(SYN_split[[1]], SYN_split[[2]])

SNCA_raw <- GetAssayData(SYN_split[[1]], layer='counts')
PLX_raw <- GetAssayData(SYN_split[[2]], layer='counts')

SNCA <- CreateSeuratObject(counts = SNCA_raw, project = "SNCA", min.cells = 3, min.features = 500)
SNCA_PLX <- CreateSeuratObject(counts = PLX_raw, project = "SNCA_PLX", min.cells = 3, min.features = 500)


# ------------------------------------------------------------------------
# Calculate Principal Component Threshold
# ------------------------------------------------------------------------

pct <- SYN_merged[["pca"]]@stdev / sum(SYN_merged[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]

co1

co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

# last point where change of % of variation is more than 0.1%.
co2

pcs <- min(co1, co2)

# ------------------------------------------------------------------------
# clustree
# ------------------------------------------------------------------------
library(clustree)

clustree(SYN_merged)


# ------------------------------------------------------------------------
# Add Module Score
# ------------------------------------------------------------------------

SYN_merged <- AddModuleScore(SYN_merged,
                             features = c('P2ry13','Cx3cr1','Maf','Ltc4s','Olfml3','Bin1','Ifngr1','Frmd4a','Alox5ap','Lpar6','Csf1r','Selplg','Slco2b1','P2ry12','Stab1','Cysltr1','Gpr34',
                                          'Fscn1','Adrb2','Sting1','Pmepa1','Marcks','Zfhx3','Ptgs1','Srgap2','Slc2a5','Ssh2','Rps6ka1','Ccr5','Cmtm6','Txnip',
                                          'Adgrg1','Otulinl','Rnase4','Sparc','Rgs2','Lrba'),
                             name='Homeostatic_Microglia')
SYN_merged <- AddModuleScore(SYN_merged,
                             features = c('Cd9','Spp1','Cd83','Ramp1','Ctsd','Ccl3','Trem2','Gpnmb','Itgax','Lpl','Cd63','Apoe','Npc2','Ctsb','Ctsz','Axl','Rplp1','Cadm1','Ank','Cst7',
                                          'Clec7a','H2-D1','Fth1','Serpine2','Ctsl','Lgals3bp','Cd52','H2-K1','Cd68','Csf1','Or5v1'),
                             name='DAMLike_Microglia')
SYN_merged <- AddModuleScore(SYN_merged,
                             features = c('Cd74','Apoe','Cd9','Capg','B2m','Rps19','Clec7a','Rps28','Rps14','Rpl32','Rps21','Rplp1','Npc2','Itgax','Rpl23','Ftl1','Tyrobp','Trem2','Lpl',
                                          'Cadm1','Ctsz','Cst7','H2-D1','Fth1','Serpine2','Ctsl','Lgals3bp','Cd52','H2-K1','Cd68','Csf1','Or5v1','Hexa','Ctsa',
                                          'Gusb','Baiap2l2','Igf1'),
                             name='HLA_Microglia')
SYN_merged <- AddModuleScore(SYN_merged,
                             features = c('Ccl4','Ccl3','Cd83','Spp1','Lpl','Plek','Itgax','Npc2','Axl','Apoe','Capg','Clec7a','Rps28','Ctsz','Rpl32','Rpl23','B2m','Pld3','Cst7','H2-D1',
                                          'Fth1','Serpine2','Tyrobp','Ctsl','Lgals3bp','Cd52','H2-K1','Cd68','Csf1','Cadm1','Or5v1','Hexa','Ctsa','Gusb','Rps14','Rplp1','Baiap2l2',
                                          'Igf1','Rpl18a'),
                             name='CRM_Microglia')
SYN_merged <- AddModuleScore(SYN_merged,
                             features = c('Ifit2','Bst2','B2m','Lgals3bp','Ccl3','Spp1','Axl','Apoe','Gnas','Clec7a','Cst7','Itgax','Serpine2','Csf1','Ctsl','Cd83','Igf1'),
                             name='IRM_Microglia')

FeaturePlot(SYN_merged,features = "IRM_Microglia1", label = TRUE)

# clean up seurat object
SYN_merged@meta.data <- SYN_merged@meta.data %>% select(-matches("Homeostatic_Microglia[0-9]+"))



# ------------------------------------------------------------------------
# Monocle3
# ------------------------------------------------------------------------
library(monocle3)

SYN_subset <- subset(SYN_merged, idents = c(0,3,4,12,16,20))

cds <- as.CellDataSet(SYN_subset)
cds <- cluster_cells(cds)
cds <- learn_graph(cds)
cds <- order_cells(cds)


# ------------------------------------------------------------------------
# Module Score
# ------------------------------------------------------------------------

SYN_merged <- AddModuleScore(SYN_merged,
                       features = c('P2ry12','Tmem119','Sall1','Fcrl2','Gpr34','Hexb','Olfml3','Clec7a','Siglech','Siglece','Ptprc'),
                       name='Microglia')

# Plot scores
FeaturePlot(SYN_merged,
            features = 'Microglia1', label = F)



# ------------------------------------------------------------------------
# Steps for Getting X vs Y marker excels
# ------------------------------------------------------------------------

clusterx_markers <- FindMarkers(SYN_merged, ident.1=c(29), ident.2=c(0,2,5,10,13,16,17,18,22,27,38))

clusterx_markers <- remove_Rp_from_markers(clusterx_markers)

clusterx_markers <- clusterx_markers[abs(clusterx_markers$pct.1-clusterx_markers$pct.2)>0.1,]

clusterx_markers$p_val <- NULL

write.csv(clusterx_markers, '/Users/tbehr/Desktop/29_v_myeloid.csv')


# ------------------------------------------------------------------------
# Import Barcodes and Define as Separate Cluster
# ------------------------------------------------------------------------
# (mostly converted into add_cluster_from_barcodes() function)

barcodes_th2 <- read.csv('data/Barcodes_Th2.csv', header = T)

new_clusters <- as.character(Idents(SYN_merged_relabeled))
new_clusters[colnames(SYN_merged_relabeled) %in% barcodes_th2$Barcode] <- 'Th2'
Idents(SYN_merged_relabeled) <- new_clusters


Idents(SYN_merged_relabeled)[colnames(SYN_merged_relabeled) %in% barcodes_th2$Barcode] <- 'Th2'





























