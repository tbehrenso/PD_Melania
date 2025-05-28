# Load required packages
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(ggplot2)
library(DOSE)

# ------------------------------
# Load DEG list
# 
deg_genes <- read.table("/Users/tbehr/Desktop/DEGlist.txt", header = FALSE, stringsAsFactors = FALSE)
gene_symbols <- deg_genes$V1

# ------------------------------
# Convert gene symbols to Entrez IDs
gene_df <- bitr(gene_symbols, fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = org.Mm.eg.db)

# Remove duplicates (if any)
gene_entrez <- unique(gene_df$ENTREZID)

# ------------------------------
# GO Enrichment Analysis (Biological Process)
go_bp <- enrichGO(gene = gene_entrez,
                  OrgDb = org.Mm.eg.db,
                  ont = "BP",
                  keyType = "ENTREZID",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2,
                  readable = TRUE)

# ------------------------------
# KEGG Pathway Enrichment
kegg_res <- enrichKEGG(gene = gene_entrez,
                       organism = 'mmu',  # mouse KEGG code
                       keyType = "kegg",
                       pvalueCutoff = 0.05)

# Convert Entrez back to gene symbols for readability
kegg_res <- setReadable(kegg_res, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")

# ------------------------------
# Visualizations

# GO Dotplot
dotplot(go_bp, showCategory = 20) + ggtitle("Clusterprolfiler: GO Biological Process Enrichment")

# KEGG Barplot
barplot(kegg_res, showCategory = 20) + ggtitle("Clusterprolfiler: KEGG Pathway Enrichment")

# Enrichment Map (optional)
emapplot(pairwise_termsim(go_bp))


# ------------------------------------------------------------
#  Gene Set Enrichment Analysis
# ------------------------------------------------------------

# logFC values
logfc_list <- clusterX.markers$avg_log2FC

# name vector
names(logfc_list) <- rownames(clusterX.markers)

# order by descending logFC (needed for GSEA)
logfc_list <- sort(logfc_list, decreasing=T)

# Run GSEA
gse <- gseGO(geneList = logfc_list,
             ont = 'ALL',
             OrgDb = 'org.Mm.eg.db',
             keyType = 'SYMBOL'
             )

dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)

emapplot(gse, showCategory = 10)

cnetplot(gse, categorySize="pvalue", foldChange=logfc_list, showCategory = 3)

ridgeplot(gse) + labs(x = "enrichment distribution")

gseaplot(gse, by = "all", title = gse$Description[1], geneSetID = 1)

# ------------------------------
# GSEA for KEGG

ids<-bitr(names(logfc_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb='org.Mm.eg.db')

dedup_ids <- ids[!duplicated(ids[c("SYMBOL")]),]

df2 <- clusterX.markers[rownames(clusterX.markers) %in% dedup_ids$SYMBOL,]

df2$Y <- dedup_ids$ENTREZID

kegg_gene_list <- df2$avg_log2FC

names(kegg_gene_list) <- df2$Y

# re,pve MA
kegg_gene_list<-na.omit(kegg_gene_list)

# order descending 
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)

# Run GSEA for KEGG
kk2 <- gseKEGG(geneList = kegg_gene_list,
               organism = 'mmu',
               keyType = 'ncbi-geneid',
               pvalueCutoff = 0.1)

dotplot(kk2, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
cnetplot(kk2, categorySize="pvalue", foldChange=kegg_gene_list)
ridgeplot(kk2) + labs(x = "enrichment distribution")
gseaplot(kk2, by = "all", title = kk2$Description[1], geneSetID = 1)

library(pathview)

dme <- pathview(gene.data=kegg_gene_list, pathway.id="mmu04915", species = 'mmu')










