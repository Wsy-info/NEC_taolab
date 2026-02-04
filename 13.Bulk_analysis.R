### bulk RNA-seq analysis
### GSEA, Scatter plot, and Deconvolution analysis

library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(fgsea)
library(GseaVis)
library(msigdbr)
library(org.Rn.eg.db)
library(clusterProfiler)
library(enrichplot)
library(biomaRt)
library(SingleCellExperiment)
library(BisqueRNA)


### load data
load("../bulk/DPI_bulk.RData")

### GSEA
### get marker genes
marker_inflammatory_epi <- read.table("round2_DEGs_epi_celltype.txt", sep = "\t")
marker_inflammatory_epi <- marker_inflammatory_epi %>% filter(cluster == "Inflammatory Epithelial cell") %>% top_n(n = 50, wt = avg_logFC)

### gene to rat gene
rat_mart <- useMart(biomart = "ensembl", dataset = "rnorvegicus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
human_mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
results <- getLDS(
    attributes = c("hgnc_symbol"),
    filters = "hgnc_symbol",
    values = marker_inflammatory_epi$gene,
    mart = human_mart,
    attributesL = c("rgd_symbol"),
    martL = rat_mart,
    uniqueRows = TRUE
  )

### rat gene to ENTREZID
results1 <- bitr(results$RGD.symbol,
                fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = "org.Rn.eg.db"
)
marker_gene_list <- list("GO:0000001|Inflammatory_Epi_markers" = results1$ENTREZID)

### get FCgenelist
FCgenelist <- read.table("../dpi/GSEA/NEC_Ctrl_all_genes_diff_test.txt")
FCgenelist_ctrl <- FCgenelist$log2FoldChange
names(FCgenelist_ctrl) <- rownames(FCgenelist)
gene_entrezid <- bitr(names(FCgenelist_ctrl),
                      fromType = "SYMBOL",
                      toType = "ENTREZID",
                      OrgDb = "org.Rn.eg.db"
)

rownames(gene_entrezid) <- gene_entrezid$SYMBOL
FCgenelist_ctrl <- sort(FCgenelist_ctrl[which(names(FCgenelist_ctrl) %in% gene_entrezid$SYMBOL)], decreasing = T)
names(FCgenelist_ctrl) <- gene_entrezid[names(FCgenelist_ctrl), "ENTREZID"]

### GSEA
fgsea_result_ctrl <- fgsea(marker_gene_list,
                            FCgenelist_ctrl, minSize = 5, maxSize = 500, nperm = 1000)

### Plot
fgseaRes_ctrl <- fgsea_result_ctrl %>%
  tidyr::separate_wider_delim(pathway, delim = "|", names = c("ID", "Description")) %>%
  rowwise() %>%
  mutate(leadingEdge =paste0(unlist(leadingEdge), collapse = "/"))

enrich_df_ctrl <- data.frame(ID = fgseaRes_ctrl$ID,
                        Description = fgseaRes_ctrl$Description,
                        setSize = fgseaRes_ctrl$size,
                        enrichmentScore = fgseaRes_ctrl$ES,
                        NES = fgseaRes_ctrl$NES,
                        pvalue = fgseaRes_ctrl$pval,
                        p.adjust = fgseaRes_ctrl$padj,
                        core_enrichment = fgseaRes_ctrl$leadingEdge)

rownames(enrich_df_ctrl) <- enrich_df_ctrl$ID

en_ctrl <- GseaVis::dfGO2gseaResult(enrich.df = enrich_df_ctrl,
                               geneList = FCgenelist_ctrl,
                               OrgDb = org.Rn.eg.db)

en_ctrl <- setReadable(x = en_ctrl, OrgDb = org.Rn.eg.db)
en_ctrl@geneSets$`GO:0000001` <- marker_gene_list$`GO:0000001|Inflammatory_Epi_markers`
enrichplot::gseaplot2(x = en_ctrl, geneSetID = enrich_df_ctrl$ID[1])

### get the position of Duox2 and Duoxa2
mygene <- c("79107", "499879")

# plot
pdf(file = "dpi_NEC_Ctrl_GSEA.pdf", width = 8, height = 4)
gseaNb(en_ctrl,
       geneSetID = enrich_df_ctrl$ID[1],
       addPval = T,
       newGsea = T,
       addGene = mygene
)
dev.off()



### Scatter plot
deg_result <- data.frame(gene = names(FCgenelist), log2FC = FCgenelist)
deg_result <- deg_result %>% dplyr::mutate(rank = -1 * rank(log2FC, ties.method = "max"))

p1 <- deg_result %>%
  ggplot() +
  geom_point(aes(x = rank, y = log2FC, color = log2FC, size = abs(log2FC))) +
  scale_x_continuous(breaks = c(-14022, -9022, -4022),
                     labels = c(0, 5000, 10000)) +
  scale_color_gradient2(low = '#F0E442', high = '#ff0000', mid = "#ffffff", midpoint = 0.5, name = "log2FC") +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(data = deg_result %>% filter(gene %in% c("Duox2", "Duoxa2")), aes(x = rank, y = log2FC), color = "black") +
  geom_text_repel(data = deg_result %>% filter(gene %in% c("Duox2", "Duoxa2")),
                  aes(x = rank + 10, y = log2FC, label = gene),
                  box.padding = 0.5, nudge_x = 10,
                  nudge_y = 0.2, segment.curvature = -0.1,
                  segment.ncp = 3, segment.angle = 20,
                  direction = "y", hjust = "left") +
  scale_size(range = c(1, 7), name = "log2FC") +
  labs(x = "Rank of DEGs",
       y = "log2FC", title = "NEC vs Ctrl") +
  theme_bw() +
  theme(
    axis.text = element_text(color = "#000000", size = 12.5),
    axis.title = element_text(color = "#000000", size = 15)
  )

deg_result_DPI <- data.frame(gene = names(FCgenelist_DPI), log2FC = FCgenelist_DPI)
deg_result_DPI <- deg_result_DPI %>% dplyr::mutate(rank = -1 * rank(log2FC, ties.method = "max"))

p2 <- deg_result_DPI %>%
  ggplot() +
  geom_point(aes(x = rank, y = log2FC, color = log2FC, size = abs(log2FC))) +
  scale_x_continuous(breaks = c(-13917, -8917, -3917),
                     labels = c(0, 5000, 10000)) +
  scale_color_gradient2(low = '#F0E442', high = '#ff0000', mid = "#ffffff", midpoint = 0.5, name = "log2FC") +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(data = deg_result_DPI %>% filter(gene %in% c("Duox2", "Duoxa2")), aes(x = rank, y = log2FC), color = "black") +
  geom_text_repel(data = deg_result_DPI %>% filter(gene %in% c("Duox2", "Duoxa2")),
                  aes(x = rank + 10, y = log2FC, label = gene),
                  box.padding = 0.5, nudge_x = 10,
                  nudge_y = 0.2, segment.curvature = -0.1,
                  segment.ncp = 3, segment.angle = 20,
                  direction = "y", hjust = "left") +
  scale_size(range = c(1, 7), name = "log2FC") +
  labs(x = "Rank of DEGs",
       y = "log2FC", title = "NEC+DPI vs NEC") +
  theme_bw() +
  theme(
    axis.text = element_text(color = "#000000", size = 12.5),
    axis.title = element_text(color = "#000000", size = 15)
  )


p1 + p2
ggsave("DEGs_rank_log2FC.pdf", width = 10, height = 6)


### Deconvolution analysis
load("../ExpressionSet_scRNA_Major.RData")
### count
load("../nec_dpi_adj_count.RData")
### ExpressionSet
bulk.eset <- Biobase::ExpressionSet(assayData = count)

load("markers_summary.RData")
markers_maj_use <- markers_maj %>% group_by(cluster) %>% top_n(n = 30, wt = avg_logFC)
res_maj <- BisqueRNA::ReferenceBasedDecomposition(bulk.eset, sc.eset_maj, markers = markers_maj_use, use.overlap = FALSE)

### plot
proportion_maj <- res_maj$bulk.props[c("Myeloid cell", "Neutrophil", "Epithelial cell"), ]
proportion_maj <- data.frame(type = rep(c("NEC-1", "NEC-2", "DPI-1", "DPI-2"), each = 3),
                             Major_clusters = c("Myeloid cell", "Neutrophil", "Epithelial cell"),
                             propor = c(proportion_maj[, 1], proportion_maj[, 2], proportion_maj[, 3], proportion_maj[, 4]))
proportion_maj$type_Maj <- paste0(proportion_maj$Major_clusters, "|", proportion_maj$type)
proportion_maj$type_Maj <- factor(proportion_maj$type_Maj, levels = rev(c("Myeloid cell|NEC-1", "Myeloid cell|NEC-2", "Myeloid cell|DPI-1", "Myeloid cell|DPI-2",
                                                                      "Neutrophil|NEC-1", "Neutrophil|NEC-2", "Neutrophil|DPI-1", "Neutrophil|DPI-2",
                                                                      "Epithelial cell|NEC-1", "Epithelial cell|NEC-2", "Epithelial cell|DPI-1", "Epithelial cell|DPI-2")))

ggplot(data = proportion_maj, aes(x = propor, y = type_Maj, fill = Major_clusters, color = Major_clusters)) +
  geom_point(size = 5) +
  scale_fill_manual(values = c("Myeloid cell" = "#395D97",
                               "Neutrophil" = "#53ACA6",
                               "Epithelial cell" = "#C91315"
)) +
  scale_color_manual(values = c("Myeloid cell" = "#395D97",
                               "Neutrophil" = "#53ACA6",
                               "Epithelial cell" = "#C91315"
)) +
  labs(title = NULL, x = 'Cell proportion', y = NULL) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 11, color = "black"),
        axis.title = element_text(size = 13)
  )
ggsave("Major_cellratio.pdf", width = 8, height = 4)