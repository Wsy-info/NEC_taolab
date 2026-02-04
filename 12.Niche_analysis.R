### Pie plot of Niche, co-localization analysis and Niche pathway enrichment analysis

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(org.Hs.eg.db)
library(clusterProfiler)



### load data
load("round2_NEC.harmony_singlet.RData")
load("all_st.RData")
load("niche_score.RData")

### get color for Major_clusters
color_cell_2 <- c("#F1A96E", ### B/Plasma cell
                  "#395D97", ### Myeloid cell
                  "#53ACA6", ### Neutrophil
                  "#E78382", ### Mast cell
                  "#865EA5", ### T/NK/ILC
                  "#C65FA4", ### Endothelial cell
                  "#C91315", ### Epithelial cell
                  "#A8B8D1", ### Stromal cell
                  "#7C3F3A", ### Glial cell
                  "#B4BC6E" ### Neuronal cell
)
names(color_cell_2) <- levels(NEC.final.harmony@meta.data$Major_clusters)

### get color for Minor_clusters
load("NEC_color.RData")

### pie plot
for(i in c("Normal", "Adjacent1", "NEC1")){
  df_pie_i <- all_st@meta.data %>% filter(Sample == i)
  for(j in as.character(unique(df_pie_i$Niche))){
    df_pie <- as.data.frame(table((all_st@meta.data %>% filter(Niche == j & Sample == i))$Minor_clusters))
    df_pie$Major_clusters <- rep(c("B/Plasma cell", "Myeloid cell", "Neutrophil", "Mast cell", "T/NK/ILC",
                                   "Epithelial cell", "Endothelial cell", "Stromal cell", "Glial cell", "Neuronal cell"),
                                 c(11, 14, 12, 1, 25, 15, 13, 22, 1, 1))
    df_pie <- df_pie %>% filter(Freq != 0 | Var1 == "no_label")
    df_pie$Major_clusters <- factor(df_pie$Major_clusters, levels = c("B/Plasma cell", "Myeloid cell", "Neutrophil", "Mast cell", "T/NK/ILC",
                                                                      "Epithelial cell", "Endothelial cell", "Stromal cell", "Glial cell", "Neuronal cell"))


    ggplot() +
      geom_col(data = df_pie %>% group_by(Major_clusters) %>% reframe(Freq = sum(Freq)),
               aes(x = 1, y = Freq, fill = Major_clusters), color = "white", width = 1.5) +
      geom_col(data = df_pie, aes(x = 2, y = Freq, fill = Var1), color = "white", width = 0.5) +
      scale_fill_manual(values = c(color_cell_2, all_color_2, add_color)) +
      coord_polar(theta = "y") +
      theme_void() +
      theme(legend.position = "none") +
      xlim(-0.5, 2.8)

    ggsave(paste0("pie/", i, "_", j, ".pdf"), width = 5, height = 5)
  }
}



### Co-localization analysis for IAF2
iaf2_niche_score <- niche_score %>% filter(Minor_clusters == "IAF2")
mean_score <- colMeans(iaf2_niche_score)

all(names(mean_score) == celltype_map_df$col_name)
mean_score_df <- data.frame(score = mean_score,
                            Major_clusters = celltype_map_df$Major_clusters,
                            Minor_clusters = celltype_map_df$Minor_clusters)

mean_score_df_use <- mean_score_df %>% filter(Major_clusters %in% c("B/Plasma cell", "T/NK/ILC", "Neutrophil", "Myeloid cell", "Mast cell"))
mean_score_df_use$Major_clusters <- factor(mean_score_df_use$Major_clusters, levels = c("B/Plasma cell", "Myeloid cell", "Neutrophil", "Mast cell", "T/NK/ILC"))

mean_score_df_use %>% group_by(Major_clusters) %>% reframe(mean(score))

ggplot(data = mean_score_df_use[6:68, ], aes(x = Major_clusters, y = score, fill = Major_clusters)) +
  geom_violin(scale = "width", adjust = 1.5, trim = F, alpha = 0.7) +
  stat_summary(fun = "mean", geom = "crossbar", width = 1, color = "black") +
  geom_quasirandom(width = 0.2, size = 1.2, alpha = 0.8) +
  scale_fill_manual(values = c("B/Plasma cell" = "#F1A96E",
                               "Myeloid cell" = "#395D97",
                               "Neutrophil" = "#53ACA6",
                               "Mast cell" = "#E78382",
                               "T/NK/ILC" = "#865EA5"
)) +
  theme_classic() +
  theme(axis.text = element_text(size = 12, color = "black"),
        legend.position = "none") +
  labs(x = "", y = "Enriched score of IAF2") +
  RotatedAxis()
ggsave("IAF2_Immune_cell_enriched_score_vlnplot2_new.pdf", width = 5, height = 6)


### pathway enrichment analysis for N9
all_st_log <- NormalizeData(all_st, normalization.method = "LogNormalize")
Idents(all_st_log) <- "Niche"

markers <- FindMarkers(all_st_log, ident.1 = "N9", only.pos = T)
N9_df <- bitr(rownames(markers),
                fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = "org.Hs.eg.db")
N9_pathway <- enrichGO(gene = unique(N9_df$ENTREZID),
                       OrgDb = org.Hs.eg.db,
                       keyType = 'ENTREZID',
                       ont = "BP", pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.1,
                       readable = TRUE)
N9_pathway <- N9_pathway@result %>% filter(pvalue < 0.05 & Count > 2)

n9_pathway_use <- N9_pathway %>% filter(Description %in% select_pathway_n9)
n9_pathway_use$Description <- factor(n9_pathway_use$Description, levels = rev(n9_pathway_use$Description))
ggplot(n9_pathway_use, aes(x = -log10(pvalue), y = Description)) +
  geom_bar(stat = "identity", color = "#e545d2", fill = "#e545d2", alpha = 0.5) +
  theme_classic() +
  theme(axis.text = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 14),
        axis.ticks = element_line(colour = "black")) +
  labs(x = '-log10(pvalue)', y = "")
ggsave(filename = "N9_pathway_barplot_new.pdf", width = 8, height = 5)