### UMAP, FeaturePlot, DotPlot for scRNA-seq

library(Seurat)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(scRNAtoolVis)

### load data
load("round2_NEC.harmony_singlet.RData")

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

### get color for Tissue
color_tissue_2 <- c("#9195D9", "#C2BBC2", "#E66EAC")
names(color_tissue_2) <- c("Normal", "Adjacent", "NEC")


### UMAP for Tissue
Idents(NEC.final.harmony) <- "Tissue"
umap_data <- DimPlot(NEC.final.harmony)$data

ggplot(data = umap_data, aes(x = UMAP_1, y = UMAP_2, color = ident)) +
  geom_point(size = 0.5, stroke = 0.1) +
  scale_color_manual(values = color_tissue_2) +
  facet_wrap(~ident) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank()) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none") +
  labs(x = "", y = "")
ggsave("nec_tissue_cell_umap.pdf", width = 12, height = 6)


### UMAP for Minor_clusters
Idents(NEC.final.harmony) <- "Minor_clusters"
umap_data <- DimPlot(NEC.final.harmony)$data
ggplot(data = umap_data, aes(x = UMAP_1, y = UMAP_2, color = ident)) +
  geom_point(size = 0.5, stroke = 0.1) +
  scale_color_manual(values = all_color_2) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank()) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none") +
  labs(x = "", y = "")
ggsave("nec_minor_cell_umap.pdf", width = 6, height = 6)


### FeaturePlot for Epithelial cells
load("round2_NEC.epi.harmony_singlet.RData")
### get colors
reds_palette <- brewer.pal(n = 9, name = "Reds")
blue_palette <- brewer.pal(n = 9, name = "Blues")

### plot
p1 <- FeaturePlot(NEC.epi.harmony, features = c("DUOX2")) +
  scale_color_gradientn(colors = c("lightgrey", reds_palette[3:9]))
p2 <- FeaturePlot(NEC.epi.harmony, features = c("DUOXA2")) +
  scale_color_gradientn(colors = c("lightgrey", reds_palette[3:9]))
p1 + p2
ggsave("Inflammatory_epithelial_cell_marker_featureplot.pdf", width = 5, height = 10)

p3 <- FeaturePlot(NEC.epi.harmony, features = c("FGA")) +
  scale_color_gradientn(colors = c("lightgrey", blue_palette[c(3:4, 7, 9)]))
p4 <- FeaturePlot(NEC.epi.harmony, features = c("FGB")) +
  scale_color_gradientn(colors = c("lightgrey", blue_palette[c(3:4, 7, 9)]))
p3 + p4
ggsave("FGB_enterocyte_marker_featureplot.pdf", width = 5, height = 10)


### DotPlot for Epithelial cells
pdf("dotplot/epi.dotplot.pdf", height = 10.5, width = 8)
jjDotPlot(object = NEC.epi.harmony,
          id = "Minor_clusters",
          x.text.angle = 45,
          x.text.vjust = 1,
          gene = gene_show,
          ytree = F,
          dot.max = 6,
          dot.col = c("#2166ac", "white", "#b2182b"),
          anno = T) +
  coord_flip() +
  theme(axis.text.x = element_text(vjust = 0.5)) +
  scale_size(limits = c(0, 100), range = c(1, 6), breaks = c(0, 25, 50, 75, 100))
dev.off()