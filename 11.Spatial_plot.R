### Spatial data plot

library(Seurat)
library(dplyr)
library(ggplot2)

### NEC1
NEC1_all <- read.csv("NEC1_all.csv", row.names = 1)

### Major_clusters
NEC1_all_plot <- data.frame(cell = NEC1_all$cell_name,
                            x = NEC1_all$x + NEC1_all$xcoord,
                            y = NEC1_all$y + NEC1_all$ycoord,
                            ident = NEC1_all$Major_clusters,
                            boundary = "segmentation")
NEC1_all_plot$ident <- factor(NEC1_all_plot$ident)
SingleImagePlot(data = NEC1_all_plot, col.by = "ident", cols = "polychrome",
                shuffle.cols = F, size = 0.5, alpha = 1,
                mols.alpha = 1, border.color = "white",
                border.size = 0.1, dark.background = T)
ggsave("NEC1_all_final.pdf", width = 8, height = 8)
