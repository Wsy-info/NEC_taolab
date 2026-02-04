### Cell type ratio

library(Seurat)
library(dplyr)
library(ggplot2)
library(ggalluvial)
library(data.table)
library(Startrac)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)

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


### Cell number of Major_clusters
df <- as.data.frame(table(NEC.final.harmony@meta.data$Major_clusters))
df <- df %>% arrange(Freq)
colnames(df) <- c("Major_clusters", "Number")
df$Major_clusters <- as.character(df$Major_clusters)
df$Major_clusters <- factor(df$Major_clusters, levels = df$Major_clusters)
df$label <- paste0(df$Major_clusters, " (", scales::comma(df$Number), ")")

### barplot
ggplot(data = df, aes(x = Number, y = Major_clusters, fill = Major_clusters)) +
  geom_col() +
  geom_text(aes(label = label), hjust = -0.1, size = 4) +
  scale_fill_manual(values = color_cell_2) +
  scale_x_continuous(breaks = c(0, 10000, 20000, 30000, 40000), labels = c(0, 1, 2, 3, 4)) +
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12, color = "black"),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none"
  ) +
  labs(x = "Number of cells(n × 10,000)", y = "")
ggsave("Major_clusters_number_nec.pdf", width = 8, height = 8)


### Ratio
data_df <- NEC.final.harmony@meta.data
cellratio_df <- as.data.frame(prop.table(table(data_df$Major_clusters, data_df$Tissue), margin = 2))
colnames(cellratio_df) <- c("Cell_Type", "Tissue", "Ratio")

### plot
ggplot(cellratio_df, aes(x = Tissue, y = Ratio, fill = Cell_Type, alluvium = Cell_Type)) +
  geom_col(width = 0.5, color = NA) +
  geom_flow(width = 0.4, alpha = 0.2) +
  scale_fill_manual(values = color_cell_2) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(title = NULL, y = 'Cell proportion', x = NULL) +
  theme_classic() +
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 13, color = "black"),
        legend.position = "none")
ggsave("nec_major_cell_cellratio_barplot.pdf", width = 10, height = 9)


### Ro/e analysis
cellInfo.tb <- data_df[, c("Patient", "Tissue", "Sample", "Major_clusters")] %>% data.table()
cellInfo.tb$Major_clusters <- as.character(cellInfo.tb$Major_clusters)
cellInfo.tb$Tissue <- cellInfo.tb$Tissue %>% as.factor()
loc.avai.vec <- levels(cellInfo.tb[["Tissue"]])

### calculate Ro/e
startrac.dist <- unclass(calTissueDist(cellInfo.tb,
                                       byPatient = F,
                                       colname.cluster = "Major_clusters",
                                       colname.patient = "Patient",
                                       colname.tissue = "Tissue"))

startrac.dist <- startrac.dist[, loc.avai.vec]
cuts <- c(0, 0.0000001, 1, 1.5, 1.9, Inf)
startrac.dist.bin.values <- factor(c("-", "+/-", "+", "++", "+++"), levels = c("-", "+/-", "+", "++", "+++"))
startrac.dist.bin <- matrix(startrac.dist.bin.values[findInterval(startrac.dist, cuts)], ncol = ncol(startrac.dist))
colnames(startrac.dist.bin) <- colnames(startrac.dist)
rownames(startrac.dist.bin) <- rownames(startrac.dist.bin)

### plot
startrac.dist <- startrac.dist[levels(NEC.final.harmony@meta.data$Major_clusters), ]
getPalette <- colorRampPalette(brewer.pal(11, "Blues")[1:6])
pdf("nec.Roe_sig_blue.pdf", width = 5)
Heatmap(startrac.dist, name ="Ro/e",
        rect_gp = gpar(col = "white", lwd = 1.5),
        cluster_rows = F,
        cluster_columns = F,
        column_names_rot = 0,
        column_names_centered = TRUE,
        row_names_side = "left",
        row_title = NULL, column_title = NULL,
        column_names_gp = gpar(fontsize = 10),
        row_names_gp = gpar(fontsize = 10),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", startrac.dist[i, j]), x, y, gp = gpar(fontsize = 10))
        },
        col = getPalette(5),
        heatmap_legend_param = list(
          title = "Ro/e",
          break_dist = 1,
          col_fun = colorRamp2(c(0, 1, 1.5, 1.9, 4), getPalette(5)),
          at = c(0, 1, 1.5, 1.9, 4),
          labels = c('0','1','1.5','1.9','max')
 ))
dev.off()