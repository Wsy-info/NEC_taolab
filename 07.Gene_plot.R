### Violin plot and heatmap plot of gene expression

library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)

### load data
load("round2_NEC.epi.harmony_singlet.RData")
load("round2_NEC.tcell.harmony_singlet.RData")


### violin plot
features <- c('DUOX2', 'DUOXA2', 'S100A9',  'S100A11', 'S100A12',
              'S100A16', 'SAA1', 'SAA2', 'ICAM1', 'ITGA2',
              'ASS1', 'STOM', 'LCN2', 'NOS2', 'PLAU',
              'ANXA2', 'SLC6A14', 'WARS', 'XDH', 'RAC1',
              'MMP1', 'MMP10', 'CCL20', 'CXCL1', 'CXCL3',
              'CXCL5', 'CXCL11', 'IL1B', 'TNFSF10', 'FGB',
              'FGA', 'REG1B', 'HMOX1', 'RETREG1'
)

vln.data <- VlnPlot(NEC.epi.harmony, features = features[1])$data
for(i in features[-1]){
  vln.data_i <- VlnPlot(NEC.epi.harmony, features = i)$data
  vln.data$new <- vln.data_i[, 1]
  colnames(vln.data)[length(colnames(vln.data))] <- i
}


vln.dat <- FetchData(NEC.epi.harmony, c(features, "Minor_clusters"), slot = "data")
vln.dat$Cell <- rownames(vln.dat)
vln.dat.melt <- reshape2::melt(vln.dat, id.vars = c("Cell", "Minor_clusters"),
                               measure.vars = features,
                               variable.name = "gene",
                               value.name = "Expr")

### violin plot
ggplot(vln.dat.melt, aes(Minor_clusters, Expr, fill = Minor_clusters)) +
  geom_violin(scale = "width", adjust = 1, trim = TRUE) +
  stat_summary(aes(group = vars(gene)), fun = mean, geom = "crossbar", size = 0.5, color = "black", position = position_dodge(0.9)) +
  scale_y_continuous(expand = c(0, 0), position = "right", labels = function(x)
    c(rep(x = "", times = length(x) - 2), x[length(x) - 1], "")) +
  facet_grid(rows = vars(gene), scales = "free", switch = "y") +
  scale_fill_manual(values = color_epi_2) +
  theme_cowplot(font_size = 12) +
  theme(legend.position = "none", panel.spacing = unit(0, "lines"),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = NA, color = "black"),
        plot.margin = margin(7, 7, 0, 7, "pt"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        strip.text.y.left = element_text(angle = 0),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black")
  )
ggsave("gene_violinplot_epi.pdf", width = 8, height = 18)


### heatmap
geneset_feature <- c("CCR7", "LEF1", "SELL", "TCF7", ### Naive markers
                     "CD69", "CXCR6", "NR4A1", "NR4A3", "RUNX3",### Resident
                     "CD274", "CTLA4", "HAVCR2", "LAG3", "PDCD1", "TIGIT", ### Inhibitory receptors
                     "BHLHE40", "CCL3", "CCL4", "CCL5", "CST7", "EOMES", "GNLY", "GZMA", "GZMB", "GZMK", "GZMM", "IFNG", "KLRK1", "NKG7", "PRF1", ### Cytotoxicity and Cytokines
                     "CD28", "ICOS", "TNFRSF4", "TNFRSF9", "TNFRSF18", ### Co-stimulatory molecules
                     "IL2RA", "IL2RB", "IL2RG", ### IL2R
                     "FOS","FOSB","JUN","JUNB", "TBX21", "ZEB2" ### Transcriptional factor(TF)
)

expression_geneset_feature <- AverageExpression(NEC.tcell.harmony, features = (geneset_feature))$RNA

### heatmap plot
pheatmap::pheatmap(expression_geneset_feature,
                   cluster_rows = F,
                   cluster_cols = F,
                   border_color = "white",
                   scale = "row",
                   gaps_row = c(4, 9, 15, 30, 35, 38),
                   cellwidth = 10, cellheight = 10,
                   # filename = "tcell_feature_heatmap.pdf",
                   width = 9, height = 10.5,
                   color = colorRampPalette(c("#2d6bb4", "white", "#e6251e"))(50)
)
