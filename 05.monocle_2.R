### monocle2

library(Seurat)
library(dplyr)
library(ggplot2)
library(monocle)

### load data
load("round2_NEC.neu.harmony_singlet.RData")

### use count matrix
expr_matrix <- GetAssayData(NEC.neu.harmony, slot = "counts")
p_data <- NEC.neu.harmony@meta.data
f_data <- data.frame(gene_short_name = rownames(NEC.neu.harmony@assays$RNA),
                      row.names = row.names(NEC.neu.harmony@assays$RNA))
pd <- new('AnnotatedDataFrame', data = p_data)
fd <- new('AnnotatedDataFrame', data = f_data)

### create cds object
cds_pre <- newCellDataSet(expr_matrix,
                          phenoData = pd,
                          featureData = fd,
                          expressionFamily = negbinomial.size())

cds_pre <- estimateSizeFactors(cds_pre)
cds_pre <- estimateDispersions(cds_pre)

NEC.neu.harmony <- FindVariableFeatures(NEC.neu.harmony, nfeatures = 2000)
gene_sle <- VariableFeatures(NEC.neu.harmony)
cds <- setOrderingFilter(cds_pre, gene_sle)
plot_ordering_genes(cds)

### reduction
cds <- reduceDimension(cds, reduction_method = 'DDRTree',
                       max_components = 2)

### order
cds <- orderCells(cds, reverse = T)

### get information for plot
data_df <- t(reducedDimS(cds)) %>%
  as.data.frame() %>%
  dplyr::select('Component 1' = 1, 'Component 2' = 2) %>%
  rownames_to_column("Cells") %>%
  mutate(pData(cds)$State,
         pData(cds)$Pseudotime,
         pData(cds)$orig.ident,
         pData(cds)$Minor_clusters)
colnames(data_df) <- c("cells", "Component_1", "Component_2", "State",
                       "Pseudotime", "Tissue", "Cell_subtype")

reduced_dim_coords <- reducedDimK(cds)
ca_space_df <- Matrix::t(reduced_dim_coords) %>%
  as.data.frame() %>%
  dplyr::select(prin_graph_dim_1 = 1, prin_graph_dim_2 = 2) %>%
  mutate(sample_name = rownames(.), sample_state = rownames(.))
dp_mst <- minSpanningTree(cds)
edge_df <- dp_mst %>%
  igraph::as_data_frame() %>%
  dplyr::select(source = "from", target = "to") %>%
  left_join(ca_space_df %>% dplyr::select(source = "sample_name", source_prin_graph_dim_1 = "prin_graph_dim_1", source_prin_graph_dim_2 = "prin_graph_dim_2"), by = "source") %>%
  left_join(ca_space_df %>% dplyr::select(target = "sample_name", target_prin_graph_dim_1 = "prin_graph_dim_1", target_prin_graph_dim_2 = "prin_graph_dim_2"), by = "target")

save(edge_df, dp_mst, ca_space_df, reduced_dim_coords, data_df, file = "monocle2_neu_info_df.RData")


### Pseudotime plot
ggplot() +
  geom_segment(aes_string(x = "source_prin_graph_dim_1",
                          y = "source_prin_graph_dim_2",
                          xend = "target_prin_graph_dim_1",
                          yend = "target_prin_graph_dim_2"),
               linewidth = 0.75, linetype = "solid", na.rm = TRUE, data = edge_df) +
  geom_point(data = data_df, aes(x = Component_1,
                                 y = Component_2,
                                 color = Pseudotime)) +
  scale_color_distiller(palette = "YlGnBu", direction = 1) +
  theme_classic() +
  theme(legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.position = "bottom") +
  labs(x = "Componet 1", y = "Componet 2")
ggsave("monocle_pseu_neu.pdf")


### for Minor_clusters
ggplot() +
  geom_segment(aes_string(x = "source_prin_graph_dim_1",
                          y = "source_prin_graph_dim_2",
                          xend = "target_prin_graph_dim_1",
                          yend = "target_prin_graph_dim_2"),
               size = 0.75, linetype = "solid", na.rm = TRUE, data = edge_df) +
  geom_point(data = data_df, aes(x = Component_1,
                                 y = Component_2,
                                 color = Cell_subtype)) +
  scale_color_manual(name = "Cell_subtype",
                     values = color_neu_2) +
  theme_classic() +
  theme(legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.position = "bottom") +
  labs(x = "Componet 1", y = "Componet 2") +
  guides(color = guide_legend(override.aes = list(size = 5)))
ggsave("monocle_min_neu.pdf")


### ridge plot
rigde_df <- pData(cds)
ggplot(rigde_df, aes(x = Pseudotime, y = Minor_clusters,
                     fill = Minor_clusters, color = Minor_clusters)) +
  geom_density_ridges(rel_min_height = 0.001, alpha = 0.5) +
  stat_density_ridges(aes(color = Minor_clusters),
                      quantile_lines = TRUE,
                      quantiles = 2,
                      linewidth = 2,
                      alpha = 0.5) +
  scale_fill_manual(values = rev(color_neu_2)) +
  scale_color_manual(values = rev(color_neu_2)) +
  scale_y_discrete("") +
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size = 12)) +
  theme_classic()
ggsave("monocle_rigde_neu.pdf", width = 10, height = 6)


### plot the expression dynamics of marker genes
marker_plot <- plot_genes_in_pseudotime(cds[selected_genes, ],
                         color_by = "Minor_clusters",
                         nrow = 4,
                         ncol = 4,
                         cell_size = 1.2, trend_formula = "~ sm.ns(Pseudotime, df=4)") +
  geom_vline(xintercept = c(5, 20.4, 19.7)) +
  theme(legend.position = "none",
        plot.title = element_text(size = 15),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 13)
  )

marker_plot
ggsave(filename = "Markers_pseudotime_plot.pdf", width = 8, height = 6)
