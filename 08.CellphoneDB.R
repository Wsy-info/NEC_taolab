### CellphoneDB result visualization

library(dplyr)
library(ggplot2)
library(pheatmap)
library(circlize)
library(fmsb)
library(scales)

load("round2_NEC.harmony_singlet.RData")

load("CellphoneDB_result.RData")

inter_heatmap <- function(source_index, target_index, compare_by, scale_lr){
  source_celltype <- levels(NEC.final.harmony@meta.data$Minor_clusters)[source_index]
  target_celltype <- levels(NEC.final.harmony@meta.data$Minor_clusters)[target_index]

  involve_pairs <- as.character(pairs$pair[which((pairs$celltype_a %in% source_celltype & pairs$celltype_b %in% target_celltype))])
  ### get all means of selected celltype
  means_wanted <- means[, c("interacting_pair", "gene_a", "gene_b", "receptor_a", "receptor_b", involve_pairs)]

  ### remove LR
  means_wanted <- means_wanted[which(!(means_wanted$gene_a == "" & means_wanted$gene_b == "")), ]

  ### get all pvalues of selected celltype
  pvalue_wanted <- pvalue[rownames(means_wanted), c("interacting_pair", "gene_a", "gene_b", "receptor_a", "receptor_b", involve_pairs)]

  ### only retain p<0.05
  for(i in 6:(dim(means_wanted)[2])){
    for(j in 1:(dim(means_wanted)[1])){
      if(pvalue_wanted[j, i] >= 0.05){
        means_wanted[j, i] <- 0
      }
    }
  }

  ### only retain LR with sum of means > 0
  means_wanted <- means_wanted[which(rowSums(means_wanted[, -c(1:5)]) > 0), ]
  rownames(means_wanted) <- means_wanted$interacting_pair

  ### only retain means
  means_wanted_new <- means_wanted[, -c(1:5)]
  means_wanted_new <- means_wanted_new[which(rowSums(means_wanted_new) > 0), ]

  if(scale_lr){
    ### z-score scale
    means_wanted_new_scale <- apply(means_wanted_new, 1, function(x){
      (x - min(x))/(max(x) - min(x))
    })
  } else {
    means_wanted_new_scale <- t(means_wanted_new)
  }

  involve_pairs_2 <- pairs[which((pairs$celltype_a %in% source_celltype & pairs$celltype_b %in% target_celltype)), ]
  rownames_order <- paste0(involve_pairs_2$celltype_a, "|", involve_pairs_2$celltype_b)
  follow_order <- c()
  if(compare_by == "target"){
    for(i in source_celltype){
      follow_order_i <- paste0(i, "|", target_celltype)
      follow_order <- c(follow_order, follow_order_i)
    }
  } else {
    if(compare_by == "source"){
      for(i in target_celltype){
      follow_order_i <- paste0(source_celltype, "|", i)
      follow_order <- c(follow_order, follow_order_i)
    }
    }
  }

  rownames(means_wanted_new_scale) <- rownames_order
  means_wanted_new_scale <- means_wanted_new_scale[follow_order, ]
  return(t(means_wanted_new_scale))
}


### Inflammatory Epi & Immune cell
matrix_1 <- inter_heatmap(67, 1:58, "target", T)

pheatmap::pheatmap(matrix_1,
                   scale = "none",
                   cluster_cols = F,
                   cluster_rows = F,
                   border_color = "grey",
                   gaps_col = c(10, 23, 34, 35),
                   color = colorRampPalette(c("white", "#F2D287", "#e35487"))(1000),
                   cellheight = 12,
                   cellwidth = 12,
                   filename = "cellphonedb_Inflammatory_epi_immune_cell.pdf"
)