### Cellchat

library(Seurat)
library(CellChat)
library(patchwork)
library(ComplexHeatmap)
library(dplyr)
library(ggplot2)
library(reshape2)
library(wesanderson)
library(circlize)

load("round2_NEC.harmony_singlet.RData")

data.input <- NEC.final.harmony@assays$RNA@data
meta <- NEC.final.harmony@meta.data
cellchat <- createCellChat(object = data.input,
                           meta = meta,
                           group.by = 'Minor_clusters')
cellchat <- setIdent(cellchat, ident.use = 'Minor_clusters')

CellChatDB <- CellChatDB.human
CellChatDB.use <- subsetDB(CellChatDB, search = 'Secreted Signaling')
cellchat@DB <- CellChatDB

cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat,
                                       only.pos = TRUE,
                                       thresh.pc = 0,
                                       thresh.fc = 0,
                                       thresh.p = 0.05)

cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat, slot.name = 'netP')
cellchat <- computeCommunProbPathway(cellchat)
df.netp <- subsetCommunication(cellchat, slot.name = 'netP')
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
groupSize <- as.numeric(table(cellchat@idents))

### save
save(cellchat, groupSize, file = "All_cell_chat_all.RData")

### LYVE1+ Mac, CD163L1+ Mac vs T/NK/ILC
netVisual_matrix <- netVisual_bubble(cellchat,
                                     sources.use = c(11, 12),
                                     targets.use = 36:58,
                                     angle.x = 45,
                                     remove.isolate = F
)
netVisual_matrix <- netVisual_matrix$data
netVisual_matrix <- netVisual_matrix %>% filter(pathway_name %in% c("MHC-II", "PD-L1", "CD80", "CD86"))
netVisual_matrix <- dcast(netVisual_matrix, interaction_name_2 ~ source.target, value.var = "prob", fill = 0)
rownames(netVisual_matrix) <- netVisual_matrix$interaction_name_2
netVisual_matrix$interaction_name_2 <- NULL
netVisual_matrix <- netVisual_matrix[c(1, 3, 5:16), ]


### ComplexHeatmap
max(netVisual_matrix)/2
col_fun <- colorRamp2(c(0, 0.1586494, 0.3172989), c("white", "#F2D287", "#e35487"))
col_label <- sapply(strsplit(colnames(netVisual_matrix), split = " -> "),"[", 2)

pdf("lyve1_cd163l1_tissue_mac_cd4t.pdf", width = 8, height = 8)
Heatmap(netVisual_matrix[, c(1, 9, 2, 10, 3, 11, 4, 12, 5, 13, 6, 14, 7, 15, 8, 16)],
        col = col_fun,
        cluster_rows = F,
        cluster_columns = F,
        column_title = NULL,
        column_gap = unit(0, "mm"),
        border = TRUE,
        rect_gp = gpar(col = "grey", lwd = 1),
        border_gp = gpar(col = "black", lwd = 0.8)
)
dev.off()