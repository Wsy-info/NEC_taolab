### for scRNA-seq

library(Seurat)
library(dplyr)
library(ggplot2)
library(harmony)

### Sample 1
nec_data_s1 <- Read10X("use_count/Sample_1_result/outs/filtered_feature_bc_matrix")
colnames(nec_data_s1) <- paste0("Nec_s1_",colnames(nec_data_s1))

doublet <- read.table("s1_Scrublet_result.txt", sep = "\t", header = T)
nec_data_s1 <- nec_data_s1[, !(colnames(nec_data_s1) %in% doublet$CB)]

### Create Seurat Object and QC
nec_data_s1 <- CreateSeuratObject(counts = nec_data_s1, min.cell = 3, min.features = 200, project = "sample_1")

nec_data_s1[["percent.mt"]] <- PercentageFeatureSet(nec_data_s1, pattern = "^MT-")
nec_data_s1[["percent.ribo"]] <- PercentageFeatureSet(nec_data_s1, pattern = "^RPS|^RPL")
nec_data_s1[["percent.hb"]] <- PercentageFeatureSet(nec_data_s1, pattern = "^HBA|^HBB")
head(nec_data_s1@meta.data,5)

nec_data_s1@meta.data$Patient <- "Patient1"
nec_data_s1@meta.data$Tissue <- "NEC"
nec_data_s1@meta.data$Sample <- "Patient1_NEC"

### qc_std_plot
nFeature_lower <- 500
nFeature_upper <- 6000
nCount_lower <- 1000
nCount_upper <- 100000
pMT_lower <- 0
pMT_upper <- 50
pHB_lower <- 0
pHB_upper <- 5

nec_data_s1 <- subset(nec_data_s1, subset = nFeature_RNA < nFeature_upper &
  nFeature_RNA > nFeature_lower &
  nCount_RNA < nCount_upper &
  nCount_RNA > nCount_lower &
  percent.mt < pMT_upper &
  percent.hb < pHB_upper)

### similar code for other samples

### merge seurat object
data <- Reduce(merge, list)
dim(data)

### Downstream Analysis
NEC.harmony <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000)
NEC.harmony <- FindVariableFeatures(NEC.harmony, selection.method = "vst", nfeatures = 3000)

mito.genes <- grep(pattern = "^MT-", rownames(NEC.harmony), value = TRUE)
ribo.genes <- grep(pattern = "^RP[LS]", rownames(NEC.harmony), value = TRUE)
NEC.harmony <- NEC.harmony[setdiff(rownames(NEC.harmony), mito.genes), ]
NEC.harmony <- NEC.harmony[setdiff(rownames(NEC.harmony), ribo.genes), ]

### regress out cell cycle
s.genes <- Seurat::cc.genes.updated.2019$s.genes
g2m.genes <- Seurat::cc.genes.updated.2019$g2m.genes
NEC.harmony <- CellCycleScoring(NEC.harmony, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
NEC.harmony <- ScaleData(NEC.harmony, vars.to.regress = c("S.Score", "G2M.Score"))
NEC.harmony <- RunPCA(NEC.harmony, features = NEC.harmony@assays$RNA@var.features)

### Remove batch effect
NEC.harmony <- RunHarmony(NEC.harmony, max.iter.harmony = 50, group.by.vars = "Patient")
NEC.harmony <- FindNeighbors(NEC.harmony, reduction = "harmony", dims = 1:30) %>% FindClusters(resolution = 1.2)
NEC.harmony <- RunUMAP(NEC.harmony, reduction = "harmony", dims = 1:30, check_duplicates = FALSE)
NEC.harmony <- RunTSNE(NEC.harmony, reduction = "harmony", dims = 1:30, check_duplicates = FALSE)