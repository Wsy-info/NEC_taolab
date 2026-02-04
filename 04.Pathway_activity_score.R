### Pathway activity score

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

### load data
load("round2_NEC.harmony_singlet.RData")
load("round2_Hea.harmony_singlet.RData")

### merge
All.final.harmony <- merge(Normal.final.harmony, NEC.final.harmony)
All.final.harmony@meta.data$Tissue <- factor(All.final.harmony@meta.data$Tissue, levels = c("Normal", "Adjacent", "NEC"))

### get msigdb pathway
load("reference/msigdb/all_msigdbr_list.RData")

### color
color_tissue_2 <- c("#9195D9", "#C2BBC2", "#E66EAC")
comparsion <- list(c("Normal", "Adjacent"), c("Normal", "NEC"), c("Adjacent", "NEC"))


### calculate pathway activity acorss NEC, Adjacent, and Normal
### BoxPlot and ANOVA test

### 1.TLR
gene <- all_msigdbr_list[["REACTOME_TOLL_LIKE_RECEPTOR_CASCADES"]]
All.final.harmony <- AddModuleScore(All.final.harmony,
                                    features = list(gene),
                                    name = 'TLR4')

### 2.NLR
gene <- all_msigdbr_list[["KEGG_NOD_LIKE_RECEPTOR_SIGNALING_PATHWAY"]]
All.final.harmony <- AddModuleScore(All.final.harmony,
                                    features = list(gene),
                                    name = 'NLR')

### 3.inflammatory
gene <- all_msigdbr_list[["HALLMARK_INFLAMMATORY_RESPONSE"]]
All.final.harmony <- AddModuleScore(All.final.harmony,
                                    features = list(gene),
                                    name = 'inflammatory')

### 4. NFkB
gene <- all_msigdbr_list[["HINATA_NFKB_IMMU_INF"]]
All.final.harmony <- AddModuleScore(All.final.harmony,
                                          features = list(gene),
                                          name = 'NFKB')

### 5.ROS
gene <- all_msigdbr_list[["GOBP_RESPONSE_TO_REACTIVE_OXYGEN_SPECIES"]]
All.final.harmony <- AddModuleScore(All.final.harmony,
                                    features = list(gene),
                                    name = 'ROS')

### 6.Apoptosis
gene <- all_msigdbr_list[["HALLMARK_APOPTOSIS"]]
All.final.harmony <- AddModuleScore(All.final.harmony,
                                    features = list(gene),
                                    name = 'Apoptosis')

### 7.Negative proliferation
gene <- all_msigdbr_list[["GOBP_NEGATIVE_REGULATION_OF_CELL_POPULATION_PROLIFERATION"]]
All.final.harmony <- AddModuleScore(All.final.harmony,
                                    features = list(gene),
                                    name = 'Proliferation')

### 8.Cell death
gene <- all_msigdbr_list[["GOBP_POSITIVE_REGULATION_OF_CELL_DEATH"]]
All.final.harmony <- AddModuleScore(All.final.harmony,
                                    features = list(gene),
                                    name = 'Death')

### boxplot
p1 <- ggplot(data = All.final.harmony@meta.data, aes(x = Tissue, y = TLR41, color = Tissue, fill = Tissue)) +
  geom_boxplot(alpha = 0.3, linewidth = 1, outlier.color = NA) +
  theme_classic() +
  scale_color_manual(values = color_tissue_2) +
  scale_fill_manual(values = color_tissue_2) +
  labs(title = "TLR", x = "", y = 'Pathway score') +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 11),
        axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(size = 14),
        legend.position = "none")

p2 <- ggplot(data = All.final.harmony@meta.data, aes(x = Tissue, y = NLR1, color = Tissue, fill = Tissue)) +
  geom_boxplot(alpha = 0.3, linewidth = 1, outlier.color = NA) +
  theme_classic() +
  scale_color_manual(values = color_tissue_2) +
  scale_fill_manual(values = color_tissue_2) +
  labs(title = "NLR", x = "", y = 'Pathway score') +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 11),
        axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(size = 14),
        legend.position = "none")

p3 <- ggplot(data = All.final.harmony@meta.data, aes(x = Tissue, y = inflammatory1, color = Tissue, fill = Tissue)) +
  geom_boxplot(alpha = 0.3, linewidth = 1, outlier.color = NA) +
  theme_classic() +
  scale_color_manual(values = color_tissue_2) +
  scale_fill_manual(values = color_tissue_2) +
  labs(title = "Inflammatory", x = "", y = 'Pathway score') +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 11),
        axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(size = 14),
        legend.position = "none")

p4 <- ggplot(data = All.final.harmony@meta.data, aes(x = Tissue, y = Proliferation1, color = Tissue, fill = Tissue)) +
  geom_boxplot(alpha = 0.3, linewidth = 1, outlier.color = NA) +
  theme_classic() +
  scale_color_manual(values = color_tissue_2) +
  scale_fill_manual(values = color_tissue_2) +
  labs(title = "Negative reg.of proliferation", x = "", y = 'Pathway score') +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 11),
        axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(size = 14),
        legend.position = "none")

p5 <- ggplot(data = All.final.harmony@meta.data, aes(x = Tissue, y = NFKB1, color = Tissue, fill = Tissue)) +
  geom_boxplot(alpha = 0.3, linewidth = 1, outlier.color = NA) +
  theme_classic() +
  scale_color_manual(values = color_tissue_2) +
  scale_fill_manual(values = color_tissue_2) +
  labs(title = "NF-kB", x = "", y = 'Pathway score') +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 11),
        axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(size = 14),
        legend.position = "none")

p6 <- ggplot(data = All.final.harmony@meta.data, aes(x = Tissue, y = Death1, color = Tissue, fill = Tissue)) +
  geom_boxplot(alpha = 0.3, linewidth = 1, outlier.color = NA) +
  theme_classic() +
  scale_color_manual(values = color_tissue_2) +
  scale_fill_manual(values = color_tissue_2) +
  labs(title = "Cell death", x = "", y = 'Pathway score') +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 11),
        axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(size = 14),
        legend.position = "none")

p7 <- ggplot(data = All.final.harmony@meta.data, aes(x = Tissue, y = Apoptosis1, color = Tissue, fill = Tissue)) +
  geom_boxplot(alpha = 0.3, linewidth = 1, outlier.color = NA) +
  theme_classic() +
  scale_color_manual(values = color_tissue_2) +
  scale_fill_manual(values = color_tissue_2) +
  labs(title = "Apoptosis", x = "", y = 'Pathway score') +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 11),
        axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(size = 14),
        legend.position = "none")

p8 <- ggplot(data = All.final.harmony@meta.data, aes(x = Tissue, y = ROS1, color = Tissue, fill = Tissue)) +
  geom_boxplot(alpha = 0.3, linewidth = 1, outlier.color = NA) +
  theme_classic() +
  scale_color_manual(values = color_tissue_2) +
  scale_fill_manual(values = color_tissue_2) +
  labs(title = "Reactive oxygen species", x = "", y = 'Pathway score') +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 11),
        axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(size = 14),
        legend.position = "none")

p1 + p2 + p3 + p5 + p8 + p7 + p4 + p6 + plot_layout(ncol = 4)
ggsave("all_pathway_8.pdf", width = 12, height = 6)

### ANOVA test
df <- All.final.harmony@meta.data[, c("TLR41", "Tissue")]
anova_model <- aov(TLR41 ~ Tissue, data = df)
bartlett.test(TLR41 ~ Tissue, data = df)
summary(anova_model)
TukeyHSD(anova_model)

df <- All.final.harmony@meta.data[, c("NLR1", "Tissue")]
anova_model <- aov(NLR1 ~ Tissue, data = df)
bartlett.test(NLR1 ~ Tissue, data = df)
summary(anova_model)
TukeyHSD(anova_model)

df <- All.final.harmony@meta.data[, c("inflammatory1", "Tissue")]
anova_model <- aov(inflammatory1 ~ Tissue, data = df)
bartlett.test(inflammatory1 ~ Tissue, data = df)
summary(anova_model)
TukeyHSD(anova_model)

df <- All.final.harmony@meta.data[, c("Proliferation1", "Tissue")]
anova_model <- aov(Proliferation1 ~ Tissue, data = df)
bartlett.test(Proliferation1 ~ Tissue, data = df)
summary(anova_model)
TukeyHSD(anova_model)

df <- All.final.harmony@meta.data[, c("Death1", "Tissue")]
anova_model <- aov(Death1 ~ Tissue, data = df)
bartlett.test(Death1 ~ Tissue, data = df)
summary(anova_model)
TukeyHSD(anova_model)

df <- All.final.harmony@meta.data[, c("Apoptosis1", "Tissue")]
anova_model <- aov(Apoptosis1 ~ Tissue, data = df)
bartlett.test(Apoptosis1 ~ Tissue, data = df)
summary(anova_model)
TukeyHSD(anova_model)

df <- All.final.harmony@meta.data[, c("ROS1", "Tissue")]
anova_model <- aov(ROS1 ~ Tissue, data = df)
bartlett.test(ROS1 ~ Tissue, data = df)
summary(anova_model)
TukeyHSD(anova_model)

df <- All.final.harmony@meta.data[, c("NFKB1", "Tissue")]
anova_model <- aov(NFKB1 ~ Tissue, data = df)
bartlett.test(NFKB1 ~ Tissue, data = df)
summary(anova_model)
TukeyHSD(anova_model)


### Heatmap for pathway activity score
### Neutrophil
load("round2_NEC.neu.harmony_singlet.RData")
all_score_mean <- NEC.neu.harmony@meta.data %>%
  group_by(Minor_clusters) %>%
  reframe("Primary granules" = mean(Pri1),
          "Secondary granules" = mean(Sec1),
          "Teriary granules" = mean(Ter1),
          "Secretory vesicles" = mean(Sv1),
          "Activation" = mean(Act1),
          "Maturation" = mean(Mat1),
          "Neutrophil aging" = mean(Aging1),
          "IFN-gamma" = mean(IFN_gamma1),
          "Leukocyte chemotaxis" = mean(Che1),
          "Chemokine activity" = mean(Chem1),
          "Glycolysis" = mean(Gly1),
          "Hypoxia" = mean(Hyp1),
          "T cell migration" = mean(T1),
          "Response to stimulus"= mean(stress1)
) %>% as.data.frame()


rownames(all_score_mean) <- all_score_mean$Minor_clusters
all_score_mean <- all_score_mean[, -1]

### plot
pheatmap::pheatmap(all_score_mean,
                   cluster_rows = F,
                   cluster_cols = F,
                   cellwidth = 28,
                   cellheight = 28,
                   border_color = NA,
                   scale = "column",
                   color = colorRampPalette(c("#2d6bb4", "white", "#e6251e"))(50),
                   # filename = "neu_heatmap.pdf"
)


### Lollipop plot for pathway activity score
### Stromal cell
load("round2_NEC.str.harmony_singlet.RData")
NEC.str.harmony.filter <- subset(NEC.str.harmony, subset = Minor_clusters %in% c("Stromal 3(C7+)", "Stromal 4", "IAF1", "IAF2", "GREM1+ Fib"))

pathway_score <- NEC.str.harmony.filter@meta.data[, 20:23]
all_means <- colMeans(pathway_score[, -1])

p <- list()
for(i in names(all_means)){
  min <- min(all_score[which(all_score$Pathway == i), "Score"])
  max <- max(all_score[which(all_score$Pathway == i), "Score"])
  p[[i]] <- ggplot(data = all_score[which(all_score$Pathway == i), ],
       aes(x = Score, y = Minor_clusters, color = Minor_clusters)) +
    geom_vline(xintercept = all_means[i], color = "black", linetype = 2, linewidth = 0.5) +
    geom_point(size = 6) +
    geom_hline(yintercept = 1, color = color_str_2[4], linewidth = 1) +
    geom_hline(yintercept = 2, color = color_str_2[5], linewidth = 1) +
    geom_hline(yintercept = 3, color = color_str_2[6], linewidth = 1) +
    geom_hline(yintercept = 4, color = color_str_2[7], linewidth = 1) +
    geom_hline(yintercept = 5, color = color_str_2[11], linewidth = 1) +
    scale_color_manual(values = color_str_2) +
    scale_x_continuous(limits = c(min - (max - min) * 0.1, max + (max - min) * 0.1)) +
    theme_bw() +
    theme(axis.text.y = element_blank(),
          axis.text.x = element_text(size = 12, color = "black"),
          panel.grid = element_blank()) +
    labs(x = "", y = "", title = i)
}

p[[1]] + p[[2]] + p[[3]] + p[[4]] + plot_layout(nrow = 1, guides = "collect")
ggsave("str_pathway.pdf", width = 12, height = 3)
