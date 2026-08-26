# =============================================================================
# ICH Brain scRNA-seq Analysis: Microglia/Myeloid Cell Clustering
# Mouse 10x Chromium | Female vs. Male Comparison
# =============================================================================
#
# Pipeline Overview:
#   1.  Install & Load Packages
#   2.  Load Raw 10X Data
#   3.  Create Seurat Objects & Add Metadata
#   4.  Calculate Mitochondrial Percentage
#   5.  QC Violin Plots (Pre-Filter)
#   6.  Preview Cell Loss at Chosen Cutoffs
#   7.  Apply Quality Filters
#   8.  Merge & Save Pre-Integration Objects
#   9.  Normalization & Preprocessing
#   10. Clustering
#   11. QC Metrics (Post-Clustering)
#   12. Cell-Type Marker Visualization
#   13. Cluster Annotation (All Cell Types)
#   14. Subset Myeloid/Microglial Clusters
#   15. Sex-Balanced Downsampling
#   16. Re-label Clusters on Balanced Object
#   17. UMAP Plots (Balanced Object)
#   18. QC Plots (Balanced Object)
#   19. Gene Lists
#   20. DotPlots
#   21. Cell Composition Barplots
#   22. Differential Expression & Heatmap
#   23. Export for Morpheus Heatmap Tool
#   24. Pheatmap from Morpheus Matrix
#   25. Cluster-Specific Marker Genes (per Cluster)
#   26. Female vs. Male DGE per Cluster
#   27. Female vs. Male DGE Across All Microglia (Combined)
#   28. Save Final Object
#   29. Module Scoring: Gene Signatures
#   30. Add All Module Scores
#   31. Filter & Clean Module Labels
#   32. Clip, Rescale & Plot Selected Modules
#   33. Save Module Score UMAP to PDF
#   34. Hallmark GSEA: Female vs. Male per Myeloid Cluster
# =============================================================================


# ── 1. Install & Load Required Packages ──────────────────────────────────────

library(dplyr)
library(ggplot2)
library(patchwork)
library(viridis)
library(openxlsx)
library(tidyr)
library(pheatmap)
library(limma)
library(dplyr)
library(ggplot2)
library(patchwork)
library(Seurat)
library(SeuratData)
library(sctransform)
library(HGNChelper)
library(multtest)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)


# ── 2. Load Raw 10X Data ──────────────────────────────────────────────────────
# Update paths below to point to your local filtered_feature_bc_matrix
# directories. $`Gene Expression` is specified because CITE-seq data is used;
# this ensures only the transcriptomic modality is passed to Seurat.

F.Brain <- Read10X(data.dir = "data/F_Brain/filtered_feature_bc_matrix")
M.Brain <- Read10X(data.dir = "data/M_Brain/filtered_feature_bc_matrix")


# ── 3. Create Seurat Objects & Add Metadata ───────────────────────────────────

FBRAIN <- CreateSeuratObject(
  counts       = F.Brain$`Gene Expression`,
  project      = "FBrain",
  min.cells    = 3,
  min.features = 200
)
FBRAIN$tissue <- "Brain"
FBRAIN$sex    <- "Female"

MBRAIN <- CreateSeuratObject(
  counts       = M.Brain$`Gene Expression`,
  project      = "MBrain",
  min.cells    = 3,
  min.features = 200
)
MBRAIN$tissue <- "Brain"
MBRAIN$sex    <- "Male"


# ── 4. Calculate Mitochondrial Percentage ─────────────────────────────────────
# Adds a 'percent.mt' column: the fraction of counts from mitochondrial genes
# (genes starting with "mt-").

FBRAIN[["percent.mt"]] <- PercentageFeatureSet(FBRAIN, pattern = "^mt-")
MBRAIN[["percent.mt"]] <- PercentageFeatureSet(MBRAIN, pattern = "^mt-")


# ── 5. QC Violin Plots (Pre-Filter) ───────────────────────────────────────────
# Visualize nFeatures, nCounts, and percent.mt before filtering to inform
# cutoff decisions.

sample_list <- list(FBRAIN, MBRAIN)

for (x in sample_list) {
  print(
    VlnPlot(x, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  )
}


# ── 6. Preview Cell Loss at Chosen Cutoffs ────────────────────────────────────
# Evaluate how many cells would be retained under the cutoffs:
#   200 < nFeature_RNA < 6000  &  percent.mt < 5
# Run this before subsetting to avoid having to re-load data.

pre  <- vector()
post <- vector()

for (x in sample_list) {
  pre  <- append(pre,  table(x@active.ident))
  y    <- subset(x, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 5)
  post <- append(post, table(y@active.ident))
}

cell_loss <- merge(
  data.frame(pre),
  data.frame(post),
  by = "row.names"
)
cell_loss$percent_lost <- ((cell_loss$pre - cell_loss$post) / cell_loss$pre) * 100

print(cell_loss)

# Optionally save the cell-loss summary:
# write.csv(cell_loss, file = "outputs/CellLossAfterFiltering.csv", row.names = FALSE)


# ── 7. Apply Quality Filters ──────────────────────────────────────────────────
# Adjust cutoffs here if the cell-loss preview (section 6) suggests a
# different threshold is more appropriate for your data.

FBRAIN <- subset(FBRAIN, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 5)
MBRAIN <- subset(MBRAIN, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 5)


# ── 8. Merge & Save Pre-Integration Objects ───────────────────────────────────
# Merge both samples into a single Seurat object for integration.
# Individual filtered objects are also saved for reference.

brains <- merge(FBRAIN, y = MBRAIN)

saveRDS(FBRAIN,  file = "outputs/ICHBrains_F.rds")
saveRDS(MBRAIN,  file = "outputs/ICHBrains_M.rds")
saveRDS(brains,  file = "outputs/ICHBrains.rds")

# Free memory from individual sample objects
rm(FBRAIN, MBRAIN, F.Brain, M.Brain, sample_list, cell_loss, pre, post); gc()


# ── 9. Normalization & Preprocessing ─────────────────────────────────────────
# SCTransform regresses out mitochondrial percentage and replaces the standard
# NormalizeData → ScaleData workflow. CCA integration corrects for batch
# effects between Female and Male samples.

obj <- brains

obj <- SCTransform(obj, vars.to.regress = "percent.mt", verbose = FALSE)
obj <- RunPCA(obj, verbose = FALSE)
obj <- IntegrateLayers(
  object               = obj,
  method               = CCAIntegration,
  normalization.method = "SCT",
  verbose              = FALSE
)

# saveRDS(obj, file = "outputs/ICHIntegratedMandFBrains.rds")

rm(brains); gc()


# ── 10. Clustering ────────────────────────────────────────────────────────────

ktest <- FindNeighbors(obj, reduction = "integrated.dr", dims = 1:20)
ktest <- FindClusters(ktest, resolution = 0.2, algorithm = 1)
ktest <- RunUMAP(ktest, dims = 1:20, reduction = "integrated.dr")

DefaultAssay(ktest) <- "RNA"
ktest <- JoinLayers(ktest)

rm(obj); gc()


# ── 11. QC Metrics (Post-Clustering) ─────────────────────────────────────────

p1 <- VlnPlot(ktest, features = "percent.mt",   group.by = "orig.ident")
p2 <- VlnPlot(ktest, features = "nFeature_RNA", group.by = "orig.ident")
p3 <- VlnPlot(ktest, features = "nCount_RNA",   group.by = "orig.ident")

# ggsave("outputs/qc_percent_mt.png",   plot = p1, dpi = 300, width = 10, height = 6)
# ggsave("outputs/qc_nFeature_RNA.png", plot = p2, dpi = 300, width = 10, height = 6)
# ggsave("outputs/qc_nCount_RNA.png",   plot = p3, dpi = 300, width = 10, height = 6)


# ── 12. Cell-Type Marker Visualization ───────────────────────────────────────

DefaultAssay(ktest) <- "RNA"

cell_type_markers <- list(
  Microglia        = c("Cx3cr1", "Sall1", "P2ry12", "Il18", "Cd74"),
  Mono_Macs        = c("Cd163", "Cd14", "Ptprc", "Cd86"),
  Dendritic_Cells  = c("Itgax", "Ptcra", "Krt71", "Ttc24", "Lilra4"),
  T_Cells          = c("Cd3e", "Cd3d", "Cd3g", "Lck", "Il7r", "Cd8a", "Trac"),
  B_Cells          = c("Cd19", "Jchain", "Ms4a1", "Igha", "Ighm", "Ighg", "Ighd"),
  Astrocytes       = c("Aldh1l1", "Gfap", "Aqp4", "Fgfr3"),
  Oligodendrocytes = c("Plp1", "Mbp", "Cnp", "Mog", "Olig1", "Olig2", "Mobp")
)

for (cell_type in names(cell_type_markers)) {
  print(VlnPlot(ktest, features = cell_type_markers[[cell_type]]))
}

# Feature plots for canonical markers
canonical_markers <- c("Cx3cr1", "P2ry12", "Fcrls", "Cd163", "Cd3e", "Ighm", "Aqp4", "Plp1")

for (gene in canonical_markers) {
  print(FeaturePlot(ktest, features = gene))
}


# ── 13. Cluster Annotation (All Cell Types) ───────────────────────────────────
# Assign biologically meaningful labels to each Seurat cluster number.
# Update new_cluster_ids if your clustering resolution produces different
# cluster numbers.

new_cluster_ids <- c(
  "Micro 1", "Micro 2",  "Micro 3",  "Micro 4",  "Micro 5",  "Mono Macs",
  "Micro 6", "B-Cells",  "Micro 7",  "Micro 8",  "Astrocytes", "Micro 9",
  "T-Cells", "Oligodendrocytes", "Micro 10"
)
names(new_cluster_ids) <- levels(ktest)
ktest <- RenameIdents(ktest, new_cluster_ids)

DimPlot(ktest, reduction = "umap", label = TRUE)
DimPlot(ktest, reduction = "umap", group.by = "orig.ident")
DimPlot(ktest, reduction = "umap", split.by = "orig.ident", label = TRUE)


# ── 14. Subset Myeloid/Microglial Clusters ───────────────────────────────────

myeloid_clusters <- c(
  "Micro 1", "Micro 2", "Micro 3", "Micro 4",  "Micro 5",
  "Mono Macs", "Micro 6", "Micro 7", "Micro 8", "Micro 9", "Micro 10"
)

sub_obj <- subset(ktest, idents = myeloid_clusters)
DimPlot(sub_obj, reduction = "umap", label = TRUE)


# ── 15. Sex-Balanced Downsampling ─────────────────────────────────────────────
# Downsample male cells proportionally to match female cell count per cluster,
# preserving cluster composition across both sexes.

table(sub_obj$orig.ident)

Idents(sub_obj) <- "seurat_clusters"

F_obj <- subset(sub_obj, subset = orig.ident == "FBrain")
M_obj <- subset(sub_obj, subset = orig.ident == "MBrain")

n_F <- nrow(F_obj@meta.data)

Idents(M_obj) <- "seurat_clusters"

counts_F <- table(F_obj$seurat_clusters)
counts_M <- table(M_obj$seurat_clusters)

common_clusters   <- intersect(names(counts_F[counts_F > 0]), names(counts_M[counts_M > 0]))
cell_counts_M     <- counts_M[common_clusters]
props             <- cell_counts_M / sum(cell_counts_M)
cells_per_cluster <- round(props * n_F)

set.seed(123)
M_keep <- unlist(lapply(names(cells_per_cluster), function(cl) {
  available <- WhichCells(M_obj, idents = cl)
  sample(available, size = min(length(available), cells_per_cluster[cl]))
}))

F_cells      <- WhichCells(sub_obj, expression = orig.ident == "FBrain")
balanced_obj <- subset(sub_obj, cells = c(F_cells, M_keep))

table(balanced_obj$orig.ident)

rm(F_obj, M_obj, ktest); gc()


# ── 16. Re-label Clusters on Balanced Object ──────────────────────────────────
# Maps numeric Seurat cluster IDs back to meaningful names on the
# downsampled balanced object.

new_labels <- c(
  "0"  = "Micro 1",
  "1"  = "Micro 2",
  "2"  = "Micro 3",
  "3"  = "Micro 4",
  "4"  = "Micro 5",
  "5"  = "Mono Macs",
  "6"  = "Micro 6",
  "8"  = "Micro 7",
  "9"  = "Micro 8",
  "11" = "Micro 9",
  "14" = "Micro 10"
)

Idents(balanced_obj) <- "seurat_clusters"
balanced_obj <- RenameIdents(balanced_obj, new_labels)


# ── 17. UMAP Plots (Balanced Object) ──────────────────────────────────────────

colnames(sub_obj@reductions$umap@cell.embeddings) <- c("UMAP 1", "UMAP 2")
sub_obj@reductions$umap@key <- "UMAP "

p_umap <- DimPlot(
  balanced_obj,
  reduction = "umap",
  label     = TRUE,
  repel     = TRUE
) +
  xlab("UMAP 1") + ylab("UMAP 2")

p_umap_split <- DimPlot(
  balanced_obj,
  reduction = "umap",
  split.by  = "orig.ident",
  label     = TRUE,
  repel     = TRUE
) +
  xlab("UMAP 1") + ylab("UMAP 2")

p_umap
p_umap_split


# ── 18. QC Plots (Balanced Object) ────────────────────────────────────────────

p_mt_clust  <- VlnPlot(balanced_obj, features = "percent.mt")                            + ggtitle("MT% by Cluster")
p_mt_sample <- VlnPlot(balanced_obj, features = "percent.mt",   group.by = "orig.ident") + ggtitle("MT% by Sample")
p_nf_clust  <- VlnPlot(balanced_obj, features = "nFeature_RNA")                          + ggtitle("nFeature by Cluster")
p_nf_sample <- VlnPlot(balanced_obj, features = "nFeature_RNA", group.by = "orig.ident") + ggtitle("nFeature by Sample")
p_nc_clust  <- VlnPlot(balanced_obj, features = "nCount_RNA")                            + ggtitle("nCount by Cluster")
p_nc_sample <- VlnPlot(balanced_obj, features = "nCount_RNA",   group.by = "orig.ident") + ggtitle("nCount by Sample")

# ggsave("outputs/qc_balanced_percent_mt.png",   plot = p_mt_clust | p_mt_sample, width = 18, height = 6, dpi = 300)
# ggsave("outputs/qc_balanced_nFeature_RNA.png", plot = p_nf_clust | p_nf_sample, width = 18, height = 6, dpi = 300)
# ggsave("outputs/qc_balanced_nCount_RNA.png",   plot = p_nc_clust | p_nc_sample, width = 18, height = 6, dpi = 300)


# ── 19. Gene Lists ────────────────────────────────────────────────────────────

core_genes <- c(
  # Homeostatic markers
  "P2ry12", "Cx3cr1", "Tmem119", "Tgfbr1",
  # cGAS-STING pathway
  "Cgas", "Tmem173", "Nfkbia",
  # Interferon response
  "Ifitm3", "Ifnar2", "Ifi27l2a", "Oasl2", "Irf7", "Stat1", "Isg15", "Cxcl10",
  # Pro-inflammatory cytokines
  "Tnf", "Il6", "Ccl2",
  # DAM / lipid metabolism
  "Apoe", "Lpl", "Fabp5", "Ch25h", "Gpnmb", "Lgals3bp", "Fth1", "Ftl1",
  # DAM receptors & phagocytosis
  "Trem2", "Clec7a", "Axl", "Itgax", "Cd9", "Cd63", "Tyrobp", "Csf1", "Csf1r",
  # Lysosomal/cathepsins
  "Ctsb", "Ctsd", "Ctsl", "Ctsz", "Cst7",
  # Damage signals
  "Lgals3", "Spp1",
  # Mitochondrial integrity
  "Tomm22", "Gpx4", "Fdx1", "Aifm1", "Immt"
)

homeostatic_microglia_genes <- c(
  "P2ry12", "P2ry13", "Cx3cr1", "Tmem119", "Siglech", "Fcrls",
  "Sall1", "Sall3", "Mef2c",
  "Hexb", "Cst3", "Csf1r", "Lpl",
  "C1qa", "C1qb", "C1qc",
  "Olfml3", "Gpr34", "Tgfbr1", "Tgfbr2", "Penk", "Itgb5"
)

DAM_genes <- c(
  "Apoe", "Lpl", "Ctsb", "Ctsd", "Cst7", "Tyrobp", "B2m", "Fth1", "Ftl1",
  "Spp1", "Lgals3", "Fabp5",
  "Trem2", "Clec7a", "Axl", "Itgax", "Cd9", "Cd63", "Csf1", "Csf1r",
  "Ch25h", "Gpnmb", "Lgals3bp", "Ctsl", "Ctsz"
)

IRM_genes <- c(
  "Ifit1", "Ifit2", "Ifit3", "Ifi27l2a", "Ifi44", "Ifi44l", "Isg15", "Irf7",
  "Oasl1", "Oasl2", "Mx1", "Mx2",
  "Stat1", "Stat2", "Ddx58", "Dhx58", "Rsad2", "Gbp2", "Gbp3", "Gbp5",
  "Usp18", "Cmpk2",
  "Cxcl10", "Ccl5", "Ccl12"
)

cgas_sting_genes <- c(
  "Mb21d1", "Trex1", "Dnase2a", "Zbp1",
  "Tmem173", "Tbk1", "Ikbke", "Irf3", "Irf7", "Nfkb1", "Nfkb2", "Rela",
  "Ifnb1", "Ifna4", "Ifna2",
  "Isg15", "Ifit1", "Ifit2", "Ifit3", "Ifi44", "Ifi44l",
  "Oasl1", "Oasl2", "Mx1", "Mx2", "Rsad2", "Usp18", "Cmpk2",
  "Cxcl10", "Ccl5", "Il6", "Tnf", "Il1b"
)

ARM_genes <- c(
  "H2-Ab1", "H2-Aa", "H2-Eb1", "Cd74",
  "Lyz2", "Lyz1", "Ctsb", "Ctsd", "Cst7", "Tyrobp", "B2m",
  "Il1b", "Tnf", "Cxcl10", "Ccl5", "Ccl2", "Ccl12",
  "Stat1", "Stat2", "Irf7", "Isg15", "Ifit1", "Ifit3", "Oasl1",
  "Lgals3", "Lpl", "Apoe", "Ctsl", "Ctsz",
  "Spp1", "Gpnmb", "Fabp5", "Ch25h"
)

stroke_SIM_genes <- c(
  "Spp1", "Lgals3", "Gpnmb", "Fabp5", "Ch25h",
  "Ccl2", "Ccl3", "Ccl4", "Ccl5", "Cxcl10",
  "Il1b", "Tnf", "Il6",
  "Apoe", "Lpl", "Ctsb", "Ctsd", "Cst7", "Tyrobp",
  "H2-Ab1", "H2-Aa", "H2-Eb1", "Cd74",
  "Ifit1", "Ifit3", "Isg15", "Irf7"
)

gene_lists <- list(
  Homeostatic = homeostatic_microglia_genes,
  DAM         = DAM_genes,
  IRM         = IRM_genes,
  cGAS_STING  = cgas_sting_genes,
  ARM         = ARM_genes,
  Stroke_SIM  = stroke_SIM_genes
)

# Filter to genes present in the object
core_genes <- core_genes[core_genes %in% rownames(balanced_obj)]


# ── 20. DotPlots ──────────────────────────────────────────────────────────────

dotplot_theme <- theme_bw(base_size = 14) +
  theme(
    axis.text.x  = element_text(size = 12, angle = 45, hjust = 1, color = "black"),
    axis.text.y  = element_text(size = 12, color = "black"),
    panel.grid   = element_blank(),
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 10),
    plot.margin  = margin(20, 20, 20, 20)
  )

# All clusters combined
DotPlot(balanced_obj, features = core_genes, group.by = "seurat_clusters") +
  RotatedAxis() +
  scale_color_viridis_c(option = "C") +
  scale_size(range = c(1, 8), breaks = c(0, 10, 20, 30, 40, 50)) +
  labs(color = "Average Expression", size = "Percent Expressed") +
  guides(color = guide_colorbar(order = 1), size = guide_legend(order = 2)) +
  dotplot_theme

# Multi-gene-list dotplot
DotPlot(balanced_obj, features = core_genes) +
  RotatedAxis() +
  scale_color_viridis_c(option = "C") +
  scale_size(range = c(1, 8)) +
  labs(color = "Average Expression", size = "Percent Expressed") +
  guides(color = guide_colorbar(order = 1), size = guide_legend(order = 2)) +
  dotplot_theme

# ── 21. Cell Composition Barplots ─────────────────────────────────────────────

composition_tbl <- as.data.frame(
  table(balanced_obj@active.ident, balanced_obj@meta.data$orig.ident)
)
names(composition_tbl) <- c("cluster", "sample", "count")

sex_labels <- c("FBrain" = "Female", "MBrain" = "Male")

# Stacked count barplot
p_count <- ggplot(composition_tbl, aes(fill = cluster, y = count, x = sample)) +
  geom_bar(position = "stack", stat = "identity", color = "white", linewidth = 0.3) +
  scale_x_discrete(labels = sex_labels) +
  labs(x = "Sample", y = "Cell Count", fill = "Cluster") +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 20),
    axis.text.x  = element_text(size = 14, color = "black"),
    axis.text.y  = element_text(size = 16, color = "black"),
    legend.title = element_text(size = 16),
    legend.text  = element_text(size = 14)
  )
p_count

# Proportional barplot
p_fraction <- ggplot(composition_tbl, aes(fill = cluster, y = count, x = sample)) +
  geom_bar(position = "fill", stat = "identity", color = "white", linewidth = 0.3) +
  scale_x_discrete(labels = sex_labels) +
  labs(x = "Sample", y = "Population Fraction", fill = "Cluster") +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 20),
    axis.text.x  = element_text(size = 14, color = "black"),
    axis.text.y  = element_text(size = 16, color = "black"),
    legend.title = element_text(size = 16),
    legend.text  = element_text(size = 14)
  )
p_fraction


# ── 22. Differential Expression & Heatmap ────────────────────────────────────

Idents(balanced_obj) <- "seurat_clusters"
balanced_obj <- RenameIdents(balanced_obj, new_labels)

DefaultAssay(balanced_obj) <- "RNA"
balanced_obj <- JoinLayers(balanced_obj)

markers <- FindAllMarkers(
  balanced_obj,
  only.pos        = TRUE,
  min.pct         = 0.1,
  logfc.threshold = 0.25,
  test.use        = "wilcox"
)

top5  <- markers %>% group_by(cluster) %>% top_n(n = 5,  wt = avg_log2FC)
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

balanced_small <- subset(balanced_obj, downsample = 200)

p_heatmap <- DoHeatmap(
  balanced_small,
  features         = top5$gene,
  raster           = FALSE,
  size             = 6,
  group.bar.height = 0.02
) +
  NoLegend() +
  theme(
    axis.text.y  = element_text(size = 8),
    axis.ticks.y = element_blank(),
    plot.margin  = margin(10, 60, 10, 10)
  )

p_heatmap
# ggsave("outputs/heatmap_top5_markers.png", plot = p_heatmap, width = 20, height = 15, dpi = 300)

rm(balanced_small); gc()


# ── 23. Export for Morpheus Heatmap Tool ──────────────────────────────────────
p_F <- DotPlot(subset(balanced_obj, subset = orig.ident == "FBrain"),
               features = core_genes) + RotatedAxis()
p_M <- DotPlot(subset(balanced_obj, subset = orig.ident == "MBrain"),
               features = core_genes) + RotatedAxis()

df_F <- p_F$data
df_M <- p_M$data

df_F$sex <- "Female"
df_M$sex <- "Male"

df_combined <- rbind(df_F, df_M)
# write.xlsx(df_combined, file = "outputs/DotPlot_core_genes_F_vs_M.xlsx", row.names = FALSE)

morpheus_mat <- df_combined %>%
  select(features.plot, id, sex, avg.exp) %>%
  mutate(cluster_sex = paste0(id, "_", sex)) %>%
  select(features.plot, cluster_sex, avg.exp) %>%
  pivot_wider(names_from = cluster_sex, values_from = avg.exp)

# write.csv(morpheus_mat, "outputs/Morpheus_core_genes_matrix_by_sex.csv", row.names = FALSE)


# ── 24. Pheatmap from Morpheus Matrix ────────────────────────────────────────
# Read back the exported matrix for a publication-ready clustered heatmap.

morpheus_mat <- as.matrix(read.csv("outputs/Morpheus_core_genes_matrix_by_sex.csv", row.names = 1))

## Convert all columns to numeric safely
df_num <- as.data.frame(lapply(df, function(x) as.numeric(as.character(x))))
rownames(df_num) <- rownames(df)

## Convert to matrix
mat <- as.matrix(df_num)

## ------------------------------------------------------------
## 2. Z-score NORMALIZATION
## ------------------------------------------------------------
mat_scaled <- scale(mat)

## ------------------------------------------------------------
## 3. Clean row names
## ------------------------------------------------------------
rownames(mat_scaled) <- gsub("_", " ", rownames(mat_scaled))

## ------------------------------------------------------------
## 4. Build ROW annotations (Cluster + Sex)
## ------------------------------------------------------------
samples <- rownames(mat_scaled)

row_annot <- data.frame(
  Cluster = sub(" [^ ]+$", "", samples),
  Sex     = sub(".* ", "", samples)
)
rownames(row_annot) <- samples

## ------------------------------------------------------------
## 5. Colors for annotations
## ------------------------------------------------------------
cluster_levels <- unique(row_annot$Cluster)
cluster_colors <- colorRampPalette(brewer.pal(12, "Set3"))(length(cluster_levels))
names(cluster_colors) <- cluster_levels

sex_colors <- c(Female = "#E41A1C", Male = "#377EB8")

ha_rows <- rowAnnotation(
  df = row_annot,
  col = list(
    Cluster = cluster_colors,
    Sex     = sex_colors
  )
)

ha_rows@anno_list$Cluster@show_legend <- FALSE
ha_rows@anno_list$Sex@show_legend     <- FALSE

## ------------------------------------------------------------
## 6. Heatmap color function
## ------------------------------------------------------------
col_fun <- colorRamp2(
  c(min(mat_scaled, na.rm = TRUE), 0, max(mat_scaled, na.rm = TRUE)),
  c("#2166AC", "white", "#B2182B")
)

## ------------------------------------------------------------
## 7. Draw heatmap WITH row splits between clusters
## ------------------------------------------------------------
morpheus_heatmap <- Heatmap(
  mat_scaled,
  name = "Z-score",
  col  = col_fun,
  
  ## ✅ KEY: split rows by Cluster — creates white gap between each cluster
  row_split = factor(row_annot$Cluster,
                     levels = cluster_levels),   # ✅ preserves your CSV order
  
  ## ✅ Controls the white gap size between clusters
  row_gap = unit(2.5, "mm"),                       # ✅ adjust mm to taste
  
  ## ✅ Remove the cluster label that appears in the gap by default
  row_title = NULL,                              # ✅ set to cluster_levels if you want labels
  
  left_annotation = ha_rows,
  cluster_rows    = FALSE,                       # ✅ preserve order within each cluster
  cluster_columns = FALSE,
  cluster_row_slices = FALSE,                    # ✅ preserve cluster order from CSV
  
  show_row_names   = TRUE,
  show_column_names = TRUE,
  row_names_gp     = gpar(fontsize = 14),
  column_names_gp  = gpar(fontsize = 15),
  
  heatmap_legend_param = list(
    title         = "Expression",
    legend_height = unit(4, "cm")
  )
)

print(morpheus_heatmap)


# ── 25. Cluster-Specific Marker Genes (per Cluster) ──────────────────────────
# Identifies genes that distinguish each myeloid cluster from all others.
# Results saved as individual CSVs, one per cluster.

for (cluster_id in myeloid_clusters) {
  dge <- FindMarkers(
    balanced_obj,
    ident.1         = cluster_id,
    verbose         = TRUE,
    logfc.threshold = 0.25
  )
  write.csv(
    dge,
    file = paste0("outputs/ClusterMarkers_", cluster_id, ".csv")
  )
}


# ── 26. Female vs. Male DGE per Cluster ───────────────────────────────────────
# For each myeloid cluster, compares Female vs. Male cells using Wilcoxon test.
# Positive log2FC = higher expression in Female.

for (cluster_id in myeloid_clusters) {
  dge <- FindMarkers(
    balanced_obj,
    ident.1         = "Female",
    ident.2         = "Male",
    group.by        = "sex",
    subset.ident    = cluster_id,
    verbose         = TRUE,
    logfc.threshold = 0.25
  )
  write.csv(
    dge,
    file = paste0("outputs/FemVsMale_", cluster_id, ".csv")
  )
}


# ── 27. Female vs. Male DGE Across All Microglia (Combined) ──────────────────
# Pools all microglial clusters (excludes Mono_Macs) and runs a single
# Female vs. Male comparison across the entire microglial population.

microglia_only <- c(
  "Micro 1", "Micro 2", "Micro 3", "Micro 4",  "Micro 5",
  "Micro 6", "Micro 7", "Micro 8", "Micro 9", "Micro 10"
)

micro_obj <- subset(balanced_obj, idents = microglia_only)

dge_all_micro <- FindMarkers(
  micro_obj,
  ident.1         = "Female",
  ident.2         = "Male",
  group.by        = "sex",
  verbose         = TRUE,
  logfc.threshold = 0.25
)

# write.csv(dge_all_micro, file = "outputs/FemVsMale_AllMicroglia.csv")

rm(micro_obj); gc()


# ── 28. Save Final Object ─────────────────────────────────────────────────────

# saveRDS(balanced_obj, file = "outputs/BrainsFinalClustering_balancedFandM_final.rds")

# =============================================================================
# MODULE SCORING
# =============================================================================


# ── 29. Module Scoring: Gene Signatures ───────────────────────────────────────
# Gene sets from MSigDB hallmark collections and curated literature sources.
# All gene symbols are normalized to mouse title case (e.g. "Tnf") before
# matching against rownames of the Seurat object.

tnfa_genes <- c(
  "Abca1","Ackr3","Areg","Atf3","Atp2b1","B4galt1","B4galt5","Bcl2a1d","Bcl3","Bcl6","Bhlhe40",
  "Birc2","Birc3","Bmp2","Btg1","Btg2","Btg3","Ccl20","Ccl5","Ccn1","Ccnd1","Ccnl1","Ccrl2",
  "Cd44","Cd69","Cd80","Cd83","Cdkn1a","Cebpb","Cebpd","Cflar","Clcf1","Csf1","Csf2","Cxcl1",
  "Cxcl10","Cxcl11","Cxcl2","Cxcl5","Dennd5a","Dnajb4","Dram1","Dusp1","Dusp2","Dusp4","Dusp5",
  "Edn1","Efna1","Egr1","Egr2","Egr3","Ehd1","Eif1","Ets2","F2rl1","F3","Fjx1","Fos","Fosb",
  "Fosl1","Fosl2","Fut4","G0s2","Gadd45a","Gadd45b","Gch1","Gem","Gfpt2","Gpr183","Hbegf","Hes1",
  "Icam1","Icosl","Id2","Ier2","Ier3","Ier5","Ifih1","Ifit2","Ifngr2","Il12b","Il15ra","Il18",
  "Il1a","Il1b","Il23a","Il6","Il6st","Il7r","Inhba","Irf1","Irs2","Jag1","Jun","Junb","Kdm6b",
  "Klf10","Klf2","Klf4","Klf6","Klf9","Kynu","Lamb3","Ldlr","Lif","Litaf","Maff","Map2k3","Map3k8",
  "Marcks","Mcl1","Msc","Mxd1","Myc","Nampt","Nfat5","Nfe2l2","Nfil3","Nfkb1","Nfkb2","Nfkbia",
  "Nfkbie","Ninj1","Nr4a1","Nr4a2","Nr4a3","Olr1","Panx1","Pde4b","Pdlim5","Per1","Pfkfb3",
  "Phlda1","Phlda2","Plau","Plaur","Plek","Plk2","Plpp3","Pmepa1","Pnrc1","Ppp1r15a","Ptger4",
  "Ptgs2","Ptpre","Ptx3","Rcan1","Rel","Rela","Relb","Rhob","Rigi","Ripk2","Rnf19b","Sat1","Sdc4",
  "Serpinb2","Serpinb8","Serpine1","Sgk1","Slc16a6","Slc2a3","Slc2a6","Smad3","Snn","Socs3","Sod2",
  "Sphk1","Spsb1","Sqstm1","Stat5a","Tank","Tap1","Tgif1","Tiparp","Tlr2","Tnc","Tnf","Tnfaip2",
  "Tnfaip3","Tnfaip6","Tnfaip8","Tnfrsf9","Tnfsf9","Tnip1","Tnip2","Traf1","Trib1","Trip10",
  "Tsc22d1","Tubb2a","Vegfa","Yrdc","Zbtb10","Zc3h12a","Zfp36"
)

interferon_alpha_genes <- c(
  "Adar","B2m","Batf2","Bst2","C1s","Casp1","Casp8","Ccrl2","Cd47","Cd74",
  "Cmpk2","Cnp","Csf1","Cxcl10","Cxcl11","Ddx60","Dhx58","Eif2ak2","Elf1",
  "Epsti1","Mvb12a","Tent5a","Cmtr1","Gbp2","Gbp4","Gmpr","Herc6","Hla-c",
  "Ifi27","Ifi30","Ifi35","Ifi44","Ifi44l","Ifih1","Ifit2","Ifit3","Ifitm1",
  "Ifitm2","Ifitm3","Il15","Il4r","Il7","Irf1","Irf2","Irf7","Irf9","Isg15",
  "Isg20","Lamp3","Lap3","Lgals3bp","Lpar6","Ly6e","Mov10","Mx1","Ncoa7","Nmi",
  "Nub1","Oas1","Oasl","Ogfr","Parp12","Parp14","Parp9","Plscr1","Pnpt1","Helz2",
  "Procr","Psma3","Psmb8","Psmb9","Psme1","Psme2","Ripk2","Rnf31","Rsad2","Rtp4",
  "Samd9","Samd9l","Sell","Slc25a28","Sp110","Stat2","Tap1","Tdrd7","Tmem140","Trafd1",
  "Trim14","Trim21","Trim25","Trim26","Trim5","Txnip","Uba7","Ube2l6","Usp18","Wars1"
)

inflammatory_response_genes <- c(
  "Abca1","Abi1","Acvr1b","Acvr2a","Adgre1","Adm","Adora2b","Adrm1","Ahr",
  "Aplnr","Aqp9","Atp2a2","Atp2b1","Atp2c1","Axl","Bdkrb1","Best1","Bst2",
  "Btg2","C3ar1","C5ar1","Calcrl","Ccl17","Ccl2","Ccl20","Ccl22","Ccl24",
  "Ccl5","Ccl7","Ccr7","Ccrl2","Cd14","Cd40","Cd48","Cd55","Cd69","Cd70",
  "Cd82","Cdkn1a","Chst2","Clec5a","Cmklr1","Csf1","Csf3","Csf3r","Cx3cl1",
  "Cxcl10","Cxcl11","Cxcl6","Cxcl8","Cxcl9","Cxcr6","Cybb","Dcbld2","Ebi3",
  "Edn1","Eif2ak2","Emp3","Ereg","F3","Ffar2","Fpr1","Fzd5","Gabbr1","Gch1",
  "Gna15","Gnai3","Gp1ba","Gpc3","Gpr132","Gpr183","Has2","Hbegf","Hif1a",
  "Hpn","Hrh1","Icam1","Icam4","Icoslg","Ifitm1","Ifnar1","Ifngr2","Il10",
  "Il10ra","Il12b","Il15","Il15ra","Il18","Il18r1","Il18rap","Il1a","Il1b",
  "Il1r1","Il2rb","Il4r","Il6","Il7r","Inhba","Irak2","Irf1","Irf7","Itga5",
  "Itgb3","Itgb8","Kcna3","Kcnj2","Kcnmb2","Kif1b","Klf6","Lamp3","Lck","Lcp2",
  "Ldlr","Lif","Lpar1","Lta","Ly6e","Lyn","Marco","Mefv","Mep1a","Met","Mmp14",
  "Msr1","Mxd1","Myc","Nampt","Ndp","Nfkb1","Nfkbia","Nlrp3","Nmi","Nmur1","Nod2",
  "Npffr2","Olr1","Oprk1","Osm","Osmr","P2rx4","P2rx7","P2ry2","Pcdh7","Pde4b","Pdpn",
  "Pik3r5","Plaur","Prok2","Psen1","Ptafr","Ptger2","Ptger4","Ptgir","Ptpre","Pvr","Raf1",
  "Rasgrp1","Rela","Rgs1","Rgs16","Rhog","Ripk2","Rnf144b","Ros1","Rtp4","Scarf1","Scn1b",
  "Sele","Selenos","Sell","Sema4d","Serpine1","Sgms2","Slamf1","Slc11a2","Slc1a2","Slc28a2",
  "Slc31a1","Slc31a2","Slc4a4","Slc7a1","Slc7a2","Sphk1","Sri","Stab1","Tacr1","Tacr3","Tapbp",
  "Timp1","Tlr1","Tlr2","Tlr3","Tnfaip6","Tnfrsf1b","Tnfrsf9","Tnfsf10","Tnfsf15","Tnfsf9","Tpbg","Vip"
)

DNA_repair_genes <- c(
  "Aaas","Ada","Adcy6","Adrm1","Ago4","Ak1","Ak3","Alyref","Aprt","Arl6ip1","Bcam",
  "Bcap31","Bola2","Brf2","Cant1","Ccno","Cda","Cetn2","Clp1","Cmpk2","Cox17","Cstf3",
  "Dad1","Dctn4","Ddb1","Ddb2","Dgcr8","Dguok","Dut","Edf1","Eif1b","Ell","Eloa",
  "Ercc1","Ercc2","Ercc3","Ercc4","Ercc5","Ercc8","Fen1","Gmpr2","Gpx4","Gsdme","Gtf2a2",
  "Gtf2b","Gtf2f1","Gtf2h1","Gtf2h3","Gtf2h5","Gtf3c5","Guk1","Hcls1","Hprt1","Impdh2",
  "Itpa","Lig1","Mpc2","Mpg","Mrpl40","Ncbp2","Nelfb","Nelfcd","Nelfe","Nfx1","Nme1",
  "Nme3","Nme4","Npr2","Nt5c","Nt5c3a","Nudt21","Nudt9","Pcna","Pde4b","Pde6g","Pnp",
  "Pola1","Pola2","Polb","Pold1","Pold3","Pold4","Pole4","Polh","Poll","Polr1c","Polr1d",
  "Polr1h","Polr2a","Polr2c","Polr2d","Polr2e","Polr2f","Polr2g","Polr2h","Polr2i","Polr2j",
  "Polr2k","Polr3c","Polr3gl","Pom121","Prim1","Rad51","Rad52","Rae1","Rala","Rbx1","Rev3l",
  "Rfc2","Rfc3","Rfc4","Rfc5","Rnmt","Rpa2","Rpa3","Rrm2b","Sac3d1","Sdcbp","Sec61a1",
  "Sf3a3","Smad5","Snapc4","Snapc5","Srsf6","Ssrp1","Stx3","Supt4h1","Supt5h","Surf1",
  "Taf10","Taf12","Taf13","Taf1c","Taf6","Taf9","Tarbp2","Tk2","Tmed2","Tp53","Tsg101",
  "Tyms","Umps","Upf3b","Usp11","Vps28","Vps37b","Vps37d","Xpc","Znf707","Zwint"
)

serotonin_genes <- c(
  "Gpm6b","Itgb3","Nos1","Slc18a1","Slc18a2","Slc18a3","Slc18b1","Slc22a1",
  "Slc22a2","Slc22a3","Slc29a4","Slc6a4","Snca"
)

oxidative_stress_genes <- c(
  "Abcc1","Atox1","Cat","Cdkn2d","Egln2","Ercc2","Fes","Ftl","G6pd",
  "Gclc","Gclm","Glrx","Glrx2","Gpx3","Gpx4","Gsr","Hhex","Hmox2",
  "Ipcef1","Junb","Lamtor5","Lsp1","Mbp","Mgst1","Mpo","Msra","Ndufa6",
  "Ndufb4","Ndufs2","Nqo1","Oxsr1","Pdlim1","Pfkp","Prdx1","Prdx2",
  "Prdx4","Prdx6","Prnp","Ptpa","Sbno2","Scaf4","Selenos","Sod1",
  "Sod2","Srxn1","Stk25","Txn","Txnrd1","Txnrd2"
)

serotonin_kegg_genes <- c(
  "Tph1","Ddc","Slc18a1","Slc18a2","Cacna1c","Cacna1d","Cacna1f",
  "Cacna1s","Htr2a","Htr2b","Htr2c","Gnaq","Plcb1","Plcb2","Plcb3",
  "Plcb4","Itpr1","Itpr2","Itpr3","Prkca","Prkcb","Prkcg","Mapk1",
  "Mapk3","Pla2g4b","Pla2g4e","Pla2g4f","Pla2g4d","Pla2g4a","Pla2g4c",
  "Cyp2c8","Cyp2c9","Cyp2c18","Cyp2c19","Cyp2d6","Cyp2j2","Cyp4x1",
  "Alox5","Alox12","Alox12b","Alox15","Alox15b","Ptgs1","Ptgs2",
  "Htr3c","Htr3d","Htr3e","Htr3a","Htr3b","Htr4","Htr6","Htr7",
  "Gnas","Adcy5","Prkaca","Prkacb","Prkacg","Kcnn2","Kcnd2",
  "Gabrb1","Gabrb2","Gabrb3","Rapgef3","App","Htr1a","Htr1b",
  "Htr1d","Htr1e","Htr1f","Htr5a","Gnai1","Gnai2","Gnai3",
  "Gnao1","Gnb1","Gnb2","Gnb3","Gnb4","Gnb5","Gng2","Gng3","Gng4",
  "Gng5","Gng7","Gng8","Gng10","Gng11","Gng12","Gng13","Gngt1","Gngt2",
  "Casp3","Dusp1","Hras","Kras","Nras","Araf","Braf","Raf1","Map2k1",
  "Cacna1a","Cacna1b","Kcnj3","Kcnj6","Kcnj9","Kcnj5","Trpc1","Slc6a4",
  "Maoa","Maob"
)

cgas_sting_module_genes <- c(
  "Aars2","Akt1","Aurkb","Banf1","Btk","Cgas","Ddx41","Irf3","Irgm","Lyplal1",
  "Map3k7","Marchf5","Parp1","Pcbp2","Ppp6c","Prkdc","Rnf39","Slc19a1",
  "Smpdl3a","Spsb3","Sting1","Tab1","Tmem173","Tbk1","Trex1","Zdhhc18","Zdhhc9"
)

ifn_a_module_genes <- c(
  "Ap3b1","Ap3d1","Chuk","D1pas1","Ddx3x","Dhx36","Dhx9","Flt3","Gbp4","Havcr2",
  "Hmgb1","Hspd1","Ifih1","Irf3","Irf7","Mavs","Mmp12","Nlrc3","Nmb","Nmbr",
  "Nmi","Ptpn22","Ptprs","Rigi","Setd2","Stat1","Tbk1","Tlr3","Tlr4","Tlr7",
  "Tlr8","Tlr9","Traf3ip3","Trex1","Trim65","Zc3hav1"
)

microglia_homeostatic_genes <- c(
  "Csmd3","Med12l","Ccr5","Cst3","Gpr155","Cx3cr1","Gcnt1","Gpr34","Gtf2h2",
  "Tmem119","Arhgap5","Mfap3","Golm1","P2ry13","Rab39","Pmepa1","Sall1",
  "Selplg","P2ry12","Sparc","Tlr3","Lrrc3","Plxdc2","Cd164","Lrba","Hexb",
  "Olfml3","Cd81","Crybb1","Tmem173","Srgap2","Txnip",
  "Cmtm6","Lpcat2","Rhob","Maf","Rgs2","Slco2b1","Glul","Siglech","Lgmn",
  "Csf1r","Marcks","Serinc3"
)


# ── 30. Add All Module Scores ─────────────────────────────────────────────────
# Gene symbols are already in mouse title case above.
# Signatures are filtered to genes present in the object before scoring.

signatures <- list(
  TNFA_SIGNALING_VIA_NFKB = tnfa_genes,
  INTERFERON_ALPHA         = interferon_alpha_genes,
  INFLAMMATORY_RESPONSE    = inflammatory_response_genes,
  DNA_REPAIR               = DNA_repair_genes,
  SEROTONIN                = serotonin_genes,
  OXIDATIVE_STRESS         = oxidative_stress_genes,
  SEROTONIN_KEGG           = serotonin_kegg_genes,
  CGAS_STING               = cgas_sting_module_genes,
  IFN_A                    = ifn_a_module_genes,
  MICROGLIA_HOMEOSTATIC    = microglia_homeostatic_genes
)

signatures_clean <- lapply(signatures, function(g) {
  intersect(g, rownames(balanced_obj))
})

# Drop any signature with zero matched genes
signatures_clean <- signatures_clean[sapply(signatures_clean, length) > 0]

cat("\nGenes matched per signature:\n")
print(sapply(signatures_clean, length))

for (sig in names(signatures_clean)) {
  balanced_obj <- AddModuleScore(
    balanced_obj,
    features = list(signatures_clean[[sig]]),
    name     = paste0(sig, "_SCORE")
  )
}

# Seurat appends "1" to the name — strip it for clean column names
colnames(balanced_obj@meta.data) <- sub("_SCORE1$", "_SCORE", colnames(balanced_obj@meta.data))

# Confirm column names
module_cols <- grep("_SCORE$", colnames(balanced_obj@meta.data), value = TRUE)
cat("\nModule score columns:\n")
print(module_cols)

# Build UMAP long table for module score plotting
umap_df <- as.data.frame(Embeddings(balanced_obj, reduction = "umap"))
colnames(umap_df)[1:2] <- c("umap_1", "umap_2")
umap_df$cell <- rownames(umap_df)

umap_long <- umap_df %>%
  left_join(
    balanced_obj@meta.data %>%
      select(all_of(module_cols)) %>%
      mutate(cell = rownames(balanced_obj@meta.data)),
    by = "cell"
  ) %>%
  pivot_longer(
    cols      = all_of(module_cols),
    names_to  = "module",
    values_to = "score"
  )

# ── 31. Filter & Clean Module Labels ──────────────────────────────────────────
# Select the biologically relevant modules for the publication figure and
# apply clean display names.

modules_to_plot <- c(
  "INFLAMMATORY_RESPONSE_SCORE",
  "INTERFERON_ALPHA_SCORE",
  "TNFA_SIGNALING_VIA_NFKB_SCORE",
  "OXIDATIVE_STRESS_SCORE",
  "DNA_REPAIR_SCORE",
  "MICROGLIA_HOMEOSTATIC_SCORE"
)

module_display_names <- c(
  "MICROGLIA_HOMEOSTATIC_SCORE"    = "Homeostatic microglia",
  "INFLAMMATORY_RESPONSE_SCORE"    = "Inflammatory response",
  "INTERFERON_ALPHA_SCORE"         = "Interferon-\u03b1",
  "TNFA_SIGNALING_VIA_NFKB_SCORE"  = "TNF\u03b1 signaling via NF-\u03baB",
  "OXIDATIVE_STRESS_SCORE"         = "Oxidative stress",
  "DNA_REPAIR_SCORE"               = "DNA repair"
)

module_level_order <- c(
  "Homeostatic microglia",
  "Inflammatory response",
  "Interferon-\u03b1",
  "TNF\u03b1 signaling via NF-\u03baB",
  "Oxidative stress",
  "DNA repair"
)

umap_long_subset <- umap_long %>%
  filter(module %in% modules_to_plot) %>%
  mutate(
    module = recode(module, !!!module_display_names),
    module = factor(module, levels = module_level_order)
  )


# ── 32. Clip, Rescale & Plot Selected Modules ─────────────────────────────────
# Clips the top 1% of scores per module to reduce the influence of outlier
# cells, then rescales to 0–1 for a consistent color axis across panels.

umap_long_subset <- umap_long_subset %>%
  group_by(module) %>%
  mutate(
    score_clipped = pmin(score, quantile(score, 0.99)),
    score_scaled  = rescale(score_clipped)
  ) %>%
  ungroup()

p_modules <- ggplot(umap_long_subset, aes(umap_1, umap_2, color = score_scaled)) +
  geom_point(size = 0.15, alpha = 0.8) +
  scale_color_viridis_c(
    option = "C",
    name   = "Module Score",
    guide  = guide_colorbar(
      barheight = unit(3,   "cm"),
      barwidth  = unit(0.3, "cm")
    )
  ) +
  facet_wrap(~ module, ncol = 2, strip.position = "top") +
  theme_minimal(base_size = 8) +
  theme(
    strip.text       = element_text(size = 15, margin = margin(b = 3)),
    strip.background = element_blank(),
    axis.title       = element_blank(),
    axis.text        = element_blank(),
    axis.ticks       = element_blank(),
    panel.grid       = element_blank(),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.3),
    legend.position  = "right",
    legend.title     = element_text(size = 11),
    legend.text      = element_text(size = 9)
  )

p_modules


# ── 33. Save Module Score UMAP to PDF ─────────────────────────────────────────

# cairo_pdf("outputs/module_scoring_selected.pdf", width = 10, height = 8.33)
# print(p_modules)
# dev.off()
# 
# ggsave(
#   "outputs/module_scoring_selected.png",
#   plot   = p_modules,
#   width  = 3000,
#   height = 2500,
#   units  = "px",
#   dpi    = 300)

# =============================================================================
# SECTION 34 — Hallmark GSEA: Female vs. Male per Myeloid Cluster
# =============================================================================
#
#   34a. Load packages & resolve conflicts
#   34b. Prepare objects & load Hallmark gene sets
#   34c. DE ranking helper function
#   34d. Run GSEA per cluster
#   34e. Summarize & export results
#   34f. NES heatmap tile plot
#   34g. Enrichment plot for top pathway in Micro 1
#   34h. DotPlot (final — filtered, title-case, ordered clusters)
# =============================================================================

# ── 34a. Load Packages & Resolve Conflicts ────────────────────────────────────

library(Seurat)
library(dplyr)
library(data.table)
library(msigdbr)
library(fgsea)
library(ggplot2)
library(purrr)
library(tidyverse)
library(conflicted)

# Install presto for fast Wilcoxon DE (only needs to be run once)
# install.packages("devtools")
# devtools::install_github("immunogenomics/presto")

conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::rename)
conflicts_prefer(dplyr::mutate)
conflicts_prefer(dplyr::slice)
conflicts_prefer(base::intersect)
conflicts_prefer(base::union)
conflicts_prefer(base::setdiff)


# ── 34b. Prepare Objects & Load Hallmark Gene Sets ────────────────────────────

balanced_obj <- JoinLayers(balanced_obj)

DefaultAssay(balanced_obj) <- "RNA"

# Pull mouse Hallmark gene sets from MSigDB
hallmark_mouse <- msigdbr(species = "Mus musculus", category = "H")

hallmark_list <- hallmark_mouse %>%
  split(.$gs_name) %>%
  map(~ unique(.x$gene_symbol))

# Ensure orig.ident is a two-level factor with FBrain as the reference
balanced_obj$orig.ident <- factor(balanced_obj$orig.ident, levels = c("FBrain", "MBrain"))

# Confirm named cluster identities are set correctly
cat("\nCluster identity counts:\n")
print(table(Idents(balanced_obj)))


# ── 34c. DE Ranking Helper Function ──────────────────────────────────────────
# Subsets to one cluster, switches identity to orig.ident, and runs
# FindMarkers (FBrain vs MBrain). Returns the full DE table with a
# rank_metric column (avg_log2FC) used to order genes for fgsea.

get_ranked_stats_one_cluster <- function(seu,
                                         cluster_label,
                                         group_var       = "orig.ident",
                                         group1          = "FBrain",
                                         group2          = "MBrain",
                                         test.use        = "wilcox",
                                         min.pct         = 0.1,
                                         logfc.threshold = 0) {
  seu_cl         <- subset(seu, idents = cluster_label)
  Idents(seu_cl) <- group_var
  
  markers <- FindMarkers(
    object          = seu_cl,
    ident.1         = group1,
    ident.2         = group2,
    test.use        = test.use,
    min.pct         = min.pct,
    logfc.threshold = logfc.threshold
  )
  
  markers$gene <- rownames(markers)
  
  # Support both Seurat v4 (avg_logFC) and v5 (avg_log2FC)
  stat_col <- if ("avg_log2FC" %in% colnames(markers)) "avg_log2FC" else "avg_logFC"
  
  markers %>% mutate(rank_metric = .data[[stat_col]])
}


# ── 34d. Run GSEA per Cluster ─────────────────────────────────────────────────
# Clusters defined explicitly to match levels(Idents(balanced_obj)).
# nperm removed so fgseaMultilevel is used automatically (more accurate).
# Micro 10 has only 26 cells and will likely be skipped by the n < 50 guard.

clusters <- c(
  "Micro 1", "Micro 2", "Micro 3", "Micro 4",  "Micro 5",
  "Mono Macs", "Micro 6", "Micro 7", "Micro 8", "Micro 9", "Micro 10"
)

gsea_results_list <- list()
de_tables_list    <- list()

set.seed(123)

for (cl in clusters) {
  message("Running DE + GSEA for cluster: ", cl)
  
  de_cl        <- get_ranked_stats_one_cluster(balanced_obj, cluster_label = cl)
  
  ranks        <- de_cl$rank_metric
  names(ranks) <- de_cl$gene
  
  # Keep only genes present in any Hallmark set
  common_genes <- base::intersect(names(ranks), unique(unlist(hallmark_list)))
  ranks        <- ranks[common_genes]
  
  if (length(ranks) < 50) {
    warning("Cluster ", cl, " has too few overlapping genes with Hallmark sets; skipping.")
    next
  }
  
  fgsea_res         <- fgsea(
    pathways = hallmark_list,
    stats    = ranks,
    minSize  = 10,
    maxSize  = 500
    # nperm removed — fgseaMultilevel runs automatically and is more accurate
  )
  fgsea_res$cluster <- cl
  
  gsea_results_list[[cl]] <- fgsea_res
  de_tables_list[[cl]]    <- de_cl
}


# ── 34e. Summarize & Export Results ──────────────────────────────────────────

gsea_results <- bind_rows(gsea_results_list) %>%
  arrange(cluster, padj)

sig_gsea <- gsea_results %>%
  dplyr::filter(padj < 0.05)

cat("\nSignificant pathway hits (padj < 0.05):\n")
print(head(sig_gsea))

# Flatten leadingEdge list-column to semicolon-separated strings for CSV export
gsea_results_export <- gsea_results %>%
  mutate(
    leadingEdge = sapply(leadingEdge, paste, collapse = ";"),
    pathway     = as.character(pathway),
    cluster     = as.character(cluster)
  )

# write.csv(
#   gsea_results_export,
#   "outputs/Hallmark_GSEA_F_vs_M_by_microglia_cluster.csv",
#   row.names = FALSE
# )

# ── 34f. NES Heatmap Tile Plot ────────────────────────────────────────────────
# Wide matrix of NES values: rows = significant pathways, columns = clusters.

nes_mat <- sig_gsea %>%
  dplyr::select(pathway, cluster, NES) %>%
  pivot_wider(names_from = cluster, values_from = NES)

pathway_names     <- nes_mat$pathway
nes_mat           <- as.matrix(nes_mat[, -1])
rownames(nes_mat) <- pathway_names

nes_df <- as.data.frame(nes_mat) %>%
  tibble::rownames_to_column("pathway") %>%
  pivot_longer(-pathway, names_to = "cluster", values_to = "NES")

ggplot(nes_df, aes(x = cluster, y = pathway, fill = NES)) +
  geom_tile() +
  scale_fill_gradient2(low = "#2266ac", mid = "white", high = "#b82c34") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title  = element_blank()
  ) +
  ggtitle("Hallmark GSEA NES (Female vs Male) by microglia cluster")


# ── 34g. Enrichment Plot: Top Pathway in Micro 1 ─────────────────────────────
# Identifies the most significant pathway in Micro 1 and plots the classic
# running-sum enrichment curve.

cl_of_interest <- "Micro 1"

res_cl      <- gsea_results %>%
  dplyr::filter(cluster == cl_of_interest) %>%
  arrange(padj)

top_pathway <- res_cl$pathway[1]

de_cl        <- de_tables_list[[cl_of_interest]]
ranks        <- de_cl$rank_metric
names(ranks) <- de_cl$gene
ranks        <- sort(ranks, decreasing = TRUE)

plotEnrichment(
  hallmark_list[[top_pathway]],
  ranks
) + ggtitle(paste(cl_of_interest, "-", top_pathway))


# ── 34h. DotPlot: Final (Filtered, Title-Case, Ordered Clusters) ─────────────
# Pathway names are cleaned (HALLMARK_ prefix and underscores removed, then
# title-cased). Six off-topic pathways are excluded. Clusters are ordered
# left-to-right by microglial subtype with Mono Macs last.

cluster_order <- c(
  "Micro 1", "Micro 2", "Micro 3", "Micro 4",  "Micro 5",
  "Micro 6", "Micro 7", "Micro 8", "Micro 9", "Micro 10", "Mono Macs"
)

dot_df <- sig_gsea %>%
  mutate(
    pathway      = as.character(pathway),
    cluster      = as.character(cluster),
    neglog10padj = -log10(padj),
    pathway      = gsub("^HALLMARK_", "", pathway),
    pathway      = gsub("_", " ", pathway),
    pathway      = tools::toTitleCase(tolower(pathway))
  ) %>%
  dplyr::filter(!pathway %in% c(
    "Xenobiotic Metabolism",
    "Notch Signaling",
    "Kras Signaling Up",
    "Hedgehog Signaling",
    "Androgen Response",
    "Estrogen Response Late"
  )) %>%
  mutate(cluster = factor(cluster, levels = cluster_order))

p_gsea_dot <- ggplot(dot_df, aes(x = cluster, y = pathway)) +
  geom_point(aes(size = neglog10padj, color = NES)) +
  scale_color_gradient2(low = "#2266ac", mid = "white", high = "#b82c34") +
  scale_size_continuous(range = c(1, 8)) +
  theme_bw() +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 14, color = "black"),
    axis.text.y  = element_text(size = 14, color = "black"),
    axis.title   = element_blank(),
    plot.title   = element_text(hjust = 0, size = 16, color = "black"),
    legend.text  = element_text(size = 14, color = "black"),
    legend.title = element_text(size = 14, color = "black")
  ) +
  ggtitle("Hallmark GSEA DotPlot (Female vs Male)")

p_gsea_dot

# ggsave(
#   "outputs/GSEA_dotplot_final.png",
#   plot   = p_gsea_dot,
#   width  = 3500,
#   height = 2500,
#   units  = "px",
#   dpi    = 300
# )

# cairo_pdf("outputs/GSEA_dotplot_final.pdf", width = 11.67, height = 8.33)
# print(p_gsea_dot)
# dev.off()