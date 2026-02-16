#!/usr/bin/env Rscript

# ============================================================
# Pan-cancer perivascular-like macrophages + CAF ligands (NicheNet)
# Dataset: GSE210347
# Inputs:  GSE210347_meta.txt, GSE210347_counts.Rds
# Outputs: ~/projects/nichenet/results/figures and ~/projects/nichenet/results/tables
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(nichenetr)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(pheatmap)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(patchwork)
})

set.seed(123)

# ----------------------------
# Paths (edit once)
# ----------------------------
base_dir   <- path.expand("~/projects/nichenet")
fig_dir    <- file.path(base_dir, "results", "figures")
table_dir  <- file.path(base_dir, "results", "tables")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

meta_path   <- "GSE210347_meta.txt"
counts_path <- "GSE210347_counts.Rds"

nichenet_path <- path.expand("~/nichenet_data")
lr_path       <- file.path(nichenet_path, "lr_network_human_21122021.rds")
ltm_path      <- file.path(nichenet_path, "ligand_target_matrix_nsga2r_final.rds")
wn_path       <- file.path(nichenet_path, "weighted_networks_nsga2r_final.rds")

# ----------------------------
# Load data
# ----------------------------
message("Loading metadata and counts ...")
meta <- read.table(meta_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
counts <- readRDS(counts_path)

stopifnot("cellname" %in% colnames(meta))
stopifnot("celltype" %in% colnames(meta))
stopifnot("tissue" %in% colnames(meta))
stopifnot(all(meta$cellname %in% colnames(counts)))

rownames(meta) <- meta$cellname
counts <- counts[, meta$cellname, drop = FALSE]

message("Cells: ", nrow(meta))
message("Genes: ", nrow(counts))

# ============================================================
# Figure 1: Global UMAP (subsampled)
# ============================================================
message("Building global UMAP (50k subsample) ...")

meta_sub <- meta
if (nrow(meta) > 50000) {
  idx <- sample(seq_len(nrow(meta)), 50000)
  meta_sub <- meta[idx, , drop = FALSE]
}

counts_sub <- counts[, rownames(meta_sub), drop = FALSE]

global_obj <- CreateSeuratObject(
  counts = counts_sub,
  meta.data = meta_sub,
  min.cells = 3,
  min.features = 200
)

global_obj <- NormalizeData(global_obj, verbose = FALSE)
global_obj <- FindVariableFeatures(global_obj, verbose = FALSE)
global_obj <- ScaleData(global_obj, features = VariableFeatures(global_obj), verbose = FALSE)
global_obj <- RunPCA(global_obj, npcs = 30, verbose = FALSE)
global_obj <- RunUMAP(global_obj, dims = 1:30, verbose = FALSE)

p_global <- DimPlot(global_obj, group.by = "celltype", label = TRUE, repel = TRUE) +
  ggtitle("Global UMAP (50k subsample)")

ggsave(file.path(fig_dir, "01_global_umap.png"), p_global, width = 12, height = 8, dpi = 300)

rm(global_obj, counts_sub, meta_sub); gc()

# ============================================================
# Subset: Myeloid -> CD68+ macrophage-enriched
# ============================================================
message("Subsetting myeloid ...")

myeloid_cells <- rownames(meta)[meta$celltype == "myeloid"]
stopifnot(length(myeloid_cells) > 0)

counts_myeloid <- counts[, myeloid_cells, drop = FALSE]
meta_myeloid <- meta[myeloid_cells, , drop = FALSE]

myeloid <- CreateSeuratObject(counts = counts_myeloid, meta.data = meta_myeloid)
myeloid <- NormalizeData(myeloid, verbose = FALSE)
myeloid <- FindVariableFeatures(myeloid, verbose = FALSE)
myeloid <- ScaleData(myeloid, features = VariableFeatures(myeloid), verbose = FALSE)
myeloid <- RunPCA(myeloid, npcs = 30, verbose = FALSE)
myeloid <- RunUMAP(myeloid, dims = 1:30, verbose = FALSE)

cd68_present <- "CD68" %in% rownames(myeloid)
if (!cd68_present) stop("CD68 not found in expression matrix rownames(). Check gene symbols.")

cd68_counts <- GetAssayData(myeloid, slot = "counts")["CD68", ]
myeloid$CD68_pos <- cd68_counts > 0

p_myeloid_cd68 <- DimPlot(myeloid, group.by = "CD68_pos") +
  ggtitle("Myeloid UMAP: CD68+ highlight")

ggsave(file.path(fig_dir, "02_myeloid_cd68.png"), p_myeloid_cd68, width = 12, height = 7, dpi = 300)

mac_cells <- colnames(myeloid)[myeloid$CD68_pos]
counts_mac <- GetAssayData(myeloid, slot = "counts")[, mac_cells, drop = FALSE]
meta_mac <- meta[mac_cells, , drop = FALSE]

mye <- CreateSeuratObject(counts = counts_mac, meta.data = meta_mac)
mye <- NormalizeData(mye, verbose = FALSE)
mye <- FindVariableFeatures(mye, verbose = FALSE)
mye <- ScaleData(mye, features = VariableFeatures(mye), verbose = FALSE)
mye <- RunPCA(mye, npcs = 30, verbose = FALSE)
mye <- RunUMAP(mye, dims = 1:30, verbose = FALSE)

p_mac_umap_celltype <- DimPlot(mye, group.by = "celltype") + ggtitle("CD68+ subset")
p_mac_umap_tissue <- DimPlot(mye, group.by = "tissue") + ggtitle("CD68+ subset by tissue")
ggsave(file.path(fig_dir, "03_macrophage_umap.png"), p_mac_umap_celltype | p_mac_umap_tissue,
       width = 14, height = 6, dpi = 300)

rm(myeloid, counts_myeloid, meta_myeloid, counts_mac, meta_mac); gc()

# ============================================================
# Figure 3-4: Marker panels + perivascular-like score
# ============================================================
message("Scoring perivascular-like program ...")

cat_a <- c("LYVE1", "FOLR2", "MRC1", "CD163", "STAB1")
cat_b <- c("VEGFA", "ANGPT2", "PDGFB", "SPP1", "TREM2", "EDNRB")
cat_c <- c("CCL18", "MMP9", "CHI3L1")

all_peri_markers <- unique(c(cat_a, cat_b, cat_c))
all_peri_markers <- all_peri_markers[all_peri_markers %in% rownames(mye)]

if (length(all_peri_markers) < 5) {
  stop("Too few perivascular markers detected in data. Detected: ",
       paste(all_peri_markers, collapse = ", "))
}

p_markers <- FeaturePlot(
  mye,
  features = all_peri_markers,
  ncol = 4,
  min.cutoff = "q10",
  max.cutoff = "q90"
) + plot_annotation(title = "Perivascular-like marker panel (feature plots)")

ggsave(file.path(fig_dir, "04_perivascular_markers.png"), p_markers, width = 16, height = 14, dpi = 300)

mye <- AddModuleScore(mye, features = list(all_peri_markers), name = "Perivascular")
thr <- as.numeric(stats::quantile(mye$Perivascular1, 0.75, na.rm = TRUE))
mye$Perivascular_high <- mye$Perivascular1 > thr

p_score <- FeaturePlot(mye, features = "Perivascular1", reduction = "umap") +
  ggtitle("Perivascular module score") +
  scale_color_viridis_c()

p_gate <- DimPlot(mye, group.by = "Perivascular_high") +
  ggtitle("PV-high gate (top quartile)")

ggsave(file.path(fig_dir, "05_perivascular_score_gate.png"), p_score | p_gate,
       width = 14, height = 6, dpi = 300)

# PV score by tissue
p_violin <- VlnPlot(mye, features = "Perivascular1", group.by = "tissue", pt.size = 0) +
  ggtitle("Perivascular module score by tissue") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(fig_dir, "06_perivascular_score_by_tissue.png"), p_violin,
       width = 14, height = 6, dpi = 300)

# ============================================================
# Fibroblasts: subset + PDGFRA/PDGFRB
# ============================================================
message("Processing fibroblasts ...")

fib_cells <- rownames(meta)[meta$celltype == "fibroblast"]
stopifnot(length(fib_cells) > 0)

counts_fib <- counts[, fib_cells, drop = FALSE]
meta_fib <- meta[fib_cells, , drop = FALSE]

fibro <- CreateSeuratObject(counts = counts_fib, meta.data = meta_fib)
fibro <- NormalizeData(fibro, verbose = FALSE)
fibro <- FindVariableFeatures(fibro, verbose = FALSE)
fibro <- ScaleData(fibro, features = VariableFeatures(fibro), verbose = FALSE)
fibro <- RunPCA(fibro, npcs = 30, verbose = FALSE)
fibro <- RunUMAP(fibro, dims = 1:30, verbose = FALSE)

pdgfra_present <- "PDGFRA" %in% rownames(fibro)
pdgfrb_present <- "PDGFRB" %in% rownames(fibro)

if (!pdgfra_present || !pdgfrb_present) {
  stop("PDGFRA and/or PDGFRB not found in fibroblast gene names. Check gene symbols.")
}

raw_fib <- GetAssayData(fibro, slot = "counts")
fibro$CAF_subtype <- "Other_fibroblast"
fibro$CAF_subtype[raw_fib["PDGFRA", ] > 0 & raw_fib["PDGFRB", ] > 0] <- "PDGFRA_PDGFRB_fibro"

p_fib_markers <- FeaturePlot(fibro, features = c("PDGFRA", "PDGFRB"), ncol = 2) +
  plot_annotation(title = "Fibroblast marker expression (PDGFRA, PDGFRB)")

ggsave(file.path(fig_dir, "07_fibroblast_markers.png"), p_fib_markers,
       width = 12, height = 6, dpi = 300)

rm(raw_fib, counts_fib, meta_fib); gc()

# ============================================================
# Load NicheNet resources
# ============================================================
message("Loading NicheNet resources ...")

lr_network <- readRDS(lr_path) %>% dplyr::distinct(from, to)
ligand_target_matrix <- readRDS(ltm_path)
weighted_networks <- readRDS(wn_path)

all_receptors <- unique(lr_network$to)

# ============================================================
# Eligibility per tissue (minimum cell counts)
# ============================================================
MIN_RECEIVER <- 250
MIN_SENDER <- 250

receiver_counts <- mye@meta.data %>%
  dplyr::group_by(tissue) %>%
  dplyr::summarise(n_peri_high = sum(Perivascular_high), .groups = "drop")

sender_counts <- fibro@meta.data %>%
  dplyr::group_by(tissue) %>%
  dplyr::summarise(n_caf = sum(CAF_subtype == "PDGFRA_PDGFRB_fibro"), .groups = "drop")

eligibility <- dplyr::full_join(receiver_counts, sender_counts, by = "tissue") %>%
  dplyr::mutate(
    n_peri_high = ifelse(is.na(n_peri_high), 0, n_peri_high),
    n_caf = ifelse(is.na(n_caf), 0, n_caf),
    qualifies = n_peri_high >= MIN_RECEIVER & n_caf >= MIN_SENDER
  ) %>%
  dplyr::arrange(dplyr::desc(n_peri_high))

write.csv(eligibility, file.path(table_dir, "percancer_eligibility.csv"), row.names = FALSE)

qualifying_tissues <- eligibility$tissue[eligibility$qualifies]
if (length(qualifying_tissues) == 0) {
  MIN_RECEIVER <- 100
  MIN_SENDER <- 100
  qualifying_tissues <- eligibility$tissue[eligibility$n_peri_high >= MIN_RECEIVER & eligibility$n_caf >= MIN_SENDER]
}

message("Qualifying tissues (NicheNet): ", paste(qualifying_tissues, collapse = ", "))

# ============================================================
# NicheNet per tissue
# ============================================================
per_cancer_results <- list()

for (tiss in qualifying_tissues) {

  receiver_cells <- colnames(mye)[mye$Perivascular_high & mye$tissue == tiss]
  sender_cells <- colnames(fibro)[fibro$CAF_subtype == "PDGFRA_PDGFRB_fibro" & fibro$tissue == tiss]

  if (length(receiver_cells) < 50 || length(sender_cells) < 50) next

  expr_receiver <- rownames(mye)[Matrix::rowSums(
    GetAssayData(mye, slot = "counts")[, receiver_cells, drop = FALSE] > 0
  ) / length(receiver_cells) >= 0.05]

  expr_sender <- rownames(fibro)[Matrix::rowSums(
    GetAssayData(fibro, slot = "counts")[, sender_cells, drop = FALSE] > 0
  ) / length(sender_cells) >= 0.05]

  expressed_receptors_tiss <- intersect(all_receptors, expr_receiver)

  potential_ligands_tiss <- lr_network %>%
    dplyr::filter(to %in% expressed_receptors_tiss) %>%
    dplyr::pull(from) %>%
    unique()

  potential_ligands_focused <- intersect(potential_ligands_tiss, expr_sender)
  if (length(potential_ligands_focused) < 5) next

  mye_tiss <- subset(mye, tissue == tiss)
  Idents(mye_tiss) <- mye_tiss$Perivascular_high

  de <- FindMarkers(
    mye_tiss,
    ident.1 = TRUE,
    ident.2 = FALSE,
    min.pct = 0.1,
    logfc.threshold = 0.25,
    verbose = FALSE
  )

  geneset <- de %>%
    as.data.frame() %>%
    dplyr::filter(p_val_adj <= 0.05, avg_log2FC >= 0.25) %>%
    rownames() %>%
    intersect(rownames(ligand_target_matrix))

  background <- intersect(expr_receiver, rownames(ligand_target_matrix))
  if (length(geneset) < 10 || length(background) < 50) next

  activities <- predict_ligand_activities(
    geneset = geneset,
    background_expressed_genes = background,
    ligand_target_matrix = ligand_target_matrix,
    potential_ligands = potential_ligands_focused
  ) %>%
    dplyr::arrange(dplyr::desc(aupr_corrected)) %>%
    dplyr::mutate(tissue = tiss, rank = dplyr::row_number())

  per_cancer_results[[tiss]] <- activities
  rm(mye_tiss); gc()
}

if (length(per_cancer_results) < 1) stop("No tissues produced NicheNet results. Check thresholds and metadata.")

all_activities <- dplyr::bind_rows(per_cancer_results)
write.csv(all_activities, file.path(table_dir, "nichenet_ligand_activities_per_cancer.csv"), row.names = FALSE)

# ============================================================
# Figure 6: Heatmap (top ligands per tissue)
# ============================================================
top_n <- 30

top_by_tissue <- all_activities %>%
  dplyr::group_by(tissue) %>%
  dplyr::slice_max(order_by = aupr_corrected, n = top_n, with_ties = FALSE) %>%
  dplyr::ungroup()

mat <- top_by_tissue %>%
  dplyr::select(tissue, test_ligand, aupr_corrected) %>%
  tidyr::pivot_wider(names_from = tissue, values_from = aupr_corrected, values_fill = 0)

ligands <- mat$test_ligand
hm <- as.matrix(mat[, setdiff(colnames(mat), "test_ligand"), drop = FALSE])
rownames(hm) <- ligands

png(file.path(fig_dir, "08_nichenet_ligand_activity_heatmap.png"), width = 1400, height = 1000, res = 150)
pheatmap::pheatmap(
  hm,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  main = "NicheNet ligand activity (per cancer; top ligands per tissue)"
)
dev.off()

# ============================================================
# Figure 7: Conserved ligands summary
# ============================================================
conserved <- all_activities %>%
  dplyr::group_by(test_ligand) %>%
  dplyr::summarise(
    n_tissues = dplyr::n_distinct(tissue),
    mean_activity = mean(aupr_corrected, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(n_tissues), dplyr::desc(mean_activity))

write.csv(conserved, file.path(table_dir, "conserved_ligands_summary.csv"), row.names = FALSE)

p_cons <- conserved %>%
  dplyr::slice_head(n = 30) %>%
  ggplot(aes(x = reorder(test_ligand, mean_activity), y = mean_activity)) +
  geom_col() +
  coord_flip() +
  ggtitle("Top ligands by mean activity (summary)") +
  xlab(NULL) + ylab("Mean AUPR (corrected)")

ggsave(file.path(fig_dir, "09_conserved_ligands_barplot.png"), p_cons, width = 10, height = 8, dpi = 300)

# ============================================================
# Figure 8: KEGG enrichment (PV-high vs PV-low DE genes, pooled)
# ============================================================
message("KEGG enrichment (pooled PV-high vs PV-low) ...")

Idents(mye) <- mye$Perivascular_high
de_all <- FindMarkers(mye, ident.1 = TRUE, ident.2 = FALSE,
                      min.pct = 0.1, logfc.threshold = 0.25, verbose = FALSE)

sig <- de_all %>%
  as.data.frame() %>%
  dplyr::filter(p_val_adj <= 0.05, avg_log2FC >= 0.25) %>%
  rownames()

if (length(sig) >= 20) {
  entrez <- bitr(sig, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) %>%
    dplyr::pull(ENTREZID) %>%
    unique()

  ekegg <- enrichKEGG(gene = entrez, organism = "hsa", pvalueCutoff = 0.05)
  kegg_df <- as.data.frame(ekegg)
  write.csv(kegg_df, file.path(table_dir, "kegg_enrichment_pvhigh.csv"), row.names = FALSE)

  p_kegg <- dotplot(ekegg, showCategory = 15) + ggtitle("KEGG enrichment (PV-high vs PV-low)")
  ggsave(file.path(fig_dir, "10_kegg_enrichment_dotplot.png"), p_kegg, width = 10, height = 7, dpi = 300)
} else {
  message("Too few significant DE genes for KEGG enrichment (n = ", length(sig), ").")
}

message("Done.")
