############################################################
## COMPLETE INTEGRATED PIPELINE
## ICH Survivor Bulk RNA-seq Analysis
## Analyses: DEG + GSEA + GSVA + ISG Score + Outcomes
############################################################

# ── Reload session ────────────────────────────────────────────────
#load("my_session_0604.RData")
# ── Save session ──────────────────────────────────────────────────
#save.image(file = "my_session_0701.RData")
############################################################
## 0. PACKAGES
############################################################

suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(EnhancedVolcano)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(msigdbr)
  library(fgsea)
  library(GSVA)
  library(pheatmap)
  library(data.table)
  library(broom)
  library(gt)
  library(gtsummary)
  library(logistf)
  library(ggrepel)
  library(conflicted)
  library(effectsize)
  library(factoextra)
})

# ── Resolve ALL namespace conflicts upfront ───────────────────────
# dplyr vs other packages
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::rename)
conflicts_prefer(dplyr::mutate)
conflicts_prefer(dplyr::arrange)
conflicts_prefer(dplyr::slice)
conflicts_prefer(dplyr::count)
conflicts_prefer(dplyr::between)
conflicts_prefer(dplyr::first)
conflicts_prefer(dplyr::last)
conflicts_prefer(dplyr::lag)
conflicts_prefer(dplyr::lead)

# base vs lubridate vs GenomicRanges
conflicts_prefer(base::intersect)
conflicts_prefer(base::union)
conflicts_prefer(base::setdiff)

# stats vs other packages
conflicts_prefer(stats::filter)
conflicts_prefer(stats::lag)

# Verify no remaining conflicts
conflict_scout()
# Should return an empty list if all conflicts are resolved

############################################################
## 1. LOAD AND PREPARE METADATA
############################################################

meta_data <- read.csv(
  "/home/ea576/project/REASSESS/Final_Reassess_data_041726_cleaned.csv"
)

names(meta_data)[names(meta_data) == "ses1"] <- "sex"

meta_data <- meta_data %>%
  dplyr::relocate(sex, .after = idrc) %>%
  dplyr::mutate(sex = as.character(sex)) %>%
  {
    is_empty_sex <- is.na(.$sex) | .$sex == ""
    .$sex[is_empty_sex & .$sex_birth %in% c(1,"1")] <- "M"
    .$sex[is_empty_sex & .$sex_birth %in% c(2,"2")] <- "F"
    .
  } %>%
  dplyr::filter(rowSums(is.na(.)) != ncol(.)) %>%
  dplyr::filter(!(is.na(.[[1]]) | .[[1]] == "" |
                    trimws(.[[1]]) == ""))

rownames(meta_data) <- meta_data$idrc

############################################################
## 2. LOAD AND PROCESS COUNT DATA — PIPELINE 2
############################################################

count_data <- readr::read_csv(
  "/home/ea576/project/REASSESS_YCGA_Pipeline_2/DESeq2/gene_count_matrix.csv"
)

ids_to_check <- c(
  "317-009","234-002","216-004","157-005","131-006","21-394094",
  "134-022","21-394199","21-394335","134-020","502-002","101-001",
  "502-004","134-025","101-002","313-04","205-002","134-023",
  "131-009","134-030","205-012","134-024","241-001","174-002",
  "134-029","306-005","134-036","100-015","223-006","101-007",
  "225-003","234-004","122-003","101-008","174-012","323-008",
  "157-018","122-002","323-009","101-006","131-024","323-005",
  "243-014","134-044"
)

missing_ids    <- ids_to_check[!ids_to_check %in% meta_data$idrc]
matched_in_ori <- missing_ids[missing_ids %in% meta_data$ori_study_id]

lookup <- meta_data %>%
  dplyr::filter(ori_study_id %in% matched_in_ori) %>%
  dplyr::select(ori_study_id, idrc)

for (i in seq_len(nrow(lookup))) {
  col_to_rename <- grep(lookup$ori_study_id[i],
                        names(count_data), value = TRUE)
  if (length(col_to_rename) == 1)
    names(count_data)[names(count_data) == col_to_rename] <-
      lookup$idrc[i]
}

col_313 <- grep("313-04", names(count_data), value = TRUE)
if (length(col_313) == 1)
  names(count_data)[names(count_data) == col_313] <- "313-010"

pattern_matches <- sapply(names(count_data), function(col) {
  any(grepl(col, meta_data$idrc, fixed = TRUE)) ||
    any(grepl(gsub("[^0-9]","",col),
              gsub("[^0-9]","",meta_data$idrc)))
})
pattern_missing       <- names(count_data)[!pattern_matches]
count_data_filtered   <- count_data[
  , c(TRUE, !names(count_data)[-1] %in% pattern_missing)
]

keep_ensembl_ids <- c(
  "ENSG00000138641.18","ENSG00000015479.20","ENSG00000182487.13",
  "ENSG00000290871.1", "ENSG00000285437.2", "ENSG00000160408.16",
  "ENSG00000183604.17"
)

filtered_df <- count_data_filtered %>%
  tibble::column_to_rownames(var = "gene_id") %>%
  { .[which(rowSums(.) > 10), ] } %>%
  as.data.frame() %>%
  cbind(RowName = rownames(.), .) %>%
  { rownames(.) <- NULL; . } %>%
  tidyr::separate(RowName,
                  into = c("Ensembl_ID","GENE_ID"), sep = "\\|")

duplicate_rows <- filtered_df[
  duplicated(filtered_df$GENE_ID) |
    duplicated(filtered_df$GENE_ID, fromLast = TRUE), ]

filtered_df <- filtered_df %>%
  dplyr::filter(!(GENE_ID %in% duplicate_rows$GENE_ID) |
                  Ensembl_ID %in% keep_ensembl_ids)

rownames(filtered_df) <- filtered_df$GENE_ID
filtered_df$GENE_ID   <- NULL
filtered_df           <- cbind(RowName = rownames(filtered_df),
                               filtered_df)
rownames(filtered_df) <- NULL
filtered_df$Ensembl_ID <- paste(filtered_df$Ensembl_ID,
                                filtered_df$RowName, sep = "|")
rownames(filtered_df)  <- filtered_df$Ensembl_ID
filtered_df$RowName    <- NULL
filtered_df$Ensembl_ID <- NULL

############################################################
## 3. LOAD AND PROCESS COUNT DATA — PIPELINE 1
############################################################

count_data_old <- readr::read_csv(
  "/gpfs/gibbs/project/sansing/ea576/REASSESS_YCGA_Pipeline/DESeq2/gene_count_matrix.csv"
)

colnames(count_data_old)[16] <- "208-001"
colnames(count_data_old)[17] <- "225-001"
colnames(count_data_old)[18] <- "225-002"
colnames(count_data_old)[27] <- "320-001"

filtered_df_2 <- count_data_old %>%
  tibble::column_to_rownames(var = "gene_id") %>%
  { .[which(rowSums(.) > 10), ] } %>%
  as.data.frame() %>%
  cbind(RowName = rownames(.), .) %>%
  { rownames(.) <- NULL; . } %>%
  tidyr::separate(RowName,
                  into = c("Ensembl_ID","GENE_ID"), sep = "\\|")

duplicate_rows_2 <- filtered_df_2[
  duplicated(filtered_df_2$GENE_ID) |
    duplicated(filtered_df_2$GENE_ID, fromLast = TRUE), ]

filtered_df_2 <- filtered_df_2 %>%
  dplyr::filter(!(GENE_ID %in% duplicate_rows_2$GENE_ID) |
                  Ensembl_ID %in% keep_ensembl_ids)

rownames(filtered_df_2) <- filtered_df_2$GENE_ID
filtered_df_2$GENE_ID   <- NULL
filtered_df_2           <- cbind(RowName = rownames(filtered_df_2),
                                 filtered_df_2)
rownames(filtered_df_2) <- NULL
filtered_df_2$Ensembl_ID <- paste(filtered_df_2$Ensembl_ID,
                                  filtered_df_2$RowName, sep = "|")
rownames(filtered_df_2)  <- filtered_df_2$Ensembl_ID
filtered_df_2$RowName    <- NULL
filtered_df_2$Ensembl_ID <- NULL

############################################################
## 4. COMBINE COUNT MATRICES
############################################################

combined_counts <- dplyr::full_join(
  filtered_df   %>% tibble::rownames_to_column("gene"),
  filtered_df_2 %>% tibble::rownames_to_column("gene"),
  by = "gene"
) %>%
  tibble::column_to_rownames("gene")

combined_counts[is.na(combined_counts)] <- 0

############################################################
## 5. BUILD MINI METADATA
############################################################

cols_keep <- c(
  "idrc","redcap_data_access_group","ori_study","ori_study_id",
  "ich_diag_date","initial_gcs","hist_mrs","initial_nihss",
  "age_ich","sex_birth","followup_mrs","mrs_date",
  "adm_date_ind_ich","disch_date_ind_ich","len_stay",
  "redcap_event_name","common_date","weight_lb_fu","adl_total",
  "tics_total","mrs","barthel_total","euroqol_mobility",
  "euroqol_selfcare","euroqol_usualact","euroqol_paindisc",
  "euroqol_anxdepr","hlthstat","iqdcode_avg_pre","iqdcode_avg_post",
  "story","dgs01_qsorres","dgs02_qsorres","digord","delay",
  "dateexam","height_ft","height_in","weight_lb","sbp","dbp",
  "nih_date","nih_total","content","fluency","aud_verbal",
  "seq_command","repetition","obj_naming","aphasia_sum",
  "bed_aphasia","apraxia","aphas_class","mmse13_qsorres",
  "age","sex","ses4","hx1","hx2","hx3","hx21","hx21d",
  "dm1","liv1","liv2","ich_hemisphere_adj","ich_location_adj"
)

mini_data <- meta_data[, intersect(cols_keep, names(meta_data))]

# ── Core recoding ─────────────────────────────────────────────────
mini_data <- mini_data %>%
  dplyr::mutate(
    hilo_mrs = dplyr::case_when(
      mrs <= 3 ~ "Good",
      mrs >  3 ~ "Poor",
      TRUE     ~ NA_character_
    ),
    hilo_mrs = factor(hilo_mrs, levels = c("Good", "Poor")),
    sex = dplyr::case_when(
      sex == "M" ~ "Male",
      sex == "F" ~ "Female",
      TRUE       ~ NA_character_
    ),
    sex = factor(sex, levels = c("Male", "Female")),
    age = dplyr::case_when(
      idrc == "134-029" ~ 72, idrc == "134-044" ~ 77,
      idrc == "243-014" ~ 53, idrc == "234-004" ~ 52,
      TRUE ~ age
    ),
    # MASTER eurodepr encoding — defined ONCE here, never changed
    eurodepr = dplyr::case_when(
      euroqol_anxdepr == 1        ~ "No",
      euroqol_anxdepr %in% c(2,3) ~ "Yes",
      TRUE                        ~ NA_character_
    ),
    eurodepr = factor(eurodepr, levels = c("No", "Yes"))
  )

rownames(mini_data) <- mini_data$idrc

############################################################
## 6. EXCLUSIONS
############################################################


############################################################
## 7. ALIGN METADATA TO COUNTS → combined_meta
############################################################

combined_meta <- data.frame(idrc = colnames(combined_counts)) %>%
  dplyr::left_join(mini_data, by = "idrc") %>%
  tibble::column_to_rownames("idrc")

stopifnot(all(rownames(combined_meta) == colnames(combined_counts)))

############################################################
## 8. FINAL COMBINED_META PREPARATION
############################################################

combined_meta <- combined_meta %>%
  dplyr::mutate(
    sex = factor(sex, levels = c("Male", "Female")),
    hilo_mrs = dplyr::case_when(
      mrs <= 3 ~ "Good",
      mrs >  3 ~ "Poor",
      TRUE     ~ NA_character_
    ),
    hilo_mrs = factor(hilo_mrs, levels = c("Good", "Poor")),
    ich_location_adj = factor(ich_location_adj,
                              levels = c("Deep", "Lobar")),
    ich_hemisphere_adj = factor(ich_hemisphere_adj,
                                levels = c("Left", "Right")),
    sbp  = as.numeric(sbp),
    dbp  = as.numeric(dbp),
    age  = as.numeric(age),
    mrs  = as.numeric(mrs),
    # time since ICH in years
    time_since_ich = as.numeric(
      as.Date(common_date) - as.Date(ich_diag_date)
    ) / 365.25
  )

cat("=== Final combined_meta ===\n")
cat("N patients:", nrow(combined_meta), "\n")
print(table(combined_meta$sex,      useNA = "always"))
print(table(combined_meta$hilo_mrs, useNA = "always"))
print(table(combined_meta$eurodepr, useNA = "always"))
###################################### Done ######################################

################################################################################
# ICH Survivor Bulk RNA-seq Analysis
# Integrated Pipeline: DEG + GSEA + GSVA + ISG Score + Outcomes
#
# This script covers:
#   Section 9  — DESeq2 base setup and low-count filtering
#   Section 10 — VST normalization
#   Section 11 — Helper functions (DESeq2, GSEA, plotting)
#   Section 12 — Gene set collections + covariate scaling
#   Section 12b — GSVA sample-level pathway activity
################################################################################


################################################################################
# SECTION 9: DESeq2 SETUP
#
# Builds a DESeqDataSet from the combined count matrix and metadata,
# then applies a minimum count filter to remove lowly expressed genes.
################################################################################

dds_full <- DESeqDataSetFromMatrix(
  countData = combined_counts,
  colData   = combined_meta,
  design    = ~ 1            # intercept-only; design is set per contrast later
)

# Retain genes with >= 10 counts in at least 3 samples
keep     <- rowSums(counts(dds_full) >= 10) >= 3
dds_full <- dds_full[keep, ]

cat("Genes after low-count filtering:", nrow(dds_full), "\n")


################################################################################
# SECTION 10: VST NORMALIZATION
#
# Variance-stabilizing transformation (VST) is computed once and reused
# throughout the pipeline to ensure consistency.
#
#   vst_mat_blind — blind = TRUE  → used for GSVA and unsupervised analyses
#   vst_mat       — alias for vst_mat_blind; used throughout unless noted
################################################################################

vst_mat_blind <- assay(vst(dds_full, blind = TRUE))
vst_mat       <- vst_mat_blind

cat("VST matrix dimensions (genes x samples):", dim(vst_mat), "\n")


################################################################################
# SECTION 11: HELPER FUNCTIONS
#
# Reusable functions for gene name cleaning, DESeq2 result formatting,
# ranked-list construction, GSEA execution, and NES bar chart plotting.
################################################################################

# ------------------------------------------------------------------------------
# clean_gene_names()
#
# Converts ENSEMBL|SYMBOL rownames to plain SYMBOL rownames.
# If duplicates remain after stripping, expression values are averaged.
#
# Args:
#   mat  — numeric matrix with rownames in "ENSG...|SYMBOL" format
#
# Returns:
#   Matrix with SYMBOL rownames; duplicates resolved by row-wise mean
# ------------------------------------------------------------------------------
clean_gene_names <- function(mat) {
  
  # Strip ENSEMBL prefix if rownames contain a pipe character
  if (any(grepl("\\|", rownames(mat)))) {
    gene_df <- tibble::tibble(raw = rownames(mat)) %>%
      tidyr::separate(raw, into = c("ENSEMBL", "SYMBOL"),
                      sep = "\\|", fill = "right")
    rownames(mat) <- gene_df$SYMBOL
  }
  
  # Remove rows with missing or empty gene symbols
  mat <- mat[!is.na(rownames(mat)) & rownames(mat) != "", ]
  
  # Resolve duplicates by averaging expression across duplicate gene symbols
  if (any(duplicated(rownames(mat)))) {
    mat <- mat %>%
      as.data.frame() %>%
      tibble::rownames_to_column("gene") %>%
      dplyr::group_by(gene) %>%
      dplyr::summarize(dplyr::across(everything(), mean), .groups = "drop") %>%
      tibble::column_to_rownames("gene") %>%
      as.matrix()
  }
  
  mat
}


# ------------------------------------------------------------------------------
# map_ensembl_to_symbol()
#
# Maps ENSEMBL IDs from DESeq2 results to HGNC gene symbols and Entrez IDs
# using the org.Hs.eg.db annotation database. One mapping per ENSEMBL ID.
#
# Args:
#   dds_obj — DESeqDataSet after running DESeq()
#
# Returns:
#   Data frame with columns: ENSEMBL, SYMBOL, ENTREZID
# ------------------------------------------------------------------------------
map_ensembl_to_symbol <- function(dds_obj) {
  clusterProfiler::bitr(
    sub("\\..*", "", rownames(results(dds_obj))),  # strip version suffix
    fromType = "ENSEMBL",
    toType   = c("SYMBOL", "ENTREZID"),
    OrgDb    = org.Hs.eg.db
  ) %>%
    dplyr::group_by(ENSEMBL) %>%
    dplyr::slice(1) %>%           # keep first match per ENSEMBL ID
    dplyr::ungroup()
}


# ------------------------------------------------------------------------------
# build_res_df()
#
# Extracts DESeq2 results and appends gene symbol annotations.
# Rownames are expected in "ENSG...|SYMBOL" format.
# Final gene symbol priority: GENE_ID (from rowname) > SYMBOL (from bitr)
#                             > ENSEMBL (fallback)
#
# Args:
#   dds_obj — DESeqDataSet after running DESeq()
#   ens2sym — data frame from map_ensembl_to_symbol()
#
# Returns:
#   Data frame with DESeq2 statistics and a `final_symbol` column
# ------------------------------------------------------------------------------
build_res_df <- function(dds_obj, ens2sym) {
  as.data.frame(results(dds_obj)) %>%
    tibble::rownames_to_column("RowName") %>%
    tidyr::separate(RowName,
                    into = c("Ensembl_ID", "GENE_ID"),
                    sep  = "\\|", fill = "right") %>%
    dplyr::mutate(Ensembl_clean = sub("\\..*", "", Ensembl_ID)) %>%
    dplyr::left_join(ens2sym, by = c("Ensembl_clean" = "ENSEMBL")) %>%
    dplyr::mutate(
      final_symbol = dplyr::case_when(
        !is.na(GENE_ID) & GENE_ID != "" ~ GENE_ID,
        !is.na(SYMBOL)                  ~ SYMBOL,
        TRUE                            ~ Ensembl_clean
      )
    )
}


# ------------------------------------------------------------------------------
# build_ranked_list()
#
# Constructs a named, sorted vector of DESeq2 Wald statistics for GSEA input.
# Genes without a test statistic or gene symbol are excluded.
# Duplicate gene symbols are removed (first occurrence retained after sorting).
#
# Args:
#   res_df — data frame from build_res_df()
#
# Returns:
#   Named numeric vector sorted descending by Wald statistic
# ------------------------------------------------------------------------------
build_ranked_list <- function(res_df) {
  res_df %>%
    dplyr::filter(!is.na(stat), !is.na(final_symbol)) %>%
    dplyr::arrange(dplyr::desc(stat)) %>%
    dplyr::distinct(final_symbol, .keep_all = TRUE) %>%
    dplyr::pull(stat, name = final_symbol)
}


# ------------------------------------------------------------------------------
# run_gsea()
#
# Runs gene set enrichment analysis using clusterProfiler::GSEA.
# Seed is fixed at 42 for reproducibility.
#
# Args:
#   ranked_list — named numeric vector from build_ranked_list()
#   term2gene   — two-column data frame: gene set name and gene symbol
#
# Returns:
#   gseaResult object
# ------------------------------------------------------------------------------
run_gsea <- function(ranked_list, term2gene) {
  clusterProfiler::GSEA(
    geneList     = ranked_list,
    TERM2GENE    = term2gene,
    pvalueCutoff = 1,     # retain all results; filter downstream
    eps          = 0,
    seed         = 42
  )
}


# ------------------------------------------------------------------------------
# plot_nes_barchart()
#
# Plots a horizontal bar chart of normalized enrichment scores (NES) for
# significantly enriched pathways (adjusted p < 0.05), annotated with
# q-values. Top `n` pathways by absolute NES are displayed.
#
# Args:
#   gsea_obj    — gseaResult object from run_gsea()
#   title_str   — plot title string
#   subtitle_str— plot subtitle string
#   pos_label   — label for positively enriched pathways (default: "Up")
#   neg_label   — label for negatively enriched pathways (default: "Down")
#   n           — maximum number of pathways to display (default: 20)
#
# Returns:
#   ggplot object
# ------------------------------------------------------------------------------
plot_nes_barchart <- function(gsea_obj,
                              title_str,
                              subtitle_str,
                              pos_label = "Up",
                              neg_label = "Down",
                              n = 20) {
  
  top <- gsea_obj@result %>%
    dplyr::filter(p.adjust < 0.05) %>%
    dplyr::mutate(
      direction   = ifelse(NES > 0, pos_label, neg_label),
      Description = Description %>%
        gsub("HALLMARK_", "", .) %>%
        gsub("_", " ", .)        %>%
        stringr::str_to_title()
    ) %>%
    dplyr::slice_max(abs(NES), n = n)
  
  ggplot(top, aes(x = reorder(Description, NES), y = NES, fill = direction)) +
    geom_col() +
    geom_text(
      aes(label = paste0("q=", round(p.adjust, 3))),
      hjust  = ifelse(top$NES > 0, -0.1, 1.1),
      size   = 4,
      colour = "black"
    ) +
    scale_fill_manual(values = setNames(
      c("#7E1811", "#456787"),
      c(pos_label, neg_label)
    )) +
    scale_y_continuous(
      limits = c(min(top$NES) * 1.1, max(top$NES) * 1.3)
    ) +
    coord_flip() +
    theme_bw() +
    theme(
      panel.grid    = element_blank(),
      axis.text     = element_text(size = 14, colour = "black"),
      axis.title.x  = element_text(size = 16, colour = "black"),
      axis.title.y  = element_text(colour = "black"),
      legend.text   = element_text(size = 12, colour = "black"),
      legend.title  = element_text(colour = "black"),
      plot.title    = element_text(size = 16, colour = "black"),
      plot.subtitle = element_text(colour = "black"),
      strip.text    = element_text(colour = "black"),
      text          = element_text(colour = "black")
    ) +
    labs(
      title    = title_str,
      subtitle = subtitle_str,
      x        = NULL,
      y        = "Normalized Enrichment Score (NES)",
      fill     = "Direction"
    )
}

# ------------------------------------------------------------------------------
# build_ranked_list_tidy()
#
# Constructs a named, sorted vector of Wald statistics for GSEA input from
# a tidy DESeq2 result (extracted with tidy = TRUE). Gene IDs are in the
# "ENSG...|SYMBOL" format produced by DESeq2 on this dataset; the SYMBOL
# portion after the pipe is used as the gene name.
#
# Genes without a test statistic or resolvable symbol are excluded.
# Duplicate symbols are removed, keeping the first occurrence after sorting.
#
# Args:
#   res_tidy — tidy DESeq2 result data frame with columns: gene_id, statistic
#
# Returns:
#   Named numeric vector sorted descending by Wald statistic
# ------------------------------------------------------------------------------
build_ranked_list_tidy <- function(res_tidy) {
  
  # Confirm expected columns are present before proceeding
  required_cols <- c("gene_id", "stat")
  missing_cols  <- setdiff(required_cols, colnames(res_tidy))
  if (length(missing_cols) > 0) {
    stop("build_ranked_list_tidy(): missing required columns: ",
         paste(missing_cols, collapse = ", "),
         "\nColumns found: ", paste(colnames(res_tidy), collapse = ", "))
  }
  
  res_tidy %>%
    dplyr::filter(!is.na(.data[["stat"]])) %>%
    dplyr::mutate(
      symbol = dplyr::case_when(
        grepl("\\|", .data[["gene_id"]]) ~ sub(".*\\|", "", .data[["gene_id"]]),
        TRUE                             ~ .data[["gene_id"]]
      )
    ) %>%
    dplyr::filter(
      !is.na(.data[["symbol"]]),
      .data[["symbol"]] != "",
      !grepl("^ENSG", .data[["symbol"]])
    ) %>%
    dplyr::arrange(dplyr::desc(.data[["stat"]])) %>%
    dplyr::distinct(.data[["symbol"]], .keep_all = TRUE) %>%
    dplyr::pull(.data[["stat"]], name = .data[["symbol"]])
}

################################################################################
# SECTION 12: GENE SET COLLECTIONS
#
# Downloads MSigDB gene sets via msigdbr for use in GSEA.
# Active collections: Hallmark (H) and KEGG (C2).
# Immunologic (C7/IMMUNESIGDB) is available but disabled by default.
################################################################################

hallmark_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, gene_symbol)

t2g_kegg <- msigdbr(species = "Homo sapiens",
                    category    = "C2",
                    subcategory = "CP:KEGG") %>%
  dplyr::select(gs_name, gene_symbol)

# t2g_immune <- msigdbr(species = "Homo sapiens",
#                       category    = "C7",
#                       subcategory = "IMMUNESIGDB") %>%
#   dplyr::select(gs_name, gene_symbol)

gene_set_collections <- list(
  Hallmark = hallmark_t2g,
  KEGG     = t2g_kegg
  # Immunologic = t2g_immune
)


################################################################################
# SECTION 12b: COVARIATE SCALING
#
# Continuous covariates are z-score standardized once here so that
# DESeq2 models throughout the pipeline use consistent scaled predictors.
# Scaling diagnostics are printed to confirm mean ≈ 0, SD ≈ 1.
################################################################################

combined_meta <- combined_meta %>%
  dplyr::mutate(
    age_scaled     = as.numeric(scale(age)),
    age_ich_scaled = as.numeric(scale(age_ich)),
    nihss_scaled   = as.numeric(scale(initial_nihss)),
    sbp_scaled     = as.numeric(scale(sbp)),
    dbp_scaled     = as.numeric(scale(dbp))
  )

# Scaling diagnostics
cat("age_scaled     — mean:", round(mean(combined_meta$age_scaled,     na.rm = TRUE), 4),
    " SD:", round(sd(combined_meta$age_scaled,     na.rm = TRUE), 4), "\n")
cat("age_ich_scaled — mean:", round(mean(combined_meta$age_ich_scaled, na.rm = TRUE), 4),
    " SD:", round(sd(combined_meta$age_ich_scaled, na.rm = TRUE), 4), "\n")
cat("nihss_scaled   — mean:", round(mean(combined_meta$nihss_scaled,   na.rm = TRUE), 4),
    " SD:", round(sd(combined_meta$nihss_scaled,   na.rm = TRUE), 4), "\n")


################################################################################
# SECTION 12c: GSVA — SAMPLE-LEVEL PATHWAY ACTIVITY
#
# Computes per-sample Hallmark pathway enrichment scores using GSVA.
# Input: blind VST matrix (gene symbols as rownames).
# Output: heatmap of pathway scores annotated by sex, with sample clusters
#         determined by hierarchical clustering of the GSVA score matrix.
#
# Cluster count k is selected empirically using elbow and silhouette plots.
# Set k below after reviewing those plots.
################################################################################

# --- Step 1: Build named gene-set list from hallmark_t2g ---------------------

pathways_list <- hallmark_t2g %>%
  dplyr::group_by(gs_name) %>%
  dplyr::summarise(genes = list(gene_symbol), .groups = "drop") %>%
  tibble::deframe()


# --- Step 2: Map VST matrix rownames to gene symbols -------------------------
# Rownames are in "ENSG...|SYMBOL" format; extract the symbol after the pipe.
# Rows that still begin with "ENSG" after splitting are unmapped and dropped.

rowname_map <- tibble::tibble(RowName = rownames(vst_mat_blind)) %>%
  dplyr::mutate(
    Ensembl_ID   = sub("\\|.*", "", RowName),
    final_symbol = sub(".*\\|", "", RowName)
  ) %>%
  dplyr::filter(
    !is.na(final_symbol),
    final_symbol != "",
    !grepl("^ENSG", final_symbol)   # drop rows where symbol was not resolved
  )

vst_df <- as.data.frame(vst_mat_blind) %>%
  tibble::rownames_to_column("RowName") %>%
  dplyr::left_join(rowname_map, by = "RowName") %>%
  dplyr::filter(!is.na(final_symbol), final_symbol != "") %>%
  dplyr::distinct(final_symbol, .keep_all = TRUE) %>%   # one row per symbol
  tibble::column_to_rownames("final_symbol") %>%
  dplyr::select(-RowName, -Ensembl_ID)

vst_matrix_gsva <- as.matrix(vst_df)
cat("VST matrix for GSVA:", nrow(vst_matrix_gsva), "genes x",
    ncol(vst_matrix_gsva), "samples\n")


# --- Step 3: Run GSVA --------------------------------------------------------

gsva_param <- GSVA::gsvaParam(
  exprData = vst_matrix_gsva,
  geneSets = pathways_list,
  minSize  = 15,
  maxSize  = 500
)

gsva_res <- GSVA::gsva(gsva_param, verbose = TRUE)
cat("GSVA result dimensions (pathways x samples):", dim(gsva_res), "\n")


# --- Step 4: Pathway activity summary ----------------------------------------

pathway_summary <- data.frame(
  pathway    = rownames(gsva_res),
  mean_score = rowMeans(gsva_res),
  sd_score   = apply(gsva_res, 1, sd),
  cv         = apply(gsva_res, 1, sd) / abs(rowMeans(gsva_res))
) %>%
  dplyr::arrange(dplyr::desc(abs(mean_score)))

cat("\nTop 10 most active pathways:\n")
print(head(pathway_summary, 10))


# --- Step 5: Build column annotation for heatmap ----------------------------

annot_df <- combined_meta %>%
  tibble::rownames_to_column("idrc") %>%
  dplyr::filter(idrc %in% colnames(gsva_res)) %>%
  dplyr::select(idrc, sex) %>%
  tibble::column_to_rownames("idrc")

annot_df <- annot_df[colnames(gsva_res), , drop = FALSE]   # match column order

col_annot <- data.frame(
  Sex       = as.character(annot_df$sex),
  row.names = rownames(annot_df)
)
col_annot[is.na(col_annot)] <- "Unknown"   # pheatmap requires no NAs

ann_colors <- list(
  Sex = c("Male" = "#364981", "Female" = "#FFC482")
)

cat("Sex distribution in annotation:\n")
print(table(col_annot$Sex, useNA = "always"))


# --- Step 6: Clean pathway names for display ---------------------------------

gsva_plot_mat <- gsva_res
rownames(gsva_plot_mat) <- rownames(gsva_plot_mat) %>%
  gsub("HALLMARK_", "", .) %>%
  gsub("_", " ", .)         %>%
  stringr::str_to_title()


# --- Step 7: Initial heatmap (without cluster borders) -----------------------

p_gsva_heatmap <- pheatmap::pheatmap(
  gsva_plot_mat,
  annotation_col           = col_annot,
  annotation_colors        = ann_colors,
  scale                    = "row",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method        = "complete",
  show_colnames            = FALSE,
  fontsize_row             = 11,
  color                    = colorRampPalette(c("#456787", "white", "#7E1811"))(100),
  main                     = "GSVA Hallmark Pathway Enrichment — All Samples"
)

# To save:
# png("GSVA_hallmark_all_samples.png", width = 14, height = 8, units = "in", res = 300)
# print(p_gsva_heatmap)
# dev.off()

# cairo_pdf("GSVA_hallmark_all_samples.pdf", width = 14, height = 8)
# print(p_gsva_heatmap)
# dev.off()


# --- Step 8: Determine optimal number of clusters (k) -----------------------
# Extract the column dendrogram used by pheatmap, then assess k using
# the elbow (WSS) and silhouette methods. Set k below after reviewing plots.

col_dend <- p_gsva_heatmap$tree_col
cat("Clustering method:", col_dend$method, "\n")    # should be "complete"
cat("Number of samples:", length(col_dend$order), "\n")

# Elbow plot
fviz_nbclust(
  t(gsva_res),
  FUNcluster = hcut,
  method     = "wss",
  k.max      = 10,
  hc_func    = "hclust",
  hc_metric  = "euclidean",
  hc_method  = "complete"
) + labs(title = "Elbow Method — Optimal k (Euclidean + Complete)")

# Silhouette plot
fviz_nbclust(
  t(gsva_res),
  FUNcluster = hcut,
  method     = "silhouette",
  k.max      = 10,
  hc_func    = "hclust",
  hc_metric  = "euclidean",
  hc_method  = "complete"
) + labs(title = "Silhouette Method — Optimal k (Euclidean + Complete)")


# --- Step 9: Cut dendrogram and assign sample clusters -----------------------
# Set k after reviewing elbow and silhouette plots above.

k               <- 2
sample_clusters <- cutree(col_dend, k = k)

cat("\nCluster sizes (k =", k, "):\n")
print(table(sample_clusters))


# --- Step 10: Redraw heatmap with cluster borders ----------------------------

p_gsva_heatmap <- pheatmap::pheatmap(
  gsva_plot_mat,
  annotation_col           = col_annot,
  annotation_colors        = ann_colors,
  scale                    = "row",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method        = "complete",
  cutree_cols              = k,
  show_colnames            = FALSE,
  fontsize_row             = 11,
  color                    = colorRampPalette(c("#456787", "white", "#7E1811"))(100),
  main                     = paste0("GSVA Hallmark Pathway Enrichment — ", k, " Clusters")
)

# To save:
# png("GSVA_hallmark_clustered.png", width = 14, height = 8, units = "in", res = 300)
# print(p_gsva_heatmap)
# dev.off()

################################################################################
# SECTION 12d: GSVA CLUSTER — METADATA ASSOCIATION TESTS
#
# Tests whether key clinical and demographic variables differ between
# GSVA-derived sample clusters using chi-square (categorical) or Wilcoxon
# rank-sum (continuous) tests. Results are visualized as a -log10(p) bar plot
# faceted by variable type.
#
# Requires: sample_clusters (from Section 12c), combined_meta
################################################################################

# Variables to test for association with GSVA cluster membership
keep_vars <- c("sex", "age", "race", "ethnicity", "mrs", "euroqol_anxdepr",
               "cesd_total")

# Variable type classification
# - categorical_cols: tested with chi-square
# - force_continuous: coerced to numeric and tested with Wilcoxon
# - euroqol_anxdepr is treated as categorical by default; move to
#   force_continuous if an ordinal/numeric treatment is preferred
categorical_cols <- c("sex", "race", "ethnicity")
force_continuous <- c("age", "mrs", "cesd_total")

# Attach cluster assignments to metadata
meta_cluster <- combined_meta %>%
  dplyr::mutate(
    sample  = rownames(.),
    cluster = factor(sample_clusters[rownames(.)])
  ) %>%
  dplyr::filter(!is.na(cluster))


# --- Run association tests ---------------------------------------------------

results_list <- list()

for (col in keep_vars) {
  
  if (!col %in% colnames(meta_cluster)) {
    message("Skipping ", col, " — column not found in metadata")
    next
  }
  
  var <- meta_cluster[[col]]
  
  if (col %in% force_continuous) {
    var <- suppressWarnings(as.numeric(var))
  }
  
  # Skip variables with no usable variation
  if (all(is.na(var)) || length(unique(na.omit(var))) < 2) {
    message("Skipping ", col, " — insufficient variation")
    next
  }
  
  tryCatch({
    
    if (col %in% categorical_cols ||
        (is.character(var) && !col %in% force_continuous)) {
      
      tbl  <- table(meta_cluster$cluster, var, useNA = "no")
      test <- chisq.test(tbl, simulate.p.value = TRUE, B = 2000)
      results_list[[col]] <- data.frame(
        variable = col,
        type     = "categorical",
        test     = "Chi-square",
        pvalue   = test$p.value,
        n_nonNA  = sum(!is.na(var))
      )
      
    } else {
      
      test <- wilcox.test(as.numeric(var) ~ meta_cluster$cluster, exact = FALSE)
      results_list[[col]] <- data.frame(
        variable = col,
        type     = "continuous",
        test     = "Wilcoxon",
        pvalue   = test$p.value,
        n_nonNA  = sum(!is.na(var))
      )
    }
    
  }, error = function(e) {
    message("Skipping ", col, ": ", e$message)
  })
}


# --- Compile and display results ---------------------------------------------

# Human-readable axis labels for the bar plot
readable_labels <- c(
  sex             = "Sex",
  age             = "Age",
  race            = "Race",
  ethnicity       = "Ethnicity",
  mrs             = "mRS",
  euroqol_anxdepr = "EQ-5D Anxiety/Depression",
  cesd_total      = "CESD Total"
)

meta_test_results <- dplyr::bind_rows(results_list) %>%
  dplyr::mutate(padj = p.adjust(pvalue, method = "BH")) %>%
  dplyr::arrange(pvalue)

print(meta_test_results)


# --- Bar plot: -log10(p) by variable, faceted by type -----------------------

p_cluster_meta <- meta_test_results %>%
  dplyr::mutate(
    sig    = dplyr::case_when(
      padj   < 0.05 ~ "padj < 0.05",
      pvalue < 0.05 ~ "nominal p < 0.05",
      TRUE          ~ "ns"
    ),
    log10p = -log10(pvalue),
    type   = factor(type, levels = c("continuous", "categorical"))
  ) %>%
  ggplot(aes(x = log10p, y = reorder(variable, log10p), fill = sig)) +
  geom_col() +
  geom_vline(xintercept = -log10(0.05),
             linetype = "dashed", colour = "red", linewidth = 0.7) +
  scale_fill_manual(values = c(
    "padj < 0.05"      = "firebrick",
    "nominal p < 0.05" = "#7E1811",
    "ns"               = "grey75"
  )) +
  scale_y_discrete(labels = readable_labels) +
  facet_grid(type ~ ., scales = "free_y", space = "free_y") +
  theme_bw() +
  theme(
    text             = element_text(family = "Arial", colour = "black"),
    plot.title       = element_text(size = 14),
    plot.subtitle    = element_text(size = 12),
    axis.title       = element_text(size = 13),
    axis.text        = element_text(size = 12),
    legend.title     = element_text(size = 12),
    legend.text      = element_text(size = 11),
    strip.text       = element_text(size = 12),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey90")
  ) +
  labs(
    title    = "Metadata variables associated with GSVA cluster membership",
    subtitle = "Red dashed line = nominal p < 0.05 | Faceted by variable type",
    x        = "-log10(p-value)",
    y        = NULL,
    fill     = "Significance"
  )

print(p_cluster_meta)

# To save:
# ggsave("GSVA_cluster_metadata_associations.png", plot = p_cluster_meta,
#        width = 12, height = 7, units = "in", dpi = 300)

# cairo_pdf("GSVA_cluster_metadata_associations.pdf", width = 12, height = 7)
# print(p_cluster_meta)
# dev.off()


################################################################################
# SECTION 13: DIFFERENTIAL EXPRESSION BY SEX (age-adjusted)
#
# Primary DESeq2 contrast: Female vs Male, adjusting for scaled age.
# Thresholds: adjusted p < 0.05, |log2FC| > 0.25.
# Outputs: DEG table, ranked gene list for GSEA, volcano plot.
################################################################################

# Retain samples with non-missing sex and scaled age
keep_sex <- !is.na(combined_meta$sex) & !is.na(combined_meta$age_scaled)

dds_sex <- DESeqDataSetFromMatrix(
  countData = combined_counts[rownames(dds_full), keep_sex],
  colData   = combined_meta[keep_sex, ],
  design    = ~ age_scaled + sex
)

dds_sex$sex <- relevel(dds_sex$sex, ref = "Male")
dds_sex     <- DESeq(dds_sex)

cat("Model coefficients:\n")
print(resultsNames(dds_sex))

# Annotate results and build ranked list for GSEA
ens2sym_sex <- map_ensembl_to_symbol(dds_sex)
res_df_sex  <- build_res_df(dds_sex, ens2sym_sex)
ranked_sex  <- build_ranked_list(res_df_sex)

cat("Total DEG rows:", nrow(res_df_sex),
    "| Genes with symbol:", sum(!is.na(res_df_sex$SYMBOL)), "\n")

# Significant DEGs
sig_genes_sex <- res_df_sex %>%
  dplyr::filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 0.25)

cat("Significant sex DEGs (padj < 0.05, |LFC| > 0.25):", nrow(sig_genes_sex), "\n")


# --- Volcano plot ------------------------------------------------------------
# Labels: top 65 significant DEGs by rank score, plus genes of interest.
# Colors: blue = higher in Male, red = higher in Female, grey = not significant.

selectLab_sex <- res_df_sex %>%
  dplyr::filter(padj < 0.05, abs(log2FoldChange) > 0.25) %>%
  dplyr::mutate(rank_score = -log10(padj) * abs(log2FoldChange)) %>%
  dplyr::arrange(dplyr::desc(rank_score)) %>%
  dplyr::slice_head(n = 65) %>%
  dplyr::filter(!grepl("^ENSG", final_symbol)) %>%
  dplyr::pull(final_symbol) %>%
  base::union(c("AHSP", "CA1", "OSBP2", "IFI27", "SRY"))

keyvals_sex <- dplyr::case_when(
  res_df_sex$log2FoldChange < -0.25 & res_df_sex$padj < 0.05 ~ "#456787",
  res_df_sex$log2FoldChange >  0.25 & res_df_sex$padj < 0.05 ~ "#7E1811",
  TRUE                                                         ~ "grey70"
)
names(keyvals_sex) <- keyvals_sex

p_volcano_sex <- EnhancedVolcano(
  res_df_sex,
  lab                   = res_df_sex$final_symbol,
  selectLab             = selectLab_sex,
  title                 = "Female vs Male DEGs (age-adjusted)",
  subtitle              = NULL,
  x                     = "log2FoldChange",
  y                     = "padj",
  pCutoff               = 0.05,
  FCcutoff              = 0.25,
  pointSize             = 4.0,
  labSize               = 5.5,
  colAlpha              = 1,
  colCustom             = keyvals_sex,
  drawConnectors        = TRUE,
  widthConnectors       = 0.5,
  colConnectors         = "grey40",
  maxoverlapsConnectors = Inf,
  min.segment.length    = 0,
  boxedLabels           = TRUE
) +
  theme_classic(base_size = 16) +
  theme(
    legend.position = "none",
    plot.title      = element_text(size = 25)
  )

print(p_volcano_sex)

# To save:
# ggsave("volcano_sex_age_adjusted.png", plot = p_volcano_sex,
#        width = 15, height = 7.5, units = "in", dpi = 300)

cairo_pdf("volcano_sex_age_adjusted.pdf", width = 15, height = 7.5)
print(p_volcano_sex)
dev.off()
################################################################################
# SECTION 14: DEPRESSION PREVALENCE BY SEX AND FUNCTIONAL OUTCOME
#
# Describes the distribution of EQ-5D-defined depression (eurodepr) by sex
# within each functional outcome group (Good / Poor mRS).
#
# Key finding: Among patients with poor outcome (mRS > 3), females had a
# significantly higher prevalence of depression than males
# (53.8% vs 13.0%; OR = 7.25, 95% CI: 1.21–57.87, p = 0.018, Fisher's test).
# This sex difference in depression within the poor outcome group may partly
# explain the divergent transcriptomic signatures observed in the DEG and
# GSEA analyses.
################################################################################
##############################################
## Examine prevalence of depression x sex and plot 
############################################
library(ggsignif)

# -------------------------
# 1. Run chi-square test first
# -------------------------
chisq_data <- combined_meta %>%
  dplyr::filter(!is.na(eurodepr), !is.na(sex)) %>%
  dplyr::mutate(eurodepr = factor(eurodepr), sex = factor(sex))

chisq_result <- chisq.test(table(chisq_data$sex, chisq_data$eurodepr))

p_val   <- chisq_result$p.value
p_label <- ifelse(p_val < 0.001, "p < 0.001",
                  ifelse(p_val < 0.01,  "p < 0.01",
                         ifelse(p_val < 0.05,  paste0("p = ", round(p_val, 3)),
                                paste0("p = ", round(p_val, 3)))))

cat("Chi-square p-value:", p_val, "\n")

# -------------------------
# 2. Plot with p-value bracket
# -------------------------
plot_data <- combined_meta %>%
  dplyr::filter(!is.na(eurodepr), !is.na(sex)) %>%
  dplyr::group_by(sex, eurodepr) %>%
  dplyr::summarise(n = n(), .groups = "drop") %>%
  dplyr::group_by(sex) %>%
  dplyr::mutate(
    pct   = round(n / sum(n) * 100, 1),
    total = sum(n)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::filter(eurodepr == "Yes")

depression_by_sex <- ggplot2::ggplot(plot_data, aes(x = sex, y = pct, fill = sex)) +
  ggplot2::geom_col(width = 0.5) +
  ggplot2::geom_text(
    aes(label = paste0(pct, "%\n(n=", n, "/", total, ")")),
    vjust    = -0.5,
    fontface = "bold",
    size     = 5
  ) +
  
  # ✅ p-value bracket
  ggsignif::geom_signif(
    comparisons      = list(c("Female", "Male")),
    annotations      = p_label,
    y_position       = 78,                        # ✅ adjust height of bracket
    tip_length       = 0.02,
    textsize         = 5,
    fontface         = "bold",
    color            = "black",
    vjust            = -0.3                        # ✅ pushes label above bracket
  ) +
  
  ggplot2::scale_fill_manual(
    values = c("Male" = "#364981", "Female" = "#FFC482")
  ) +
  ggplot2::scale_y_continuous(
    limits = c(0, 95),
    expand = expansion(mult = c(0, 0.05)),
    labels = scales::percent_format(scale = 1)
  ) +
  ggplot2::labs(
    title   = "Prevalence of Depression and Anxiety by Sex",
    caption = "EQ-5D item 5; Yes = moderate or severe",
    x       = NULL,
    y       = "Percentage Depressed/Anxious (%)"
  ) +
  ggplot2::theme_classic(base_size = 14) +
  ggplot2::theme(
    plot.title      = ggplot2::element_text(face = "bold", size = 15),
    plot.caption    = ggplot2::element_text(colour = "grey50", size = 10),
    axis.text       = ggplot2::element_text(colour = "black", size = 13),
    legend.position = "none"
  )

print(depression_by_sex)

# ggsave(
#   "depression_by_sex.png",
#   plot   = n,
#   width  = 8,
#   height = 6,
#   units  = "in"
# )

## save
# cairo_pdf("depression_by_sex.pdf", width = 8, height = 6)
# print(depression_by_sex)
# dev.off()

# --- Depression prevalence by sex: Poor mRS ----------------------------------

cat("Depression by sex — Poor mRS (mRS > 3):\n")
combined_meta %>%
  dplyr::filter(hilo_mrs == "Poor") %>%
  dplyr::group_by(sex) %>%
  dplyr::summarise(
    n_total       = dplyr::n(),
    n_depressed   = sum(eurodepr == "Yes", na.rm = TRUE),
    n_not_dep     = sum(eurodepr == "No",  na.rm = TRUE),
    n_NA          = sum(is.na(eurodepr)),
    pct_depressed = round(n_depressed / n_total * 100, 1),
    .groups       = "drop"
  ) %>%
  print()


# --- Depression prevalence by sex: Good mRS ----------------------------------

cat("Depression by sex — Good mRS (mRS <= 3):\n")
combined_meta %>%
  dplyr::filter(hilo_mrs == "Good") %>%
  dplyr::group_by(sex) %>%
  dplyr::summarise(
    n_total       = dplyr::n(),
    n_depressed   = sum(eurodepr == "Yes", na.rm = TRUE),
    n_not_dep     = sum(eurodepr == "No",  na.rm = TRUE),
    n_NA          = sum(is.na(eurodepr)),
    pct_depressed = round(n_depressed / n_total * 100, 1),
    .groups       = "drop"
  ) %>%
  print()


# --- Fisher's exact test: depression ~ sex within Poor mRS -------------------

poor_tab <- combined_meta %>%
  dplyr::filter(hilo_mrs == "Poor", !is.na(eurodepr)) %>%
  dplyr::select(sex, eurodepr) %>%
  table()

print(poor_tab)
fisher.test(poor_tab)


# --- Stacked bar plot: depression by sex in Poor mRS -------------------------
## Fisher p-value
ft       <- fisher.test(poor_tab)
pval_txt <- paste0("Fisher's exact p = ", format.pval(ft$p.value, digits = 3))

## Plot data
plot_data <- combined_meta %>%
  dplyr::filter(hilo_mrs == "Poor", !is.na(eurodepr)) %>%
  dplyr::group_by(sex, eurodepr) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::group_by(sex) %>%
  dplyr::mutate(pct = round(n / sum(n) * 100, 1))

## Plot
p_depr_poormrs <- ggplot2::ggplot(plot_data, ggplot2::aes(x = sex, y = pct, fill = eurodepr)) +
  ggplot2::geom_bar(stat = "identity", position = "stack", width = 0.5) +
  ggplot2::geom_text(
    ggplot2::aes(label = paste0(pct, "%\n(n=", n, ")")),
    position = ggplot2::position_stack(vjust = 0.5),
    size = 6, fontface = "bold", colour = "white"
  ) +
  ggplot2::scale_fill_manual(
    values = c("No" = "#456787", "Yes" = "#7E1811"),
    labels = c("No" = "No Depression/Anxiety", "Yes" = "Depression/Anxiety"),
    name   = "EQ-5D Depression/Anxiety"
  ) +
  ggplot2::scale_y_continuous(limits = c(0, 120), breaks = seq(0, 100, 20)) +
  ggplot2::labs(
    title = "Depression and Anxiety by Sex in Patients with Poor Functional Outcome (mRS 4-6)",
    x = "Sex", y = "Percentage (%)"
  ) +
  ## Prism-style bracket
  ggplot2::annotate("segment", x = 1, xend = 2, y = 110, yend = 110) +
  ggplot2::annotate("segment", x = 1, xend = 1, y = 106, yend = 110) +
  ggplot2::annotate("segment", x = 2, xend = 2, y = 106, yend = 110) +
  ggplot2::annotate("text", x = 1.5, y = 113, label = pval_txt,
                    hjust = 0.5, fontface = "bold", size = 5) +
  ggplot2::theme_classic(base_size = 18) +
  ggplot2::theme(
    plot.title      = ggplot2::element_text(face = "bold", size = 15),
    axis.text       = ggplot2::element_text(colour = "black"),
    legend.title    = ggplot2::element_text(face = "bold", size = 12),
    legend.position = "right"
  )
print(p_depr_poormrs)

## Save as PDF
# pdf("depression_by_sex_poor_outcome.pdf", width = 11, height = 7, useDingbats = FALSE)
# print(p_depr_poormrs)
# dev.off()

# To save:
# ggsave("depression_by_sex_poor_outcome.png", plot = p_depr_poormrs,
#        width = 11, height = 7, units = "in", dpi = 300)


################################################################################
# SECTION 15: PCA — SEX COHORT
#
# PCA is computed on the sex-contrast VST matrix (blind = FALSE).
# Two visualizations are produced:
#   1. Points colored by IFI27 expression (continuous, plasma palette)
#   2. Points colored by sex
#
# IFI27 is an interferon-stimulated gene of interest identified in the DEG
# analysis; its expression gradient across PC space is assessed here.
################################################################################

vsd_vis <- vst(dds_sex, blind = FALSE)

pca_data    <- plotPCA(vsd_vis, intgroup = "sex", returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))
pca_data$idrc <- rownames(pca_data)

# Append clinical metadata to PCA data frame
pca_data <- pca_data %>%
  dplyr::left_join(
    combined_meta %>%
      tibble::rownames_to_column("idrc") %>%
      dplyr::select(idrc, eurodepr, hilo_mrs, ich_hemisphere_adj,
                    ich_location_adj),
    by = "idrc"
  )

# Append IFI27 VST expression
ifi27_id       <- grep("IFI27$", rownames(vst_mat), value = TRUE)[1]
pca_data$IFI27 <- as.numeric(vst_mat[ifi27_id, pca_data$idrc])


# --- PCA colored by IFI27 expression -----------------------------------------
# Points are drawn in ascending IFI27 order so high-expressing samples
# appear on top.
p_pca_ifi27 <- pca_data %>%
  dplyr::arrange(IFI27) %>%
  ggplot(aes(PC1, PC2, color = IFI27, shape = sex)) +
  geom_point(
    position = position_jitter(width = 0.3, height = 0.3, seed = 42),
    size     = 4
  ) +
  scale_color_viridis_c(option = "plasma") +
  scale_shape_manual(values = c("Female" = 16, "Male" = 17)) +
  labs(
    x     = paste0("PC1: ", percent_var[1], "% variance"),
    y     = paste0("PC2: ", percent_var[2], "% variance"),
    title = "PCA — colored by IFI27 expression"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text  = element_text(size = 14),
    axis.title = element_text(size = 16)
  )

print(p_pca_ifi27)

# To save:
# ggsave("PCA_IFI27_expression.png", plot = p_pca_ifi27,
#        width = 10, height = 6, units = "in", dpi = 300)

# pdf("PCA_IFI27_expression.pdf", width = 10, height = 6, useDingbats = FALSE)
# print(p_pca_ifi27)
# dev.off()

# --- PCA colored by sex ------------------------------------------------------

p_pca_sex <- ggplot(pca_data, aes(PC1, PC2, color = sex)) +
  geom_point(size = 4) +
  scale_color_manual(values = c("Male" = "#364981", "Female" = "#FFC482")) +
  labs(
    x     = paste0("PC1: ", percent_var[1], "% variance"),
    y     = paste0("PC2: ", percent_var[2], "% variance"),
    title = "PCA — colored by sex"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text  = element_text(size = 14),
    axis.title = element_text(size = 16)
  )

print(p_pca_sex)


################################################################################
# SECTION 16: CORE ISG GENE HEATMAP — Sex × Outcome × Depression
#
# Visualizes z-scored VST expression of a curated set of immune and
# neuroinflammatory genes across all samples, annotated by sex, functional
# outcome (mRS Good/Poor), and depression status.
#
# Two layouts are produced:
#   A. Fully clustered (ward.D2) — reveals transcriptomic groupings
#   B. Manually ordered by Sex → Outcome → Depression — facilitates
#      direct group comparisons
#
# Expression values are capped at ±2 SD before plotting.
# Note: pheatmap objects must be saved with png() + dev.off(), not ggsave().
################################################################################

# --- Define curated gene panel -----------------------------------------------
# Genes were selected based on relevance to ICH neuroinflammation, interferon
# signaling, microglial activation (DAM), and ferroptosis pathways.

core_genes <- c(
  # NF-kB / innate immune signaling
  "CGAS", "NFKBIA", "STING1",
  # Interferon-stimulated gene (ISG) signature
  "IFITM3", "IFNAR2", "IFI27", "OASL", "IRF7", "STAT1", "ISG15", "CXCL10",
  # Inflammatory cytokines and chemokines
  "TNF", "IL6", "CCL2",
  # Lipid / secreted
  "LGALS3BP",
  # Disease-associated microglia (DAM) markers (blood-detectable)
  "SPP1", "ITGAX", "CD9", "CD63", "TYROBP", "CSF1", "CSF1R",
  # Cathepsins
  "CTSB", "CTSD", "CTSL", "CTSZ",
  # Iron homeostasis / ferritin
  "FTH1", "FTL",
  # Galectin
  "LGALS3",
  # Mitochondrial / ferroptosis
  "TOMM22", "GPX4", "FDX1", "AIFM1", "IMMT"
)


# --- Match core genes to VST matrix ------------------------------------------
# VST rownames are in "ENSG...|SYMBOL" format; extract symbol after the pipe.

vst_symbols   <- gsub(".*\\|", "", rownames(vst_mat))
genes_found   <- core_genes[core_genes %in% vst_symbols]
genes_missing <- setdiff(core_genes, vst_symbols)

cat("Core genes found:  ", length(genes_found),   "\n")
cat("Core genes missing:", length(genes_missing), "\n")
if (length(genes_missing) > 0) cat("Missing:", paste(genes_missing, collapse = ", "), "\n")

row_ids       <- rownames(vst_mat)[vst_symbols %in% genes_found]
mat           <- vst_mat[row_ids, ]
rownames(mat) <- gsub(".*\\|", "", rownames(mat))


# --- Z-score and cap expression values ---------------------------------------

mat_scaled <- t(scale(t(mat)))
mat_scaled <- pmin(pmax(mat_scaled, -2), 2)   # cap at ±2 SD


# --- Build column annotation -------------------------------------------------

anno_col <- combined_meta[colnames(mat), c("sex", "hilo_mrs", "eurodepr")] %>%
  dplyr::rename(Sex = sex, Outcome = hilo_mrs, Depression = eurodepr) %>%
  dplyr::mutate(across(everything(), as.character)) %>%
  as.data.frame()

rownames(anno_col)       <- colnames(mat)
anno_col[is.na(anno_col)] <- "Unknown"

anno_colors <- list(
  Sex        = c("Male"    = "#364981",  "Female"  = "#FFC482"),   # kept as requested
  Outcome    = c("Good"    = "#009E73",  "Poor"    = "#D55E00",    "Unknown" = "grey80"),
  Depression  = c("No"      = "#0072B2",  "Yes"     = "#E69F00",    "Unknown" = "grey80")
)


# --- Option A: Fully clustered heatmap ---------------------------------------

p_heatmap_clustered <- pheatmap::pheatmap(
  mat_scaled,
  annotation_col    = anno_col,
  annotation_colors = anno_colors,
  cluster_rows      = TRUE,
  cluster_cols      = TRUE,
  clustering_method = "ward.D2",
  color             = colorRampPalette(c("#456787", "white", "#7E1811"))(100),
  breaks            = seq(-2, 2, length.out = 101),
  show_colnames     = FALSE,
  show_rownames     = TRUE,
  fontsize_row      = 9,
  fontsize          = 11,
  border_color      = NA,
  main              = "Core ISG Genes — All Samples (Clustered)"
)

#To save:
# png("core_genes_heatmap_clustered.png", width = 18, height = 8, units = "in", res = 300)
# print(p_heatmap_clustered)
# dev.off()

# --- Option B: Manually ordered heatmap (Sex → Outcome → Depression) ---------

sample_order <- combined_meta %>%
  tibble::rownames_to_column("idrc") %>%
  dplyr::filter(idrc %in% colnames(mat)) %>%
  dplyr::arrange(
    factor(sex,      levels = c("Male",   "Female")),
    factor(hilo_mrs, levels = c("Good",   "Poor")),
    factor(eurodepr, levels = c("No",     "Yes"))
  ) %>%
  dplyr::pull(idrc)

# Insert column gaps at group boundaries (Male-Good | Male-Poor | Female-Good)
gap_positions <- cumsum(
  table(paste0(combined_meta[sample_order, "sex"], "_",
               combined_meta[sample_order, "hilo_mrs"]))
)[-4]

p_heatmap_split <- pheatmap::pheatmap(
  mat_scaled[, sample_order],
  annotation_col    = anno_col[sample_order, , drop = FALSE],
  annotation_colors = anno_colors,
  cluster_rows      = TRUE,
  cluster_cols      = FALSE,
  clustering_method = "ward.D2",
  gaps_col          = as.numeric(gap_positions),
  color             = colorRampPalette(c("#456787", "white", "#7E1811"))(100),
  breaks            = seq(-2, 2, length.out = 101),
  show_colnames     = FALSE,
  show_rownames     = TRUE,
  fontsize_row      = 12,
  fontsize          = 11,
  border_color      = NA,
  main              = "Core ISG Genes - Grouped by Sex, Outcome, and Depression"
)

# To save:
# png("core_genes_heatmap_grouped.png", width = 19.5, height = 9, units = "in", res = 300)
# print(p_heatmap_split)
# dev.off()

cairo_pdf("core_genes_heatmap_grouped.pdf", width = 19.5, height = 9)
print(p_heatmap_split)
dev.off()

################################################################################
# SECTION 17: MULTI-VARIABLE GSEA PIPELINE
#
# Fits a single DESeq2 model adjusting simultaneously for age, mRS,
# depression (eurodepr), and sex, then extracts three independent contrasts:
#
#   1. Sex       — Female vs Male      (adjusted for age, mRS, eurodepr)
#   2. mRS       — continuous effect   (adjusted for age, sex, eurodepr)
#   3. Depression — Yes vs No          (adjusted for age, sex, mRS)
#
# GSEA is run for each contrast against all gene set collections defined
# in Section 12. Results are summarized and saved as dot plots.
#
# Requires: combined_counts, combined_meta (with age_scaled, sex, mrs,
#           eurodepr columns), gene_set_collections, run_gsea()
################################################################################


# --- Step 1: Select complete-case samples ------------------------------------
# Retain only samples with non-missing values for all model covariates.

keep_full <- combined_meta %>%
  tibble::rownames_to_column("idrc") %>%
  dplyr::filter(
    !is.na(age_scaled),
    !is.na(sex),
    !is.na(mrs),
    !is.na(eurodepr)
  ) %>%
  dplyr::pull(idrc)

cat("Samples entering full model:", length(keep_full), "\n")


# --- Step 2: Prepare metadata ------------------------------------------------

meta_full <- combined_meta[keep_full, ] %>%
  dplyr::mutate(
    sex      = factor(sex,      levels = c("Male",  "Female")),
    eurodepr = factor(eurodepr, levels = c("No",    "Yes")),
    mrs      = as.numeric(mrs)
  )

# Sanity checks
cat("Model covariate summary:\n")
cat("  age_scaled range:", range(meta_full$age_scaled, na.rm = TRUE), "\n")
cat("  sex levels:      ", levels(meta_full$sex),                     "\n")
cat("  mrs range:       ", range(meta_full$mrs,        na.rm = TRUE), "\n")
cat("  eurodepr levels: ", levels(meta_full$eurodepr),                "\n")
cat("  eurodepr counts:\n")
print(table(meta_full$eurodepr))


# --- Step 3: Build and filter DESeq2 object ----------------------------------
# Minimum count filter: >= 10 counts in at least 5 samples.

dds_full_model <- DESeqDataSetFromMatrix(
  countData = combined_counts[, keep_full],
  colData   = meta_full,
  design    = ~ age_scaled + mrs + eurodepr + sex
)

dds_full_model <- dds_full_model[
  rowSums(counts(dds_full_model) >= 10) >= 5, ]

cat("Genes after filtering:", nrow(dds_full_model), "\n")


# --- Step 4: Run DESeq2 ------------------------------------------------------

dds_full_model <- DESeq(dds_full_model)

cat("Model coefficients:\n")
print(resultsNames(dds_full_model))


# --- Step 5: Extract results for each contrast -------------------------------

# Sex: Female vs Male (adjusted for age, mRS, depression)
res_sex_full <- results(
  dds_full_model,
  contrast = c("sex", "Female", "Male"),
  tidy     = TRUE
) %>% dplyr::rename(gene_id = row)

# mRS: continuous effect (adjusted for age, sex, depression)
res_mrs_full <- results(
  dds_full_model,
  name = "mrs",
  tidy = TRUE
) %>% dplyr::rename(gene_id = row)

# Depression: Yes vs No (adjusted for age, sex, mRS)
res_eurodepr <- results(
  dds_full_model,
  contrast = c("eurodepr", "Yes", "No"),
  tidy     = TRUE
) %>% dplyr::rename(gene_id = row)


# --- Step 6: Build ranked gene lists for GSEA --------------------------------
# Uses make_ranked() — expects a tidy DESeq2 result with a `stat` column
# and gene IDs resolvable to HGNC symbols.

ranked_sex_full <- build_ranked_list_tidy(res_sex_full)
ranked_mrs_full <- build_ranked_list_tidy(res_mrs_full)
ranked_eurodepr <- build_ranked_list_tidy(res_eurodepr)

cat("Ranked list lengths:\n")
cat("  Sex (full model):  ", length(ranked_sex_full), "\n")
cat("  mRS (full model):  ", length(ranked_mrs_full), "\n")
cat("  Depression:        ", length(ranked_eurodepr), "\n")

# --- Step 7: Run GSEA for each contrast and collection ----------------------
# run_gsea() is defined in Section 11. All gene set collections from
# Section 12 are tested for each contrast.

run_gsea_all_collections <- function(ranked_list, label) {
  cat("\nRunning GSEA —", label, "\n")
  result <- lapply(names(gene_set_collections), function(coll) {
    cat("  Collection:", coll, "\n")
    run_gsea(ranked_list, gene_set_collections[[coll]])
  })
  names(result) <- names(gene_set_collections)
  result
}

gsea_sex_full <- run_gsea_all_collections(ranked_sex_full, "Sex (full model)")
gsea_mrs_full <- run_gsea_all_collections(ranked_mrs_full, "mRS (full model)")
gsea_eurodepr <- run_gsea_all_collections(ranked_eurodepr, "Depression")


# --- Step 8: Print significant Hallmark pathways per contrast ----------------

print_sig_pathways <- function(gsea_list, label) {
  cat("\n===", label, "— Hallmark significant pathways ===\n")
  gsea_list[["Hallmark"]]@result %>%
    dplyr::filter(p.adjust < 0.05) %>%
    dplyr::arrange(NES) %>%
    dplyr::select(Description, NES, pvalue, p.adjust) %>%
    print()
}

print_sig_pathways(gsea_sex_full, "SEX (full model)")
print_sig_pathways(gsea_mrs_full, "mRS (full model)")
print_sig_pathways(gsea_eurodepr, "DEPRESSION")


# --- Step 9: GSEA dot plot function ------------------------------------------
# Plots significantly enriched pathways (adjusted p < 0.05) as a dot plot
# with NES on the x-axis and gene set size mapped to point size.
# Returns NULL with a message if no pathways reach significance.
#
# Args:
#   gsea_obj  — gseaResult object
#   title     — plot title string
#   pos_label — legend label for positively enriched pathways
#   neg_label — legend label for negatively enriched pathways
#
# Returns:
#   ggplot object or NULL

plot_gsea_dotplot <- function(gsea_obj,
                              title,
                              pos_label = "↑ Up",
                              neg_label = "↓ Down") {
  
  sig <- gsea_obj@result %>%
    dplyr::filter(!is.na(p.adjust), p.adjust < 0.05) %>%
    dplyr::mutate(
      Description = stringr::str_remove(Description, "^HALLMARK_|^KEGG_"),
      Description = stringr::str_replace_all(Description, "_", " "),
      Description = stringr::str_to_title(Description),
      Direction   = ifelse(NES > 0, pos_label, neg_label)
    ) %>%
    dplyr::arrange(NES)
  
  if (nrow(sig) == 0) {
    message("No significant pathways for: ", title)
    return(NULL)
  }
  
  sig$Description <- factor(sig$Description, levels = sig$Description)
  
  ggplot2::ggplot(sig, ggplot2::aes(
    x     = NES,
    y     = Description,
    size  = setSize,
    color = Direction
  )) +
    ggplot2::geom_point() +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::scale_color_manual(values = c(
      setNames("#7E1811", pos_label),
      setNames("#456787", neg_label)
    )) +
    ggplot2::scale_size_continuous(range = c(3, 8)) +
    ggplot2::labs(
      title = title,
      x     = "Normalised Enrichment Score (NES)",
      y     = NULL,
      size  = "Gene set size",
      color = "Direction"
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(face = "bold", hjust = 0.5),
      axis.text.y      = ggplot2::element_text(size = 11, color = "black"),
      panel.grid.minor = ggplot2::element_blank(),
      legend.text      = ggplot2::element_text(size = 11, color = "black")
    )
}


# --- Step 10: Define contrasts for plotting ----------------------------------
# Each entry specifies the GSEA result list, plot title, file prefix,
# and directional labels that will appear in the plot legend.

contrasts_plot <- list(
  
  list(
    gsea      = gsea_sex_full,
    label     = "Sex — Female vs Male\n(age + mRS + depression adjusted)",
    prefix    = "sex_full",
    pos_label = "↑ Higher in Female",
    neg_label = "↓ Higher in Male"
  ),
  
  list(
    gsea      = gsea_mrs_full,
    label     = "mRS — continuous\n(age + sex + depression adjusted)",
    prefix    = "mrs_full",
    pos_label = "↑ Higher mRS",
    neg_label = "↓ Lower mRS"
  ),
  
  list(
    gsea      = gsea_eurodepr,
    label     = "Depression — Yes vs No\n(age + sex + mRS adjusted)",
    prefix    = "eurodepr",
    pos_label = "↑ Higher in Depressed",
    neg_label = "↓ Higher in Non-depressed"
  )
)


# --- Step 11: Generate plots (print only) ------------------------------------

for (contrast in contrasts_plot) {
  for (coll in names(gene_set_collections)) {
    
    p <- plot_gsea_dotplot(
      gsea_obj  = contrast$gsea[[coll]],
      title     = paste0(contrast$label, " — ", coll),
      pos_label = contrast$pos_label,
      neg_label = contrast$neg_label
    )
    
    if (!is.null(p)) print(p)
  }
}


# --- Step 12: Save each plot individually ------------------------------------

# ## sex_full
# p_sex_full_Hallmark <- plot_gsea_dotplot(
#   gsea_obj  = gsea_sex_full[["Hallmark"]],
#   title     = "Sex — Female vs Male\n(age + mRS + depression adjusted) — Hallmark",
#   pos_label = "↑ Higher in Female",
#   neg_label = "↓ Higher in Male"
# )
# cairo_pdf("gsea_dotplot_sex_full_Hallmark.pdf", width = 9, height = 6)
# print(p_sex_full_Hallmark)
# dev.off()
# 
# p_sex_full_KEGG <- plot_gsea_dotplot(
#   gsea_obj  = gsea_sex_full[["KEGG"]],
#   title     = "Sex — Female vs Male\n(age + mRS + depression adjusted) — KEGG",
#   pos_label = "↑ Higher in Female",
#   neg_label = "↓ Higher in Male"
# )
# cairo_pdf("gsea_dotplot_sex_full_KEGG.pdf", width = 9, height = 12)
# print(p_sex_full_KEGG)
# dev.off()
# 
# ## mrs_full
# p_mrs_full_Hallmark <- plot_gsea_dotplot(
#   gsea_obj  = gsea_mrs_full[["Hallmark"]],
#   title     = "mRS — continuous\n(age + sex + depression adjusted) — Hallmark",
#   pos_label = "↑ Higher mRS",
#   neg_label = "↓ Lower mRS"
# )
# cairo_pdf("gsea_dotplot_mrs_full_Hallmark.pdf", width = 8, height = 7.5)
# print(p_mrs_full_Hallmark)
# dev.off()
# 
# p_mrs_full_KEGG <- plot_gsea_dotplot(
#   gsea_obj  = gsea_mrs_full[["KEGG"]],
#   title     = "mRS — continuous\n(age + sex + depression adjusted) — KEGG",
#   pos_label = "↑ Higher mRS",
#   neg_label = "↓ Lower mRS"
# )
# cairo_pdf("gsea_dotplot_mrs_full_KEGG.pdf", width = 9, height = 6)
# print(p_mrs_full_KEGG)
# dev.off()
# 
# ## eurodepr
# p_eurodepr_Hallmark <- plot_gsea_dotplot(
#   gsea_obj  = gsea_eurodepr[["Hallmark"]],
#   title     = "Depression — Yes vs No\n(age + sex + mRS adjusted) — Hallmark",
#   pos_label = "↑ Higher in Depressed",
#   neg_label = "↓ Higher in Non-depressed"
# )
# cairo_pdf("gsea_dotplot_eurodepr_Hallmark.pdf", width = 9, height = 6.8)
# print(p_eurodepr_Hallmark)
# dev.off()
# 
# p_eurodepr_KEGG <- plot_gsea_dotplot(
#   gsea_obj  = gsea_eurodepr[["KEGG"]],
#   title     = "Depression — Yes vs No\n(age + sex + mRS adjusted) — KEGG",
#   pos_label = "↑ Higher in Depressed",
#   neg_label = "↓ Higher in Non-depressed"
# )
# cairo_pdf("gsea_dotplot_eurodepr_KEGG.pdf", width = 9, height = 6)
# print(p_eurodepr_KEGG)
# dev.off()

################################################################################
# SECTION 18: GSEA ENRICHMENT PLOTS — TOP PATHWAYS PER CONTRAST
#
# For each contrast and gene set collection, generates enrichment running-score
# plots (via enrichplot::gseaplot2) for the top 3 most significant pathways
# (ranked by adjusted p-value). Plots are printed and saved automatically.
#
# Requires: contrasts_plot (from Section 17), gene_set_collections
################################################################################

for (contrast in contrasts_plot) {
  for (coll in names(gene_set_collections)) {
    
    # Select top 3 significant pathways by adjusted p-value
    top_paths <- contrast$gsea[[coll]]@result %>%
      dplyr::filter(!is.na(p.adjust), p.adjust < 0.05) %>%
      dplyr::arrange(p.adjust) %>%
      dplyr::slice_head(n = 3) %>%
      dplyr::pull(ID)
    
    for (pw in top_paths) {
      
      p <- enrichplot::gseaplot2(
        contrast$gsea[[coll]],
        geneSetID = pw,
        title     = paste0(contrast$label, "\n", pw)
      )
      
      print(p)
      
      # Sanitize pathway name for use in filename
      clean_pw <- pw %>%
        stringr::str_remove("^HALLMARK_|^KEGG_") %>%
        stringr::str_replace_all("[^A-Za-z0-9]", "_")
      
      filename <- paste0("gsea_enrichment_", contrast$prefix,
                         "_", coll, "_", clean_pw, ".png")
      
      ggplot2::ggsave(filename = filename, plot = p,
                      width = 8, height = 6, dpi = 300)
      cat("Saved:", filename, "\n")
    }
  }
}


################################################################################
# SECTION 19: IFI27 EXPRESSION — SEX DIFFERENCES
#
# Extracts VST-normalized IFI27 expression for all samples and tests for
# sex differences using a linear model adjusted for age, mRS, and depression.
# Results are displayed as a violin + sina plot.
#
# IFI27 is an interferon-stimulated gene elevated in innate immune activation;
# sex differences in its expression may reflect divergent immune profiles
# following ICH.
#
# Requires: vst_mat, combined_counts, mini_data
################################################################################

library(ggforce)

# --- Build gene ID lookup for IFI27 ------------------------------------------
# VST rownames are in "ENSG...|SYMBOL" format; grep for exact symbol match.

gene_ids <- c(IFI27 = grep("IFI27$", rownames(vst_mat), value = TRUE)[1])


# --- Build expression metadata frame -----------------------------------------
# Aligns sample-level expression with clinical metadata from mini_data.

expr_meta <- data.frame(idrc = colnames(combined_counts)) %>%
  dplyr::left_join(mini_data, by = "idrc")

rownames(expr_meta) <- expr_meta$idrc

# Verify alignment
cat("Required columns present in expr_meta:\n")
cat("  idrc:     ", "idrc"     %in% colnames(expr_meta), "\n")
cat("  eurodepr: ", "eurodepr" %in% colnames(expr_meta), "\n")
cat("  mrs:      ", "mrs"      %in% colnames(expr_meta), "\n")

unmatched <- setdiff(colnames(vst_mat), rownames(expr_meta))
if (length(unmatched) > 0) {
  cat("Samples in VST matrix not found in expr_meta:", unmatched, "\n")
} else {
  cat("All VST samples matched in expr_meta.\n")
}


# --- Assemble per-gene expression data frame ---------------------------------

expr_df_master <- lapply(names(gene_ids), function(g) {
  row_id <- gene_ids[[g]]
  data.frame(
    idrc               = colnames(vst_mat),
    gene               = g,
    expr               = as.numeric(vst_mat[row_id, ]),
    sex                = expr_meta[colnames(vst_mat), "sex"],
    age                = expr_meta[colnames(vst_mat), "age"],
    eurodepr           = expr_meta[colnames(vst_mat), "eurodepr"],
    hilo_mrs           = expr_meta[colnames(vst_mat), "hilo_mrs"],
    mrs                = expr_meta[colnames(vst_mat), "mrs"],
    tics_total         = expr_meta[colnames(vst_mat), "tics_total"],
    ich_location_adj   = expr_meta[colnames(vst_mat), "ich_location_adj"],
    ich_hemisphere_adj = expr_meta[colnames(vst_mat), "ich_hemisphere_adj"],
    sbp                = expr_meta[colnames(vst_mat), "sbp"]
  )
}) %>%
  dplyr::bind_rows()

cat("Expression data frame:", nrow(expr_df_master), "rows,",
    ncol(expr_df_master), "columns\n")


# --- Clean and factor variables ----------------------------------------------

cat("eurodepr distribution (pre-factoring):\n")
print(table(expr_df_master$eurodepr, useNA = "always"))

expr_df_clean <- expr_df_master %>%
  tidyr::drop_na(expr, sex, age) %>%
  dplyr::mutate(
    sex      = factor(sex,      levels = c("Female", "Male")),
    age      = as.numeric(age),
    eurodepr = factor(eurodepr, levels = c("No",     "Yes")),
    hilo_mrs = factor(hilo_mrs, levels = c("Good",   "Poor")),
    mrs      = as.numeric(mrs)
  )

cat("eurodepr levels (reference = first):", levels(expr_df_clean$eurodepr), "\n")


# --- Adjusted sex difference test per gene -----------------------------------
# Linear model: expr ~ sex + age + mRS + depression.
# The sex coefficient p-value is extracted and converted to a significance label.

sex_age_tests <- expr_df_clean %>%
  dplyr::group_by(gene) %>%
  dplyr::group_modify(~ {
    fit     <- lm(expr ~ sex + age + mrs + eurodepr, data = .x)
    out     <- summary(fit)$coefficients
    sex_row <- grep("^sex", rownames(out), value = TRUE)
    data.frame(p_sex = out[sex_row, "Pr(>|t|)"])
  }) %>%
  dplyr::mutate(
    signif = dplyr::case_when(
      p_sex < 0.001 ~ "***",
      p_sex < 0.01  ~ "**",
      p_sex < 0.05  ~ "*",
      TRUE          ~ "ns"
    )
  )

# Append test results and compute y-position for significance label
expr_df_plot <- expr_df_clean %>%
  dplyr::left_join(sex_age_tests, by = "gene") %>%
  dplyr::group_by(gene) %>%
  dplyr::mutate(y_pos = max(expr) * 0.95) %>%
  dplyr::ungroup()


# --- Violin + sina plot ------------------------------------------------------

p_ifi27_sex <- ggplot(expr_df_plot, aes(x = sex, y = expr, fill = sex)) +
  geom_violin(trim = FALSE, alpha = 0.6, scale = "width") +
  ggforce::geom_sina(size = 1.8, alpha = 0.8, colour = "black", seed = 123) +
  geom_text(aes(x = 1.5, y = y_pos, label = signif),
            color = "black", size = 10) +
  scale_fill_manual(values = c("Female" = "#FFC482", "Male" = "#1B2F6E")) +
  labs(
    title    = "Sex Differences in IFI27 Expression",
    subtitle = "Adjusted for Age, mRS, and Depression Status",
    x        = "Sex",
    y        = "Expression (VST)"
  ) +
  theme_bw() +
  theme(
    panel.grid      = element_blank(),
    strip.text      = element_text(size = 16, color = "black"),
    axis.text       = element_text(size = 14, color = "black"),
    axis.title      = element_text(size = 16, color = "black"),
    plot.title      = element_text(size = 18, color = "black"),
    legend.position = "none"
  )

print(p_ifi27_sex)

# To save:
# ggsave("IFI27_sex_differences.png", plot = p_ifi27_sex,
#        width = 6.9, height = 5.6, units = "in", dpi = 300)

# pdf("IFI27_sex_differences.pdf", width = 6.9, height = 5.6, useDingbats = FALSE)
# print(p_ifi27_sex)
# dev.off()

################################################################################
# SECTION 20: SEROTONIN GENE ANALYSIS — SEX × DEPRESSION
#
# Tests whether serotonin pathway gene expression differs by depression status,
# and whether that relationship is moderated by sex.
#
# Part A: Sex × Depression interaction test (LRT) for each gene
# Part B: Overall Wilcoxon test, depressed vs non-depressed (all patients)
#
# Effect size: rank-biserial correlation (r)
# Multiple testing correction: Benjamini-Hochberg (BH)
#
# Requires: vst_mat_blind, combined_meta, clean_gene_names()
################################################################################

# --- Step 1: Define candidate serotonin gene panel ---------------------------

serotonin_set <- c("TPH1", "SLC6A4", "MAOA", "MAOB",
                   "HTR1B", "HTR2B", "HTR3A", "HTR6", "HTR7")


# --- Step 2: Identify genes present in the VST matrix -----------------------

cleaned_mat  <- clean_gene_names(vst_mat_blind)
found_sero   <- serotonin_set[serotonin_set %in% rownames(cleaned_mat)]
missing_sero <- setdiff(serotonin_set, rownames(cleaned_mat))

cat("Serotonin genes found:  ", length(found_sero),   "\n")
cat("Serotonin genes missing:", length(missing_sero), "\n")
if (length(missing_sero) > 0) cat("Missing:", paste(missing_sero, collapse = ", "), "\n")

serotonin_set <- found_sero   # restrict to available genes


# --- Step 3: Build analysis data frame ---------------------------------------
# Transpose the expression matrix and join with combined_meta.

sero_expr <- cleaned_mat[serotonin_set, ] %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column("idrc")

sero_meta <- combined_meta %>%
  tibble::rownames_to_column("idrc") %>%
  dplyr::left_join(sero_expr, by = "idrc")

cat("sero_meta dimensions:", dim(sero_meta), "\n")
cat("Serotonin genes in sero_meta:",
    paste(serotonin_set[serotonin_set %in% colnames(sero_meta)],
          collapse = ", "), "\n")


# --- Part A: Sex × Depression interaction test -------------------------------
# For each gene: fits interaction model (eurodepr * sex) and main-effects model
# (eurodepr + sex), then uses a likelihood ratio test (LRT via anova()) to
# assess whether the interaction term improves model fit.
# Marginal group means (Female/Male × Depressed/Not) are also reported.

cat("\n=== Part A: Sex × Depression Interaction ===\n")

interaction_results <- purrr::map_dfr(serotonin_set, function(gene) {
  
  df <- sero_meta %>%
    dplyr::filter(!is.na(eurodepr), !is.na(sex), !is.na(.data[[gene]])) %>%
    dplyr::mutate(
      eurodepr = factor(eurodepr, levels = c("No",     "Yes")),
      sex      = factor(sex,      levels = c("Female", "Male"))
    )
  
  mod_full <- lm(as.formula(paste(gene, "~ eurodepr * sex")), data = df)
  mod_main <- lm(as.formula(paste(gene, "~ eurodepr + sex")), data = df)
  
  coefs   <- summary(mod_full)$coefficients
  int_row  <- grep("eurodepr.*sex|sex.*eurodepr", rownames(coefs), value = TRUE)
  depr_row <- grep("^eurodeprYes", rownames(coefs), value = TRUE)
  
  # Print coefficient names for the first gene as a diagnostic check
  if (gene == serotonin_set[1]) {
    cat("Coefficient rownames (", gene, "):", rownames(coefs), "\n")
  }
  
  lrt_p  <- anova(mod_main, mod_full)$`Pr(>F)`[2]
  int_est <- if (length(int_row)  > 0) coefs[int_row[1],  "Estimate"]   else NA
  int_p   <- if (length(int_row)  > 0) coefs[int_row[1],  "Pr(>|t|)"]  else NA
  
  # Marginal group means
  safe <- function(x) if (length(x) == 0) NA_real_ else x
  
  mm <- df %>%
    dplyr::group_by(sex, eurodepr) %>%
    dplyr::summarise(mean = mean(.data[[gene]], na.rm = TRUE), .groups = "drop")
  
  female_no  <- safe(mm$mean[mm$sex == "Female" & mm$eurodepr == "No"])
  female_yes <- safe(mm$mean[mm$sex == "Female" & mm$eurodepr == "Yes"])
  male_no    <- safe(mm$mean[mm$sex == "Male"   & mm$eurodepr == "No"])
  male_yes   <- safe(mm$mean[mm$sex == "Male"   & mm$eurodepr == "Yes"])
  
  data.frame(
    gene          = gene,
    female_no     = round(female_no,              3),
    female_yes    = round(female_yes,             3),
    male_no       = round(male_no,                3),
    male_yes      = round(male_yes,               3),
    diff_female   = round(female_yes - female_no, 3),
    diff_male     = round(male_yes   - male_no,   3),
    interaction_b = round(int_est,                4),
    interaction_p = round(int_p,                  4),
    lrt_p         = round(lrt_p,                  4),
    r2_full       = round(summary(mod_full)$r.squared, 3)
  )
}) %>%
  dplyr::mutate(
    lrt_p_adj         = round(p.adjust(lrt_p,         method = "BH"), 4),
    interaction_p_adj = round(p.adjust(interaction_p, method = "BH"), 4),
    signif            = dplyr::case_when(
      lrt_p_adj < 0.05 ~ "*",
      lrt_p     < 0.05 ~ "nom",
      lrt_p     < 0.10 ~ "†",
      TRUE             ~ "ns"
    )
  ) %>%
  dplyr::arrange(lrt_p)

cat("Interaction results (sorted by LRT p-value):\n")
print(interaction_results)


# --- Part B: Overall Wilcoxon test (all patients) ----------------------------
# Tests whether expression differs between depressed and non-depressed patients
# across the full sample, irrespective of sex.

cat("\n=== Part B: Overall Depressed vs Non-depressed ===\n")

overall_results <- purrr::map_dfr(serotonin_set, function(gene) {
  
  df <- sero_meta %>%
    dplyr::filter(!is.na(eurodepr), !is.na(.data[[gene]]))
  
  wt      <- wilcox.test(as.formula(paste(gene, "~ eurodepr")),
                         data = df, exact = FALSE)
  r_val <- -(effectsize::rank_biserial(df[[gene]],
                                       df$eurodepr)$r_rank_biserial)
  med_no  <- median(df[[gene]][df$eurodepr == "No"],  na.rm = TRUE)
  med_yes <- median(df[[gene]][df$eurodepr == "Yes"], na.rm = TRUE)
  
  data.frame(
    gene     = gene,
    n_no     = sum(df$eurodepr == "No"),
    n_yes    = sum(df$eurodepr == "Yes"),
    med_no   = round(med_no,           3),
    med_yes  = round(med_yes,          3),
    med_diff = round(med_yes - med_no, 3),
    effect_r = round(r_val,            3),
    W        = wt$statistic,
    p.value  = round(wt$p.value,       4)
  )
}) %>%
  dplyr::mutate(
    p_adj     = round(p.adjust(p.value, method = "BH"), 4),
    direction = ifelse(med_diff > 0, "↑ depressed", "↓ depressed"),
    signif    = dplyr::case_when(
      p_adj   < 0.05 ~ "*",
      p.value < 0.05 ~ "nom",
      p.value < 0.10 ~ "†",
      TRUE           ~ "ns"
    )
  ) %>%
  dplyr::arrange(p.value)

cat("Overall results (sorted by p-value):\n")
print(overall_results)


################################################################################
# SECTION 21: SEROTONIN GENE VISUALIZATIONS
#
# Visualization A: Interaction plot — mean expression per Sex × Depression group
#                  Faceted by gene; crossing lines indicate a possible interaction
#
# Visualization B: Volcano-style dot plot — rank-biserial effect size vs
#                  -log10(p) for the overall depressed vs non-depressed test
################################################################################

# --- Visualization A: Sex × Depression interaction plot ----------------------

int_plot_df <- interaction_results %>%
  tidyr::pivot_longer(
    cols      = c(female_no, female_yes, male_no, male_yes),
    names_to  = "group",
    values_to = "mean_expr"
  ) %>%
  dplyr::mutate(
    sex      = ifelse(stringr::str_starts(group, "female"), "Female", "Male"),
    eurodepr = ifelse(stringr::str_ends(group,  "yes"),    "Yes",    "No"),
    gene     = factor(gene, levels = interaction_results$gene),
    label    = dplyr::case_when(
      lrt_p < 0.05 ~ paste0(gene, "*"),
      lrt_p < 0.10 ~ paste0(gene, "†"),
      TRUE         ~ as.character(gene)
    )
  )

p_interaction <- ggplot(
  int_plot_df,
  aes(x = eurodepr, y = mean_expr, color = sex, group = sex)
) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  facet_wrap(~ gene, scales = "free_y", ncol = 3) +
  scale_color_manual(values = c("Female" = "#FFC482", "Male" = "#1B2F6E")) +
  labs(
    title    = "Sex × Depression Interaction — Serotonin Genes",
    subtitle = "* LRT p < 0.05 | † LRT p < 0.10 | Crossing lines suggest interaction",
    x        = "Depression status",
    y        = "Mean VST expression",
    color    = "Sex"
  ) +
  theme_bw() +
  theme(
    panel.grid       = element_blank(),
    strip.background = element_rect(fill = "grey90"),
    strip.text       = element_text(size = 10, face = "bold"),
    axis.text        = element_text(size = 10),
    legend.position  = "bottom"
  )

print(p_interaction)

# To save:
# ggsave("serotonin_sex_x_depression_interaction.png", plot = p_interaction,
#        width = 10, height = 10, dpi = 300)

# pdf("serotonin_sex_x_depression_interaction.pdf", width = 10, height = 10, useDingbats = FALSE)
# print(p_interaction)
# dev.off()

# --- Visualization B: Volcano-style effect size plot (overall test) ----------

p_overall_sero <- ggplot(
  overall_results,
  aes(x = effect_r, y = -log10(p.value), color = signif, label = gene)
) +
  geom_point(size = 4) +
  ggrepel::geom_text_repel(size = 3.5, max.overlaps = 20) +
  geom_vline(xintercept = 0,           linetype = "dashed", colour = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", colour = "red",
             alpha = 0.6) +
  scale_color_manual(values = c(
    "*"   = "#c0392b",
    "nom" = "#e67e22",
    "†"   = "#f1c40f",
    "ns"  = "grey60"
  )) +
  labs(
    title    = "Serotonin Genes: Depressed vs Non-depressed (all patients)",
    subtitle = "Rank-biserial effect size | positive = higher in depressed",
    x        = "Effect size r (rank-biserial)",
    y        = "-log10(p-value)",
    color    = "Significance"
  ) +
  theme_bw() +
  theme(
    panel.grid      = element_blank(),
    legend.position = "right"
  )

print(p_overall_sero)

# To save:
# ggsave("serotonin_overall_depression_volcano.png", plot = p_overall_sero,
#        width = 8, height = 6, dpi = 300)

################################################################################
# SECTION 22: SEX-STRATIFIED SEROTONIN GENE ANALYSIS
#
# Runs Wilcoxon rank-sum tests separately within Female and Male patients to
# assess whether depression-associated expression differences in serotonin
# genes are sex-specific. This tests directionality and magnitude of the
# depression effect within each sex.
#
# Multiple testing correction is applied within each sex group (BH method).
# Genes with fewer than 3 observations per depression group are skipped.
#
# Requires: sero_meta, serotonin_set, rank_biserial()
################################################################################

# --- Helper: run sex-stratified Wilcoxon for one gene × sex combination ------

run_stratified_wilcoxon <- function(gene, sx, data) {
  
  df <- data %>%
    dplyr::filter(sex == sx, !is.na(eurodepr), !is.na(.data[[gene]]))
  
  n_no  <- sum(df$eurodepr == "No",  na.rm = TRUE)
  n_yes <- sum(df$eurodepr == "Yes", na.rm = TRUE)
  
  if (n_no < 3 || n_yes < 3) {
    return(data.frame(
      sex      = sx, gene     = gene,
      n_no     = n_no,  n_yes    = n_yes,
      med_no   = NA,    med_yes  = NA,
      med_diff = NA,    effect_r = NA,
      W        = NA,    p.value  = NA
    ))
  }
  
  wt      <- wilcox.test(as.formula(paste(gene, "~ eurodepr")),
                         data = df, exact = FALSE)
  r_val   <- -(effectsize::rank_biserial(df[[gene]],
                                         df$eurodepr)$r_rank_biserial)
  med_no  <- median(df[[gene]][df$eurodepr == "No"],  na.rm = TRUE)
  med_yes <- median(df[[gene]][df$eurodepr == "Yes"], na.rm = TRUE)
  
  data.frame(
    sex      = sx,
    gene     = gene,
    n_no     = n_no,
    n_yes    = n_yes,
    med_no   = round(med_no,           3),
    med_yes  = round(med_yes,          3),
    med_diff = round(med_yes - med_no, 3),
    effect_r = round(r_val,            3),
    W        = wt$statistic,
    p.value  = round(wt$p.value,       4)
  )
}


# --- Run stratified analysis -------------------------------------------------

stratified_results <- purrr::map_dfr(
  c("Female", "Male"),
  function(sx) purrr::map_dfr(serotonin_set, run_stratified_wilcoxon,
                              sx = sx, data = sero_meta)
) %>%
  dplyr::group_by(sex) %>%
  dplyr::mutate(
    p_adj     = round(p.adjust(p.value, method = "BH"), 4),
    direction = dplyr::case_when(
      med_diff >  0 ~ "↑ in depressed",
      med_diff <  0 ~ "↓ in depressed",
      med_diff == 0 ~ "no change"
    ),
    signif = dplyr::case_when(
      p_adj   < 0.05 ~ "*",
      p.value < 0.05 ~ "nom",
      p.value < 0.10 ~ "†",
      TRUE           ~ "ns"
    )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(sex, p.value)

cat("=== Female patients ===\n")
print(dplyr::filter(stratified_results, sex == "Female"))

cat("\n=== Male patients ===\n")
print(dplyr::filter(stratified_results, sex == "Male"))


################################################################################
# SECTION 23: SEX-STRATIFIED SEROTONIN VISUALIZATIONS
#
# Visualization C: Dot plot — effect size by gene, colored by sex,
#                  shape encodes significance. Useful for comparing direction
#                  and magnitude across sexes side by side.
#
# Visualization D: Diverging bar chart — same data as C but displayed as bars;
#                  left = downregulated in depressed, right = upregulated.
#                  Bar alpha encodes significance level.
################################################################################

# --- Visualization C: Dot plot (effect size × gene, by sex) ------------------

p_strat_dot <- ggplot(
  dplyr::filter(stratified_results, !is.na(effect_r)),
  aes(x = effect_r, y = reorder(gene, effect_r),
      color = sex, shape = signif)
) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(size = 4, position = position_dodge(width = 0.5)) +
  scale_color_manual(values = c("Female" = "#FFC482", "Male" = "#1B2F6E")) +
  scale_shape_manual(values = c("*" = 8, "nom" = 17, "†" = 2, "ns" = 16)) +
  scale_x_continuous(limits = c(-0.75, 0.80),
                     breaks = seq(-0.6, 0.6, by = 0.2),
                     labels = scales::label_number())+
  labs(
    title    = "Depression Effect on Serotonin Genes — by Sex",
    subtitle = "Positive = higher in depressed | * p_adj < 0.05 | nom = p < 0.05 uncorrected",
    x        = "Effect size r (rank-biserial)",
    y        = "Gene",
    color    = "Sex",
    shape    = "Significance"
  ) +
  theme_bw() +
  theme(
    panel.grid.major.y = element_line(colour = "grey90"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position    = "right",
    axis.text          = element_text(size = 11)
  )

print(p_strat_dot)

# To save:
# ggsave("serotonin_stratified_dotplot.png", plot = p_strat_dot,
#        width = 8, height = 5, dpi = 300)


# --- Visualization D: Diverging bar chart (effect size, alpha = significance)-

p_strat_bar <- stratified_results %>%
  dplyr::filter(!is.na(effect_r)) %>%
  dplyr::mutate(gene = factor(gene, levels = serotonin_set)) %>%
  ggplot(aes(x = effect_r, y = gene, fill = sex, alpha = signif)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_vline(xintercept =  0.0, colour = "black",  linewidth = 0.6) +
  geom_vline(xintercept =  0.3, colour = "grey60", linewidth = 0.4,
             linetype = "dotted") +
  geom_vline(xintercept = -0.3, colour = "grey60", linewidth = 0.4,
             linetype = "dotted") +
  # Significance annotations in the right margin
  geom_text(
    aes(
      x     = 0.72,
      label = dplyr::case_when(
        p.value < 0.05 ~ "*",
        p.value < 0.10 ~ "†",
        TRUE           ~ ""
      ),
      color = sex
    ),
    position   = position_dodge(width = 0.7),
    size       = 5,
    fontface   = "bold",
    show.legend = FALSE
  ) +
  scale_fill_manual(values  = c("Female" = "#FFC482", "Male" = "#1B2F6E")) +
  scale_color_manual(values = c("Female" = "#FFC482", "Male" = "#1B2F6E")) +
  scale_alpha_manual(values = c("*" = 1.00, "nom" = 0.75,
                                "†" = 0.75, "ns"  = 0.40),
                     guide  = "none") +
  scale_x_continuous(limits = c(-0.75, 0.80),
                     breaks = seq(-0.6, 0.6, by = 0.2),
                     labels = scales::label_number()) +
  labs(
    title    = "Serotonin Gene Expression: Depressed vs Non-Depressed by Sex",
    subtitle = paste("Rank-biserial effect size | positive = higher in depressed",
                     "| * p < 0.05 | † p < 0.10 | faded = ns"),
    x        = "Effect size r (rank-biserial)",
    y        = NULL,
    fill     = "Sex"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major.y = element_line(colour = "grey93"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position    = "bottom",
    axis.text.y        = element_text(size = 11, face = "bold"),
    axis.text.x        = element_text(size = 10),
    plot.title         = element_text(face = "bold"),
    plot.subtitle      = element_text(size = 9, colour = "grey40")
  )

print(p_strat_bar)

# To save:
# ggsave("serotonin_stratified_diverging_bar.png", plot = p_strat_bar,
#        width = 9, height = 5, dpi = 300)


################################################################################
# SECTION 24: HTR6 EXPRESSION — DEPRESSED VS NON-DEPRESSED (FULL SAMPLE)
#
# Focused single-gene analysis of HTR6, which encodes the serotonin receptor 6.
# Tests expression differences between depressed and non-depressed patients
# across the full sample using a Wilcoxon rank-sum test.
# Effect size is reported as rank-biserial correlation (r).
################################################################################

htr6_df <- sero_meta %>%
  dplyr::filter(!is.na(eurodepr), !is.na(HTR6))

cat("HTR6 sample counts:\n")
cat("  Total:          ", nrow(htr6_df), "\n")
cat("  Non-depressed:  ", sum(htr6_df$eurodepr == "No"),  "\n")
cat("  Depressed:      ", sum(htr6_df$eurodepr == "Yes"), "\n")

# Wilcoxon test and effect size
htr6_wt    <- wilcox.test(HTR6 ~ eurodepr, data = htr6_df, exact = FALSE)
htr6_r <- effectsize::rank_biserial(htr6_df$HTR6, htr6_df$eurodepr)$r_rank_biserial
cat("effect r =", round(htr6_r, 3), "\n")

# Descriptive statistics by group
htr6_df %>%
  dplyr::group_by(eurodepr) %>%
  dplyr::summarise(
    n      = dplyr::n(),
    median = round(median(HTR6, na.rm = TRUE), 3),
    mean   = round(mean(HTR6,   na.rm = TRUE), 3),
    sd     = round(sd(HTR6,     na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  print()

cat("\nWilcoxon W =", htr6_wt$statistic,          "\n")
cat("p-value    =", round(htr6_wt$p.value, 4),    "\n")
cat("effect r   =", round(htr6_r, 3),              "\n")

# Boxplot + jitter
p_htr6 <- ggplot(htr6_df, aes(x = eurodepr, y = HTR6, fill = eurodepr)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.5) +
  geom_jitter(aes(colour = eurodepr), width = 0.12, size = 2.5, alpha = 0.8) +
  ggpubr::stat_compare_means(
    method  = "wilcox.test",
    label   = "p.format",
    size    = 4,
    label.x = 1.4
  ) +
  scale_fill_manual(values   = c("No" = "#B0C4DE", "Yes" = "#E07B7B")) +
  scale_colour_manual(values = c("No" = "#456787", "Yes" = "#7E1811")) +
  labs(
    title    = "HTR6 Expression: Depressed vs Non-Depressed",
    subtitle = paste0("Full sample | Wilcoxon rank-sum | n = ", nrow(htr6_df)),
    x        = "Depression (EQ-5D)",
    y        = "HTR6 expression (VST-normalized)"
  ) +
  theme_bw() +
  theme(
    panel.grid      = element_blank(),
    axis.text       = element_text(size = 12, colour = "black"),
    axis.title      = element_text(size = 13),
    plot.title      = element_text(size = 14),
    legend.position = "none"
  )

print(p_htr6)

# To save:
# ggsave("HTR6_depression_expression.png", plot = p_htr6,
#        width = 7, height = 5, units = "in", dpi = 300)

cairo_pdf("HTR6_depression_expression.pdf", width = 7, height = 5)
print(p_htr6)
dev.off()

################################################################################
# SECTION 25: IFI27 HIGH-EXPRESSERS — DEPRESSION PREVALENCE BY SEX
#
# Tests whether the sex difference in depression prevalence is concentrated
# among patients with high IFI27 expression, using Fisher's exact test.
#
# IFI27 high/low is defined relative to the sex-specific median, so that
# the threshold is calibrated separately for males and females.
#
# Requires: expr_df_clean (from Section 19)
################################################################################

# --- Compute sex-specific IFI27 medians --------------------------------------

ifi27_median_male <- median(
  expr_df_clean$expr[expr_df_clean$gene == "IFI27" &
                       expr_df_clean$sex == "Male"],
  na.rm = TRUE
)

ifi27_median_female <- median(
  expr_df_clean$expr[expr_df_clean$gene == "IFI27" &
                       expr_df_clean$sex == "Female"],
  na.rm = TRUE
)

cat("Sex-specific IFI27 medians:\n")
cat("  Male:  ", ifi27_median_male,   "\n")
cat("  Female:", ifi27_median_female, "\n")


# --- Classify samples as IFI27 High/Low (sex-specific threshold) -------------

ifi27_depr <- expr_df_clean %>%
  dplyr::filter(gene == "IFI27") %>%
  tidyr::drop_na(eurodepr) %>%
  dplyr::mutate(
    ifi27_group = dplyr::case_when(
      sex == "Male"   & expr >  ifi27_median_male   ~ "High",
      sex == "Male"   & expr <= ifi27_median_male   ~ "Low",
      sex == "Female" & expr >  ifi27_median_female ~ "High",
      sex == "Female" & expr <= ifi27_median_female ~ "Low"
    )
  )

cat("\nIFI27 group distribution by sex:\n")
print(table(ifi27_depr$sex, ifi27_depr$ifi27_group, useNA = "always"))


# --- Fisher's exact test: sex ~ depression in High-IFI27 patients only -------

elevated_ifi27 <- dplyr::filter(ifi27_depr, ifi27_group == "High")
fisher_result  <- fisher.test(table(elevated_ifi27$sex, elevated_ifi27$eurodepr))

cat("\nFisher's exact test — sex × depression in IFI27-high patients:\n")
print(fisher_result)


# --- Stacked bar plot: depression proportion by sex in IFI27-high patients ---

p_ifi27_depr <- ggplot(elevated_ifi27, aes(x = sex, fill = eurodepr)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("No" = "#456787", "Yes" = "#7E1811")) +
  annotate("text", x = 1.5, y = 1.08,
           label = paste0("Fisher's p = ", round(fisher_result$p.value, 4)),
           size = 5, hjust = 0.5, colour = "grey20", lineheight = 1.4) +
  coord_cartesian(clip = "off") +
  labs(
    title = "Depression Prevalence in Elevated IFI27 Patients\nby Sex",
    x     = "Sex",
    y     = "Proportion",
    fill  = "Depression"
  ) +
  theme_bw() +
  theme(
    panel.grid   = element_blank(),
    axis.title.x = element_text(size = 18, color = "black"),
    axis.title.y = element_text(size = 18, color = "black"),
    axis.text.x  = element_text(size = 18, color = "black"),
    axis.text.y  = element_text(size = 18, color = "black"),
    plot.title   = element_text(size = 20, color = "black", hjust = 0.5),
    legend.title = element_text(size = 14, color = "black"),
    legend.text  = element_text(size = 13, color = "black"),
    plot.margin  = margin(t = 40, r = 10, b = 10, l = 30)
  )

print(p_ifi27_depr)


# To save:
# ggsave("IFI27_high_depression_by_sex.png", plot = p_ifi27_depr,
#        width = 7.2, height = 5, units = "in", dpi = 300)

# cairo_pdf("IFI27_high_depression_by_sex.pdf", width = 7.2, height = 5)
# print(p_ifi27_depr)
# dev.off()

################################################################################
# SECTION 26: CORRELATION ANALYSES
#
# Three parts:
#   A: Mini correlation — Depression & Demographics (race excluded)
#   B: Mini correlation — Depression & Demographics (ethnicity included)
#   C: Race association table (Kruskal-Wallis + Chi-square)
#   D: Full correlation matrix — all patient metadata
#
# Requires: combined_meta (from Section X)
################################################################################

# ── Step 1: Define helper function — masked corrplot ─────────────────────────

draw_masked_corrplot <- function(cor_mat, p_mat, coef_disp,
                                 tl.cex     = 0.9,
                                 number.cex = 0.9,
                                 cl.cex     = 1.0,
                                 title      = "",
                                 mar        = c(0, 0, 2, 0)) {
  corrplot::corrplot(
    cor_mat,
    method      = "color",
    type        = "upper",
    order       = "original",
    tl.col      = "black",
    tl.cex      = tl.cex,
    col         = colorRampPalette(c("#456787", "white", "#7E1811"))(200),
    p.mat       = p_mat,
    sig.level   = 0.05,
    insig       = "pch",
    pch         = 4,
    pch.cex     = 2,
    pch.col     = "black",
    addCoef.col = NULL,
    cl.cex      = cl.cex,
    title       = title,
    mar         = mar
  )
  
  # Overlay rho values for significant cells only
  n <- nrow(cor_mat)
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (j > i && !is.na(coef_disp[i, j])) {
        text(
          x      = j,
          y      = n + 1 - i,
          labels = formatC(coef_disp[i, j], digits = 2, format = "f"),
          col    = "black",
          cex    = number.cex
        )
      }
    }
  }
}


# ── Step 2: Mini correlation — Depression & Demographics (race excluded) ──────

cols_subset_a <- c(
  "age", "sex",
  "mrs", "euroqol_anxdepr", "cesd_total"
  # race      -> removed: nominal, 5 groups, not suitable for Spearman
  # ethnicity -> removed: collinear with race (chi-sq p = 0.0005), n = 9 Hispanic
)

num_meta_a <- combined_meta %>%
  dplyr::select(any_of(cols_subset_a)) %>%
  dplyr::mutate(
    sex = dplyr::case_when(
      sex == "Male"   ~ 0,   # Male = 0, Female = 1
      sex == "Female" ~ 1,
      TRUE            ~ NA_real_
    )
  ) %>%
  dplyr::select(where(is.numeric))

cat("Part A variables:\n")
print(colnames(num_meta_a))
cat("\nN rows:", nrow(num_meta_a), "\n")

cat("\nMissingness:\n")
num_meta_a %>%
  dplyr::summarise(across(everything(), ~ sum(!is.na(.x)))) %>%
  tidyr::pivot_longer(everything(),
                      names_to  = "variable",
                      values_to = "n_nonNA") %>%
  dplyr::mutate(pct_complete = round(n_nonNA / nrow(num_meta_a) * 100, 1)) %>%
  print(n = Inf)

readable_labels_a <- c(
  age             = "Age",
  sex             = "Sex",
  mrs             = "mRS",
  euroqol_anxdepr = "EQ-5D\nAnxiety/Depression",
  cesd_total      = "CESD Total"
)

cor_mat_a <- cor(num_meta_a, use = "pairwise.complete.obs", method = "spearman")

n_vars_a    <- ncol(num_meta_a)
var_names_a <- colnames(num_meta_a)

p_mat_a <- matrix(NA_real_, nrow = n_vars_a, ncol = n_vars_a,
                  dimnames = list(var_names_a, var_names_a))

for (i in seq_len(n_vars_a)) {
  for (j in seq_len(n_vars_a)) {
    if (i == j) {
      p_mat_a[i, j] <- 0
    } else {
      complete_idx <- complete.cases(num_meta_a[, i], num_meta_a[, j])
      n_complete   <- sum(complete_idx)
      if (n_complete >= 5) {
        test_result   <- cor.test(
          num_meta_a[complete_idx, i],
          num_meta_a[complete_idx, j],
          method = "spearman",
          exact  = FALSE
        )
        p_mat_a[i, j] <- test_result$p.value
      }
    }
  }
}

ordered_vars_a <- names(readable_labels_a)
ordered_vars_a <- ordered_vars_a[ordered_vars_a %in% rownames(cor_mat_a)]

cor_mat_a <- cor_mat_a[ordered_vars_a, ordered_vars_a]
p_mat_a   <- p_mat_a[ordered_vars_a,   ordered_vars_a]

stopifnot(identical(rownames(cor_mat_a), rownames(p_mat_a)))

new_labels_a <- dplyr::coalesce(readable_labels_a[colnames(cor_mat_a)],
                                colnames(cor_mat_a))

colnames(cor_mat_a) <- rownames(cor_mat_a) <- new_labels_a
colnames(p_mat_a)   <- rownames(p_mat_a)   <- new_labels_a

stopifnot(identical(colnames(cor_mat_a), colnames(p_mat_a)))

coef_display_a <- cor_mat_a
coef_display_a[p_mat_a >= 0.05] <- NA
diag(coef_display_a) <- NA

draw_masked_corrplot(
  cor_mat    = cor_mat_a,
  p_mat      = p_mat_a,
  coef_disp  = coef_display_a,
  tl.cex     = 0.9,
  number.cex = 0.9,
  cl.cex     = 1.5,
  title      = "Spearman Correlation - Depression & Demographics",
  mar        = c(0, 0, 2, 0)
)

# png("correlation_matrix_depression_subset_2.png",
#     width = 10, height = 10, units = "in", res = 300)

draw_masked_corrplot(
  cor_mat    = cor_mat_a,
  p_mat      = p_mat_a,
  coef_disp  = coef_display_a,
  tl.cex     = 1.4,
  number.cex = 1.2,
  cl.cex     = 1.2,
  title      = "",
  mar        = c(0, 0, 0, 0)
)
dev.off()
cat("Saved: correlation_matrix_depression_subset_2.png\n")

# Preview
draw_masked_corrplot(
  cor_mat    = cor_mat_a,
  p_mat      = p_mat_a,
  coef_disp  = coef_display_a,
  tl.cex     = 0.9,
  number.cex = 0.9,
  cl.cex     = 1.5,
  title      = "Spearman Correlation - Depression & Demographics",
  mar        = c(0, 0, 2, 0)
)

# Save
cairo_pdf("correlation_matrix_depression_subset_2.pdf", width = 10, height = 10)
draw_masked_corrplot(
  cor_mat    = cor_mat_a,
  p_mat      = p_mat_a,
  coef_disp  = coef_display_a,
  tl.cex     = 1.4,
  number.cex = 1.2,
  cl.cex     = 1.2,
  title      = "",
  mar        = c(0, 0, 0, 0)
)
dev.off()
cat("Saved: correlation_matrix_depression_subset_2.pdf\n")


# ── Step 3: Mini correlation — Depression & Demographics (ethnicity included) ─

num_meta_b <- combined_meta %>%
  dplyr::select(any_of(c("age", "sex", "ethnicity",
                         "mrs", "euroqol_anxdepr", "cesd_total"))) %>%
  dplyr::mutate(
    sex = dplyr::case_when(
      sex == "Male"   ~ 0,
      sex == "Female" ~ 1,
      TRUE            ~ NA_real_
    ),
    ethnicity = dplyr::case_when(
      ethnicity == 1 ~ 1,   # Hispanic = 1
      ethnicity == 2 ~ 0,   # Non-Hispanic = 0
      TRUE           ~ NA_real_
    )
  ) %>%
  dplyr::select(where(is.numeric))

cat("Part B variables:\n")
print(colnames(num_meta_b))
cat("\nN rows:", nrow(num_meta_b), "\n")

cat("\nMissingness:\n")
num_meta_b %>%
  dplyr::summarise(across(everything(), ~ sum(!is.na(.x)))) %>%
  tidyr::pivot_longer(everything(),
                      names_to  = "variable",
                      values_to = "n_nonNA") %>%
  dplyr::mutate(pct_complete = round(n_nonNA / nrow(num_meta_b) * 100, 1)) %>%
  print(n = Inf)

readable_labels_b <- c(
  age             = "Age",
  sex             = "Sex",
  ethnicity       = "Ethnicity\n(Hispanic)",
  mrs             = "mRS",
  euroqol_anxdepr = "EQ-5D\nAnxiety/Depression",
  cesd_total      = "CESD Total"
)

cor_mat_b <- cor(num_meta_b, use = "pairwise.complete.obs", method = "spearman")

n_vars_b    <- ncol(num_meta_b)
var_names_b <- colnames(num_meta_b)

p_mat_b <- matrix(NA_real_, nrow = n_vars_b, ncol = n_vars_b,
                  dimnames = list(var_names_b, var_names_b))

for (i in seq_len(n_vars_b)) {
  for (j in seq_len(n_vars_b)) {
    if (i == j) {
      p_mat_b[i, j] <- 0
    } else {
      complete_idx <- complete.cases(num_meta_b[, i], num_meta_b[, j])
      n_complete   <- sum(complete_idx)
      if (n_complete >= 5) {
        test_result   <- cor.test(
          num_meta_b[complete_idx, i],
          num_meta_b[complete_idx, j],
          method = "spearman",
          exact  = FALSE
        )
        p_mat_b[i, j] <- test_result$p.value
      }
    }
  }
}

ordered_vars_b <- names(readable_labels_b)
ordered_vars_b <- ordered_vars_b[ordered_vars_b %in% rownames(cor_mat_b)]

cor_mat_b <- cor_mat_b[ordered_vars_b, ordered_vars_b]
p_mat_b   <- p_mat_b[ordered_vars_b,   ordered_vars_b]

stopifnot(identical(rownames(cor_mat_b), rownames(p_mat_b)))

new_labels_b <- dplyr::coalesce(readable_labels_b[colnames(cor_mat_b)],
                                colnames(cor_mat_b))

colnames(cor_mat_b) <- rownames(cor_mat_b) <- new_labels_b
colnames(p_mat_b)   <- rownames(p_mat_b)   <- new_labels_b

stopifnot(identical(colnames(cor_mat_b), colnames(p_mat_b)))

coef_display_b <- cor_mat_b
coef_display_b[p_mat_b >= 0.05] <- NA
diag(coef_display_b) <- NA

png("correlation_matrix_depression_subset.png",
    width = 10, height = 10, units = "in", res = 300)
draw_masked_corrplot(
  cor_mat    = cor_mat_b,
  p_mat      = p_mat_b,
  coef_disp  = coef_display_b,
  tl.cex     = 1.4,
  number.cex = 1.2,
  cl.cex     = 1.2,
  title      = "",
  mar        = c(0, 0, 0, 0)
)
dev.off()
cat("Saved: correlation_matrix_depression_subset.png\n")


# ── Step 4: Race association table ───────────────────────────────────────────

race_labels <- c(
  "1" = "American Indian/Alaska Native",
  "2" = "Asian",
  "3" = "Black or African American",
  "4" = "Native Hawaiian/Pacific Islander",
  "5" = "White or Caucasian"
)

race_df <- combined_meta %>%
  dplyr::mutate(
    race_label = dplyr::recode(as.character(race), !!!race_labels),
    race_label = factor(race_label)
  ) %>%
  dplyr::filter(!is.na(race_label))

cat("Race associations:\n")

kw_vars <- c("age", "mrs", "euroqol_anxdepr", "cesd_total")

kw_results <- purrr::map_dfr(kw_vars, function(v) {
  df_v <- race_df %>% dplyr::filter(!is.na(.data[[v]]))
  kt   <- kruskal.test(as.formula(paste(v, "~ race_label")), data = df_v)
  data.frame(
    variable = v,
    n        = nrow(df_v),
    H        = round(kt$statistic, 2),
    df       = kt$parameter,
    p.value  = round(kt$p.value, 4)
  )
})

cat("Kruskal-Wallis: continuous/ordinal variables by race\n")
print(kw_results)

cat_vars <- c("sex", "ethnicity")

chi_results <- purrr::map_dfr(cat_vars, function(v) {
  df_v <- race_df %>% dplyr::filter(!is.na(.data[[v]]))
  ct   <- chisq.test(table(df_v[[v]], df_v$race_label), simulate.p.value = TRUE)
  data.frame(
    variable = v,
    n        = nrow(df_v),
    X2       = round(ct$statistic, 2),
    p.value  = round(ct$p.value, 4)
  )
})

cat("\nChi-square: categorical variables by race\n")
print(chi_results)


# ── Step 5: Full correlation matrix — all patient metadata ────────────────────

cols_for_cor <- c(
  "age", "sex", "mrs", "hilo_mrs", "initial_gcs", "initial_nihss",
  "hist_mrs", "len_stay", "nih_total", "mmse13_qsorres", "barthel_total",
  "euroqol_mobility", "euroqol_selfcare", "euroqol_usualact",
  "euroqol_paindisc", "euroqol_anxdepr", "hlthstat",
  "iqdcode_avg_pre", "iqdcode_avg_post",
  "story", "dgs01_qsorres", "dgs02_qsorres", "digord", "delay",
  "content", "fluency", "aud_verbal", "seq_command", "repetition",
  "obj_naming", "aphasia_sum", "apraxia",
  "adl_total", "tics_total", "dejong_total", "cesd_total",
  "sbp", "dbp", "weight_lb", "age_ich", "ses4",
  "hx1", "hx2", "hx21", "hx21d", "hx3", "dm1", "race", "ethnicity",
  "ich_hemisphere_adj", "ich_location_adj"
)

num_meta <- combined_meta %>%
  dplyr::select(any_of(cols_for_cor)) %>%
  dplyr::mutate(
    sex = dplyr::case_when(
      sex == "Male"   ~ 1,
      sex == "Female" ~ 2,
      TRUE            ~ NA_real_
    ),
    hilo_mrs = dplyr::case_when(
      hilo_mrs == "Low"  ~ 1,
      hilo_mrs == "High" ~ 2,
      TRUE               ~ NA_real_
    ),
    ich_hemisphere_adj = dplyr::case_when(
      ich_hemisphere_adj == "Left"  ~ 1,
      ich_hemisphere_adj == "Right" ~ 2,
      TRUE                          ~ NA_real_
    ),
    ich_location_adj = dplyr::case_when(
      ich_location_adj == "Deep"  ~ 1,
      ich_location_adj == "Lobar" ~ 2,
      TRUE                        ~ NA_real_
    )
  ) %>%
  dplyr::select(where(is.numeric))

always_keep <- c("ich_hemisphere_adj", "ich_location_adj")
n_total     <- nrow(num_meta)

num_meta <- num_meta %>%
  dplyr::select(
    where(~ sum(!is.na(.x)) > n_total / 2) | any_of(always_keep)
  )

num_meta <- num_meta %>%
  dplyr::select(
    where(~ isTRUE(sd(.x, na.rm = TRUE) > 0)) | any_of(always_keep)
  )

cat("Variables kept:\n")
print(colnames(num_meta))
cat("\nN variables:", ncol(num_meta), "\n")

num_meta %>%
  dplyr::summarise(across(everything(), ~ sum(!is.na(.x)))) %>%
  tidyr::pivot_longer(everything(),
                      names_to  = "variable",
                      values_to = "n_nonNA") %>%
  dplyr::mutate(pct_complete = round(n_nonNA / n_total * 100, 1)) %>%
  dplyr::arrange(pct_complete) %>%
  print(n = Inf)

cat("\nVerifying combined_meta is unchanged:\n")
print(table(combined_meta$sex,                useNA = "always"))
print(table(combined_meta$hilo_mrs,           useNA = "always"))
print(table(combined_meta$ich_hemisphere_adj, useNA = "always"))
print(table(combined_meta$ich_location_adj,   useNA = "always"))

readable_labels <- c(
  age                = "Age",
  sex                = "Sex",
  mrs                = "mRS",
  hilo_mrs           = "mRS (Hi/Lo)",
  initial_gcs        = "Initial GCS",
  initial_nihss      = "Initial NIHSS",
  hist_mrs           = "Historical mRS",
  len_stay           = "Length of Stay",
  nih_total          = "NIH Total",
  mmse13_qsorres     = "MMSE",
  barthel_total      = "Barthel Index",
  euroqol_mobility   = "EQ5D Mobility",
  euroqol_selfcare   = "EQ5D Self-Care",
  euroqol_usualact   = "EQ5D Usual Activity",
  euroqol_paindisc   = "EQ5D Pain",
  euroqol_anxdepr    = "EQ5D Anxiety/Depression",
  hlthstat           = "Health Status",
  iqdcode_avg_pre    = "IQCODE Pre",
  iqdcode_avg_post   = "IQCODE Post",
  story              = "Story Recall",
  dgs01_qsorres      = "Digit Span Fwd",
  dgs02_qsorres      = "Digit Span Bwd",
  digord             = "Digit Ordering",
  delay              = "Delayed Recall",
  content            = "WAB Content",
  fluency            = "WAB Fluency",
  aud_verbal         = "WAB Auditory",
  seq_command        = "WAB Sequential",
  repetition         = "WAB Repetition",
  obj_naming         = "WAB Naming",
  aphasia_sum        = "Aphasia Sum",
  apraxia            = "Apraxia",
  adl_total          = "ADL Total",
  tics_total         = "TICS Total",
  dejong_total       = "De Jong Loneliness",
  cesd_total         = "CESD Total",
  sbp                = "Systolic BP",
  dbp                = "Diastolic BP",
  weight_lb          = "Weight (lbs)",
  age_ich            = "Age at ICH",
  ses4               = "Marital Status",
  hx1                = "Hx Hypertension",
  hx2                = "Hx Diabetes",
  hx3                = "Hx Stroke",
  hx21               = "Hx Depression",
  hx21d              = "Hx Suicide",
  dm1                = "Diabetes",
  race               = "Race",
  ethnicity          = "Ethnicity",
  ich_hemisphere_adj = "ICH Hemisphere",
  ich_location_adj   = "ICH Location"
)

cor_mat_raw <- cor(num_meta, use = "pairwise.complete.obs", method = "spearman")

n_vars_d    <- ncol(num_meta)
var_names_d <- colnames(num_meta)

p_mat_raw <- matrix(NA_real_, nrow = n_vars_d, ncol = n_vars_d,
                    dimnames = list(var_names_d, var_names_d))

for (i in seq_len(n_vars_d)) {
  for (j in seq_len(n_vars_d)) {
    if (i == j) {
      p_mat_raw[i, j] <- 0
    } else {
      complete_idx <- complete.cases(num_meta[, i], num_meta[, j])
      n_complete   <- sum(complete_idx)
      if (n_complete >= 5) {
        test_result     <- cor.test(
          num_meta[complete_idx, i],
          num_meta[complete_idx, j],
          method = "spearman",
          exact  = FALSE
        )
        p_mat_raw[i, j] <- test_result$p.value
      }
    }
  }
}

stopifnot(identical(colnames(cor_mat_raw), colnames(p_mat_raw)))

ordered_vars_d <- names(readable_labels)
ordered_vars_d <- ordered_vars_d[ordered_vars_d %in% rownames(cor_mat_raw)]

cor_mat_raw <- cor_mat_raw[ordered_vars_d, ordered_vars_d]
p_mat_raw   <- p_mat_raw[ordered_vars_d,   ordered_vars_d]

stopifnot(identical(rownames(cor_mat_raw), rownames(p_mat_raw)))

new_labels_d <- dplyr::coalesce(readable_labels[colnames(cor_mat_raw)],
                                colnames(cor_mat_raw))

colnames(cor_mat_raw) <- rownames(cor_mat_raw) <- new_labels_d
colnames(p_mat_raw)   <- rownames(p_mat_raw)   <- new_labels_d

cor_mat <- cor_mat_raw
p_mat   <- p_mat_raw

stopifnot(identical(colnames(cor_mat), colnames(p_mat)))
cat("Labels applied. Variables:\n")
print(colnames(cor_mat))

corrplot::corrplot(
  cor_mat,
  method      = "color",
  type        = "upper",
  order       = "original",
  tl.col      = "black",
  tl.cex      = 0.7,
  col         = colorRampPalette(c("#456787", "white", "#7E1811"))(200),
  p.mat       = p_mat,
  sig.level   = 0.05,
  insig       = "pch",
  pch         = 4,
  pch.cex     = 0.8,
  pch.col     = "grey60",
  addCoef.col = NULL,
  title       = "Spearman Correlation Matrix - Patient Metadata",
  mar         = c(0, 0, 2, 0)
)

sex_cor_report <- data.frame(
  variable = colnames(cor_mat),
  rho      = cor_mat["Sex", ],
  p_value  = p_mat["Sex", ]
) %>%
  dplyr::filter(variable != "Sex") %>%
  dplyr::mutate(
    p_adj     = p.adjust(p_value, method = "BH"),
    rho       = round(rho,     3),
    p_value   = round(p_value, 4),
    p_adj     = round(p_adj,   4),
    direction = ifelse(rho > 0, "higher in Female", "higher in Male")
  ) %>%
  dplyr::arrange(p_value)

cat("\nSex correlations (ranked by p-value):\n")
print(sex_cor_report)

depr_label <- colnames(cor_mat)[
  grepl("Anxiety|anxdepr|EQ5D Anx", colnames(cor_mat), ignore.case = TRUE)
][1]

cat("\nUsing depression label:", depr_label, "\n")

depr_cor_report <- data.frame(
  variable = colnames(cor_mat),
  rho      = cor_mat[depr_label, ],
  p_value  = p_mat[depr_label, ]
) %>%
  dplyr::filter(variable != depr_label) %>%
  dplyr::mutate(
    p_adj     = p.adjust(p_value, method = "BH"),
    rho       = round(rho,     3),
    p_value   = round(p_value, 4),
    p_adj     = round(p_adj,   4),
    direction = ifelse(rho > 0, "higher in Depressed", "lower in Depressed")
  ) %>%
  dplyr::arrange(p_value)

cat("\nDepression (EQ5D) correlations (ranked by p-value):\n")
print(depr_cor_report)

png("correlation_matrix_final_new3.png",
    width = 20, height = 20, units = "in", res = 300)

corrplot::corrplot(
  cor_mat,
  method      = "color",
  type        = "upper",
  order       = "original",
  tl.col      = "black",
  tl.cex      = 2,
  cl.cex      = 2.5,
  cl.ratio    = 0.1,
  cl.length   = 5,
  col         = colorRampPalette(c("#456787", "white", "#7E1811"))(200),
  p.mat       = p_mat,
  sig.level   = 0.05,
  insig       = "pch",
  pch         = 4,
  pch.cex     = 0.8,
  pch.col     = "grey60",
  addCoef.col = NULL,
  mar         = c(0, 0, 0, 0)
)

dev.off()
cat("Saved: correlation_matrix_final_new3.png\n")

cairo_pdf("correlation_matrix_final_new3.pdf", width = 20, height = 20)

corrplot::corrplot(
  cor_mat,
  method      = "color",
  type        = "upper",
  order       = "original",
  tl.col      = "black",
  tl.cex      = 2,
  cl.cex      = 2.5,
  cl.ratio    = 0.1,
  cl.length   = 5,
  col         = colorRampPalette(c("#456787", "white", "#7E1811"))(200),
  p.mat       = p_mat,
  sig.level   = 0.05,
  insig       = "pch",
  pch         = 4,
  pch.cex     = 0.8,
  pch.col     = "grey60",
  addCoef.col = NULL,
  mar         = c(0, 0, 0, 0)
)

dev.off()
cat("Saved: correlation_matrix_final_new3.pdf\n")
