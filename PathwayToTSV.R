# Combines all lineage CSVs from pathway_results/
# into one navigable TSV for Excel, and also combines all
# per-dataset marker CSVs into a second TSV.

library(dplyr)
library(readr)
library(stringr)

# Settings
PROJECT_ROOT <- "/Users/nicholastong/Wavelet part 2 march 2026/New Dataset Wavelet Project"

PATHWAY_DIR  <- file.path(PROJECT_ROOT, "pathway_results_improved")

OUTPUT_TSV_PATHWAYS <- file.path(PROJECT_ROOT, "combined_pathway_annotations.tsv")
OUTPUT_TSV_MARKERS  <- file.path(PROJECT_ROOT, "combined_marker_table.tsv")

# Marker files expected in project root, produced in Main
# Example:
#   without_dwt_markers.csv
#   band3_view0_markers.csv
#   band3_view1_markers.csv
#   ...
MARKER_DIR <- PROJECT_ROOT

# ------------------------------------------------------------
# Desired dataset order
# ------------------------------------------------------------
dataset_order <- c(
  "without_dwt",
  "3band_view0",
  "3band_view1",
  "3band_view2",
  "4band_view0",
  "4band_view1",
  "4band_view2",
  "4band_view3"
)

dataset_display <- c(
  "without_dwt" = "without DWT",
  "3band_view0" = "3_0",
  "3band_view1" = "3_1",
  "3band_view2" = "3_2",
  "4band_view0" = "4_0",
  "4band_view1" = "4_1",
  "4band_view2" = "4_2",
  "4band_view3" = "4_3"
)

# Marker thresholds used in OliviaMain.
# Update these if you used different thresholds for different runs.
marker_threshold_map <- c(
  "without_dwt" = 0.1,
  "3band_view0" = 0.1,
  "3band_view1" = 0.1,
  "3band_view2" = 0.1,
  "4band_view0" = 0.1,
  "4band_view1" = 0.1,
  "4band_view2" = 0.1,
  "4band_view3" = 0.1
)

# ------------------------------------------------------------
# Helper: safely convert cluster to numeric for ordering
# ------------------------------------------------------------
safe_cluster_num <- function(x) {
  out <- suppressWarnings(as.numeric(as.character(x)))
  out[is.na(out)] <- Inf
  out
}

# ============================================================
# PART 1: COMBINE PATHWAY RESULT CSVs
# ============================================================

pathway_csv_files <- list.files(
  PATHWAY_DIR,
  pattern = "_top_enriched_pathways\\.csv$",
  full.names = TRUE
)

if (length(pathway_csv_files) == 0) {
  stop("No pathway CSV files found in: ", PATHWAY_DIR)
}

pathway_tables <- lapply(pathway_csv_files, function(f) {
  df <- read.csv(f, stringsAsFactors = FALSE)
  
  # infer dataset from file if missing
  if (!"dataset" %in% colnames(df)) {
    dataset_name <- sub("_top_enriched_pathways\\.csv$", "", basename(f))
    df$dataset <- dataset_name
  }
  
  df
})

combined_pathways <- bind_rows(pathway_tables)

required_pathway_cols <- c(
  "cluster",
  "pathway_name",
  "top_marker_genes",
  "adjusted_p_value",
  "significance_cutoff_used",
  "dataset"
)

missing_pathway_cols <- setdiff(required_pathway_cols, colnames(combined_pathways))
if (length(missing_pathway_cols) > 0) {
  stop(
    "Missing required columns in combined pathway data: ",
    paste(missing_pathway_cols, collapse = ", ")
  )
}

combined_pathways <- combined_pathways %>%
  mutate(
    dataset = as.character(dataset),
    cluster = as.character(cluster),
    dataset = factor(dataset, levels = dataset_order)
  ) %>%
  arrange(dataset, safe_cluster_num(cluster)) %>%
  mutate(
    navigation_id = paste(dataset_display[as.character(dataset)], cluster),
    final_annotations = ""
  ) %>%
  select(
    navigation_id,
    cluster,
    pathway_name,
    top_marker_genes,
    adjusted_p_value,
    significance_cutoff_used,
    dataset,
    final_annotations
  )

write.table(
  combined_pathways,
  file = OUTPUT_TSV_PATHWAYS,
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = TRUE
)

cat("Combined pathway TSV written to:\n", OUTPUT_TSV_PATHWAYS, "\n")

# ============================================================
# PART 2: COMBINE MARKER CSVs
# ============================================================
  
  marker_csv_files <- list.files(
    MARKER_DIR,
    pattern = "_markers\\.csv$",
    full.names = TRUE
  )
  
  if (length(marker_csv_files) == 0) {
    stop("No marker CSV files found in: ", MARKER_DIR)
  }
  
  # keep only datasets we care about
  marker_csv_files <- marker_csv_files[
    basename(marker_csv_files) %in% paste0(dataset_order, "_markers.csv")
  ]
  
  if (length(marker_csv_files) == 0) {
    stop("No matching marker CSV files found for the expected dataset names.")
  }
  
  marker_tables <- lapply(marker_csv_files, function(f) {
    df <- read.csv(f, stringsAsFactors = FALSE)
    
    dataset_name <- sub("_markers\\.csv$", "", basename(f))
    
    # Determine marker gene column
    gene_col_candidates <- c("gene", "Gene", "feature", "Feature")
    gene_col <- gene_col_candidates[gene_col_candidates %in% colnames(df)][1]
    if (is.na(gene_col)) {
      stop("No gene column found in marker file: ", basename(f))
    }
    
    # Determine cluster column
    cluster_col_candidates <- c("cluster", "Cluster")
    cluster_col <- cluster_col_candidates[cluster_col_candidates %in% colnames(df)][1]
    if (is.na(cluster_col)) {
      stop("No cluster column found in marker file: ", basename(f))
    }
    
    # Determine p-value column
    pval_col_candidates <- c("p_val_adj", "p_val", "pvalue", "PValue", "p_val_adj.1")
    pval_col <- pval_col_candidates[pval_col_candidates %in% colnames(df)][1]
    if (is.na(pval_col)) {
      stop("No p-value column found in marker file: ", basename(f))
    }
    
    out <- data.frame(
      marker_gene_name = df[[gene_col]],
      cluster = as.character(df[[cluster_col]]),
      threshold_used = marker_threshold_map[[dataset_name]],
      p_value = df[[pval_col]],
      dataset = dataset_name,
      stringsAsFactors = FALSE
    )
    
    out
  })
  
  combined_markers <- bind_rows(marker_tables) %>%
    mutate(
      dataset = factor(dataset, levels = dataset_order)
    ) %>%
    arrange(dataset, safe_cluster_num(cluster), p_value, marker_gene_name) %>%
    select(
      marker_gene_name,
      cluster,
      threshold_used,
      p_value,
      dataset
    )
  
  write.table(
    combined_markers,
    file = OUTPUT_TSV_MARKERS,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = TRUE
  )
  
  cat("Combined marker TSV written to:\n", OUTPUT_TSV_MARKERS, "\n")