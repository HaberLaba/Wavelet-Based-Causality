files <- c(
  "/Users/oliviaxu/Desktop/wavelet_project/without_dwt/without_dwt_annotations.csv",
  "/Users/oliviaxu/Desktop/wavelet_project/3band_2reg/wavelet_3band_2reg_view0_annotations.csv",
  "/Users/oliviaxu/Desktop/wavelet_project/3band_2reg/wavelet_3band_2reg_view1_annotations.csv",
  "/Users/oliviaxu/Desktop/wavelet_project/3band_2reg/wavelet_3band_2reg_view2_annotations.csv",
  "/Users/oliviaxu/Desktop/wavelet_project/4band_2reg/wavelet_4band_2reg_view0_annotations.csv",
  "/Users/oliviaxu/Desktop/wavelet_project/4band_2reg/wavelet_4band_2reg_view1_annotations.csv",
  "/Users/oliviaxu/Desktop/wavelet_project/4band_2reg/wavelet_4band_2reg_view2_annotations.csv",
  "/Users/oliviaxu/Desktop/wavelet_project/4band_2reg/wavelet_4band_2reg_view3_annotations.csv"
)

common_names <- NULL

for (i in seq_along(files)) {
  df <- read.csv(files[i], stringsAsFactors = FALSE)
  names_in_file <- df[[1]]
  
  if (is.null(common_names)) {
    common_names <- names_in_file
  } else {
    common_names <- intersect(common_names, names_in_file)
  }
}

print(common_names)

marker_files <- c(
  "/Users/oliviaxu/Desktop/wavelet_project/without_dwt/without_dwt_markers.csv",
  "/Users/oliviaxu/Desktop/wavelet_project/3band_2reg/wavelet_3band_2reg_view0_markers.csv",
  "/Users/oliviaxu/Desktop/wavelet_project/3band_2reg/wavelet_3band_2reg_view1_markers.csv",
  "/Users/oliviaxu/Desktop/wavelet_project/3band_2reg/wavelet_3band_2reg_view2_markers.csv",
  "/Users/oliviaxu/Desktop/wavelet_project/4band_2reg/wavelet_4band_2reg_view0_markers.csv",
  "/Users/oliviaxu/Desktop/wavelet_project/4band_2reg/wavelet_4band_2reg_view1_markers.csv",
  "/Users/oliviaxu/Desktop/wavelet_project/4band_2reg/wavelet_4band_2reg_view2_markers.csv",
  "/Users/oliviaxu/Desktop/wavelet_project/4band_2reg/wavelet_4band_2reg_view3_markers.csv"
)

top_k <- 10

view_gene_set <- function(path, k = top_k) {
  read_csv(path, show_col_types = FALSE) %>%
    mutate(cluster = as.character(cluster)) %>%
    group_by(cluster) %>%
    slice_max(order_by = avg_log2FC, n = k, with_ties = FALSE) %>%
    ungroup() %>%
    pull(gene) %>%
    unique()
}

per_view_sets <- lapply(marker_files, view_gene_set)
names(per_view_sets) <- basename(marker_files)

sizes <- sapply(per_view_sets, length)
cat("Per-view set sizes (union of per-cluster top", top_k, "):\n")
print(sizes)

common_genes <- Reduce(intersect, per_view_sets)
