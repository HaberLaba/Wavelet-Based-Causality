# --- Gene–Cell
library(igraph)
library(ggraph)
library(ggplot2)
library(dplyr)
library(cowplot)
library(tidyr)
library(grid)

setwd("/Users/oliviaxu/Desktop/wavelet_project/without_dwt/without_dwt_matrices")

edges <- read.csv("without_dwt_1_5.csv", stringsAsFactors = FALSE, check.names = FALSE) %>%
  transmute(
    source = trimws(as.character(Source)),
    target = trimws(as.character(Target)),
    weight = suppressWarnings(as.numeric(Signed_weight))
  )

gc_edges <- edges %>%
  filter(!grepl("^Celltype_", source, ignore.case = TRUE) &
           grepl("^Celltype_", target, ignore.case = TRUE)) %>%
  mutate(
    target = sub("^Celltype_\\s*", "", target, ignore.case = TRUE),
    target = trimws(gsub("\\s+", " ", target))
  )

gc_edges <- gc_edges %>%
  group_by(source, target) %>%
  summarise(weight = mean(weight, na.rm = TRUE), .groups = "drop")

gene_names <- sort(unique(gc_edges$source))
cell_names <- sort(unique(gc_edges$target))

cell_x_pos <- 0.9
gene_x_pos <- 1.6

center_seq <- function(n, step = 1) {
  if (n <= 0) return(numeric(0))
  if (n == 1) return(0)
  offset <- (n - 1) / 2
  (seq_len(n) - 1 - offset) * step
}

col_span_gene <- 18
col_span_cell <- 18

dy_gene <- if (length(gene_names) <= 1) 0 else col_span_gene / (length(gene_names) - 1)
dy_cell <- if (length(cell_names) <= 1) 0 else col_span_cell / (length(cell_names) - 1)

y_gene_vals <- center_seq(length(gene_names), step = dy_gene)
y_cell_vals <- center_seq(length(cell_names), step = dy_cell)

cells_df <- if (length(cell_names) > 0) {
  data.frame(
    name       = cell_names,
    x          = rep(cell_x_pos, length(cell_names)),
    y          = y_cell_vals,
    full_label = cell_names,
    type       = rep("Cell", length(cell_names)),
    color      = rep("#FFD580", length(cell_names)),
    stringsAsFactors = FALSE
  )
} else {
  data.frame(
    name = character(0), x = numeric(0), y = numeric(0),
    full_label = character(0), type = character(0), color = character(0),
    stringsAsFactors = FALSE
  )
}

genes_df <- if (length(gene_names) > 0) {
  data.frame(
    name       = gene_names,
    x          = rep(gene_x_pos, length(gene_names)),
    y          = y_gene_vals,
    full_label = gene_names,
    type       = rep("Gene", length(gene_names)),
    color      = rep("#56B4E9", length(gene_names)),
    stringsAsFactors = FALSE
  )
} else {
  data.frame(
    name = character(0), x = numeric(0), y = numeric(0),
    full_label = character(0), type = character(0), color = character(0),
    stringsAsFactors = FALSE
  )
}

layout_df <- rbind(cells_df, genes_df)

lab_dx <- 0.05
layout_df <- layout_df %>%
  mutate(
    x_label = ifelse(type == "Cell", x - lab_dx, x + lab_dx),
    hjust   = ifelse(type == "Cell", 1, 0)
  )

pad_x <- 0.5
pad_y <- 0.8

x_min <- cell_x_pos - pad_x
x_max <- gene_x_pos + pad_x

edge_df <- gc_edges %>%
  inner_join(layout_df %>% select(name, x_source = x, y_source = y), by = c("source" = "name")) %>%
  inner_join(layout_df %>% select(name, x_target = x, y_target = y), by = c("target" = "name"))

col_sep   <- abs(gene_x_pos - cell_x_pos)
gap_start <- 0.03 * col_sep
gap_end   <- 0.08 * col_sep

edge_df <- edge_df %>%
  mutate(
    sx = sign(x_target - x_source),
    x_start = x_source + sx * gap_start,
    y_start = y_source,
    x_end   = x_target - sx * gap_end,
    y_end   = y_target,
    edge_color = ifelse(weight > 0, "#009851", "#C52C32")
  )

min_width <- 0.25; max_width <- 1.25
w_in <- abs(edge_df$weight)
edge_df$edge_size <- if (length(w_in) > 1 && diff(range(w_in)) > 0) {
  min_width + (w_in - min(w_in)) * (max_width - min_width) / (max(w_in) - min(w_in))
} else {
  rep((min_width + max_width) / 2, length(w_in))
}

node_strength_df <- bind_rows(
  edge_df %>% transmute(name = source, w = abs(weight)),
  edge_df %>% transmute(name = target, w = abs(weight))
) %>%
  group_by(name) %>%
  summarise(strength = sum(w, na.rm = TRUE), .groups = "drop")

layout_df <- layout_df %>%
  left_join(node_strength_df, by = "name") %>%
  mutate(strength = ifelse(is.na(strength), 0, strength))

plot_obj <- ggplot() +
  geom_segment(
    data = edge_df,
    aes(x = x_start, y = y_start, xend = x_end, yend = y_end,
        color = edge_color, linewidth = edge_size),
    arrow = arrow(length = unit(0.4, "cm"), type = "closed", angle = 20),
    show.legend = FALSE
  ) +
  geom_point(
    data = layout_df,
    aes(x = x, y = y, color = color, size = strength)
  ) +
  scale_size_continuous(range = c(8, 18), guide = "none") +
  geom_text(
    data = layout_df,
    aes(x = x_label, y = y, label = full_label, hjust = hjust),
    size = 5, color = "black", fontface = "bold"
  ) +
  scale_color_identity() +
  scale_linewidth_identity() +
  scale_x_continuous(limits = c(x_min, x_max), expand = expansion(mult = 0)) +
  scale_y_continuous(limits = c(y_min, y_max), expand = expansion(mult = 0)) +
  coord_cartesian(clip = "off") +
  theme_void() +
  theme(
    plot.margin = margin(t = 20, r = 30, b = 20, l = 30),
    plot.title = element_text(
      family = "Helvetica", size = 40, face = "bold", hjust = 0.5, vjust = 1
    )
  ) +
  ggtitle("3-Band and 4-Band Acinar")

ggsave(
  filename = "/Users/oliviaxu/Desktop/gene_cell_Acinar.png",
  plot = plot_obj,
  width = 15, height = 25, dpi = 600, bg = "white"
)











y0 <- 0.82
dy <- 0.08
ys <- function(i) y0 - dy * (i - 1)

sign_df <- data.frame(
  x_start = c(0, 0),
  x_end   = c(1.0, 1.0),
  y       = c(ys(1), ys(2)),
  color   = c("#009851", "#C52C32"),
  label   = c("positive causal link", "negative causal link")
)

abs_w  <- abs(edge_df$weight)
w_q    <- round(as.numeric(quantile(abs_w, probs = c(0, 0.5, 1), na.rm = TRUE)), 3)
names(w_q) <- c("min","med","max")

lw_min <- min(edge_df$edge_size, na.rm = TRUE)
lw_max <- max(edge_df$edge_size, na.rm = TRUE)
lw_mid <- (lw_min + lw_max) / 2

w_df <- data.frame(
  x_start = 0,
  x_end   = 1.0,
  y       = c(ys(3), ys(4), ys(5)),
  size    = c(lw_max, lw_mid, lw_min),
  label   = paste0("| weight | = ", c(w_q["max"], w_q["med"], w_q["min"]))
)

legend_plot <- ggplot() +
  geom_segment(
    data = sign_df,
    aes(x = x_start, y = y, xend = x_end, yend = y, color = color),
    arrow = arrow(length = unit(0.24, "cm"), type = "closed", angle = 20),
    linewidth = 1.4, lineend = "round", show.legend = FALSE
  ) +
  geom_text(
    data = sign_df,
    aes(x = 1.2, y = y, label = label),
    hjust = 0, size = 6.0, fontface = "bold", color = "black"
  ) +
  geom_segment(
    data = w_df,
    aes(x = x_start, y = y, xend = x_end, yend = y),
    color = "grey20",
    linewidth = w_df$size,
    lineend = "round",
    show.legend = FALSE
  ) +
  geom_text(
    data = w_df,
    aes(x = 1.2, y = y, label = label),
    hjust = 0, size = 6.0, fontface = "bold", color = "black"
  ) +
  scale_color_identity() +
  coord_cartesian(
    xlim = c(0, 3.0),
    ylim = c(ys(5) - 0.02, y0 + 0.02),
    expand = FALSE, clip = "off"
  ) +
  theme_void() +
  ggtitle("Legend") +
  theme(
    plot.title  = element_text(size = 30, face = "bold", hjust = 0.5),
    plot.margin = margin(2, 6, 2, 6),
    plot.background = element_rect(fill = "white", color = NA)
  )

legend_x <- 0.70
legend_y <- 0.35
legend_w <- 0.32
legend_h <- 0.28

final_plot <- ggdraw() +
  draw_plot(plot_obj, x = 0.08, y = 0, width = 1, height = 1) +
  draw_plot(legend_plot, x = legend_x, y = legend_y,
            width = legend_w, height = legend_h)
final_plot

ggsave(
  filename = "/Users/oliviaxu/Desktop/gene_cell_Acinar.png",
  plot = final_plot,
  width = 25,
  height = 20,
  dpi = 600,
  bg = "white"
)
