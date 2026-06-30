# --- Gene–Gene
library(ggplot2)
library(dplyr)
library(grid)
library(cowplot)

setwd("/Users/nicholastong/Desktop/wavelet_project/4band_2reg/wavelet_4band_2reg_view2_matrices")

edges <- read.csv("4band_2reg_view2_2_5.csv", stringsAsFactors = FALSE, check.names = FALSE) %>%
  transmute(
    source = trimws(as.character(Source)),
    target = trimws(as.character(Target)),
    weight = suppressWarnings(as.numeric(Signed_weight))
  )

gg_edges <- edges %>%
  filter(!grepl("^Celltype_", source, ignore.case = TRUE) &
           !grepl("^Celltype_", target, ignore.case = TRUE)) %>%
  mutate(
    source = trimws(source),
    target = trimws(target),
    weight = suppressWarnings(as.numeric(weight))
  ) %>%
  group_by(source, target) %>%
  summarise(weight = mean(weight, na.rm = TRUE), .groups = "drop")

gene_names <- sort(unique(c(gg_edges$source, gg_edges$target)))
n_genes <- length(gene_names)

if (!exists("circle_radius")) circle_radius <- 5

if (n_genes == 0) {
  layout_df <- data.frame(
    name = character(0), x = numeric(0), y = numeric(0),
    label = character(0), full_label = character(0),
    stringsAsFactors = FALSE
  )
} else {
  theta <- if (n_genes == 1) 0 else seq(0, 2*pi, length.out = n_genes + 1)[-(n_genes+1)]
  layout_df <- data.frame(
    name       = gene_names,
    x          = circle_radius * cos(theta),
    y          = circle_radius * sin(theta),
    label      = paste0("G", seq_along(gene_names)),
    full_label = gene_names,
    stringsAsFactors = FALSE
  )
}

edge_df <- gg_edges %>%
  inner_join(layout_df %>% select(name, x_source = x, y_source = y),
             by = c("source" = "name")) %>%
  inner_join(layout_df %>% select(name, x_target = x, y_target = y),
             by = c("target" = "name")) %>%
  mutate(
    dx = x_target - x_source, dy = y_target - y_source,
    L  = pmax(sqrt(dx^2 + dy^2), 1e-9),
    ux = dx / L, uy = dy / L,
    off = pmin(0.04 * circle_radius, 0.05 * L),
    x_start = x_source + ux * off,  y_start = y_source + uy * off,
    x_end   = x_target - ux * off,  y_end   = y_target - uy * off,
    edge_color = ifelse(weight > 0, "#009851", "#C52C32")
  )

min_width <- 0.25; max_width <- 1.25
w_abs <- abs(edge_df$weight)
edge_df$edge_size <- if (length(w_abs) > 1 && diff(range(w_abs)) > 0) {
  min_width + (w_abs - min(w_abs)) * (max_width - min_width) / (max(w_abs) - min(w_abs))
} else {
  rep((min_width + max_width) / 2, length(w_abs))
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

g1_row <- subset(layout_df, label == "G1")
others <- subset(layout_df, label != "G1")

plot_gene_gene <- ggplot() +
  geom_segment(
    data = edge_df,
    aes(x = x_start, y = y_start, xend = x_end, yend = y_end,
        color = edge_color, linewidth = edge_size),
    arrow = arrow(length = unit(0.25, "cm"), type = "closed", angle = 20),
    lineend = "round",
    show.legend = FALSE
  ) +
  geom_point(
    data = others,
    aes(x = x, y = y, size = strength),
    shape = 21, color = "#56B4E9", fill = "#56B4E9", stroke = 0.7
  ) +
  geom_point(
    data = g1_row,
    aes(x = x, y = y, size = strength),
    shape = 24, color = "#7E57C2", fill = "#7E57C2", stroke = 0.9
  ) +
  scale_size_continuous(range = c(10, 20), guide = "none") +
  geom_text(data = layout_df, aes(x = x, y = y, label = label),
            size = 3.6, color = "black", fontface = "bold") +
  scale_color_identity() +
  scale_linewidth_identity() +
  coord_equal(clip = "off") +
  theme_void() +
  theme(
    plot.margin = margin(20, 20, 20, 20),
    plot.title  = element_text(family = "Helvetica", size = 40, face = "bold",
                               hjust = 0.5, vjust = 1)
  )

n_items <- nrow(layout_df)
per_col <- ceiling(n_items / 4)

leg_cols <- lapply(1:4, function(k) {
  idx_start <- (k - 1) * per_col + 1
  idx_end   <- min(k * per_col, n_items)
  if (idx_start > n_items) return(layout_df[0, , drop = FALSE])
  df <- layout_df[idx_start:idx_end, , drop = FALSE]
  df$legend_y <- rev(seq_len(nrow(df)))
  df
})

legend_x_point <- 0.5
legend_x_eq    <- legend_x_point + 0.5
legend_x_real  <- legend_x_eq + 0.3

make_node_legend_col <- function(df) {
  ggplot(df) +
    geom_text(aes(x = legend_x_point, y = legend_y, label = label),
              size = 6, color = "black", fontface = "bold", hjust = 0.5) +
    geom_text(aes(x = legend_x_eq, y = legend_y, label = "="),
              size = 6, color = "black", fontface = "bold", hjust = 0) +
    geom_text(aes(x = legend_x_real, y = legend_y, label = full_label),
              size = 6, color = "black", fontface = "bold", hjust = 0) +
    coord_cartesian(xlim = c(0, 4),
                    ylim = c(0, ifelse(nrow(df) > 0, nrow(df), 1) + 1),
                    expand = FALSE, clip = "off") +
    theme_void()
}

legend_cols_plots <- lapply(leg_cols, make_node_legend_col)

node_legend_combined <- cowplot::plot_grid(
  plotlist = legend_cols_plots,
  ncol = 4, align = "h"
)

legend_with_title <- cowplot::plot_grid(
  ggdraw() + draw_label("Legend", fontface = "bold", size = 30, hjust = 0.5),
  node_legend_combined,
  ggdraw() + draw_label(
    "Note: Purple triangle marks G1 (starting node); remaining genes arranged counterclockwise.",
    size = 14, fontface = "italic", hjust = 0.5
  ),
  ncol = 1,
  rel_heights = c(0.12, 1, 0.12)
)

final_gene_gene <- ggdraw() +
  draw_label(
    "High 3 (Lineage 5)",
    x = 0.5, y = 0.98,
    hjust = 0.5, vjust = 1,
    fontface = "bold", size = 40
  ) +
  draw_plot(plot_gene_gene, x = -0.2, y = 0, width = 1.05, height = 1) +
  draw_plot(legend_with_title, x = 0.72, y = 0.05, width = 0.28, height = 0.90)

final_gene_gene

#ggsave(
#  filename = "/Users/oliviaxu/Desktop/gene_gene_weights_scaled.png",
#  plot = final_gene_gene,
#  width = 30, height = 20, dpi = 600, bg = "white"
#)
