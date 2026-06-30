setwd("/users/nicholastong/desktop/wavelet_project/without_dwt/without_dwt_matrices")


library(ggplot2)
library(dplyr)
library(tibble)
set.seed(123)

#Simulate 15 genes, 6 cells
num_genes <- 15
num_cells <- 6
gene_names <- paste0("GENE_", seq_len(num_genes))
cell_names <- paste0("CELL_", seq_len(num_cells))
gene_labels <- paste0("G", seq_len(num_genes))
cell_labels <- paste0("C", seq_len(num_cells))
min_width <- 0.3
max_width <- 1.25

#Simulate gene-gene edges
gene_gene_edges <- expand.grid(source = gene_names, target = gene_names, stringsAsFactors = FALSE)
gene_gene_edges <- gene_gene_edges[gene_gene_edges$source != gene_gene_edges$target, ]
gene_gene_edges <- gene_gene_edges[sample(nrow(gene_gene_edges), 30), ]
gene_gene_edges$weight <- round(rnorm(nrow(gene_gene_edges), 0, 0.9), 2)
gene_gene_edges <- gene_gene_edges[abs(gene_gene_edges$weight) > 0.18, ]
gene_gene_edges <- gene_gene_edges[!duplicated(gene_gene_edges[, c("source", "target")]), ]

cell_cell_edges <- expand.grid(source = cell_names, target = cell_names, stringsAsFactors = FALSE)
cell_cell_edges <- cell_cell_edges[cell_cell_edges$source != cell_cell_edges$target, ]
cell_cell_edges <- cell_cell_edges[sample(nrow(cell_cell_edges), 9), ]
cell_cell_edges$weight <- round(rnorm(nrow(cell_cell_edges), 0, 0.9), 2)
cell_cell_edges <- cell_cell_edges[abs(cell_cell_edges$weight) > 0.24, ]
cell_cell_edges <- cell_cell_edges[!duplicated(cell_cell_edges[, c("source", "target")]), ]

gene_cell_edges <- expand.grid(gene = gene_names, cell = cell_names, stringsAsFactors=FALSE)
gene_cell_edges <- gene_cell_edges[sample(nrow(gene_cell_edges), 20), ]
gene_cell_edges <- gene_cell_edges[!duplicated(gene_cell_edges[,c("gene","cell")]), ]
gene_cell_edges$weight <- round(rnorm(nrow(gene_cell_edges), 0, 1.15), 2)
gene_cell_edges <- gene_cell_edges[abs(gene_cell_edges$weight) > 0.22, ]
gene_cell_edges <- gene_cell_edges[!duplicated(gene_cell_edges[, c("gene", "cell")]), ]
gene_cell_edges <- gene_cell_edges %>% rename(source = gene, target = cell)

#Layout
gap_wide <- 0.17
panel_width <- 2.2
panel_left <- c(0, panel_width + gap_wide, 2*(panel_width + gap_wide))
x_panels <- c(panel_left[1], panel_left[2], panel_left[3])

# Update box y-limits allowing space for internal titles above graphs
box_ymin <- -1.65
box_ymax <- 1.7   # slightly taller box to fit title inside top
title_y <- box_ymax - 0.18  # title inside top, 0.18 units below top border

#Gene-gene graph (left panel)
angle_g <- seq(0, 2*pi, length.out = num_genes + 1)[1:num_genes]
x_gene <- cos(angle_g) * 0.85 + x_panels[1] + panel_width/2
y_gene <- sin(angle_g)
layout_gene <- tibble(
  name = gene_names, label = gene_labels,
  x = x_gene, y = y_gene, color = "#56B4E9"
)

#Cell-cell graph (right panel)
angle_c <- seq(0, 2*pi, length.out = num_cells + 1)[1:num_cells]
x_cell <- cos(angle_c) * 0.85 + x_panels[3] + panel_width/2
y_cell <- sin(angle_c)
layout_cell <- tibble(
  name = cell_names, label = cell_labels,
  x = x_cell, y = y_cell, color = "#FFD580"
)

#Bipartite gene-cell (center panel)
center_panel_left <- x_panels[2]
center_panel_right <- x_panels[2] + panel_width
x_gene_gc <- center_panel_left + 0.55
x_cell_gc <- center_panel_right - 0.55
y_gene_gc <- seq(1.2, -1.2, length.out = num_genes)
y_cell_gc <- seq(1.2, -1.2, length.out = num_cells)
layout_gene_gc <- tibble(
  name = gene_names, label = gene_labels, x = x_gene_gc, y = y_gene_gc, color = "#56B4E9"
)
layout_cell_gc <- tibble(
  name = cell_names, label = cell_labels, x = x_cell_gc, y = y_cell_gc, color = "#FFD580"
)
layout_gc <- bind_rows(layout_gene_gc, layout_cell_gc)

boxes <- tibble(
  xmin = x_panels,
  xmax = x_panels + panel_width,
  ymin = box_ymin,
  ymax = box_ymax
)

#Split segments
offset <- 0.13
make_split_segment <- function(layout, from, to, weight, offset) {
  x1 <- layout$x[layout$name == from]; y1 <- layout$y[layout$name == from]
  x2 <- layout$x[layout$name == to]; y2 <- layout$y[layout$name == to]
  dx <- x2 - x1; dy <- y2 - y1
  dist <- sqrt(dx^2 + dy^2)
  ux <- dx / dist; uy <- dy / dist
  x_start <- x1 + ux * offset; y_start <- y1 + uy * offset
  x_end <- x2 - ux * offset; y_end <- y2 - uy * offset
  edge_col <- ifelse(weight > 0, "#009851", "#C52C32")
  edge_size <- min(max_width, max(min_width, abs(weight)))
  data.frame(x_start = x_start, y_start = y_start, x_end = x_end, y_end = y_end,
             color = edge_col, size = edge_size, weight = weight)
}
gene_gene_edges_df <- bind_rows(
  mapply(function(f, t, w) make_split_segment(layout_gene, f, t, w, offset),
         gene_gene_edges$source, gene_gene_edges$target, gene_gene_edges$weight, SIMPLIFY = FALSE)
)
cell_cell_edges_df <- bind_rows(
  mapply(function(f, t, w) make_split_segment(layout_cell, f, t, w, offset),
         cell_cell_edges$source, cell_cell_edges$target, cell_cell_edges$weight, SIMPLIFY = FALSE)
)
gene_cell_edges_df <- gene_cell_edges %>%
  left_join(layout_gc %>% select(name, x_source = x, y_source = y), by = c("source" = "name")) %>%
  left_join(layout_gc %>% select(name, x_target = x, y_target = y), by = c("target" = "name")) %>%
  rowwise() %>% mutate(
    dx = x_target - x_source, dy = y_target - y_source, dist = sqrt(dx^2 + dy^2),
    ux = dx / dist, uy = dy / dist,
    x_start = x_source + ux * offset, y_start = y_source + uy * offset,
    x_end = x_target - ux * offset, y_end = y_target - uy * offset
  ) %>% ungroup()
gene_cell_edges_df$color <- ifelse(gene_cell_edges_df$weight > 0, "#009851", "#C52C32")
gene_cell_edges_df$size <- min_width + (abs(gene_cell_edges_df$weight) - min(abs(gene_cell_edges_df$weight))) *
  (max_width - min_width) / (max(abs(gene_cell_edges_df$weight)) - min(abs(gene_cell_edges_df$weight)) + 1e-9)

#Plot
p <- ggplot() +
  #Boxes
  geom_rect(data = boxes, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = NA, color = "black", linewidth = 1.35) +
  
  #Graph Titles
  geom_text(aes(x = x_panels + panel_width / 2, y = title_y,
                label = c("Gene-Gene", "Gene-Cell", "Cell-Cell")),
            fontface = "bold", family = "Helvetica", size = 6, color = "black", vjust = 1) +
  
  #Gene-gene panel
  geom_segment(data = gene_gene_edges_df,
               aes(x = x_start, y = y_start, xend = x_end, yend = y_end, color = color, size = size),
               arrow = arrow(length = unit(0.28, "cm"), type = "closed"),
               lineend = "round", show.legend = FALSE) +
  geom_point(data = layout_gene, aes(x = x, y = y, color = color), size = 7) +
  geom_text(data = layout_gene, aes(x = x, y = y, label = label),
            color = "black", family = "Helvetica", fontface = "bold", size = 3.2) +
  
  #Cell-cell panel
  geom_segment(data = cell_cell_edges_df,
               aes(x = x_start, y = y_start, xend = x_end, yend = y_end, color = color, size = size),
               arrow = arrow(length = unit(0.28, "cm"), type = "closed"),
               lineend = "round", show.legend = FALSE) +
  geom_point(data = layout_cell, aes(x = x, y = y, color = color), size = 7) +
  geom_text(data = layout_cell, aes(x = x, y = y, label = label),
            color = "black", family = "Helvetica", fontface = "bold", size = 3.2) +
  
  #Gene-cell center panel
  geom_segment(data = gene_cell_edges_df,
               aes(x = x_start, y = y_start, xend = x_end, yend = y_end, color = color, size = size),
               arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
               lineend = "round", show.legend = FALSE) +
  geom_point(data = layout_gene_gc, aes(x = x, y = y, color = color), size = 7) +
  geom_text(data = layout_gene_gc, aes(x = x, y = y, label = label),
            color = "black", family = "Helvetica", fontface = "bold", size = 3.2) +
  geom_point(data = layout_cell_gc, aes(x = x, y = y, color = color), size = 7) +
  geom_text(data = layout_cell_gc, aes(x = x, y = y, label = label),
            color = "black", family = "Helvetica", fontface = "bold", size = 3.2) +
  
  scale_color_identity() +
  scale_size(range = c(min_width, max_width)) +
  
  scale_x_continuous(limits = c(-0.3, max(x_panels + panel_width) + 0.6), expand = c(0,0)) +
  scale_y_continuous(limits = c(box_ymin - 0.35, box_ymax + 0.18), expand = c(0,0)) +
  
  theme_void() +
  ggtitle("Multilayer Causal Graph") +
  theme(
    plot.title = element_text(family = "Helvetica", size = 22, face = "bold", hjust = 0.5, vjust = -1),
    plot.margin = margin(15, 15, 15, 15)
  )

ggsave("simulatedmultilayergraph.png", p, width = 16, height = 7, dpi = 300)
print(p)
