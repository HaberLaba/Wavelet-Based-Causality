library(plotly)
library(dplyr)
library(reticulate)
library(plotly)

#use_condaenv("r-reticulate", required = TRUE)

#use_python("/usr/bin/python3", required = TRUE)

setwd("/users/nicholastong/test/")


edges <- read.csv("3band_2reg_view0_1_5.csv", stringsAsFactors = FALSE)
colnames(edges) <- tolower(colnames(edges))
names(edges) <- gsub("signed_weight", "weight", names(edges), ignore.case = TRUE)
edges$source <- trimws(as.character(edges$source))
edges$target <- trimws(as.character(edges$target))
edges$weight <- suppressWarnings(as.numeric(edges$weight))

all_nodes <- unique(c(edges$source, edges$target))
cell_types <- all_nodes[grepl("celltype", all_nodes, ignore.case = TRUE)]
genes <- setdiff(all_nodes, cell_types)

gene_labels <- if (length(genes) > 0) paste0("G", seq_along(genes)) else character(0)
cell_labels <- if (length(cell_types) > 0) paste0("C", seq_along(cell_types)) else character(0)

rad_genes <- 220; rad_cells <- 65
z_genes <- 0; z_cells <- 80
label_offset <- 10

tg <- seq(0, 2*pi, length.out=length(genes)+1)[-1]
genes_coords <- data.frame(
  name=genes, label=gene_labels,
  x=rad_genes*cos(tg), y=rad_genes*sin(tg),
  z=z_genes, type="Gene"
)
tc <- seq(0, 2*pi, length.out=length(cell_types)+1)[-1]
cells_coords <- data.frame(
  name=cell_types, label=cell_labels,
  x=rad_cells*cos(tc), y=rad_cells*sin(tc),
  z=z_cells + label_offset, type="Cell"
)

nodes <- rbind(genes_coords, cells_coords)
rownames(nodes) <- nodes$name

node_degree_table <- table(c(edges$source, edges$target))
nodes$degree <- as.numeric(node_degree_table[nodes$name])
nodes$degree[is.na(nodes$degree)] <- 0

size_scale <- function(deg, base, max) {
  ifelse(deg <= 1, base, pmin(base + (deg - 1) * ((max - base) / 10), max))
}
nodes$size <- ifelse(nodes$type == "Gene",
                     size_scale(nodes$degree, 16, 32),
                     size_scale(nodes$degree, 24, 52))

edge_and_cone <- function(from, to, weight) {
  fx <- nodes[from, "x"]; fy <- nodes[from, "y"]; fz <- nodes[from, "z"]
  tx <- nodes[to, "x"]; ty <- nodes[to, "y"]; tz <- nodes[to, "z"]
  dx <- tx - fx; dy <- ty - fy; dz <- tz - fz
  dist <- sqrt(dx^2 + dy^2 + dz^2)
  if (dist == 0) return(NULL)
  ux <- dx / dist; uy <- dy / dist; uz <- dz / dist
  offset_start <- (nodes[from, "size"] / 2) * 0.7
  cone_length <- 6
  line_end_x <- tx - ux * cone_length
  line_end_y <- ty - uy * cone_length
  line_end_z <- tz - uz * cone_length
  fx <- fx + ux * offset_start
  fy <- fy + uy * offset_start
  fz <- fz + uz * offset_start
  raw_width <- (abs(weight)^0.9) / 2
  wid <- max(0.25, min(raw_width, 1.25))
  col <- if (weight >= 0) "#009851" else "#C52C32"
  line_trace <- list(
    x = c(fx, line_end_x), y = c(fy, line_end_y), z = c(fz, line_end_z),
    mode = "lines", type = "scatter3d",
    line = list(color = col, width = wid),
    hoverinfo = "text",
    text = sprintf("%s → %s (%.2f)", from, to, weight),
    showlegend = FALSE
  )
  cone_size <- if (nodes[to, "type"] == "Gene") 16 else 10
  cone_trace <- list(
    type = "cone",
    x = line_end_x, y = line_end_y, z = line_end_z,
    u = ux, v = uy, w = uz,
    sizemode = "absolute",
    sizeref = cone_size,
    anchor = "tip",
    showscale = FALSE,
    colorscale = list(c(0, col), c(1, col)),
    showlegend = FALSE
  )
  list(line = line_trace, cone = cone_trace)
}

p <- plot_ly(width = 1800, height = 900)

for(i in seq_len(nrow(edges))) {
  e <- edges[i,]
  if (e$source %in% nodes$name && e$target %in% nodes$name) {
    ec <- edge_and_cone(e$source, e$target, e$weight)
    if (!is.null(ec)) {
      p <- do.call(add_trace, c(list(p), ec$line))
      p <- do.call(add_trace, c(list(p), ec$cone))
    }
  }
}

p <- add_trace(p,
               x = genes_coords$x, y = genes_coords$y, z = genes_coords$z,
               type = "scatter3d", mode = "markers+text",
               text = genes_coords$label, textposition = "middle center",
               textfont = list(family = "Helvetica", size = 7, color = "black"),
               marker = list(size = nodes$size[match(genes_coords$name, nodes$name)],
                             color = "#0072B2", symbol = "circle",
                             line = list(color = "black", width = 1.7)),
               hoverinfo = "text", showlegend = FALSE)

p <- add_trace(p,
               x = cells_coords$x, y = cells_coords$y, z = cells_coords$z,
               type = "scatter3d", mode = "markers+text",
               text = cells_coords$label, textposition = "middle center",
               textfont = list(family = "Helvetica", size = 10, color = "black"),
               marker = list(size = nodes$size[match(cells_coords$name, nodes$name)],
                             color = "#FFD580", symbol = "circle",
                             line = list(color = "black", width = 2)),
               hoverinfo = "text", showlegend = FALSE)

mapping <- rbind(
  data.frame(label = genes_coords$label, name = genes_coords$name, type = "Gene"),
  data.frame(label = cells_coords$label, name = cells_coords$name, type = "Cell")
)

total <- nrow(mapping)
chunk_size <- ceiling(total / 4)
chunks <- split(1:total, ceiling(seq_along(1:total) / chunk_size))

for(chunk in chunks) {
  for(i in chunk) {
    display_name <- paste(mapping$label[i], "=", mapping$name[i])
    color <- if(mapping$type[i] == "Gene") "#0072B2" else "#FFD580"
    p <- add_trace(p,
                   x = 0, y = 0, z = 0,
                   type = "scatter3d",
                   mode = "markers",
                   opacity = 0,
                   marker = list(size = 8, color = color, symbol = "circle"),
                   name = display_name,
                   showlegend = TRUE,
                   hoverinfo = "none",
                   inherit = FALSE)
  }
}

scene_cfg <- list(
  xaxis = list(title = "", showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE, range = c(-250, 250)),
  yaxis = list(title = "", showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE, range = c(-250, 250)),
    zaxis = list(title = "", showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE, range = c(-50, 120))
)

p <- layout(
  p,
  scene = scene_cfg,
  title = list(
    text = "<b>3D Multilayer Causal Graph (created by finalists with Rstudio)</b>",
    font = list(family = "Helvetica", size = 20, color = "black"),
    x = 0.4,
    xanchor = "center",
    y = 0.95
  ),
  font = list(family = "Helvetica"),
  margin = list(l = 0, r = 300, t = 80, b = 150),
  legend = list(
    orientation = "h",
    x = 0.5,
    y = -0.2,
    xanchor = "center",
    yanchor = "top",
    traceorder = "normal"
  )
)

p <- config(p, displayModeBar = TRUE)
print(p)

# Save using Kaleido backend
#save_plot_png_kaleido <- function(plot_obj, filename = "new_withoutdwt_1_3d.png",
#                                  width = 1800, height = 900, scale = 3) {
#  save_image(plot_obj, file = filename, width = width, height = height, scale = scale)
#  message("High resolution PNG saved to: ", normalizePath(filename))
#}


#save_plot_png_kaleido(p, "new_withoutdwt_1_3d.png")