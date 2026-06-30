# --- Celltype–Celltype --------------------------------------------------------
library(ggplot2)
library(dplyr)
library(grid)
library(cowplot)

setwd("/Users/oliviaxu/Desktop/wavelet_project/3band_2reg/wavelet_3band_2reg_view1_matrices")

edges <- read.csv("3band_2reg_view1_4_5.csv", stringsAsFactors = FALSE, check.names = FALSE) %>%
  transmute(
    source = trimws(as.character(Source)),
    target = trimws(as.character(Target)),
    weight = suppressWarnings(as.numeric(Signed_weight))
  )

cc_edges <- edges %>%
  filter(grepl("^Celltype_", source, ignore.case = TRUE) &
           grepl("^Celltype_", target, ignore.case = TRUE)) %>%
  mutate(
    source = trimws(gsub("\\s+", " ", sub("^Celltype_\\s*", "", source, ignore.case = TRUE))),
    target = trimws(gsub("\\s+", " ", sub("^Celltype_\\s*", "", target, ignore.case = TRUE)))
  ) %>%
  group_by(source, target) %>%
  summarise(weight = mean(weight, na.rm = TRUE), .groups = "drop")

cell_names <- sort(unique(c(cc_edges$source, cc_edges$target)))
n_cells <- length(cell_names)
circle_radius <- 5

if (n_cells == 0) {
  layout_df <- data.frame(name=character(0), x=numeric(0), y=numeric(0),
                          full_label=character(0), color=character(0),
                          stringsAsFactors=FALSE)
} else {
  theta <- if (n_cells == 1) 0 else seq(0, 2*pi, length.out = n_cells + 1)[-(n_cells+1)]
  layout_df <- data.frame(
    name       = cell_names,
    x          = circle_radius * cos(theta),
    y          = circle_radius * sin(theta),
    full_label = cell_names,
    color      = "#FFD580",
    stringsAsFactors = FALSE
  )
}

extra <- 0.5

edge_df <- cc_edges %>%
  inner_join(layout_df %>% select(name, x_source = x, y_source = y),
             by = c("source" = "name")) %>%
  inner_join(layout_df %>% select(name, x_target = x, y_target = y),
             by = c("target" = "name")) %>%
  mutate(
    dx = x_target - x_source, dy = y_target - y_source,
    L  = sqrt(dx^2 + dy^2), L = ifelse(L==0, 1e-9, L),
    ux = dx / L, uy = dy / L,
    x_start = x_source + ux * extra,
    y_start = y_source + uy * extra,
    x_end   = x_target - ux * extra,
    y_end   = y_target - uy * extra,
    edge_color = ifelse(weight > 0, "#009851", "#C52C32")
  )

min_width <- 0.25; max_width <- 1.25
w <- abs(edge_df$weight)
edge_df$edge_size <- if (length(w)>1 && diff(range(w))>0) {
  min_width + (w - min(w))*(max_width-min_width)/(max(w)-min(w))
} else rep((min_width+max_width)/2, length(w))

node_strength_df <- bind_rows(
  edge_df %>% transmute(name = source, w = abs(weight)),
  edge_df %>% transmute(name = target, w = abs(weight))
) %>%
  group_by(name) %>%
  summarise(strength = sum(w, na.rm = TRUE), .groups = "drop")

layout_df <- layout_df %>%
  left_join(node_strength_df, by = "name") %>%
  mutate(strength = ifelse(is.na(strength), 0, strength))

label_pad <- 0.4
layout_df <- layout_df %>%
  mutate(
    theta = atan2(y, x),
    x_label = x + cos(theta) * label_pad,
    y_label = y + sin(theta) * label_pad,
    hjust = ifelse(cos(theta) >  0.2, 0,
                   ifelse(cos(theta) < -0.2, 1, 0.5)),
    vjust = ifelse(sin(theta) >  0.2, 0,
                   ifelse(sin(theta) < -0.2, 1, 0.5))
  )

plot_cell_cell <- ggplot() +
  geom_segment(
    data=edge_df,
    aes(x=x_start, y=y_start, xend=x_end, yend=y_end,
        color=edge_color, linewidth=edge_size),
    arrow=arrow(length=unit(0.25,"cm"), type="closed", angle=20),
    show.legend=FALSE, lineend="round"
  ) +
  geom_point(
    data=layout_df,
    aes(x=x,y=y,size=strength),
    color="#FFD580", fill="#FFD580", shape=21, stroke=0.6
  ) +
  scale_size_continuous(range=c(8,18), guide="none") +
  geom_text(
    data=layout_df,
    aes(x=x_label, y=y_label, label=full_label, hjust=hjust, vjust=vjust),
    size=6, fontface="bold"
  ) +
  scale_color_identity() +
  scale_linewidth_identity() +
  scale_x_continuous(limits=c(-circle_radius-0.8, circle_radius+0.8),
                     expand=expansion(mult=0)) +
  scale_y_continuous(limits=c(-circle_radius-0.8, circle_radius+0.8),
                     expand=expansion(mult=0)) +
  coord_equal(clip="off") +
  theme_void() +
  ggtitle("3-Band and 4-Band Acinar") +
  theme(
    plot.title = element_text(size = 40, face = "bold", hjust = 0.5),
    plot.title.position = "plot",
    plot.margin=margin(t=20,r=20,b=20,l=20)
  )

ggsave("cell_cell.png", plot = plot_cell_cell,
       width = 12, height = 8, dpi = 600, bg = "white")





































y0 <- 0.8; dy <- 0.02; ys <- function(i) y0 - dy*(i-1)

sign_df <- data.frame(
  x_start=c(0,0), x_end=c(1,1),
  y=c(ys(1),ys(2)),
  color=c("#009851","#C52C32"),
  label=c("positive causal link","negative causal link")
)

abs_w <- abs(edge_df$weight)
w_q <- round(as.numeric(quantile(abs_w, probs=c(1,0.5,0), na.rm=TRUE)),3)
names(w_q) <- c("max","med","min")
lw_min <- min(edge_df$edge_size,na.rm=TRUE)
lw_max <- max(edge_df$edge_size,na.rm=TRUE)
lw_mid <- (lw_min+lw_max)/2

w_df <- data.frame(
  x_start=0, x_end=1,
  y=c(ys(3),ys(4),ys(5)),
  size=c(lw_max,lw_mid,lw_min),
  label=paste0("| weight | = ", c(w_q["max"], w_q["med"], w_q["min"]))
)

legend_plot_cc <- ggplot() +
  geom_segment(data=sign_df,
               aes(x=x_start,y=y,xend=x_end,yend=y,color=color),
               arrow=arrow(length=unit(0.20,"cm"), type="closed", angle=20),
               linewidth=1.2, lineend="round") +
  geom_text(data=sign_df, aes(x=1.15,y=y,label=label),
            hjust=0, size=6, fontface="bold") +
  geom_segment(data=w_df,
               aes(x=x_start,y=y,xend=x_end,yend=y),
               color="grey20", linewidth=w_df$size, lineend="round") +
  geom_text(data=w_df, aes(x=1.15,y=y,label=label),
            hjust=0, size=6, fontface="bold") +
  scale_color_identity() +
  coord_cartesian(xlim=c(0,3), ylim=c(ys(5)-0.01,y0+0.01),
                  expand=FALSE, clip="off") +
  theme_void() +
  ggtitle("Legend") +
  theme(
    plot.title=element_text(size=30, face="bold", hjust=0.5),
    plot.margin=margin(2,6,2,6),
    plot.background=element_rect(fill="white", color=NA)
  )

base_x <- -0.07; base_w <- 0.90
legend_w <- 0.23; legend_h <- 0.5
legend_x <- 0.75; legend_y <- (1 - legend_h)/2

final_cell_cell <- ggdraw() +
  draw_label("3-Band and 4-Band Acinar", x=0.5, y=0.98,
             hjust=0.5, vjust=1, fontface="bold", size=40) +
  draw_plot(plot_cell_cell, x=base_x, y=0, width=base_w, height=1) +
  draw_plot(legend_plot_cc, x=legend_x, y=legend_y,
            width=legend_w, height=legend_h)
final_cell_cell

ggsave(
  filename = "/Users/oliviaxu/Desktop/cell_cell_Acinar.png",
  plot = final_cell_cell,
  width = 29,
  height = 20,
  dpi = 600,
  bg = "white"
)
