# Figures for the intro slides: one Sen2 double-induction plate (Ara + 0.1 mM IPTG),
# the same plate shown in the paper (Fig. 7A).
# Run from the repo root: Rscript slides/make_figures.R
library(tidyverse)

plate_file <- "data/iris/SpecGentAraIPTG-10050201-3A_1.jpg.iris"

plate_map <- read_csv("data/AlphaDesign_384wp.identifiers.csv", show_col_types = FALSE) %>%
  select(row = row.384, column = column.384, genename)

plate <- read_tsv(plate_file, skip = 6, show_col_types = FALSE) %>%
  left_join(plate_map, by = c("row", "column"))

# box around the 4 colonies (2 clones x 2 technical replicates) of one gene
gene_box <- function(gene) {
  plate %>% filter(genename == gene) %>%
    summarise(xmin = min(column) - 0.5, xmax = max(column) + 0.5,
              ymin = min(row) - 0.5,    ymax = max(row) + 0.5, genename = gene)
}

plate_plot <- ggplot(plate, aes(column, row, fill = opacity)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  scale_fill_gradient(low = "black", high = "white", name = "Iris opacity") +
  scale_y_reverse(breaks = 1:16, labels = LETTERS[1:16], expand = c(0, 0)) +
  scale_x_continuous(breaks = 1:24, position = "top", expand = c(0, 0)) +
  coord_fixed() +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank())

ggsave("slides/figures/plate_heatmap.png", plate_plot, width = 9, height = 5.4, dpi = 150, bg = "white")

boxes <- bind_rows(gene_box("50aa_1_55"), gene_box("50aa_2_1"))
box_colours <- c("50aa_1_55" = "red", "50aa_2_1" = "orange")

plate_plot_boxes <- function(genes) {
  plate_plot +
    geom_rect(data = filter(boxes, genename %in% genes),
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, colour = genename),
              inherit.aes = FALSE, fill = NA, linewidth = 1.2) +
    scale_colour_manual(values = box_colours, name = NULL)
}

ggsave("slides/figures/plate_heatmap_hit.png", plate_plot_boxes("50aa_1_55"),
       width = 9, height = 5.4, dpi = 150, bg = "white")
ggsave("slides/figures/plate_heatmap_hit_mixed.png", plate_plot_boxes(c("50aa_1_55", "50aa_2_1")),
       width = 9, height = 5.4, dpi = 150, bg = "white")

# values for the slide text: opacity of each highlighted colony relative to the plate's GFP mean
gfp_mean <- mean(plate$opacity[plate$genename == "gfp"])
plate %>% filter(genename %in% names(box_colours)) %>%
  mutate(rel_to_gfp = round(opacity / gfp_mean, 2)) %>%
  select(genename, row, column, opacity, rel_to_gfp) %>%
  print()
