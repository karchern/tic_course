# Figures for the intro slides: Sen2 plates from plate replicate 3A_1. The double-induction plate
# (Ara + 0.1 mM IPTG) is the one shown in the paper (Fig. 7A).
# Run from the repo root: Rscript slides/make_figures.R
library(tidyverse)

plate_map <- read_csv("data/AlphaDesign_384wp.identifiers.csv", show_col_types = FALSE) %>%
  select(row = row.384, column = column.384, genename)

read_plate <- function(file) {
  read_tsv(file, skip = 6, show_col_types = FALSE) %>%
    left_join(plate_map, by = c("row", "column"))
}

araiptg <- read_plate("data/iris/SpecGentAraIPTG-10050201-3A_1.jpg.iris")
iptg    <- read_plate("data/iris/SpecGentIPTG-100501-3A_1.jpg.iris")
none    <- read_plate("data/iris/SpecGent-1005-3A_1.jpg.iris")

# box around the 4 colonies (2 clones x 2 technical replicates) of one gene
gene_box <- function(gene) {
  plate_map %>% filter(genename == gene) %>%
    summarise(xmin = min(column) - 0.5, xmax = max(column) + 0.5,
              ymin = min(row) - 0.5,    ymax = max(row) + 0.5, genename = gene)
}

box_colours <- c("50aa_1_55" = "red", "50aa_1_14" = "#1f77b4", "50aa_2_5" = "#ff7f0e")
boxes <- map_dfr(names(box_colours), gene_box)

plate_heatmap <- function(plate, genes = character(0), title = NULL) {
  ggplot(plate, aes(column, row, fill = opacity)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    geom_rect(data = filter(boxes, genename %in% genes),
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, colour = genename),
              inherit.aes = FALSE, fill = NA, linewidth = 1.2) +
    scale_fill_gradient(low = "black", high = "white", name = "Iris opacity") +
    scale_colour_manual(values = box_colours, name = NULL) +
    scale_y_reverse(breaks = 1:16, labels = LETTERS[1:16], expand = c(0, 0)) +
    scale_x_continuous(breaks = 1:24, position = "top", expand = c(0, 0)) +
    coord_fixed() +
    labs(x = NULL, y = NULL, title = title) +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank())
}

save_fig <- function(p, name) ggsave(file.path("slides/figures", name), p, width = 9, height = 5.4, dpi = 150, bg = "white")

save_fig(plate_heatmap(araiptg), "plate_heatmap.png")
save_fig(plate_heatmap(araiptg, "50aa_1_55"), "plate_heatmap_hit.png")
save_fig(plate_heatmap(araiptg, c("50aa_1_55", "50aa_1_14")), "plate_heatmap_hit_toxic.png")
save_fig(plate_heatmap(araiptg, "50aa_1_14", "Ara + IPTG (toxin + library gene)"), "plate_heatmap_toxic_araiptg.png")
save_fig(plate_heatmap(iptg, "50aa_1_14", "IPTG only (library gene)"), "plate_heatmap_toxic_iptg.png")
# dark spots: colonies with opacity 0 on IPTG plates are present on the uninduced plate
save_fig(plate_heatmap(araiptg, "50aa_2_5", "Ara + IPTG (toxin + library gene)"), "plate_heatmap_dark_araiptg.png")
save_fig(plate_heatmap(none, "50aa_2_5", "No induction"), "plate_heatmap_dark_none.png")

# values for the slide text: opacity of each highlighted colony relative to its plate's GFP mean
bind_rows(`Ara + IPTG` = araiptg, IPTG = iptg, none = none, .id = "plate") %>%
  group_by(plate) %>%
  mutate(rel_to_gfp = round(opacity / mean(opacity[genename == "gfp"]), 2)) %>%
  filter(genename %in% names(box_colours)) %>%
  select(plate, genename, row, column, opacity, rel_to_gfp) %>%
  print(n = Inf)
