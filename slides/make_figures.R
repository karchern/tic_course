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

# Highlight colours avoid the red of the tetrad grid and the green of the GFP
# block below, so a highlighted strain never reads as "just another tetrad".
box_colours <- c("50aa_1_55" = "#0072B2", "50aa_1_14" = "#7B3FA0", "50aa_2_5" = "#E69F00")
boxes <- map_dfr(names(box_colours), gene_box) %>%
  # keep the legend in reveal order, so an entry never shifts when a later
  # figure adds another gene below it
  mutate(genename = factor(genename, levels = names(box_colours)))

# The same tetrad structure drawn on the plate photo: every 96-well source well
# becomes one 2x2 block in 384 format. The GFP wells form a single run, so they
# get one box rather than eight (they are averaged into one per-plate baseline).
tetrads <- plate_map %>%
  group_by(tx = ceiling(column / 2), ty = ceiling(row / 2)) %>%
  summarise(gfp = all(genename == "gfp"),
            xmin = min(column) - 0.5, xmax = max(column) + 0.5,
            ymin = min(row) - 0.5,    ymax = max(row) + 0.5, .groups = "drop")
library_tetrads <- filter(tetrads, !gfp)
gfp_block <- tetrads %>% filter(gfp) %>%
  summarise(xmin = min(xmin), xmax = max(xmax), ymin = min(ymin), ymax = max(ymax))

plate_heatmap <- function(plate, genes = character(0), title = NULL, show_tetrads = FALSE) {
  rect <- function(d, colour, linewidth)
    geom_rect(data = d, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              inherit.aes = FALSE, fill = NA, colour = colour, linewidth = linewidth)
  ggplot(plate, aes(column, row, fill = opacity)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    (if (show_tetrads) list(rect(library_tetrads, "#cc1a1a", 0.35),
                            rect(gfp_block, "#1a9940", 1.1)) else NULL) +
    geom_rect(data = filter(boxes, genename %in% genes),
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, colour = genename),
              inherit.aes = FALSE, fill = NA, linewidth = 1.6) +
    scale_fill_gradient(low = "black", high = "white", name = "Iris opacity") +
    scale_colour_manual(values = box_colours, breaks = names(box_colours), name = NULL) +
    scale_y_reverse(breaks = 1:16, labels = LETTERS[1:16], expand = c(0, 0)) +
    scale_x_continuous(breaks = 1:24, position = "top", expand = c(0, 0)) +
    coord_fixed() +
    labs(x = NULL, y = NULL, title = title,
         subtitle = if (show_tetrads) "red: one library strain (4 colonies)    green: the 8 GFP control wells" else NULL) +
    theme_minimal(base_size = 12) +
    # Pin the legend block to the top and keep the opacity bar above the gene
    # boxes, so the legend sits in the same place on every figure. Otherwise a
    # figure with two boxed genes has a taller legend, which is centred, and the
    # whole block shifts between slides.
    theme(panel.grid = element_blank(), legend.justification = "top") +
    guides(fill = guide_colourbar(order = 1), colour = guide_legend(order = 2))
}

save_fig <- function(p, name) ggsave(file.path("slides/figures", name), p, width = 9, height = 5.4, dpi = 150, bg = "white")

save_fig(plate_heatmap(araiptg, show_tetrads = TRUE), "plate_heatmap.png")
save_fig(plate_heatmap(araiptg, "50aa_1_55", show_tetrads = TRUE), "plate_heatmap_hit.png")
save_fig(plate_heatmap(araiptg, c("50aa_1_55", "50aa_1_14"), show_tetrads = TRUE), "plate_heatmap_hit_toxic.png")
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

# Alignment check: the heatmap's row/column must be the real position on the plate.
# photo_vs_iris.csv holds the brightness measured from the plate photo at each
# fitted colony position; plotting it beside the Iris values in the same layout
# shows whether the two agree cell by cell.
if (file.exists("slides/figures/photo_vs_iris.csv")) {
  pv <- read_csv("slides/figures/photo_vs_iris.csv", show_col_types = FALSE)
  align <- pv %>%
    pivot_longer(c(photo, opacity), names_to = "source", values_to = "v") %>%
    group_by(source) %>% mutate(v = (v - min(v)) / (max(v) - min(v))) %>% ungroup() %>%
    mutate(source = factor(source, levels = c("photo", "opacity"),
                           labels = c("read off the plate photo", "Iris opacity")))
  p <- ggplot(align, aes(column, row, fill = v)) +
    geom_tile(colour = "white", linewidth = 0.2) +
    facet_wrap(~ source) +
    scale_fill_gradient(low = "black", high = "white", guide = "none") +
    scale_y_reverse(breaks = 1:16, labels = LETTERS[1:16], expand = c(0, 0)) +
    scale_x_continuous(breaks = seq(2, 24, 2), position = "top", expand = c(0, 0)) +
    coord_fixed() +
    labs(x = NULL, y = NULL, title = "Alignment check: two independent readings of the same plate",
         subtitle = sprintf("r = %.3f between photo brightness and Iris opacity at the same row and column",
                            cor(pv$photo, pv$opacity))) +
    theme_minimal(base_size = 11) +
    theme(panel.grid = element_blank())
  ggsave("slides/figures/alignment_check.png", p, width = 11, height = 4.4, dpi = 150, bg = "white")
}
