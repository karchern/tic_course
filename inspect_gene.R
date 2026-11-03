# A zoom into the raw colony measurements behind one gene.
#
# Left:  what each colony does when the toxin is switched on. A line joins each
#        physical colony position across the two conditions.
# Right: the same colonies as they sit on the plate. Every well is pinned as a
#        2x2 quartet (two clones, each pinned twice), with its opacity printed.
#
# A GFP quartet from the same plates is shown alongside as the reference: GFP
# does not block the toxin, so it is what "no effect" looks like.
#
#   source("inspect_gene.R")
#   inspect_gene("50aa_1_55")
library(tidyverse)
library(patchwork)

inspect_gene <- function(gene, colonies = NULL) {
  if (is.null(colonies)) colonies <- read_tsv("data/colonies.tsv", show_col_types = FALSE)
  if (!gene %in% colonies$genename) stop(gene, " is not in the data")

  # One GFP well, near the middle of the GFP block, as the reference quartet.
  gfp_well <- colonies %>%
    filter(genename == "gfp") %>%
    distinct(well_row = (row + 1) %/% 2, well_col = (column + 1) %/% 2) %>%
    arrange(well_row, well_col) %>%
    slice(ceiling(n() / 2))

  shown <- colonies %>%
    filter(genename == gene |
           (genename == "gfp" &
            (row + 1) %/% 2 == gfp_well$well_row &
            (column + 1) %/% 2 == gfp_well$well_col)) %>%
    mutate(what      = factor(if_else(genename == "gfp", "GFP control", gene),
                              levels = c(gene, "GFP control")),
           position  = paste0(LETTERS[row], column),
           condition = factor(induced, levels = c("IPTG", "AraIPTG"),
                              labels = c("IPTG\ngene only", "Ara + IPTG\ngene + toxin")),
           # where each colony sits inside its own 2x2 quartet
           quartet_row = row - 2 * ((row + 1) %/% 2) + 2,
           quartet_col = column - 2 * ((column + 1) %/% 2) + 2)

  strain_colours <- setNames(c("firebrick", "grey40"), c(gene, "GFP control"))
  # text on a tile is dark on the pale half of the scale, pale on the dark half
  mid_opacity <- max(shown$opacity) / 2

  summary_panel <- ggplot(shown, aes(condition, opacity, colour = what, group = position)) +
    geom_line(alpha = 0.4) +
    geom_point(size = 3, alpha = 0.85) +
    facet_wrap(~ paste("plate copy", plate_rep), ncol = 1) +
    scale_y_continuous(limits = c(0, NA), labels = scales::label_number(scale = 1e-3, suffix = "k")) +
    scale_colour_manual(values = strain_colours) +
    labs(x = NULL, y = "Iris opacity", colour = NULL) +
    theme_bw(base_size = 12) +
    theme(legend.position = "top", panel.grid.minor = element_blank())

  quartet_panel <- ggplot(shown, aes(quartet_col, quartet_row)) +
    geom_tile(aes(fill = opacity), colour = "white", linewidth = 1.5) +
    geom_text(aes(label = round(opacity / 1000), colour = opacity > mid_opacity),
              size = 3.4, show.legend = FALSE) +
    facet_grid(paste("plate copy", plate_rep) ~ what + condition, switch = "y") +
    scale_y_reverse() +
    scale_fill_gradient(low = "black", high = "white", guide = "none") +
    scale_colour_manual(values = c(`TRUE` = "grey10", `FALSE` = "grey95")) +
    coord_fixed() +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 11) +
    theme(axis.text = element_blank(), panel.grid = element_blank())

  summary_panel + quartet_panel +
    plot_layout(widths = c(1, 1.3)) +
    plot_annotation(title = paste("Raw colonies:", gene),
                    subtitle = "every well is pinned as a 2x2 quartet; GFP shown as the no-effect reference",
                    caption = "squares are the colonies as they sit on the plate, labelled with opacity in thousands",
                    theme = theme(plot.title = element_text(size = 15, face = "bold")))
}
