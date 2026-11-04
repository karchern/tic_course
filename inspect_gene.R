# A zoom into the raw data behind one gene.
#
# Left:  what each colony does when the toxin is switched on. A line joins each
#        physical colony position across the two conditions. These are the Iris
#        opacity numbers the analysis actually uses.
# Right: the colonies themselves, cut straight out of the plate photographs.
#        Every well is pinned as a 2x2 quartet (two clones, each pinned twice).
#
# A GFP quartet from the same plates is shown alongside as the reference: GFP
# does not block the toxin, so it is what "no effect" looks like.
#
# Only plates we have a photograph of can appear on the right. Photos go in
# data/photos, named after the plate like the .iris files; see plate_photos.R.
#
#   source("inspect_gene.R")
#   inspect_gene("50aa_1_55")
library(tidyverse)
library(patchwork)
source("plate_photos.R")

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
           plate     = paste(drugs, inducers, copy, sep = "-"),
           well_row  = (row + 1) %/% 2,
           well_col  = (column + 1) %/% 2)

  strain_colours <- setNames(c("firebrick", "grey40"), c(gene, "GFP control"))

  summary_panel <- ggplot(shown, aes(condition, opacity, colour = what, group = position)) +
    geom_line(alpha = 0.4) +
    geom_point(size = 3, alpha = 0.85) +
    facet_wrap(~ paste("plate copy", plate_rep), ncol = 1) +
    scale_y_continuous(limits = c(0, NA), labels = scales::label_number(scale = 1e-3, suffix = "k")) +
    scale_colour_manual(values = strain_colours) +
    labs(x = NULL, y = "Iris opacity", colour = NULL,
         subtitle = "what Iris measured") +
    theme_bw(base_size = 12) +
    theme(legend.position = "top", panel.grid.minor = element_blank())

  # The raw colonies, for the plates we have a photograph of.
  have <- intersect(unique(shown$plate), plates_with_photos())
  missing <- setdiff(unique(shown$plate), plates_with_photos())
  if (length(have) == 0)
    stop("no plate photographs found in ", photo_dir, " for the plates in the data")

  pixels <- shown %>%
    filter(plate %in% have) %>%
    distinct(plate, plate_rep, what, condition, well_row, well_col) %>%
    rowwise() %>%
    reframe(across(everything()), quartet_pixels(plate, well_row, well_col))

  photo_panel <- ggplot(pixels, aes(x, y, fill = grey)) +
    geom_raster() +
    facet_grid(paste("plate copy", plate_rep) ~ what + condition, switch = "y") +
    scale_y_reverse() +
    scale_fill_gradient(low = "black", high = "white", guide = "none") +
    coord_fixed() +
    labs(x = NULL, y = NULL, subtitle = "the colonies themselves") +
    theme_minimal(base_size = 11) +
    theme(axis.text = element_blank(), panel.grid = element_blank(),
          panel.spacing = unit(6, "pt"))

  note <- if (length(missing) == 0) {
    "right: every quartet cut out of its plate photograph"
  } else {
    message("no photograph for: ", paste(missing, collapse = ", "),
            "\n  (drop one in ", photo_dir, " named after the plate to include it)")
    sprintf("right: real colonies, for the %d of %d plates here we have a photograph of",
            length(have), length(have) + length(missing))
  }

  summary_panel + photo_panel +
    plot_layout(widths = c(1, 1.3)) +
    plot_annotation(title = paste("Raw data:", gene),
                    subtitle = "every well is pinned as a 2x2 quartet; GFP shown as the no-effect reference",
                    caption = note,
                    theme = theme(plot.title = element_text(size = 15, face = "bold")))
}
