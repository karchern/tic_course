# A zoom into the raw colony measurements behind one gene.
#
# Shows every colony of that gene on every plate copy, under IPTG alone (gene
# on, toxin off) and under Ara + IPTG (both on), with a few GFP colonies from
# the same plates as a reference. A line joins each physical colony position
# across the two conditions, so the drop caused by the toxin is visible.
#
#   source("inspect_gene.R")
#   inspect_gene("50aa_1_55")
library(tidyverse)

inspect_gene <- function(gene, colonies = NULL, n_gfp = 4) {
  if (is.null(colonies)) colonies <- read_tsv("data/colonies.tsv", show_col_types = FALSE)
  if (!gene %in% colonies$genename) stop(gene, " is not in the data")

  # The same few GFP positions on every plate, spread across the plate, so the
  # reference is the same colonies in both conditions.
  gfp_positions <- colonies %>%
    filter(genename == "gfp") %>%
    distinct(row, column) %>%
    arrange(row, column) %>%
    slice(round(seq(1, n(), length.out = n_gfp)))

  shown <- colonies %>%
    filter(genename == gene |
           (genename == "gfp" & paste(row, column) %in% paste(gfp_positions$row, gfp_positions$column))) %>%
    mutate(what     = if_else(genename == "gfp", "GFP control", gene),
           what     = factor(what, levels = c(gene, "GFP control")),
           position = paste0(LETTERS[row], column),
           condition = factor(induced, levels = c("IPTG", "AraIPTG"),
                              labels = c("IPTG\ngene only", "Ara + IPTG\ngene + toxin")))

  # the numbers behind the picture
  shown %>%
    group_by(plate_rep, what, condition) %>%
    summarise(colonies = n(), opacity = paste(sort(opacity), collapse = ", "), .groups = "drop") %>%
    arrange(plate_rep, what, condition) %>%
    print(n = Inf, width = Inf)

  ggplot(shown, aes(condition, opacity, colour = what, group = position)) +
    geom_line(alpha = 0.4) +
    geom_point(size = 3, alpha = 0.85) +
    facet_wrap(~ paste("plate copy", plate_rep)) +
    scale_y_continuous(limits = c(0, NA), labels = scales::label_number(scale = 1e-3, suffix = "k")) +
    scale_colour_manual(values = setNames(c("firebrick", "grey55"), c(gene, "GFP control"))) +
    labs(x = NULL, y = "Iris opacity", colour = NULL,
         title = paste("Raw colonies:", gene),
         subtitle = paste("each line is one colony position;", n_gfp, "GFP colonies per plate for reference")) +
    theme_bw(base_size = 12) +
    theme(legend.position = "top", panel.grid.minor = element_blank())
}
