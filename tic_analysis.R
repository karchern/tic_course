# Toxin Inhibition Conjugation (TIC) screen: from colony measurements to a
# volcano plot of RcaT blockers.
#
# Each colony is one E. coli strain carrying two plasmids: the RcaT toxin
# (arabinose-inducible) and a designed blocker from the library (IPTG-inducible).
# A blocker works if its colonies still grow when both are induced.
#
# Run from the repo root:  Rscript preprocess.R && Rscript tic_analysis.R
library(tidyverse)

out_dir  <- "output"
p_cutoff <- 0.01


# 1. Load the preprocessed colonies ------------------------------------------
# preprocess.R built this table: every .iris file read and joined to the plate
# map, so each colony carries its strain and gene name, and the plate name
# decoded into plasmid, retron, induction and plate copy. It kept the pir116
# plasmid, RcaT-Sen2, low IPTG, plate copies A_1 and B_1, and the two
# conditions compared below: IPTG alone and Ara + IPTG.
colonies <- read_tsv("data/colonies.tsv", show_col_types = FALSE)


# 3. Correct the plate edges -------------------------------------------------
colonies <- colonies %>%
  group_by(plate_rep, induced) %>%
  mutate(interior_median_opacity = median(opacity[row %in% 2:15 & column %in% 2:23])) %>%
  mutate(exterior_median_opacity = median(opacity[!(row %in% 2:15 & column %in% 2:23)])) %>%
  mutate(correction_factor = interior_median_opacity / exterior_median_opacity) %>%
  mutate(measurement_location = ifelse(
    row %in% 2:15 & column %in% 2:23, "plate_interior", "plate_edge"
  )) %>%
  mutate(opacity_edge_corrected = if_else(measurement_location == "plate_edge", opacity * correction_factor, opacity))

# 4. Normalise to the GFP controls on the same plate -------------------------
fitness <- colonies %>%
  group_by(plate_rep, induced) %>%
  mutate(normalized_fitness = opacity_edge_corrected / mean(opacity_edge_corrected[genename == "gfp"])) %>%
  ungroup()

# 5. Effect size -------------------------------------------------------------
# How much better does a strain grow with toxin AND blocker than with the toxin
# alone? One ratio per clone per plate copy, so 6 values for each blocker.

ratios <- fitness %>%
  select(
    plate_rep, induced, strain, plate_96, genename, normalized_fitness
  ) %>%
  pivot_wider(names_from = induced, values_from = normalized_fitness) %>%
  mutate(ratio = AraIPTG / IPTG)

# Check for undefined ratios (IPTG values very low -> pinning errors)
ratios <-  ratios %>%
  filter(is.finite(ratio))

gfp_ratios <- ratios %>% filter(genename == 'gfp')
library_gene_ratios <- ratios %>% filter(genename != 'gfp')

# We will run a naive t-test between ratios of library genes against the ratio of GFP values. 
# Since we already normalized by gfp (per plate), we can do this over all data points

library_gene_ratios_listed <- library_gene_ratios %>%
  group_by(genename) %>%
  summarize(library_gene_ratio = list(ratio))

gfp_ratios_listed <- gfp_ratios %>%
  group_by(genename) %>%
  summarize(gfp_ratio = list(ratio)) %>%
  select(-genename)

results <- library_gene_ratios_listed %>%
  ungroup() %>%
  mutate(
    t.test.one.sided = 
      map(library_gene_ratio, function(lib_gene_rat) {
#        print(c(length(lib_gene), length(gfp_gene)))
        if (length(lib_gene_rat) > 1) {
          t.test(
            lib_gene_rat, mu = 1
          )          
        } else {
          return(NA)
        }
      }
    )
  ) %>%
  mutate(t.test.pvalue = map_dbl(t.test.one.sided, \(x) x$p.value)) %>%
  mutate(t.test.pvalue.corrected = p.adjust(t.test.pvalue, method = 'BH')) %>%
  mutate(med_ratio = map_dbl(library_gene_ratio, median)) %>%
  arrange(t.test.pvalue.corrected) %>%
  rename(
    effect = med_ratio, p_value = t.test.pvalue.corrected
  )

# controls <- ratios %>%
#   filter(genename == "gfp") %>%
#   group_by(strain) %>%
#   summarise(ratio = mean(ratio), .groups = "drop")


# # 6. Test each blocker against the GFP controls ------------------------------
# # Effect size 1 means "behaves like GFP"; above 1 means it rescues growth.
# compare_with_gfp <- function(ratio, retron) {
#   control <- controls$ratio[controls$retron == retron[1]]
#   tibble(effect  = mean(ratio) / mean(control),
#          p_value = t.test(ratio, control)$p.value,
#          n = length(ratio), n_control = length(control))
# }

# results <- ratios %>%
#   filter(genename != "gfp") %>%
#   group_by(retron, genename) %>%
#   reframe(compare_with_gfp(ratio, retron)) %>%
#   left_join(cutoffs, by = "retron") %>%
#   mutate(hit = p_value < p_cutoff & effect > cutoff)


# 7. Volcano plot ------------------------------------------------------------
volcano <- ggplot(results, aes(effect, -log10(p_value))) +
  #geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed", colour = "grey50") +
  #geom_vline(data = cutoffs, aes(xintercept = cutoff), linetype = "dashed", colour = "grey50") +
  #geom_point(aes(colour = hit), alpha = 0.7) +
  geom_point(alpha = 0.7)
#  geom_text(aes(label = if_else(hit, genename, "")), size = 2.5, vjust = -0.9, check_overlap = TRUE) +
#  geom_text(aes(label = if_else(genename, "")), size = 2.5, vjust = -0.9, check_overlap = TRUE) +
#  facet_wrap(~ retron, scales = "free_x") +
  # scale_x_continuous(expand = expansion(mult = 0.1)) +   # room for the labels
  # scale_colour_manual(values = c(`FALSE` = "grey60", `TRUE` = "firebrick"), guide = "none") +
  # labs(x = "effect size:  (toxin + blocker) / toxin,  relative to GFP",
  #      y = "-log10(p)",
  #      title = "Designed blockers of RcaT toxicity",
  #      subtitle = "dashed lines: the significance and effect-size cutoffs") +
  # theme_bw()

dir.create(out_dir, showWarnings = FALSE)
ggsave(file.path(out_dir, "volcano.png"), volcano, width = 9, height = 4.5, dpi = 150, bg = "white")
results %>% select(genename, effect, p_value, t.test.pvalue) %>%
  write_tsv(file.path(out_dir, "hits.tsv"))

results %>% select(genename, effect, p_value) %>% slice_head(n = 10) %>% print()
