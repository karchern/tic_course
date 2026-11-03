# Toxin Inhibition Conjugation (TIC) screen: from Iris colony measurements to a
# volcano plot of RcaT blockers.
#
# Each colony is one E. coli strain carrying two plasmids: the RcaT toxin
# (arabinose-inducible) and a designed blocker from the library (IPTG-inducible).
# A blocker works if its colonies still grow when both are induced.
#
# Method follows Jendrusch et al. 2025, Mol Syst Biol, Fig. 7B,C.
# Run from the repo root:  Rscript tic_analysis.R
library(tidyverse)

iris_dir <- "data/iris"
map_file <- "data/AlphaDesign_384wp.identifiers.csv"
out_dir  <- "output"

p_cutoff <- 0.01
cutoffs  <- tibble(retron = c("Sen2", "Eco9"), cutoff = c(1.1, 1.3))


# 1. Load --------------------------------------------------------------------
# The map says which strain sits at each of the 384 positions. One 96-well
# source plate is pinned into 384 format, so each well becomes a 2x2 block:
# two clones, each pinned twice. "strain" identifies one clone of one well.
plate_map <- read_csv(map_file, show_col_types = FALSE) %>%
  select(row = row.384, column = column.384, strain, genename)

read_iris <- function(path) {
  read_tsv(path, skip = 6, show_col_types = FALSE) %>%
    mutate(plate = str_remove(basename(path), "\\.jpg\\.iris$"))
}

colonies <- list.files(iris_dir, "\\.iris$", full.names = TRUE) %>%
  map(read_iris) %>%
  list_rbind() %>%
  left_join(plate_map, by = c("row", "column"))

# The plate name carries the design, e.g. SpecGentAraIPTG-10050201-3A_1:
# antibiotics, inducer concentrations, then which copy of the plate it is.
colonies <- colonies %>%
  separate_wider_delim(plate, "-", names = c("drugs", "inducers", "copy")) %>%
  mutate(
    block     = as.integer(str_sub(copy, 1, 1)),
    plate_rep = str_sub(copy, 2),
    plasmid   = if_else(block %in% c(1, 2), "p-RcaT", "pir116"),
    retron    = if_else(block %in% c(1, 3), "Sen2", "Eco9"),
    induced   = str_remove(str_remove(drugs, "Spec"), "Gent|tet"),
    induced   = if_else(induced == "", "none", induced),
    iptg      = case_when(inducers %in% c("100501", "10050201", "1001001", "100100201") ~ "low",
                          inducers %in% c("10051", "1005021", "100101", "10010021")     ~ "high",
                          .default = "none"))


# 2. Filter ------------------------------------------------------------------
colonies <- colonies %>%
  filter(plasmid == "pir116",   # the high-copy toxin plasmid
         iptg != "high",        # low blocker induction only
         !is.na(genename))      # positions the map leaves empty

# Colonies measured as opacity 0 are KEPT. The paper's methods call them pinning
# misses and drop them, but they sit at the same few genes on every plate and
# only ever on IPTG plates, so they are strains killed by their own blocker
# rather than misses. Keeping them also reproduces the published numbers more
# closely. To drop them instead, uncomment:
# colonies <- filter(colonies, opacity > 0)


# 3. Correct the plate edges -------------------------------------------------
# Colonies on the outer row and column have no neighbours on one side, so they
# grow bigger. Each outer line is scaled to the median of the plate's interior.
colonies <- colonies %>%
  group_by(retron, plate_rep, induced) %>%
  mutate(interior = median(opacity[row %in% 2:15 & column %in% 2:23])) %>%
  group_by(retron, plate_rep, induced, row) %>%
  mutate(opacity = if_else(row %in% c(1, 16), opacity * interior / median(opacity), opacity)) %>%
  group_by(retron, plate_rep, induced, column) %>%
  mutate(opacity = if_else(column %in% c(1, 24), opacity * interior / median(opacity), opacity)) %>%
  ungroup()


# 4. Normalise to the GFP controls on the same plate -------------------------
# GFP does not block the toxin, so its colonies show what "no effect" looks
# like on that plate. Dividing by their mean puts every plate on one footing.
fitness <- colonies %>%
  group_by(retron, plate_rep, induced) %>%
  mutate(fitness = opacity / mean(opacity[genename == "gfp"])) %>%
  ungroup() %>%
  # average the two technical replicates of each clone
  group_by(retron, plate_rep, induced, strain, genename) %>%
  summarise(fitness = mean(fitness), .groups = "drop")


# 5. Effect size -------------------------------------------------------------
# How much better does a strain grow with toxin AND blocker than with the toxin
# alone? One ratio per clone per plate copy, so 6 values for each blocker.
ratios <- fitness %>%
  filter(induced %in% c("AraIPTG", "Ara")) %>%
  pivot_wider(names_from = induced, values_from = fitness) %>%
  mutate(ratio = AraIPTG / Ara) %>%
  filter(is.finite(ratio))

# The paper compares each blocker (n = 6) with 16 GFP values, one per GFP
# strain, so the GFP ratios are averaged over the three plate copies first.
controls <- ratios %>%
  filter(genename == "gfp") %>%
  group_by(retron, strain) %>%
  summarise(ratio = mean(ratio), .groups = "drop")


# 6. Test each blocker against the GFP controls ------------------------------
# Effect size 1 means "behaves like GFP"; above 1 means it rescues growth.
compare_with_gfp <- function(ratio, retron) {
  control <- controls$ratio[controls$retron == retron[1]]
  tibble(effect  = mean(ratio) / mean(control),
         p_value = t.test(ratio, control)$p.value,
         n = length(ratio), n_control = length(control))
}

results <- ratios %>%
  filter(genename != "gfp") %>%
  group_by(retron, genename) %>%
  reframe(compare_with_gfp(ratio, retron)) %>%
  left_join(cutoffs, by = "retron") %>%
  mutate(hit = p_value < p_cutoff & effect > cutoff)


# 7. Volcano plot ------------------------------------------------------------
volcano <- ggplot(results, aes(effect, -log10(p_value))) +
  geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed", colour = "grey50") +
  geom_vline(data = cutoffs, aes(xintercept = cutoff), linetype = "dashed", colour = "grey50") +
  geom_point(aes(colour = hit), alpha = 0.7) +
  geom_text(aes(label = if_else(hit, genename, "")), size = 2.5, vjust = -0.9, check_overlap = TRUE) +
  facet_wrap(~ retron, scales = "free_x") +
  scale_x_continuous(expand = expansion(mult = 0.1)) +   # room for the labels
  scale_colour_manual(values = c(`FALSE` = "grey60", `TRUE` = "firebrick"), guide = "none") +
  labs(x = "effect size:  (toxin + blocker) / toxin,  relative to GFP",
       y = "-log10(p)",
       title = "Designed blockers of RcaT toxicity",
       subtitle = "dashed lines: the significance and effect-size cutoffs") +
  theme_bw()

dir.create(out_dir, showWarnings = FALSE)
ggsave(file.path(out_dir, "volcano.png"), volcano, width = 9, height = 4.5, dpi = 150, bg = "white")
write_tsv(results, file.path(out_dir, "hits.tsv"))

results %>% count(retron, hit) %>% print()
