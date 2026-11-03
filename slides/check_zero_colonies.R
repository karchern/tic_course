# Where do opacity-0 colonies occur on the Sen2 (pir116) plates, and are they
# pinning misses or induced death?
# Pairs every colony across the four conditions within one plate replicate and
# one IPTG level, so each row is the same physical position in all conditions.
# Rscript slides/check_zero_colonies.R
suppressMessages(library(tidyverse))

plate_map <- read_csv("data/AlphaDesign_384wp.identifiers.csv", show_col_types = FALSE) %>%
  select(row = row.384, column = column.384, genename)

sen2 <- list.files("data/iris", "^SpecGent.*-3[AB]_", full.names = TRUE)

plates <- map_dfr(sen2, function(f) {
  id <- sub("\\.jpg\\.iris$", "", basename(f))
  parts <- str_split(id, "-")[[1]]          # e.g. SpecGentAraIPTG / 10050201 / 3A_1
  read_tsv(f, skip = 6, show_col_types = FALSE) %>%
    mutate(replicate = parts[3],
           condition = recode(parts[1], SpecGent = "none", SpecGentAra = "Ara",
                              SpecGentIPTG = "IPTG", SpecGentAraIPTG = "AraIPTG"),
           iptg = case_when(parts[2] %in% c("100501", "10050201") ~ "0.1mM",
                            parts[2] %in% c("10051", "1005021")   ~ "1mM",
                            TRUE ~ "none"))
}) %>% left_join(plate_map, by = c("row", "column"))

cat("== opacity-0 colonies per condition ==\n")
plates %>% group_by(condition) %>%
  summarise(plates = n_distinct(paste(iptg, replicate)), colonies = n(), zeros = sum(opacity == 0)) %>%
  print()

# the uninduced and Ara-only plates have one IPTG level ("none"), so they are
# joined on position and replicate alone
paired <- plates %>%
  filter(condition %in% c("IPTG", "AraIPTG")) %>%
  select(row, column, genename, replicate, iptg, condition, opacity) %>%
  pivot_wider(names_from = condition, values_from = opacity) %>%
  left_join(plates %>% filter(condition == "none") %>% select(row, column, replicate, none = opacity),
            by = c("row", "column", "replicate")) %>%
  left_join(plates %>% filter(condition == "Ara") %>% select(row, column, replicate, Ara = opacity),
            by = c("row", "column", "replicate"))

cat("\n== the question: IPTG-only zero, Ara+IPTG growing? ==\n")
paired %>%
  summarise(`IPTG=0 and AraIPTG>0` = sum(IPTG == 0 & AraIPTG > 0),
            `AraIPTG=0 and IPTG>0` = sum(AraIPTG == 0 & IPTG > 0),
            `both zero`            = sum(IPTG == 0 & AraIPTG == 0),
            `pairs compared`       = n()) %>%
  print()

cat("\n== every position that is zero in at least one induced condition ==\n")
paired %>%
  filter(IPTG == 0 | AraIPTG == 0) %>%
  arrange(genename, replicate, iptg) %>%
  select(genename, row, column, replicate, iptg, none, Ara, IPTG, AraIPTG) %>%
  print(n = Inf)

# The two genes that account for most of the zeros, next to GFP for reference:
# every colony of theirs, in all four conditions, across all three plate replicates.
shown <- c("50aa_3_46", "50aa_2_5", "gfp")
fig <- plates %>%
  filter(genename %in% shown) %>%
  mutate(genename  = factor(genename, levels = shown),
         condition = factor(condition, levels = c("none", "Ara", "IPTG", "AraIPTG"),
                            labels = c("no\ninduction", "Ara\n(toxin)", "IPTG\n(library gene)", "Ara + IPTG\n(both)"))) %>%
  ggplot(aes(condition, opacity)) +
  geom_hline(yintercept = 0, colour = "grey70") +
  geom_jitter(aes(colour = opacity == 0), width = 0.18, height = 0, alpha = 0.8, size = 2) +
  facet_wrap(~ genename) +
  scale_colour_manual(values = c(`FALSE` = "grey30", `TRUE` = "firebrick"),
                      labels = c(`FALSE` = "measured", `TRUE` = "opacity 0"), name = NULL) +
  scale_y_continuous(labels = scales::label_number(scale = 1e-3, suffix = "k")) +
  labs(x = NULL, y = "Iris opacity",
       title = "Every colony of these genes, all Sen2 plates",
       subtitle = "Zeros appear only once IPTG induces the library gene; the same positions grow without induction") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), legend.position = "top")

ggsave("slides/figures/zero_colonies_by_condition.png", fig, width = 10, height = 4.6, dpi = 150, bg = "white")
cat("\nwrote slides/figures/zero_colonies_by_condition.png\n")
