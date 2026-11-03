# Rebuild the published (Sept 2023) analysis from the paper methods and compare
# it to reference/published_hits_pir116_low.txt.
# Input: output/repro_dirty/output/alpha_df.txt (written by alphafold_tactic_patched.Rmd)
# Run from the repo root: Rscript repro/compare_published.R
suppressMessages(library(tidyverse))

df  <- read_tsv("output/repro_dirty/output/alpha_df.txt", show_col_types = FALSE)
pub <- read_tsv("reference/published_hits_pir116_low.txt", show_col_types = FALSE) %>%
  # the table's own "hit" column is an older, looser call (25 Sen2 / 19 Eco9 hits);
  # the paper (Fig. 7B/C: 17 Sen2 / 5 Eco9) applies p < 0.01 & effect size > 1.1 (Sen2) / 1.3 (Eco9)
  mutate(hit = ifelse(pval < 0.01 & efsize > ifelse(retron == "Sen2", 1.1, 1.3), "hit", "neutral"))

# pir116 (Gent) plates, low blocker induction only
d <- df %>%
  filter(grepl("^Gent", plate_id)) %>%
  separate(plate_id, c("cond", "conc", "rep"), sep = "-", remove = FALSE) %>%
  filter(conc %in% c("1005", "100502", "100501", "10050201")) %>%
  mutate(retron   = ifelse(substr(rep, 1, 1) == "3", "Sen2", "Eco9"),
         platerep = substr(rep, 3, 3),
         cond     = recode(cond, Gent = "none", GentAra = "Ara", GentIPTG = "IPTG", GentAraIPTG = "AraIPTG"),
         well     = paste0(row.96, column.96),
         # clone (biological replicate) as in data/AlphaDesign_384wp.identifiers.csv: odd 384 rows = clone 1
         clone    = ifelse(row %% 2 == 1, 1, 2))

# normalise to the per-plate GFP mean; average the 2 technical replicates (zeros kept)
fit <- d %>%
  group_by(plate_id) %>%
  mutate(normop = opacity.corrected / mean(opacity.corrected[genename == "gfp"])) %>%
  group_by(retron, cond, platerep, well, clone, genename) %>%
  summarise(fit = mean(normop), .groups = "drop")

# effect = (toxin + blocker) / toxin, per clone and plate replicate (n = 6 per blocker)
ratios <- fit %>%
  filter(cond %in% c("Ara", "AraIPTG"), genename != "empty") %>%
  pivot_wider(names_from = cond, values_from = fit) %>%
  mutate(ratio = AraIPTG / Ara) %>%
  filter(is.finite(ratio))

# GFP reference: one value per GFP strain (averaged over plate replicates, n = 16)
gfp <- ratios %>%
  filter(genename == "gfp") %>%
  group_by(retron, well, clone) %>%
  summarise(ratio = mean(ratio), .groups = "drop")

res <- ratios %>%
  filter(genename != "gfp") %>%
  group_by(retron, genename) %>%
  summarise(efsize = mean(ratio) / mean(gfp$ratio[gfp$retron == first(retron)]),
            pval   = t.test(ratio, gfp$ratio[gfp$retron == first(retron)])$p.value,
            .groups = "drop") %>%
  mutate(hit = ifelse(pval < 0.01 & efsize > ifelse(retron == "Sen2", 1.1, 1.3), "hit", "neutral"))

cmp <- inner_join(res, pub, by = c("genename", "retron"), suffix = c("_rebuilt", "_published"))

stats <- cmp %>%
  group_by(retron) %>%
  summarise(n                = n(),
            efsize_pearson   = cor(efsize_rebuilt, efsize_published),
            efsize_spearman  = cor(efsize_rebuilt, efsize_published, method = "spearman"),
            log10p_pearson   = cor(log10(pval_rebuilt), log10(pval_published)),
            log10p_spearman  = cor(pval_rebuilt, pval_published, method = "spearman"),
            hits_rebuilt     = sum(hit_rebuilt == "hit"),
            hits_published   = sum(hit_published == "hit"),
            hits_shared      = sum(hit_rebuilt == "hit" & hit_published == "hit"))
print(stats, width = Inf)

dir.create("output/compare", showWarnings = FALSE, recursive = TRUE)
write_tsv(cmp, "output/compare/rebuilt_vs_published.tsv")

long <- cmp %>%
  transmute(retron, genename, hit_published,
            `effect size`  = paste(efsize_rebuilt, efsize_published),
            `-log10(p)`    = paste(-log10(pval_rebuilt), -log10(pval_published))) %>%
  pivot_longer(c(`effect size`, `-log10(p)`), names_to = "metric") %>%
  separate(value, c("rebuilt", "published"), sep = " ", convert = TRUE)
labs_df <- stats %>%
  transmute(retron,
            `effect size` = sprintf("r = %.2f, rho = %.2f", efsize_pearson, efsize_spearman),
            `-log10(p)`   = sprintf("r = %.2f, rho = %.2f", log10p_pearson, log10p_spearman)) %>%
  pivot_longer(-retron, names_to = "metric", values_to = "label")

p <- ggplot(long, aes(published, rebuilt)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey60") +
  geom_point(aes(colour = hit_published), alpha = 0.7) +
  geom_text(data = labs_df, aes(label = label), x = -Inf, y = Inf, hjust = -0.05, vjust = 1.3, size = 3.2) +
  facet_wrap(metric ~ retron, scales = "free") +
  scale_colour_manual(values = c(hit = "firebrick", neutral = "grey40")) +
  labs(x = "published (230908_hits_onlyLowInduction_pir116.txt)", y = "rebuilt from paper methods",
       colour = "paper hit call", title = "Rebuilt vs published: pir116, low blocker induction") +
  theme_bw()
ggsave("output/compare/rebuilt_vs_published.png", p, width = 9, height = 8, dpi = 120)
