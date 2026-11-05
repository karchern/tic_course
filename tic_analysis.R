library(tidyverse)
library(ggrepel)

out_dir  <- "output"
p_cutoff <- 0.01

# Load data
colonies <- read_tsv("data/colonies.tsv", show_col_types = FALSE)
#colonies <- read_tsv("data/colonies_predocs.tsv", show_col_types = FALSE)

# 2. Sanity checks: does the experiment behave as it should? 
ggplot(
    colonies %>%
      mutate(is_gfp = ifelse(genename == "gfp", "GFP", "non-GFP"))
      , aes(induction, opacity, fill = is_gfp)) +
  geom_boxplot(outlier.size = 0.4, linewidth = 0.3) +
  scale_fill_manual(values = c("GFP" = 'green', 'non-GFP' = 'grey')) +
  theme_bw()

# 3. Correct the plate edges -------------------------------------------------
colonies <- colonies %>%
  group_by(plate_rep, induction) %>%
  mutate(interior_median_opacity = median(opacity[row %in% 2:15 & column %in% 2:23])) %>%
  mutate(exterior_median_opacity = median(opacity[!(row %in% 2:15 & column %in% 2:23)])) %>%
  mutate(correction_factor = interior_median_opacity / exterior_median_opacity) %>%
  mutate(measurement_location = ifelse(
    row %in% 2:15 & column %in% 2:23, "plate_interior", "plate_edge"
  )) %>%
  mutate(opacity_edge_corrected = if_else(measurement_location == "plate_edge", opacity * correction_factor, opacity))

# 4. Normalise to the GFP controls on the same plate -------------------------
fitness <- colonies %>%
  group_by(plate_rep, induction) %>%
  mutate(normalized_fitness = opacity_edge_corrected / mean(opacity_edge_corrected[genename == "gfp"])) %>%
  ungroup()

# 5. Effect size -------------------------------------------------------------
# How much better does a strain grow with toxin AND blocker than with the toxin
# alone? One ratio per clone per plate copy, so 6 values for each blocker.

ratios <- fitness %>%
  select(
    plate_rep, induction, strain, plate_96, genename, normalized_fitness
  ) %>%
  pivot_wider(names_from = induction, values_from = normalized_fitness) %>%
  mutate(ratio = AraIPTG / IPTG)

# Check for undefined ratios (IPTG values very low -> pinning errors)
ratios %>%
  filter(!is.finite(ratio))

# Go on...

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


# 7. Volcano plot ------------------------------------------------------------
volcano <- ggplot(results, aes(effect, -log10(p_value))) +
  geom_point(alpha = 0.7) +
  theme_bw()

results %>% select(genename, effect, p_value, t.test.pvalue) %>%
  write_tsv(file.path(out_dir, "hits.tsv"))

source("inspect_gene.R")
inspect_gene("50aa_1_55")

