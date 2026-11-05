# 1. Load libs and data
library(tidyverse)
library(ggrepel)

out_dir  <- "output"
p_cutoff <- 0.01

# Load data
#colonies <- read_tsv("data/colonies.tsv", show_col_types = FALSE)
colonies <- read_tsv("data/colonies_predocs.tsv", show_col_types = FALSE)

# 2. Sanity checks: plot raw opacities by condition
p <- ggplot(
    colonies, aes(
        x = induction, y = opacity
    )
) + 
geom_boxplot() +
theme_bw()

# 3. Correct the plate edges -------------------------------------------------
# Do it next time :)

# 4. Normalise to the GFP controls on the same plate -------------------------
fitness <- colonies %>%
    group_by(
        plate_rep, induction
    ) %>%
    mutate(
        normalized_fitness = opacity /  mean(opacity[genename == 'gfp'])
    ) %>%
    ungroup()

# 5. Effect size -------------------------------------------------------------
# How much better does a strain grow with toxin AND blocker than with the toxin
# alone? One ratio per clone per plate copy, so 6 values for each blocker.

fitness <- fitness %>%
    select(
        plate_rep, induction, strain, plate_96, genename, normalized_fitness
    ) %>%
    pivot_wider(
        names_from = induction, values_from = normalized_fitness
    ) %>%
    mutate(
        ratio = AraIPTG / IPTG
    )
    # mutate(
    #     ratio = blocker_plus_toxicity/toxicity
    # )

# Check for undefined ratios (IPTG values very low -> pinning errors)
fitness <- fitness %>% 
    filter(is.finite(ratio))

# We will run a naive t-test between ratios of library genes against the ratio
# of GFP values. Since we already normalized by gfp (per plate), we can do this
# over all data points

t_test_results <- fitness %>%
    group_by(genename) %>% 
    nest() %>%
    mutate(
        p_val = map_dbl(data, \(x) {
            t.test(x$ratio, mu = 1)$p.value
        }),
        median_ratio = map_dbl(data, \(x) {
            median(x$ratio)
        })
    )

# 7. Volcano plot ------------------------------------------------------------
volcano_plot <- ggplot(
    t_test_results, aes(x = median_ratio, y = -log10(p_val)) 
) +
    geom_point() +
    theme_bw() +
    geom_text_repel(
        data = t_test_results %>% filter(median_ratio > 1) %>% arrange(p_val) %>% head(10),
        aes(x = median_ratio, y = -log10(p_val), label = genename)
    ) +
    NULL



# Inspect some strong hits
#inspect_gene("50aa_1_55", "iris_predocs")

results_students <- t_test_results %>%
    mutate(
        p_val_students = p_val,
        median_ratio_students = median_ratio
    ) %>%
    select(
        genename,
        p_val_students,
        median_ratio_students
    ) 

top_10_students <- results_students %>%
    arrange(p_val_students) %>%
    filter(median_ratio_students > 1) %>%
    head(10) %>%
    select(genename) %>%
    mutate(top10 = TRUE)
results_students <- left_join(results_students, top_10_students)

results_paper <- read_tsv('/Users/karcher/Downloads/HitList_Paper.txt') %>%
    mutate(
        p_val_paper = pval,
        median_ratio_paper = efsize
    ) %>%
    select(
        genename,
        p_val_paper,
        median_ratio_paper
    ) %>%
    mutate(top10 = TRUE)

top_10_paper <- results_paper %>%
    arrange(p_val_paper) %>%
    filter(median_ratio_paper > 1) %>%
    head(10) %>%
    select(genename)
results_paper <- left_join(results_paper, top_10_paper)

a <- inner_join(results_paper, results_students, by = 'genename') %>%
    mutate(top10.x = ifelse(is.na(top10.x), FALSE, TRUE)) %>%
    mutate(top10.y = ifelse(is.na(top10.y), FALSE, TRUE)) %>%
    mutate(
        hit_group = case_when(
            top10.x & top10.y ~ 'strong_hit_both',
            top10.x | top10.y ~ 'strong_hit_one\nbut not other',
            TRUE ~ 'else'
        )
    )

ggplot(data = a, aes(x = median_ratio_paper, y = median_ratio_students))