# 1. Load libs and data
library(tidyverse)
library(ggrepel)

out_dir  <- "output"
p_cutoff <- 0.01

# Load data
#colonies <- read_tsv("data/colonies.tsv", show_col_types = FALSE)
colonies <- read_tsv("data/colonies_predocs.tsv", show_col_types = FALSE)



# 2. Sanity checks: does the experiment behave as it should?




# 3. Correct the plate edges -------------------------------------------------




# 4. Normalise to the GFP controls on the same plate -------------------------




# 5. Effect size -------------------------------------------------------------
# How much better does a strain grow with toxin AND blocker than with the toxin
# alone? One ratio per clone per plate copy, so 6 values for each blocker.




# Check for undefined ratios (IPTG values very low -> pinning errors)




# We will run a naive t-test between ratios of library genes against the ratio
# of GFP values. Since we already normalized by gfp (per plate), we can do this
# over all data points




# 7. Volcano plot ------------------------------------------------------------




# Inspect some strong hits


