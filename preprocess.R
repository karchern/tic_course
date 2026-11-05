# Preprocessing for the TIC screen.
#
# Reads the raw Iris measurements, attaches the plate map so every colony knows
# which strain it is, decodes the plate names into experimental conditions, and
# keeps only the plates the analysis uses. Writes one tidy table of colonies.
#
# Run from the repo root:  Rscript preprocess.R
library(tidyverse)

iris_dir <- "data/iris_predocs"
#iris_dir <- "data/iris"

map_file <- "data/AlphaDesign_384wp.identifiers.csv"

out_file <- "data/colonies_predocs.tsv"
#out_file <- "data/colonies.tsv"


# 1. Load --------------------------------------------------------------------
# The map says which strain sits at each of the 384 positions. One 96-well
# source plate is pinned into 384 format, so each well becomes a 2x2 block:
# two clones, each pinned twice. "strain" identifies one clone of one well.
plate_map <- read_csv(map_file, show_col_types = FALSE) %>%
  select(row = row.384, column = column.384, strain, genename, plate_96 = plate)

read_iris <- function(path) {
  read_tsv(path, skip = 6, show_col_types = FALSE) %>%
    mutate(plate = str_remove(basename(path), "\\.jpg\\.iris$"))
}

# FOR THURSDAY( CLAUDE IGNORE):
# - Ensure that 'plate' column has exactly this format: SpecGent-1005-3A_1
colonies <- list.files(iris_dir, "\\.iris$", full.names = TRUE) %>%
  map(read_iris) %>%
  list_rbind() %>%
  left_join(plate_map, by = c("row", "column"))

# The plate name carries the design, e.g. SpecGentAraIPTG-10050201-3A_1:
# antibiotics, inducer concentrations, then which copy of the plate it is.
colonies <- colonies %>%
  separate_wider_delim(plate, "-", names = c("drugs", "inducers", "copy"), cols_remove = FALSE) %>%
  mutate(
    #block     = as.integer(str_sub(copy, 1, 1)),
    plate_rep = str_sub(copy, 2),
#    plasmid   = if_else(block %in% c(1, 2), "p-RcaT", "pir116"),
    plasmid   = 'pir116',
#    retron    = if_else(block %in% c(1, 3), "Sen2", "Eco9"),
    retron    = 'Sen2',
    induced   = str_remove(str_remove(drugs, "Spec"), "Gent|tet"),
    induced   = if_else(induced == "", "none", induced),
    iptg      = case_when(inducers %in% c("100501", "10050201", "1001001", "100100201") ~ "low",
                          inducers %in% c("10051", "1005021", "100101", "10010021")     ~ "high",
                          .default = "none")) %>%
    filter(plate_rep != "A_2")


# 2. Filter ------------------------------------------------------------------
colonies <- colonies %>%
  filter(plasmid == "pir116",   # the high-copy toxin plasmid
         iptg != "high",        # low blocker induction only
         retron == "Sen2",
         !is.na(genename) # positions the map leaves empty
         # all four inductions are kept: the analysis uses IPTG and AraIPTG,
         # the sanity checks need the uninduced and Ara-only plates too
         ) %>%
  group_by(plate_rep, induced, strain) %>%
  ungroup()

# 3. Keep only columns we need
colonies <- colonies %>%
  select(
    row, column, opacity, plate, copy, strain, genename, plate_96, plate_rep, induced
  ) %>%
  rename(induction = induced)

write_tsv(colonies, out_file)
cat("wrote", out_file, "-", nrow(colonies), "colonies,",
    n_distinct(colonies$plate_rep), "plate copies,",
    n_distinct(colonies$genename), "strains\n")
