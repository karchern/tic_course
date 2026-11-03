# What the original script actually does with zero-opacity colonies.
# The effect-size and p-value sections of alphafold_tactic.Rmd are reproduced
# here verbatim (extracted from repro/alphafold_tactic_patched.Rmd), then
# instrumented. Run the patched Rmd first so alpha_df.txt exists.
# Rscript repro/trace_zero_handling.R
suppressMessages(library(tidyverse)); library(data.table)
df <- fread("output/repro_dirty/output/alpha_df.txt", sep = "	", header = TRUE)
df = df %>%
  mutate(meanop = ave(opacity.corrected, plate_id, strain,FUN= mean))

#norm by gfp per plate
normgfp <- function (df) {
  
  meangfp <- dplyr::filter (df,  genename == "gfp") %>% dplyr::group_by (plate_id) %>%
    dplyr::reframe(meangfp = mean(meanop))
  
  test <- inner_join (df, meangfp, by = "plate_id") %>%
    dplyr::mutate(normop = opacity.corrected / meangfp)
  
  return(test)
}

df <- normgfp(df)

df = df %>%
  mutate(copynum_ind = sapply(strsplit(as.character(plate_id), "-"),'[',1),
         concentrations = sapply(strsplit(as.character(plate_id), "-"),'[',2),
         condition = sapply(strsplit(as.character(plate_id), "-"),'[',3)) %>%
  mutate(replicate_rcat = sapply(strsplit(as.character(plate_id), "_"),'[',2),
        condition = sapply(strsplit(as.character(condition), "_"),'[',1))


df = df %>%
  mutate(copynum = case_when(condition %in% c(1,2) ~ "p-RcaT",
                               TRUE ~ "pir116"),
         condition_replicate_rcat = interaction(condition, replicate_rcat, sep = "_"),
         rcat_biol_rep = case_when(condition_replicate_rcat %in% c("1_1", "1_2", "2_1", "2_2",  "3_1", "3_2", "4_1", "4_2") ~ 1,
                                 TRUE ~ 2),
         rcat_tech_rep = case_when(condition_replicate_rcat %in% c("1_1", "1_3", "2_1", "2_3",  "3_1", "3_3", "4_1", "4_3") ~ 1,
                                   TRUE ~ 2),
         retron = case_when(condition %in% c(1,3) ~ "Sen2",
                            condition %in% c(2,4) ~ "Eco9"))

df =df %>%
  mutate(antibiotic = ifelse(grepl("Gent", copynum_ind), "Gent", 
                             ifelse(grepl("tet", copynum_ind), "tet", NA)))

remove_patterns <- function(text, pattern1, pattern2) {
  result <- ifelse(text == pattern1 | text == pattern2, "no induction", gsub(paste0(pattern1, "|", pattern2), "", text))
  return(result)
}

df =df %>%
  mutate(induction = remove_patterns(copynum_ind, pattern1="Gent", pattern2="tet"))



df = df %>%
  mutate(IPTGconc = case_when(concentrations %in% c("10050201", "100501", "100100201","1001001") ~ "0.1mM",
                              concentrations %in% c("1005021", "10051", "10010021","100101") ~ "1mM",
                              TRUE ~ "0mM"),
         induction = interaction(induction, IPTGconc, sep = "_"))


# effect size calculation
test = df %>%
  filter(normop > 0) %>%
  dplyr::select(normop, genename, strain, replicate, rcat_biol_rep, rcat_tech_rep, induction, plate_id, copynum, retron) %>%
  pivot_wider(id_cols = c(strain,genename, replicate, rcat_biol_rep, rcat_tech_rep, copynum, retron), names_from = induction, values_from = normop) %>%
  mutate(efsize_low = case_when(IPTG_0.1mM != 0 ~ AraIPTG_0.1mM - (IPTG_0.1mM * Ara_0mM),
                                TRUE ~ AraIPTG_0.1mM),
         efsize_high = case_when(IPTG_1mM != 0 ~ AraIPTG_1mM - (IPTG_1mM * Ara_0mM),
                                TRUE ~ AraIPTG_1mM))
  # mutate(efsize_low = case_when(IPTG_0.1mM != 0 ~ AraIPTG_0.1mM/IPTG_0.1mM,
  #                               TRUE ~ AraIPTG_0.1mM),
  #        efsize_high = case_when(IPTG_1mM != 0 ~ AraIPTG_1mM/IPTG_1mM,
  #                               TRUE ~ AraIPTG_1mM))



is.nan.data.frame <- function(x)
do.call(cbind, lapply(x, is.nan))
test[is.nan(test)] <- 0


cat("\nA. rows per condition in df (before any filter)\n")
df %>% group_by(induction) %>% summarise(rows = n(), zeros = sum(normop == 0)) %>% as.data.frame() %>% print()

cat("\nB. why is a condition NA in the wide table? missing plate vs zero\n")
wide_all <- df %>%
  dplyr::select(normop, genename, strain, replicate, rcat_biol_rep, rcat_tech_rep, induction, copynum, retron) %>%
  pivot_wider(id_cols = c(strain, genename, replicate, rcat_biol_rep, rcat_tech_rep, copynum, retron),
              names_from = induction, values_from = normop)
for (cond in c("Ara_0mM", "IPTG_0.1mM", "AraIPTG_0.1mM")) {
  cat(sprintf("   %-14s NA after filter: %4d   (never measured: %4d, measured as zero: %4d)\n", cond,
              sum(is.na(test[[cond]])), sum(is.na(wide_all[[cond]])), sum(wide_all[[cond]] == 0, na.rm = TRUE)))
}

cat("\nC. cause of a missing efsize_low, per colony position\n")
test %>%
  mutate(cause = case_when(is.na(Ara_0mM)       ~ "no Ara value (missing plate or zero)",
                           is.na(AraIPTG_0.1mM) ~ "AraIPTG zero",
                           is.na(IPTG_0.1mM)    ~ "IPTG zero only",
                           TRUE                 ~ "complete")) %>%
  count(cause, efsize_low_is_NA = is.na(efsize_low)) %>% as.data.frame() %>% print()

cat("\nD. IPTG-zero-only rows: efsize equals AraIPTG alone?\n")
test %>% filter(is.na(IPTG_0.1mM), !is.na(Ara_0mM), !is.na(AraIPTG_0.1mM)) %>%
  transmute(genename, AraIPTG = round(AraIPTG_0.1mM, 3), efsize_low = round(efsize_low, 3)) %>%
  head(4) %>% as.data.frame() %>% print()

cat("\nE. per-gene consequence (pir116, low induction)\n")
test %>% filter(copynum == "pir116") %>% group_by(genename, retron) %>%
  summarise(reps = n(), usable = sum(!is.na(efsize_low)), .groups = "drop") %>%
  summarise(groups = n(), lose_some = sum(usable < reps), lose_all = sum(usable == 0),
            min_usable = min(usable), median_usable = median(usable)) %>% as.data.frame() %>% print()

cat("\nF. the p-value path (signif) is NOT filtered\n")
lowrows <- df %>% filter(induction %in% c("IPTG_0.1mM", "AraIPTG_0.1mM", "Ara_0mM", "no induction_0mM"))
cat("   rows entering lm/anova:", nrow(lowrows), " of which normop == 0:", sum(lowrows$normop == 0), "\n")

cat("\nG. the 768 rows with no Ara measurement at all\n")
wide_all %>% filter(is.na(Ara_0mM)) %>% count(copynum, retron) %>% as.data.frame() %>% print()
cat("   Ara plates present:\n")
df %>% filter(induction == "Ara_0mM") %>% distinct(plate_id) %>% arrange(plate_id) %>% pull(plate_id) %>% print()
