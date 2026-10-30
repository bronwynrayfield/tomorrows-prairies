
################################################################################
##         ANALYSIS: ACCURACY ASSESSMENT (BINARY GRASS VS NON-GRASS)          ##
################################################################################

##----------------------- Alison Long - June 2026 -------------------------

# STEP 8: Assess accuracy of datasets mapping grasslands by comparing them
#         to ground-truthed point data for 2015 to present. This script combines 
#         the compiled results from script 8.1 and the crop insurance-derived
#         results from script 8.2.

# Ran separate scripts for the two types of point data - will merge here and
# calculate metrics

## Goal:
##   Assess the accuracy of multiple grassland rasters vs. ground-truth points
##   using time-aligned rules that depend on class (NATIVE or TAME).
##
##   1) Handles PLI as 4 separate rasters: pli_ap, pli_cu, pli_mg, pli_mm
##   2) Allows MULTIPLE valid metadata years per raster (e.g., pli_mm = 2020|2021)
##   3) Filters points by class-based time rules using ANY allowed year per raster
##   4) Computes a matched_year per record for per-year summaries
##   5) Produces both "grouped" (legacy: pli combined) and "atomic" (each raster) outputs
##
## Time Rules (per record):
##   - must have non-NA raster prediction (0 or 1)
##   - class-based time rule using allowed years for that raster
##   - rules differ for the compiled vs crop insurance datasets
##
## Exclusions:
##   - Source used for training: exclude if Training_norm matches either atomic raster code or grouped alias
##
# ----------------------------
# 1) Setup
# ----------------------------

# load packages
library(readxl)
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(purrr)
library(tibble)
library(sf)
library(scales)

# Small utility: "x %||% y" returns y if x is NULL
`%||%` <- function(x, y) if (is.null(x)) y else x

# ----------------------------
# 2) Directories and file paths 
# ----------------------------
acc_cm_dir <- "R:/CanadianGrasslands/WG_Classification/results/accuracy_assessment/binary/compiled_derived"
acc_ci_dir <- "R:/CanadianGrasslands/WG_Classification/results/accuracy_assessment/binary/crop_ins_derived/v2_of_dataset"
out_dir <- "R:/CanadianGrasslands/WG_Classification/results/accuracy_assessment/binary/combined_CI_and_compiled"

acc_cm_file <- file.path(acc_cm_dir, "27Jul2026/accuracy_assessment_data_binary_pli_sep_compiled_27Jul2026.csv")
acc_ci_file <- file.path(acc_ci_dir, "23Jul2026/lenient/accuracy_assessment_data_binary_pli_sep_cropinsV2_lenient_23Jul2026.csv")
#acc_ci_file <- file.path(acc_ci_dir, "23Jul2026/standard/accuracy_assessment_data_binary_pli_sep_cropinsV2_23Jul2026.csv")
metadata_file <- "R:/CanadianGrasslands/WG_Classification/github/tomorrows-prairies/processeddatasets/Grass/metadata_full_extent_binary_pli_sep_noplow_aci2024and2025.csv"

# ----------------------------
# 3) Configuration
# ----------------------------
# IMPORTANT: Handle PLI as four separate rasters in the points file and metadata.
pli_variants <- c("pli_ap", "pli_cu", "pli_mg", "pli_mm")

# Grouping helper: collapse pli_* -> "pli"; all others remain themselves
raster_to_group <- function(x) ifelse(x %in% pli_variants, "pli", x)

# ----------------------------
# 4) Load points and clean
# ----------------------------
# Load
acc_cm <- readr::read_csv(acc_cm_file, show_col_types = FALSE)
acc_ci <- readr::read_csv(acc_ci_file, show_col_types = FALSE)

# Inspect structure
str(acc_cm)
str(acc_ci)

# Count ACI_field_2024 point/raster combos in acc_cm (12,602)
acc_cm %>%
  filter(Sorc_ID == "ACI_field_2024") %>% 
  tally()

# Count ACI_field_2024 point/raster combos in acc_cm (22,618) - larger due to lenient year matching for tame/other
acc_ci %>%
  filter(Sorc_ID == "ACI_field_2024") %>% 
  tally()

# Exclude ACI_field_2024 points from acc_cm since they are present in acc_ci also
acc_cm <- acc_cm %>%
  filter(Sorc_ID != "ACI_field_2024")

# Add missing columns
acc_ci <- acc_ci %>%
  mutate(Ecrg_ID = NA_real_)

# Clean up ecoregion names
acc_ci <- acc_ci %>%
  mutate(Ecrgn_n = case_when(
    Ecrgn_n == "Okanogan Highland"~ "Okanagan Highland", 
    Ecrgn_n == "Okanogan Range"~ "Okanagan Range",
    Ecrgn_n == "Southwest Manotoba Uplands"~ "Southwest Manitoba Uplands",
    Ecrgn_n == "St-Laurent Lowlands"~ "St. Lawrence Lowlands",
    TRUE ~ Ecrgn_n))

# exclude two rows where pli variant is in wrong ecodistrict (see code row 1030 below)
# View the rows that will be excluded
excluded_rows <- acc_ci %>%
  filter(Ecds_ID %in% c(745, 830) & raster == "pli_mm")

print(excluded_rows, width = Inf)

#    Pt_ID Ecds_ID Ecrgn_n         Ecprv Training Sorc_ID     LC_grass Yer_Cll LC_Cl_N raster  pred
# 1 298163     745 Aspen Parkland   10.2 NA       ACI_CI_2020        0    2020 OTHER   pli_mm     0
# 2 303882     830 Mixed Grassland  10.3 NA       ACI_CI_2020        0    2020 OTHER   pli_mm     0

#   raster_group allowed_years matched_year quad  Ecrg_ID
# 1 pli          2020;2021             2020 tn         NA
# 2 pli          2020;2021             2020 tn         NA

# Remove them from the dataset
acc_ci <- acc_ci %>%
  filter(!(Ecds_ID %in% c(745, 830) & raster == "pli_mm"))

# Make sure aci raster has same name for both
# acc_cm <- acc_cm %>%
#   mutate(raster = if_else(raster == "aci", "aci2024", raster))

# acc_cm <- acc_cm %>%
#   mutate(raster_group = if_else(raster_group == "aci", "aci2024", raster_group))

# Reorder to match
acc_ci <- acc_ci %>%
  select(names(acc_cm))

# ----------------------------
# 5) Create unique point IDs
# ----------------------------
# Add column for dataset 
acc_ci <- acc_ci %>%
  mutate(Dataset = "CI")

acc_cm <- acc_cm %>%
  mutate(Dataset = "CM")

# Create column for concatenation
acc_ci <- acc_ci %>%
  mutate(Pt_Dataset = stringr::str_c(Pt_ID, Dataset, sep = "_"))

acc_cm <- acc_cm %>%
  mutate(Pt_Dataset = stringr::str_c(Pt_ID, Dataset, sep = "_"))

# Remove original point ID columns
acc_ci <- acc_ci %>%
  select(-Pt_ID)

acc_cm <- acc_cm %>%
  select(-Pt_ID)

# ----------------------------
# 6) Merge into one dataset
# ----------------------------
combined <- bind_rows(acc_ci, acc_cm)
nrow(combined) == nrow(acc_ci) + nrow(acc_cm)

# ----------------------------
# 7) Export combined dataset
# ----------------------------
# Collapse lists to strings 
combined_out <- combined %>%
  mutate(allowed_years = map_chr(allowed_years, ~ paste(.x, collapse = ";")))

# Export
readr::write_csv(combined_out, file.path(out_dir, "accuracy_assessment_data_binary_pli_sep_CIandCM_28Jul2026.csv"))

# ----------------------------
# 8) Metrics helpers
# ----------------------------
safe_div <- function(n, d) ifelse(d == 0 | is.na(d), NA_real_, n / d)

compute_counts <- function(df) {
  # df must have columns: LC_grass (0/1), pred (0/1)
  tp <- sum(df$LC_grass == 1 & df$pred == 1)
  tn <- sum(df$LC_grass == 0 & df$pred == 0)
  fp <- sum(df$LC_grass == 0 & df$pred == 1)
  fn <- sum(df$LC_grass == 1 & df$pred == 0)
  tibble(TP = tp, TN = tn, FP = fp, FN = fn, N_evaluated = tp + tn + fp + fn)
}

compute_metrics <- function(tp, tn, fp, fn) {
  n <- tp + tn + fp + fn
  
  acc  <- safe_div(tp + tn, n)
  prec <- safe_div(tp, tp + fp)
  rec  <- safe_div(tp, tp + fn)
  spec <- safe_div(tn, tn + fp)
  npv  <- safe_div(tn, tn + fn)
  
  f1 <- ifelse(is.na(prec) | is.na(rec) | (prec + rec) == 0,
               NA_real_, 2 * prec * rec / (prec + rec))
  
  bal <- if (is.na(rec) || is.na(spec)) NA_real_ else mean(c(rec, spec)) # require recall and specificity
  
  # Cohen's Kappa
  pe <- safe_div(
    (as.numeric(tp + fp) * as.numeric(tp + fn) +
       as.numeric(fn + tn) * as.numeric(fp + tn)),
    as.numeric(n)^2)
  kappa <- safe_div(acc - pe, 1 - pe)
  
  # MCC
  mcc <- safe_div(
    (as.numeric(tp) * as.numeric(tn) - as.numeric(fp) * as.numeric(fn)),
    sqrt((as.numeric(tp + fp) * as.numeric(tp + fn) *
            as.numeric(tn + fp) * as.numeric(tn + fn)))
  )
  
  # Error rates
  fpr <- safe_div(fp, fp + tn)
  fnr <- safe_div(fn, fn + tp)
  
  # Prevalence
  prev <- safe_div(tp + fn, n)
  
  tibble(
    Overall_Accuracy     = acc,
    Precision_grass      = prec,
    Recall_grass         = rec,
    Specificity_notgrass = spec,
    NPV_notgrass         = npv,
    F1_grass             = f1,
    Balanced_Accuracy    = bal,
    Kappa                = kappa,
    MCC                  = mcc,
    False_Positive_Rate  = fpr,
    False_Negative_Rate  = fnr,
    Prevalence_grass     = prev
  )
}

summarize_block <- function(df, key_cols) {
  # Generic summarizer: group by key_cols, then compute counts + metrics
  df %>%
    group_by(across(all_of(key_cols))) %>%
    group_modify(~ {
      counts  <- compute_counts(.x)
      metrics <- compute_metrics(counts$TP, counts$TN, counts$FP, counts$FN)
      bind_cols(counts, metrics)
    }) %>%
    ungroup()
}

# ----------------------------
# 9) Per-raster summaries (ATOMIC and GROUPED)
# ----------------------------

# ATOMIC: one row per (raster, matched_year)
summary_by_raster_atomic <- summarize_block(
  combined,
  key_cols = c("raster", "matched_year")
) %>%
  rename(raster_year = matched_year) %>%
  arrange(desc(Overall_Accuracy))

# GROUPED: collapse pli_* -> "pli"; combine counts across all matched_year within group
summary_by_raster_grouped <- combined %>%
  mutate(raster = raster_group) %>%
  group_by(raster) %>%
  summarise(
    TP = sum(LC_grass == 1 & pred == 1),
    TN = sum(LC_grass == 0 & pred == 0),
    FP = sum(LC_grass == 0 & pred == 1),
    FN = sum(LC_grass == 1 & pred == 0),
    N_evaluated = TP + TN + FP + FN,
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(metrics = list(compute_metrics(TP, TN, FP, FN))) %>%
  ungroup() %>%
  tidyr::unnest_wider(metrics) %>%
  mutate(raster_year = NA_integer_) %>%
  relocate(raster, raster_year) %>%
  arrange(desc(Overall_Accuracy))

# ----------------------------
# 10) Confusion counts tables
# ----------------------------
confusion_counts_atomic <- summary_by_raster_atomic %>%
  select(raster, raster_year, TP, TN, FP, FN, N_evaluated)

confusion_counts_grouped <- summary_by_raster_grouped %>%
  select(raster, raster_year, TP, TN, FP, FN, N_evaluated)

# ----------------------------
# 11) Per ground-truth year breakdowns
# ----------------------------
# ATOMIC: per raster × matched_year × Yer_Cll
per_year_summary_atomic <- summarize_block(
  combined,
  key_cols = c("raster", "matched_year", "Yer_Cll")
) %>%
  rename(raster_year = matched_year) %>%
  arrange(raster, Yer_Cll)

# GROUPED: per grouped raster × Yer_Cll (mix of multiple matched years inside)
per_year_summary_grouped <- summarize_block(
  combined %>% mutate(raster = raster_group),
  key_cols = c("raster", "Yer_Cll")
) %>%
  mutate(raster_year = NA_integer_) %>%
  relocate(raster, raster_year) %>%
  arrange(raster, Yer_Cll)

# ----------------------------
# 12) Save outputs (CSV)
# ----------------------------
# GROUPED keep legacy filenames
readr::write_csv(summary_by_raster_grouped, file.path(out_dir, "raster_accuracy_summary_combined.csv"))
readr::write_csv(confusion_counts_grouped,  file.path(out_dir, "raster_confusion_matrices_combined.csv"))
readr::write_csv(per_year_summary_grouped,  file.path(out_dir, "raster_per_year_accuracy_combined.csv"))

# ATOMIC with _atomic suffix
readr::write_csv(summary_by_raster_atomic, file.path(out_dir, "raster_accuracy_summary_atomic_combined.csv"))
readr::write_csv(confusion_counts_atomic,  file.path(out_dir, "raster_confusion_matrices_atomic_combined.csv"))
readr::write_csv(per_year_summary_atomic,  file.path(out_dir, "raster_per_year_accuracy_atomic_combined.csv"))

# ----------------------------
# 13) Accuracy by Ecodistrict
# ----------------------------

# Ensure ecodistrict column exists
if (!"Ecds_ID" %in% names(combined)) {
  stop("Ecds_ID column not found in combined.")
}

# ---- ATOMIC: raster × ecodistrict ----
accuracy_by_ecodistrict_atomic <- summarize_block(
  combined,
  key_cols = c("raster", "Ecds_ID")
) %>%
  arrange(raster, Ecds_ID)

# ---- GROUPED: collapse pli_* -> pli ----
accuracy_by_ecodistrict_grouped <- combined %>%
  mutate(raster = raster_group) %>%
  group_by(raster, Ecds_ID) %>%
  summarise(
    TP = sum(LC_grass == 1 & pred == 1),
    TN = sum(LC_grass == 0 & pred == 0),
    FP = sum(LC_grass == 0 & pred == 1),
    FN = sum(LC_grass == 1 & pred == 0),
    N_evaluated = TP + TN + FP + FN,
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(metrics = list(compute_metrics(TP, TN, FP, FN))) %>%
  ungroup() %>%
  tidyr::unnest_wider(metrics) %>%
  arrange(raster, Ecds_ID)

# ---- OPTIONAL: include matched raster year (atomic finer breakdown) ----
accuracy_by_ecodistrict_year_atomic <- summarize_block(
  combined,
  key_cols = c("raster", "matched_year", "Ecds_ID")
) %>%
  rename(raster_year = matched_year) %>%
  arrange(raster, raster_year, Ecds_ID)

#------ Save CSVs -----

# GROUPED
readr::write_csv(
  accuracy_by_ecodistrict_grouped,
  file.path(out_dir, "accuracy_by_ecodistrict_grouped_combinedCIandCM.csv")
)

# ATOMIC
readr::write_csv(
  accuracy_by_ecodistrict_atomic,
  file.path(out_dir, "accuracy_by_ecodistrict_atomic_combinedCIandCM.csv")
)

# ATOMIC + YEAR (detailed output)
readr::write_csv(
  accuracy_by_ecodistrict_year_atomic,
  file.path(out_dir, "accuracy_by_ecodistrict_year_atomic_combinedCIandCM.csv")
)

# ----------------------------
# 14) Create hierarchical fallbacks for ecodistricts with small N
# ----------------------------

#---- Define threshold ----
min_n <- 20

#---- Compute ecoregion-level metrics ----
accuracy_by_ecoregion_atomic <- combined %>%
  group_by(raster, Ecrgn_n) %>%
  summarise(
    TP = sum(LC_grass == 1 & pred == 1),
    TN = sum(LC_grass == 0 & pred == 0),
    FP = sum(LC_grass == 0 & pred == 1),
    FN = sum(LC_grass == 1 & pred == 0),
    N_evaluated = TP + TN + FP + FN,
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(metrics = list(compute_metrics(TP, TN, FP, FN))) %>%
  ungroup() %>%
  tidyr::unnest_wider(metrics)

accuracy_by_ecoregion_grouped <- combined %>%
  mutate(raster = raster_group) %>%
  group_by(raster, Ecrgn_n) %>%
  summarise(
    TP = sum(LC_grass == 1 & pred == 1),
    TN = sum(LC_grass == 0 & pred == 0),
    FP = sum(LC_grass == 0 & pred == 1),
    FN = sum(LC_grass == 1 & pred == 0),
    N_evaluated = TP + TN + FP + FN,
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(metrics = list(compute_metrics(TP, TN, FP, FN))) %>%
  ungroup() %>%
  tidyr::unnest_wider(metrics)

#---- Compute ecoprovince-level metrics ----
accuracy_by_ecoprovince_atomic <- combined %>%
  group_by(raster, Ecprv) %>%
  summarise(
    TP = sum(LC_grass == 1 & pred == 1),
    TN = sum(LC_grass == 0 & pred == 0),
    FP = sum(LC_grass == 0 & pred == 1),
    FN = sum(LC_grass == 1 & pred == 0),
    N_evaluated = TP + TN + FP + FN,
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(metrics = list(compute_metrics(TP, TN, FP, FN))) %>%
  ungroup() %>%
  tidyr::unnest_wider(metrics)

accuracy_by_ecoprovince_grouped <- combined %>%
  mutate(raster = raster_group) %>%
  group_by(raster, Ecprv) %>%
  summarise(
    TP = sum(LC_grass == 1 & pred == 1),
    TN = sum(LC_grass == 0 & pred == 0),
    FP = sum(LC_grass == 0 & pred == 1),
    FN = sum(LC_grass == 1 & pred == 0),
    N_evaluated = TP + TN + FP + FN,
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(metrics = list(compute_metrics(TP, TN, FP, FN))) %>%
  ungroup() %>%
  tidyr::unnest_wider(metrics)

#accuracy_by_ecoprovince <- accuracy_by_ecoprovince %>%
#  rename(
#    TP_ecoP = TP,
#    TN_ecoP = TN,
#    FP_ecoP = FP,
#    FN_ecoP = FN,
#    N_evaluated_ecoP = N_evaluated
#  )

accuracy_by_ecoprovince_atomic <- accuracy_by_ecoprovince_atomic %>%
  rename_with(~ paste0(.x, "_ecoP"), -c(raster, Ecprv))

accuracy_by_ecoprovince_grouped <- accuracy_by_ecoprovince_grouped %>%
  rename_with(~ paste0(.x, "_ecoP"), -c(raster, Ecprv))

#---- Define mapping ----

# Get mapping (one row per ecodistrict)
ecodist_to_regions <- combined %>%
  distinct(Ecds_ID, Ecrgn_n, Ecprv)

# examime ecodistricts associated with more than one ecoregon/prov
combined %>%
  distinct(Ecds_ID, Ecrgn_n, Ecprv) %>%
  count(Ecds_ID) %>%
  filter(n > 1)

combined %>%
  distinct(Ecds_ID, Ecrgn_n, Ecprv) %>%
  filter(Ecds_ID %in%
           (combined %>%
              distinct(Ecds_ID, Ecrgn_n, Ecprv) %>%
              count(Ecds_ID) %>%
              filter(n > 1) %>%
              pull(Ecds_ID))) %>%
  print(n = Inf)

# Ecds_ID Ecrgn_n                    Ecprv
# 1    1010 Okanogan Highland           14.3
# 2    1009 Okanogan Range              14.3
# 3     543 St-Laurent Lowlands          8.1
# 4     544 St-Laurent Lowlands          8.1
# 5     546 St-Laurent Lowlands          8.1
# 6     541 St-Laurent Lowlands          8.1
# 7     547 St-Laurent Lowlands          8.1
# 8     545 St-Laurent Lowlands          8.1
# 9     542 St-Laurent Lowlands          8.1
# 10    1011 Okanogan Highland           14.3
# 11     855 Southwest Manotoba Uplands  10.2
# 12    1010 Okanagan Highland           14.3
# 13    1009 Okanagan Range              14.3
# 14    1011 Okanagan Highland           14.3
# 15     543 St. Lawrence Lowlands        8.1
# 16     544 St. Lawrence Lowlands        8.1
# 17     546 St. Lawrence Lowlands        8.1
# 18     541 St. Lawrence Lowlands        8.1
# 19     547 St. Lawrence Lowlands        8.1
# 20     545 St. Lawrence Lowlands        8.1
# 21     542 St. Lawrence Lowlands        8.1
# 22     855 Southwest Manitoba Uplands  10.2

# Join
accuracy_by_ecodistrict_atomic2 <- accuracy_by_ecodistrict_atomic %>%
  left_join(ecodist_to_regions, by = "Ecds_ID")

accuracy_by_ecodistrict_grouped2 <- accuracy_by_ecodistrict_grouped %>%
  left_join(ecodist_to_regions, by = "Ecds_ID")

# Check whether an ecodistrict belongs to multiple regions/provinces
ecodist_to_regions %>%
  count(Ecds_ID) %>%
  filter(n > 1)

#------- Join ----------
fallback_joined_atomic <- accuracy_by_ecodistrict_atomic2 %>%
  left_join(
    accuracy_by_ecoregion_atomic,
    by = c("raster", "Ecrgn_n"),
    suffix = c("_ecoD", "_ecoR")
  ) %>%
  left_join(
    accuracy_by_ecoprovince_atomic,
    by = c("raster", "Ecprv")
    #suffix = c("", "_ecoP")
  )

# check names
names(fallback_joined_atomic)

fallback_joined_grouped <- accuracy_by_ecodistrict_grouped2 %>%
  left_join(
    accuracy_by_ecoregion_grouped,
    by = c("raster", "Ecrgn_n"),
    suffix = c("_ecoD", "_ecoR")
  ) %>%
  left_join(
    accuracy_by_ecoprovince_grouped,
    by = c("raster", "Ecprv")
    #suffix = c("", "_ecoP")
  )

names(fallback_joined_grouped) %>%
  grep("N_evaluated", ., value = TRUE)

#---- Apply fallback logic ----
# Replace metrics when N_evaluated is too small
accuracy_with_fallback_atomic <- fallback_joined_atomic %>%
  mutate(
    source_level = case_when(
      N_evaluated_ecoD >= min_n ~ "ecodistrict",
      N_evaluated_ecoR >= min_n ~ "ecoregion",
      TRUE ~ "ecoprovince"),
    TP = case_when(
      source_level == "ecodistrict" ~ TP_ecoD,
      source_level == "ecoregion" ~ TP_ecoR,
      TRUE ~ TP_ecoP),
    TN = case_when(
      source_level == "ecodistrict" ~ TN_ecoD,
      source_level == "ecoregion" ~ TN_ecoR,
      TRUE ~ TN_ecoP),
    FP = case_when(
      source_level == "ecodistrict" ~ FP_ecoD,
      source_level == "ecoregion" ~ FP_ecoR,
      TRUE ~ FP_ecoP),
    FN = case_when(
      source_level == "ecodistrict" ~ FN_ecoD,
      source_level == "ecoregion" ~ FN_ecoR,
      TRUE ~ FN_ecoP),
    N_evaluated = case_when(
      source_level == "ecodistrict" ~ N_evaluated_ecoD,
      source_level == "ecoregion" ~ N_evaluated_ecoR,
      TRUE ~ N_evaluated_ecoP),
    # recompute metrics from chosen counts
    metrics = pmap(
      list(TP, TN, FP, FN),
      ~ compute_metrics(..1, ..2, ..3, ..4))) %>%
  tidyr::unnest_wider(metrics) %>%
  select(
    raster, Ecds_ID, Ecrgn_n, Ecprv,
    source_level, everything())

accuracy_with_fallback_grouped <- fallback_joined_grouped %>%
  mutate(
    source_level = case_when(
      !is.na(N_evaluated_ecoD) & N_evaluated_ecoD >= min_n ~ "ecodistrict",
      !is.na(N_evaluated_ecoR) & N_evaluated_ecoR >= min_n ~ "ecoregion",
      TRUE ~ "ecoprovince"),
    TP = case_when(
      source_level == "ecodistrict" ~ TP_ecoD,
      source_level == "ecoregion" ~ TP_ecoR,
      TRUE ~ TP_ecoP),
    TN = case_when(
      source_level == "ecodistrict" ~ TN_ecoD,
      source_level == "ecoregion" ~ TN_ecoR,
      TRUE ~ TN_ecoP),
    FP = case_when(
      source_level == "ecodistrict" ~ FP_ecoD,
      source_level == "ecoregion" ~ FP_ecoR,
      TRUE ~ FP_ecoP),
    FN = case_when(
      source_level == "ecodistrict" ~ FN_ecoD,
      source_level == "ecoregion" ~ FN_ecoR,
      TRUE ~ FN_ecoP),
    N_evaluated = case_when(
      source_level == "ecodistrict" ~ N_evaluated_ecoD,
      source_level == "ecoregion" ~ N_evaluated_ecoR,
      TRUE ~ N_evaluated_ecoP),
    # recompute metrics from chosen counts
    metrics = pmap(
      list(TP, TN, FP, FN),
      ~ compute_metrics(..1, ..2, ..3, ..4))) %>%
  tidyr::unnest_wider(metrics) %>%
  select(
    raster, Ecds_ID, Ecrgn_n, Ecprv,
    source_level, everything())

#---- Export ----

# GROUPED
readr::write_csv(
  accuracy_with_fallback_grouped,
  file.path(out_dir, "binary_accuracy_ecoregprovfallback_grouped_combinedCIandCM_28Jul2026.csv")
)

# ATOMIC
readr::write_csv(
  accuracy_with_fallback_atomic,
  file.path(out_dir, "binary_accuracy_ecoregprovfallback_atomic_combinedCIandCM_28Jul2026.csv")
)

# ----------------------------
# 15) Bootstrapping 
# ----------------------------
# Define bootstrap function
boot_f1 <- function(df, n_boot = 1000) {
  
  vals <- replicate(n_boot, {
    
    samp <- df[sample(
      nrow(df),
      replace = TRUE
    ), ]
    
    counts <- compute_counts(samp)
    
    compute_metrics(
      counts$TP,
      counts$TN,
      counts$FP,
      counts$FN
    )$F1_grass
    
  })
  
  vals <- vals[!is.na(vals)]
  
  if (length(vals) == 0) {
    return(
      tibble(
        F1_boot = NA_real_,
        F1_lwr = NA_real_,
        F1_upr = NA_real_,
        F1_width = NA_real_
      )
    )
  }
  
  tibble(
    F1_boot = mean(vals),
    F1_lwr = quantile(vals, 0.025),
    F1_upr = quantile(vals, 0.975),
    F1_width = quantile(vals, 0.975) -
      quantile(vals, 0.025)
  )
}

# Create bootstrap dataset corresponding to the 
# level actually used by each row
boot_results <- accuracy_with_fallback_grouped %>%
  select(
    raster,
    Ecds_ID,
    Ecrgn_n,
    Ecprv,
    source_level
  ) %>%
  rowwise() %>%
  mutate(
    boot = list({
      
      if (source_level == "ecodistrict") {
        
        dat <- combined %>%
          mutate(raster = raster_group) %>%
          filter(
            raster == cur_data()$raster,
            Ecds_ID == cur_data()$Ecds_ID
          )
        
      } else if (source_level == "ecoregion") {
        
        dat <- combined %>%
          mutate(raster = raster_group) %>%
          filter(
            raster == cur_data()$raster,
            Ecrgn_n == cur_data()$Ecrgn_n
          )
        
      } else {
        
        dat <- combined %>%
          mutate(raster = raster_group) %>%
          filter(
            raster == cur_data()$raster,
            Ecprv == cur_data()$Ecprv
          )
        
      }
      
      boot_f1(dat, n_boot = 1000)
      
    })
  ) %>%
  unnest_wider(boot)

# Join results back to table
accuracy_with_fallback_grouped <- accuracy_with_fallback_grouped %>%
  left_join(
    boot_results %>%
      select(
        raster,
        Ecds_ID,
        F1_boot,
        F1_lwr,
        F1_upr,
        F1_width
      ),
    by = c("raster", "Ecds_ID")
  )

# ----------------------------
# 16) Diagnostics
# ----------------------------

#-------- Count number of NAs in each metic --------#

# How often is MCC NA at each level?

fallback_joined_grouped %>%
  summarise(
    n_rows = n(),
    MCC_ecoD_NA = sum(is.na(MCC_ecoD)),
    MCC_ecoR_NA = sum(is.na(MCC_ecoR)),
    MCC_ecoP_NA = sum(is.na(MCC_ecoP)),
    pct_MCC_ecoD_NA = mean(is.na(MCC_ecoD)) * 100,
    pct_MCC_ecoR_NA = mean(is.na(MCC_ecoR)) * 100,
    pct_MCC_ecoP_NA = mean(is.na(MCC_ecoP)) * 100
  )

# without ACI2025:
#n_rows MCC_ecoD_NA MCC_ecoR_NA MCC_ecoP_NA pct_MCC_ecoD_NA pct_MCC_ecoR_NA pct_MCC_ecoP_NA
#  3927        1614        1035         885            41.1            26.4            22.5

# With ACI2025:
# n_rows MCC_ecoD_NA MCC_ecoR_NA MCC_ecoP_NA pct_MCC_ecoD_NA pct_MCC_ecoR_NA pct_MCC_ecoP_NA
#   4015        1731        1149         998            43.1            28.6            24.9

# exact match (not lenient):
# n_rows MCC_ecoD_NA MCC_ecoR_NA MCC_ecoP_NA pct_MCC_ecoD_NA pct_MCC_ecoR_NA pct_MCC_ecoP_NA
#   3213        1996        1431        1258            62.1            44.5            39.2

fallback_joined_atomic %>%
  summarise(
    n_rows = n(),
    MCC_ecoD_NA = sum(is.na(MCC_ecoD)),
    MCC_ecoR_NA = sum(is.na(MCC_ecoR)),
    MCC_ecoP_NA = sum(is.na(MCC_ecoP)),
    pct_MCC_ecoD_NA = mean(is.na(MCC_ecoD)) * 100,
    pct_MCC_ecoR_NA = mean(is.na(MCC_ecoR)) * 100,
    pct_MCC_ecoP_NA = mean(is.na(MCC_ecoP)) * 100
  )

# without ACI2025:
#n_rows MCC_ecoD_NA MCC_ecoR_NA MCC_ecoP_NA pct_MCC_ecoD_NA pct_MCC_ecoR_NA pct_MCC_ecoP_NA
#  3929        1617        1057         908            41.2            26.9            23.1

# with ACI2025:
# n_rows MCC_ecoD_NA MCC_ecoR_NA MCC_ecoP_NA pct_MCC_ecoD_NA pct_MCC_ecoR_NA pct_MCC_ecoP_NA
#   4015        1731        1149        1000            43.1            28.6            24.9

# exact match (not lenient):
# n_rows MCC_ecoD_NA MCC_ecoR_NA MCC_ecoP_NA pct_MCC_ecoD_NA pct_MCC_ecoR_NA pct_MCC_ecoP_NA
#   3213        1996        1431        1285            62.1            44.5            40.0

# How often is overall accuracy NA at each level?

fallback_joined_grouped %>%
  summarise(
    n_rows = n(),
    OA_ecoD_NA = sum(is.na(Overall_Accuracy_ecoD)),
    OA_ecoR_NA = sum(is.na(Overall_Accuracy_ecoR)),
    OA_ecoP_NA = sum(is.na(Overall_Accuracy_ecoP)),
    pct_OA_ecoD_NA = mean(is.na(Overall_Accuracy_ecoD)) * 100,
    pct_OA_ecoR_NA = mean(is.na(Overall_Accuracy_ecoR)) * 100,
    pct_OA_ecoP_NA = mean(is.na(Overall_Accuracy_ecoP)) * 100
  )

# without ACI2025:
#n_rows OA_ecoD_NA OA_ecoR_NA OA_ecoP_NA pct_OA_ecoD_NA pct_OA_ecoR_NA pct_OA_ecoP_NA
#  3927          0          0          0              0              0              0

# with ACI2025:
# n_rows OA_ecoD_NA OA_ecoR_NA OA_ecoP_NA pct_OA_ecoD_NA pct_OA_ecoR_NA pct_OA_ecoP_NA
#   4015          0          0          0              0              0              0

# exact match (not lenient):
# n_rows OA_ecoD_NA OA_ecoR_NA OA_ecoP_NA pct_OA_ecoD_NA pct_OA_ecoR_NA pct_OA_ecoP_NA
#   3213          0          0          0              0              0              0

fallback_joined_atomic %>%
  summarise(
    n_rows = n(),
    OA_ecoD_NA = sum(is.na(Overall_Accuracy_ecoD)),
    OA_ecoR_NA = sum(is.na(Overall_Accuracy_ecoR)),
    OA_ecoP_NA = sum(is.na(Overall_Accuracy_ecoP)),
    pct_OA_ecoD_NA = mean(is.na(Overall_Accuracy_ecoD)) * 100,
    pct_OA_ecoR_NA = mean(is.na(Overall_Accuracy_ecoR)) * 100,
    pct_OA_ecoP_NA = mean(is.na(Overall_Accuracy_ecoP)) * 100
  )

# without ACI2025:
#n_rows OA_ecoD_NA OA_ecoR_NA OA_ecoP_NA pct_OA_ecoD_NA pct_OA_ecoR_NA pct_OA_ecoP_NA
#  3929          0          0          0              0              0              0

# with ACI2025:
# n_rows OA_ecoD_NA OA_ecoR_NA OA_ecoP_NA pct_OA_ecoD_NA pct_OA_ecoR_NA pct_OA_ecoP_NA
#   4015          0          0          0              0              0              0

# How often is balanced accuracy NA at each level?

fallback_joined_grouped %>%
  summarise(
    n_rows = n(),
    BA_ecoD_NA = sum(is.na(Balanced_Accuracy_ecoD)),
    BA_ecoR_NA = sum(is.na(Balanced_Accuracy_ecoR)),
    BA_ecoP_NA = sum(is.na(Balanced_Accuracy_ecoP)),
    pct_BA_ecoD_NA = mean(is.na(Balanced_Accuracy_ecoD)) * 100,
    pct_BA_ecoR_NA = mean(is.na(Balanced_Accuracy_ecoR)) * 100,
    pct_BA_ecoP_NA = mean(is.na(Balanced_Accuracy_ecoP)) * 100
  )

# without ACI2025:
#n_rows BA_ecoD_NA BA_ecoR_NA BA_ecoP_NA pct_BA_ecoD_NA pct_BA_ecoR_NA pct_BA_ecoP_NA
#  3927       1399        976        882           35.6           24.9           22.5

# with ACI2025:
# n_rows BA_ecoD_NA BA_ecoR_NA BA_ecoP_NA pct_BA_ecoD_NA pct_BA_ecoR_NA pct_BA_ecoP_NA
#   4015       1515       1090        995           37.7           27.1           24.8

# exact match (not lenient):
# n_rows BA_ecoD_NA BA_ecoR_NA BA_ecoP_NA pct_BA_ecoD_NA pct_BA_ecoR_NA pct_BA_ecoP_NA
#   3213       1801       1381       1258           56.1           43.0           39.2

fallback_joined_atomic %>%
  summarise(
    n_rows = n(),
    BA_ecoD_NA = sum(is.na(Balanced_Accuracy_ecoD)),
    BA_ecoR_NA = sum(is.na(Balanced_Accuracy_ecoR)),
    BA_ecoP_NA = sum(is.na(Balanced_Accuracy_ecoP)),
    pct_BA_ecoD_NA = mean(is.na(Balanced_Accuracy_ecoD)) * 100,
    pct_BA_ecoR_NA = mean(is.na(Balanced_Accuracy_ecoR)) * 100,
    pct_BA_ecoP_NA = mean(is.na(Balanced_Accuracy_ecoP)) * 100
  )

# without ACI2025:
#n_rows BA_ecoD_NA BA_ecoR_NA BA_ecoP_NA pct_BA_ecoD_NA pct_BA_ecoR_NA pct_BA_ecoP_NA
#  3929       1402        998        905           35.7           25.4           23.0

# with ACI2025:
# n_rows BA_ecoD_NA BA_ecoR_NA BA_ecoP_NA pct_BA_ecoD_NA pct_BA_ecoR_NA pct_BA_ecoP_NA
#   4015       1515       1090        997           37.7           27.1           24.8

# Using fallbacks, how often are MCC and OA na?
accuracy_with_fallback_grouped %>%
  summarise(
    n_rows = n(),
    OA_NA  = sum(is.na(Overall_Accuracy)),
    OA_pct_NA = 100 * mean(is.na(Overall_Accuracy)),
    F1_NA = sum(is.na(F1_grass)),
    F1_pct_NA = 100 * mean(is.na(F1_grass)),
    BA_NA  = sum(is.na(Balanced_Accuracy)),
    BA_pct_NA = 100 * mean(is.na(Balanced_Accuracy)),
    MCC_NA = sum(is.na(MCC)),
    MCC_pct_NA = 100 * mean(is.na(MCC))
  )

# without ACI2025:
#n_rows OA_NA OA_pct_NA F1_NA F1_pct_NA BA_NA BA_pct_NA MCC_NA MCC_pct_NA
#  3927     0         0   376      9.57  1114      28.4   1312       33.4

# with ACI2025:
# n_rows OA_NA OA_pct_NA F1_NA F1_pct_NA BA_NA BA_pct_NA MCC_NA MCC_pct_NA
#   4015     0         0   423      10.5  1223      30.5   1421       35.4

# exact match (not lenient):
# n_rows OA_NA OA_pct_NA F1_NA F1_pct_NA BA_NA BA_pct_NA MCC_NA MCC_pct_NA
#   3213     0         0   272      8.47  1575      49.0   1754       54.6

accuracy_with_fallback_atomic %>%
  summarise(
    n_rows = n(),
    OA_NA  = sum(is.na(Overall_Accuracy)),
    OA_pct_NA = 100 * mean(is.na(Overall_Accuracy)),
    F1_NA = sum(is.na(F1_grass)),
    F1_pct_NA = 100 * mean(is.na(F1_grass)),
    BA_NA  = sum(is.na(Balanced_Accuracy)),
    BA_pct_NA = 100 * mean(is.na(Balanced_Accuracy)),
    MCC_NA = sum(is.na(MCC)),
    MCC_pct_NA = 100 * mean(is.na(MCC))
  )

# without ACI2025:
#n_rows OA_NA OA_pct_NA F1_NA F1_pct_NA BA_NA BA_pct_NA MCC_NA MCC_pct_NA
#  3929     0         0   377      9.60  1127      28.7   1325       33.7

# with ACI2025:
# n_rows OA_NA OA_pct_NA F1_NA F1_pct_NA BA_NA BA_pct_NA MCC_NA MCC_pct_NA
#   4015     0         0   423      10.5  1223      30.5   1421       35.4

#------ Why are there more rows for grouped than atomic? ------#

# Check row counts after every step
list(
  accuracy_by_ecodistrict_atomic  = nrow(accuracy_by_ecodistrict_atomic),
  accuracy_by_ecodistrict_grouped = nrow(accuracy_by_ecodistrict_grouped),
  
  accuracy_by_ecodistrict_atomic2  = nrow(accuracy_by_ecodistrict_atomic2),
  accuracy_by_ecodistrict_grouped2 = nrow(accuracy_by_ecodistrict_grouped2),
  
  fallback_joined_atomic  = nrow(fallback_joined_atomic),
  fallback_joined_grouped = nrow(fallback_joined_grouped)
)

# Check whether grouped contains duplicate raster × ecodistrict combinations
accuracy_by_ecodistrict_grouped %>%
  count(raster, Ecds_ID) %>%
  filter(n > 1)

# Check whether duplicates are introduced after joining ecoregions/ecoprovinces
dup_keys <- fallback_joined_grouped %>%
  count(raster, Ecds_ID) %>%
  filter(n > 1)

fallback_joined_grouped %>%
  semi_join(dup_keys, by = c("raster", "Ecds_ID")) %>%
  arrange(raster, Ecds_ID)

# Check for duplicated keys in ecodistrict mapping
ecodist_to_regions %>%
  count(Ecds_ID) %>%
  filter(n > 1)

ecodist_to_regions %>%
  count(Ecds_ID, Ecrgn_n, Ecprv) %>%
  filter(n > 1)

# Should be equal:
nrow(ecodist_to_regions)
n_distinct(ecodist_to_regions$Ecds_ID)

# Check ecoregion table for duplicated join keys (should be 0)
accuracy_by_ecoregion_grouped %>%
  count(raster, Ecrgn_n) %>%
  filter(n > 1)

accuracy_by_ecoregion_atomic %>%
  count(raster, Ecrgn_n) %>%
  filter(n > 1)

accuracy_by_ecoprovince_grouped %>%
  count(raster, Ecprv) %>%
  filter(n > 1)

accuracy_by_ecoprovince_grouped %>%
  count(raster, Ecprv) %>%
  filter(n > 1)

# Use anti-joins to find rows present only in grouped
group_keys <- fallback_joined_grouped %>%
  distinct(raster, Ecds_ID)

atomic_keys <- fallback_joined_atomic %>%
  distinct(raster, Ecds_ID)

group_keys %>%
  anti_join(atomic_keys,
            by = c("raster", "Ecds_ID")) %>%
  nrow()

# Compare unique key counts
# If rows > unique keys then duplicates exist
# If rows == unique keys then no duplicates (grouped genuinely contains more raster–ecodistrict combinations)
cat(
  "Grouped rows:", nrow(fallback_joined_grouped), "\n",
  "Grouped unique raster-Ecds_ID:",
  fallback_joined_grouped %>%
    distinct(raster, Ecds_ID) %>%
    nrow(), "\n\n",
  
  "Atomic rows:", nrow(fallback_joined_atomic), "\n",
  "Atomic unique raster-Ecds_ID:",
  fallback_joined_atomic %>%
    distinct(raster, Ecds_ID) %>%
    nrow(), "\n"
)

# compare ecodistrict keys directly
atomic_keys <- accuracy_by_ecodistrict_atomic %>%
  distinct(raster, Ecds_ID)

group_keys <- accuracy_by_ecodistrict_grouped %>%
  distinct(raster, Ecds_ID)

group_keys %>%
  anti_join(
    atomic_keys,
    by = c("raster", "Ecds_ID")
  ) %>%
  arrange(raster, Ecds_ID) %>%
  print(n = Inf)

# Check to see if 67 "extra" grouped keys are  comparison artifact
atomic_keys_grouped <- accuracy_by_ecodistrict_atomic %>%
  mutate(
    raster = if_else(
      raster %in% c("pli_ap","pli_cu","pli_mg","pli_mm"),
      "pli",
      raster
    )
  ) %>%
  distinct(raster, Ecds_ID)

group_keys <- accuracy_by_ecodistrict_grouped %>%
  distinct(raster, Ecds_ID)

anti_join(
  group_keys,
  atomic_keys_grouped,
  by = c("raster", "Ecds_ID")
)

# check raster to group
table(combined$raster_group, useNA = "ifany")

# old
# aci2024     cec     cgr    cnvc      dw    eccc     esa    esri     gcc     glc     gpw     lcc     lcr 
#  271060   74117   28144   39971  268285   60363   62429  276256    1400   38962  275777   74117   74117 
# luts10  luts30     mgi     ngi    pglc     pli    plvi    pudl 
#  74117   74117    7336   24479   50713   30589     602   37736 

# new - added ACI2025, removed duplicate points
# aci2024 aci2025     cec     cgr    cnvc      dw    eccc     esa    esri     gcc     glc     gpw     lcc 
#  265807    5884   73879   28144   39733  268285   60363   62191  271016    1167   38724  270552   73879 
#   lcr  luts10  luts30     mgi     ngi    pglc     pli    plvi    pudl 
# 73879   73879   73879    7336   24479   50713   30587     602   37736 

# exact matching (not lenient)
# aci2024 aci2025     cec     cgr    cnvc      dw    eccc     esa    esri     gcc     glc     gpw     lcc 
# 610     869   73879   28144   39733  263283   27340   29109    5850    1167   31153    5808   73879 
# lcr  luts10  luts30     mgi     ngi    pglc     pli    plvi    pudl 
# 73879   73879   73879    7336   24479   50713   23660     602   37736 

# check atomic
table(combined$raster, useNA = "ifany")

# aci2024 aci2025     cec     cgr    cnvc      dw    eccc     esa    esri     gcc     glc     gpw     lcc 
# 265807    5884   73879   28144   39733  268285   60363   62191  271016    1167   38724  270552   73879 
# lcr  luts10  luts30     mgi     ngi    pglc  pli_ap  pli_cu  pli_mg  pli_mm    plvi    pudl 
# 73879   73879   73879    7336   24479   50713     710     920   20638    8319     602   37736 

# exact match:
# aci2024 aci2025     cec     cgr    cnvc      dw    eccc     esa    esri     gcc     glc     gpw     lcc 
# 610     869   73879   28144   39733  263283   27340   29109    5850    1167   31153    5808   73879 
# lcr  luts10  luts30     mgi     ngi    pglc  pli_ap  pli_cu  pli_mg  pli_mm    plvi    pudl 
# 73879   73879   73879    7336   24479   50713     710     920   13711    8319     602   37736 

combined %>%
  count(raster_group) %>%
  ggplot(aes(x = reorder(raster_group, n), y = n)) +
  geom_col() +
  coord_flip() +
  scale_y_continuous(labels = scales::label_comma()) +
  labs(
    x = "Dataset",
    y = "Count of validation points",
    title = "Number of validation points by dataset"
  ) +
  theme_bw()

unique(combined$raster_group)

combined %>%
  filter(raster_group == "aci") %>%
  count(raster, sort = TRUE)

combined %>%
  filter(raster_group == "aci") %>%
  distinct(raster)

# Where are missing values?
accuracy_with_fallback_grouped %>%
  summarise(
    no_grass =
      sum(TP + FN == 0, na.rm = TRUE),
    
    no_notgrass =
      sum(TN + FP == 0, na.rm = TRUE),
    
    no_positive_predictions =
      sum(TP + FP == 0, na.rm = TRUE),
    
    no_negative_predictions =
      sum(TN + FN == 0, na.rm = TRUE),
    
    BA_NA =
      sum(is.na(Balanced_Accuracy)),
    
    MCC_NA =
      sum(is.na(MCC))
  )

# without ACI2025:
#no_grass no_notgrass no_positive_predictions no_negative_predictions BA_NA MCC_NA
#      87        1027                     283                      18  1114   1312

# with ACI2025:
# no_grass no_notgrass no_positive_predictions no_negative_predictions BA_NA MCC_NA
#      132        1091                    287                      19  1223   1421

# exact match (not lenient):
# no_grass no_notgrass no_positive_predictions no_negative_predictions BA_NA MCC_NA
#       12        1563                     259                      28  1575   1754

# Why different numbers in PLI grouped vs atomic? - this is fixed
atomic_keys_grouped <- accuracy_by_ecodistrict_atomic %>%
  mutate(
    raster = if_else(
      raster %in% c("pli_ap", "pli_cu", "pli_mg", "pli_mm"),
      "pli",
      raster
    )
  ) %>%
  distinct(raster, Ecds_ID)

group_keys <- accuracy_by_ecodistrict_grouped %>%
  distinct(raster, Ecds_ID)

atomic_keys_grouped %>%
  anti_join(
    group_keys,
    by = c("raster", "Ecds_ID")
  )

nrow(accuracy_by_ecodistrict_atomic)
nrow(accuracy_by_ecodistrict_grouped)

nrow(atomic_keys_grouped)
nrow(group_keys)

# Difference is occurring because two atomic rows collapse onto existing rows when recoding
# pli_ap, pli_cu, etc to pli

accuracy_by_ecodistrict_atomic %>%
  mutate(
    raster = if_else(
      raster %in% c("pli_ap","pli_cu","pli_mg","pli_mm"),
      "pli",
      raster
    )
  ) %>%
  count(raster, Ecds_ID) %>%
  filter(n > 1)

# in atomic table, there are 2 ecodistricts where more than one PLI variant occurs:
# raster Ecds_ID     n
# pli        745     2
# pli        830     2

# Ecodistrict 745 should be in AP
# Ecodistrict 830 should be in MG

accuracy_by_ecodistrict_atomic %>%
  filter(Ecds_ID %in% c(745, 830)) %>%
  arrange(Ecds_ID, raster) %>%
  select(raster, Ecds_ID, N_evaluated)

combined %>%
  filter(Ecds_ID %in% c(745, 830)) %>%
  count(Ecds_ID, raster)

accuracy_by_ecodistrict_atomic %>%
  filter(Ecds_ID %in% c(745, 830)) %>%
  arrange(Ecds_ID, raster) %>%
  print(n=Inf)

#   raster  Ecds_ID    TP    TN    FP    FN N_evaluated Overall_Accuracy Precision_grass Recall_grass
# 1 aci2024     745   191  7398    57    94        7740            0.980           0.770        0.670
# 2 cec         745   192   205     0   142         539            0.737           1            0.575
# 3 cgr         745    59    60    17    91         227            0.524           0.776        0.393
# 4 cnvc        745   145     0     0    80         225            0.644           1            0.644
# 5 dw          745   203  6554   901    81        7739            0.873           0.184        0.715
# 6 eccc        745   134   205     1    49         389            0.871           0.993        0.732
# 7 esa         745   144   198     8    39         389            0.879           0.947        0.787
# 8 esri        745   114  7157   298   171        7740            0.939           0.277        0.4  
# 9 glc         745   117    64    13    33         227            0.797           0.9          0.78 
# 10 gpw        745   204  7333   122    81        7740            0.974           0.626        0.716
# 11 lcc        745   191   205     0   143         539            0.735           1            0.572
# 12 lcr        745   149   205     0   185         539            0.657           1            0.446
# 13 luts10     745   140   205     0   194         539            0.640           1            0.419
# 14 luts30     745   143   205     0   191         539            0.646           1            0.428
# 15 ngi        745    74     0     1    12          87            0.851           0.987        0.860
# 16 pglc       745     0     0     0     3           3            0              NA            0    
# 17 pli_ap     745    63     0     0    19          82            0.768           1            0.768
# 18 pli_mm     745     0     1     0     0           1            1              NA           NA    
# 19 pudl       745   102     0     0   123         225            0.453           1            0.453
# 20 aci2024    830   534  1900    54    74        2562            0.950           0.908        0.878
# 21 aci2025    830    65     0     0    12          77            0.844           1            0.844
# 22 cec        830   803   109     1   362        1275            0.715           0.999        0.689
# 23 cgr        830   517     0     0    20         537            0.963           1            0.963
# 24 cnvc       830   644     0     0   148         792            0.813           1            0.813
# 25 dw         830   473  1810   144   165        2592            0.881           0.767        0.741
# 26 eccc       830   776   106     4    49         935            0.943           0.995        0.941
# 27 esa        830   785   102     8    40         935            0.949           0.990        0.952
# 28 esri       830   461  1591   363   147        2562            0.801           0.559        0.758
# 29 glc        830   499     0     0    38         537            0.929           1            0.929
# 30 gpw        830   522  1798   156    86        2562            0.906           0.770        0.859
# 31 lcc        830   798   110     0   367        1275            0.712           1            0.685
# 32 lcr        830   912   110     0   253        1275            0.802           1            0.783
# 33 luts10     830   863   110     0   302        1275            0.763           1            0.741
# 34 luts30     830   877   110     0   288        1275            0.774           1            0.753
# 35 ngi        830   508     0     0    29         537            0.946           1            0.946
# 36 pglc       830  1036   105     5   129        1275            0.895           0.995        0.889
# 37 pli_mg     830   875    88    21   204        1188            0.811           0.977        0.811
# 38 pli_mm     830     0     1     0     0           1            1              NA           NA    
# 39 pudl       830   652     0     0   140         792            0.823           1            0.823

#------ Check correlations among metrics ------#

accuracy_with_fallback_grouped %>%
  select(
    Overall_Accuracy,
    Balanced_Accuracy,
    F1_grass,
    MCC
  ) %>%
  cor(use = "pairwise.complete.obs")

# without ACI2025:
#                   Overall_Accuracy Balanced_Accuracy  F1_grass       MCC
# Overall_Accuracy         1.0000000         0.6874857 0.7184568 0.7073522
# Balanced_Accuracy        0.6874857         1.0000000 0.7667255 0.8115127
# F1_grass                 0.7184568         0.7667255 1.0000000 0.6402146
# MCC                      0.7073522         0.8115127 0.6402146 1.0000000

# with ACI2025:
#                   Overall_Accuracy Balanced_Accuracy  F1_grass       MCC
# Overall_Accuracy         1.0000000         0.6732306 0.7075173 0.6917296
# Balanced_Accuracy        0.6732306         1.0000000 0.7481058 0.8216997
# F1_grass                 0.7075173         0.7481058 1.0000000 0.6261291
# MCC                      0.6917296         0.8216997 0.6261291 1.0000000

# exact match (not lenient):
#                   Overall_Accuracy Balanced_Accuracy  F1_grass       MCC
# Overall_Accuracy         1.0000000         0.6613476 0.8386212 0.6708037
# Balanced_Accuracy        0.6613476         1.0000000 0.6780332 0.8011216
# F1_grass                 0.8386212         0.6780332 1.0000000 0.5205585
# MCC                      0.6708037         0.8011216 0.5205585 1.0000000

#------ Check distribution of number of points ------#

accuracy_with_fallback_grouped %>%
  ggplot(aes(x = N_evaluated)) +
  geom_histogram(binwidth = 5, color = "black", fill = "steelblue") +
  geom_vline(
    xintercept = 20,
    linetype = "dashed",
    color = "red",
    linewidth = 1
  ) +
  labs(
    title = "Distribution of Validation Sample Sizes",
    x = "n evaluated",
    y = "Count of ecodistricts"
  ) +
  theme_bw()

# highly right-skewed - trying log-scaled x-axis
accuracy_with_fallback_grouped %>%
  ggplot(aes(x = N_evaluated)) +
  geom_histogram(bins = 50, color = "black", fill = "steelblue") +
  scale_x_log10() +
  labs(
    title = "Distribution of Validation Sample Sizes",
    x = "n evaluated (log10 scale)",
    y = "Count of ecodistricts"
  ) +
  theme_bw()

# Number of rows with N less than 20
accuracy_with_fallback_grouped %>%
  summarise(
    n_below_20 = sum(N_evaluated < 20),
    pct_below_20 = 100 * mean(N_evaluated < 20)
  )

# without ACI2025:
# n_below_20 pct_below_20
#         25        0.637

# with ACI2025:
# n_below_20 pct_below_20
#         40        0.996

# exact match (not lenient):
# n_below_20 pct_below_20
#         49         1.53

accuracy_with_fallback_atomic %>%
  summarise(
    n_below_20 = sum(N_evaluated < 20),
    pct_below_20 = 100 * mean(N_evaluated < 20)
  )

# without ACI2025:
# n_below_20 pct_below_20
#         26        0.662

# with ACI2025:
# n_below_20 pct_below_20
#         40        0.996

# More info on rows with N less than 20 (atomic)
accuracy_with_fallback_atomic %>%
  filter(N_evaluated < 20) %>%
  select(
    raster,
    Ecds_ID,
    Ecrgn_n,
    Ecprv,
    source_level,
    N_evaluated,
    F1_grass,
    Overall_Accuracy,
    TP, TN, FP, FN
  ) %>%
  arrange(N_evaluated, raster, Ecds_ID) %>%
  print(n = Inf)

# without ACI2025:
#   raster Ecds_ID Ecrgn_n   Ecprv source_level N_evaluated F1_grass Overall_Accuracy    TP    TN    FP    FN
# 1 glc        374 Lac Seul~   6.2 ecoprovince            1   NA               1          0     1     0     0
# 2 ngi        672 Mid-Bore~   9.3 ecoprovince            1   NA               0          0     0     0     1
# 3 ngi        850 Lake Man~  10.1 ecoprovince            1    1               1          1     0     0     0
# 4 pli_mm     745 Aspen Pa~  10.2 ecoprovince            1   NA               1          0     1     0     0
# 5 cec        960 Cascade ~  13.2 ecoprovince            2   NA               0          0     0     0     2
# 6 cec        982 Fraser B~  14.1 ecoprovince            2   NA               0          0     0     0     2
# 7 cnvc       960 Cascade ~  13.2 ecoprovince            2   NA               0          0     0     0     2
# 8 esa        960 Cascade ~  13.2 ecoprovince            2    0.667           0.5        1     0     0     1
# 9 esa        982 Fraser B~  14.1 ecoprovince            2    1               1          2     0     0     0
# 10 glc        982 Fraser B~  14.1 ecoprovince            2   NA               0          0     0     0     2
# 11 lcc        960 Cascade ~  13.2 ecoprovince            2   NA               0          0     0     0     2
# 12 lcc        982 Fraser B~  14.1 ecoprovince            2    0.667           0.5        1     0     0     1
# 13 lcr        960 Cascade ~  13.2 ecoprovince            2   NA               0          0     0     0     2
# 14 lcr        982 Fraser B~  14.1 ecoprovince            2   NA               0          0     0     0     2
# 15 luts10     960 Cascade ~  13.2 ecoprovince            2   NA               0          0     0     0     2
# 16 luts10     982 Fraser B~  14.1 ecoprovince            2   NA               0          0     0     0     2
# 17 luts30     960 Cascade ~  13.2 ecoprovince            2   NA               0          0     0     0     2
# 18 luts30     982 Fraser B~  14.1 ecoprovince            2   NA               0          0     0     0     2
# 19 gcc       1015 Southern~  14.4 ecoprovince            6    1               1          6     0     0     0
# 20 cgr        375 Lake of ~   6.5 ecoprovince           11    0.167           0.0909     1     0     0    10
# 21 cgr        687 Boreal T~   9.2 ecoprovince           11   NA               0          0     0     0    11
# 22 cgr        688 Boreal T~   9.2 ecoprovince           11   NA               0          0     0     0    11
# 23 cgr        695 Boreal T~   9.2 ecoprovince           11   NA               0          0     0     0    11
# 24 cgr        701 Boreal T~   9.2 ecoprovince           11   NA               0          0     0     0    11
# 25 cgr        709 Boreal T~   9.2 ecoprovince           11   NA               0          0     0     0    11
# 26 pudl       375 Lake of ~   6.5 ecoprovince           11    0.533           0.364      4     0     0     7

# with ACI2025:
#   raster   Ecds_ID Ecrgn_n  Ecprv source_level N_evaluated F1_grass Overall_Accuracy    TP    TN    FP    FN
# 1 aci2024      374 Lac Seu~   6.2 ecoprovince            1   NA               0          0     0     1     0
# 2 aci2025      960 Cascade~  13.2 ecoprovince            1   NA               0          0     0     0     1
# 3 glc          374 Lac Seu~   6.2 ecoprovince            1   NA               1          0     1     0     0
# 4 ngi          672 Mid-Bor~   9.3 ecoprovince            1   NA               0          0     0     0     1
# 5 ngi          850 Lake Ma~  10.1 ecoprovince            1    1               1          1     0     0     0
# 6 cec          960 Cascade~  13.2 ecoprovince            2   NA               0          0     0     0     2
# 7 cec          982 Fraser ~  14.1 ecoprovince            2   NA               0          0     0     0     2
# 8 cnvc         960 Cascade~  13.2 ecoprovince            2   NA               0          0     0     0     2
# 9 esa          960 Cascade~  13.2 ecoprovince            2    0.667           0.5        1     0     0     1
# 10 esa         982 Fraser ~  14.1 ecoprovince            2    1               1          2     0     0     0
# 11 glc         982 Fraser ~  14.1 ecoprovince            2   NA               0          0     0     0     2
# 12 lcc         960 Cascade~  13.2 ecoprovince            2   NA               0          0     0     0     2
# 13 lcc         982 Fraser ~  14.1 ecoprovince            2    0.667           0.5        1     0     0     1
# 14 lcr         960 Cascade~  13.2 ecoprovince            2   NA               0          0     0     0     2
# 15 lcr         982 Fraser ~  14.1 ecoprovince            2   NA               0          0     0     0     2
# 16 luts10      960 Cascade~  13.2 ecoprovince            2   NA               0          0     0     0     2
# 17 luts10      982 Fraser ~  14.1 ecoprovince            2   NA               0          0     0     0     2
# 18 luts30      960 Cascade~  13.2 ecoprovince            2   NA               0          0     0     0     2
# 19 luts30      982 Fraser ~  14.1 ecoprovince            2   NA               0          0     0     0     2
# 20 aci2025     591 Peace L~   9.2 ecoprovince            3   NA               0.667      0     2     1     0
# 21 aci2025     617 Mid-Bor~   9.2 ecoprovince            3   NA               0.667      0     2     1     0
# 22 aci2025     982 Fraser ~  14.1 ecoprovince            4   NA               0.25       0     1     3     0
# 23 aci2025     727 Aspen P~  10.2 ecoprovince            6    0.8             0.667      4     0     0     2
# 24 aci2025     739 Aspen P~  10.2 ecoprovince            6    0.8             0.667      4     0     0     2
# 25 aci2025     750 Aspen P~  10.2 ecoprovince            6    0.8             0.667      4     0     0     2
# 26 aci2025     757 Aspen P~  10.2 ecoprovince            6    0.8             0.667      4     0     0     2
# 27 aci2025     955 Eastern~  13.1 ecoprovince            6   NA               1          0     6     0     0
# 28 aci2025     956 Eastern~  13.1 ecoprovince            6   NA               1          0     6     0     0
# 29 gcc        1015 Souther~  14.4 ecoprovince            6    1               1          6     0     0     0
# 30 aci2025    1017 Norther~  14.4 ecoprovince            7    0.444           0.286      2     0     0     5
# 31 aci2025    1018 Norther~  14.4 ecoprovince            7    0.444           0.286      2     0     0     5
# 32 cgr         375 Lake of~   6.5 ecoprovince           11    0.167           0.0909     1     0     0    10
# 33 cgr         687 Boreal ~   9.2 ecoprovince           11   NA               0          0     0     0    11
# 34 cgr         688 Boreal ~   9.2 ecoprovince           11   NA               0          0     0     0    11
# 35 cgr         695 Boreal ~   9.2 ecoprovince           11   NA               0          0     0     0    11
# 36 cgr         701 Boreal ~   9.2 ecoprovince           11   NA               0          0     0     0    11
# 37 cgr         709 Boreal ~   9.2 ecoprovince           11   NA               0          0     0     0    11
# 38 pudl        375 Lake of~   6.5 ecoprovince           11    0.533           0.364      4     0     0     7
# 39 aci2025     618 Western~   9.1 ecoprovince           13    0.588           0.462      5     1     0     7
# 40 aci2025     631 Western~   9.1 ecoprovince           13    0.588           0.462      5     1     0     7

# exact match (not lenient):
#   raster  Ecds_ID Ecrgn_n  Ecprv source_level N_evaluated F1_grass Overall_Accuracy    TP    TN    FP    FN
# 1 aci2024     850 Lake Ma~  10.1 ecoprovince            1   NA               0          0     0     0     1
# 2 aci2025     960 Cascade~  13.2 ecoprovince            1   NA               0          0     0     0     1
# 3 esri        850 Lake Ma~  10.1 ecoprovince            1   NA               0          0     0     0     1
# 4 gpw         850 Lake Ma~  10.1 ecoprovince            1    1               1          1     0     0     0
# 5 gpw         960 Cascade~  13.2 ecoprovince            1   NA               0          0     0     0     1
# 6 ngi         672 Mid-Bor~   9.3 ecoprovince            1   NA               0          0     0     0     1
# 7 ngi         850 Lake Ma~  10.1 ecoprovince            1    1               1          1     0     0     0
# 8 aci2024     960 Cascade~  13.2 ecoprovince            2   NA               0          0     0     0     2
# 9 aci2024     982 Fraser ~  14.1 ecoprovince            2   NA               0          0     0     0     2
# 10 cec         960 Cascade~  13.2 ecoprovince            2   NA               0          0     0     0     2
# 11 cec         982 Fraser ~  14.1 ecoprovince            2   NA               0          0     0     0     2
# 12 cnvc        960 Cascade~  13.2 ecoprovince            2   NA               0          0     0     0     2
# 13 esa         960 Cascade~  13.2 ecoprovince            2    0.667           0.5        1     0     0     1
# 14 esa         982 Fraser ~  14.1 ecoprovince            2    1               1          2     0     0     0
# 15 esri        960 Cascade~  13.2 ecoprovince            2    1               1          2     0     0     0
# 16 esri        982 Fraser ~  14.1 ecoprovince            2   NA               0          0     0     0     2
# 17 glc         960 Cascade~  13.2 ecoprovince            2   NA               0          0     0     0     2
# 18 glc         982 Fraser ~  14.1 ecoprovince            2   NA               0          0     0     0     2
# 19 gpw         982 Fraser ~  14.1 ecoprovince            2    1               1          2     0     0     0
# 20 lcc         960 Cascade~  13.2 ecoprovince            2   NA               0          0     0     0     2
# 21 lcc         982 Fraser ~  14.1 ecoprovince            2    0.667           0.5        1     0     0     1
# 22 lcr         960 Cascade~  13.2 ecoprovince            2   NA               0          0     0     0     2
# 23 lcr         982 Fraser ~  14.1 ecoprovince            2   NA               0          0     0     0     2
# 24 luts10      960 Cascade~  13.2 ecoprovince            2   NA               0          0     0     0     2
# 25 luts10      982 Fraser ~  14.1 ecoprovince            2   NA               0          0     0     0     2
# 26 luts30      960 Cascade~  13.2 ecoprovince            2   NA               0          0     0     0     2
# 27 luts30      982 Fraser ~  14.1 ecoprovince            2   NA               0          0     0     0     2
# 28 aci2025     591 Peace L~   9.2 ecoprovince            3   NA               0.667      0     2     1     0
# 29 aci2025     617 Mid-Bor~   9.2 ecoprovince            3   NA               0.667      0     2     1     0
# 30 aci2025     982 Fraser ~  14.1 ecoprovince            4   NA               0.25       0     1     3     0
# 31 aci2025     727 Aspen P~  10.2 ecoprovince            6    0.8             0.667      4     0     0     2
# 32 aci2025     739 Aspen P~  10.2 ecoprovince            6    0.8             0.667      4     0     0     2
# 33 aci2025     750 Aspen P~  10.2 ecoprovince            6    0.8             0.667      4     0     0     2
# 34 aci2025     757 Aspen P~  10.2 ecoprovince            6    0.8             0.667      4     0     0     2
# 35 aci2025     955 Eastern~  13.1 ecoprovince            6   NA               1          0     6     0     0
# 36 aci2025     956 Eastern~  13.1 ecoprovince            6   NA               1          0     6     0     0
# 37 gcc        1015 Souther~  14.4 ecoprovince            6    1               1          6     0     0     0
# 38 aci2025    1017 Norther~  14.4 ecoprovince            7    0.444           0.286      2     0     0     5
# 39 aci2025    1018 Norther~  14.4 ecoprovince            7    0.444           0.286      2     0     0     5
# 40 cgr         375 Lake of~   6.5 ecoprovince           11    0.167           0.0909     1     0     0    10
# 41 cgr         687 Boreal ~   9.2 ecoprovince           11   NA               0          0     0     0    11
# 42 cgr         688 Boreal ~   9.2 ecoprovince           11   NA               0          0     0     0    11
# 43 cgr         695 Boreal ~   9.2 ecoprovince           11   NA               0          0     0     0    11
# 44 cgr         701 Boreal ~   9.2 ecoprovince           11   NA               0          0     0     0    11
# 45 cgr         709 Boreal ~   9.2 ecoprovince           11   NA               0          0     0     0    11
# 46 pudl        375 Lake of~   6.5 ecoprovince           11    0.533           0.364      4     0     0     7
# 47 aci2025     618 Western~   9.1 ecoprovince           13    0.588           0.462      5     1     0     7
# 48 aci2025     631 Western~   9.1 ecoprovince           13    0.588           0.462      5     1     0     7
# 49 gpw         631 Western~   9.1 ecoprovince           17    0.692           0.529      9     0     0     8

# More info on rows with only 1 or 2 points (grouped)
accuracy_with_fallback_grouped %>%
  filter(N_evaluated %in% c(1, 2)) %>%
  select(
    raster,
    Ecds_ID,
    Ecrgn_n,
    Ecprv,
    source_level,
    N_evaluated,
    F1_grass,
    Overall_Accuracy
  ) %>%
  arrange(N_evaluated, raster, Ecds_ID) %>%
  print(n = Inf)

# without ACI2025:
#   raster Ecds_ID Ecrgn_n             Ecprv source_level N_evaluated F1_grass Overall_Accuracy
# 1 glc        374 Lac Seul Upland       6.2 ecoprovince            1   NA                  1  
# 2 ngi        672 Mid-Boreal Lowland    9.3 ecoprovince            1   NA                  0  
# 3 ngi        850 Lake Manitoba Plain  10.1 ecoprovince            1    1                  1  
# 4 cec        960 Cascade Ranges       13.2 ecoprovince            2   NA                  0  
# 5 cec        982 Fraser Basin         14.1 ecoprovince            2   NA                  0  
# 6 cnvc       960 Cascade Ranges       13.2 ecoprovince            2   NA                  0  
# 7 esa        960 Cascade Ranges       13.2 ecoprovince            2    0.667              0.5
# 8 esa        982 Fraser Basin         14.1 ecoprovince            2    1                  1  
# 9 glc        982 Fraser Basin         14.1 ecoprovince            2   NA                  0  
# 10 lcc        960 Cascade Ranges       13.2 ecoprovince            2   NA                  0  
# 11 lcc        982 Fraser Basin         14.1 ecoprovince            2    0.667              0.5
# 12 lcr        960 Cascade Ranges       13.2 ecoprovince            2   NA                  0  
# 13 lcr        982 Fraser Basin         14.1 ecoprovince            2   NA                  0  
# 14 luts10     960 Cascade Ranges       13.2 ecoprovince            2   NA                  0  
# 15 luts10     982 Fraser Basin         14.1 ecoprovince            2   NA                  0  
# 16 luts30     960 Cascade Ranges       13.2 ecoprovince            2   NA                  0  
# 17 luts30     982 Fraser Basin         14.1 ecoprovince            2   NA                  0 

# with ACI2025:
#    raster  Ecds_ID Ecrgn_n             Ecprv source_level N_evaluated F1_grass Overall_Accuracy
# 1 aci2024      374 Lac Seul Upland       6.2 ecoprovince            1   NA                  0  
# 2 aci2025      960 Cascade Ranges       13.2 ecoprovince            1   NA                  0  
# 3 glc          374 Lac Seul Upland       6.2 ecoprovince            1   NA                  1  
# 4 ngi          672 Mid-Boreal Lowland    9.3 ecoprovince            1   NA                  0  
# 5 ngi          850 Lake Manitoba Plain  10.1 ecoprovince            1    1                  1  
# 6 cec          960 Cascade Ranges       13.2 ecoprovince            2   NA                  0  
# 7 cec          982 Fraser Basin         14.1 ecoprovince            2   NA                  0  
# 8 cnvc         960 Cascade Ranges       13.2 ecoprovince            2   NA                  0  
# 9 esa          960 Cascade Ranges       13.2 ecoprovince            2    0.667              0.5
# 10 esa         982 Fraser Basin         14.1 ecoprovince            2    1                  1  
# 11 glc         982 Fraser Basin         14.1 ecoprovince            2   NA                  0  
# 12 lcc         960 Cascade Ranges       13.2 ecoprovince            2   NA                  0  
# 13 lcc         982 Fraser Basin         14.1 ecoprovince            2    0.667              0.5
# 14 lcr         960 Cascade Ranges       13.2 ecoprovince            2   NA                  0  
# 15 lcr         982 Fraser Basin         14.1 ecoprovince            2   NA                  0  
# 16 luts10      960 Cascade Ranges       13.2 ecoprovince            2   NA                  0  
# 17 luts10      982 Fraser Basin         14.1 ecoprovince            2   NA                  0  
# 18 luts30      960 Cascade Ranges       13.2 ecoprovince            2   NA                  0  
# 19 luts30      982 Fraser Basin         14.1 ecoprovince            2   NA                  0 

# exact match (not lenient):
#   raster  Ecds_ID Ecrgn_n             Ecprv source_level N_evaluated F1_grass Overall_Accuracy
# 1 aci2024     850 Lake Manitoba Plain  10.1 ecoprovince            1   NA                  0  
# 2 aci2025     960 Cascade Ranges       13.2 ecoprovince            1   NA                  0  
# 3 esri        850 Lake Manitoba Plain  10.1 ecoprovince            1   NA                  0  
# 4 gpw         850 Lake Manitoba Plain  10.1 ecoprovince            1    1                  1  
# 5 gpw         960 Cascade Ranges       13.2 ecoprovince            1   NA                  0  
# 6 ngi         672 Mid-Boreal Lowland    9.3 ecoprovince            1   NA                  0  
# 7 ngi         850 Lake Manitoba Plain  10.1 ecoprovince            1    1                  1  
# 8 aci2024     960 Cascade Ranges       13.2 ecoprovince            2   NA                  0  
# 9 aci2024     982 Fraser Basin         14.1 ecoprovince            2   NA                  0  
# 10 cec         960 Cascade Ranges       13.2 ecoprovince            2   NA                  0  
# 11 cec         982 Fraser Basin         14.1 ecoprovince            2   NA                  0  
# 12 cnvc        960 Cascade Ranges       13.2 ecoprovince            2   NA                  0  
# 13 esa         960 Cascade Ranges       13.2 ecoprovince            2    0.667              0.5
# 14 esa         982 Fraser Basin         14.1 ecoprovince            2    1                  1  
# 15 esri        960 Cascade Ranges       13.2 ecoprovince            2    1                  1  
# 16 esri        982 Fraser Basin         14.1 ecoprovince            2   NA                  0  
# 17 glc         960 Cascade Ranges       13.2 ecoprovince            2   NA                  0  
# 18 glc         982 Fraser Basin         14.1 ecoprovince            2   NA                  0  
# 19 gpw         982 Fraser Basin         14.1 ecoprovince            2    1                  1  
# 20 lcc         960 Cascade Ranges       13.2 ecoprovince            2   NA                  0  
# 21 lcc         982 Fraser Basin         14.1 ecoprovince            2    0.667              0.5
# 22 lcr         960 Cascade Ranges       13.2 ecoprovince            2   NA                  0  
# 23 lcr         982 Fraser Basin         14.1 ecoprovince            2   NA                  0  
# 24 luts10      960 Cascade Ranges       13.2 ecoprovince            2   NA                  0  
# 25 luts10      982 Fraser Basin         14.1 ecoprovince            2   NA                  0  
# 26 luts30      960 Cascade Ranges       13.2 ecoprovince            2   NA                  0  
# 27 luts30      982 Fraser Basin         14.1 ecoprovince            2   NA                  0  

accuracy_with_fallback_atomic %>%
  filter(N_evaluated %in% c(1, 2)) %>%
  select(
    raster,
    Ecds_ID,
    Ecrgn_n,
    Ecprv,
    source_level,
    N_evaluated,
    F1_grass,
    Overall_Accuracy
  ) %>%
  arrange(N_evaluated, raster, Ecds_ID) %>%
  print(n = Inf)

#   raster Ecds_ID Ecrgn_n             Ecprv source_level N_evaluated F1_grass Overall_Accuracy
# 1 glc        374 Lac Seul Upland       6.2 ecoprovince            1   NA                  1  
# 2 ngi        672 Mid-Boreal Lowland    9.3 ecoprovince            1   NA                  0  
# 3 ngi        850 Lake Manitoba Plain  10.1 ecoprovince            1    1                  1  
# 4 pli_mm     745 Aspen Parkland       10.2 ecoprovince            1   NA                  1  
# 5 cec        960 Cascade Ranges       13.2 ecoprovince            2   NA                  0  
# 6 cec        982 Fraser Basin         14.1 ecoprovince            2   NA                  0  
# 7 cnvc       960 Cascade Ranges       13.2 ecoprovince            2   NA                  0  
# 8 esa        960 Cascade Ranges       13.2 ecoprovince            2    0.667              0.5
# 9 esa        982 Fraser Basin         14.1 ecoprovince            2    1                  1  
# 10 glc        982 Fraser Basin         14.1 ecoprovince            2   NA                  0  
# 11 lcc        960 Cascade Ranges       13.2 ecoprovince            2   NA                  0  
# 12 lcc        982 Fraser Basin         14.1 ecoprovince            2    0.667              0.5
# 13 lcr        960 Cascade Ranges       13.2 ecoprovince            2   NA                  0  
# 14 lcr        982 Fraser Basin         14.1 ecoprovince            2   NA                  0  
# 15 luts10     960 Cascade Ranges       13.2 ecoprovince            2   NA                  0  
# 16 luts10     982 Fraser Basin         14.1 ecoprovince            2   NA                  0  
# 17 luts30     960 Cascade Ranges       13.2 ecoprovince            2   NA                  0  
# 18 luts30     982 Fraser Basin         14.1 ecoprovince            2   NA                  0

#------ How often is OA (as acc source) 0 or 1? ------#

accuracy_with_fallback_grouped_cons <- accuracy_with_fallback_grouped %>%
  mutate(
    acc = coalesce(F1_grass, Overall_Accuracy),
    acc_source = if_else(
      !is.na(F1_grass),
      "F1_grass",
      "Overall_Accuracy"
    )
  )

accuracy_with_fallback_grouped_cons %>%
  filter(
    acc_source == "Overall_Accuracy",
    acc %in% c(0, 1)
  ) %>%
  count(acc)

# without ACI2025:
#     acc     n
# 1     0    60
# 2     1    26

# with ACI2025:
#     acc     n
# 1     0    62
# 2     1    29

# exact match (not lenient):
#     acc     n
# 1     0    77
# 2     1     4

accuracy_with_fallback_atomic_cons <- accuracy_with_fallback_atomic %>%
  mutate(
    acc = coalesce(F1_grass, Overall_Accuracy),
    acc_source = if_else(
      !is.na(F1_grass),
      "F1_grass",
      "Overall_Accuracy"
    )
  )

accuracy_with_fallback_atomic_cons %>%
  filter(
    acc_source == "Overall_Accuracy",
    acc %in% c(0, 1)
  ) %>%
  count(acc)

#     acc     n
# 1     0    60
# 2     1    27

# Are these concentrated in raster/eco unit combos with small N?
accuracy_with_fallback_grouped_cons %>%
  filter(
    acc_source == "Overall_Accuracy",
    acc %in% c(0, 1)
  ) %>%
  summarise(
    n = n(),
    min_n = min(N_evaluated),
    median_n = median(N_evaluated),
    max_n = max(N_evaluated)
  )

# without ACI2025:
#  n min_n median_n max_n
# 86     1       40   329

# with ACI2025:
#  n min_n median_n max_n
# 91     1       40   329

# exact match (not lenient):
#   n min_n median_n max_n
#  81     1       27   247

accuracy_with_fallback_grouped_cons %>%
  filter(
    acc_source == "Overall_Accuracy",
    acc %in% c(0, 1)
  ) %>%
  count(N_evaluated, sort = TRUE) %>%
  print(n = Inf)

# without ACI2025:
#   N_evaluated     n
# 1          40    12
# 2           2    11
# 3         133    10
# 4          11     5
# 5          27     4
# 6          25     3
# 7          29     3
# 8           1     2
# 9          23     2
# 10         35     2
# 11         134     2
# 12         135     2
# 13          20     1
# 14          21     1
# 15          22     1
# 16          28     1
# 17          31     1
# 18          32     1
# 19          34     1
# 20          38     1
# 21          43     1
# 22          48     1
# 23          51     1
# 24          53     1
# 25          67     1
# 26          71     1
# 27          76     1
# 28          83     1
# 29          84     1
# 30          98     1
# 31         104     1
# 32         114     1
# 33         154     1
# 34         156     1
# 35         163     1
# 36         171     1
# 37         186     1
# 38         187     1
# 39         211     1
# 40         329     1

# with ACI2025:
#   N_evaluated     n
# 1          40    12
# 2           2    11
# 3         133    10
# 4          11     5
# 5           1     4
# 6          27     4
# 7          25     3
# 8          29     3
# 9           6     2
# 10          23     2
# 11          35     2
# 12         134     2
# 13         135     2
# 14          20     1
# 15          21     1
# 16          22     1
# 17          28     1
# 18          31     1
# 19          32     1
# 20          34     1
# 21          38     1
# 22          39     1
# 23          43     1
# 24          48     1
# 25          51     1
# 26          53     1
# 27          67     1
# 28          71     1
# 29          76     1
# 30          83     1
# 31          84     1
# 32          98     1
# 33         104     1
# 34         114     1
# 35         154     1
# 36         156     1
# 37         163     1
# 38         171     1
# 39         186     1
# 40         187     1
# 41         211     1
# 42         329     1

# exact match (not lenient):
#   N_evaluated     n
# 1           2    15
# 2         133     9
# 3          40     7
# 4           1     5
# 5          11     5
# 6          23     4
# 7          27     4
# 8          25     3
# 9           6     2
# 10          29     2
# 11          35     2
# 12         134     2
# 13         168     2
# 14          20     1
# 15          21     1
# 16          22     1
# 17          28     1
# 18          31     1
# 19          32     1
# 20          34     1
# 21          38     1
# 22          41     1
# 23          43     1
# 24          51     1
# 25          65     1
# 26          67     1
# 27          98     1
# 28         114     1
# 29         127     1
# 30         135     1
# 31         161     1
# 32         247     1

accuracy_with_fallback_grouped_cons %>%
  filter(
    acc_source == "Overall_Accuracy",
    acc %in% c(0, 1)
  ) %>%
  select(
    raster,
    Ecds_ID,
    Ecrgn_n,
    Ecprv,
    source_level,
    N_evaluated,
    acc,
    acc_source,
    TP, TN, FP, FN
  ) %>%
  arrange(
    desc(N_evaluated),
    raster,
    Ecds_ID
  ) %>%
  print(n = Inf)

# without ACI2025:
#   raster  Ecds_ID Ecrgn_n           Ecprv source_level N_evaluated   acc acc_source    TP    TN    FP    FN
# 1 glc         544 St. Lawrence Low~   8.1 ecodistrict          329     1 Overall_A~     0   329     0     0
# 2 glc         566 Lake Erie Lowland   8.2 ecodistrict          211     1 Overall_A~     0   211     0     0
# 3 glc         542 St. Lawrence Low~   8.1 ecodistrict          187     1 Overall_A~     0   187     0     0
# 4 glc         543 St. Lawrence Low~   8.1 ecodistrict          186     1 Overall_A~     0   186     0     0
# 5 aci2024     587 Peace Lowland       9.2 ecodistrict          171     1 Overall_A~     0   171     0     0
# 6 glc         411 Algonquin-Lake N~   6.6 ecodistrict          163     1 Overall_A~     0   163     0     0
# 7 glc         955 Eastern Vancouve~  13.1 ecoregion            156     1 Overall_A~     0   156     0     0
# 8 glc         956 Eastern Vancouve~  13.1 ecodistrict          154     1 Overall_A~     0   154     0     0
# 9 aci2024     607 Wabasca Lowland     9.2 ecodistrict          135     1 Overall_A~     0   135     0     0
# 10 dw          607 Wabasca Lowland     9.2 ecodistrict          135     1 Overall_A~     0   135     0     0
# 11 ngi        1017 Northern Contine~  14.4 ecoregion            134     0 Overall_A~     0     0     0   134
# 12 ngi        1019 Northern Contine~  14.4 ecoregion            134     0 Overall_A~     0     0     0   134
# 13 glc         410 Lac Temiscamingu~   6.6 ecoregion            133     1 Overall_A~     0   133     0     0
# 14 glc         587 Peace Lowland       9.2 ecoregion            133     0 Overall_A~     0     0     0   133
# 15 glc         588 Peace Lowland       9.2 ecoregion            133     0 Overall_A~     0     0     0   133
# 16 glc         590 Peace Lowland       9.2 ecoregion            133     0 Overall_A~     0     0     0   133
# 17 glc         592 Peace Lowland       9.2 ecoregion            133     0 Overall_A~     0     0     0   133
# 18 glc         594 Peace Lowland       9.2 ecoregion            133     0 Overall_A~     0     0     0   133
# 19 glc         596 Peace Lowland       9.2 ecoregion            133     0 Overall_A~     0     0     0   133
# 20 glc         597 Peace Lowland       9.2 ecoregion            133     0 Overall_A~     0     0     0   133
# 21 glc         598 Peace Lowland       9.2 ecoregion            133     0 Overall_A~     0     0     0   133
# 22 glc         600 Peace Lowland       9.2 ecoregion            133     0 Overall_A~     0     0     0   133
# 23 glc         979 Fraser Plateau     14.2 ecodistrict          114     0 Overall_A~     0     0     0   114
# 24 glc         408 Lac Temiscamingu~   6.6 ecodistrict          104     1 Overall_A~     0   104     0     0
# 25 ngi        1018 Northern Contine~  14.4 ecodistrict           98     0 Overall_A~     0     0     0    98
# 26 glc         546 St. Lawrence Low~   8.1 ecodistrict           84     1 Overall_A~     0    84     0     0
# 27 glc         545 St. Lawrence Low~   8.1 ecodistrict           83     1 Overall_A~     0    83     0     0
# 28 glc         383 Thunder Bay-Quet~   6.5 ecodistrict           76     1 Overall_A~     0    76     0     0
# 29 glc         412 Algonquin-Lake N~   6.6 ecodistrict           71     1 Overall_A~     0    71     0     0
# 30 glc         623 Western Alberta ~   9.1 ecodistrict           67     0 Overall_A~     0     0     0    67
# 31 glc         951 Pacific Ranges     13.2 ecoregion             53     1 Overall_A~     0    53     0     0
# 32 pudl        727 Aspen Parkland     10.2 ecodistrict           51     0 Overall_A~     0     0     0    51
# 33 glc         950 Pacific Ranges     13.2 ecodistrict           48     1 Overall_A~     0    48     0     0
# 34 pudl        684 Boreal Transition   9.2 ecodistrict           43     0 Overall_A~     0     0     0    43
# 35 aci2024     589 Peace Lowland       9.2 ecodistrict           40     1 Overall_A~     0    40     0     0
# 36 aci2024     601 Mid-Boreal Uplan~   9.2 ecodistrict           40     1 Overall_A~     0    40     0     0
# 37 cnvc        622 Western Alberta ~   9.1 ecoregion             40     0 Overall_A~     0     0     0    40
# 38 cnvc        623 Western Alberta ~   9.1 ecoregion             40     0 Overall_A~     0     0     0    40
# 39 cnvc        628 Western Alberta ~   9.1 ecoregion             40     0 Overall_A~     0     0     0    40
# 40 cnvc        629 Western Alberta ~   9.1 ecoregion             40     0 Overall_A~     0     0     0    40
# 41 cnvc        630 Western Alberta ~   9.1 ecoregion             40     0 Overall_A~     0     0     0    40
# 42 cnvc        631 Western Alberta ~   9.1 ecoregion             40     0 Overall_A~     0     0     0    40
# 43 dw          601 Mid-Boreal Uplan~   9.2 ecodistrict           40     1 Overall_A~     0    40     0     0
# 44 esri        589 Peace Lowland       9.2 ecodistrict           40     1 Overall_A~     0    40     0     0
# 45 esri        601 Mid-Boreal Uplan~   9.2 ecodistrict           40     1 Overall_A~     0    40     0     0
# 46 gpw         601 Mid-Boreal Uplan~   9.2 ecodistrict           40     0 Overall_A~     0     0    40     0
# 47 eccc       1017 Northern Contine~  14.4 ecodistrict           38     0 Overall_A~     0     0     0    38
# 48 glc         686 Boreal Transition   9.2 ecodistrict           35     0 Overall_A~     0     0     0    35
# 49 ngi         631 Western Alberta ~   9.1 ecodistrict           35     0 Overall_A~     0     0     0    35
# 50 glc         978 Fraser Plateau     14.2 ecodistrict           34     0 Overall_A~     0     0     0    34
# 51 glc         591 Peace Lowland       9.2 ecodistrict           32     0 Overall_A~     0     0     0    32
# 52 pudl        692 Boreal Transition   9.2 ecodistrict           31     0 Overall_A~     0     0     0    31
# 53 glc         572 Lake Erie Lowland   8.2 ecodistrict           29     1 Overall_A~     0    29     0     0
# 54 glc         692 Boreal Transition   9.2 ecodistrict           29     0 Overall_A~     0     0     0    29
# 55 ngi        1016 Northern Contine~  14.4 ecodistrict           29     0 Overall_A~     0     0     0    29
# 56 ngi         686 Boreal Transition   9.2 ecodistrict           28     0 Overall_A~     0     0     0    28
# 57 cnvc        678 Boreal Transition   9.2 ecodistrict           27     0 Overall_A~     0     0     0    27
# 58 glc         586 Peace Lowland       9.2 ecodistrict           27     0 Overall_A~     0     0     0    27
# 59 ngi         703 Boreal Transition   9.2 ecodistrict           27     0 Overall_A~     0     0     0    27
# 60 pudl        755 Aspen Parkland     10.2 ecodistrict           27     0 Overall_A~     0     0     0    27
# 61 lcr         586 Peace Lowland       9.2 ecodistrict           25     0 Overall_A~     0     0     0    25
# 62 luts10      586 Peace Lowland       9.2 ecodistrict           25     0 Overall_A~     0     0     0    25
# 63 luts30      586 Peace Lowland       9.2 ecodistrict           25     0 Overall_A~     0     0     0    25
# 64 glc         409 Lac Temiscamingu~   6.6 ecodistrict           23     1 Overall_A~     0    23     0     0
# 65 glc         593 Peace Lowland       9.2 ecodistrict           23     0 Overall_A~     0     0     0    23
# 66 cnvc        684 Boreal Transition   9.2 ecodistrict           22     0 Overall_A~     0     0     0    22
# 67 pudl        590 Peace Lowland       9.2 ecodistrict           21     0 Overall_A~     0     0     0    21
# 68 pudl        839 Aspen Parkland     10.2 ecodistrict           20     0 Overall_A~     0     0     0    20
# 69 cgr         687 Boreal Transition   9.2 ecoprovince           11     0 Overall_A~     0     0     0    11
# 70 cgr         688 Boreal Transition   9.2 ecoprovince           11     0 Overall_A~     0     0     0    11
# 71 cgr         695 Boreal Transition   9.2 ecoprovince           11     0 Overall_A~     0     0     0    11
# 72 cgr         701 Boreal Transition   9.2 ecoprovince           11     0 Overall_A~     0     0     0    11
# 73 cgr         709 Boreal Transition   9.2 ecoprovince           11     0 Overall_A~     0     0     0    11
# 74 cec         960 Cascade Ranges     13.2 ecoprovince            2     0 Overall_A~     0     0     0     2
# 75 cec         982 Fraser Basin       14.1 ecoprovince            2     0 Overall_A~     0     0     0     2
# 76 cnvc        960 Cascade Ranges     13.2 ecoprovince            2     0 Overall_A~     0     0     0     2
# 77 glc         982 Fraser Basin       14.1 ecoprovince            2     0 Overall_A~     0     0     0     2
# 78 lcc         960 Cascade Ranges     13.2 ecoprovince            2     0 Overall_A~     0     0     0     2
# 79 lcr         960 Cascade Ranges     13.2 ecoprovince            2     0 Overall_A~     0     0     0     2
# 80 lcr         982 Fraser Basin       14.1 ecoprovince            2     0 Overall_A~     0     0     0     2
# 81 luts10      960 Cascade Ranges     13.2 ecoprovince            2     0 Overall_A~     0     0     0     2
# 82 luts10      982 Fraser Basin       14.1 ecoprovince            2     0 Overall_A~     0     0     0     2
# 83 luts30      960 Cascade Ranges     13.2 ecoprovince            2     0 Overall_A~     0     0     0     2
# 84 luts30      982 Fraser Basin       14.1 ecoprovince            2     0 Overall_A~     0     0     0     2
# 85 glc         374 Lac Seul Upland     6.2 ecoprovince            1     1 Overall_A~     0     1     0     0
# 86 ngi         672 Mid-Boreal Lowla~   9.3 ecoprovince            1     0 Overall_A~     0     0     0     1

#with ACI2025:
#   raster  Ecds_ID Ecrgn_n           Ecprv source_level N_evaluated   acc acc_source    TP    TN    FP    FN
# 1 glc         544 St. Lawrence Low~   8.1 ecodistrict          329     1 Overall_A~     0   329     0     0
# 2 glc         566 Lake Erie Lowland   8.2 ecodistrict          211     1 Overall_A~     0   211     0     0
# 3 glc         542 St. Lawrence Low~   8.1 ecodistrict          187     1 Overall_A~     0   187     0     0
# 4 glc         543 St. Lawrence Low~   8.1 ecodistrict          186     1 Overall_A~     0   186     0     0
# 5 aci2024     587 Peace Lowland       9.2 ecodistrict          171     1 Overall_A~     0   171     0     0
# 6 glc         411 Algonquin-Lake N~   6.6 ecodistrict          163     1 Overall_A~     0   163     0     0
# 7 glc         955 Eastern Vancouve~  13.1 ecoregion            156     1 Overall_A~     0   156     0     0
# 8 glc         956 Eastern Vancouve~  13.1 ecodistrict          154     1 Overall_A~     0   154     0     0
# 9 aci2024     607 Wabasca Lowland     9.2 ecodistrict          135     1 Overall_A~     0   135     0     0
# 10 dw          607 Wabasca Lowland     9.2 ecodistrict          135     1 Overall_A~     0   135     0     0
# 11 ngi        1017 Northern Contine~  14.4 ecoregion            134     0 Overall_A~     0     0     0   134
# 12 ngi        1019 Northern Contine~  14.4 ecoregion            134     0 Overall_A~     0     0     0   134
# 13 glc         410 Lac Temiscamingu~   6.6 ecoregion            133     1 Overall_A~     0   133     0     0
# 14 glc         587 Peace Lowland       9.2 ecoregion            133     0 Overall_A~     0     0     0   133
# 15 glc         588 Peace Lowland       9.2 ecoregion            133     0 Overall_A~     0     0     0   133
# 16 glc         590 Peace Lowland       9.2 ecoregion            133     0 Overall_A~     0     0     0   133
# 17 glc         592 Peace Lowland       9.2 ecoregion            133     0 Overall_A~     0     0     0   133
# 18 glc         594 Peace Lowland       9.2 ecoregion            133     0 Overall_A~     0     0     0   133
# 19 glc         596 Peace Lowland       9.2 ecoregion            133     0 Overall_A~     0     0     0   133
# 20 glc         597 Peace Lowland       9.2 ecoregion            133     0 Overall_A~     0     0     0   133
# 21 glc         598 Peace Lowland       9.2 ecoregion            133     0 Overall_A~     0     0     0   133
# 22 glc         600 Peace Lowland       9.2 ecoregion            133     0 Overall_A~     0     0     0   133
# 23 glc         979 Fraser Plateau     14.2 ecodistrict          114     0 Overall_A~     0     0     0   114
# 24 glc         408 Lac Temiscamingu~   6.6 ecodistrict          104     1 Overall_A~     0   104     0     0
# 25 ngi        1018 Northern Contine~  14.4 ecodistrict           98     0 Overall_A~     0     0     0    98
# 26 glc         546 St. Lawrence Low~   8.1 ecodistrict           84     1 Overall_A~     0    84     0     0
# 27 glc         545 St. Lawrence Low~   8.1 ecodistrict           83     1 Overall_A~     0    83     0     0
# 28 glc         383 Thunder Bay-Quet~   6.5 ecodistrict           76     1 Overall_A~     0    76     0     0
# 29 glc         412 Algonquin-Lake N~   6.6 ecodistrict           71     1 Overall_A~     0    71     0     0
# 30 glc         623 Western Alberta ~   9.1 ecodistrict           67     0 Overall_A~     0     0     0    67
# 31 glc         951 Pacific Ranges     13.2 ecoregion             53     1 Overall_A~     0    53     0     0
# 32 pudl        727 Aspen Parkland     10.2 ecodistrict           51     0 Overall_A~     0     0     0    51
# 33 glc         950 Pacific Ranges     13.2 ecodistrict           48     1 Overall_A~     0    48     0     0
# 34 pudl        684 Boreal Transition   9.2 ecodistrict           43     0 Overall_A~     0     0     0    43
# 35 aci2024     589 Peace Lowland       9.2 ecodistrict           40     1 Overall_A~     0    40     0     0
# 36 aci2024     601 Mid-Boreal Uplan~   9.2 ecodistrict           40     1 Overall_A~     0    40     0     0
# 37 cnvc        622 Western Alberta ~   9.1 ecoregion             40     0 Overall_A~     0     0     0    40
# 38 cnvc        623 Western Alberta ~   9.1 ecoregion             40     0 Overall_A~     0     0     0    40
# 39 cnvc        628 Western Alberta ~   9.1 ecoregion             40     0 Overall_A~     0     0     0    40
# 40 cnvc        629 Western Alberta ~   9.1 ecoregion             40     0 Overall_A~     0     0     0    40
# 41 cnvc        630 Western Alberta ~   9.1 ecoregion             40     0 Overall_A~     0     0     0    40
# 42 cnvc        631 Western Alberta ~   9.1 ecoregion             40     0 Overall_A~     0     0     0    40
# 43 dw          601 Mid-Boreal Uplan~   9.2 ecodistrict           40     1 Overall_A~     0    40     0     0
# 44 esri        589 Peace Lowland       9.2 ecodistrict           40     1 Overall_A~     0    40     0     0
# 45 esri        601 Mid-Boreal Uplan~   9.2 ecodistrict           40     1 Overall_A~     0    40     0     0
# 46 gpw         601 Mid-Boreal Uplan~   9.2 ecodistrict           40     0 Overall_A~     0     0    40     0
# 47 aci2024     571 Lake Erie Lowland   8.2 ecodistrict           39     1 Overall_A~     0    39     0     0
# 48 eccc       1017 Northern Contine~  14.4 ecodistrict           38     0 Overall_A~     0     0     0    38
# 49 glc         686 Boreal Transition   9.2 ecodistrict           35     0 Overall_A~     0     0     0    35
# 50 ngi         631 Western Alberta ~   9.1 ecodistrict           35     0 Overall_A~     0     0     0    35
# 51 glc         978 Fraser Plateau     14.2 ecodistrict           34     0 Overall_A~     0     0     0    34
# 52 glc         591 Peace Lowland       9.2 ecodistrict           32     0 Overall_A~     0     0     0    32
# 53 pudl        692 Boreal Transition   9.2 ecodistrict           31     0 Overall_A~     0     0     0    31
# 54 glc         572 Lake Erie Lowland   8.2 ecodistrict           29     1 Overall_A~     0    29     0     0
# 55 glc         692 Boreal Transition   9.2 ecodistrict           29     0 Overall_A~     0     0     0    29
# 56 ngi        1016 Northern Contine~  14.4 ecodistrict           29     0 Overall_A~     0     0     0    29
# 57 ngi         686 Boreal Transition   9.2 ecodistrict           28     0 Overall_A~     0     0     0    28
# 58 cnvc        678 Boreal Transition   9.2 ecodistrict           27     0 Overall_A~     0     0     0    27
# 59 glc         586 Peace Lowland       9.2 ecodistrict           27     0 Overall_A~     0     0     0    27
# 60 ngi         703 Boreal Transition   9.2 ecodistrict           27     0 Overall_A~     0     0     0    27
# 61 pudl        755 Aspen Parkland     10.2 ecodistrict           27     0 Overall_A~     0     0     0    27
# 62 lcr         586 Peace Lowland       9.2 ecodistrict           25     0 Overall_A~     0     0     0    25
# 63 luts10      586 Peace Lowland       9.2 ecodistrict           25     0 Overall_A~     0     0     0    25
# 64 luts30      586 Peace Lowland       9.2 ecodistrict           25     0 Overall_A~     0     0     0    25
# 65 glc         409 Lac Temiscamingu~   6.6 ecodistrict           23     1 Overall_A~     0    23     0     0
# 66 glc         593 Peace Lowland       9.2 ecodistrict           23     0 Overall_A~     0     0     0    23
# 67 cnvc        684 Boreal Transition   9.2 ecodistrict           22     0 Overall_A~     0     0     0    22
# 68 pudl        590 Peace Lowland       9.2 ecodistrict           21     0 Overall_A~     0     0     0    21
# 69 pudl        839 Aspen Parkland     10.2 ecodistrict           20     0 Overall_A~     0     0     0    20
# 70 cgr         687 Boreal Transition   9.2 ecoprovince           11     0 Overall_A~     0     0     0    11
# 71 cgr         688 Boreal Transition   9.2 ecoprovince           11     0 Overall_A~     0     0     0    11
# 72 cgr         695 Boreal Transition   9.2 ecoprovince           11     0 Overall_A~     0     0     0    11
# 73 cgr         701 Boreal Transition   9.2 ecoprovince           11     0 Overall_A~     0     0     0    11
# 74 cgr         709 Boreal Transition   9.2 ecoprovince           11     0 Overall_A~     0     0     0    11
# 75 aci2025     955 Eastern Vancouve~  13.1 ecoprovince            6     1 Overall_A~     0     6     0     0
# 76 aci2025     956 Eastern Vancouve~  13.1 ecoprovince            6     1 Overall_A~     0     6     0     0
# 77 cec         960 Cascade Ranges     13.2 ecoprovince            2     0 Overall_A~     0     0     0     2
# 78 cec         982 Fraser Basin       14.1 ecoprovince            2     0 Overall_A~     0     0     0     2
# 79 cnvc        960 Cascade Ranges     13.2 ecoprovince            2     0 Overall_A~     0     0     0     2
# 80 glc         982 Fraser Basin       14.1 ecoprovince            2     0 Overall_A~     0     0     0     2
# 81 lcc         960 Cascade Ranges     13.2 ecoprovince            2     0 Overall_A~     0     0     0     2
# 82 lcr         960 Cascade Ranges     13.2 ecoprovince            2     0 Overall_A~     0     0     0     2
# 83 lcr         982 Fraser Basin       14.1 ecoprovince            2     0 Overall_A~     0     0     0     2
# 84 luts10      960 Cascade Ranges     13.2 ecoprovince            2     0 Overall_A~     0     0     0     2
# 85 luts10      982 Fraser Basin       14.1 ecoprovince            2     0 Overall_A~     0     0     0     2
# 86 luts30      960 Cascade Ranges     13.2 ecoprovince            2     0 Overall_A~     0     0     0     2
# 87 luts30      982 Fraser Basin       14.1 ecoprovince            2     0 Overall_A~     0     0     0     2
# 88 aci2024     374 Lac Seul Upland     6.2 ecoprovince            1     0 Overall_A~     0     0     1     0
# 89 aci2025     960 Cascade Ranges     13.2 ecoprovince            1     0 Overall_A~     0     0     0     1
# 90 glc         374 Lac Seul Upland     6.2 ecoprovince            1     1 Overall_A~     0     1     0     0
# 91 ngi         672 Mid-Boreal Lowla~   9.3 ecoprovince            1     0 Overall_A~     0     0     0     1

# exact match (not lenient):
#   raster  Ecds_ID Ecrgn_n           Ecprv source_level N_evaluated   acc acc_source    TP    TN    FP    FN
# 1 glc        1006 Thompson-Okanaga~  14.3 ecodistrict          247     0 Overall_A~     0     0     0   247
# 2 glc        1002 Interior Transit~  14.3 ecoregion            168     0 Overall_A~     0     0     0   168
# 3 glc        1003 Interior Transit~  14.3 ecoregion            168     0 Overall_A~     0     0     0   168
# 4 glc        1001 Interior Transit~  14.3 ecodistrict          161     0 Overall_A~     0     0     0   161
# 5 dw          607 Wabasca Lowland     9.2 ecodistrict          135     1 Overall_A~     0   135     0     0
# 6 ngi        1017 Northern Contine~  14.4 ecoregion            134     0 Overall_A~     0     0     0   134
# 7 ngi        1019 Northern Contine~  14.4 ecoregion            134     0 Overall_A~     0     0     0   134
# 8 glc         587 Peace Lowland       9.2 ecoregion            133     0 Overall_A~     0     0     0   133
# 9 glc         588 Peace Lowland       9.2 ecoregion            133     0 Overall_A~     0     0     0   133
# 10 glc         590 Peace Lowland       9.2 ecoregion            133     0 Overall_A~     0     0     0   133
# 11 glc         592 Peace Lowland       9.2 ecoregion            133     0 Overall_A~     0     0     0   133
# 12 glc         594 Peace Lowland       9.2 ecoregion            133     0 Overall_A~     0     0     0   133
# 13 glc         596 Peace Lowland       9.2 ecoregion            133     0 Overall_A~     0     0     0   133
# 14 glc         597 Peace Lowland       9.2 ecoregion            133     0 Overall_A~     0     0     0   133
# 15 glc         598 Peace Lowland       9.2 ecoregion            133     0 Overall_A~     0     0     0   133
# 16 glc         600 Peace Lowland       9.2 ecoregion            133     0 Overall_A~     0     0     0   133
# 17 esri        567 Lake Erie Lowland   8.2 ecodistrict          127     0 Overall_A~     0     0     0   127
# 18 glc         979 Fraser Plateau     14.2 ecodistrict          114     0 Overall_A~     0     0     0   114
# 19 ngi        1018 Northern Contine~  14.4 ecodistrict           98     0 Overall_A~     0     0     0    98
# 20 glc         623 Western Alberta ~   9.1 ecodistrict           67     0 Overall_A~     0     0     0    67
# 21 esri        568 Lake Erie Lowland   8.2 ecodistrict           65     0 Overall_A~     0     0     0    65
# 22 pudl        727 Aspen Parkland     10.2 ecodistrict           51     0 Overall_A~     0     0     0    51
# 23 pudl        684 Boreal Transition   9.2 ecodistrict           43     0 Overall_A~     0     0     0    43
# 24 esri        570 Lake Erie Lowland   8.2 ecodistrict           41     0 Overall_A~     0     0     0    41
# 25 cnvc        622 Western Alberta ~   9.1 ecoregion             40     0 Overall_A~     0     0     0    40
# 26 cnvc        623 Western Alberta ~   9.1 ecoregion             40     0 Overall_A~     0     0     0    40
# 27 cnvc        628 Western Alberta ~   9.1 ecoregion             40     0 Overall_A~     0     0     0    40
# 28 cnvc        629 Western Alberta ~   9.1 ecoregion             40     0 Overall_A~     0     0     0    40
# 29 cnvc        630 Western Alberta ~   9.1 ecoregion             40     0 Overall_A~     0     0     0    40
# 30 cnvc        631 Western Alberta ~   9.1 ecoregion             40     0 Overall_A~     0     0     0    40
# 31 dw          601 Mid-Boreal Uplan~   9.2 ecodistrict           40     1 Overall_A~     0    40     0     0
# 32 eccc       1017 Northern Contine~  14.4 ecodistrict           38     0 Overall_A~     0     0     0    38
# 33 glc         686 Boreal Transition   9.2 ecodistrict           35     0 Overall_A~     0     0     0    35
# 34 ngi         631 Western Alberta ~   9.1 ecodistrict           35     0 Overall_A~     0     0     0    35
# 35 glc         978 Fraser Plateau     14.2 ecodistrict           34     0 Overall_A~     0     0     0    34
# 36 glc         591 Peace Lowland       9.2 ecodistrict           32     0 Overall_A~     0     0     0    32
# 37 pudl        692 Boreal Transition   9.2 ecodistrict           31     0 Overall_A~     0     0     0    31
# 38 glc         692 Boreal Transition   9.2 ecodistrict           29     0 Overall_A~     0     0     0    29
# 39 ngi        1016 Northern Contine~  14.4 ecodistrict           29     0 Overall_A~     0     0     0    29
# 40 ngi         686 Boreal Transition   9.2 ecodistrict           28     0 Overall_A~     0     0     0    28
# 41 cnvc        678 Boreal Transition   9.2 ecodistrict           27     0 Overall_A~     0     0     0    27
# 42 glc         586 Peace Lowland       9.2 ecodistrict           27     0 Overall_A~     0     0     0    27
# 43 ngi         703 Boreal Transition   9.2 ecodistrict           27     0 Overall_A~     0     0     0    27
# 44 pudl        755 Aspen Parkland     10.2 ecodistrict           27     0 Overall_A~     0     0     0    27
# 45 lcr         586 Peace Lowland       9.2 ecodistrict           25     0 Overall_A~     0     0     0    25
# 46 luts10      586 Peace Lowland       9.2 ecodistrict           25     0 Overall_A~     0     0     0    25
# 47 luts30      586 Peace Lowland       9.2 ecodistrict           25     0 Overall_A~     0     0     0    25
# 48 glc         593 Peace Lowland       9.2 ecodistrict           23     0 Overall_A~     0     0     0    23
# 49 glc         988 Columbia Mountai~  14.4 ecoregion             23     0 Overall_A~     0     0     0    23
# 50 glc         990 Columbia Mountai~  14.4 ecoregion             23     0 Overall_A~     0     0     0    23
# 51 glc         991 Columbia Mountai~  14.4 ecoregion             23     0 Overall_A~     0     0     0    23
# 52 cnvc        684 Boreal Transition   9.2 ecodistrict           22     0 Overall_A~     0     0     0    22
# 53 pudl        590 Peace Lowland       9.2 ecodistrict           21     0 Overall_A~     0     0     0    21
# 54 pudl        839 Aspen Parkland     10.2 ecodistrict           20     0 Overall_A~     0     0     0    20
# 55 cgr         687 Boreal Transition   9.2 ecoprovince           11     0 Overall_A~     0     0     0    11
# 56 cgr         688 Boreal Transition   9.2 ecoprovince           11     0 Overall_A~     0     0     0    11
# 57 cgr         695 Boreal Transition   9.2 ecoprovince           11     0 Overall_A~     0     0     0    11
# 58 cgr         701 Boreal Transition   9.2 ecoprovince           11     0 Overall_A~     0     0     0    11
# 59 cgr         709 Boreal Transition   9.2 ecoprovince           11     0 Overall_A~     0     0     0    11
# 60 aci2025     955 Eastern Vancouve~  13.1 ecoprovince            6     1 Overall_A~     0     6     0     0
# 61 aci2025     956 Eastern Vancouve~  13.1 ecoprovince            6     1 Overall_A~     0     6     0     0
# 62 aci2024     960 Cascade Ranges     13.2 ecoprovince            2     0 Overall_A~     0     0     0     2
# 63 aci2024     982 Fraser Basin       14.1 ecoprovince            2     0 Overall_A~     0     0     0     2
# 64 cec         960 Cascade Ranges     13.2 ecoprovince            2     0 Overall_A~     0     0     0     2
# 65 cec         982 Fraser Basin       14.1 ecoprovince            2     0 Overall_A~     0     0     0     2
# 66 cnvc        960 Cascade Ranges     13.2 ecoprovince            2     0 Overall_A~     0     0     0     2
# 67 esri        982 Fraser Basin       14.1 ecoprovince            2     0 Overall_A~     0     0     0     2
# 68 glc         960 Cascade Ranges     13.2 ecoprovince            2     0 Overall_A~     0     0     0     2
# 69 glc         982 Fraser Basin       14.1 ecoprovince            2     0 Overall_A~     0     0     0     2
# 70 lcc         960 Cascade Ranges     13.2 ecoprovince            2     0 Overall_A~     0     0     0     2
# 71 lcr         960 Cascade Ranges     13.2 ecoprovince            2     0 Overall_A~     0     0     0     2
# 72 lcr         982 Fraser Basin       14.1 ecoprovince            2     0 Overall_A~     0     0     0     2
# 73 luts10      960 Cascade Ranges     13.2 ecoprovince            2     0 Overall_A~     0     0     0     2
# 74 luts10      982 Fraser Basin       14.1 ecoprovince            2     0 Overall_A~     0     0     0     2
# 75 luts30      960 Cascade Ranges     13.2 ecoprovince            2     0 Overall_A~     0     0     0     2
# 76 luts30      982 Fraser Basin       14.1 ecoprovince            2     0 Overall_A~     0     0     0     2
# 77 aci2024     850 Lake Manitoba Pl~  10.1 ecoprovince            1     0 Overall_A~     0     0     0     1
# 78 aci2025     960 Cascade Ranges     13.2 ecoprovince            1     0 Overall_A~     0     0     0     1
# 79 esri        850 Lake Manitoba Pl~  10.1 ecoprovince            1     0 Overall_A~     0     0     0     1
# 80 gpw         960 Cascade Ranges     13.2 ecoprovince            1     0 Overall_A~     0     0     0     1
# 81 ngi         672 Mid-Boreal Lowla~   9.3 ecoprovince            1     0 Overall_A~     0     0     0     1

# ----------------------------
# 16) Export file for consensus mapping
# ----------------------------

# Primary metric: F1_grass 
# Fallback metric: Overall Accuracy 

# Grouped
accuracy_with_fallback_grouped_cons <- accuracy_with_fallback_grouped %>%
  mutate(
    acc = coalesce(F1_grass, Overall_Accuracy),
    acc_source = if_else(
      !is.na(F1_grass),
      "F1_grass",
      "Overall_Accuracy"
    )
  )

# Add quality flag for rows with N less than 20
accuracy_with_fallback_grouped_cons <- accuracy_with_fallback_grouped_cons %>%
  mutate(
    low_confidence_acc =
      #acc_source == "Overall_Accuracy" &
      N_evaluated < 20
  )

# Atomic
accuracy_with_fallback_atomic_cons <- accuracy_with_fallback_atomic %>%
  mutate(
    acc = coalesce(F1_grass, Overall_Accuracy),
    acc_source = if_else(
      !is.na(F1_grass),
      "F1_grass",
      "Overall_Accuracy"
    )
  )

# Add quality flag for rows with N less than 20
accuracy_with_fallback_atomic_cons <- accuracy_with_fallback_atomic_cons %>%
  mutate(
    low_confidence_acc =
      #acc_source == "Overall_Accuracy" &
      N_evaluated < 20
  )

#---- Export ----

# GROUPED
readr::write_csv(
  accuracy_with_fallback_grouped_cons,
  file.path(out_dir, "binary_accuracy_grouped_combinedCIandCM_forconsensus_28Jul2026.csv")
)

# ATOMIC
readr::write_csv(
  accuracy_with_fallback_atomic_cons,
  file.path(out_dir, "binary_accuracy_atomic_combinedCIandCM_forconsensus_28Jul2026.csv")
)

#----- Plot ----

# Grouped
accuracy_with_fallback_grouped_cons %>%
  group_by(raster) %>%
  summarise(acc = mean(acc, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = reorder(raster, acc), y = acc)) +
  geom_col() +
  coord_flip() +
  scale_y_continuous(labels = scales::label_percent()) +
  labs(
    x = "Dataset",
    y = "Accuracy (%)",
    title = "Accuracy by dataset (grass vs non-grass)"
  ) +
  theme_bw()

# with number of points
# Accuracy values
acc_df <- accuracy_with_fallback_grouped_cons %>%
  group_by(raster) %>%
  summarise(acc = mean(acc, na.rm = TRUE), .groups = "drop")

# Validation point counts
count_df <- combined %>%
  count(raster_group, name = "n") %>%
  rename(raster = raster_group)

# Join counts to accuracy table
plot_df <- acc_df %>%
  left_join(count_df, by = "raster")

p <- ggplot(plot_df,
       aes(x = reorder(raster, acc), y = acc)) +
  geom_col() +
  geom_text(
    aes(y = acc * 0.75, 
        label = paste0(scales::percent(acc, accuracy = 0.1),
                       " (n=", scales::comma(n), ")")),
    color = "white",
    size = 4
  ) +
  coord_flip() +
  scale_y_continuous(labels = scales::label_percent()) +
  labs(
    x = "Dataset",
    y = "Accuracy (%)",
    title = "Accuracy by dataset (grass vs non-grass)"
  ) +
  theme_bw()

p

ggsave(file.path(out_dir, "accuracy_by_dataset_binary.png"),
  plot = p,
  width = 8,
  height = 6,
  dpi = 300
)
