
####################################################################################
##            PRELIMINARY ACCURACY ASSESSMENT: BINARY GRASS VS NON-GRASS          ##
####################################################################################

##----------------------- Alison Long - June 2026 -------------------------

# step 8: assess accuracy of datasets mapping grasslands by comparing them
#         to ground-truthed point data for 2015 to present

## THIS VERSION COMBINES COMBINES COMPILED AND CROP-INSURANCE DERIVED POINT DATA
## Ran separate scripts for the two types of point data - will merge here and
## calculate metrics

## NEED TO REMOVE MORE TRAINING DATA FOR ACI??

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
# 2) File paths 
# ----------------------------
acc_cm_file <- "R:/CanadianGrasslands/WG_Classification/results/accuracy_assessment/binary/compiled_derived/accuracy_assessment_data_binary_pli_sep_compiled_23Jun2026.csv"
acc_ci_file <- "R:/CanadianGrasslands/WG_Classification/results/accuracy_assessment/binary/crop_ins_derived/accuracy_assessment_data_binary_pli_sep_cropinsV2_lenient_23Jun2026.csv"
metadata_file <- "R:/CanadianGrasslands/WG_Classification/github/tomorrows-prairies/processeddatasets/Grass/metadata_full_extent_binary_pli_sep.csv"

out_dir <- "R:/CanadianGrasslands/WG_Classification/results/accuracy_assessment/binary/combined_CI_and_compiled"

# ----------------------------
# 3) Configuration
# ----------------------------
# IMPORTANT: Handle PLI as four separate rasters in the points file and metadata.
pli_variants <- c("pli_ap", "pli_cu", "pli_mg", "pli_mm")

# Grouping helper: collapse pli_* -> "pli"; all others remain themselves
raster_to_group <- function(x) ifelse(x %in% pli_variants, "pli", x)

# ----------------------------
# 4) Load points and normalize order of fields
# ----------------------------
# Load
acc_cm <- readr::read_csv(acc_cm_file, show_col_types = FALSE)
acc_ci <- readr::read_csv(acc_ci_file, show_col_types = FALSE)

# Inspect structure
str(acc_cm)
str(acc_ci)

# Add missing columns
#acc_ci <- acc_ci %>%
#  mutate(Training = "")

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
readr::write_csv(combined_out, file.path(out_dir, "accuracy_assessment_data_binary_pli_sep_CIandCM_23Jun2026.csv"))

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
accuracy_by_ecoregion <- eval_data %>%
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
accuracy_by_ecoprovince <- eval_data %>%
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

accuracy_by_ecoprovince <- accuracy_by_ecoprovince %>%
  rename_with(~ paste0(.x, "_ecoP"), -c(raster, Ecprv))

#---- Define mapping ----

# Get mapping (one row per ecodistrict)
ecodist_to_regions <- eval_data %>%
  distinct(Ecds_ID, Ecrgn_n, Ecprv)

accuracy_by_ecodistrict_atomic2 <- accuracy_by_ecodistrict_atomic %>%
  left_join(ecodist_to_regions, by = "Ecds_ID")

accuracy_by_ecodistrict_grouped2 <- accuracy_by_ecodistrict_grouped %>%
  left_join(ecodist_to_regions, by = "Ecds_ID")

#------- Join ----------
fallback_joined_atomic <- accuracy_by_ecodistrict_atomic2 %>%
  left_join(
    accuracy_by_ecoregion,
    by = c("raster", "Ecrgn_n"),
    suffix = c("_ecoD", "_ecoR")
  ) %>%
  left_join(
    accuracy_by_ecoprovince,
    by = c("raster", "Ecprv")
    #suffix = c("", "_ecoP")
  )

# check names
names(fallback_joined_atomic)

fallback_joined_grouped <- accuracy_by_ecodistrict_grouped2 %>%
  left_join(
    accuracy_by_ecoregion,
    by = c("raster", "Ecrgn_n"),
    suffix = c("_ecoD", "_ecoR")
  ) %>%
  left_join(
    accuracy_by_ecoprovince,
    by = c("raster", "Ecprv")
    #suffix = c("", "_ecoP")
  )

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
  file.path(out_dir, "binary_accuracy_ecoregprovfallback_grouped_combinedCIandCM_date.csv")
)

# ATOMIC
readr::write_csv(
  accuracy_with_fallback_atomic,
  file.path(out_dir, "binary_accuracy_ecoregprovfallback_atomic_combinedCIandCM_date.csv")
)

