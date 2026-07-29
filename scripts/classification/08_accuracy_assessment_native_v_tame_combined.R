
####################################################################################
##                 PRELIMINARY ACCURACY ASSESSMENT: NATIVE VS TAME                ##
####################################################################################

##----------------------- Alison Long - June 2026 -------------------------

# step 8 (multiclass): assess accuracy of datasets mapping grasslands
# against ground-truthed point data (2015 to present), with 3 classes:
#   0 = not grass, 1 = tame grass, 2 = native grass

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
##   5) Produces both "grouped" (pli combined) and "atomic" (each raster) outputs
##
## Time Rules (per record):
##   - must have non-NA raster prediction (0 or 1)
##   - class-based time rule using allowed years for that raster:
##       NATIVE: any( allowed_year <= Yer_Cll ); matched_year = max(allowed_year ≤ Yer_Cll)
##       TAME:   any( allowed_year == Yer_Cll ); matched_year = Yer_Cll (if allowed)
##       OTHER:  any( allowed_year == Yer_Cll ); matched_year = Yer_Cll (if allowed)
##
## Exclusions:
##   - Source used for training: exclude if Source_norm matches either atomic raster code or grouped alias
##
# Key changes from binary version:
#   - Uses LC_class3 (0/1/2) instead of LC_grass (0/1)
#   - Rasters expected as 0/1/2
#   - Time-aligned rules:
#         native (2): latest allowed_year <= Yer_Cll
#         tame   (1): allowed_year == Yer_Cll
#         other  (0): allowed_year == Yer_Cll
#   - 3×3 confusion matrices and per-class metrics
#   - More accuracy metrics (not just recall)

# ----------------------------
# 1) Setup
# ----------------------------

# load packages
library(readxl)
library(writexl)
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
acc_cm_dir <- "R:/CanadianGrasslands/WG_Classification/results/accuracy_assessment/native_v_tame/compiled_derived"
acc_ci_dir <- "R:/CanadianGrasslands/WG_Classification/results/accuracy_assessment/native_v_tame/crop_ins_derived"
out_dir    <- "R:/CanadianGrasslands/WG_Classification/results/accuracy_assessment/native_v_tame/combined_CI_and_compiled"

acc_cm_file <- file.path(acc_cm_dir, "18Jul2026/eval_data_18Jul2026_natvtam_withnongrass_compiled.csv")
acc_ci_file <- file.path(acc_ci_dir, "18Jul2026/eval_data_18Jul2026_natvtam_withnongrass_lenient.csv")
metadata_file <- "R:/CanadianGrasslands/WG_Classification/github/tomorrows-prairies/processeddatasets/Grass/metadata_full_extent_binary_pli_sep.csv"

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

# Count ACI_field_2024 point/raster combos in acc_cm (10,716)
acc_cm %>%
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

# Make sure aci raster has same name for both
# acc_cm <- acc_cm %>%
#   mutate(raster = if_else(raster == "aci", "aci2024", raster))
# 
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
readr::write_csv(combined_out, file.path(out_dir, "accuracy_assessment_data_natvtam_pli_sep_CIandCM_CIcombined_18Jul2026.csv"))

# ----------------------------
# 8) Multiclass metrics helpers
# ----------------------------

safe_div <- function(n, d) ifelse(d == 0 | is.na(d), NA_real_, n / d)

# ----------------------------
# MCC helpers
# ----------------------------
binary_mcc <- function(tp, tn, fp, fn) {
  
  # force numeric (double) to avoid integer overflow
  tp <- as.numeric(tp)
  tn <- as.numeric(tn)
  fp <- as.numeric(fp)
  fn <- as.numeric(fn)
  
  num <- tp * tn - fp * fn
  
  den <- sqrt(
    (tp + fp) *
      (tp + fn) *
      (tn + fp) *
      (tn + fn)
  )
  
  # robust guard
  if (is.na(den) || den == 0) return(NA_real_)
  
  num / den
}

multiclass_mcc <- function(cm) {
  
  cm <- as.matrix(cm)
  if (!is.numeric(cm)) stop("Confusion matrix must be numeric")
  
  n <- as.numeric(sum(cm))
  if (n == 0) return(NA_real_)
  
  r <- as.numeric(rowSums(cm))
  c <- as.numeric(colSums(cm))
  
  sum_diag <- sum(diag(cm))
  sum_row_col <- sum(r * c)
  
  numerator <- sum_diag * n - sum_row_col
  
  denom_row <- n^2 - sum(r^2)
  denom_col <- n^2 - sum(c^2)
  
  if (denom_row <= 0 || denom_col <= 0) return(NA_real_)
  
  numerator / sqrt(denom_row * denom_col)
}

# ----------------------------
# Build 3x3 confusion matrix
# ----------------------------
compute_confusion_3 <- function(df,
                                truth_col = "LC_class3",
                                pred_col  = "pred") {
  
  lvls <- c(0L, 1L, 2L)
  
  truth <- factor(df[[truth_col]], levels = lvls)
  pred  <- factor(df[[pred_col]],  levels = lvls)
  
  tab <- table(truth, pred)
  tab <- as.matrix(tab)[as.character(lvls), as.character(lvls), drop = FALSE]
  
  counts_long <- as_tibble(as.data.frame(as.table(tab))) %>%
    transmute(
      Truth = as.integer(as.character(truth)),
      Pred  = as.integer(as.character(pred)),
      N     = as.integer(Freq)
    )
  
  list(tab = tab, n = sum(tab), counts_long = counts_long)
}

# ----------------------------
# Compute metrics from 3x3 CM
# ----------------------------
compute_metrics_3 <- function(tab) {
  
  n     <- as.numeric(sum(tab))
  diagv <- as.numeric(diag(tab))
  rows  <- as.numeric(rowSums(tab))   # support (truth)
  cols  <- as.numeric(colSums(tab))   # predicted totals
  
  # Number of classes actually present
  n_truth_classes <- sum(rows > 0)
  n_pred_classes <- sum(cols > 0)
  
  # ---- Accuracy ----
  overall <- safe_div(sum(diagv), n)
  
  # ---- Precision / Recall / F1 ----
  prec <- safe_div(diagv, cols)
  rec  <- safe_div(diagv, rows)
  f1   <- ifelse(
    is.na(prec) | is.na(rec) | (prec + rec) == 0,
    NA_real_,
    2 * prec * rec / (prec + rec)
  )
  
  macro_prec <- mean(prec, na.rm = TRUE)
  macro_rec  <- mean(rec,  na.rm = TRUE)
  macro_f1   <- mean(f1,   na.rm = TRUE)
  
  balanced_acc <- macro_rec
  
  weights     <- safe_div(rows, sum(rows))
  weighted_f1 <- sum(weights * f1, na.rm = TRUE)
  
  # ---- Cohen's kappa ----
  pe     <- safe_div(sum(rows * cols), n * n)
  kappa <- ifelse(is.na(pe) | (1 - pe) == 0, NA_real_, (overall - pe) / (1 - pe))
  
  # ---- Multiclass MCC ----
  mcc_multiclass <- multiclass_mcc(tab)
  
  # ---- Degenerate truth check ----
  n_truth_classes <- sum(rows > 0)
  
  if (n_truth_classes < 2) {
    kappa <- NA_real_           # make sure NA not 0
    mcc_multiclass <- NA_real_
  }
  
  # ---- One-vs-rest MCC ----
  classes <- seq_len(nrow(tab))
  
  ovr_mcc <- sapply(classes, function(i) {
    tp <- tab[i, i]
    fn <- sum(tab[i, -i])
    fp <- sum(tab[-i, i])
    tn <- sum(tab[-i, -i])
    binary_mcc(tp, tn, fp, fn)
  })
  
  mcc_ovr_macro <- mean(ovr_mcc, na.rm = TRUE)
  
  tibble(
    N_evaluated        = as.integer(n),
    
    Truth_Classes_Present = n_truth_classes,
    Pred_Classes_Present  = n_pred_classes,
    
    Overall_Accuracy  = overall,
    Balanced_Accuracy = balanced_acc,
    Kappa             = kappa,
    
    MCC_Multiclass    = mcc_multiclass,
    MCC_OvR_Macro     = mcc_ovr_macro,
    MCC_OvR_0         = ovr_mcc[1],
    MCC_OvR_1         = ovr_mcc[2],
    MCC_OvR_2         = ovr_mcc[3],
    
    Macro_Precision   = macro_prec,
    Macro_Recall      = macro_rec,
    Macro_F1          = macro_f1,
    Weighted_F1       = weighted_f1,
    
    Precision_0 = prec[1], Recall_0 = rec[1], F1_0 = f1[1], Support_0 = rows[1],
    Precision_1 = prec[2], Recall_1 = rec[2], F1_1 = f1[2], Support_1 = rows[2],
    Precision_2 = prec[3], Recall_2 = rec[3], F1_2 = f1[3], Support_2 = rows[3]
  )
}

# ----------------------------
# Summarize metrics per group
# ----------------------------
summarize_block_multiclass <- function(df, key_cols) {
  
  df %>%
    group_by(across(all_of(key_cols))) %>%
    group_modify(~ {
      c3  <- compute_confusion_3(.x)
      met <- compute_metrics_3(c3$tab)
      met
    }) %>%
    ungroup()
}

# ----------------------------
# Return long-form confusion
# ----------------------------
confusion_block_multiclass <- function(df, key_cols) {
  
  df %>%
    group_by(across(all_of(key_cols))) %>%
    group_modify(~ compute_confusion_3(.x)$counts_long) %>%
    ungroup()
}

# ----------------------------
# Create output tables
# ----------------------------
metrics_by_dataset <- summarize_block_multiclass(
  combined,
  key_cols = "raster_group"
)

confusions_by_dataset <- confusion_block_multiclass(
  combined,
  key_cols = "raster_group"
)

metrics_by_dataset

# ----------------------------
# Export output tables
# ----------------------------
out_xlsx <- file.path(out_dir, "nat_v_tame_v_nongrass_accuracy_assessment_results_combined_18Jul2026.xlsx")

write_xlsx(
  list(
    "Metrics_by_Dataset"   = metrics_by_dataset,
    "Confusion_Long_Form"  = confusions_by_dataset
  ),
  path = out_xlsx
)

#####################################################################################
#####################################################################################

# ----------------------------
# 13) Accuracy by Ecodistrict (TRICLASS)
# ----------------------------

if (!"Ecds_ID" %in% names(combined)) {
  stop("Ecds_ID column not found in combined.")
}

# ---- ATOMIC: raster × ecodistrict ----
accuracy_by_ecodistrict_atomic <- summarize_block_multiclass(
  combined,
  key_cols = c("raster", "Ecds_ID")
) %>%
  arrange(raster, Ecds_ID)

# ---- GROUPED: collapse pli_* -> pli ----
accuracy_by_ecodistrict_grouped <- combined %>%
  mutate(raster = raster_group) %>%
  summarize_block_multiclass(
    key_cols = c("raster", "Ecds_ID")
  ) %>%
  arrange(raster, Ecds_ID)

# ---- OPTIONAL: atomic with year ----
accuracy_by_ecodistrict_year_atomic <- summarize_block_multiclass(
  combined,
  key_cols = c("raster", "matched_year", "Ecds_ID")
) %>%
  rename(raster_year = matched_year) %>%
  arrange(raster, raster_year, Ecds_ID)

# ---- Save ----
readr::write_csv(accuracy_by_ecodistrict_grouped,
                 file.path(out_dir, "accuracy_by_ecodistrict_grouped_combined.csv"))

readr::write_csv(accuracy_by_ecodistrict_atomic,
                 file.path(out_dir, "accuracy_by_ecodistrict_atomic_combined.csv"))

readr::write_csv(accuracy_by_ecodistrict_year_atomic,
                 file.path(out_dir, "accuracy_by_ecodistrict_year_atomic_combined.csv"))

# ----------------------------
# 14) Multiclass fallback hierarchy
# ecodistrict -> ecoregion -> ecoprovince
# ----------------------------

min_n <- 20

# ---- Ecoregion metrics ----
accuracy_by_ecoregion_atomic <- summarize_block_multiclass(
  combined,
  key_cols = c("raster", "Ecrgn_n")
)

accuracy_by_ecoregion_grouped <- combined %>%
mutate(raster = raster_group) %>%
summarize_block_multiclass(
  key_cols = c("raster", "Ecrgn_n")
)

# ---- Ecoprovince metrics ----
accuracy_by_ecoprovince_atomic <- summarize_block_multiclass(
  combined,
  key_cols = c("raster", "Ecprv")
) %>%
  rename_with(
    ~ paste0(.x, "_ecoP"),
    -c(raster, Ecprv)
  )

accuracy_by_ecoprovince_grouped <- combined %>%
  mutate(raster = raster_group) %>%
  summarize_block_multiclass(
    key_cols = c("raster", "Ecprv")
  ) %>%
  rename_with(
    ~ paste0(.x, "_ecoP"),
    -c(raster, Ecprv)
  )

# ---- Mapping table ----
ecodist_to_regions <- combined %>%
  distinct(Ecds_ID, Ecrgn_n, Ecprv)

# --------------------------------------------------
# ATOMIC
# --------------------------------------------------

accuracy_by_ecodistrict_atomic2 <- accuracy_by_ecodistrict_atomic %>%
  left_join(ecodist_to_regions, by = "Ecds_ID")

fallback_joined_atomic <- accuracy_by_ecodistrict_atomic2 %>%
  left_join(
    accuracy_by_ecoregion_atomic,
    by = c("raster", "Ecrgn_n"),
    suffix = c("_ecoD", "_ecoR")
  ) %>%
  left_join(
    accuracy_by_ecoprovince_atomic,
    by = c("raster", "Ecprv")
  )

# Make sure suffixes are present 
names(fallback_joined_atomic)

metric_cols <- setdiff(
  names(accuracy_by_ecodistrict_atomic),
  c("raster", "Ecds_ID")
)

accuracy_with_fallback_atomic <- fallback_joined_atomic %>%
  mutate(
    source_level = case_when(
      N_evaluated_ecoD >= min_n ~ "ecodistrict",
      N_evaluated_ecoR >= min_n ~ "ecoregion",
      TRUE ~ "ecoprovince"
    )
  )

for (m in metric_cols) {
  
  accuracy_with_fallback_atomic[[m]] <- case_when(
    accuracy_with_fallback_atomic$source_level == "ecodistrict" ~
      accuracy_with_fallback_atomic[[paste0(m, "_ecoD")]],
    
    accuracy_with_fallback_atomic$source_level == "ecoregion" ~
      accuracy_with_fallback_atomic[[paste0(m, "_ecoR")]],
    
    TRUE ~
      accuracy_with_fallback_atomic[[paste0(m, "_ecoP")]]
  )
}

# --------------------------------------------------
# GROUPED
# --------------------------------------------------

accuracy_by_ecodistrict_grouped2 <- accuracy_by_ecodistrict_grouped %>%
  left_join(ecodist_to_regions, by = "Ecds_ID")

fallback_joined_grouped <- accuracy_by_ecodistrict_grouped2 %>%
  left_join(
    accuracy_by_ecoregion_grouped,
    by = c("raster", "Ecrgn_n"),
    suffix = c("_ecoD", "_ecoR")
  ) %>%
  left_join(
    accuracy_by_ecoprovince_grouped,
    by = c("raster", "Ecprv")
  )

# Make sure suffixes are present 
names(fallback_joined_grouped)

accuracy_with_fallback_grouped <- fallback_joined_grouped %>%
  mutate(
    source_level = case_when(
      N_evaluated_ecoD >= min_n ~ "ecodistrict",
      N_evaluated_ecoR >= min_n ~ "ecoregion",
      TRUE ~ "ecoprovince"
    )
  )

for (m in metric_cols) {
  
  accuracy_with_fallback_grouped[[m]] <- case_when(
    accuracy_with_fallback_grouped$source_level == "ecodistrict" ~
      accuracy_with_fallback_grouped[[paste0(m, "_ecoD")]],
    
    accuracy_with_fallback_grouped$source_level == "ecoregion" ~
      accuracy_with_fallback_grouped[[paste0(m, "_ecoR")]],
    
    TRUE ~
      accuracy_with_fallback_grouped[[paste0(m, "_ecoP")]]
  )
}

# ----------------------------
# Export multiclass fallback results
# ----------------------------

readr::write_csv(
  accuracy_with_fallback_atomic,
  file.path(
    out_dir,
    "nattam_accuracy_by_ecodistrict_atomic_ecoregprovfallback_combined_18Jul2026.csv"
  )
)

readr::write_csv(
  accuracy_with_fallback_grouped,
  file.path(
    out_dir,
    "nattam_accuracy_by_ecodistrict_grouped_ecoregprovfallback_combined_18Jul2026.csv"
  )
)

# ----------------------------
# 15) Diagnostics
# ----------------------------
table(combined$raster)

table(combined$raster_group)

sum(is.na(combined$raster_group))

combined %>%
  count(raster, raster_group) %>%
  filter(raster %in% c("pli_ap","pli_cu","pli_mg","pli_mm"))

fallback_joined_grouped %>%
  summarise(
    ecoR_missing =
      sum(is.na(N_evaluated_ecoR)),
    ecoP_missing =
      sum(is.na(N_evaluated_ecoP))
  )

fallback_joined_grouped %>%
  filter(
    is.na(N_evaluated_ecoR) |
      is.na(N_evaluated_ecoP)
  ) %>%
  count(raster)

#------ NA rates for major metrics at each level ------#

fallback_joined_grouped %>%
  summarise(
    n_rows = n(),
    
    OA_ecoD_NA = sum(is.na(Overall_Accuracy_ecoD)),
    OA_ecoR_NA = sum(is.na(Overall_Accuracy_ecoR)),
    OA_ecoP_NA = sum(is.na(Overall_Accuracy_ecoP)),
    
    BA_ecoD_NA = sum(is.na(Balanced_Accuracy_ecoD)),
    BA_ecoR_NA = sum(is.na(Balanced_Accuracy_ecoR)),
    BA_ecoP_NA = sum(is.na(Balanced_Accuracy_ecoP)),
    
    Kappa_ecoD_NA = sum(is.na(Kappa_ecoD)),
    Kappa_ecoR_NA = sum(is.na(Kappa_ecoR)),
    Kappa_ecoP_NA = sum(is.na(Kappa_ecoP)),
    
    MCC_ecoD_NA = sum(is.na(MCC_Multiclass_ecoD)),
    MCC_ecoR_NA = sum(is.na(MCC_Multiclass_ecoR)),
    MCC_ecoP_NA = sum(is.na(MCC_Multiclass_ecoP))
  ) %>%
  print(width = Inf)

# n_rows OA_ecoD_NA OA_ecoR_NA OA_ecoP_NA BA_ecoD_NA BA_ecoR_NA BA_ecoP_NA Kappa_ecoD_NA Kappa_ecoR_NA
#   1605          0          0          0          0          0          0           524           102
# Kappa_ecoP_NA MCC_ecoD_NA MCC_ecoR_NA MCC_ecoP_NA
#            43         553         117          49

fallback_joined_atomic %>%
  summarise(
    n_rows = n(),
    
    OA_ecoD_NA = sum(is.na(Overall_Accuracy_ecoD)),
    OA_ecoR_NA = sum(is.na(Overall_Accuracy_ecoR)),
    OA_ecoP_NA = sum(is.na(Overall_Accuracy_ecoP)),
    
    BA_ecoD_NA = sum(is.na(Balanced_Accuracy_ecoD)),
    BA_ecoR_NA = sum(is.na(Balanced_Accuracy_ecoR)),
    BA_ecoP_NA = sum(is.na(Balanced_Accuracy_ecoP)),
    
    Kappa_ecoD_NA = sum(is.na(Kappa_ecoD)),
    Kappa_ecoR_NA = sum(is.na(Kappa_ecoR)),
    Kappa_ecoP_NA = sum(is.na(Kappa_ecoP)),
    
    MCC_ecoD_NA = sum(is.na(MCC_Multiclass_ecoD)),
    MCC_ecoR_NA = sum(is.na(MCC_Multiclass_ecoR)),
    MCC_ecoP_NA = sum(is.na(MCC_Multiclass_ecoP))
  ) %>%
  print(width = Inf)

# n_rows OA_ecoD_NA OA_ecoR_NA OA_ecoP_NA BA_ecoD_NA BA_ecoR_NA BA_ecoP_NA Kappa_ecoD_NA Kappa_ecoR_NA
#   1607          0          0          0          0          0          0           527           124
# Kappa_ecoP_NA MCC_ecoD_NA MCC_ecoR_NA MCC_ecoP_NA
#            66         556         139          72

#------ Final fallback NA rates ------#

accuracy_with_fallback_grouped %>%
summarise(
  n_rows = n(),
  
  OA_NA = sum(is.na(Overall_Accuracy)),
  OA_pct_NA = 100 * mean(is.na(Overall_Accuracy)),
  
  Macro_F1_NA = sum(is.na(Macro_F1)),
  Macro_F1_pct_NA = 100 * mean(is.na(Macro_F1)),
  
  Weighted_F1_NA = sum(is.na(Weighted_F1)),
  Weighted_F1_pct_NA = 100 * mean(is.na(Weighted_F1)),
  
  BA_NA = sum(is.na(Balanced_Accuracy)),
  BA_pct_NA = 100 * mean(is.na(Balanced_Accuracy)),
  
  Kappa_NA = sum(is.na(Kappa)),
  Kappa_pct_NA = 100 * mean(is.na(Kappa)),
  
  MCC_NA = sum(is.na(MCC_Multiclass)),
  MCC_pct_NA = 100 * mean(is.na(MCC_Multiclass))
) %>%
  print(width = Inf)

# n_rows OA_NA OA_pct_NA Macro_F1_NA Macro_F1_pct_NA Weighted_F1_NA Weighted_F1_pct_NA BA_NA BA_pct_NA
#   1605     0         0          51            3.18              0                  0     0         0
# Kappa_NA Kappa_pct_NA MCC_NA MCC_pct_NA
#      346         21.6    361       22.5

accuracy_with_fallback_atomic %>%
  summarise(
    n_rows = n(),
    
    OA_NA = sum(is.na(Overall_Accuracy)),
    OA_pct_NA = 100 * mean(is.na(Overall_Accuracy)),
    
    Macro_F1_NA = sum(is.na(Macro_F1)),
    Macro_F1_pct_NA = 100 * mean(is.na(Macro_F1)),
    
    Weighted_F1_NA = sum(is.na(Weighted_F1)),
    Weighted_F1_pct_NA = 100 * mean(is.na(Weighted_F1)),
    
    BA_NA = sum(is.na(Balanced_Accuracy)),
    BA_pct_NA = 100 * mean(is.na(Balanced_Accuracy)),
    
    Kappa_NA = sum(is.na(Kappa)),
    Kappa_pct_NA = 100 * mean(is.na(Kappa)),
    
    MCC_NA = sum(is.na(MCC_Multiclass)),
    MCC_pct_NA = 100 * mean(is.na(MCC_Multiclass))
  ) %>%
  print(width = Inf)

# n_rows OA_NA OA_pct_NA Macro_F1_NA Macro_F1_pct_NA Weighted_F1_NA Weighted_F1_pct_NA BA_NA BA_pct_NA
#   1607     0         0          51            3.17              0                  0     0         0
# Kappa_NA Kappa_pct_NA MCC_NA MCC_pct_NA
#      359         22.3    374       23.3

#------ Where are MCC/Kappa missing? ------#

accuracy_with_fallback_grouped %>%
  count(Truth_Classes_Present)

#  Truth_Classes_Present     n
#1                     1   346
#2                     2   728
#3                     3   531

accuracy_with_fallback_grouped %>%
  count(Pred_Classes_Present)

#   Pred_Classes_Present     n
# 1                    1    33
# 2                    2   315
# 3                    3  1257

accuracy_with_fallback_grouped %>%
  summarise(
    one_truth_class =
      sum(Truth_Classes_Present == 1),
    
    two_truth_classes =
      sum(Truth_Classes_Present == 2),
    
    three_truth_classes =
      sum(Truth_Classes_Present == 3),
    
    one_pred_class =
      sum(Pred_Classes_Present == 1),
    
    two_pred_classes =
      sum(Pred_Classes_Present == 2),
    
    three_pred_classes =
      sum(Pred_Classes_Present == 3),
    
    Kappa_NA =
      sum(is.na(Kappa)),
    
    MCC_NA =
      sum(is.na(MCC_Multiclass))
  ) %>%
  print(width = Inf)

# one_truth_class two_truth_classes three_truth_classes one_pred_class two_pred_classes three_pred_classes
#             346               728                 531             33              315               1257
# Kappa_NA MCC_NA
#      346    361

# Which class is driving the 346 single-class units?
accuracy_with_fallback_grouped %>%
  summarise(
    only_class0 =
      sum(
        Support_0 > 0 &
          Support_1 == 0 &
          Support_2 == 0
      ),
    
    only_class1 =
      sum(
        Support_0 == 0 &
          Support_1 > 0 &
          Support_2 == 0
      ),
    
    only_class2 =
      sum(
        Support_0 == 0 &
          Support_1 == 0 &
          Support_2 > 0
      )
  )

# only_class0 only_class1 only_class2
#          24          11         311

# Which class pairings dominate the 728 two-class units?
accuracy_with_fallback_grouped %>%
  mutate(
    combo = case_when(
      Support_0 > 0 & Support_1 > 0 & Support_2 == 0 ~ "0+1",
      Support_0 > 0 & Support_1 == 0 & Support_2 > 0 ~ "0+2",
      Support_0 == 0 & Support_1 > 0 & Support_2 > 0 ~ "1+2",
      TRUE ~ NA_character_
    )
  ) %>%
  count(combo)

#   combo     n
# 1 0+1     441
# 2 0+2      33
# 3 1+2     254
# 4 NA      877

#------ Which classes are absent? ------#

accuracy_with_fallback_grouped %>%
summarise(
  no_class0 = sum(Support_0 == 0, na.rm = TRUE),
  no_class1 = sum(Support_1 == 0, na.rm = TRUE),
  no_class2 = sum(Support_2 == 0, na.rm = TRUE)
)

# no_class0 no_class1 no_class2
#       576       368       476

accuracy_with_fallback_atomic %>%
  summarise(
    no_class0 = sum(Support_0 == 0, na.rm = TRUE),
    no_class1 = sum(Support_1 == 0, na.rm = TRUE),
    no_class2 = sum(Support_2 == 0, na.rm = TRUE)
  )

# no_class0 no_class1 no_class2
#       588       369       477

# ----- Missing one-vs-rest MCCs ------#

accuracy_with_fallback_grouped %>%
  summarise(
    MCC0_NA = sum(is.na(MCC_OvR_0)),
    MCC1_NA = sum(is.na(MCC_OvR_1)),
    MCC2_NA = sum(is.na(MCC_OvR_2))
  )

# MCC0_NA MCC1_NA MCC2_NA
#     603     440     868 - highest number of NAs for native

accuracy_with_fallback_grouped %>%
  summarise(
    no_class0_truth = sum(Support_0 == 0),
    no_class1_truth = sum(Support_1 == 0),
    no_class2_truth = sum(Support_2 == 0),
    
    no_class0_pred = sum(Precision_0 %>% is.na()),
    no_class1_pred = sum(Precision_1 %>% is.na()),
    no_class2_pred = sum(Precision_2 %>% is.na())
  )

# Native MCC problems are because the class is absent from truth and pred:
# no_class0_truth no_class1_truth no_class2_truth no_class0_pred no_class1_pred no_class2_pred
#             576             368             476             10            106            265

accuracy_with_fallback_atomic %>%
  summarise(
    MCC0_NA = sum(is.na(MCC_OvR_0)),
    MCC1_NA = sum(is.na(MCC_OvR_1)),
    MCC2_NA = sum(is.na(MCC_OvR_2))
  )

# MCC0_NA MCC1_NA MCC2_NA
#     616     441     881

#------ Grouped vs atomic results ------#

accuracy_with_fallback_grouped %>%
  summarise(
    OA_NA = sum(is.na(Overall_Accuracy)),
    MCC_NA = sum(is.na(MCC_Multiclass))
  )

# OA_NA MCC_NA
#     0    361

accuracy_with_fallback_atomic %>%
  summarise(
    OA_NA = sum(is.na(Overall_Accuracy)),
    MCC_NA = sum(is.na(MCC_Multiclass))
  )

# OA_NA MCC_NA
#     0    374

#------ Check for grouped/atomic key mismatches ------#

atomic_keys <- accuracy_by_ecodistrict_atomic %>%
  distinct(raster, Ecds_ID)

group_keys <- accuracy_by_ecodistrict_grouped %>%
  distinct(raster, Ecds_ID)

atomic_keys_grouped <- accuracy_by_ecodistrict_atomic %>%
  mutate(
    raster = if_else(
      raster %in% c(
        "pli_ap",
        "pli_cu",
        "pli_mg",
        "pli_mm"
      ),
      "pli",
      raster
    )
  ) %>%
  distinct(raster, Ecds_ID)

anti_join(
  group_keys,
  atomic_keys_grouped,
  by = c("raster", "Ecds_ID")
)

anti_join(
  atomic_keys_grouped,
  group_keys,
  by = c("raster", "Ecds_ID")
)

#------ Which classes are often missing ------#

accuracy_with_fallback_grouped %>%
  count(
    Truth_Classes_Present,
    Pred_Classes_Present
  ) %>%
  arrange(
    Truth_Classes_Present,
    Pred_Classes_Present
  )

#   Truth_Classes_Present Pred_Classes_Present     n
# 1                     1                    1    18
# 2                     1                    2    60
# 3                     1                    3   268
# 4                     2                    1    15
# 5                     2                    2   204
# 6                     2                    3   509
# 7                     3                    2    51
# 8                     3                    3   480

ggplot(accuracy_with_fallback_grouped,
  aes(
    factor(Truth_Classes_Present),
    factor(Pred_Classes_Present))) +
  geom_count()

accuracy_with_fallback_atomic %>%
  count(
    Truth_Classes_Present,
    Pred_Classes_Present
  ) %>%
  arrange(
    Truth_Classes_Present,
    Pred_Classes_Present
  )

#   Truth_Classes_Present Pred_Classes_Present     n
# 1                     1                    1    19
# 2                     1                    2    60
# 3                     1                    3   280
# 4                     2                    1    15
# 5                     2                    2   204
# 6                     2                    3   497
# 7                     3                    2    51
# 8                     3                    3   481

ggplot(accuracy_with_fallback_atomic,
       aes(
         factor(Truth_Classes_Present),
         factor(Pred_Classes_Present))) +
  geom_count()

#------ Check correlations among metrics ------#

accuracy_with_fallback_grouped %>%
  select(
    Overall_Accuracy,
    Balanced_Accuracy,
    Weighted_F1,
    Macro_F1
  ) %>%
  cor(use = "pairwise.complete.obs")

#                   Overall_Accuracy Balanced_Accuracy Weighted_F1  Macro_F1
# Overall_Accuracy         1.0000000         0.8568095   0.6588147 0.8817312
# Balanced_Accuracy        0.8568095         1.0000000   0.4776229 0.7680445
# Weighted_F1              0.6588147         0.4776229   1.0000000 0.4552734
# Macro_F1                 0.8817312         0.7680445   0.4552734 1.0000000

#------ Check distribution of sample sizes --------#

accuracy_with_fallback_grouped %>%
  ggplot(aes(x = N_evaluated)) +
  geom_histogram(
    binwidth = 5,
    color = "black",
    fill = "steelblue"
  ) +
  geom_vline(
    xintercept = 20,
    linetype = "dashed",
    color = "red",
    linewidth = 1
  ) +
  labs(
    title = "Distribution of Validation Sample Sizes",
    x = "N evaluated",
    y = "Count of ecological units"
  ) +
  theme_bw()

# log-scaled version to address right skew
accuracy_with_fallback_grouped %>%
  ggplot(aes(x = N_evaluated)) +
  geom_histogram(
    bins = 50,
    color = "black",
    fill = "steelblue"
  ) +
  scale_x_log10() +
  labs(
    title = "Distribution of Validation Sample Sizes",
    x = "N evaluated (log10 scale)",
    y = "Count of ecological units"
  ) +
  theme_bw()

#------ Number of rows with n less than 20 ------#

accuracy_with_fallback_grouped %>%
  summarise(
    n_below_20 = sum(N_evaluated < 20),
    pct_below_20 = 100 * mean(N_evaluated < 20)
  )

# n_below_20 pct_below_20
#         10        0.623

# Examine low-sample rows
accuracy_with_fallback_grouped %>%
  filter(N_evaluated < 20) %>%
  select(
    raster,
    Ecds_ID,
    Ecrgn_n,
    Ecprv,
    source_level,
    N_evaluated,
    Truth_Classes_Present,
    Pred_Classes_Present,
    Macro_F1,
    Balanced_Accuracy,
    Weighted_F1,
    Overall_Accuracy,
    Support_0,
    Support_1,
    Support_2
  ) %>%
  arrange(N_evaluated, raster, Ecds_ID) %>%
  print(width = Inf)

#   raster Ecds_ID Ecrgn_n                Ecprv source_level N_evaluated Truth_Classes_Present
# 1 ngi        623 Western Alberta Upland   9.1 ecoprovince            1                     1
# 2 ngi        850 Lake Manitoba Plain     10.1 ecoprovince            1                     1
# 3 cnvc       960 Cascade Ranges          13.2 ecoprovince            2                     1
# 4 cgr        375 Lake of the Woods        6.5 ecoprovince           11                     2
# 5 cgr        687 Boreal Transition        9.2 ecoprovince           11                     1
# 6 cgr        688 Boreal Transition        9.2 ecoprovince           11                     1
# 7 cgr        695 Boreal Transition        9.2 ecoprovince           11                     1
# 8 cgr        701 Boreal Transition        9.2 ecoprovince           11                     1
# 9 cgr        709 Boreal Transition        9.2 ecoprovince           11                     1
# 10 pudl      375 Lake of the Woods        6.5 ecoprovince           11                     1

#   Pred_Classes_Present Macro_F1 Balanced_Accuracy Weighted_F1 Overall_Accuracy Support_0 Support_1
# 1                    1  NaN                 0               0            0             0         0
# 2                    1  NaN                 0               0            0             0         0
# 3                    1  NaN                 0               0            0             0         0
# 4                    2  NaN                 0               0            0             0         4
# 5                    1  NaN                 0               0            0             0         0
# 6                    1  NaN                 0               0            0             0         0
# 7                    1  NaN                 0               0            0             0         0
# 8                    1  NaN                 0               0            0             0         0
# 9                    1  NaN                 0               0            0             0         0
# 10                   3    0.429             0.273           0            0.273         0         0

# Support_2
# 1         1
# 2         1
# 3         2
# 4         7
# 5        11
# 6        11
# 7        11
# 8        11
# 9        11
# 10       11

#------ Examine rows with only 1 or 2 points ------#

accuracy_with_fallback_grouped %>%
  filter(N_evaluated %in% c(1, 2)) %>%
  select(
    raster,
    Ecds_ID,
    Ecrgn_n,
    Ecprv,
    source_level,
    N_evaluated,
    Truth_Classes_Present,
    Pred_Classes_Present,
    Macro_F1,
    Balanced_Accuracy,
    Weighted_F1,
    Overall_Accuracy
  ) %>%
  arrange(N_evaluated, raster, Ecds_ID) %>%
  print(width = Inf)

#   raster Ecds_ID Ecrgn_n                Ecprv source_level N_evaluated Truth_Classes_Present
# 1 ngi        623 Western Alberta Upland   9.1 ecoprovince            1                     1
# 2 ngi        850 Lake Manitoba Plain     10.1 ecoprovince            1                     1
# 3 cnvc       960 Cascade Ranges          13.2 ecoprovince            2                     1

#   Pred_Classes_Present Macro_F1 Balanced_Accuracy Weighted_F1 Overall_Accuracy
# 1                    1      NaN                 0           0                0
# 2                    1      NaN                 0           0                0
# 3                    1      NaN                 0           0                0

#------ How often is fallback used? ------#

accuracy_with_fallback_grouped_cons %>%
  count(acc_source)

#   acc_source            n
# 1 Balanced_Accuracy    51
# 2 Macro_F1           1554

accuracy_with_fallback_grouped_cons %>%
  filter(acc_source == "Balanced_Accuracy") %>%
  summarise(
    n = n(),
    pct = 100 * n() / nrow(accuracy_with_fallback_grouped_cons)
  )

#  n   pct
# 51  3.18

#------ Why is Macro_F1 undefined? ------#

accuracy_with_fallback_grouped %>%
  filter(is.na(Macro_F1)) %>%
  summarise(
    n_rows = n(),
    one_truth_class = sum(Truth_Classes_Present == 1),
    two_truth_classes = sum(Truth_Classes_Present == 2),
    three_truth_classes = sum(Truth_Classes_Present == 3),
    
    one_pred_class = sum(Pred_Classes_Present == 1),
    two_pred_classes = sum(Pred_Classes_Present == 2),
    three_pred_classes = sum(Pred_Classes_Present == 3)
  ) %>%
  print(width = Inf)

# n_rows one_truth_class two_truth_classes three_truth_classes one_pred_class two_pred_classes
#     51              28                23                   0             25               26

# three_pred_classes
#                  0

#------ How many fallback rows have extreme BA values? ------#

accuracy_with_fallback_grouped_cons %>%
  filter(
    acc_source == "Balanced_Accuracy",
    acc %in% c(0, 1)
  ) %>%
  count(acc)

#   acc     n
#     0    51

#------ Are extreme fallback values concentrated in small n? ------#

accuracy_with_fallback_grouped_cons %>%
  filter(
    acc_source == "Balanced_Accuracy",
    acc %in% c(0, 1)
  ) %>%
  summarise(
    n = n(),
    min_n = min(N_evaluated),
    median_n = median(N_evaluated),
    max_n = max(N_evaluated)
  )

#   n min_n median_n max_n
#  51     1       32   239

#------ Distribution of n among fallback rows ------#

accuracy_with_fallback_grouped_cons %>%
  filter(
    acc_source == "Balanced_Accuracy",
    acc %in% c(0, 1)
  ) %>%
  count(N_evaluated, sort = TRUE) %>%
  print(n = Inf)

#   N_evaluated     n
# 1          40     8
# 2          11     6
# 3          27     4
# 4          32     4
# 5          45     4
# 6          21     3
# 7           1     2
# 8          20     2
# 9          34     2
# 10          43     2
# 11           2     1
# 12          22     1
# 13          24     1
# 14          26     1
# 15          30     1
# 16          31     1
# 17          35     1
# 18          38     1
# 19          41     1
# 20          51     1
# 21          60     1
# 22          74     1
# 23         202     1
# 24         239     1

#------ Inspect all fallback rows ------#

accuracy_with_fallback_grouped_cons %>%
  filter(acc_source == "Balanced_Accuracy") %>%
  select(
    raster,
    Ecds_ID,
    Ecrgn_n,
    Ecprv,
    source_level,
    N_evaluated,
    Truth_Classes_Present,
    Pred_Classes_Present,
    acc,
    acc_source,
    Macro_F1,
    Balanced_Accuracy,
    Overall_Accuracy,
    Weighted_F1,
    Support_0,
    Support_1,
    Support_2,
    F1_0,
    F1_1,
    F1_2
  ) %>%
  arrange(
    desc(N_evaluated),
    raster,
    Ecds_ID
  ) %>%
  print(n = Inf, width = Inf)

#   raster Ecds_ID Ecrgn_n                     Ecprv source_level N_evaluated Truth_Classes_Present
# 1 ngi        731 Aspen Parkland               10.2 ecodistrict          239                     1
# 2 ngi        740 Aspen Parkland               10.2 ecodistrict          202                     1
# 3 cnvc      1017 Northern Continental Divide  14.4 ecodistrict           74                     2
# 4 pudl       678 Boreal Transition             9.2 ecodistrict           60                     2
# 5 pudl       727 Aspen Parkland               10.2 ecodistrict           51                     2
# 6 ngi        683 Boreal Transition             9.2 ecoregion             45                     2
# 7 ngi        692 Boreal Transition             9.2 ecoregion             45                     2
# 8 ngi        701 Boreal Transition             9.2 ecoregion             45                     2
# 9 ngi        708 Boreal Transition             9.2 ecoregion             45                     2
# 10 pudl       684 Boreal Transition             9.2 ecodistrict           43                     2
# 11 pudl       760 Aspen Parkland               10.2 ecodistrict           43                     1
# 12 cnvc       380 Lake of the Woods             6.5 ecoregion             41                     1
# 13 cnvc       622 Western Alberta Upland        9.1 ecoregion             40                     2
# 14 cnvc       623 Western Alberta Upland        9.1 ecoregion             40                     2
# 15 cnvc       628 Western Alberta Upland        9.1 ecoregion             40                     2
# 16 cnvc       629 Western Alberta Upland        9.1 ecoregion             40                     2
# 17 cnvc       630 Western Alberta Upland        9.1 ecoregion             40                     2
# 18 cnvc       631 Western Alberta Upland        9.1 ecoregion             40                     2
# 19 gpw        601 Mid-Boreal Uplands            9.2 ecodistrict           40                     1
# 20 pudl       592 Peace Lowland                 9.2 ecodistrict           40                     2
# 21 eccc      1017 Northern Continental Divide  14.4 ecodistrict           38                     1
# 22 pudl       703 Boreal Transition             9.2 ecodistrict           35                     2
# 23 ngi        727 Aspen Parkland               10.2 ecodistrict           34                     1
# 24 ngi        744 Aspen Parkland               10.2 ecodistrict           34                     1
# 25 cnvc       375 Lake of the Woods             6.5 ecodistrict           32                     1
# 26 pudl       686 Boreal Transition             9.2 ecodistrict           32                     2
# 27 pudl       690 Boreal Transition             9.2 ecodistrict           32                     1
# 28 pudl       728 Aspen Parkland               10.2 ecodistrict           32                     2
# 29 pudl       692 Boreal Transition             9.2 ecodistrict           31                     2
# 30 gpw        628 Western Alberta Upland        9.1 ecodistrict           30                     1
# 31 cnvc       678 Boreal Transition             9.2 ecodistrict           27                     2
# 32 ngi        684 Boreal Transition             9.2 ecodistrict           27                     1
# 33 ngi        746 Aspen Parkland               10.2 ecodistrict           27                     1
# 34 pudl       755 Aspen Parkland               10.2 ecodistrict           27                     1
# 35 gpw        382 Thunder Bay-Quetico           6.5 ecodistrict           26                     1
# 36 cnvc       847 Lake Manitoba Plain          10.1 ecodistrict           24                     1
# 37 cnvc       684 Boreal Transition             9.2 ecodistrict           22                     2
# 38 cnvc       997 Eastern Continental Ranges   14.4 ecoregion             21                     1
# 39 cnvc       999 Eastern Continental Ranges   14.4 ecoregion             21                     1
# 40 pudl       590 Peace Lowland                 9.2 ecodistrict           21                     2
# 41 cnvc       839 Aspen Parkland               10.2 ecodistrict           20                     1
# 42 pudl       839 Aspen Parkland               10.2 ecodistrict           20                     1
# 43 cgr        375 Lake of the Woods             6.5 ecoprovince           11                     2
# 44 cgr        687 Boreal Transition             9.2 ecoprovince           11                     1
# 45 cgr        688 Boreal Transition             9.2 ecoprovince           11                     1
# 46 cgr        695 Boreal Transition             9.2 ecoprovince           11                     1
# 47 cgr        701 Boreal Transition             9.2 ecoprovince           11                     1
# 48 cgr        709 Boreal Transition             9.2 ecoprovince           11                     1
# 49 cnvc       960 Cascade Ranges               13.2 ecoprovince            2                     1
# 50 ngi        623 Western Alberta Upland        9.1 ecoprovince            1                     1
# 51 ngi        850 Lake Manitoba Plain          10.1 ecoprovince            1                     1

#   Pred_Classes_Present   acc acc_source        Macro_F1 Balanced_Accuracy Overall_Accuracy Weighted_F1
# 1                    2     0 Balanced_Accuracy      NaN                 0                0           0
# 2                    2     0 Balanced_Accuracy      NaN                 0                0           0
# 3                    2     0 Balanced_Accuracy      NaN                 0                0           0
# 4                    2     0 Balanced_Accuracy      NaN                 0                0           0
# 5                    1     0 Balanced_Accuracy      NaN                 0                0           0
# 6                    2     0 Balanced_Accuracy      NaN                 0                0           0
# 7                    2     0 Balanced_Accuracy      NaN                 0                0           0
# 8                    2     0 Balanced_Accuracy      NaN                 0                0           0
# 9                    2     0 Balanced_Accuracy      NaN                 0                0           0
# 10                    1     0 Balanced_Accuracy      NaN                 0                0           0
# 11                    2     0 Balanced_Accuracy      NaN                 0                0           0
# 12                    2     0 Balanced_Accuracy      NaN                 0                0           0
# 13                    1     0 Balanced_Accuracy      NaN                 0                0           0
# 14                    1     0 Balanced_Accuracy      NaN                 0                0           0
# 15                    1     0 Balanced_Accuracy      NaN                 0                0           0
# 16                    1     0 Balanced_Accuracy      NaN                 0                0           0
# 17                    1     0 Balanced_Accuracy      NaN                 0                0           0
# 18                    1     0 Balanced_Accuracy      NaN                 0                0           0
# 19                    1     0 Balanced_Accuracy      NaN                 0                0           0
# 20                    2     0 Balanced_Accuracy      NaN                 0                0           0
# 21                    1     0 Balanced_Accuracy      NaN                 0                0           0
# 22                    2     0 Balanced_Accuracy      NaN                 0                0           0
# 23                    2     0 Balanced_Accuracy      NaN                 0                0           0
# 24                    2     0 Balanced_Accuracy      NaN                 0                0           0
# 25                    2     0 Balanced_Accuracy      NaN                 0                0           0
# 26                    2     0 Balanced_Accuracy      NaN                 0                0           0
# 27                    2     0 Balanced_Accuracy      NaN                 0                0           0
# 28                    2     0 Balanced_Accuracy      NaN                 0                0           0
# 29                    1     0 Balanced_Accuracy      NaN                 0                0           0
# 30                    1     0 Balanced_Accuracy      NaN                 0                0           0
# 31                    1     0 Balanced_Accuracy      NaN                 0                0           0
# 32                    2     0 Balanced_Accuracy      NaN                 0                0           0
# 33                    2     0 Balanced_Accuracy      NaN                 0                0           0
# 34                    1     0 Balanced_Accuracy      NaN                 0                0           0
# 35                    2     0 Balanced_Accuracy      NaN                 0                0           0
# 36                    2     0 Balanced_Accuracy      NaN                 0                0           0
# 37                    1     0 Balanced_Accuracy      NaN                 0                0           0
# 38                    2     0 Balanced_Accuracy      NaN                 0                0           0
# 39                    2     0 Balanced_Accuracy      NaN                 0                0           0
# 40                    1     0 Balanced_Accuracy      NaN                 0                0           0
# 41                    2     0 Balanced_Accuracy      NaN                 0                0           0
# 42                    1     0 Balanced_Accuracy      NaN                 0                0           0
# 43                    2     0 Balanced_Accuracy      NaN                 0                0           0
# 44                    1     0 Balanced_Accuracy      NaN                 0                0           0
# 45                    1     0 Balanced_Accuracy      NaN                 0                0           0
# 46                    1     0 Balanced_Accuracy      NaN                 0                0           0
# 47                    1     0 Balanced_Accuracy      NaN                 0                0           0
# 48                    1     0 Balanced_Accuracy      NaN                 0                0           0
# 49                    1     0 Balanced_Accuracy      NaN                 0                0           0
# 50                    1     0 Balanced_Accuracy      NaN                 0                0           0
# 51                    1     0 Balanced_Accuracy      NaN                 0                0           0

#   Support_0 Support_1 Support_2  F1_0  F1_1  F1_2
# 1         0         0       239    NA    NA    NA
# 2         0         0       202    NA    NA    NA
# 3         0         2        72    NA    NA    NA
# 4         0        53         7    NA    NA    NA
# 5         0        17        34    NA    NA    NA
# 6         1         0        44    NA    NA    NA
# 7         1         0        44    NA    NA    NA
# 8         1         0        44    NA    NA    NA
# 9         1         0        44    NA    NA    NA
# 10         0        14        29    NA    NA    NA
# 11         0         0        43    NA    NA    NA
# 12         0         0        41    NA    NA    NA
# 13         0         3        37    NA    NA    NA
# 14         0         3        37    NA    NA    NA
# 15         0         3        37    NA    NA    NA
# 16         0         3        37    NA    NA    NA
# 17         0         3        37    NA    NA    NA
# 18         0         3        37    NA    NA    NA
# 19        40         0         0    NA    NA    NA
# 20         0        36         4    NA    NA    NA
# 21         0         0        38    NA    NA    NA
# 22         0         8        27    NA    NA    NA
# 23         0         0        34    NA    NA    NA
# 24         0         0        34    NA    NA    NA
# 25         0         0        32    NA    NA    NA
# 26         0         1        31    NA    NA    NA
# 27         0         0        32    NA    NA    NA
# 28         0         3        29    NA    NA    NA
# 29         0        16        15    NA    NA    NA
# 30         0        30         0    NA    NA    NA
# 31         0        24         3    NA    NA    NA
# 32         0         0        27    NA    NA    NA
# 33         0         0        27    NA    NA    NA
# 34         0         0        27    NA    NA    NA
# 35         0        26         0    NA    NA    NA
# 36         0         0        24    NA    NA    NA
# 37         0         3        19    NA    NA    NA
# 38         0         0        21    NA    NA    NA
# 39         0         0        21    NA    NA    NA
# 40         0        20         1    NA    NA    NA
# 41         0         0        20    NA    NA    NA
# 42         0         0        20    NA    NA    NA
# 43         0         4         7    NA    NA    NA
# 44         0         0        11    NA    NA    NA
# 45         0         0        11    NA    NA    NA
# 46         0         0        11    NA    NA    NA
# 47         0         0        11    NA    NA    NA
# 48         0         0        11    NA    NA    NA
# 49         0         0         2    NA    NA    NA
# 50         0         0         1    NA    NA    NA
# 51         0         0         1    NA    NA    NA

accuracy_with_fallback_grouped_cons %>%
  filter(acc_source == "Balanced_Accuracy") %>%
  summarise(
    no_class0 = sum(Support_0 == 0),
    no_class1 = sum(Support_1 == 0),
    no_class2 = sum(Support_2 == 0)
  )

# no_class0 no_class1 no_class2
#        46        30         3

#------ CLASS-SPECIFIC METRICS ------#

# Summary by class
class_metrics_long <- accuracy_with_fallback_grouped %>%
  select(
    raster,
    Ecds_ID,
    N_evaluated,
    starts_with("Precision_"),
    starts_with("Recall_"),
    starts_with("F1_"),
    starts_with("MCC_OvR_")
  ) %>%
  pivot_longer(
    cols = -c(raster, Ecds_ID, N_evaluated),
    names_to = c("metric", "class"),
    names_pattern = "(.*)_(\\d)$"
  )

class_metrics_long %>%
  group_by(class, metric) %>%
  summarise(
    n = sum(!is.na(value)),
    NAs = sum(is.na(value)),
    mean = mean(value, na.rm = TRUE),
    median = median(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    q25 = quantile(value, 0.25, na.rm = TRUE),
    q75 = quantile(value, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

# class metric        n  mean median    sd    q25   q75
# 1 0     F1         1016 0.722  0.837 0.280 0.563  0.943
# 2 0     MCC_OvR    1002 0.512  0.558 0.260 0.356  0.711
# 3 0     Precision  1595 0.446  0.434 0.419 0      0.914
# 4 0     Recall     1029 0.860  0.949 0.205 0.833  0.984
# 5 1     F1         1093 0.499  0.559 0.273 0.272  0.736
# 6 1     MCC_OvR    1165 0.351  0.359 0.283 0.101  0.590
# 7 1     Precision  1499 0.458  0.535 0.377 0      0.816
# 8 1     Recall     1237 0.454  0.469 0.312 0.16   0.733
# 9 2     F1          966 0.502  0.516 0.294 0.240  0.755
# 10 2     MCC_OvR     737 0.210  0.122 0.238 0.0333 0.392
# 11 2     Precision  1340 0.588  0.854 0.445 0      1    
# 12 2     Recall     1129 0.409  0.382 0.324 0.0954 0.694
# 13 NA    NA        50230 0.474  0.515 0.348 0.101  0.769

# Metric distributions among classes
class_metrics_long %>%
  filter(metric == "F1") %>%
  ggplot(aes(class, value)) +
  geom_boxplot() +
  labs(
    x = "Class",
    y = "F1 score",
    title = "Distribution of F1 by class"
  ) +
  theme_bw()

class_metrics_long %>%
  filter(metric == "Recall") %>%
  ggplot(aes(class, value)) +
  geom_boxplot() +
  theme_bw()

class_metrics_long %>%
  filter(metric == "Precision") %>%
  ggplot(aes(class, value)) +
  geom_boxplot() +
  theme_bw()

# How much support does each class have?
accuracy_with_fallback_grouped %>%
  summarise(
    median_support0 = median(Support_0),
    median_support1 = median(Support_1),
    median_support2 = median(Support_2),
    
    mean_support0 = mean(Support_0),
    mean_support1 = mean(Support_1),
    mean_support2 = mean(Support_2)
  )

# median_support0 median_support1 median_support2 mean_support0 mean_support1 mean_support2
#              33              47              33          330.          123.          351.

accuracy_with_fallback_grouped %>%
  select(
    Support_0,
    Support_1,
    Support_2
  ) %>%
  pivot_longer(
    everything(),
    names_to = "class",
    values_to = "support"
  ) %>%
  ggplot(aes(support)) +
  geom_histogram(bins = 40) +
  facet_wrap(~ class, scales = "free_y") +
  theme_bw()

# How often does each class drive Macro F1 failures?
accuracy_with_fallback_grouped %>%
  filter(is.na(Macro_F1)) %>%
  summarise(
    no_class0 = sum(Support_0 == 0),
    no_class1 = sum(Support_1 == 0),
    no_class2 = sum(Support_2 == 0)
  )

# no_class0 no_class1 no_class2
#        46        30         3

accuracy_with_fallback_grouped %>%
  filter(is.na(Macro_F1)) %>%
  mutate(
    missing_pattern = paste0(
      ifelse(Support_0 == 0, "0", ""),
      ifelse(Support_1 == 0, "1", ""),
      ifelse(Support_2 == 0, "2", "")
    )
  ) %>%
  count(missing_pattern, sort = TRUE)

#   missing_pattern     n
# 1 01                 25
# 2 0                  19
# 3 1                   4
# 4 02                  2
# 5 12                  1

# Which rasters struggle with native grassland?

# Average native-class F1 by source
accuracy_with_fallback_grouped %>%
  group_by(raster) %>%
  summarise(
    mean_F1_native = mean(F1_2, na.rm = TRUE),
    median_F1_native = median(F1_2, na.rm = TRUE),
    n_valid = sum(!is.na(F1_2))
  ) %>%
  arrange(mean_F1_native)

#   raster mean_F1_native median_F1_native n_valid
# 1 gpw             0.197            0.133      46
# 2 aci             0.246            0.136      52
# 3 plvi            0.316            0.347      13
# 4 mgi             0.439            0.426      32
# 5 cnvc            0.484            0.541     175
# 6 pudl            0.498            0.499     148
# 7 pli             0.505            0.495      67
# 8 eccc            0.555            0.561     131
# 9 cgr             0.578            0.660     118
# 10 pglc           0.606            0.621      83
# 11 ngi            0.607            0.702     101

# Average tame-class F1 by source
accuracy_with_fallback_grouped %>%
  group_by(raster) %>%
  summarise(
    mean_F1_tame = mean(F1_1, na.rm = TRUE),
    median_F1_tame = median(F1_1, na.rm = TRUE),
    n_valid = sum(!is.na(F1_1))
  ) %>%
  arrange(mean_F1_tame)

#   raster mean_F1_tame median_F1_tame n_valid
# 1 pudl         0.0547         0.0401      50
# 2 ngi          0.0952         0.0547      30
# 3 cnvc         0.119          0.127       42
# 4 cgr          0.154          0.166       54
# 5 plvi         0.224          0.207        7
# 6 mgi          0.397          0.377       39
# 7 pli          0.411          0.412       45
# 8 gpw          0.432          0.420      271
# 9 aci          0.625          0.692      275
# 10 pglc        0.679          0.680       86
# 11 eccc        0.742          0.770      194

# Identify systematic omission vs commission errors
# Recall < Precision → class is being omitted.
# Precision < Recall → class is being over-predicted.
accuracy_with_fallback_grouped %>%
  summarise(
    native_precision = mean(Precision_2, na.rm = TRUE),
    native_recall = mean(Recall_2, na.rm = TRUE),
    
    tame_precision = mean(Precision_1, na.rm = TRUE),
    tame_recall = mean(Recall_1, na.rm = TRUE),
    
    not_precision = mean(Precision_0, na.rm = TRUE),
    not_recall = mean(Recall_0, na.rm = TRUE)
  )

# native_precision native_recall tame_precision tame_recall not_precision not_recall
#            0.588         0.409          0.458       0.454         0.446      0.860

# Relationship between sample size and class-specific performance
accuracy_with_fallback_grouped %>%
  ggplot(
    aes(
      N_evaluated,
      F1_2
    )
  ) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "loess") +
  scale_x_log10() +
  labs(
    x = "N evaluated",
    y = "F1 Native",
    title = "Native F1 versus sample size"
  ) +
  theme_bw()

# Quantify whether native grasslands are intrinsically harder
# If class 2 remains substantially lower even after excluding NA values and accounting 
# for sample size, that would support a conclusion that native grassland is genuinely 
# the most difficult class for the underlying datasets.
accuracy_with_fallback_grouped %>%
  summarise(
    mean_F1_0 = mean(F1_0, na.rm = TRUE),
    mean_F1_1 = mean(F1_1, na.rm = TRUE),
    mean_F1_2 = mean(F1_2, na.rm = TRUE),
    
    median_F1_0 = median(F1_0, na.rm = TRUE),
    median_F1_1 = median(F1_1, na.rm = TRUE),
    median_F1_2 = median(F1_2, na.rm = TRUE)
  )

# Sample size plot
accuracy_with_fallback_grouped %>%
  summarise(
    cor_native =
      cor(log10(N_evaluated), F1_2,
          use = "complete.obs"),
    
    cor_tame =
      cor(log10(N_evaluated), F1_1,
          use = "complete.obs"),
    
    cor_not =
      cor(log10(N_evaluated), F1_0,
          use = "complete.obs")
  )

# cor_native cor_tame cor_not
#      0.108   -0.202  0.0477

# ----------------------------
# 16) Export file for consensus mapping
# ----------------------------

# Primary metric: Macro F1 
# Fallback metric: Balanced Accuracy 

# Grouped
accuracy_with_fallback_grouped_cons <- accuracy_with_fallback_grouped %>%
  mutate(
    acc = coalesce(Macro_F1, Balanced_Accuracy),
    acc_source = if_else(
      !is.na(Macro_F1),
      "Macro_F1",
      "Balanced_Accuracy"
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
    acc = coalesce(Macro_F1, Balanced_Accuracy),
    acc_source = if_else(
      !is.na(Macro_F1),
      "Macro_F1",
      "Balanced_Accuracy"
    )
  )

# Add quality flag for rows with N less than 20
accuracy_with_fallback_atomic_cons <- accuracy_with_fallback_atomic_cons %>%
  mutate(
    low_confidence_acc =
      #acc_source == "Overall_Accuracy" &
      N_evaluated < 20
  )

# Class-specific metrics (class-specific F1 or fallbacks if NA)
accuracy_with_fallback_grouped_cons <- accuracy_with_fallback_grouped_cons %>%
  mutate(
    native_metric = coalesce(F1_2, Macro_F1, Balanced_Accuracy),
    tame_metric   = coalesce(F1_1, Macro_F1, Balanced_Accuracy),
    not_metric    = coalesce(F1_0, Macro_F1, Balanced_Accuracy))

accuracy_with_fallback_atomic_cons <- accuracy_with_fallback_atomic_cons %>%
  mutate(
    native_metric = coalesce(F1_2, Macro_F1, Balanced_Accuracy),
    tame_metric   = coalesce(F1_1, Macro_F1, Balanced_Accuracy),
    not_metric    = coalesce(F1_0, Macro_F1, Balanced_Accuracy))

#---- Export ----

# GROUPED
readr::write_csv(
  accuracy_with_fallback_grouped_cons,
  file.path(out_dir, "nattam_accuracy_grouped_combinedCIandCM_forconsensus_22Jul2026.csv")
)

# ATOMIC
readr::write_csv(
  accuracy_with_fallback_atomic_cons,
  file.path(out_dir, "nattam_accuracy_atomic_combinedCIandCM_forconsensus_22Jul2026.csv")
)

# ----------------------------
# 17) Create rasters of accuracy by ecodistrict
# ----------------------------

summary_by_ecodistrict <- summarize_block_multiclass(
  combined,
  key_cols = c("raster", "Ecds_ID")
)

# ---- global MCC fallback (FIXED) ----
global_mcc <- metrics_by_dataset %>%
  select(raster_group, MCC_global = MCC_Multiclass) %>%
  rename(raster = raster_group)

# ---- weights ----
weights_ecodist_mcc <- summary_by_ecodistrict %>%
  select(raster, Ecds_ID, MCC_Multiclass, N_evaluated) %>%
  
  rename(
    MCC = MCC_Multiclass,
    N_evaluated = N_evaluated
  ) %>%
  
  mutate(
    MCC = ifelse(N_evaluated < 30, NA, MCC)
  ) %>%
  
  left_join(global_mcc, by = "raster") %>%
  
  mutate(
    MCC_final = case_when(
      is.na(MCC) ~ MCC_global,
      TRUE ~ 0.7 * MCC + 0.3 * MCC_global
    ),
    
    MCC_final = (MCC_final + 1) / 2,
    MCC_final = pmin(pmax(MCC_final, 0), 1)
  ) %>%
  
  select(raster, Ecds_ID, MCC_final)

library(terra)

eco_vect <- vect("R:/CanadianGrasslands/processeddatasets/boundaries/national-ecological-framework-for-canada/Ecodistricts/ecodistricts_clean_class_sa.shp")
eco_vect$Ecds_ID <- as.character(eco_vect$ECODISTRIC)

weights_ecodist_mcc$Ecds_ID <- as.character(weights_ecodist_mcc$Ecds_ID)

template <- rast("R:/CanadianGrasslands/WG_Classification/processed/aafc_aci/aligned/aafc_aci_merged_reproj_algn_binary_noshrub.tif")

raster_names <- unique(weights_ecodist_mcc$raster)

mcc_rasters <- vector("list", length(raster_names))
names(mcc_rasters) <- raster_names

for (r_name in raster_names) {
  
  w_df <- weights_ecodist_mcc %>% filter(raster == r_name)
  
  eco_joined <- merge(eco_vect, w_df, by = "Ecds_ID", all.x = TRUE)
  
  fallback_val <- global_mcc$MCC_global[global_mcc$raster == r_name]
  fallback_val <- (fallback_val + 1) / 2
  
  eco_joined$MCC_final[is.na(eco_joined$MCC_final)] <- fallback_val
  
  r_mcc <- rasterize(
    eco_joined,
    template,
    field = "MCC_final",
    touches = TRUE
  )
  
  names(r_mcc) <- r_name
  
  mcc_rasters[[r_name]] <- r_mcc
}

mcc_stack <- rast(mcc_rasters)

writeRaster(
  mcc_stack,
  file.path(out_dir, "mcc_weight_stack_triclass.tif"),
  overwrite = TRUE,
  wopt = list(datatype = "FLT4S",
              gdal = c("COMPRESS=ZSTD", "TILED=YES", "BIGTIFF=YES"))
)

# ----------------------------
# 15) Plots: Map by ecodistrict
# ----------------------------

# ---- Match CRS to template raster ----
ecodistricts <- terra::project(eco_vect, template)

# ---- Convert to sf for ggplot ----
ecodistricts_sf <- sf::st_as_sf(ecodistricts)

# ---- Ensure join key is character ----
ecodistricts_sf <- ecodistricts_sf %>%
  mutate(Ecds_ID = as.character(Ecds_ID))

accuracy_by_ecodistrict_grouped <- accuracy_by_ecodistrict_grouped %>%
  mutate(Ecds_ID = as.character(Ecds_ID))

# ----------------------------
# Join spatial + accuracy data
# ----------------------------
map_data_grouped <- ecodistricts_sf %>%
  left_join(accuracy_by_ecodistrict_grouped, by = "Ecds_ID")

# ----------------------------
# Metrics to map
# ----------------------------
metrics_to_map <- c(
  "Overall_Accuracy",
  "Balanced_Accuracy",
  "Kappa",
  "MCC_Multiclass",
  "N_evaluated"
)

# ----------------------------
# Function to make maps
# ----------------------------
make_map <- function(data, raster_name, metric_name) {
  
  df <- data %>%
    filter(raster == raster_name)
  
  # Define fill scale 
  if (metric_name == "N_evaluated") {
    fill_scale <- scale_fill_viridis_c(
      na.value = "grey90",
      option = "C",
      name = "N evaluated"
    )
  } else {
    fill_scale <- scale_fill_viridis_c(
      na.value = "grey90",
      option = "C",
      labels = scales::percent_format(accuracy = 1),
      name = metric_name
    )
  }
  
  # Build plot
  p <- ggplot(df) +
    geom_sf(aes(fill = .data[[metric_name]]), color = NA) +
    fill_scale +
    ggtitle(paste0("Raster: ", raster_name, " | ", metric_name)) +
    theme_minimal()
  
  return(p)
}

# ----------------------------
# Generate and save maps
# ----------------------------
unique_rasters <- sort(unique(map_data_grouped$raster))

for (r in unique_rasters) {
  
  for (m in metrics_to_map) {
    
    p <- make_map(map_data_grouped, r, m)
    
    outfile <- file.path(
      out_dir,
      paste0("map_", r, "_", m, ".png")
    )
    
    ggsave(
      filename = outfile,
      plot = p,
      width = 8,
      height = 6,
      dpi = 300
    )
  }
}

# ------ Mean N_evaluated per ecodistrict (across rasters)------

# ---- Compute mean N across rasters ----
mean_N_by_ecodistrict <- accuracy_by_ecodistrict_atomic %>%
  group_by(Ecds_ID) %>%
  summarise(
    mean_N_evaluated = mean(N_evaluated, na.rm = TRUE),
    n_rasters = sum(!is.na(N_evaluated)),
    .groups = "drop"
  )

# ---- Join to spatial data ----
mean_N_by_ecodistrict <- mean_N_by_ecodistrict %>%
  mutate(Ecds_ID = as.character(Ecds_ID))

mean_N_map <- ecodistricts_sf %>%
  left_join(mean_N_by_ecodistrict, by = "Ecds_ID")

# --- Plot map ------------
p_mean_N <- ggplot(mean_N_map) +
  geom_sf(aes(fill = mean_N_evaluated), color = NA) +
  scale_fill_viridis_c(
    option = "C",
    na.value = "grey90",
    name = "Mean N evaluated"
  ) +
  ggtitle("Mean Sample Size (N_evaluated) Across Rasters by Ecodistrict") +
  theme_minimal()

# ---- Save ----
ggsave(
  filename = file.path(out_dir, "map_mean_N_evaluated_ecodistrict.png"),
  plot = p_mean_N,
  width = 8,
  height = 6,
  dpi = 300
)

# ----------------------------
# 16) Plots: Overall accuracy (GROUPED and ATOMIC)
# ----------------------------

# grouped summary 
summary_by_raster_grouped <- metrics_by_dataset

# atomic summary
summary_by_raster_atomic <- summarize_block_multiclass(
  combined,
  key_cols = "raster"
)

# GROUPED plot
p_grouped <- ggplot(summary_by_raster_grouped,
                    aes(x = reorder(raster_group, Overall_Accuracy), y = Overall_Accuracy)) +
  geom_col(fill = "#2d7f5e") +
  geom_text(aes(label = paste0("n = ", N_evaluated)),
            hjust = 1.1,
            color = "white",
            size = 3.5) +
  coord_flip() +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.02))) +
  labs(title = "Overall Accuracy by raster",
       x = "Raster", y = "Overall Accuracy") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(out_dir, "raster_accuracy_plot_nattam_cropsins.png"),
       p_grouped, width = 7, height = 5, dpi = 200)

# ATOMIC plot
p_atomic <- ggplot(summary_by_raster_atomic,
                   aes(x = reorder(raster, Overall_Accuracy), y = Overall_Accuracy)) +
  geom_col(fill = "#2d7f5e") +
  geom_text(aes(label = paste0("n = ", N_evaluated)),
            hjust = 1.1,
            color = "white",
            size = 3.5) +
  coord_flip() +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.02))) +
  labs(title = "Overall Accuracy by Raster",
       x = "Raster", y = "Overall Accuracy") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(out_dir, "raster_accuracy_plot_atomic_nattam_cropins.png"),
       p_atomic, width = 7, height = 5, dpi = 200)

#################################################################
#################################################################

# Compare macro metrics
summary_long <- summary_by_raster_grouped %>%
  select(raster_group,
         Overall_Accuracy,
         Balanced_Accuracy,
         Macro_F1,
         MCC_Multiclass) %>%
  pivot_longer(-raster_group, names_to = "metric", values_to = "value")

ggplot(summary_long,
       aes(x = reorder(raster_group, value), y = value, fill = metric)) +
  geom_col(position = "dodge") +
  coord_flip() +
  theme_minimal()