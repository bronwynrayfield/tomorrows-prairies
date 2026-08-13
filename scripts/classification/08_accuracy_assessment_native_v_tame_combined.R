
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

acc_cm_file <- file.path(acc_cm_dir, "27Jul2026/eval_data_27Jul2026_natvtam_withnongrass_compiled.csv")
acc_ci_file <- file.path(acc_ci_dir, "23Jul2026/eval_data_23Jul2026_natvtam_withnongrass_lenient.csv")
#acc_ci_file <- file.path(acc_ci_dir, "23Jul2026/eval_data_23Jul2026_natvtam_withnongrass.csv")
metadata_file <- "R:/CanadianGrasslands/WG_Classification/github/tomorrows-prairies/processeddatasets/Grass/metadata_full_extent_nattam_pli_sep.csv"

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

# Count ACI_field_2024 point/raster combos in acc_cm (5,463)
acc_cm %>%
  filter(Sorc_ID == "ACI_field_2024") %>% 
  tally()

# Count ACI_field_2024 point/raster combos in acc_cm (10,477) - larger due to lenient year matching for tame/other
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

#   Pt_ID Ecds_ID Ecrgn_n         Ecprv Training Sorc_ID     LC_class3 Yer_Cll LC_Cl_N raster  pred
# 1 691703 298163     745 Aspen Parkland   10.2 NA       ACI_CI_2020         0    2020 OTHER   pli_mm     0
# 2 703523 303882     830 Mixed Grassland  10.3 NA       ACI_CI_2020         0    2020 OTHER   pli_mm     0

#   raster_group allowed_years matched_year classpair Ecrg_ID
# 1 pli          2020;2021             2020 T0_P0          NA
# 2 pli          2020;2021             2020 T0_P0          NA

# Remove them from the dataset
acc_ci <- acc_ci %>%
  filter(!(Ecds_ID %in% c(745, 830) & raster == "pli_mm"))

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
readr::write_csv(combined_out, file.path(out_dir, "accuracy_assessment_data_natvtam_pli_sep_CIandCM_CIcombined_29Jul2026.csv"))

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
out_xlsx <- file.path(out_dir, "nat_v_tame_v_nongrass_accuracy_assessment_results_combined_29Jul2026.xlsx")

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
    "nattam_accuracy_by_ecodistrict_atomic_ecoregprovfallback_combined_29Jul2026.csv"
  )
)

readr::write_csv(
  accuracy_with_fallback_grouped,
  file.path(
    out_dir,
    "nattam_accuracy_by_ecodistrict_grouped_ecoregprovfallback_combined_29Jul2026.csv"
  )
)

# ----------------------------
# 15) Diagnostics
# ----------------------------
table(combined$raster)

# lenient:
# aci2024 aci2025     cgr    cnvc    eccc     gpw     mgi     ngi    pglc  pli_ap  pli_cu  pli_mg  pli_mm 
#  265807    5884   28144   39733   60363  270552    7336   22415   50713     710     920   20638    8319 
# plvi    pudl 
#  602   37736 

# exact match:
# aci2024 aci2025     cgr    cnvc    eccc     gpw     mgi     ngi    pglc  pli_ap  pli_cu  pli_mg  pli_mm 
# 610     869   28144   39733   27340    5808    7336   22415   50713     710     920   13711    8319 
# plvi    pudl 
# 602   37736 

table(combined$raster_group)

# lenient:
# aci2024 aci2025     cgr    cnvc    eccc     gpw     mgi     ngi    pglc     pli    plvi    pudl 
#  265807    5884   28144   39733   60363  270552    7336   22415   50713   30587     602   37736

# exact match:
# aci2024 aci2025     cgr    cnvc    eccc     gpw     mgi     ngi    pglc     pli    plvi    pudl 
# 610     869   28144   39733   27340    5808    7336   22415   50713   23660     602   37736 

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

# without ACI2025:
# n_rows OA_ecoD_NA OA_ecoR_NA OA_ecoP_NA BA_ecoD_NA BA_ecoR_NA BA_ecoP_NA Kappa_ecoD_NA Kappa_ecoR_NA
#   1605          0          0          0          0          0          0           524           102
# Kappa_ecoP_NA MCC_ecoD_NA MCC_ecoR_NA MCC_ecoP_NA
#            43         553         117          49

# with ACI2025:
# n_rows OA_ecoD_NA OA_ecoR_NA OA_ecoP_NA BA_ecoD_NA BA_ecoR_NA BA_ecoP_NA Kappa_ecoD_NA Kappa_ecoR_NA
#   1693          0          0          0          0          0          0           632           216
# Kappa_ecoP_NA MCC_ecoD_NA MCC_ecoR_NA MCC_ecoP_NA
#           156         661         231         162

# exact match (not lenient):
# n_rows OA_ecoD_NA OA_ecoR_NA OA_ecoP_NA BA_ecoD_NA BA_ecoR_NA BA_ecoP_NA Kappa_ecoD_NA Kappa_ecoR_NA
#   1200          0          0          0          0          0          0           754           311
# Kappa_ecoP_NA MCC_ecoD_NA MCC_ecoR_NA MCC_ecoP_NA
#           250         777         326         256

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

# without ACI2025:
# n_rows OA_ecoD_NA OA_ecoR_NA OA_ecoP_NA BA_ecoD_NA BA_ecoR_NA BA_ecoP_NA Kappa_ecoD_NA Kappa_ecoR_NA
#   1607          0          0          0          0          0          0           527           124
# Kappa_ecoP_NA MCC_ecoD_NA MCC_ecoR_NA MCC_ecoP_NA
#            66         556         139          72

# with ACI2025:
# n_rows OA_ecoD_NA OA_ecoR_NA OA_ecoP_NA BA_ecoD_NA BA_ecoR_NA BA_ecoP_NA Kappa_ecoD_NA Kappa_ecoR_NA
#   1693          0          0          0          0          0          0           632           216
# Kappa_ecoP_NA MCC_ecoD_NA MCC_ecoR_NA MCC_ecoP_NA
#           158         661         231         164

# exact match (not lenient):
# n_rows OA_ecoD_NA OA_ecoR_NA OA_ecoP_NA BA_ecoD_NA BA_ecoR_NA BA_ecoP_NA Kappa_ecoD_NA Kappa_ecoR_NA
#   1200          0          0          0          0          0          0           754           311
# Kappa_ecoP_NA MCC_ecoD_NA MCC_ecoR_NA MCC_ecoP_NA
#           277         777         326         283

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

# without ACI2025:
# n_rows OA_NA OA_pct_NA Macro_F1_NA Macro_F1_pct_NA Weighted_F1_NA Weighted_F1_pct_NA BA_NA BA_pct_NA
#   1605     0         0          51            3.18              0                  0     0         0
# Kappa_NA Kappa_pct_NA MCC_NA MCC_pct_NA
#      346         21.6    361       22.5

# with ACI2025:
# n_rows OA_NA OA_pct_NA Macro_F1_NA Macro_F1_pct_NA Weighted_F1_NA Weighted_F1_pct_NA BA_NA BA_pct_NA
#   1693     0         0          52            3.07              0                  0     0         0
# Kappa_NA Kappa_pct_NA MCC_NA MCC_pct_NA
#      451         26.6    466       27.5

# exact match (not lenient):
# n_rows OA_NA OA_pct_NA Macro_F1_NA Macro_F1_pct_NA Weighted_F1_NA Weighted_F1_pct_NA BA_NA BA_pct_NA
#   1200     0         0          58            4.83              0                  0     0         0
# Kappa_NA Kappa_pct_NA MCC_NA MCC_pct_NA
#      587         48.9    599       49.9

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

# without ACI2025:
# n_rows OA_NA OA_pct_NA Macro_F1_NA Macro_F1_pct_NA Weighted_F1_NA Weighted_F1_pct_NA BA_NA BA_pct_NA
#   1607     0         0          51            3.17              0                  0     0         0
# Kappa_NA Kappa_pct_NA MCC_NA MCC_pct_NA
#      359         22.3    374       23.3

# with ACI2025:
# n_rows OA_NA OA_pct_NA Macro_F1_NA Macro_F1_pct_NA Weighted_F1_NA Weighted_F1_pct_NA BA_NA BA_pct_NA
#   1693     0         0          52            3.07              0                  0     0         0
# Kappa_NA Kappa_pct_NA MCC_NA MCC_pct_NA
#      451         26.6    466       27.5

#------ Where are MCC/Kappa missing? ------#

accuracy_with_fallback_grouped %>%
  count(Truth_Classes_Present)

# without ACI2025:
#  Truth_Classes_Present     n
#1                     1   346
#2                     2   728
#3                     3   531

# with ACI2025:
#   Truth_Classes_Present     n
# 1                     1   451
# 2                     2   692
# 3                     3   550

# exact match (not lenient):
#   Truth_Classes_Present     n
# 1                     1   587
# 2                     2   343
# 3                     3   270

accuracy_with_fallback_grouped %>%
  count(Pred_Classes_Present)

# without ACI2025:
#   Pred_Classes_Present     n
# 1                    1    33
# 2                    2   315
# 3                    3  1257

# with ACI2025:
# Pred_Classes_Present     n
# 1                    1    38
# 2                    2   357
# 3                    3  1298

# exact match (not lenient):
#   Pred_Classes_Present     n
# 1                    1    35
# 2                    2   205
# 3                    3   960

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

# without ACI2025:
# one_truth_class two_truth_classes three_truth_classes one_pred_class two_pred_classes three_pred_classes
#             346               728                 531             33              315               1257
# Kappa_NA MCC_NA
#      346    361

# with ACI2025:
# one_truth_class two_truth_classes three_truth_classes one_pred_class two_pred_classes three_pred_classes
#             451               692                 550             38              357               1298
# Kappa_NA MCC_NA
#      451    466

# exact match (not lenient):
# one_truth_class two_truth_classes three_truth_classes one_pred_class two_pred_classes three_pred_classes
#             587               343                 270             35              205                960
# Kappa_NA MCC_NA
#      587    599

# Which class is driving the 346/451 single-class units?
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

# without ACI2025:
# only_class0 only_class1 only_class2
#          24          11         311

# with ACI2025:
# only_class0 only_class1 only_class2
#          69          52         330

# exact match (not lenient):
# only_class0 only_class1 only_class2
#           6          44         537

# Which class pairings dominate the 728/692 two-class units?
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

# without ACI2025:
#   combo     n
# 1 0+1     441
# 2 0+2      33
# 3 1+2     254
# 4 NA      877

# with ACI2025:
#   combo     n
# 1 0+1     402
# 2 0+2      31
# 3 1+2     259
# 4 NA     1001

# exact match (not lenient):
#   combo     n
# 1 0+1       5
# 2 0+2      19
# 3 1+2     319
# 4 NA      857

#------ Which classes are absent? ------#

accuracy_with_fallback_grouped %>%
summarise(
  no_class0 = sum(Support_0 == 0, na.rm = TRUE),
  no_class1 = sum(Support_1 == 0, na.rm = TRUE),
  no_class2 = sum(Support_2 == 0, na.rm = TRUE)
)

# without ACI2025:
# no_class0 no_class1 no_class2
#       576       368       476

# with ACI2025:
# no_class0 no_class1 no_class2
#       641       430       523

# exact match (not lenient):
# no_class0 no_class1 no_class2
#       900       562        55

accuracy_with_fallback_atomic %>%
  summarise(
    no_class0 = sum(Support_0 == 0, na.rm = TRUE),
    no_class1 = sum(Support_1 == 0, na.rm = TRUE),
    no_class2 = sum(Support_2 == 0, na.rm = TRUE)
  )

# without ACI2025:
# no_class0 no_class1 no_class2
#       588       369       477

# with ACI2025:
# no_class0 no_class1 no_class2
#       641       430       523

# ----- Missing one-vs-rest MCCs ------#

accuracy_with_fallback_grouped %>%
  summarise(
    MCC0_NA = sum(is.na(MCC_OvR_0)),
    MCC1_NA = sum(is.na(MCC_OvR_1)),
    MCC2_NA = sum(is.na(MCC_OvR_2))
  )

# without ACI2025:
# MCC0_NA MCC1_NA MCC2_NA
#     603     440     868 - highest number of NAs for native

# with ACI2025:
# MCC0_NA MCC1_NA MCC2_NA
#     713     543     935

# exact match (not lenient):
# MCC0_NA MCC1_NA MCC2_NA
#     906     659     649 - highest number of NAs for non-grass

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

# without ACI2025:
# no_class0_truth no_class1_truth no_class2_truth no_class0_pred no_class1_pred no_class2_pred
#             576             368             476             10            106            265

# with ACI2025:
# no_class0_truth no_class1_truth no_class2_truth no_class0_pred no_class1_pred no_class2_pred
#             641             430             523             11            118            304

# exact match (not lenient):
# no_class0_truth no_class1_truth no_class2_truth no_class0_pred no_class1_pred no_class2_pred
#             900             562              55             16            132            127

accuracy_with_fallback_atomic %>%
  summarise(
    MCC0_NA = sum(is.na(MCC_OvR_0)),
    MCC1_NA = sum(is.na(MCC_OvR_1)),
    MCC2_NA = sum(is.na(MCC_OvR_2))
  )

# without ACI2025:
# MCC0_NA MCC1_NA MCC2_NA
#     616     441     881

# with ACI2025:
# MCC0_NA MCC1_NA MCC2_NA
#     713     543     935

# exact match (not lenient):
# MCC0_NA MCC1_NA MCC2_NA
#     906     659     649

#------ Grouped vs atomic results ------#

accuracy_with_fallback_grouped %>%
  summarise(
    OA_NA = sum(is.na(Overall_Accuracy)),
    MCC_NA = sum(is.na(MCC_Multiclass))
  )

# without ACI2025:
# OA_NA MCC_NA
#     0    361

# with ACI2025:
# OA_NA MCC_NA
#     0    466

# exact match (not lenient):
# OA_NA MCC_NA
#     0    599

accuracy_with_fallback_atomic %>%
  summarise(
    OA_NA = sum(is.na(Overall_Accuracy)),
    MCC_NA = sum(is.na(MCC_Multiclass))
  )

# without ACI2025:
# OA_NA MCC_NA
#     0    374

# with ACI2025:
# OA_NA MCC_NA
#     0    466

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

# without ACI2025:
#   Truth_Classes_Present Pred_Classes_Present     n
# 1                     1                    1    18
# 2                     1                    2    60
# 3                     1                    3   268
# 4                     2                    1    15
# 5                     2                    2   204
# 6                     2                    3   509
# 7                     3                    2    51
# 8                     3                    3   480

# with ACI2025:
#   Truth_Classes_Present Pred_Classes_Present     n
# 1                     1                    1    23
# 2                     1                    2   141
# 3                     1                    3   287
# 4                     2                    1    15
# 5                     2                    2   165
# 6                     2                    3   512
# 7                     3                    2    51
# 8                     3                    3   499

# exact match (not lenient):
#   Truth_Classes_Present Pred_Classes_Present     n
# 1                     1                    1    23
# 2                     1                    2   114
# 3                     1                    3   450
# 4                     2                    1    12
# 5                     2                    2    82
# 6                     2                    3   249
# 7                     3                    2     9
# 8                     3                    3   261

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

# without ACI2025:
#   Truth_Classes_Present Pred_Classes_Present     n
# 1                     1                    1    19
# 2                     1                    2    60
# 3                     1                    3   280
# 4                     2                    1    15
# 5                     2                    2   204
# 6                     2                    3   497
# 7                     3                    2    51
# 8                     3                    3   481

# with ACI2025:
#   Truth_Classes_Present Pred_Classes_Present     n
# 1                     1                    1    23
# 2                     1                    2   141
# 3                     1                    3   287
# 4                     2                    1    15
# 5                     2                    2   165
# 6                     2                    3   512
# 7                     3                    2    51
# 8                     3                    3   499

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

# without ACI2025:
#                   Overall_Accuracy Balanced_Accuracy Weighted_F1  Macro_F1
# Overall_Accuracy         1.0000000         0.8568095   0.6588147 0.8817312
# Balanced_Accuracy        0.8568095         1.0000000   0.4776229 0.7680445
# Weighted_F1              0.6588147         0.4776229   1.0000000 0.4552734
# Macro_F1                 0.8817312         0.7680445   0.4552734 1.0000000

# with ACI2025:
#                   Overall_Accuracy Balanced_Accuracy Weighted_F1  Macro_F1
# Overall_Accuracy         1.0000000         0.8710355   0.6313420 0.8748702
# Balanced_Accuracy        0.8710355         1.0000000   0.4460394 0.7913268
# Weighted_F1              0.6313420         0.4460394   1.0000000 0.4163951
# Macro_F1                 0.8748702         0.7913268   0.4163951 1.0000000

# exact match (not lenient):
#                   Overall_Accuracy Balanced_Accuracy Weighted_F1  Macro_F1
# Overall_Accuracy         1.0000000         0.9135399   0.2617869 0.9041072
# Balanced_Accuracy        0.9135399         1.0000000   0.2124356 0.8202068
# Weighted_F1              0.2617869         0.2124356   1.0000000 0.1018033
# Macro_F1                 0.9041072         0.8202068   0.1018033 1.0000000

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
    x = "n evaluated",
    y = "Count of ecodistricts"
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
    x = "n evaluated (log10 scale)",
    y = "Count of ecodistricts"
  ) +
  theme_bw()

#------ Number of rows with n less than 20 ------#

accuracy_with_fallback_grouped %>%
  summarise(
    n_below_20 = sum(N_evaluated < 20),
    pct_below_20 = 100 * mean(N_evaluated < 20)
  )

# without ACI2025:
# n_below_20 pct_below_20
#         10        0.623

# with ACI2025:
# n_below_20 pct_below_20
#         25         1.48

# exact match (not lenient):
# n_below_20 pct_below_20
#         31         2.58

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
  print(width = Inf, n = Inf)

# without ACI2025:
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

# with ACI2025:
#   raster  Ecds_ID Ecrgn_n                     Ecprv source_level N_evaluated Truth_Classes_Present
# 1 aci2024     374 Lac Seul Upland               6.2 ecoprovince            1                     1
# 2 aci2025     960 Cascade Ranges               13.2 ecoprovince            1                     1
# 3 ngi         623 Western Alberta Upland        9.1 ecoprovince            1                     1
# 4 ngi         850 Lake Manitoba Plain          10.1 ecoprovince            1                     1
# 5 cnvc        960 Cascade Ranges               13.2 ecoprovince            2                     1
# 6 aci2025     591 Peace Lowland                 9.2 ecoprovince            3                     1
# 7 aci2025     617 Mid-Boreal Uplands            9.2 ecoprovince            3                     1
# 8 aci2025     982 Fraser Basin                 14.1 ecoprovince            4                     1
# 9 aci2025     727 Aspen Parkland               10.2 ecoprovince            6                     1
# 10 aci2025     739 Aspen Parkland               10.2 ecoprovince            6                     1
# 11 aci2025     750 Aspen Parkland               10.2 ecoprovince            6                     1
# 12 aci2025     757 Aspen Parkland               10.2 ecoprovince            6                     1
# 13 aci2025     955 Eastern Vancouver Island     13.1 ecoprovince            6                     1
# 14 aci2025     956 Eastern Vancouver Island     13.1 ecoprovince            6                     1
# 15 aci2025    1017 Northern Continental Divide  14.4 ecoprovince            7                     1
# 16 aci2025    1018 Northern Continental Divide  14.4 ecoprovince            7                     1
# 17 cgr         375 Lake of the Woods             6.5 ecoprovince           11                     2
# 18 cgr         687 Boreal Transition             9.2 ecoprovince           11                     1
# 19 cgr         688 Boreal Transition             9.2 ecoprovince           11                     1
# 20 cgr         695 Boreal Transition             9.2 ecoprovince           11                     1
# 21 cgr         701 Boreal Transition             9.2 ecoprovince           11                     1
# 22 cgr         709 Boreal Transition             9.2 ecoprovince           11                     1
# 23 pudl        375 Lake of the Woods             6.5 ecoprovince           11                     1
# 24 aci2025     618 Western Alberta Upland        9.1 ecoprovince           13                     2
# 25 aci2025     631 Western Alberta Upland        9.1 ecoprovince           13                     2

#   Pred_Classes_Present Macro_F1 Balanced_Accuracy Weighted_F1 Overall_Accuracy Support_0 Support_1
# 1                    1  NaN                 0          0                 0             1         0
# 2                    1  NaN                 0          0                 0             0         0
# 3                    1  NaN                 0          0                 0             0         0
# 4                    1  NaN                 0          0                 0             0         0
# 5                    1  NaN                 0          0                 0             0         0
# 6                    2    0.8               0.667      0.8               0.667         3         0
# 7                    2    0.8               0.667      0.8               0.667         3         0
# 8                    2    0.4               0.25       0.4               0.25          4         0
# 9                    2    0.8               0.667      0                 0.667         0         0
# 10                    2    0.8               0.667      0                 0.667         0         0
# 11                    2    0.8               0.667      0                 0.667         0         0
# 12                    2    0.8               0.667      0                 0.667         0         0
# 13                    1    1                 1          1                 1             6         0
# 14                    1    1                 1          1                 1             6         0
# 15                    2    0.444             0.286      0                 0.286         0         0
# 16                    2    0.444             0.286      0                 0.286         0         0
# 17                    2  NaN                 0          0                 0             0         4
# 18                    1  NaN                 0          0                 0             0         0
# 19                    1  NaN                 0          0                 0             0         0
# 20                    1  NaN                 0          0                 0             0         0
# 21                    1  NaN                 0          0                 0             0         0
# 22                    1  NaN                 0          0                 0             0         0
# 23                    3    0.429             0.273      0                 0.273         0         0
# 24                    2    0.405             0.708      0.0623            0.462         1         0
# 25                    2    0.405             0.708      0.0623            0.462         1         0

#   Support_2
# 1         0
# 2         1
# 3         1
# 4         1
# 5         2
# 6         0
# 7         0
# 8         0
# 9         6
# 10         6
# 11         6
# 12         6
# 13         0
# 14         0
# 15         7
# 16         7
# 17         7
# 18        11
# 19        11
# 20        11
# 21        11
# 22        11
# 23        11
# 24        12
# 25        12

# exact match (not lenient):
#   raster  Ecds_ID Ecrgn_n                     Ecprv source_level N_evaluated Truth_Classes_Present
# 1 aci2024     850 Lake Manitoba Plain          10.1 ecoprovince            1                     1
# 2 aci2025     960 Cascade Ranges               13.2 ecoprovince            1                     1
# 3 gpw         850 Lake Manitoba Plain          10.1 ecoprovince            1                     1
# 4 gpw         960 Cascade Ranges               13.2 ecoprovince            1                     1
# 5 ngi         623 Western Alberta Upland        9.1 ecoprovince            1                     1
# 6 ngi         850 Lake Manitoba Plain          10.1 ecoprovince            1                     1
# 7 aci2024     960 Cascade Ranges               13.2 ecoprovince            2                     1
# 8 aci2024     982 Fraser Basin                 14.1 ecoprovince            2                     1
# 9 cnvc        960 Cascade Ranges               13.2 ecoprovince            2                     1
# 10 gpw         982 Fraser Basin                 14.1 ecoprovince            2                     1
# 11 aci2025     591 Peace Lowland                 9.2 ecoprovince            3                     1
# 12 aci2025     617 Mid-Boreal Uplands            9.2 ecoprovince            3                     1
# 13 aci2025     982 Fraser Basin                 14.1 ecoprovince            4                     1
# 14 aci2025     727 Aspen Parkland               10.2 ecoprovince            6                     1
# 15 aci2025     739 Aspen Parkland               10.2 ecoprovince            6                     1
# 16 aci2025     750 Aspen Parkland               10.2 ecoprovince            6                     1
# 17 aci2025     757 Aspen Parkland               10.2 ecoprovince            6                     1
# 18 aci2025     955 Eastern Vancouver Island     13.1 ecoprovince            6                     1
# 19 aci2025     956 Eastern Vancouver Island     13.1 ecoprovince            6                     1
# 20 aci2025    1017 Northern Continental Divide  14.4 ecoprovince            7                     1
# 21 aci2025    1018 Northern Continental Divide  14.4 ecoprovince            7                     1
# 22 cgr         375 Lake of the Woods             6.5 ecoprovince           11                     2
# 23 cgr         687 Boreal Transition             9.2 ecoprovince           11                     1
# 24 cgr         688 Boreal Transition             9.2 ecoprovince           11                     1
# 25 cgr         695 Boreal Transition             9.2 ecoprovince           11                     1
# 26 cgr         701 Boreal Transition             9.2 ecoprovince           11                     1
# 27 cgr         709 Boreal Transition             9.2 ecoprovince           11                     1
# 28 pudl        375 Lake of the Woods             6.5 ecoprovince           11                     1
# 29 aci2025     618 Western Alberta Upland        9.1 ecoprovince           13                     2
# 30 aci2025     631 Western Alberta Upland        9.1 ecoprovince           13                     2
# 31 gpw         631 Western Alberta Upland        9.1 ecoprovince           17                     1

#   Pred_Classes_Present Macro_F1 Balanced_Accuracy Weighted_F1 Overall_Accuracy Support_0 Support_1
# 1                    1  NaN                 0          0                 0             0         0
# 2                    1  NaN                 0          0                 0             0         0
# 3                    1  NaN                 0          0                 0             0         0
# 4                    1  NaN                 0          0                 0             0         0
# 5                    1  NaN                 0          0                 0             0         0
# 6                    1  NaN                 0          0                 0             0         0
# 7                    1  NaN                 0          0                 0             0         0
# 8                    1  NaN                 0          0                 0             0         0
# 9                    1  NaN                 0          0                 0             0         0
# 10                    1  NaN                 0          0                 0             0         0
# 11                    2    0.8               0.667      0.8               0.667         3         0
# 12                    2    0.8               0.667      0.8               0.667         3         0
# 13                    2    0.4               0.25       0.4               0.25          4         0
# 14                    2    0.8               0.667      0                 0.667         0         0
# 15                    2    0.8               0.667      0                 0.667         0         0
# 16                    2    0.8               0.667      0                 0.667         0         0
# 17                    2    0.8               0.667      0                 0.667         0         0
# 18                    1    1                 1          1                 1             6         0
# 19                    1    1                 1          1                 1             6         0
# 20                    2    0.444             0.286      0                 0.286         0         0
# 21                    2    0.444             0.286      0                 0.286         0         0
# 22                    2  NaN                 0          0                 0             0         4
# 23                    1  NaN                 0          0                 0             0         0
# 24                    1  NaN                 0          0                 0             0         0
# 25                    1  NaN                 0          0                 0             0         0
# 26                    1  NaN                 0          0                 0             0         0
# 27                    1  NaN                 0          0                 0             0         0
# 28                    3    0.429             0.273      0                 0.273         0         0
# 29                    2    0.405             0.708      0.0623            0.462         1         0
# 30                    2    0.405             0.708      0.0623            0.462         1         0
# 31                    2  NaN                 0          0                 0             0         0

#   Support_2
# 1         1
# 2         1
# 3         1
# 4         1
# 5         1
# 6         1
# 7         2
# 8         2
# 9         2
# 10         2
# 11         0
# 12         0
# 13         0
# 14         6
# 15         6
# 16         6
# 17         6
# 18         0
# 19         0
# 20         7
# 21         7
# 22         7
# 23        11
# 24        11
# 25        11
# 26        11
# 27        11
# 28        11
# 29        12
# 30        12
# 31        17

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

# without ACI2025:
#   raster Ecds_ID Ecrgn_n                Ecprv source_level N_evaluated Truth_Classes_Present
# 1 ngi        623 Western Alberta Upland   9.1 ecoprovince            1                     1
# 2 ngi        850 Lake Manitoba Plain     10.1 ecoprovince            1                     1
# 3 cnvc       960 Cascade Ranges          13.2 ecoprovince            2                     1

#   Pred_Classes_Present Macro_F1 Balanced_Accuracy Weighted_F1 Overall_Accuracy
# 1                    1      NaN                 0           0                0
# 2                    1      NaN                 0           0                0
# 3                    1      NaN                 0           0                0

# with ACI2025:
#   raster  Ecds_ID Ecrgn_n                Ecprv source_level N_evaluated Truth_Classes_Present
# 1 aci2024     374 Lac Seul Upland          6.2 ecoprovince            1                     1
# 2 aci2025     960 Cascade Ranges          13.2 ecoprovince            1                     1
# 3 ngi         623 Western Alberta Upland   9.1 ecoprovince            1                     1
# 4 ngi         850 Lake Manitoba Plain     10.1 ecoprovince            1                     1
# 5 cnvc        960 Cascade Ranges          13.2 ecoprovince            2                     1

#   Pred_Classes_Present Macro_F1 Balanced_Accuracy Weighted_F1 Overall_Accuracy
# 1                    1      NaN                 0           0                0
# 2                    1      NaN                 0           0                0
# 3                    1      NaN                 0           0                0
# 4                    1      NaN                 0           0                0
# 5                    1      NaN                 0           0                0

# exact match (not lenient):
#   raster  Ecds_ID Ecrgn_n                Ecprv source_level N_evaluated Truth_Classes_Present
# 1 aci2024     850 Lake Manitoba Plain     10.1 ecoprovince            1                     1
# 2 aci2025     960 Cascade Ranges          13.2 ecoprovince            1                     1
# 3 gpw         850 Lake Manitoba Plain     10.1 ecoprovince            1                     1
# 4 gpw         960 Cascade Ranges          13.2 ecoprovince            1                     1
# 5 ngi         623 Western Alberta Upland   9.1 ecoprovince            1                     1
# 6 ngi         850 Lake Manitoba Plain     10.1 ecoprovince            1                     1
# 7 aci2024     960 Cascade Ranges          13.2 ecoprovince            2                     1
# 8 aci2024     982 Fraser Basin            14.1 ecoprovince            2                     1
# 9 cnvc        960 Cascade Ranges          13.2 ecoprovince            2                     1
# 10 gpw         982 Fraser Basin            14.1 ecoprovince            2                     1

#   Pred_Classes_Present Macro_F1 Balanced_Accuracy Weighted_F1 Overall_Accuracy
# 1                    1      NaN                 0           0                0
# 2                    1      NaN                 0           0                0
# 3                    1      NaN                 0           0                0
# 4                    1      NaN                 0           0                0
# 5                    1      NaN                 0           0                0
# 6                    1      NaN                 0           0                0
# 7                    1      NaN                 0           0                0
# 8                    1      NaN                 0           0                0
# 9                    1      NaN                 0           0                0
# 10                    1      NaN                 0           0                0

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

#---- Export ----

# GROUPED
readr::write_csv(
  accuracy_with_fallback_grouped_cons,
  file.path(out_dir, "nattam_accuracy_grouped_combinedCIandCM_forconsensus_29Jul2026.csv")
)

# ATOMIC
readr::write_csv(
  accuracy_with_fallback_atomic_cons,
  file.path(out_dir, "nattam_accuracy_atomic_combinedCIandCM_forconsensus_29Jul2026.csv")
)

#------ How often is fallback used? ------#

accuracy_with_fallback_grouped_cons %>%
  count(acc_source)

# without ACI2025:
#   acc_source            n
# 1 Balanced_Accuracy    51
# 2 Macro_F1           1554

# with ACI2025:
#   acc_source            n
# 1 Balanced_Accuracy    52
# 2 Macro_F1           1641

# exact match (not lenient):
#   acc_source            n
# 1 Balanced_Accuracy    58
# 2 Macro_F1           1142

accuracy_with_fallback_grouped_cons %>%
  filter(acc_source == "Balanced_Accuracy") %>%
  summarise(
    n = n(),
    pct = 100 * n() / nrow(accuracy_with_fallback_grouped_cons)
  )

# without ACI2025:
#  n   pct
# 51  3.18

# with ACI2025:
# n   pct
# 52  3.07

# exact match (not lenient):
#  n   pct
# 58  4.83

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

# without ACI2025:
# n_rows one_truth_class two_truth_classes three_truth_classes one_pred_class two_pred_classes
#     51              28                23                   0             25               26
# three_pred_classes
#                  0

# with ACI2025:
# n_rows one_truth_class two_truth_classes three_truth_classes one_pred_class two_pred_classes
#     52              29                23                   0             27               25
# three_pred_classes
#                  0

# exact match (not lenient):
# n_rows one_truth_class two_truth_classes three_truth_classes one_pred_class two_pred_classes
#     58              34                24                   0             30               28
# three_pred_classes
#                  0

#------ How many fallback rows have extreme BA values? ------#

accuracy_with_fallback_grouped_cons %>%
  filter(
    acc_source == "Balanced_Accuracy",
    acc %in% c(0, 1)
  ) %>%
  count(acc)

# without ACI2025:
#   acc     n
#     0    51

# with ACI2025:
#  acc     n
#    0    52

# exact match (not lenient):
#  acc     n
#    0    58

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

# without ACI2025:
#   n min_n median_n max_n
#  51     1       32   239

# with ACI2025:
#   n min_n median_n max_n
#  52     1       32   239

# exact match (not lenient):
#   n min_n median_n max_n
#  58     1     31.5   239

#------ Distribution of n among fallback rows ------#

accuracy_with_fallback_grouped_cons %>%
  filter(
    acc_source == "Balanced_Accuracy",
    acc %in% c(0, 1)
  ) %>%
  count(N_evaluated, sort = TRUE) %>%
  print(n = Inf)

# without ACI2025:
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

# with ACI2025:
#   N_evaluated     n
# 1          40     8
# 2          11     6
# 3           1     4
# 4          27     4
# 5          32     4
# 6          45     4
# 7          21     3
# 8          20     2
# 9          34     2
# 10          43     2
# 11           2     1
# 12          22     1
# 13          24     1
# 14          30     1
# 15          31     1
# 16          35     1
# 17          38     1
# 18          41     1
# 19          51     1
# 20          60     1
# 21          74     1
# 22         202     1
# 23         239     1

# exact match (not lenient):
#   N_evaluated     n
# 1          40     7
# 2           1     6
# 3          11     6
# 4          32     5
# 5           2     4
# 6          27     4
# 7          45     4
# 8          21     3
# 9          20     2
# 10          34     2
# 11          43     2
# 12          17     1
# 13          22     1
# 14          24     1
# 15          31     1
# 16          35     1
# 17          38     1
# 18          41     1
# 19          44     1
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

# without ACI2025:
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


# with ACI2025:
#   raster  Ecds_ID Ecrgn_n                     Ecprv source_level N_evaluated Truth_Classes_Present
# 1 ngi         731 Aspen Parkland               10.2 ecodistrict          239                     1
# 2 ngi         740 Aspen Parkland               10.2 ecodistrict          202                     1
# 3 cnvc       1017 Northern Continental Divide  14.4 ecodistrict           74                     2
# 4 pudl        678 Boreal Transition             9.2 ecodistrict           60                     2
# 5 pudl        727 Aspen Parkland               10.2 ecodistrict           51                     2
# 6 ngi         683 Boreal Transition             9.2 ecoregion             45                     2
# 7 ngi         692 Boreal Transition             9.2 ecoregion             45                     2
# 8 ngi         701 Boreal Transition             9.2 ecoregion             45                     2
# 9 ngi         708 Boreal Transition             9.2 ecoregion             45                     2
# 10 pudl        684 Boreal Transition             9.2 ecodistrict           43                     2
# 11 pudl        760 Aspen Parkland               10.2 ecodistrict           43                     1
# 12 cnvc        380 Lake of the Woods             6.5 ecoregion             41                     1
# 13 cnvc        622 Western Alberta Upland        9.1 ecoregion             40                     2
# 14 cnvc        623 Western Alberta Upland        9.1 ecoregion             40                     2
# 15 cnvc        628 Western Alberta Upland        9.1 ecoregion             40                     2
# 16 cnvc        629 Western Alberta Upland        9.1 ecoregion             40                     2
# 17 cnvc        630 Western Alberta Upland        9.1 ecoregion             40                     2
# 18 cnvc        631 Western Alberta Upland        9.1 ecoregion             40                     2
# 19 gpw         601 Mid-Boreal Uplands            9.2 ecodistrict           40                     1
# 20 pudl        592 Peace Lowland                 9.2 ecodistrict           40                     2
# 21 eccc       1017 Northern Continental Divide  14.4 ecodistrict           38                     1
# 22 pudl        703 Boreal Transition             9.2 ecodistrict           35                     2
# 23 ngi         727 Aspen Parkland               10.2 ecodistrict           34                     1
# 24 ngi         744 Aspen Parkland               10.2 ecodistrict           34                     1
# 25 cnvc        375 Lake of the Woods             6.5 ecodistrict           32                     1
# 26 pudl        686 Boreal Transition             9.2 ecodistrict           32                     2
# 27 pudl        690 Boreal Transition             9.2 ecodistrict           32                     1
# 28 pudl        728 Aspen Parkland               10.2 ecodistrict           32                     2
# 29 pudl        692 Boreal Transition             9.2 ecodistrict           31                     2
# 30 gpw         628 Western Alberta Upland        9.1 ecodistrict           30                     1
# 31 cnvc        678 Boreal Transition             9.2 ecodistrict           27                     2
# 32 ngi         684 Boreal Transition             9.2 ecodistrict           27                     1
# 33 ngi         746 Aspen Parkland               10.2 ecodistrict           27                     1
# 34 pudl        755 Aspen Parkland               10.2 ecodistrict           27                     1
# 35 cnvc        847 Lake Manitoba Plain          10.1 ecodistrict           24                     1
# 36 cnvc        684 Boreal Transition             9.2 ecodistrict           22                     2
# 37 cnvc        997 Eastern Continental Ranges   14.4 ecoregion             21                     1
# 38 cnvc        999 Eastern Continental Ranges   14.4 ecoregion             21                     1
# 39 pudl        590 Peace Lowland                 9.2 ecodistrict           21                     2
# 40 cnvc        839 Aspen Parkland               10.2 ecodistrict           20                     1
# 41 pudl        839 Aspen Parkland               10.2 ecodistrict           20                     1
# 42 cgr         375 Lake of the Woods             6.5 ecoprovince           11                     2
# 43 cgr         687 Boreal Transition             9.2 ecoprovince           11                     1
# 44 cgr         688 Boreal Transition             9.2 ecoprovince           11                     1
# 45 cgr         695 Boreal Transition             9.2 ecoprovince           11                     1
# 46 cgr         701 Boreal Transition             9.2 ecoprovince           11                     1
# 47 cgr         709 Boreal Transition             9.2 ecoprovince           11                     1
# 48 cnvc        960 Cascade Ranges               13.2 ecoprovince            2                     1
# 49 aci2024     374 Lac Seul Upland               6.2 ecoprovince            1                     1
# 50 aci2025     960 Cascade Ranges               13.2 ecoprovince            1                     1
# 51 ngi         623 Western Alberta Upland        9.1 ecoprovince            1                     1
# 52 ngi         850 Lake Manitoba Plain          10.1 ecoprovince            1                     1

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
# 36                    1     0 Balanced_Accuracy      NaN                 0                0           0
# 37                    2     0 Balanced_Accuracy      NaN                 0                0           0
# 38                    2     0 Balanced_Accuracy      NaN                 0                0           0
# 39                    1     0 Balanced_Accuracy      NaN                 0                0           0
# 40                    2     0 Balanced_Accuracy      NaN                 0                0           0
# 41                    1     0 Balanced_Accuracy      NaN                 0                0           0
# 42                    2     0 Balanced_Accuracy      NaN                 0                0           0
# 43                    1     0 Balanced_Accuracy      NaN                 0                0           0
# 44                    1     0 Balanced_Accuracy      NaN                 0                0           0
# 45                    1     0 Balanced_Accuracy      NaN                 0                0           0
# 46                    1     0 Balanced_Accuracy      NaN                 0                0           0
# 47                    1     0 Balanced_Accuracy      NaN                 0                0           0
# 48                    1     0 Balanced_Accuracy      NaN                 0                0           0
# 49                    1     0 Balanced_Accuracy      NaN                 0                0           0
# 50                    1     0 Balanced_Accuracy      NaN                 0                0           0
# 51                    1     0 Balanced_Accuracy      NaN                 0                0           0
# 52                    1     0 Balanced_Accuracy      NaN                 0                0           0

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
# 35         0         0        24    NA    NA    NA
# 36         0         3        19    NA    NA    NA
# 37         0         0        21    NA    NA    NA
# 38         0         0        21    NA    NA    NA
# 39         0        20         1    NA    NA    NA
# 40         0         0        20    NA    NA    NA
# 41         0         0        20    NA    NA    NA
# 42         0         4         7    NA    NA    NA
# 43         0         0        11    NA    NA    NA
# 44         0         0        11    NA    NA    NA
# 45         0         0        11    NA    NA    NA
# 46         0         0        11    NA    NA    NA
# 47         0         0        11    NA    NA    NA
# 48         0         0         2    NA    NA    NA
# 49         1         0         0    NA    NA    NA
# 50         0         0         1    NA    NA    NA
# 51         0         0         1    NA    NA    NA
# 52         0         0         1    NA    NA    NA


# exact match (not lenient):
#   raster  Ecds_ID Ecrgn_n                     Ecprv source_level N_evaluated Truth_Classes_Present
# 1 ngi         731 Aspen Parkland               10.2 ecodistrict          239                     1
# 2 ngi         740 Aspen Parkland               10.2 ecodistrict          202                     1
# 3 cnvc       1017 Northern Continental Divide  14.4 ecodistrict           74                     2
# 4 pudl        678 Boreal Transition             9.2 ecodistrict           60                     2
# 5 pudl        727 Aspen Parkland               10.2 ecodistrict           51                     2
# 6 ngi         683 Boreal Transition             9.2 ecoregion             45                     2
# 7 ngi         692 Boreal Transition             9.2 ecoregion             45                     2
# 8 ngi         701 Boreal Transition             9.2 ecoregion             45                     2
# 9 ngi         708 Boreal Transition             9.2 ecoregion             45                     2
# 10 eccc        760 Aspen Parkland               10.2 ecodistrict           44                     2
# 11 pudl        684 Boreal Transition             9.2 ecodistrict           43                     2
# 12 pudl        760 Aspen Parkland               10.2 ecodistrict           43                     1
# 13 cnvc        380 Lake of the Woods             6.5 ecoregion             41                     1
# 14 cnvc        622 Western Alberta Upland        9.1 ecoregion             40                     2
# 15 cnvc        623 Western Alberta Upland        9.1 ecoregion             40                     2
# 16 cnvc        628 Western Alberta Upland        9.1 ecoregion             40                     2
# 17 cnvc        629 Western Alberta Upland        9.1 ecoregion             40                     2
# 18 cnvc        630 Western Alberta Upland        9.1 ecoregion             40                     2
# 19 cnvc        631 Western Alberta Upland        9.1 ecoregion             40                     2
# 20 pudl        592 Peace Lowland                 9.2 ecodistrict           40                     2
# 21 eccc       1017 Northern Continental Divide  14.4 ecodistrict           38                     1
# 22 pudl        703 Boreal Transition             9.2 ecodistrict           35                     2
# 23 ngi         727 Aspen Parkland               10.2 ecodistrict           34                     1
# 24 ngi         744 Aspen Parkland               10.2 ecodistrict           34                     1
# 25 cnvc        375 Lake of the Woods             6.5 ecodistrict           32                     1
# 26 eccc        690 Boreal Transition             9.2 ecodistrict           32                     1
# 27 pudl        686 Boreal Transition             9.2 ecodistrict           32                     2
# 28 pudl        690 Boreal Transition             9.2 ecodistrict           32                     1
# 29 pudl        728 Aspen Parkland               10.2 ecodistrict           32                     2
# 30 pudl        692 Boreal Transition             9.2 ecodistrict           31                     2
# 31 cnvc        678 Boreal Transition             9.2 ecodistrict           27                     2
# 32 ngi         684 Boreal Transition             9.2 ecodistrict           27                     1
# 33 ngi         746 Aspen Parkland               10.2 ecodistrict           27                     1
# 34 pudl        755 Aspen Parkland               10.2 ecodistrict           27                     1
# 35 cnvc        847 Lake Manitoba Plain          10.1 ecodistrict           24                     1
# 36 cnvc        684 Boreal Transition             9.2 ecodistrict           22                     2
# 37 cnvc        997 Eastern Continental Ranges   14.4 ecoregion             21                     1
# 38 cnvc        999 Eastern Continental Ranges   14.4 ecoregion             21                     1
# 39 pudl        590 Peace Lowland                 9.2 ecodistrict           21                     2
# 40 cnvc        839 Aspen Parkland               10.2 ecodistrict           20                     1
# 41 pudl        839 Aspen Parkland               10.2 ecodistrict           20                     1
# 42 gpw         631 Western Alberta Upland        9.1 ecoprovince           17                     1
# 43 cgr         375 Lake of the Woods             6.5 ecoprovince           11                     2
# 44 cgr         687 Boreal Transition             9.2 ecoprovince           11                     1
# 45 cgr         688 Boreal Transition             9.2 ecoprovince           11                     1
# 46 cgr         695 Boreal Transition             9.2 ecoprovince           11                     1
# 47 cgr         701 Boreal Transition             9.2 ecoprovince           11                     1
# 48 cgr         709 Boreal Transition             9.2 ecoprovince           11                     1
# 49 aci2024     960 Cascade Ranges               13.2 ecoprovince            2                     1
# 50 aci2024     982 Fraser Basin                 14.1 ecoprovince            2                     1
# 51 cnvc        960 Cascade Ranges               13.2 ecoprovince            2                     1
# 52 gpw         982 Fraser Basin                 14.1 ecoprovince            2                     1
# 53 aci2024     850 Lake Manitoba Plain          10.1 ecoprovince            1                     1
# 54 aci2025     960 Cascade Ranges               13.2 ecoprovince            1                     1
# 55 gpw         850 Lake Manitoba Plain          10.1 ecoprovince            1                     1
# 56 gpw         960 Cascade Ranges               13.2 ecoprovince            1                     1
# 57 ngi         623 Western Alberta Upland        9.1 ecoprovince            1                     1
# 58 ngi         850 Lake Manitoba Plain          10.1 ecoprovince            1                     1

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
# 10                    2     0 Balanced_Accuracy      NaN                 0                0           0
# 11                    1     0 Balanced_Accuracy      NaN                 0                0           0
# 12                    2     0 Balanced_Accuracy      NaN                 0                0           0
# 13                    2     0 Balanced_Accuracy      NaN                 0                0           0
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
# 29                    2     0 Balanced_Accuracy      NaN                 0                0           0
# 30                    1     0 Balanced_Accuracy      NaN                 0                0           0
# 31                    1     0 Balanced_Accuracy      NaN                 0                0           0
# 32                    2     0 Balanced_Accuracy      NaN                 0                0           0
# 33                    2     0 Balanced_Accuracy      NaN                 0                0           0
# 34                    1     0 Balanced_Accuracy      NaN                 0                0           0
# 35                    2     0 Balanced_Accuracy      NaN                 0                0           0
# 36                    1     0 Balanced_Accuracy      NaN                 0                0           0
# 37                    2     0 Balanced_Accuracy      NaN                 0                0           0
# 38                    2     0 Balanced_Accuracy      NaN                 0                0           0
# 39                    1     0 Balanced_Accuracy      NaN                 0                0           0
# 40                    2     0 Balanced_Accuracy      NaN                 0                0           0
# 41                    1     0 Balanced_Accuracy      NaN                 0                0           0
# 42                    2     0 Balanced_Accuracy      NaN                 0                0           0
# 43                    2     0 Balanced_Accuracy      NaN                 0                0           0
# 44                    1     0 Balanced_Accuracy      NaN                 0                0           0
# 45                    1     0 Balanced_Accuracy      NaN                 0                0           0
# 46                    1     0 Balanced_Accuracy      NaN                 0                0           0
# 47                    1     0 Balanced_Accuracy      NaN                 0                0           0
# 48                    1     0 Balanced_Accuracy      NaN                 0                0           0
# 49                    1     0 Balanced_Accuracy      NaN                 0                0           0
# 50                    1     0 Balanced_Accuracy      NaN                 0                0           0
# 51                    1     0 Balanced_Accuracy      NaN                 0                0           0
# 52                    1     0 Balanced_Accuracy      NaN                 0                0           0
# 53                    1     0 Balanced_Accuracy      NaN                 0                0           0
# 54                    1     0 Balanced_Accuracy      NaN                 0                0           0
# 55                    1     0 Balanced_Accuracy      NaN                 0                0           0
# 56                    1     0 Balanced_Accuracy      NaN                 0                0           0
# 57                    1     0 Balanced_Accuracy      NaN                 0                0           0
# 58                    1     0 Balanced_Accuracy      NaN                 0                0           0

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
# 10         0         1        43    NA    NA    NA
# 11         0        14        29    NA    NA    NA
# 12         0         0        43    NA    NA    NA
# 13         0         0        41    NA    NA    NA
# 14         0         3        37    NA    NA    NA
# 15         0         3        37    NA    NA    NA
# 16         0         3        37    NA    NA    NA
# 17         0         3        37    NA    NA    NA
# 18         0         3        37    NA    NA    NA
# 19         0         3        37    NA    NA    NA
# 20         0        36         4    NA    NA    NA
# 21         0         0        38    NA    NA    NA
# 22         0         8        27    NA    NA    NA
# 23         0         0        34    NA    NA    NA
# 24         0         0        34    NA    NA    NA
# 25         0         0        32    NA    NA    NA
# 26         0         0        32    NA    NA    NA
# 27         0         1        31    NA    NA    NA
# 28         0         0        32    NA    NA    NA
# 29         0         3        29    NA    NA    NA
# 30         0        16        15    NA    NA    NA
# 31         0        24         3    NA    NA    NA
# 32         0         0        27    NA    NA    NA
# 33         0         0        27    NA    NA    NA
# 34         0         0        27    NA    NA    NA
# 35         0         0        24    NA    NA    NA
# 36         0         3        19    NA    NA    NA
# 37         0         0        21    NA    NA    NA
# 38         0         0        21    NA    NA    NA
# 39         0        20         1    NA    NA    NA
# 40         0         0        20    NA    NA    NA
# 41         0         0        20    NA    NA    NA
# 42         0         0        17    NA    NA    NA
# 43         0         4         7    NA    NA    NA
# 44         0         0        11    NA    NA    NA
# 45         0         0        11    NA    NA    NA
# 46         0         0        11    NA    NA    NA
# 47         0         0        11    NA    NA    NA
# 48         0         0        11    NA    NA    NA
# 49         0         0         2    NA    NA    NA
# 50         0         0         2    NA    NA    NA
# 51         0         0         2    NA    NA    NA
# 52         0         0         2    NA    NA    NA
# 53         0         0         1    NA    NA    NA
# 54         0         0         1    NA    NA    NA
# 55         0         0         1    NA    NA    NA
# 56         0         0         1    NA    NA    NA
# 57         0         0         1    NA    NA    NA
# 58         0         0         1    NA    NA    NA


accuracy_with_fallback_grouped_cons %>%
  filter(acc_source == "Balanced_Accuracy") %>%
  summarise(
    no_class0 = sum(Support_0 == 0),
    no_class1 = sum(Support_1 == 0),
    no_class2 = sum(Support_2 == 0)
  )

# without ACI2025:
# no_class0 no_class1 no_class2
#        46        30         3

# with ACI2025:
# no_class0 no_class1 no_class2
#        46        32         3

# exact match (not lenient):
# no_class0 no_class1 no_class2
#        54        38         0

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

# without ACI2025:
#   class metric        n  mean median    sd    q25   q75
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

# with ACI2025:
#   class metric        n   NAs  mean median    sd    q25   q75
# 1 0     F1         1039   654 0.720  0.836 0.283 0.560  0.944
# 2 0     MCC_OvR     980   713 0.496  0.550 0.271 0.319  0.704
# 3 0     Precision  1682    11 0.435  0.395 0.423 0      0.913
# 4 0     Recall     1052   641 0.842  0.943 0.228 0.812  0.984
# 5 1     F1         1120   573 0.498  0.534 0.274 0.269  0.733
# 6 1     MCC_OvR    1150   543 0.343  0.353 0.277 0.101  0.574
# 7 1     Precision  1575   118 0.455  0.509 0.388 0      0.830
# 8 1     Recall     1263   430 0.446  0.455 0.304 0.166  0.718
# 9 2     F1         1004   689 0.506  0.527 0.293 0.243  0.755
# 10 2     MCC_OvR     758   935 0.217  0.139 0.239 0.0382 0.397
# 11 2     Precision  1389   304 0.585  0.833 0.443 0      1    
# 12 2     Recall     1170   523 0.417  0.394 0.326 0.105  0.708
# 13 NA    NA        51147 16573 0.470  0.503 0.350 0.0910 0.769

# exact match (not lenient):
#   class metric        n   NAs   mean median    sd     q25   q75
# 1 0     F1          287   913 0.469  0.48   0.319  0.156  0.775
# 2 0     MCC_OvR     294   906 0.347  0.325  0.321  0.0367 0.639
# 3 0     Precision  1184    16 0.0987 0      0.235  0      0    
# 4 0     Recall      300   900 0.732  0.923  0.331  0.488  0.993
# 5 1     F1          496   704 0.369  0.376  0.261  0.111  0.598
# 6 1     MCC_OvR     541   659 0.163  0.0914 0.250 -0.0314 0.353
# 7 1     Precision  1068   132 0.213  0      0.329  0      0.352
# 8 1     Recall      638   562 0.365  0.326  0.321  0.0252 0.663
# 9 2     F1         1040   160 0.565  0.605  0.284  0.345  0.821
# 10 2     MCC_OvR     551   649 0.173  0.114  0.212  0.0396 0.315
# 11 2     Precision  1073   127 0.890  1      0.244  0.941  1    
# 12 2     Recall     1145    55 0.431  0.417  0.307  0.138  0.708
# 13 NA    NA        31161 16839 0.363  0.25   0.359  0.0255 0.694

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
  labs(
    x = "Class",
    y = "Recall",
    title = "Distribution of recall by class"
  ) +
  theme_bw()

class_metrics_long %>%
  filter(metric == "Precision") %>%
  ggplot(aes(class, value)) +
  geom_boxplot() +
  labs(
    x = "Class",
    y = "Precision",
    title = "Distribution of precision by class"
  ) +
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

# without ACI2025:
# median_support0 median_support1 median_support2 mean_support0 mean_support1 mean_support2
#              33              47              33          330.          123.          351.

# with ACI2025:
# median_support0 median_support1 median_support2 mean_support0 mean_support1 mean_support2
#              26              46              32          317.          114.          334.

# exact match (not lenient):
# median_support0 median_support1 median_support2 mean_support0 mean_support1 mean_support2
#               0               2            122.          23.1          63.3          543.

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

# without ACI2025:
# no_class0 no_class1 no_class2
#        46        30         3

# with ACI2025:
# no_class0 no_class1 no_class2
#        46        32         3

# exact match (not lenient):
# no_class0 no_class1 no_class2
#        54        38         0

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

# without ACI2025:
#   missing_pattern     n
# 1 01                 25
# 2 0                  19
# 3 1                   4
# 4 02                  2
# 5 12                  1

# with ACI2025:
#   missing_pattern     n
# 1 01                 26
# 2 0                  19
# 3 1                   4
# 4 12                  2
# 5 02                  1

# exact match (not lenient):
#   missing_pattern     n
# 1 01                 34
# 2 0                  20
# 3 1                   4

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

# without ACI2025:
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

# with ACI2025:
#   raster  mean_F1_native median_F1_native n_valid
# 1 gpw              0.197            0.133      46
# 2 aci2024          0.224            0.136      51
# 3 plvi             0.316            0.347      13
# 4 mgi              0.439            0.426      32
# 5 cnvc             0.482            0.536     175
# 6 pudl             0.498            0.499     148
# 7 pli              0.505            0.495      67
# 8 eccc             0.555            0.561     131
# 9 cgr              0.578            0.660     118
# 10 pglc             0.606            0.621      83
# 11 ngi              0.607            0.702     101
# 12 aci2025          0.646            0.741      39

# exact match (not lenient):
#   raster  mean_F1_native median_F1_native n_valid
# 1 plvi             0.316            0.347      13
# 2 mgi              0.439            0.426      32
# 3 cnvc             0.482            0.536     175
# 4 pudl             0.498            0.499     148
# 5 pli              0.560            0.5        67
# 6 cgr              0.578            0.660     118
# 7 eccc             0.584            0.588     143
# 8 pglc             0.606            0.621      83
# 9 ngi              0.607            0.702     101
# 10 aci2025          0.646            0.741      39
# 11 aci2024          0.680            0.730      63
# 12 gpw              0.740            0.738      58

# Average tame-class F1 by source
accuracy_with_fallback_grouped %>%
  group_by(raster) %>%
  summarise(
    mean_F1_tame = mean(F1_1, na.rm = TRUE),
    median_F1_tame = median(F1_1, na.rm = TRUE),
    n_valid = sum(!is.na(F1_1))
  ) %>%
  arrange(mean_F1_tame)

# without ACI2025:
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

# with ACI2025:
#   raster  mean_F1_tame median_F1_tame n_valid
# 1 pudl          0.0547         0.0401      50
# 2 ngi           0.0952         0.0547      30
# 3 cnvc          0.119          0.127       42
# 4 cgr           0.154          0.166       54
# 5 plvi          0.224          0.207        7
# 6 mgi           0.397          0.377       39
# 7 pli           0.411          0.412       45
# 8 gpw           0.427          0.412      272
# 9 aci2024       0.576          0.640      225
# 10 pglc          0.679          0.680       86
# 11 aci2025       0.722          0.807       76
# 12 eccc          0.742          0.770      194

# exact match (not lenient):
#   raster  mean_F1_tame median_F1_tame n_valid
# 1 pudl          0.0547         0.0401      50
# 2 ngi           0.0952         0.0547      30
# 3 cnvc          0.119          0.127       42
# 4 cgr           0.154          0.166       54
# 5 plvi          0.224          0.207        7
# 6 eccc          0.376          0.380       92
# 7 mgi           0.397          0.377       39
# 8 aci2025       0.465          0.508       24
# 9 pli           0.560          0.557       20
# 10 gpw           0.611          0.629       52
# 11 pglc          0.679          0.680       86
# 12 aci2024     NaN             NA            0

# Identify systematic omission vs commission errors
# Recall < Precision = class is being omitted
# Precision < Recall = class is being over-predicted
accuracy_with_fallback_grouped %>%
  summarise(
    native_precision = mean(Precision_2, na.rm = TRUE),
    native_recall = mean(Recall_2, na.rm = TRUE),
    
    tame_precision = mean(Precision_1, na.rm = TRUE),
    tame_recall = mean(Recall_1, na.rm = TRUE),
    
    not_precision = mean(Precision_0, na.rm = TRUE),
    not_recall = mean(Recall_0, na.rm = TRUE)
  )

# without ACI2025:
# native_precision native_recall tame_precision tame_recall not_precision not_recall
#            0.588         0.409          0.458       0.454         0.446      0.860

# with ACI2025:
# native_precision native_recall tame_precision tame_recall not_precision not_recall
#            0.585         0.417          0.455       0.446         0.435      0.842

# exact match (not lenient):
# native_precision native_recall tame_precision tame_recall not_precision not_recall
#            0.890         0.431          0.213       0.365        0.0987      0.732

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

# with ACI2025:
# mean_F1_0 mean_F1_1 mean_F1_2 median_F1_0 median_F1_1 median_F1_2
#     0.720     0.498     0.506       0.836       0.534       0.527
# Non-grassland consistently the easiest class to classify (mean F1 = 0.72, median F1 = 0.84). 
# In contrast, tame and native grassland show substantially lower and nearly identical performance 
# (mean F1 ≈ 0.50; median F1 ≈ 0.53). These results suggest that the primary classification challenge 
# is distinguishing between grassland subtypes.

# exact match (not lenient):
# mean_F1_0 mean_F1_1 mean_F1_2 median_F1_0 median_F1_1 median_F1_2
#     0.469     0.369     0.565        0.48       0.376       0.605

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

# without ACI2025:
# cor_native cor_tame cor_not
#      0.108   -0.202  0.0477

# with ACI2025:
# cor_native cor_tame cor_not
#     0.0844   -0.212  0.0827 - sample size not major driver of F1 values

# exact match (not lenient):
# cor_native cor_tame cor_not
#      0.112   -0.348  -0.373

#---- Class-specific metrics (class-specific F1 or fallbacks if NA) ----
accuracy_with_fallback_grouped_cons2 <- accuracy_with_fallback_grouped %>%
  mutate(
    native_metric = coalesce(F1_2, Macro_F1, Balanced_Accuracy),
    tame_metric   = coalesce(F1_1, Macro_F1, Balanced_Accuracy),
    not_metric    = coalesce(F1_0, Macro_F1, Balanced_Accuracy))

accuracy_with_fallback_atomic_cons2 <- accuracy_with_fallback_atomic %>%
  mutate(
    native_metric = coalesce(F1_2, Macro_F1, Balanced_Accuracy),
    tame_metric   = coalesce(F1_1, Macro_F1, Balanced_Accuracy),
    not_metric    = coalesce(F1_0, Macro_F1, Balanced_Accuracy))

#---- Export ----

# GROUPED
readr::write_csv(
  accuracy_with_fallback_grouped_cons2,
  file.path(out_dir, "nattam_accuracy_grouped_combinedCIandCM_forconsensus_29Jul2026_final.csv")
)

# ATOMIC
readr::write_csv(
  accuracy_with_fallback_atomic_cons2,
  file.path(out_dir, "nattam_accuracy_atomic_combinedCIandCM_forconsensus_29Jul2026_final.csv")
)

#---- Diagnostics for class-specific metrics with fallbacks ----#

# Summarize distributions 
# Does one class systematically have higher/lower weights?
accuracy_with_fallback_grouped_cons2 %>%
  summarise(
    mean_native = mean(native_metric, na.rm = TRUE),
    mean_tame   = mean(tame_metric, na.rm = TRUE),
    mean_not    = mean(not_metric, na.rm = TRUE),
    
    median_native = median(native_metric, na.rm = TRUE),
    median_tame   = median(tame_metric, na.rm = TRUE),
    median_not    = median(not_metric, na.rm = TRUE)
  )

# with ACI2025:
# mean_native mean_tame mean_not median_native median_tame median_not
#       0.567     0.513    0.628         0.627       0.563      0.738

# exact match (not lenient):
# mean_native mean_tame mean_not median_native median_tame median_not
#       0.537     0.482    0.507         0.581       0.508      0.526

accuracy_with_fallback_grouped_cons2 %>%
  select(native_metric, tame_metric, not_metric) %>%
  pivot_longer(everything(),
               names_to = "class",
               values_to = "acc") %>%
  ggplot(aes(class, acc)) +
  geom_boxplot() +
  theme_bw()

accuracy_with_fallback_grouped_cons2 %>%
  group_by(raster) %>%
  summarise(
    native = mean(native_metric, na.rm = TRUE),
    tame   = mean(tame_metric, na.rm = TRUE),
    not    = mean(not_metric, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(
    c(native, tame, not),
    names_to = "class",
    values_to = "metric"
  ) %>%
  mutate(
    class = factor(
      class,
      levels = c("not", "tame", "native")
    )
  ) %>%
  ggplot(aes(
    x = reorder(raster, metric),
    y = metric,
    fill = class
  )) +
  geom_col(position = "dodge") +
  coord_flip() +
  scale_fill_manual(
    values = c(
      native = "forestgreen",
      tame   = "goldenrod",
      not    = "grey50"
    ),
    breaks = c("native", "tame", "not"),
    #guide = guide_legend(reverse = TRUE)
  ) +
  labs(
    x = "Dataset",
    y = "Mean accuracy metric",
    fill = "Class",
    title = "Class-specific accuracy metrics by dataset"
  ) +
  theme_bw()

# with number of data points

# Counts per raster
count_df <- combined %>%
  count(raster_group, name = "n") %>%
  rename(raster = raster_group)

# Data for plotting
plot_df <- accuracy_with_fallback_grouped_cons2 %>%
  group_by(raster) %>%
  summarise(
    native = mean(native_metric, na.rm = TRUE),
    tame   = mean(tame_metric, na.rm = TRUE),
    not    = mean(not_metric, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(count_df, by = "raster") %>%
  pivot_longer(
    c(native, tame, not),
    names_to = "class",
    values_to = "metric"
  ) %>%
  mutate(
    class = factor(class,
                   levels = c("not", "tame", "native")
    )
  )

# One label position per raster
label_df <- plot_df %>%
  group_by(raster, n) %>%
  summarise(y = max(metric, na.rm = TRUE) + 0.06, .groups = "drop")

p <- ggplot(plot_df, aes(
  x = reorder(raster, metric),
  y = metric,
  fill = class
)) +
  geom_col(position = position_dodge()) +
  geom_text(
    data = label_df,
    aes(x = raster, y = y, 
        label = paste0("n=",scales::comma(n))),
    inherit.aes = FALSE,
    size = 4
  ) +
  coord_flip() +
  scale_fill_manual(
    values = c(
      native = "forestgreen",
      tame   = "goldenrod",
      not    = "grey50"
    ),
    breaks = c("native", "tame", "not")
  ) +
  labs(
    x = "Dataset",
    y = "Mean accuracy metric",
    fill = "Class",
    title = "Class-specific accuracy metrics by dataset"
  ) +
  theme_bw()

p

ggsave(file.path(out_dir, "accuracy_by_dataset_triclass.png"),
  plot = p,
  width = 10,
  height = 8,
  dpi = 600
)

# How often are fallbacks used?
accuracy_with_fallback_grouped %>%
  mutate(
    native_source = case_when(
      !is.na(F1_2) ~ "F1_2",
      !is.na(Macro_F1) ~ "Macro_F1",
      TRUE ~ "Balanced_Accuracy"
    ),
    tame_source = case_when(
      !is.na(F1_1) ~ "F1_1",
      !is.na(Macro_F1) ~ "Macro_F1",
      TRUE ~ "Balanced_Accuracy"
    ),
    not_source = case_when(
      !is.na(F1_0) ~ "F1_0",
      !is.na(Macro_F1) ~ "Macro_F1",
      TRUE ~ "Balanced_Accuracy"
    )
  ) %>%
  summarise(
    native_F1 = sum(native_source == "F1_2"),
    native_macro = sum(native_source == "Macro_F1"),
    native_BA = sum(native_source == "Balanced_Accuracy"),
    
    tame_F1 = sum(tame_source == "F1_1"),
    tame_macro = sum(tame_source == "Macro_F1"),
    tame_BA = sum(tame_source == "Balanced_Accuracy"),
    
    not_F1 = sum(not_source == "F1_0"),
    not_macro = sum(not_source == "Macro_F1"),
    not_BA = sum(not_source == "Balanced_Accuracy")
  )

# with ACI2025:
# native_F1 native_macro native_BA tame_F1 tame_macro tame_BA not_F1 not_macro not_BA
#      1004          637        52    1120        521      52   1039       602     52
# 60–66% of rows use true class-specific information
# ~30–40% fall back to Macro F1
# Only ~3% fall all the way back to Balanced Accuracy

# exact match (not lenient):
# native_F1 native_macro native_BA tame_F1 tame_macro tame_BA not_F1 not_macro not_BA
#      1040          102        58     496        646      58    287       855     58

# Where is fallback occurring?
accuracy_with_fallback_grouped %>%
  filter(is.na(F1_2)) %>%
  count(raster, sort = TRUE)

# with ACI2025:
#   raster      n
# 1 gpw       254
# 2 aci2024   244
# 3 eccc       65
# 4 aci2025    54
# 5 cnvc       19
# 6 cgr        15
# 7 pudl       14
# 8 ngi        12
# 9 mgi         9
# 10 pglc        3

# exact match (not lenient):
#   raster      n
# 1 gpw        50
# 2 eccc       27
# 3 cnvc       19
# 4 cgr        15
# 5 pudl       14
# 6 ngi        12
# 7 mgi         9
# 8 aci2025     8
# 9 aci2024     3
# 10 pglc        3

accuracy_with_fallback_grouped %>%
  filter(is.na(F1_2)) %>%
  count(Ecrgn_n, sort = TRUE) %>%
  print(n = Inf)

# with ACI2025:
#   Ecrgn_n                              n
# 1 Aspen Parkland                     108
# 2 Boreal Transition                  108
# 3 Moist Mixed Grassland               55
# 4 Mid-Boreal Uplands                  53
# 5 Lake Manitoba Plain                 39
# 6 Lake Erie Lowland                   33
# 7 Peace Lowland                       33
# 8 Manitoulin-Lake Simcoe              32
# 9 Mixed Grassland                     31
# 10 Western Alberta Upland              25
# 11 St. Lawrence Lowlands               21
# 12 Interlake Plain                     16
# 13 Lake of the Woods                   15
# 14 Algonquin-Lake Nipissing             9
# 15 Fescue Grassland                     9
# 16 Lac Temiscamingue Lowland            9
# 17 Mid-Boreal Lowland                   8
# 18 Southwest Manitoba Uplands           7
# 19 Thompson-Okanagan Plateau            7
# 20 Abitibi Plains                       6
# 21 Eastern Vancouver Island             6
# 22 Columbia Mountains and Highlands     5
# 23 Thunder Bay-Quetico                  5
# 24 Cascade Ranges                       4
# 25 Clear Hills Upland                   4
# 26 Northern Continental Divide          4
# 27 Pacific Ranges                       4
# 28 Western Boreal                       4
# 29 Fraser Basin                         3
# 30 Frontenac Axis                       3
# 31 Rainy River                          3
# 32 Central Canadian Rocky Mountains     2
# 33 Eastern Continental Ranges           2
# 34 Fraser Plateau                       2
# 35 Lac Seul Upland                      2
# 36 Lake Nipigon                         2
# 37 Lower Mainland                       2
# 38 Okanagan Highland                    2
# 39 Skeena Mountains                     2
# 40 Wabasca Lowland                      2
# 41 Okanagan Range                       1
# 42 Southern Rocky Mountain Trench       1

# exact match (not lenient):
#   Ecrgn_n                         n
# 1 Boreal Transition              28
# 2 Aspen Parkland                 25
# 3 Lake Manitoba Plain            13
# 4 Mid-Boreal Uplands             12
# 5 Lake Erie Lowland              11
# 6 Western Alberta Upland         11
# 7 Manitoulin-Lake Simcoe         10
# 8 St. Lawrence Lowlands           7
# 9 Lake of the Woods               5
# 10 Cascade Ranges                  4
# 11 Abitibi Plains                  3
# 12 Algonquin-Lake Nipissing        3
# 13 Fraser Basin                    3
# 14 Lac Temiscamingue Lowland       3
# 15 Peace Lowland                   3
# 16 Thompson-Okanagan Plateau       3
# 17 Eastern Continental Ranges      2
# 18 Eastern Vancouver Island        2
# 19 Moist Mixed Grassland           2
# 20 Northern Continental Divide     2
# 21 Thunder Bay-Quetico             2
# 22 Frontenac Axis                  1
# 23 Interlake Plain                 1
# 24 Lake Nipigon                    1
# 25 Okanagan Range                  1
# 26 Rainy River                     1
# 27 Southwest Manitoba Uplands      1

# Compare class-specific metrics against macro F1
accuracy_with_fallback_grouped_cons2 %>%
  summarise(
    cor_native =
      cor(native_metric,
          Macro_F1,
          use = "complete.obs"),
    
    cor_tame =
      cor(tame_metric,
          Macro_F1,
          use = "complete.obs"),
    
    cor_not =
      cor(not_metric,
          Macro_F1,
          use = "complete.obs")
  )

# with ACI2025:
# cor_native cor_tame cor_not
#      0.857    0.863   0.845
# Class-specific metrics and Macro F1 generally tell the same story, but not
# so much so that they are redundant - class-specific weighting is introducing 
# meaningful additional information

# exact match (not lenient):
# cor_native cor_tame cor_not
#      0.896    0.941   0.940

# How often do native and tame disagree?
accuracy_with_fallback_grouped_cons2 %>%
  summarise(
    mean_diff =
      mean(native_metric - tame_metric,
           na.rm = TRUE),
    
    median_diff =
      median(native_metric - tame_metric,
             na.rm = TRUE)
  )

# with ACI2025:
# mean_diff median_diff
#    0.0536           0
# Native and tame differ more than the overall F1 summaries suggested
# Median of 0 indicates that in at least half of ecodistricts, native 
# and tame metrics are identitical - probably because both are falling
# back to same macro F1 value
# Positive mean indicates tat when they differ, native tends to be higher

# exact match (not lenient):
# mean_diff median_diff
#    0.0548           0

accuracy_with_fallback_grouped_cons2 %>%
  ggplot(
    aes(native_metric, tame_metric)
  ) +
  geom_point(alpha = 0.3) +
  geom_abline(
    slope = 1,
    intercept = 0,
    colour = "red"
  ) +
  theme_bw()

# Relationship with sample size 
accuracy_with_fallback_grouped_cons2 %>%
  summarise(
    native =
      cor(log10(N_evaluated),
          native_metric,
          use = "complete.obs"),
    
    tame =
      cor(log10(N_evaluated),
          tame_metric,
          use = "complete.obs"),
    
    not =
      cor(log10(N_evaluated),
          not_metric,
          use = "complete.obs")
  )

# with ACI2025:
# native    tame   not
#  0.183 -0.0277 0.208 - not driving class-specific metrics

# exact match (not lenient):
# native    tame    not
#  0.218 -0.0150 0.0283

# Do datasets tend to be conservative about predicting native grass?
# If native metrics are frequently higher despite lower recall, the
# consensus weighting scheme could favor datasets that are very
# selective about calling pixels native
accuracy_with_fallback_grouped_cons2 %>%
  summarise(
    native_gt_tame =
      mean(native_metric > tame_metric,
           na.rm = TRUE),
    
    tame_gt_native =
      mean(tame_metric > native_metric,
           na.rm = TRUE)
  )

# with ACI2025:
# native_gt_tame tame_gt_native
#          0.419          0.215
# Native greater than tame in 41.9% of cases 
# Tame greater than native in 21.5% of cases
# Remainder (36.6%) are ties?
# Suggests a systematic tendency for datasets to achieve higher native scores than tame scores
# Even though mean F1_native ≈ 0.506 and mean F1_tame ≈ 0.498 (almost identical)
# Indicates that similarity in means is hiding a directional tendency?
# Native classifications more often scored as better than tame classifications, 
# although the advantage is usually small?

# exact match (not lenient):
# native_gt_tame tame_gt_native
#          0.238          0.125

#----- Ecodistrict diagnostics ----

# Which ecodistricts consistently have low accuracy?
accuracy_with_fallback_grouped_cons2 %>%
  group_by(Ecds_ID) %>%
  summarise(
    native = mean(native_metric, na.rm = TRUE),
    tame   = mean(tame_metric, na.rm = TRUE),
    not    = mean(not_metric, na.rm = TRUE),
    n_rows = n(),
    .groups = "drop"
  ) %>%
  arrange(native) %>%
  print(n = 40)

# with ACI2025:
#   Ecds_ID native  tame    not n_rows
# 1     619 0.0192 0.701 0.118       1
# 2     628 0.0209 0.192 0.0455      4
# 3     985 0.118  0.118 0.118       1
# 4     997 0.128  0.146 0.325       3
# 5    1017 0.159  0.196 0.177       5
# 6     999 0.171  0.234 0.343       4
# 7     740 0.204  0.329 0.361       9
# 8    1012 0.217  0.234 0.413       3
# 9     995 0.237  0.301 0.410       4
# 10    1019 0.245  0.342 0.294       4
# 11     629 0.265  0.313 0.217       4
# 12     993 0.265  0.265 0.265       1
# 13    1013 0.265  0.265 0.265       1
# 14     960 0.271  0.271 0.271       4
# 15    1018 0.275  0.314 0.294       5
# 16     679 0.275  0.390 0.141       5
# 17     687 0.278  0.324 0.409       7
# 18     701 0.282  0.336 0.348       7
# 19     631 0.297  0.390 0.288       5
# 20     688 0.297  0.378 0.418       6
# 21     766 0.299  0.311 0.563       7
# 22    1006 0.312  0.378 0.208       4
# 23    1005 0.318  0.486 0.415       4
# 24     592 0.324  0.192 0.456       4
# 25     374 0.325  0.325 0.325       2
# 26     752 0.333  0.417 0.538       9
# 27     686 0.333  0.397 0.269       6
# 28     726 0.333  0.303 0.319       7
# 29     844 0.340  0.362 0.604       7
# 30     850 0.341  0.306 0.479       8
# 31     723 0.344  0.370 0.343       6
# 32     841 0.348  0.346 0.557       7
# 33     684 0.349  0.383 0.316       7
# 34     683 0.353  0.326 0.366       7
# 35     729 0.356  0.394 0.470      10
# 36     692 0.357  0.405 0.310       6
# 37     718 0.358  0.402 0.343       7
# 38     623 0.358  0.434 0.281       5
# 39     847 0.359  0.340 0.639       7
# 40     695 0.366  0.373 0.442       6

# exact match (not lenient):
#   Ecds_ID  native   tame    not n_rows
# 1     590 0       0      0           1
# 2     960 0       0      0           4
# 3     622 0.00962 0.207  0.108       2
# 4     628 0.00962 0.207  0.108       2
# 5     690 0.0135  0.0135 0.0135      3
# 6     592 0.0169  0.0169 0.0169      2
# 7     619 0.0192  0.414  0.217       1
# 8     610 0.0339  0.0339 0.0339      1
# 9     611 0.0339  0.0339 0.0339      1
# 10     701 0.0606  0.0846 0.0726      5
# 11     688 0.0629  0.161  0.112       4
# 12     839 0.0870  0.105  0.279       3
# 13     687 0.0873  0.0772 0.0823      5
# 14     623 0.0909  0.0909 0.0909      3
# 15     683 0.0986  0.123  0.111       5
# 16     681 0.101   0.141  0.121       3
# 17     685 0.101   0.141  0.121       3
# 18     689 0.101   0.141  0.121       3
# 19     696 0.101   0.141  0.121       3
# 20     698 0.101   0.141  0.121       3
# 21     706 0.101   0.141  0.121       3
# 22     629 0.105   0.105  0.105       2
# 23     684 0.113   0.113  0.113       5
# 24     708 0.113   0.0977 0.106       4
# 25     703 0.114   0.114  0.114       3
# 26     985 0.118   0.118  0.118       1
# 27     686 0.120   0.120  0.120       4
# 28     682 0.131   0.191  0.161       2
# 29     982 0.133   0.133  0.133       3
# 30     695 0.143   0.143  0.143       4
# 31     692 0.149   0.149  0.149       4
# 32     766 0.181   0.213  0.227       5
# 33     680 0.184   0.184  0.184       3
# 34     850 0.193   0.189  0.221       8
# 35     752 0.194   0.292  0.330       7
# 36     988 0.209   0.209  0.209       2
# 37     760 0.209   0.223  0.315       8
# 38    1016 0.209   0.209  0.209       2
# 39     599 0.217   0.0366 0.127       2
# 40     600 0.217   0.0366 0.127       2

accuracy_with_fallback_grouped_cons2 %>%
  group_by(Ecds_ID) %>%
  summarise(
    native = mean(native_metric, na.rm = TRUE),
    tame   = mean(tame_metric, na.rm = TRUE),
    not    = mean(not_metric, na.rm = TRUE),
    n_rows = n(),
    .groups = "drop"
  ) %>%
  arrange(Ecds_ID) %>%
  print(n = Inf)

# with ACI2025:
#   Ecds_ID native  tame    not n_rows
# 1     374 0.325  0.325 0.325       2
# 2     375 0.388  0.291 0.540       7
# 3     377 0.456  0.385 0.526       3
# 4     378 0.555  0.628 0.483       3
# 5     379 0.574  0.431 0.780       4
# 6     380 0.541  0.534 0.598       5
# 7     381 0.638  0.667 0.609       3
# 8     382 0.691  0.777 0.606       2
# 9     383 0.700  0.757 0.643       3
# 10     391 0.625  0.625 0.625       2
# 11     400 0.641  0.641 0.641       2
# 12     401 0.676  0.676 0.676       2
# 13     402 0.641  0.641 0.641       2
# 14     408 0.778  0.743 0.813       3
# 15     409 0.783  0.765 0.801       3
# 16     410 0.773  0.740 0.805       3
# 17     411 0.678  0.634 0.723       3
# 18     412 0.618  0.601 0.635       3
# 19     413 0.739  0.760 0.717       3
# 20     541 0.839  0.750 0.927       3
# 21     542 0.797  0.705 0.890       3
# 22     543 0.839  0.782 0.897       3
# 23     544 0.838  0.783 0.893       3
# 24     545 0.818  0.719 0.917       3
# 25     546 0.854  0.818 0.890       3
# 26     547 0.802  0.766 0.838       3
# 27     549 0.781  0.749 0.813       3
# 28     550 0.693  0.692 0.693       3
# 29     551 0.819  0.776 0.863       3
# 30     552 0.798  0.754 0.842       3
# 31     553 0.856  0.825 0.887       3
# 32     554 0.853  0.859 0.848       3
# 33     555 0.792  0.766 0.818       3
# 34     556 0.782  0.764 0.800       3
# 35     557 0.772  0.727 0.816       3
# 36     558 0.743  0.704 0.783       3
# 37     559 0.733  0.681 0.785       2
# 38     560 0.793  0.742 0.844       3
# 39     561 0.812  0.725 0.898       3
# 40     562 0.840  0.822 0.858       3
# 41     564 0.808  0.743 0.873       3
# 42     565 0.807  0.707 0.906       3
# 43     566 0.733  0.654 0.812       3
# 44     567 0.795  0.701 0.888       3
# 45     568 0.798  0.708 0.888       3
# 46     569 0.782  0.711 0.852       3
# 47     570 0.848  0.747 0.950       3
# 48     571 0.830  0.727 0.932       3
# 49     572 0.912  0.912 0.912       3
# 50     583 0.752  0.616 0.889       2
# 51     585 0.748  0.619 0.877       2
# 52     586 0.569  0.443 0.627       3
# 53     587 0.642  0.522 0.582       3
# 54     588 0.585  0.356 0.814       2
# 55     589 0.524  0.524 0.524       2
# 56     590 0.393  0.265 0.520       3
# 57     591 0.595  0.493 0.697       4
# 58     592 0.324  0.192 0.456       4
# 59     593 0.622  0.476 0.768       3
# 60     594 0.560  0.304 0.636       3
# 61     595 0.503  0.103 0.723       3
# 62     596 0.481  0.230 0.547       4
# 63     597 0.627  0.432 0.642       3
# 64     598 0.660  0.440 0.881       2
# 65     599 0.478  0.282 0.539       4
# 66     600 0.473  0.363 0.448       4
# 67     601 0.5    0.5   0.5         2
# 68     607 0.818  0.818 0.818       2
# 69     610 0.461  0.400 0.521       3
# 70     611 0.600  0.600 0.600       3
# 71     612 0.391  0.457 0.326       2
# 72     615 0.450  0.436 0.464       4
# 73     616 0.516  0.505 0.527       4
# 74     617 0.728  0.690 0.766       4
# 75     618 0.693  0.577 0.625       3
# 76     619 0.0192 0.701 0.118       1
# 77     622 0.536  0.594 0.477       4
# 78     623 0.358  0.434 0.281       5
# 79     626 0.382  0.661 0.651       2
# 80     628 0.0209 0.192 0.0455      4
# 81     629 0.265  0.313 0.217       4
# 82     630 0.378  0.552 0.205       4
# 83     631 0.297  0.390 0.288       5
# 84     647 0.735  0.757 0.713       3
# 85     650 0.369  0.487 0.251       3
# 86     652 0.692  0.590 0.794       3
# 87     657 0.687  0.595 0.780       3
# 88     658 0.695  0.500 0.890       2
# 89     659 0.695  0.500 0.890       2
# 90     660 0.519  0.466 0.571       4
# 91     661 0.805  0.639 0.971       2
# 92     662 0.601  0.25  0.951       2
# 93     669 0.705  0.663 0.810       3
# 94     672 0.546  0.466 0.843       3
# 95     677 0.468  0.407 0.738       4
# 96     678 0.425  0.383 0.466       5
# 97     679 0.275  0.390 0.141       5
# 98     680 0.476  0.474 0.479       5
# 99     681 0.524  0.448 0.581       5
# 100     682 0.540  0.493 0.562       4
# 101     683 0.353  0.326 0.366       7
# 102     684 0.349  0.383 0.316       7
# 103     685 0.477  0.352 0.584       5
# 104     686 0.333  0.397 0.269       6
# 105     687 0.278  0.324 0.409       7
# 106     688 0.297  0.378 0.418       6
# 107     689 0.507  0.407 0.588       5
# 108     690 0.449  0.388 0.510       5
# 109     691 0.493  0.543 0.849       3
# 110     692 0.357  0.405 0.310       6
# 111     693 0.744  0.521 0.968       3
# 112     694 0.704  0.481 0.926       3
# 113     695 0.366  0.373 0.442       6
# 114     696 0.531  0.452 0.592       5
# 115     697 0.661  0.431 0.891       3
# 116     698 0.503  0.396 0.591       5
# 117     699 0.531  0.618 0.850       3
# 118     700 0.684  0.485 0.884       3
# 119     701 0.282  0.336 0.348       7
# 120     702 0.762  0.542 0.983       3
# 121     703 0.458  0.454 0.463       5
# 122     704 0.788  0.620 0.957       3
# 123     705 0.719  0.457 0.981       3
# 124     706 0.525  0.433 0.599       5
# 125     707 0.734  0.525 0.944       3
# 126     708 0.389  0.364 0.398       6
# 127     709 0.402  0.352 0.624       7
# 128     710 0.699  0.605 0.794       3
# 129     711 0.569  0.465 0.631       4
# 130     714 0.523  0.380 0.666       5
# 131     715 0.653  0.527 0.779       4
# 132     716 0.479  0.319 0.612       6
# 133     717 0.737  0.550 0.924       4
# 134     718 0.358  0.402 0.343       7
# 135     720 0.613  0.659 0.773       4
# 136     723 0.344  0.370 0.343       6
# 137     724 0.484  0.346 0.644       6
# 138     726 0.333  0.303 0.319       7
# 139     727 0.408  0.368 0.471       9
# 140     728 0.382  0.301 0.441       8
# 141     729 0.356  0.394 0.470      10
# 142     730 0.395  0.416 0.409       9
# 143     731 0.449  0.345 0.457       9
# 144     732 0.498  0.361 0.429       8
# 145     733 0.435  0.449 0.622       8
# 146     734 0.571  0.569 0.627       8
# 147     735 0.432  0.300 0.539       8
# 148     736 0.573  0.471 0.680       8
# 149     737 0.506  0.386 0.525       7
# 150     738 0.403  0.342 0.503       9
# 151     739 0.385  0.470 0.467      10
# 152     740 0.204  0.329 0.361       9
# 153     741 0.570  0.349 0.566       9
# 154     742 0.445  0.490 0.668       8
# 155     743 0.519  0.487 0.608       9
# 156     744 0.405  0.423 0.506       8
# 157     745 0.445  0.429 0.565       9
# 158     746 0.464  0.322 0.518       7
# 159     747 0.567  0.421 0.505       8
# 160     748 0.629  0.448 0.607       8
# 161     749 0.502  0.318 0.554       8
# 162     750 0.407  0.388 0.438       8
# 163     751 0.434  0.315 0.413       9
# 164     752 0.333  0.417 0.538       9
# 165     753 0.561  0.429 0.650       7
# 166     754 0.620  0.556 0.744       6
# 167     755 0.477  0.414 0.548       9
# 168     756 0.390  0.408 0.512       9
# 169     757 0.529  0.394 0.569       9
# 170     758 0.506  0.337 0.542       8
# 171     759 0.678  0.500 0.694       7
# 172     760 0.405  0.365 0.542      10
# 173     761 0.575  0.477 0.548       8
# 174     762 0.652  0.529 0.622       9
# 175     763 0.402  0.421 0.493       9
# 176     764 0.520  0.242 0.657       7
# 177     765 0.372  0.348 0.644       7
# 178     766 0.299  0.311 0.563       7
# 179     767 0.584  0.582 0.779       9
# 180     768 0.684  0.230 0.655       9
# 181     769 0.606  0.531 0.570       9
# 182     770 0.643  0.655 0.734       9
# 183     771 0.684  0.544 0.625       9
# 184     772 0.523  0.351 0.620       9
# 185     773 0.581  0.625 0.718       9
# 186     774 0.628  0.401 0.643       9
# 187     775 0.670  0.655 0.737       9
# 188     776 0.527  0.533 0.438      10
# 189     777 0.638  0.535 0.522       8
# 190     778 0.491  0.601 0.522       9
# 191     779 0.402  0.543 0.567       8
# 192     780 0.553  0.598 0.890       9
# 193     781 0.377  0.513 0.578       9
# 194     782 0.692  0.725 0.850       9
# 195     783 0.653  0.735 0.823       9
# 196     784 0.693  0.726 0.765       9
# 197     785 0.555  0.651 0.673       9
# 198     786 0.675  0.617 0.726       8
# 199     787 0.772  0.638 0.806       8
# 200     788 0.709  0.348 0.618       8
# 201     789 0.554  0.622 0.722      10
# 202     790 0.665  0.450 0.595       8
# 203     791 0.696  0.515 0.708       8
# 204     792 0.679  0.662 0.848       9
# 205     793 0.669  0.627 0.752       8
# 206     794 0.705  0.557 0.735       9
# 207     795 0.836  0.795 0.927       9
# 208     796 0.730  0.630 0.686       9
# 209     797 0.720  0.654 0.657       8
# 210     798 0.482  0.470 0.561       9
# 211     799 0.640  0.559 0.582       8
# 212     800 0.681  0.614 0.688       9
# 213     801 0.607  0.561 0.557       7
# 214     802 0.714  0.652 0.778       8
# 215     803 0.707  0.664 0.772      10
# 216     804 0.731  0.460 0.609       8
# 217     805 0.595  0.313 0.618      10
# 218     806 0.666  0.407 0.466       9
# 219     807 0.776  0.683 0.830       9
# 220     808 0.730  0.595 0.744      10
# 221     809 0.682  0.545 0.736      10
# 222     810 0.765  0.697 0.809       9
# 223     811 0.839  0.739 0.858       9
# 224     812 0.733  0.513 0.635       8
# 225     813 0.779  0.627 0.759       9
# 226     814 0.639  0.497 0.593       8
# 227     815 0.655  0.532 0.589      10
# 228     816 0.579  0.714 0.837       9
# 229     817 0.785  0.638 0.758      10
# 230     818 0.711  0.649 0.742       8
# 231     819 0.750  0.639 0.534       9
# 232     820 0.771  0.656 0.782      10
# 233     821 0.645  0.570 0.809      10
# 234     822 0.582  0.682 0.749       9
# 235     823 0.718  0.612 0.829       8
# 236     824 0.651  0.650 0.675      10
# 237     825 0.549  0.666 0.794      10
# 238     826 0.624  0.705 0.884       9
# 239     827 0.695  0.678 0.568      10
# 240     828 0.777  0.682 0.865       8
# 241     829 0.699  0.594 0.706       9
# 242     830 0.743  0.714 0.808      10
# 243     831 0.589  0.669 0.784      10
# 244     832 0.758  0.682 0.750      10
# 245     833 0.724  0.606 0.612      10
# 246     834 0.727  0.725 0.826      10
# 247     835 0.871  0.713 0.845       9
# 248     836 0.650  0.609 0.721       8
# 249     837 0.640  0.617 0.698       9
# 250     838 0.627  0.587 0.535      10
# 251     839 0.378  0.237 0.624       6
# 252     840 0.600  0.368 0.831       4
# 253     841 0.348  0.346 0.557       7
# 254     843 0.477  0.374 0.588       7
# 255     844 0.340  0.362 0.604       7
# 256     846 0.452  0.347 0.400       7
# 257     847 0.359  0.340 0.639       7
# 258     848 0.570  0.480 0.752       5
# 259     849 0.431  0.321 0.629       7
# 260     850 0.341  0.306 0.479       8
# 261     851 0.676  0.522 0.685       7
# 262     852 0.389  0.277 0.648       7
# 263     853 0.990  0.990 0.990       4
# 264     854 0.699  0.434 0.965       4
# 265     855 0.622  0.573 0.696       7
# 266     950 0.504  0.504 0.504       2
# 267     951 0.552  0.552 0.552       2
# 268     955 0.700  0.700 0.700       3
# 269     956 0.703  0.703 0.703       3
# 270     959 0.847  0.847 0.847       2
# 271     960 0.271  0.271 0.271       4
# 272     962 0.678  0.814 0.543       2
# 273     969 0.678  0.814 0.543       2
# 274     971 0.515  0.527 0.513       2
# 275     972 0.662  0.585 0.494       3
# 276     976 0.515  0.527 0.514       2
# 277     978 0.604  0.601 0.540       3
# 278     979 0.866  0.469 0.420       3
# 279     982 0.588  0.680 0.497       3
# 280     985 0.118  0.118 0.118       1
# 281     988 0.405  0.475 0.544       3
# 282     990 0.405  0.475 0.544       3
# 283     991 0.574  0.574 0.574       3
# 284     992 0.549  0.653 0.757       2
# 285     993 0.265  0.265 0.265       1
# 286     995 0.237  0.301 0.410       4
# 287     997 0.128  0.146 0.325       3
# 288     999 0.171  0.234 0.343       4
# 289    1001 0.404  0.529 0.539       4
# 290    1002 0.536  0.536 0.536       1
# 291    1003 0.536  0.536 0.536       1
# 292    1004 0.526  0.565 0.488       4
# 293    1005 0.318  0.486 0.415       4
# 294    1006 0.312  0.378 0.208       4
# 295    1007 0.561  0.637 0.484       4
# 296    1008 0.394  0.492 0.433       4
# 297    1009 0.413  0.592 0.531       4
# 298    1010 0.446  0.573 0.646       4
# 299    1011 0.505  0.505 0.505       3
# 300    1012 0.217  0.234 0.413       3
# 301    1013 0.265  0.265 0.265       1
# 302    1015 0.380  0.403 0.427       3
# 303    1016 0.395  0.503 0.449       4
# 304    1017 0.159  0.196 0.177       5
# 305    1018 0.275  0.314 0.294       5
# 306    1019 0.245  0.342 0.294       4

# exact match (not lenient):
#   Ecds_ID  native   tame    not n_rows
# 1     375 0.277   0.217  0.314       5
# 2     377 0.711   0.711  0.711       1
# 3     378 0.848   0.848  0.848       1
# 4     379 0.459   0.258  0.593       1
# 5     380 0.335   0.268  0.379       3
# 6     381 0.662   0.662  0.662       1
# 7     382 0.797   0.797  0.797       1
# 8     383 0.833   0.833  0.833       1
# 9     391 0.651   0.651  0.651       1
# 10     400 0.658   0.658  0.658       1
# 11     401 0.676   0.676  0.676       1
# 12     402 0.658   0.658  0.658       1
# 13     408 0.688   0.688  0.688       1
# 14     409 0.688   0.688  0.688       1
# 15     410 0.688   0.688  0.688       1
# 16     411 0.576   0.576  0.576       1
# 17     412 0.651   0.651  0.651       1
# 18     413 0.776   0.776  0.776       1
# 19     541 0.62    0.62   0.62        1
# 20     542 0.62    0.62   0.62        1
# 21     543 0.625   0.625  0.625       1
# 22     544 0.553   0.553  0.553       1
# 23     545 0.62    0.62   0.62        1
# 24     546 0.634   0.634  0.634       1
# 25     547 0.701   0.701  0.701       1
# 26     549 0.6     0.6    0.6         1
# 27     550 0.608   0.608  0.608       1
# 28     551 0.592   0.592  0.592       1
# 29     552 0.614   0.614  0.614       1
# 30     553 0.659   0.659  0.659       1
# 31     554 0.829   0.829  0.829       1
# 32     555 0.633   0.633  0.633       1
# 33     556 0.575   0.575  0.575       1
# 34     557 0.394   0.394  0.394       1
# 35     558 0.370   0.370  0.370       1
# 36     560 0.464   0.464  0.464       1
# 37     561 0.411   0.411  0.411       1
# 38     562 0.411   0.411  0.411       1
# 39     564 0.463   0.463  0.463       1
# 40     565 0.329   0.329  0.329       1
# 41     566 0.606   0.606  0.606       1
# 42     567 0.284   0.284  0.284       1
# 43     568 0.415   0.415  0.415       1
# 44     569 0.512   0.512  0.512       1
# 45     570 0.392   0.392  0.392       1
# 46     571 0.411   0.411  0.411       1
# 47     572 0.411   0.411  0.411       1
# 48     586 0.333   0.2    0.267       1
# 49     587 0.4     0.0392 0.220       1
# 50     590 0       0      0           1
# 51     591 0.525   0.525  0.525       2
# 52     592 0.0169  0.0169 0.0169      2
# 53     593 0.417   0.417  0.417       1
# 54     594 0.4     0.0392 0.220       1
# 55     595 0.4     0.0392 0.220       1
# 56     596 0.298   0.0527 0.175       2
# 57     597 0.4     0.0392 0.220       1
# 58     599 0.217   0.0366 0.127       2
# 59     600 0.217   0.0366 0.127       2
# 60     610 0.0339  0.0339 0.0339      1
# 61     611 0.0339  0.0339 0.0339      1
# 62     615 0.307   0.307  0.307       2
# 63     616 0.307   0.307  0.307       2
# 64     617 0.690   0.690  0.690       2
# 65     618 0.588   0.405  0.222       1
# 66     619 0.0192  0.414  0.217       1
# 67     622 0.00962 0.207  0.108       2
# 68     623 0.0909  0.0909 0.0909      3
# 69     628 0.00962 0.207  0.108       2
# 70     629 0.105   0.105  0.105       2
# 71     630 0.333   0.333  0.333       2
# 72     631 0.217   0.180  0.144       5
# 73     650 0.581   0.581  0.581       1
# 74     660 0.307   0.307  0.307       2
# 75     669 0.400   0.439  0.548       1
# 76     672 0.437   0.437  0.437       1
# 77     677 0.400   0.439  0.548       1
# 78     678 0.230   0.230  0.230       3
# 79     679 0.263   0.242  0.253       3
# 80     680 0.184   0.184  0.184       3
# 81     681 0.101   0.141  0.121       3
# 82     682 0.131   0.191  0.161       2
# 83     683 0.0986  0.123  0.111       5
# 84     684 0.113   0.113  0.113       5
# 85     685 0.101   0.141  0.121       3
# 86     686 0.120   0.120  0.120       4
# 87     687 0.0873  0.0772 0.0823      5
# 88     688 0.0629  0.161  0.112       4
# 89     689 0.101   0.141  0.121       3
# 90     690 0.0135  0.0135 0.0135      3
# 91     692 0.149   0.149  0.149       4
# 92     695 0.143   0.143  0.143       4
# 93     696 0.101   0.141  0.121       3
# 94     698 0.101   0.141  0.121       3
# 95     701 0.0606  0.0846 0.0726      5
# 96     703 0.114   0.114  0.114       3
# 97     706 0.101   0.141  0.121       3
# 98     708 0.113   0.0977 0.106       4
# 99     709 0.295   0.325  0.419       5
# 100     710 0.581   0.581  0.581       1
# 101     711 0.366   0.311  0.338       2
# 102     714 0.379   0.317  0.440       3
# 103     715 0.448   0.356  0.541       1
# 104     716 0.249   0.132  0.310       3
# 105     717 0.676   0.466  0.886       1
# 106     718 0.239   0.232  0.287       5
# 107     720 0.4     0.438  0.540       1
# 108     723 0.261   0.258  0.201       4
# 109     724 0.276   0.220  0.376       3
# 110     726 0.251   0.241  0.216       5
# 111     727 0.311   0.326  0.319       9
# 112     728 0.225   0.205  0.215       6
# 113     729 0.241   0.243  0.274       8
# 114     730 0.320   0.328  0.290       9
# 115     731 0.360   0.227  0.275       7
# 116     732 0.375   0.167  0.222       6
# 117     733 0.319   0.349  0.436       6
# 118     734 0.511   0.529  0.585       6
# 119     735 0.343   0.198  0.245       6
# 120     736 0.498   0.450  0.552       6
# 121     737 0.327   0.203  0.206       5
# 122     738 0.337   0.293  0.316       7
# 123     739 0.464   0.394  0.397      10
# 124     740 0.295   0.258  0.271       9
# 125     741 0.450   0.271  0.347       9
# 126     742 0.347   0.466  0.489       6
# 127     743 0.477   0.488  0.488       7
# 128     744 0.287   0.317  0.295       6
# 129     745 0.365   0.343  0.352       9
# 130     746 0.352   0.187  0.269       5
# 131     747 0.485   0.310  0.384       6
# 132     748 0.495   0.330  0.389       6
# 133     749 0.360   0.190  0.270       6
# 134     750 0.434   0.316  0.375       8
# 135     751 0.355   0.201  0.269       7
# 136     752 0.194   0.292  0.330       7
# 137     753 0.423   0.343  0.340       5
# 138     754 0.519   0.456  0.522       4
# 139     755 0.371   0.279  0.342       7
# 140     756 0.261   0.311  0.285       7
# 141     757 0.556   0.420  0.482       9
# 142     758 0.451   0.337  0.404       8
# 143     759 0.542   0.324  0.478       4
# 144     760 0.209   0.223  0.315       8
# 145     761 0.509   0.427  0.423       6
# 146     762 0.598   0.456  0.497       7
# 147     763 0.305   0.300  0.297       7
# 148     764 0.431   0.177  0.368       5
# 149     765 0.326   0.293  0.334       5
# 150     766 0.181   0.213  0.227       5
# 151     767 0.543   0.583  0.668       7
# 152     768 0.672   0.176  0.464       7
# 153     769 0.594   0.473  0.511       7
# 154     770 0.593   0.605  0.604       7
# 155     771 0.698   0.487  0.555       7
# 156     772 0.531   0.406  0.486       9
# 157     773 0.529   0.551  0.626       9
# 158     774 0.560   0.273  0.499       7
# 159     775 0.634   0.580  0.586       7
# 160     776 0.566   0.495  0.400      10
# 161     777 0.713   0.524  0.551       6
# 162     778 0.410   0.518  0.402       7
# 163     779 0.561   0.503  0.522       8
# 164     780 0.625   0.668  0.773       9
# 165     781 0.447   0.500  0.462       9
# 166     782 0.659   0.664  0.720       9
# 167     783 0.588   0.622  0.732       7
# 168     784 0.656   0.690  0.683       7
# 169     785 0.454   0.539  0.523       7
# 170     786 0.683   0.692  0.695       6
# 171     787 0.822   0.737  0.780       6
# 172     788 0.728   0.320  0.461       8
# 173     789 0.627   0.576  0.628      10
# 174     790 0.646   0.218  0.400       6
# 175     791 0.695   0.533  0.615       6
# 176     792 0.658   0.639  0.791       7
# 177     793 0.663   0.539  0.633       6
# 178     794 0.669   0.466  0.566       7
# 179     795 0.880   0.814  0.864       7
# 180     796 0.735   0.594  0.570       7
# 181     797 0.829   0.840  0.746       6
# 182     798 0.586   0.511  0.492       9
# 183     799 0.655   0.575  0.615       6
# 184     800 0.742   0.739  0.693       9
# 185     801 0.613   0.601  0.607       5
# 186     802 0.765   0.792  0.795       6
# 187     803 0.795   0.801  0.733       8
# 188     804 0.827   0.456  0.593       6
# 189     805 0.845   0.456  0.544      10
# 190     806 0.868   0.553  0.584       9
# 191     807 0.827   0.833  0.828       7
# 192     808 0.786   0.767  0.727      10
# 193     809 0.856   0.740  0.780      10
# 194     810 0.784   0.799  0.807       7
# 195     811 0.869   0.850  0.876       7
# 196     812 0.807   0.488  0.585       6
# 197     813 0.803   0.761  0.779       7
# 198     814 0.869   0.578  0.656       8
# 199     815 0.841   0.742  0.695      10
# 200     816 0.779   0.809  0.816       9
# 201     817 0.818   0.763  0.740       8
# 202     818 0.783   0.813  0.758       6
# 203     819 0.784   0.742  0.664       7
# 204     820 0.790   0.754  0.723       8
# 205     821 0.830   0.792  0.819      10
# 206     822 0.779   0.784  0.778       9
# 207     823 0.808   0.807  0.818       6
# 208     824 0.820   0.789  0.712      10
# 209     825 0.740   0.762  0.743      10
# 210     826 0.853   0.854  0.871       9
# 211     827 0.853   0.819  0.724      10
# 212     828 0.854   0.735  0.800       7
# 213     829 0.834   0.754  0.723       9
# 214     830 0.880   0.855  0.856      10
# 215     831 0.773   0.766  0.755      10
# 216     832 0.919   0.861  0.866      10
# 217     833 0.901   0.769  0.730      10
# 218     834 0.917   0.880  0.842      10
# 219     835 0.939   0.884  0.914       7
# 220     836 0.843   0.776  0.797       8
# 221     837 0.823   0.759  0.748       9
# 222     838 0.781   0.697  0.615      10
# 223     839 0.0870  0.105  0.279       3
# 224     840 0.538   0.377  0.699       1
# 225     841 0.355   0.351  0.399       5
# 226     843 0.378   0.374  0.393       5
# 227     844 0.287   0.349  0.402       5
# 228     846 0.398   0.328  0.229       5
# 229     847 0.353   0.378  0.430       5
# 230     848 0.425   0.317  0.493       2
# 231     849 0.296   0.301  0.413       5
# 232     850 0.193   0.189  0.221       8
# 233     851 0.498   0.326  0.493       5
# 234     852 0.323   0.271  0.438       5
# 235     853 0.987   0.987  0.987       1
# 236     854 0.616   0.308  0.925       1
# 237     855 0.389   0.311  0.401       5
# 238     955 1       1      1           1
# 239     956 1       1      1           1
# 240     960 0       0      0           4
# 241     972 0.951   0.706  0.462       1
# 242     978 0.980   0.898  0.817       3
# 243     979 0.978   0.862  0.747       3
# 244     982 0.133   0.133  0.133       3
# 245     985 0.118   0.118  0.118       1
# 246     988 0.209   0.209  0.209       2
# 247     990 0.385   0.385  0.385       3
# 248     991 0.385   0.385  0.385       3
# 249     993 0.265   0.265  0.265       1
# 250     995 0.287   0.287  0.287       2
# 251     997 0.346   0.346  0.346       3
# 252     999 0.337   0.337  0.337       4
# 253    1001 0.609   0.661  0.635       3
# 254    1002 0.536   0.536  0.536       1
# 255    1003 0.536   0.536  0.536       1
# 256    1004 0.577   0.577  0.577       3
# 257    1005 0.427   0.585  0.605       4
# 258    1006 0.479   0.489  0.484       4
# 259    1007 0.496   0.549  0.602       4
# 260    1008 0.372   0.431  0.482       4
# 261    1009 0.498   0.580  0.539       4
# 262    1010 0.498   0.476  0.522       4
# 263    1011 0.516   0.516  0.516       1
# 264    1012 0.265   0.265  0.265       1
# 265    1013 0.265   0.265  0.265       1
# 266    1015 0.435   0.435  0.435       3
# 267    1016 0.209   0.209  0.209       2
# 268    1017 0.299   0.299  0.299       5
# 269    1018 0.449   0.449  0.449       5
# 270    1019 0.424   0.424  0.424       4

# Which ecodistricts show the greatest disagreement among datasets?
accuracy_with_fallback_grouped_cons2 %>%
  group_by(Ecds_ID) %>%
  summarise(
    native_sd = sd(native_metric, na.rm = TRUE),
    tame_sd   = sd(tame_metric, na.rm = TRUE),
    not_sd    = sd(not_metric, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(native_sd)) %>%
  print(n = 20)

# with ACI2025:
#   Ecds_ID native_sd tame_sd not_sd
# 1     601     0.707   0.707 0.707 
# 2     589     0.673   0.673 0.673 
# 3     611     0.503   0.503 0.503 
# 4     992     0.494   0.347 0.200 
# 5     374     0.460   0.460 0.460 
# 6     988     0.429   0.395 0.395 
# 7     990     0.429   0.395 0.395 
# 8     991     0.427   0.427 0.427 
# 9     696     0.420   0.382 0.489 
# 10     681     0.418   0.400 0.479 
# 11     378     0.416   0.426 0.443 
# 12     706     0.414   0.359 0.495 
# 13     689     0.411   0.396 0.485 
# 14     698     0.406   0.385 0.488 
# 15     703     0.402   0.398 0.423 
# 16     814     0.402   0.390 0.303 
# 17     690     0.398   0.361 0.462 
# 18     701     0.396   0.396 0.426 
# 19     971     0.392   0.459 0.0896
# 20     976     0.392   0.459 0.0886

# exact match (not lenient):
#   Ecds_ID native_sd tame_sd not_sd
# 1     630     0.471  0.471   0.471
# 2     678     0.398  0.398   0.398
# 3     591     0.389  0.389   0.389
# 4    1005     0.388  0.368   0.337
# 5     615     0.387  0.387   0.387
# 6     616     0.387  0.387   0.387
# 7     660     0.387  0.387   0.387
# 8    1004     0.381  0.381   0.381
# 9     596     0.374  0.0265  0.200
# 10     997     0.371  0.371   0.371
# 11     734     0.338  0.310   0.255
# 12     759     0.336  0.428   0.327
# 13     679     0.331  0.347   0.338
# 14     990     0.319  0.319   0.319
# 15     991     0.319  0.319   0.319
# 16     801     0.315  0.341   0.328
# 17    1017     0.307  0.307   0.307
# 18     711     0.304  0.382   0.343
# 19     999     0.304  0.304   0.304
# 20     714     0.300  0.274   0.357

accuracy_with_fallback_grouped_cons2 %>%
  group_by(Ecds_ID) %>%
  summarise(
    native_sd = sd(native_metric, na.rm = TRUE),
    tame_sd   = sd(tame_metric, na.rm = TRUE),
    not_sd    = sd(not_metric, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Ecds_ID) %>%
  print(n = Inf)

# with ACI2025:
#   Ecds_ID native_sd  tame_sd   not_sd
# 1     374   0.460    0.460    0.460  
# 2     375   0.282    0.252    0.413  
# 3     377   0.0649   0.270    0.349  
# 4     378   0.416    0.426    0.443  
# 5     379   0.138    0.291    0.163  
# 6     380   0.352    0.387    0.353  
# 7     381   0.295    0.272    0.324  
# 8     382   0.227    0.107    0.348  
# 9     383   0.170    0.119    0.251  
# 10     391   0.0358   0.0358   0.0358 
# 11     400   0.0246   0.0246   0.0246 
# 12     401   0        0        0      
# 13     402   0.0246   0.0246   0.0246 
# 14     408   0.117    0.162    0.0907 
# 15     409   0.107    0.126    0.0953 
# 16     410   0.117    0.158    0.0939 
# 17     411   0.168    0.221    0.140  
# 18     412   0.247    0.256    0.240  
# 19     413   0.210    0.196    0.229  
# 20     541   0.181    0.335    0.0278 
# 21     542   0.159    0.315    0.0562 
# 22     543   0.108    0.205    0.0367 
# 23     544   0.136    0.232    0.0412 
# 24     545   0.165    0.335    0.0277 
# 25     546   0.117    0.179    0.0546 
# 26     547   0.155    0.208    0.113  
# 27     549   0.193    0.247    0.142  
# 28     550   0.219    0.220    0.219  
# 29     551   0.138    0.212    0.0673 
# 30     552   0.169    0.242    0.104  
# 31     553   0.132    0.186    0.0787 
# 32     554   0.100    0.0951   0.106  
# 33     555   0.161    0.203    0.122  
# 34     556   0.160    0.192    0.129  
# 35     557   0.222    0.299    0.146  
# 36     558   0.230    0.298    0.162  
# 37     559   0.185    0.258    0.111  
# 38     560   0.183    0.271    0.0961 
# 39     561   0.161    0.307    0.0524 
# 40     562   0.117    0.145    0.0902 
# 41     564   0.168    0.280    0.0550 
# 42     565   0.177    0.345    0.0504 
# 43     566   0.228    0.275    0.256  
# 44     567   0.215    0.376    0.0566 
# 45     568   0.155    0.294    0.101  
# 46     569   0.125    0.239    0.0658 
# 47     570   0.177    0.350    0.0412 
# 48     571   0.180    0.349    0.0775 
# 49     572   0.0676   0.0676   0.0676 
# 50     583   0.135    0.323    0.0536 
# 51     585   0.0918   0.309    0.126  
# 52     586   0.210    0.222    0.356  
# 53     587   0.316    0.480    0.393  
# 54     588   0.272    0.359    0.184  
# 55     589   0.673    0.673    0.673  
# 56     590   0.373    0.331    0.456  
# 57     591   0.282    0.355    0.301  
# 58     592   0.361    0.231    0.507  
# 59     593   0.206    0.191    0.305  
# 60     594   0.232    0.389    0.365  
# 61     595   0.0986   0.0855   0.437  
# 62     596   0.303    0.214    0.444  
# 63     597   0.201    0.345    0.367  
# 64     598   0.150    0.262    0.0380 
# 65     599   0.349    0.324    0.482  
# 66     600   0.331    0.381    0.380  
# 67     601   0.707    0.707    0.707  
# 68     607   0.257    0.257    0.257  
# 69     610   0.384    0.356    0.425  
# 70     611   0.503    0.503    0.503  
# 71     612   0.0607   0.289    0.168  
# 72     615   0.280    0.293    0.338  
# 73     616   0.348    0.341    0.422  
# 74     617   0.0808   0.172    0.110  
# 75     618   0.157    0.259    0.351  
# 76     619  NA       NA       NA      
# 77     622   0.368    0.398    0.365  
# 78     623   0.340    0.397    0.328  
# 79     626   0.361    0.0652   0.0118 
# 80     628   0.0304   0.341    0.0569 
# 81     629   0.242    0.295    0.222  
# 82     630   0.263    0.369    0.189  
# 83     631   0.304    0.258    0.236  
# 84     647   0.0562   0.250    0.200  
# 85     650   0.278    0.269    0.305  
# 86     652   0.0982   0.240    0.167  
# 87     657   0.183    0.351    0.162  
# 88     658   0.139    0.257    0.0201 
# 89     659   0.139    0.257    0.0201 
# 90     660   0.343    0.408    0.375  
# 91     661   0.105    0.196    0.0141 
# 92     662   0.00494  0.0505   0.0406 
# 93     669   0.266    0.223    0.234  
# 94     672   0.167    0.323    0.217  
# 95     677   0.0963   0.287    0.196  
# 96     678   0.392    0.366    0.449  
# 97     679   0.189    0.326    0.124  
# 98     680   0.289    0.312    0.315  
# 99     681   0.418    0.400    0.479  
# 100     682   0.320    0.379    0.400  
# 101     683   0.341    0.333    0.378  
# 102     684   0.390    0.425    0.385  
# 103     685   0.388    0.367    0.481  
# 104     686   0.330    0.367    0.321  
# 105     687   0.332    0.324    0.427  
# 106     688   0.365    0.389    0.421  
# 107     689   0.411    0.396    0.485  
# 108     690   0.398    0.361    0.462  
# 109     691   0.294    0.260    0.124  
# 110     692   0.382    0.432    0.366  
# 111     693   0.137    0.265    0.0121 
# 112     694   0.0751   0.165    0.0543 
# 113     695   0.376    0.359    0.443  
# 114     696   0.420    0.382    0.489  
# 115     697   0.0926   0.139    0.0595 
# 116     698   0.406    0.385    0.488  
# 117     699   0.349    0.290    0.124  
# 118     700   0.0852   0.167    0.0833 
# 119     701   0.396    0.396    0.426  
# 120     702   0.0645   0.125    0.00502
# 121     703   0.402    0.398    0.423  
# 122     704   0.131    0.262    0.00994
# 123     705   0.0483   0.0906   0.0145 
# 124     706   0.414    0.359    0.495  
# 125     707   0.107    0.226    0.0398 
# 126     708   0.376    0.378    0.423  
# 127     709   0.296    0.264    0.442  
# 128     710   0.115    0.247    0.168  
# 129     711   0.304    0.369    0.387  
# 130     714   0.299    0.248    0.373  
# 131     715   0.205    0.332    0.161  
# 132     716   0.318    0.291    0.428  
# 133     717   0.141    0.266    0.0331 
# 134     718   0.274    0.349    0.281  
# 135     720   0.293    0.290    0.242  
# 136     723   0.204    0.275    0.323  
# 137     724   0.277    0.258    0.354  
# 138     726   0.208    0.265    0.325  
# 139     727   0.369    0.328    0.424  
# 140     728   0.343    0.288    0.426  
# 141     729   0.293    0.282    0.331  
# 142     730   0.266    0.322    0.329  
# 143     731   0.272    0.307    0.358  
# 144     732   0.246    0.346    0.284  
# 145     733   0.288    0.235    0.310  
# 146     734   0.307    0.281    0.322  
# 147     735   0.259    0.226    0.378  
# 148     736   0.218    0.206    0.249  
# 149     737   0.288    0.325    0.365  
# 150     738   0.267    0.299    0.384  
# 151     739   0.288    0.271    0.319  
# 152     740   0.235    0.292    0.322  
# 153     741   0.236    0.298    0.362  
# 154     742   0.273    0.102    0.227  
# 155     743   0.164    0.203    0.308  
# 156     744   0.279    0.280    0.372  
# 157     745   0.309    0.285    0.372  
# 158     746   0.289    0.222    0.390  
# 159     747   0.246    0.341    0.341  
# 160     748   0.239    0.347    0.360  
# 161     749   0.271    0.252    0.377  
# 162     750   0.276    0.296    0.294  
# 163     751   0.270    0.247    0.310  
# 164     752   0.290    0.273    0.380  
# 165     753   0.252    0.225    0.376  
# 166     754   0.180    0.146    0.217  
# 167     755   0.294    0.279    0.402  
# 168     756   0.293    0.273    0.353  
# 169     757   0.307    0.301    0.316  
# 170     758   0.190    0.235    0.380  
# 171     759   0.318    0.450    0.355  
# 172     760   0.313    0.285    0.377  
# 173     761   0.207    0.227    0.307  
# 174     762   0.267    0.320    0.344  
# 175     763   0.231    0.284    0.356  
# 176     764   0.187    0.181    0.325  
# 177     765   0.260    0.243    0.404  
# 178     766   0.270    0.253    0.444  
# 179     767   0.160    0.215    0.244  
# 180     768   0.139    0.233    0.374  
# 181     769   0.114    0.215    0.210  
# 182     770   0.185    0.143    0.198  
# 183     771   0.141    0.306    0.280  
# 184     772   0.170    0.216    0.332  
# 185     773   0.227    0.230    0.262  
# 186     774   0.239    0.379    0.364  
# 187     775   0.167    0.127    0.173  
# 188     776   0.187    0.258    0.362  
# 189     777   0.185    0.313    0.294  
# 190     778   0.271    0.263    0.319  
# 191     779   0.312    0.284    0.242  
# 192     780   0.281    0.216    0.132  
# 193     781   0.208    0.242    0.343  
# 194     782   0.167    0.0935   0.130  
# 195     783   0.243    0.148    0.152  
# 196     784   0.210    0.162    0.190  
# 197     785   0.217    0.185    0.217  
# 198     786   0.133    0.194    0.156  
# 199     787   0.134    0.351    0.212  
# 200     788   0.153    0.405    0.380  
# 201     789   0.185    0.213    0.200  
# 202     790   0.279    0.481    0.364  
# 203     791   0.149    0.356    0.355  
# 204     792   0.152    0.142    0.163  
# 205     793   0.131    0.207    0.162  
# 206     794   0.165    0.333    0.336  
# 207     795   0.145    0.149    0.0526 
# 208     796   0.160    0.208    0.265  
# 209     797   0.195    0.337    0.322  
# 210     798   0.257    0.247    0.368  
# 211     799   0.0731   0.245    0.246  
# 212     800   0.0919   0.215    0.261  
# 213     801   0.248    0.268    0.278  
# 214     802   0.165    0.281    0.171  
# 215     803   0.219    0.247    0.280  
# 216     804   0.213    0.396    0.267  
# 217     805   0.314    0.308    0.395  
# 218     806   0.382    0.368    0.300  
# 219     807   0.148    0.274    0.161  
# 220     808   0.0803   0.219    0.264  
# 221     809   0.346    0.388    0.255  
# 222     810   0.170    0.240    0.177  
# 223     811   0.141    0.233    0.143  
# 224     812   0.164    0.343    0.270  
# 225     813   0.122    0.221    0.200  
# 226     814   0.402    0.390    0.303  
# 227     815   0.349    0.329    0.335  
# 228     816   0.289    0.199    0.124  
# 229     817   0.159    0.281    0.294  
# 230     818   0.157    0.299    0.209  
# 231     819   0.253    0.307    0.435  
# 232     820   0.130    0.230    0.284  
# 233     821   0.371    0.341    0.226  
# 234     822   0.345    0.244    0.222  
# 235     823   0.142    0.279    0.117  
# 236     824   0.230    0.266    0.282  
# 237     825   0.286    0.219    0.226  
# 238     826   0.360    0.273    0.112  
# 239     827   0.291    0.257    0.376  
# 240     828   0.148    0.300    0.0908 
# 241     829   0.285    0.325    0.340  
# 242     830   0.238    0.242    0.166  
# 243     831   0.309    0.250    0.293  
# 244     832   0.326    0.313    0.247  
# 245     833   0.341    0.339    0.346  
# 246     834   0.368    0.316    0.285  
# 247     835   0.182    0.330    0.162  
# 248     836   0.365    0.250    0.202  
# 249     837   0.363    0.266    0.235  
# 250     838   0.310    0.255    0.339  
# 251     839   0.353    0.298    0.486  
# 252     840   0.0797   0.213    0.118  
# 253     841   0.217    0.223    0.303  
# 254     843   0.280    0.282    0.367  
# 255     844   0.209    0.248    0.394  
# 256     846   0.223    0.220    0.343  
# 257     847   0.218    0.279    0.430  
# 258     848   0.208    0.242    0.341  
# 259     849   0.282    0.217    0.426  
# 260     850   0.250    0.268    0.382  
# 261     851   0.339    0.348    0.371  
# 262     852   0.202    0.201    0.404  
# 263     853   0.00297  0.00297  0.00297
# 264     854   0.160    0.311    0.0269 
# 265     855   0.352    0.437    0.358  
# 266     950   0.148    0.148    0.148  
# 267     951   0.125    0.125    0.125  
# 268     955   0.261    0.261    0.261  
# 269     956   0.258    0.258    0.258  
# 270     959   0.0630   0.0630   0.0630 
# 271     960   0.322    0.322    0.322  
# 272     962   0.0702   0.123    0.0171 
# 273     969   0.0702   0.123    0.0171 
# 274     971   0.392    0.459    0.0896 
# 275     972   0.327    0.386    0.0571 
# 276     976   0.392    0.459    0.0886 
# 277     978   0.317    0.131    0.0887 
# 278     979   0.0785   0.271    0.151  
# 279     982   0.169    0.255    0.0845 
# 280     985  NA       NA       NA      
# 281     988   0.429    0.395    0.395  
# 282     990   0.429    0.395    0.395  
# 283     991   0.427    0.427    0.427  
# 284     992   0.494    0.347    0.200  
# 285     993  NA       NA       NA      
# 286     995   0.0701   0.148    0.137  
# 287     997   0.125    0.147    0.302  
# 288     999   0.132    0.214    0.249  
# 289    1001   0.232    0.238    0.140  
# 290    1002  NA       NA       NA      
# 291    1003  NA       NA       NA      
# 292    1004   0.296    0.312    0.298  
# 293    1005   0.262    0.240    0.200  
# 294    1006   0.266    0.336    0.0958 
# 295    1007   0.203    0.230    0.299  
# 296    1008   0.212    0.266    0.192  
# 297    1009   0.114    0.150    0.207  
# 298    1010   0.269    0.154    0.332  
# 299    1011   0.0370   0.0370   0.0370 
# 300    1012   0.0701   0.0793   0.168  
# 301    1013  NA       NA       NA      
# 302    1015   0.344    0.317    0.294  
# 303    1016   0.278    0.214    0.223  
# 304    1017   0.186    0.207    0.194  
# 305    1018   0.162    0.172    0.164  
# 306    1019   0.0935   0.141    0.110 

# exact match (not lenient):
#   Ecds_ID native_sd  tame_sd  not_sd
# 1     375    0.254   0.225    0.302 
# 2     377   NA      NA       NA     
# 3     378   NA      NA       NA     
# 4     379   NA      NA       NA     
# 5     380    0.293   0.273    0.329 
# 6     381   NA      NA       NA     
# 7     382   NA      NA       NA     
# 8     383   NA      NA       NA     
# 9     391   NA      NA       NA     
# 10     400   NA      NA       NA     
# 11     401   NA      NA       NA     
# 12     402   NA      NA       NA     
# 13     408   NA      NA       NA     
# 14     409   NA      NA       NA     
# 15     410   NA      NA       NA     
# 16     411   NA      NA       NA     
# 17     412   NA      NA       NA     
# 18     413   NA      NA       NA     
# 19     541   NA      NA       NA     
# 20     542   NA      NA       NA     
# 21     543   NA      NA       NA     
# 22     544   NA      NA       NA     
# 23     545   NA      NA       NA     
# 24     546   NA      NA       NA     
# 25     547   NA      NA       NA     
# 26     549   NA      NA       NA     
# 27     550   NA      NA       NA     
# 28     551   NA      NA       NA     
# 29     552   NA      NA       NA     
# 30     553   NA      NA       NA     
# 31     554   NA      NA       NA     
# 32     555   NA      NA       NA     
# 33     556   NA      NA       NA     
# 34     557   NA      NA       NA     
# 35     558   NA      NA       NA     
# 36     560   NA      NA       NA     
# 37     561   NA      NA       NA     
# 38     562   NA      NA       NA     
# 39     564   NA      NA       NA     
# 40     565   NA      NA       NA     
# 41     566   NA      NA       NA     
# 42     567   NA      NA       NA     
# 43     568   NA      NA       NA     
# 44     569   NA      NA       NA     
# 45     570   NA      NA       NA     
# 46     571   NA      NA       NA     
# 47     572   NA      NA       NA     
# 48     586   NA      NA       NA     
# 49     587   NA      NA       NA     
# 50     590   NA      NA       NA     
# 51     591    0.389   0.389    0.389 
# 52     592    0.0240  0.0240   0.0240
# 53     593   NA      NA       NA     
# 54     594   NA      NA       NA     
# 55     595   NA      NA       NA     
# 56     596    0.374   0.0265   0.200 
# 57     597   NA      NA       NA     
# 58     599    0.259   0.00376  0.131 
# 59     600    0.259   0.00376  0.131 
# 60     610   NA      NA       NA     
# 61     611   NA      NA       NA     
# 62     615    0.387   0.387    0.387 
# 63     616    0.387   0.387    0.387 
# 64     617    0.155   0.155    0.155 
# 65     618   NA      NA       NA     
# 66     619   NA      NA       NA     
# 67     622    0.0136  0.293    0.153 
# 68     623    0.157   0.157    0.157 
# 69     628    0.0136  0.293    0.153 
# 70     629    0.149   0.149    0.149 
# 71     630    0.471   0.471    0.471 
# 72     631    0.279   0.225    0.191 
# 73     650   NA      NA       NA     
# 74     660    0.387   0.387    0.387 
# 75     669   NA      NA       NA     
# 76     672   NA      NA       NA     
# 77     677   NA      NA       NA     
# 78     678    0.398   0.398    0.398 
# 79     679    0.331   0.347    0.338 
# 80     680    0.147   0.147    0.147 
# 81     681    0.0582  0.172    0.111 
# 82     682    0.0360  0.210    0.123 
# 83     683    0.0790  0.141    0.105 
# 84     684    0.168   0.168    0.168 
# 85     685    0.0582  0.172    0.111 
# 86     686    0.119   0.119    0.119 
# 87     687    0.0760  0.0712   0.0727
# 88     688    0.0519  0.267    0.150 
# 89     689    0.0582  0.172    0.111 
# 90     690    0.0233  0.0233   0.0233
# 91     692    0.272   0.272    0.272 
# 92     695    0.193   0.193    0.193 
# 93     696    0.0582  0.172    0.111 
# 94     698    0.0582  0.172    0.111 
# 95     701    0.0689  0.144    0.103 
# 96     703    0.164   0.164    0.164 
# 97     706    0.0582  0.172    0.111 
# 98     708    0.137   0.141    0.138 
# 99     709    0.279   0.286    0.392 
# 100     710   NA      NA       NA     
# 101     711    0.304   0.382    0.343 
# 102     714    0.300   0.274    0.357 
# 103     715   NA      NA       NA     
# 104     716    0.277   0.165    0.426 
# 105     717   NA      NA       NA     
# 106     718    0.171   0.197    0.197 
# 107     720   NA      NA       NA     
# 108     723    0.141   0.193    0.0946
# 109     724    0.223   0.201    0.307 
# 110     726    0.189   0.210    0.190 
# 111     727    0.280   0.279    0.274 
# 112     728    0.221   0.215    0.216 
# 113     729    0.184   0.184    0.173 
# 114     730    0.168   0.240    0.161 
# 115     731    0.236   0.237    0.189 
# 116     732    0.173   0.136    0.142 
# 117     733    0.228   0.200    0.215 
# 118     734    0.338   0.310    0.255 
# 119     735    0.204   0.152    0.184 
# 120     736    0.198   0.251    0.143 
# 121     737    0.139   0.172    0.150 
# 122     738    0.260   0.271    0.260 
# 123     739    0.231   0.287    0.247 
# 124     740    0.236   0.199    0.178 
# 125     741    0.150   0.218    0.230 
# 126     742    0.241   0.116    0.116 
# 127     743    0.162   0.202    0.258 
# 128     744    0.204   0.234    0.207 
# 129     745    0.235   0.245    0.257 
# 130     746    0.261   0.136    0.172 
# 131     747    0.225   0.314    0.269 
# 132     748    0.154   0.312    0.258 
# 133     749    0.178   0.145    0.231 
# 134     750    0.222   0.250    0.209 
# 135     751    0.234   0.109    0.117 
# 136     752    0.100   0.181    0.265 
# 137     753    0.171   0.0818   0.323 
# 138     754    0.113   0.176    0.105 
# 139     755    0.234   0.207    0.298 
# 140     756    0.163   0.226    0.226 
# 141     757    0.273   0.272    0.248 
# 142     758    0.168   0.226    0.271 
# 143     759    0.336   0.428    0.327 
# 144     760    0.185   0.213    0.297 
# 145     761    0.191   0.230    0.237 
# 146     762    0.281   0.337    0.328 
# 147     763    0.115   0.160    0.203 
# 148     764    0.169   0.140    0.291 
# 149     765    0.219   0.200    0.368 
# 150     766    0.176   0.170    0.326 
# 151     767    0.154   0.223    0.225 
# 152     768    0.146   0.233    0.333 
# 153     769    0.128   0.253    0.144 
# 154     770    0.178   0.177    0.226 
# 155     771    0.147   0.348    0.252 
# 156     772    0.154   0.250    0.291 
# 157     773    0.163   0.212    0.201 
# 158     774    0.228   0.339    0.331 
# 159     775    0.172   0.213    0.237 
# 160     776    0.202   0.284    0.314 
# 161     777    0.143   0.351    0.287 
# 162     778    0.247   0.274    0.249 
# 163     779    0.210   0.284    0.189 
# 164     780    0.228   0.135    0.150 
# 165     781    0.193   0.209    0.260 
# 166     782    0.153   0.154    0.154 
# 167     783    0.225   0.263    0.151 
# 168     784    0.217   0.196    0.191 
# 169     785    0.0954  0.130    0.118 
# 170     786    0.150   0.134    0.129 
# 171     787    0.131   0.359    0.240 
# 172     788    0.0954  0.369    0.286 
# 173     789    0.121   0.231    0.171 
# 174     790    0.181   0.385    0.315 
# 175     791    0.179   0.416    0.368 
# 176     792    0.134   0.167    0.151 
# 177     793    0.140   0.294    0.158 
# 178     794    0.138   0.337    0.358 
# 179     795    0.0772  0.161    0.119 
# 180     796    0.166   0.293    0.297 
# 181     797    0.0993  0.0861   0.276 
# 182     798    0.269   0.281    0.321 
# 183     799    0.0873  0.262    0.138 
# 184     800    0.102   0.134    0.260 
# 185     801    0.315   0.341    0.328 
# 186     802    0.180   0.169    0.170 
# 187     803    0.124   0.132    0.286 
# 188     804    0.154   0.422    0.270 
# 189     805    0.0936  0.417    0.352 
# 190     806    0.0933  0.418    0.347 
# 191     807    0.155   0.154    0.155 
# 192     808    0.0874  0.123    0.252 
# 193     809    0.163   0.378    0.275 
# 194     810    0.180   0.161    0.154 
# 195     811    0.164   0.182    0.162 
# 196     812    0.135   0.407    0.265 
# 197     813    0.130   0.179    0.150 
# 198     814    0.114   0.425    0.302 
# 199     815    0.0851  0.290    0.322 
# 200     816    0.147   0.106    0.106 
# 201     817    0.148   0.194    0.316 
# 202     818    0.142   0.109    0.189 
# 203     819    0.258   0.259    0.366 
# 204     820    0.123   0.164    0.297 
# 205     821    0.200   0.264    0.228 
# 206     822    0.191   0.186    0.193 
# 207     823    0.139   0.140    0.134 
# 208     824    0.108   0.178    0.286 
# 209     825    0.151   0.120    0.210 
# 210     826    0.108   0.106    0.0960
# 211     827    0.0779  0.134    0.315 
# 212     828    0.0895  0.317    0.165 
# 213     829    0.169   0.312    0.345 
# 214     830    0.0631  0.0904   0.102 
# 215     831    0.181   0.173    0.282 
# 216     832    0.0422  0.132    0.135 
# 217     833    0.0740  0.288    0.326 
# 218     834    0.0722  0.146    0.284 
# 219     835    0.0562  0.198    0.119 
# 220     836    0.178   0.252    0.214 
# 221     837    0.0778  0.183    0.204 
# 222     838    0.0950  0.197    0.315 
# 223     839    0.151   0.182    0.484 
# 224     840   NA      NA       NA     
# 225     841    0.196   0.217    0.203 
# 226     843    0.255   0.227    0.267 
# 227     844    0.188   0.217    0.313 
# 228     846    0.225   0.161    0.218 
# 229     847    0.226   0.278    0.375 
# 230     848    0.211   0.183    0.476 
# 231     849    0.187   0.193    0.352 
# 232     850    0.220   0.226    0.251 
# 233     851    0.282   0.163    0.317 
# 234     852    0.184   0.200    0.328 
# 235     853   NA      NA       NA     
# 236     854   NA      NA       NA     
# 237     855    0.185   0.313    0.358 
# 238     955   NA      NA       NA     
# 239     956   NA      NA       NA     
# 240     960    0       0        0     
# 241     972   NA      NA       NA     
# 242     978    0.0256  0.166    0.307 
# 243     979    0.0310  0.231    0.430 
# 244     982    0.231   0.231    0.231 
# 245     985   NA      NA       NA     
# 246     988    0.129   0.129    0.129 
# 247     990    0.319   0.319    0.319 
# 248     991    0.319   0.319    0.319 
# 249     993   NA      NA       NA     
# 250     995    0.0315  0.0315   0.0315
# 251     997    0.371   0.371    0.371 
# 252     999    0.304   0.304    0.304 
# 253    1001    0.0781  0.0998   0.0773
# 254    1002   NA      NA       NA     
# 255    1003   NA      NA       NA     
# 256    1004    0.381   0.381    0.381 
# 257    1005    0.388   0.368    0.337 
# 258    1006    0.287   0.294    0.290 
# 259    1007    0.282   0.218    0.196 
# 260    1008    0.272   0.254    0.265 
# 261    1009    0.137   0.201    0.151 
# 262    1010    0.286   0.219    0.183 
# 263    1011   NA      NA       NA     
# 264    1012   NA      NA       NA     
# 265    1013   NA      NA       NA     
# 266    1015    0.264   0.264    0.264 
# 267    1016    0.119   0.119    0.119 
# 268    1017    0.307   0.307    0.307 
# 269    1018    0.157   0.157    0.157 
# 270    1019    0.196   0.196    0.196 

# Are accuracy estimates dominated by fallback values?
accuracy_with_fallback_grouped %>%
  mutate(
    native_is_fallback = is.na(F1_2),
    tame_is_fallback   = is.na(F1_1),
    not_is_fallback    = is.na(F1_0)
  ) %>%
  group_by(Ecds_ID) %>%
  summarise(
    native_fallback_pct =
      mean(native_is_fallback) * 100,
    tame_fallback_pct =
      mean(tame_is_fallback) * 100,
    not_fallback_pct =
      mean(not_is_fallback) * 100) %>%
  print(n = Inf)

# with ACI2025:
#   Ecds_ID native_fallback_pct tame_fallback_pct not_fallback_pct
# 1     374               100                50              100  
# 2     375                57.1              42.9             42.9
# 3     377               100                 0               33.3
# 4     378               100                33.3             33.3
# 5     379                50                 0                0  
# 6     380                60                60               20  
# 7     381               100                33.3             33.3
# 8     382               100                 0               50  
# 9     383               100                33.3             33.3
# 10     391               100                 0              100  
# 11     400               100                 0              100  
# 12     401               100                 0              100  
# 13     402               100                 0              100  
# 14     408               100                33.3             33.3
# 15     409               100                33.3             33.3
# 16     410               100                33.3             33.3
# 17     411               100                33.3             33.3
# 18     412               100                33.3             33.3
# 19     413               100                33.3             33.3
# 20     541               100                33.3             33.3
# 21     542               100                33.3             33.3
# 22     543               100                33.3             33.3
# 23     544               100                33.3             33.3
# 24     545               100                33.3             33.3
# 25     546               100                33.3             33.3
# 26     547               100                33.3             33.3
# 27     549               100                33.3             33.3
# 28     550               100                33.3             33.3
# 29     551               100                33.3             33.3
# 30     552               100                33.3             33.3
# 31     553               100                33.3             33.3
# 32     554               100                33.3             33.3
# 33     555               100                33.3             33.3
# 34     556               100                33.3             33.3
# 35     557               100                33.3             33.3
# 36     558               100                33.3             33.3
# 37     559               100                50                0  
# 38     560               100                33.3             33.3
# 39     561               100                33.3             33.3
# 40     562               100                33.3             33.3
# 41     564               100                33.3             33.3
# 42     565               100                33.3             33.3
# 43     566               100                33.3             33.3
# 44     567               100                33.3             33.3
# 45     568               100                33.3             33.3
# 46     569               100                33.3             33.3
# 47     570               100                33.3             33.3
# 48     571               100                33.3             33.3
# 49     572               100                66.7             33.3
# 50     583               100                 0                0  
# 51     585               100                 0                0  
# 52     586                66.7              33.3             33.3
# 53     587                66.7              66.7             33.3
# 54     588               100                 0                0  
# 55     589               100               100                0  
# 56     590               100                33.3             33.3
# 57     591                75                50               25  
# 58     592                75                50               50  
# 59     593                66.7              33.3             33.3
# 60     594                66.7               0               33.3
# 61     595                66.7               0               33.3
# 62     596                50                25               50  
# 63     597                66.7               0               33.3
# 64     598               100                 0                0  
# 65     599                50                25               50  
# 66     600                50                25               50  
# 67     601               100               100               50  
# 68     607               100               100                0  
# 69     610                66.7              33.3             33.3
# 70     611                66.7             100               33.3
# 71     612               100                 0                0  
# 72     615                75                25               25  
# 73     616                75                25               25  
# 74     617               100                25                0  
# 75     618                66.7              33.3              0  
# 76     619                 0                 0                0  
# 77     622               100                25               25  
# 78     623               100                40               40  
# 79     626                50                 0                0  
# 80     628                75                50               75  
# 81     629               100                25               50  
# 82     630               100                25               25  
# 83     631                40                40               80  
# 84     647               100                 0                0  
# 85     650               100                 0                0  
# 86     652               100                 0                0  
# 87     657               100                 0                0  
# 88     658               100                 0                0  
# 89     659               100                 0                0  
# 90     660                75                25               25  
# 91     661               100                 0                0  
# 92     662               100                 0                0  
# 93     669                66.7              33.3              0  
# 94     672                66.7               0                0  
# 95     677                50                 0                0  
# 96     678               100                40               40  
# 97     679                60                20               40  
# 98     680                80                20               40  
# 99     681                60                20               40  
# 100     682                75                 0               25  
# 101     683                57.1              42.9             57.1
# 102     684                85.7              57.1             57.1
# 103     685                60                20               40  
# 104     686                66.7              50               66.7
# 105     687                42.9              57.1             57.1
# 106     688                50                33.3             50  
# 107     689                60                20               40  
# 108     690                80                40               40  
# 109     691                66.7               0                0  
# 110     692                83.3              50               50  
# 111     693               100                 0                0  
# 112     694               100                 0                0  
# 113     695                50                50               50  
# 114     696                60                20               40  
# 115     697               100                 0                0  
# 116     698                60                20               40  
# 117     699                66.7               0                0  
# 118     700               100                 0                0  
# 119     701                57.1              42.9             57.1
# 120     702               100                 0                0  
# 121     703                80                40               40  
# 122     704               100                 0                0  
# 123     705               100                 0                0  
# 124     706                60                20               40  
# 125     707               100                 0                0  
# 126     708                66.7              33.3             50  
# 127     709                42.9              42.9             42.9
# 128     710               100                 0                0  
# 129     711                75                 0               25  
# 130     714                80                20               20  
# 131     715               100                 0                0  
# 132     716                66.7              16.7             33.3
# 133     717               100                 0                0  
# 134     718                42.9              28.6             28.6
# 135     720                50                25                0  
# 136     723                33.3              16.7             16.7
# 137     724                50                16.7             16.7
# 138     726                28.6              28.6             28.6
# 139     727                55.6              44.4             66.7
# 140     728                50                37.5             62.5
# 141     729                20                40               50  
# 142     730                22.2              44.4             55.6
# 143     731                44.4              11.1             55.6
# 144     732                37.5               0               37.5
# 145     733                25                37.5             37.5
# 146     734                25                50               50  
# 147     735                25                12.5             37.5
# 148     736                25                37.5             37.5
# 149     737                57.1               0               28.6
# 150     738                33.3              33.3             55.6
# 151     739                 0                50               50  
# 152     740                11.1              11.1             55.6
# 153     741                33.3              11.1             33.3
# 154     742                25                37.5             37.5
# 155     743                22.2              44.4             44.4
# 156     744                37.5              50               50  
# 157     745                22.2              44.4             44.4
# 158     746                42.9              42.9             57.1
# 159     747                50                25               37.5
# 160     748                50                 0               25  
# 161     749                37.5              12.5             37.5
# 162     750                12.5              50               62.5
# 163     751                22.2              22.2             33.3
# 164     752                22.2              22.2             33.3
# 165     753                42.9              28.6             28.6
# 166     754                50                33.3             33.3
# 167     755                33.3              33.3             33.3
# 168     756                22.2              33.3             33.3
# 169     757                11.1              22.2             33.3
# 170     758                25                12.5             25  
# 171     759                57.1              28.6             28.6
# 172     760                40                40               40  
# 173     761                37.5              62.5             50  
# 174     762                22.2              33.3             33.3
# 175     763                22.2              33.3             33.3
# 176     764                28.6               0               28.6
# 177     765                28.6              42.9             42.9
# 178     766                42.9              28.6             42.9
# 179     767                22.2              33.3             44.4
# 180     768                22.2               0               33.3
# 181     769                22.2              44.4             55.6
# 182     770                22.2              44.4             44.4
# 183     771                22.2              33.3             44.4
# 184     772                22.2              22.2             33.3
# 185     773                22.2              33.3             33.3
# 186     774                22.2               0               33.3
# 187     775                22.2              44.4             44.4
# 188     776                 0                30               30  
# 189     777                25                37.5             50  
# 190     778                22.2              44.4             44.4
# 191     779                 0                37.5             50  
# 192     780                11.1              44.4             44.4
# 193     781                11.1              44.4             44.4
# 194     782                22.2              44.4             44.4
# 195     783                22.2              44.4             44.4
# 196     784                22.2              44.4             44.4
# 197     785                22.2              44.4             44.4
# 198     786                25                50               50  
# 199     787                25                37.5             50  
# 200     788                37.5               0               37.5
# 201     789                 0                40               50  
# 202     790                25                25               37.5
# 203     791                25                25               37.5
# 204     792                22.2              44.4             44.4
# 205     793                25                37.5             50  
# 206     794                22.2              22.2             33.3
# 207     795                22.2              44.4             44.4
# 208     796                22.2              44.4             44.4
# 209     797                25                50               50  
# 210     798                11.1              33.3             44.4
# 211     799                25                37.5             62.5
# 212     800                22.2              44.4             44.4
# 213     801                28.6              42.9             57.1
# 214     802                25                50               50  
# 215     803                20                40               40  
# 216     804                25                25               50  
# 217     805                 0                10               30  
# 218     806                 0                11.1             44.4
# 219     807                22.2              44.4             44.4
# 220     808                20                40               40  
# 221     809                 0                30               50  
# 222     810                22.2              44.4             44.4
# 223     811                22.2              55.6             44.4
# 224     812                25                25               50  
# 225     813                22.2              44.4             44.4
# 226     814                 0                25               50  
# 227     815                 0                30               40  
# 228     816                 0                44.4             44.4
# 229     817                20                40               40  
# 230     818                25                50               50  
# 231     819                22.2              44.4             44.4
# 232     820                20                40               40  
# 233     821                 0                40               50  
# 234     822                 0                44.4             44.4
# 235     823                25                50               50  
# 236     824                 0                40               30  
# 237     825                 0                40               40  
# 238     826                 0                44.4             44.4
# 239     827                 0                40               40  
# 240     828                25                50               50  
# 241     829                11.1              33.3             44.4
# 242     830                 0                40               50  
# 243     831                 0                40               40  
# 244     832                 0                40               50  
# 245     833                 0                30               40  
# 246     834                 0                40               40  
# 247     835                22.2              44.4             44.4
# 248     836                 0                75               50  
# 249     837                 0                44.4             44.4
# 250     838                 0                40               40  
# 251     839                83.3              33.3             33.3
# 252     840               100                 0                0  
# 253     841                28.6              28.6             28.6
# 254     843                71.4              14.3             42.9
# 255     844                42.9              28.6             28.6
# 256     846                42.9              14.3             42.9
# 257     847                42.9              28.6             28.6
# 258     848                40                 0                0  
# 259     849                42.9              28.6             28.6
# 260     850                50                37.5             37.5
# 261     851                42.9              42.9             42.9
# 262     852                42.9              14.3             42.9
# 263     853               100               100                0  
# 264     854               100                 0                0  
# 265     855                42.9              28.6             28.6
# 266     950               100               100                0  
# 267     951               100               100                0  
# 268     955               100               100                0  
# 269     956               100               100                0  
# 270     959               100               100                0  
# 271     960               100               100               50  
# 272     962               100                 0                0  
# 273     969               100                 0                0  
# 274     971                 0                 0                0  
# 275     972                66.7              33.3              0  
# 276     976                 0                 0                0  
# 277     978                 0                33.3              0  
# 278     979                 0                33.3             33.3
# 279     982               100                33.3              0  
# 280     985                 0               100              100  
# 281     988                33.3             100               33.3
# 282     990                33.3             100               33.3
# 283     991                66.7             100               33.3
# 284     992                50               100                0  
# 285     993                 0               100              100  
# 286     995                 0                25               50  
# 287     997                33.3              33.3             33.3
# 288     999                25                25               50  
# 289    1001                 0                25               25  
# 290    1002                 0               100              100  
# 291    1003                 0               100              100  
# 292    1004                75                50               50  
# 293    1005                 0                25               25  
# 294    1006                25                25               25  
# 295    1007                75                50               25  
# 296    1008                 0                50               50  
# 297    1009                25                50               25  
# 298    1010                 0                50               25  
# 299    1011                66.7             100               33.3
# 300    1012                 0                33.3             33.3
# 301    1013                 0               100              100  
# 302    1015                33.3             100               33.3
# 303    1016                50                25              100  
# 304    1017                40                60              100  
# 305    1018                 0                40              100  
# 306    1019                 0                25              100  

# exact match (not lenient):
#   Ecds_ID native_fallback_pct tame_fallback_pct not_fallback_pct
# 1     375                40                80               80  
# 2     377               100                 0              100  
# 3     378               100                 0              100  
# 4     379                 0                 0                0  
# 5     380                33.3              66.7             66.7
# 6     381               100                 0              100  
# 7     382               100                 0              100  
# 8     383               100                 0              100  
# 9     391               100                 0              100  
# 10     400               100                 0              100  
# 11     401               100                 0              100  
# 12     402               100                 0              100  
# 13     408               100                 0              100  
# 14     409               100                 0              100  
# 15     410               100                 0              100  
# 16     411               100                 0              100  
# 17     412               100                 0              100  
# 18     413               100                 0              100  
# 19     541               100                 0              100  
# 20     542               100                 0              100  
# 21     543               100                 0              100  
# 22     544               100                 0              100  
# 23     545               100                 0              100  
# 24     546               100                 0              100  
# 25     547               100                 0              100  
# 26     549               100                 0              100  
# 27     550               100                 0              100  
# 28     551               100                 0              100  
# 29     552               100                 0              100  
# 30     553               100                 0              100  
# 31     554               100                 0              100  
# 32     555               100                 0              100  
# 33     556               100                 0              100  
# 34     557               100                 0              100  
# 35     558               100                 0              100  
# 36     560               100                 0              100  
# 37     561               100                 0              100  
# 38     562               100                 0              100  
# 39     564               100                 0              100  
# 40     565               100                 0              100  
# 41     566               100                 0              100  
# 42     567               100                 0              100  
# 43     568               100                 0              100  
# 44     569               100                 0              100  
# 45     570               100                 0              100  
# 46     571               100                 0              100  
# 47     572               100                 0              100  
# 48     586                 0                 0              100  
# 49     587                 0                 0              100  
# 50     590               100               100              100  
# 51     591                50               100               50  
# 52     592                50               100              100  
# 53     593                 0               100              100  
# 54     594                 0                 0              100  
# 55     595                 0                 0              100  
# 56     596                 0                50              100  
# 57     597                 0                 0              100  
# 58     599                 0                50              100  
# 59     600                 0                50              100  
# 60     610                 0               100              100  
# 61     611                 0               100              100  
# 62     615                50                50              100  
# 63     616                50                50              100  
# 64     617               100                50               50  
# 65     618                 0               100                0  
# 66     619                 0                 0              100  
# 67     622                50                50              100  
# 68     623               100                66.7            100  
# 69     628                50                50              100  
# 70     629               100                50              100  
# 71     630               100                50              100  
# 72     631                40               100               80  
# 73     650               100                 0              100  
# 74     660                50                50              100  
# 75     669                 0                 0                0  
# 76     672                 0               100              100  
# 77     677                 0                 0                0  
# 78     678               100                66.7            100  
# 79     679                33.3              33.3            100  
# 80     680                66.7              33.3            100  
# 81     681                 0                33.3            100  
# 82     682                 0                 0              100  
# 83     683                20                60              100  
# 84     684                80                80              100  
# 85     685                 0                33.3            100  
# 86     686                50                75              100  
# 87     687                20                80              100  
# 88     688                25                50              100  
# 89     689                 0                33.3            100  
# 90     690                66.7             100              100  
# 91     692                75                75              100  
# 92     695                25               100              100  
# 93     696                 0                33.3            100  
# 94     698                 0                33.3            100  
# 95     701                40                60              100  
# 96     703                66.7              66.7            100  
# 97     706                 0                33.3            100  
# 98     708                50                50              100  
# 99     709                20                80               80  
# 100     710               100                 0              100  
# 101     711                50                 0              100  
# 102     714                66.7              33.3             66.7
# 103     715               100                 0                0  
# 104     716                33.3              33.3             66.7
# 105     717               100                 0                0  
# 106     718                 0                60               60  
# 107     720                 0                 0                0  
# 108     723                 0                50               50  
# 109     724                 0                33.3             33.3
# 110     726                 0                60               60  
# 111     727                33.3              66.7            100  
# 112     728                33.3              50              100  
# 113     729                 0                50               75  
# 114     730                 0                66.7             88.9
# 115     731                28.6              14.3             85.7
# 116     732                 0                 0               50  
# 117     733                 0                66.7             66.7
# 118     734                 0                83.3             83.3
# 119     735                 0                16.7             50  
# 120     736                 0                66.7             66.7
# 121     737                20                 0               40  
# 122     738                14.3              42.9             85.7
# 123     739                 0                70               80  
# 124     740                11.1              33.3             88.9
# 125     741                 0                33.3             55.6
# 126     742                 0                66.7             66.7
# 127     743                 0                71.4             71.4
# 128     744                16.7              66.7             83.3
# 129     745                 0                77.8             66.7
# 130     746                20                60              100  
# 131     747                33.3              33.3             66.7
# 132     748                33.3               0               33.3
# 133     749                 0                16.7             50  
# 134     750                 0                75              100  
# 135     751                 0                28.6             57.1
# 136     752                 0                28.6             57.1
# 137     753                 0                40               40  
# 138     754                25                50               75  
# 139     755                14.3              42.9             42.9
# 140     756                 0                42.9             42.9
# 141     757                 0                55.6             66.7
# 142     758                 0                50               62.5
# 143     759                25                25               50  
# 144     760                25                62.5             62.5
# 145     761                16.7              83.3             66.7
# 146     762                 0                42.9             42.9
# 147     763                 0                57.1             57.1
# 148     764                 0                 0               40  
# 149     765                 0                60               60  
# 150     766                20                40               60  
# 151     767                 0                42.9             57.1
# 152     768                 0                 0               42.9
# 153     769                 0                57.1             85.7
# 154     770                 0                57.1             57.1
# 155     771                 0                42.9             71.4
# 156     772                11.1              44.4             55.6
# 157     773                 0                55.6             55.6
# 158     774                 0                 0               42.9
# 159     775                 0                57.1             57.1
# 160     776                 0                50               50  
# 161     777                 0                50               83.3
# 162     778                 0                71.4             71.4
# 163     779                 0                62.5             87.5
# 164     780                 0                77.8             77.8
# 165     781                 0                66.7             77.8
# 166     782                 0                66.7             66.7
# 167     783                 0                57.1             57.1
# 168     784                 0                71.4             71.4
# 169     785                 0                71.4             57.1
# 170     786                 0                83.3             83.3
# 171     787                 0                66.7             83.3
# 172     788                12.5              25               62.5
# 173     789                 0                60               70  
# 174     790                 0                 0               50  
# 175     791                 0                33.3             66.7
# 176     792                 0                57.1             57.1
# 177     793                 0                50               83.3
# 178     794                 0                28.6             42.9
# 179     795                 0                71.4             57.1
# 180     796                 0                57.1             57.1
# 181     797                 0                83.3             83.3
# 182     798                 0                55.6             77.8
# 183     799                 0                66.7            100  
# 184     800                 0                77.8             77.8
# 185     801                 0                80              100  
# 186     802                 0                83.3             83.3
# 187     803                 0                75               75  
# 188     804                 0                33.3             83.3
# 189     805                 0                40               70  
# 190     806                 0                44.4             77.8
# 191     807                 0                85.7             85.7
# 192     808                 0                80               80  
# 193     809                 0                70               90  
# 194     810                 0                85.7             85.7
# 195     811                 0                85.7             85.7
# 196     812                 0                33.3             83.3
# 197     813                 0                85.7             85.7
# 198     814                 0                50               87.5
# 199     815                 0                70               80  
# 200     816                 0                88.9             88.9
# 201     817                 0                75               75  
# 202     818                 0                83.3             83.3
# 203     819                 0                85.7             85.7
# 204     820                 0                75               75  
# 205     821                 0                80               90  
# 206     822                 0                88.9             88.9
# 207     823                 0                83.3             83.3
# 208     824                 0                80               70  
# 209     825                 0                80               80  
# 210     826                 0                88.9             88.9
# 211     827                 0                80               80  
# 212     828                 0                71.4             85.7
# 213     829                 0                66.7             77.8
# 214     830                 0                80               90  
# 215     831                 0                80               80  
# 216     832                 0                80               90  
# 217     833                 0                70               80  
# 218     834                 0                80               80  
# 219     835                 0                85.7             85.7
# 220     836                 0                87.5             87.5
# 221     837                 0                77.8             77.8
# 222     838                 0                70               70  
# 223     839                66.7              66.7             66.7
# 224     840               100                 0                0  
# 225     841                 0                60               60  
# 226     843                40                40               80  
# 227     844                20                60               60  
# 228     846                20                40               60  
# 229     847                20                60               60  
# 230     848                 0                 0                0  
# 231     849                 0                60               60  
# 232     850                50                75               75  
# 233     851                20                40               80  
# 234     852                20                40               80  
# 235     853               100               100                0  
# 236     854               100                 0                0  
# 237     855                 0                 0               40  
# 238     955               100               100                0  
# 239     956               100               100                0  
# 240     960               100               100              100  
# 241     972                 0               100                0  
# 242     978                 0               100               66.7
# 243     979                 0               100               66.7
# 244     982               100               100               66.7
# 245     985                 0               100              100  
# 246     988                 0               100              100  
# 247     990                 0               100              100  
# 248     991                 0               100              100  
# 249     993                 0               100              100  
# 250     995                 0               100              100  
# 251     997                33.3             100              100  
# 252     999                25               100              100  
# 253    1001                 0                66.7            100  
# 254    1002                 0               100              100  
# 255    1003                 0               100              100  
# 256    1004                33.3              66.7            100  
# 257    1005                 0                75               75  
# 258    1006                25                75               75  
# 259    1007                25                75               75  
# 260    1008                 0                75               75  
# 261    1009                25                75               75  
# 262    1010                 0                75               75  
# 263    1011                 0               100              100  
# 264    1012                 0               100              100  
# 265    1013                 0               100              100  
# 266    1015                 0               100              100  
# 267    1016                 0               100              100  
# 268    1017                40               100              100  
# 269    1018                 0               100              100  
# 270    1019                 0               100              100 

# Ecodistricts with poor class representation
accuracy_with_fallback_grouped %>%
  group_by(Ecds_ID) %>%
  summarise(
    mean_classes =
      mean(Truth_Classes_Present)) %>%
  arrange(mean_classes) %>%
  print(n = Inf)

# with ACI2025:
#   Ecds_ID mean_classes
# 1     400         1   
# 2     401         1   
# 3     402         1   
# 4     589         1   
# 5     601         1   
# 6     607         1   
# 7     950         1   
# 8     951         1   
# 9     955         1   
# 10     956         1   
# 11     959         1   
# 12     985         1   
# 13    1002         1   
# 14    1003         1   
# 15    1011         1   
# 16     378         1.33
# 17     381         1.33
# 18     383         1.33
# 19     408         1.33
# 20     409         1.33
# 21     410         1.33
# 22     411         1.33
# 23     412         1.33
# 24     413         1.33
# 25     541         1.33
# 26     542         1.33
# 27     543         1.33
# 28     544         1.33
# 29     545         1.33
# 30     546         1.33
# 31     547         1.33
# 32     549         1.33
# 33     550         1.33
# 34     551         1.33
# 35     552         1.33
# 36     553         1.33
# 37     554         1.33
# 38     555         1.33
# 39     556         1.33
# 40     557         1.33
# 41     558         1.33
# 42     560         1.33
# 43     561         1.33
# 44     562         1.33
# 45     564         1.33
# 46     565         1.33
# 47     566         1.33
# 48     567         1.33
# 49     568         1.33
# 50     569         1.33
# 51     570         1.33
# 52     571         1.33
# 53     572         1.33
# 54     587         1.33
# 55     611         1.33
# 56     374         1.5 
# 57     382         1.5 
# 58     391         1.5 
# 59     559         1.5 
# 60     960         1.5 
# 61    1016         1.5 
# 62    1017         1.6 
# 63     377         1.67
# 64     988         1.67
# 65     990         1.67
# 66     991         1.67
# 67     591         1.75
# 68     628         1.75
# 69     734         1.75
# 70     761         1.75
# 71     786         1.75
# 72     799         1.75
# 73     802         1.75
# 74     823         1.75
# 75     836         1.75
# 76    1004         1.75
# 77    1019         1.75
# 78    1018         1.8 
# 79     695         1.83
# 80     839         1.83
# 81     801         1.86
# 82     818         1.88
# 83     770         1.89
# 84     775         1.89
# 85     778         1.89
# 86     783         1.89
# 87     784         1.89
# 88     785         1.89
# 89     792         1.89
# 90     795         1.89
# 91     796         1.89
# 92     807         1.89
# 93     810         1.89
# 94     813         1.89
# 95     819         1.89
# 96     835         1.89
# 97     375         2   
# 98     583         2   
# 99     585         2   
# 100     586         2   
# 101     588         2   
# 102     590         2   
# 103     592         2   
# 104     593         2   
# 105     594         2   
# 106     595         2   
# 107     596         2   
# 108     597         2   
# 109     598         2   
# 110     599         2   
# 111     600         2   
# 112     610         2   
# 113     612         2   
# 114     617         2   
# 115     618         2   
# 116     622         2   
# 117     623         2   
# 118     629         2   
# 119     631         2   
# 120     658         2   
# 121     659         2   
# 122     661         2   
# 123     662         2   
# 124     684         2   
# 125     686         2   
# 126     687         2   
# 127     688         2   
# 128     690         2   
# 129     693         2   
# 130     694         2   
# 131     697         2   
# 132     700         2   
# 133     701         2   
# 134     702         2   
# 135     704         2   
# 136     705         2   
# 137     707         2   
# 138     709         2   
# 139     715         2   
# 140     717         2   
# 141     728         2   
# 142     729         2   
# 143     733         2   
# 144     736         2   
# 145     738         2   
# 146     742         2   
# 147     746         2   
# 148     754         2   
# 149     760         2   
# 150     765         2   
# 151     766         2   
# 152     767         2   
# 153     787         2   
# 154     797         2   
# 155     803         2   
# 156     804         2   
# 157     811         2   
# 158     817         2   
# 159     820         2   
# 160     840         2   
# 161     853         2   
# 162     854         2   
# 163     972         2   
# 164     992         2   
# 165     993         2   
# 166    1007         2   
# 167    1008         2   
# 168    1009         2   
# 169    1013         2   
# 170    1015         2   
# 171     830         2.1 
# 172     832         2.1 
# 173     727         2.11
# 174     731         2.11
# 175     743         2.11
# 176     755         2.11
# 177     756         2.11
# 178     762         2.11
# 179     763         2.11
# 180     769         2.11
# 181     780         2.11
# 182     782         2.11
# 183     816         2.11
# 184     822         2.11
# 185     826         2.11
# 186     744         2.12
# 187     750         2.12
# 188     777         2.12
# 189     790         2.12
# 190     793         2.12
# 191     812         2.12
# 192     683         2.14
# 193     718         2.14
# 194     726         2.14
# 195     753         2.14
# 196     841         2.14
# 197     843         2.14
# 198     847         2.14
# 199     852         2.14
# 200     692         2.17
# 201     708         2.17
# 202     716         2.17
# 203     724         2.17
# 204     380         2.2 
# 205     678         2.2 
# 206     679         2.2 
# 207     680         2.2 
# 208     681         2.2 
# 209     685         2.2 
# 210     689         2.2 
# 211     696         2.2 
# 212     698         2.2 
# 213     703         2.2 
# 214     706         2.2 
# 215     739         2.2 
# 216     789         2.2 
# 217     808         2.2 
# 218     821         2.2 
# 219     825         2.2 
# 220     827         2.2 
# 221     831         2.2 
# 222     834         2.2 
# 223     730         2.22
# 224     745         2.22
# 225     751         2.22
# 226     752         2.22
# 227     771         2.22
# 228     794         2.22
# 229     829         2.22
# 230     615         2.25
# 231     616         2.25
# 232     630         2.25
# 233     660         2.25
# 234     682         2.25
# 235     711         2.25
# 236     735         2.25
# 237     747         2.25
# 238     749         2.25
# 239     779         2.25
# 240     791         2.25
# 241     814         2.25
# 242     828         2.25
# 243     999         2.25
# 244    1010         2.25
# 245     759         2.29
# 246     844         2.29
# 247     846         2.29
# 248     849         2.29
# 249     851         2.29
# 250     824         2.3 
# 251     647         2.33
# 252     650         2.33
# 253     652         2.33
# 254     657         2.33
# 255     669         2.33
# 256     672         2.33
# 257     691         2.33
# 258     699         2.33
# 259     710         2.33
# 260     723         2.33
# 261     740         2.33
# 262     773         2.33
# 263     781         2.33
# 264     837         2.33
# 265     982         2.33
# 266     997         2.33
# 267     732         2.38
# 268     850         2.38
# 269     714         2.4 
# 270     776         2.4 
# 271     809         2.4 
# 272     838         2.4 
# 273     737         2.43
# 274     764         2.43
# 275     855         2.43
# 276     757         2.44
# 277     768         2.44
# 278     772         2.44
# 279     774         2.44
# 280     798         2.44
# 281     800         2.44
# 282     806         2.44
# 283     379         2.5 
# 284     677         2.5 
# 285     720         2.5 
# 286     748         2.5 
# 287     815         2.5 
# 288     833         2.5 
# 289     995         2.5 
# 290    1001         2.5 
# 291    1005         2.5 
# 292    1006         2.5 
# 293     741         2.56
# 294     805         2.6 
# 295     848         2.6 
# 296     758         2.62
# 297     788         2.62
# 298     978         2.67
# 299     979         2.67
# 300    1012         2.67
# 301     619         3   
# 302     626         3   
# 303     962         3   
# 304     969         3   
# 305     971         3   
# 306     976         3 

# exact match (not lenient):
#   Ecds_ID mean_classes
# 1     377         1   
# 2     378         1   
# 3     381         1   
# 4     382         1   
# 5     383         1   
# 6     391         1   
# 7     400         1   
# 8     401         1   
# 9     402         1   
# 10     408         1   
# 11     409         1   
# 12     410         1   
# 13     411         1   
# 14     412         1   
# 15     413         1   
# 16     541         1   
# 17     542         1   
# 18     543         1   
# 19     544         1   
# 20     545         1   
# 21     546         1   
# 22     547         1   
# 23     549         1   
# 24     550         1   
# 25     551         1   
# 26     552         1   
# 27     553         1   
# 28     554         1   
# 29     555         1   
# 30     556         1   
# 31     557         1   
# 32     558         1   
# 33     560         1   
# 34     561         1   
# 35     562         1   
# 36     564         1   
# 37     565         1   
# 38     566         1   
# 39     567         1   
# 40     568         1   
# 41     569         1   
# 42     570         1   
# 43     571         1   
# 44     572         1   
# 45     672         1   
# 46     955         1   
# 47     956         1   
# 48     960         1   
# 49     982         1   
# 50     985         1   
# 51     988         1   
# 52     990         1   
# 53     991         1   
# 54     997         1   
# 55    1002         1   
# 56    1003         1   
# 57    1011         1   
# 58    1017         1.2 
# 59     816         1.22
# 60     822         1.22
# 61     826         1.22
# 62     695         1.25
# 63     836         1.25
# 64     999         1.25
# 65    1009         1.25
# 66    1019         1.25
# 67     807         1.29
# 68     810         1.29
# 69     813         1.29
# 70     819         1.29
# 71     835         1.29
# 72     830         1.3 
# 73     832         1.3 
# 74     690         1.33
# 75     786         1.33
# 76     802         1.33
# 77     823         1.33
# 78     978         1.33
# 79     979         1.33
# 80    1001         1.33
# 81    1004         1.33
# 82    1015         1.33
# 83     801         1.4 
# 84     808         1.4 
# 85     821         1.4 
# 86     825         1.4 
# 87     827         1.4 
# 88     831         1.4 
# 89     834         1.4 
# 90    1018         1.4 
# 91     811         1.43
# 92     591         1.5 
# 93     617         1.5 
# 94     734         1.5 
# 95     750         1.5 
# 96     799         1.5 
# 97     803         1.5 
# 98     817         1.5 
# 99     818         1.5 
# 100     820         1.5 
# 101     824         1.5 
# 102    1005         1.5 
# 103    1006         1.5 
# 104    1007         1.5 
# 105    1008         1.5 
# 106    1010         1.5 
# 107    1016         1.5 
# 108     727         1.56
# 109     780         1.56
# 110     829         1.56
# 111     778         1.57
# 112     784         1.57
# 113     375         1.6 
# 114     631         1.6 
# 115     709         1.6 
# 116     779         1.62
# 117     814         1.62
# 118     850         1.62
# 119     380         1.67
# 120     623         1.67
# 121     733         1.67
# 122     742         1.67
# 123     745         1.67
# 124     761         1.67
# 125     782         1.67
# 126     787         1.67
# 127     797         1.67
# 128     837         1.67
# 129     839         1.67
# 130     739         1.7 
# 131     809         1.7 
# 132     833         1.7 
# 133     785         1.71
# 134     688         1.75
# 135     730         1.78
# 136     757         1.78
# 137     781         1.78
# 138     800         1.78
# 139     684         1.8 
# 140     687         1.8 
# 141     701         1.8 
# 142     718         1.8 
# 143     726         1.8 
# 144     746         1.8 
# 145     789         1.8 
# 146     815         1.8 
# 147     838         1.8 
# 148     841         1.8 
# 149     843         1.8 
# 150     847         1.8 
# 151     852         1.8 
# 152     728         1.83
# 153     736         1.83
# 154     804         1.83
# 155     738         1.86
# 156     763         1.86
# 157     770         1.86
# 158     775         1.86
# 159     783         1.86
# 160     792         1.86
# 161     795         1.86
# 162     796         1.86
# 163     828         1.86
# 164     729         1.88
# 165     758         1.88
# 166     760         1.88
# 167     740         1.89
# 168     773         1.89
# 169     806         1.89
# 170     805         1.9 
# 171     586         2   
# 172     587         2   
# 173     590         2   
# 174     592         2   
# 175     593         2   
# 176     594         2   
# 177     595         2   
# 178     596         2   
# 179     597         2   
# 180     599         2   
# 181     600         2   
# 182     610         2   
# 183     611         2   
# 184     615         2   
# 185     616         2   
# 186     618         2   
# 187     619         2   
# 188     622         2   
# 189     628         2   
# 190     629         2   
# 191     630         2   
# 192     650         2   
# 193     660         2   
# 194     678         2   
# 195     679         2   
# 196     680         2   
# 197     681         2   
# 198     682         2   
# 199     683         2   
# 200     685         2   
# 201     686         2   
# 202     689         2   
# 203     692         2   
# 204     696         2   
# 205     698         2   
# 206     703         2   
# 207     706         2   
# 208     708         2   
# 209     710         2   
# 210     711         2   
# 211     715         2   
# 212     717         2   
# 213     723         2   
# 214     731         2   
# 215     743         2   
# 216     744         2   
# 217     754         2   
# 218     765         2   
# 219     766         2   
# 220     767         2   
# 221     769         2   
# 222     772         2   
# 223     776         2   
# 224     777         2   
# 225     793         2   
# 226     798         2   
# 227     812         2   
# 228     840         2   
# 229     844         2   
# 230     846         2   
# 231     849         2   
# 232     851         2   
# 233     853         2   
# 234     854         2   
# 235     972         2   
# 236     993         2   
# 237     995         2   
# 238    1012         2   
# 239    1013         2   
# 240     741         2.11
# 241     788         2.12
# 242     752         2.14
# 243     755         2.14
# 244     756         2.14
# 245     762         2.14
# 246     771         2.14
# 247     791         2.17
# 248     753         2.2 
# 249     751         2.29
# 250     794         2.29
# 251     714         2.33
# 252     716         2.33
# 253     724         2.33
# 254     735         2.33
# 255     747         2.33
# 256     749         2.33
# 257     732         2.5 
# 258     759         2.5 
# 259     790         2.5 
# 260     768         2.57
# 261     774         2.57
# 262     737         2.6 
# 263     764         2.6 
# 264     855         2.6 
# 265     748         2.67
# 266     379         3   
# 267     669         3   
# 268     677         3   
# 269     720         3   
# 270     848         3  

# Support by ecodistrict 
accuracy_with_fallback_grouped %>%
  group_by(Ecds_ID) %>%
  summarise(
    support0 = sum(Support_0),
    support1 = sum(Support_1),
    support2 = sum(Support_2)
  ) %>%
  print(n = Inf)

# with ACI2025:
#   Ecds_ID support0 support1 support2
# 1     374        2       56        0
# 2     375     1031      110      105
# 3     377      619       93        0
# 4     378       94      242        0
# 5     379     1270      294       73
# 6     380      172       86      114
# 7     381      366      617        0
# 8     382       73      370        0
# 9     383      149      344        0
# 10     391        1      112        0
# 11     400        0      106        0
# 12     401        0       94        0
# 13     402        0      106        0
# 14     408      206      122        0
# 15     409       43       99        0
# 16     410      261      164        0
# 17     411      325      191        0
# 18     412      142      116        0
# 19     413       59      104        0
# 20     541       61      424        0
# 21     542      373      433        0
# 22     543      372      112        0
# 23     544      657      314        0
# 24     545      166      424        0
# 25     546      167      112        0
# 26     547      425      226        0
# 27     549       58       42        0
# 28     550      174      287        0
# 29     551      870      553        0
# 30     552      223      122        0
# 31     553      482      354        0
# 32     554       62       97        0
# 33     555      313      272        0
# 34     556      335      466        0
# 35     557      384      499        0
# 36     558      438      698        0
# 37     559     3607     1799        0
# 38     560      318      258        0
# 39     561     4731     1246        0
# 40     562     2393      643        0
# 41     564      432      226        0
# 42     565     1018      235        0
# 43     566      420       46        0
# 44     567      575      254        0
# 45     568      698      130        0
# 46     569      589      193        0
# 47     570      810       82        0
# 48     571       76      629        0
# 49     572       58      630        0
# 50     583       44       48        0
# 51     585       30       46        0
# 52     586     2029      146        4
# 53     587      342      140       63
# 54     588     1327      220        0
# 55     589       80        0        0
# 56     590      792      499        1
# 57     591     2575     1379        7
# 58     592     1840      495      120
# 59     593     2862      360       19
# 60     594      705      464       63
# 61     595     2000      142       63
# 62     596     1866      545      139
# 63     597      348      311       63
# 64     598      799      184        0
# 65     599      809      241      179
# 66     600      106      382      179
# 67     601       80        0        0
# 68     607      270        0        0
# 69     610      268      338      116
# 70     611      163       59      116
# 71     612       12      181        0
# 72     615      189      359      119
# 73     616      193      206      117
# 74     617      120      237       32
# 75     618      157      172       12
# 76     619       10      226      102
# 77     622       77      206       37
# 78     623      154      409       84
# 79     626      406     1040       38
# 80     628       10      289      139
# 81     629        6      121       42
# 82     630       17      181       51
# 83     631        1      124      123
# 84     647      101      187       32
# 85     650       65      306       32
# 86     652     4230     1813       32
# 87     657      135      221       32
# 88     658     4167     1628        0
# 89     659     4167     1628        0
# 90     660      139      382      148
# 91     661      212       18        0
# 92     662      946       50        0
# 93     669      569      719      443
# 94     672      640      741      261
# 95     677      461     1454      704
# 96     678     2556     1296       17
# 97     679       61     1291      515
# 98     680     1128     1468      159
# 99     681     2263      477      511
# 100     682      331      239      427
# 101     683      365      412      578
# 102     684     1250     1863      123
# 103     685      994      316      515
# 104     686       38      674      178
# 105     687     5023     1384      237
# 106     688      667      637      538
# 107     689     1388      457      517
# 108     690     1542      433      161
# 109     691     1248     2615      412
# 110     692      513      918      172
# 111     693     2349      135        0
# 112     694     2687      277        0
# 113     695      823      332      212
# 114     696     3534      651      513
# 115     697     1928      319        0
# 116     698     1120      321      510
# 117     699     1384     2683      412
# 118     700     1597      247        0
# 119     701      695      445      581
# 120     702     3339      107        0
# 121     703      849      719      151
# 122     704     2108      211        0
# 123     705     4809      105        0
# 124     706     2609      372      513
# 125     707     4691      424        0
# 126     708      485      452      565
# 127     709     8093      412      176
# 128     710      347      263       32
# 129     711      343      502      509
# 130     714     1471      342      136
# 131     715      648      356        0
# 132     716      636      414      600
# 133     717     3758      567        0
# 134     718      273     1139     1043
# 135     720      407     1424      703
# 136     723      425     1252      668
# 137     724     2551      678      426
# 138     726      862     1199      900
# 139     727     5880      563      671
# 140     728     3971      284      635
# 141     729     9892     1995     2467
# 142     730     6762     2662     2067
# 143     731     6973     1129     1478
# 144     732     1434     2453    15786
# 145     733     3661      408      259
# 146     734      944      267     1141
# 147     735     3600      552    13139
# 148     736     5447      330      494
# 149     737     3309     2422    10087
# 150     738     2914      253      669
# 151     739      794      983     2305
# 152     740     2092     1100     1202
# 153     741     3967      925    14021
# 154     742      723      151      136
# 155     743     1311      217     1021
# 156     744     3210      908     1052
# 157     745    15477     1109     1696
# 158     746     1913      145      139
# 159     747     1178      704    13120
# 160     748    11590     2008    13296
# 161     749     7380      835    13121
# 162     750      313      475      949
# 163     751      845      634     5021
# 164     752     9818     1180     3060
# 165     753     8719      387      341
# 166     754     3924      228      181
# 167     755     2916      231     2911
# 168     756     6027     1206     3905
# 169     757     1174      484     3960
# 170     758     2802      347     2781
# 171     759     1796     2226    13224
# 172     760     8375     1287     1097
# 173     761      329      106      289
# 174     762     1365      438     3823
# 175     763     3834     1899     3337
# 176     764     2466     2266    13234
# 177     765     8449      323      288
# 178     766     3149      253      157
# 179     767    10534      481      333
# 180     768      786      160    20611
# 181     769     1840     1446     2425
# 182     770     5380     1606     2751
# 183     771     1940     1476     4160
# 184     772     1100      163     7237
# 185     773     6875      889      827
# 186     774     1405      248    20587
# 187     775     1659      384      788
# 188     776      835      959     2505
# 189     777      280      608     1414
# 190     778      374      196      154
# 191     779      422      680      810
# 192     780     6626      202      250
# 193     781     3153      307      471
# 194     782     9598     1247     1123
# 195     783     5163      655      501
# 196     784     2238     1079     2063
# 197     785     1174      364      210
# 198     786      978      242     1117
# 199     787     1734      648     1790
# 200     788     1005      164    20581
# 201     789     2700      802     1449
# 202     790      441      203    20580
# 203     791     1963      326     7236
# 204     792    10877      740      665
# 205     793     3989     1270      710
# 206     794     8340      945     7745
# 207     795     4224     1789     2979
# 208     796      760      554     2749
# 209     797      184      542     1106
# 210     798     4797      594      404
# 211     799      201      729     1477
# 212     800     1356      603     1043
# 213     801      191      285      728
# 214     802      311      359      680
# 215     803     2939     1820     1810
# 216     804      638      952     1632
# 217     805     1472      563    71369
# 218     806      747     3616     9088
# 219     807     1724      696     1334
# 220     808    10169     1507     5330
# 221     809     2446     2560     7904
# 222     810     1418      454      658
# 223     811      530      167      713
# 224     812      556      852     2211
# 225     813     1935      770     4027
# 226     814      128      846     2422
# 227     815      901     1329     4457
# 228     816     4837     1821     1685
# 229     817     2050      617     2160
# 230     818      378      608     1217
# 231     819      136      339     7609
# 232     820     2722      677     1824
# 233     821     2031      821     2303
# 234     822     1788      724     1669
# 235     823     1700      376      920
# 236     824      817      976      954
# 237     825     7161     1134      624
# 238     826    10160     2509     2329
# 239     827      406     1859     6159
# 240     828     6018     1912     6751
# 241     829      359      284      753
# 242     830     4237     2101     4919
# 243     831     2234      443      412
# 244     832     2885     2610    17009
# 245     833      706     2710    12700
# 246     834     2923     1678     5384
# 247     835      407      307     1561
# 248     836       88        3      320
# 249     837      796      491     2327
# 250     838      821     1876     9840
# 251     839     2629       42       60
# 252     840     1516      197        0
# 253     841     2239      685      311
# 254     843      722      262     1542
# 255     844     1189      329     1549
# 256     846      341     1107     1957
# 257     847     1490      760      439
# 258     848     2054     1956     1496
# 259     849    11449      816      444
# 260     850     2325      544      393
# 261     851     1097       21      133
# 262     852     2595      144     1566
# 263     853      964        4        0
# 264     854     1431       46        0
# 265     855     1400     2258    13231
# 266     950       96        0        0
# 267     951      106        0        0
# 268     955      322        0        0
# 269     956      318        0        0
# 270     959     1002        0        0
# 271     960      106        0        6
# 272     962      185      922        4
# 273     969      185      922        4
# 274     971      383     1842      254
# 275     972      281     1465       63
# 276     976      383     1840      254
# 277     978       91      239      101
# 278     979       19      136      270
# 279     982      187      891        4
# 280     985        0        0       32
# 281     988      156        0       45
# 282     990      156        0       45
# 283     991      130        0       40
# 284     992      156        0       13
# 285     993        0        8      482
# 286     995      288      457      825
# 287     997      288      320      113
# 288     999      288      449      364
# 289    1001     1693     2435      633
# 290    1002        0        0      172
# 291    1003        0        0      172
# 292    1004       92      174      939
# 293    1005      362      644      760
# 294    1006      154      789      406
# 295    1007      472      428       37
# 296    1008      118      134      139
# 297    1009      138      182       96
# 298    1010      248       74       69
# 299    1011      105        0       46
# 300    1012      288      328      574
# 301    1013        0        8      482
# 302    1015      110        8      488
# 303    1016        0       94      441
# 304    1017        0      322      172
# 305    1018        0      369      499
# 306    1019        0      449      318

# exact match (not lenient):
#   Ecds_ID support0 support1 support2
# 1     375       24       29      105
# 2     377        0       29        0
# 3     378        0      121        0
# 4     379       24       37       41
# 5     380       24       37      114
# 6     381        0      303        0
# 7     382        0      184        0
# 8     383        0      171        0
# 9     391        0       56        0
# 10     400        0       53        0
# 11     401        0       47        0
# 12     402        0       53        0
# 13     408        0       61        0
# 14     409        0       82        0
# 15     410        0       82        0
# 16     411        0       94        0
# 17     412        0       58        0
# 18     413        0       52        0
# 19     541        0      414        0
# 20     542        0      414        0
# 21     543        0       55        0
# 22     544        0      157        0
# 23     545        0      414        0
# 24     546        0       56        0
# 25     547        0      113        0
# 26     549        0       21        0
# 27     550        0      142        0
# 28     551        0      276        0
# 29     552        0       61        0
# 30     553        0      177        0
# 31     554        0       48        0
# 32     555        0      136        0
# 33     556        0      233        0
# 34     557        0      249        0
# 35     558        0      348        0
# 36     560        0      129        0
# 37     561        0      622        0
# 38     562        0      622        0
# 39     564        0      113        0
# 40     565        0      117        0
# 41     566        0       23        0
# 42     567        0      127        0
# 43     568        0       65        0
# 44     569        0       96        0
# 45     570        0       41        0
# 46     571        0      622        0
# 47     572        0      622        0
# 48     586        0       18        4
# 49     587        0      140       63
# 50     590        0       20        1
# 51     591        3       15        7
# 52     592        0       95      120
# 53     593        0        8       19
# 54     594        0      140       63
# 55     595        0      140       63
# 56     596        0       84      139
# 57     597        0      140       63
# 58     599        0      199      179
# 59     600        0      199      179
# 60     610        0       59      116
# 61     611        0       59      116
# 62     615        0       69      148
# 63     616        0       69      148
# 64     617        3       10       32
# 65     618        1        0       12
# 66     619        0       53      102
# 67     622        0       56      139
# 68     623        0       19       84
# 69     628        0       56      139
# 70     629        0       19       42
# 71     630        0       15       51
# 72     631        1        4      123
# 73     650        0       10       32
# 74     660        0       69      148
# 75     669      183      715      443
# 76     672        0        0      261
# 77     677      183      715      443
# 78     678        0      111       17
# 79     679        0      185      515
# 80     680        0       68      159
# 81     681        0      374      921
# 82     682        0      327      824
# 83     683        1      377      984
# 84     684        0       36      123
# 85     685        0      374      921
# 86     686        0       54      178
# 87     687        0       85      237
# 88     688        0      176      538
# 89     689        0      374      921
# 90     690        0       47      161
# 91     692        1       96      172
# 92     695        0       47      212
# 93     696        0      374      921
# 94     698        0      374      921
# 95     701        1      374      976
# 96     703        0       81      151
# 97     706        0      374      921
# 98     708        1      175      565
# 99     709       68       70      176
# 100     710        0       10       32
# 101     711        0      269      509
# 102     714       34      125      155
# 103     715       12       36        0
# 104     716       19      337      600
# 105     717       93       56        0
# 106     718      238     1063     1287
# 107     720      175      707      443
# 108     723       67      637      668
# 109     724      134      442      426
# 110     726       65      600      900
# 111     727        0      109      720
# 112     728        0       88      635
# 113     729      395      573     2467
# 114     730        8      174     2118
# 115     731       12      196     1478
# 116     732     1437     2867    18861
# 117     733       85       93      259
# 118     734       73       33     1141
# 119     735      282     1019    16205
# 120     736       86       63      494
# 121     737     1437     2586    13157
# 122     738       22       84      669
# 123     739       40      151     2350
# 124     740        8       85     1255
# 125     741      556     1385    17140
# 126     742       32       43      136
# 127     743       55       49     1021
# 128     744      283      440     1052
# 129     745      362      456     1749
# 130     746        0       12      139
# 131     747       67      431    13120
# 132     748      419      875    13296
# 133     749      303     1009    16196
# 134     750        0       38      988
# 135     751      138      317     5021
# 136     752      167      233     3060
# 137     753      413      655     3416
# 138     754       36       44      181
# 139     755      111       73     2911
# 140     756      390      557     3905
# 141     757      114      386     4004
# 142     758      106      162     2834
# 143     759     1231     2220    13224
# 144     760      396      593     1097
# 145     761       33       34      289
# 146     762      321      239     3823
# 147     763      120      682     3337
# 148     764     1422     2790    16309
# 149     765      410      632     3355
# 150     766      271      672     3231
# 151     767      637      140      333
# 152     768      418      555    24550
# 153     769       43      339     2425
# 154     770      381      616     2751
# 155     771       86      362     4160
# 156     772       63       84     7391
# 157     773      366      250      981
# 158     774       70      129    20587
# 159     775       80      107      788
# 160     776      128      555     2505
# 161     777        4      140     1414
# 162     778        8       48      154
# 163     779       18      130      964
# 164     780      440       58      404
# 165     781      140      137      597
# 166     782      666      562     1277
# 167     783      176      153      501
# 168     784      150      605     2063
# 169     785       72      115      210
# 170     786       36       53     1117
# 171     787       87      206     1790
# 172     788      425      564    24683
# 173     789      211      349     1571
# 174     790      409      600    24527
# 175     791      114      126     7236
# 176     792      696      312      665
# 177     793      177      348      710
# 178     794      545      416     7745
# 179     795      353     1108     2979
# 180     796      111      310     2749
# 181     797        4      232     1106
# 182     798       78      270      943
# 183     799        0       71     1477
# 184     800       81      331     1598
# 185     801        0       16      728
# 186     802       26       89      680
# 187     803       98      888     1810
# 188     804       14      219     1632
# 189     805       73      367    99815
# 190     806       32     1179     9439
# 191     807       60      250     1334
# 192     808      429      561     5683
# 193     809      100      697     7904
# 194     810       45       70      658
# 195     811       39       60      713
# 196     812       21      262     2211
# 197     813       75      172     4027
# 198     814        6      208     2779
# 199     815       44      411     4812
# 200     816      213      349     2034
# 201     817      100      221     2160
# 202     818       16      178     1217
# 203     819        6       85     7609
# 204     820      124      243     1824
# 205     821       94      147     2650
# 206     822       52       84     2026
# 207     823       67      120      920
# 208     824       31      316     1281
# 209     825      249      138      967
# 210     826      386      429     2686
# 211     827       22      569     6490
# 212     828      347      627     6930
# 213     829       28      205     1110
# 214     830      110      404     4919
# 215     831      100      178      769
# 216     832       87      503    17009
# 217     833       35     1213    13039
# 218     834      125      680     5741
# 219     835       16       63     1561
# 220     836        2        2      677
# 221     837       98      146     2882
# 222     838      191      782    10377
# 223     839       59       10       60
# 224     840       39       41        0
# 225     841       46      296      662
# 226     843       21      183     1908
# 227     844       43      164     1913
# 228     846       34      398     1957
# 229     847      127      651      802
# 230     848      836     1300     1129
# 231     849      424      653      808
# 232     850       76      227      745
# 233     851       46       15      481
# 234     852       92      100     1927
# 235     853       39        1        0
# 236     854       44       10        0
# 237     855     1396     2796    16310
# 238     955        6        0        0
# 239     956        6        0        0
# 240     960        0        0        6
# 241     972        5        0       63
# 242     978        5        0      317
# 243     979        2        0      270
# 244     982        4        0        4
# 245     985        0        0       32
# 246     988        0        0       83
# 247     990        0        0      124
# 248     991        0        0      124
# 249     993        0        8      482
# 250     995        0       11      733
# 251     997        0        0      113
# 252     999        0        3      364
# 253    1001        0     1151      591
# 254    1002        0        0      172
# 255    1003        0        0      172
# 256    1004        0       87     1021
# 257    1005      132      258      760
# 258    1006       36      392      480
# 259    1007      175      214      158
# 260    1008      175       67      162
# 261    1009       40       91      184
# 262    1010       26       37      164
# 263    1011        0        0       46
# 264    1012        0        8      482
# 265    1013        0        8      482
# 266    1015        0        8      574
# 267    1016        0        8      441
# 268    1017        0        2      172
# 269    1018        0        9      538
# 270    1019        0        3      318

# Native vs tame by ecodistrict 
accuracy_with_fallback_grouped_cons2 %>%
  group_by(Ecds_ID) %>%
  summarise(
    mean_diff =
      mean(native_metric - tame_metric,
           na.rm = TRUE)) %>%
  print(n = Inf)

# with ACI2025:
#   Ecds_ID mean_diff
# 1     374  0       
# 2     375  0.0970  
# 3     377  0.0706  
# 4     378 -0.0723  
# 5     379  0.144   
# 6     380  0.00682 
# 7     381 -0.0286  
# 8     382 -0.0852  
# 9     383 -0.0569  
# 10     391  0       
# 11     400  0       
# 12     401  0       
# 13     402  0       
# 14     408  0.0348  
# 15     409  0.0177  
# 16     410  0.0325  
# 17     411  0.0447  
# 18     412  0.0167  
# 19     413 -0.0213  
# 20     541  0.0886  
# 21     542  0.0929  
# 22     543  0.0573  
# 23     544  0.0553  
# 24     545  0.0989  
# 25     546  0.0362  
# 26     547  0.0360  
# 27     549  0.0322  
# 28     550  0.000410
# 29     551  0.0434  
# 30     552  0.0445  
# 31     553  0.0312  
# 32     554 -0.00548 
# 33     555  0.0262  
# 34     556  0.0181  
# 35     557  0.0448  
# 36     558  0.0393  
# 37     559  0.0519  
# 38     560  0.0512  
# 39     561  0.0866  
# 40     562  0.0177  
# 41     564  0.0650  
# 42     565  0.0992  
# 43     566  0.0786  
# 44     567  0.0935  
# 45     568  0.0901  
# 46     569  0.0706  
# 47     570  0.102   
# 48     571  0.103   
# 49     572  0       
# 50     583  0.137   
# 51     585  0.129   
# 52     586  0.125   
# 53     587  0.120   
# 54     588  0.229   
# 55     589  0       
# 56     590  0.128   
# 57     591  0.102   
# 58     592  0.132   
# 59     593  0.146   
# 60     594  0.256   
# 61     595  0.400   
# 62     596  0.250   
# 63     597  0.195   
# 64     598  0.220   
# 65     599  0.196   
# 66     600  0.110   
# 67     601  0       
# 68     607  0       
# 69     610  0.0609  
# 70     611  0       
# 71     612 -0.0657  
# 72     615  0.0139  
# 73     616  0.0109  
# 74     617  0.0382  
# 75     618  0.115   
# 76     619 -0.682   
# 77     622 -0.0581  
# 78     623 -0.0766  
# 79     626 -0.280   
# 80     628 -0.171   
# 81     629 -0.0480  
# 82     630 -0.173   
# 83     631 -0.0930  
# 84     647 -0.0219  
# 85     650 -0.118   
# 86     652  0.102   
# 87     657  0.0925  
# 88     658  0.195   
# 89     659  0.195   
# 90     660  0.0526  
# 91     661  0.166   
# 92     662  0.351   
# 93     669  0.0421  
# 94     672  0.0800  
# 95     677  0.0610  
# 96     678  0.0414  
# 97     679 -0.115   
# 98     680  0.00223 
# 99     681  0.0758  
# 100     682  0.0462  
# 101     683  0.0265  
# 102     684 -0.0334  
# 103     685  0.125   
# 104     686 -0.0639  
# 105     687 -0.0454  
# 106     688 -0.0811  
# 107     689  0.100   
# 108     690  0.0611  
# 109     691 -0.0501  
# 110     692 -0.0478  
# 111     693  0.224   
# 112     694  0.222   
# 113     695 -0.00740 
# 114     696  0.0795  
# 115     697  0.230   
# 116     698  0.107   
# 117     699 -0.0872  
# 118     700  0.200   
# 119     701 -0.0544  
# 120     702  0.221   
# 121     703  0.00451 
# 122     704  0.168   
# 123     705  0.262   
# 124     706  0.0924  
# 125     707  0.209   
# 126     708  0.0246  
# 127     709  0.0502  
# 128     710  0.0944  
# 129     711  0.103   
# 130     714  0.143   
# 131     715  0.126   
# 132     716  0.160   
# 133     717  0.187   
# 134     718 -0.0445  
# 135     720 -0.0462  
# 136     723 -0.0253  
# 137     724  0.138   
# 138     726  0.0300  
# 139     727  0.0398  
# 140     728  0.0812  
# 141     729 -0.0384  
# 142     730 -0.0210  
# 143     731  0.104   
# 144     732  0.138   
# 145     733 -0.0146  
# 146     734  0.00250 
# 147     735  0.132   
# 148     736  0.102   
# 149     737  0.121   
# 150     738  0.0615  
# 151     739 -0.0853  
# 152     740 -0.125   
# 153     741  0.221   
# 154     742 -0.0454  
# 155     743  0.0319  
# 156     744 -0.0174  
# 157     745  0.0152  
# 158     746  0.142   
# 159     747  0.145   
# 160     748  0.181   
# 161     749  0.184   
# 162     750  0.0195  
# 163     751  0.119   
# 164     752 -0.0847  
# 165     753  0.132   
# 166     754  0.0644  
# 167     755  0.0633  
# 168     756 -0.0180  
# 169     757  0.135   
# 170     758  0.169   
# 171     759  0.178   
# 172     760  0.0404  
# 173     761  0.0986  
# 174     762  0.123   
# 175     763 -0.0181  
# 176     764  0.278   
# 177     765  0.0235  
# 178     766 -0.0120  
# 179     767  0.00226 
# 180     768  0.455   
# 181     769  0.0746  
# 182     770 -0.0121  
# 183     771  0.140   
# 184     772  0.172   
# 185     773 -0.0435  
# 186     774  0.227   
# 187     775  0.0151  
# 188     776 -0.00627 
# 189     777  0.103   
# 190     778 -0.110   
# 191     779 -0.141   
# 192     780 -0.0448  
# 193     781 -0.137   
# 194     782 -0.0333  
# 195     783 -0.0828  
# 196     784 -0.0330  
# 197     785 -0.0962  
# 198     786  0.0582  
# 199     787  0.134   
# 200     788  0.362   
# 201     789 -0.0681  
# 202     790  0.216   
# 203     791  0.181   
# 204     792  0.0179  
# 205     793  0.0412  
# 206     794  0.148   
# 207     795  0.0407  
# 208     796  0.0997  
# 209     797  0.0659  
# 210     798  0.0127  
# 211     799  0.0811  
# 212     800  0.0662  
# 213     801  0.0459  
# 214     802  0.0616  
# 215     803  0.0434  
# 216     804  0.271   
# 217     805  0.282   
# 218     806  0.259   
# 219     807  0.0933  
# 220     808  0.135   
# 221     809  0.137   
# 222     810  0.0675  
# 223     811  0.0998  
# 224     812  0.220   
# 225     813  0.152   
# 226     814  0.142   
# 227     815  0.123   
# 228     816 -0.134   
# 229     817  0.146   
# 230     818  0.0616  
# 231     819  0.111   
# 232     820  0.116   
# 233     821  0.0749  
# 234     822 -0.101   
# 235     823  0.106   
# 236     824  0.00164 
# 237     825 -0.117   
# 238     826 -0.0805  
# 239     827  0.0178  
# 240     828  0.0955  
# 241     829  0.105   
# 242     830  0.0285  
# 243     831 -0.0799  
# 244     832  0.0766  
# 245     833  0.118   
# 246     834  0.00212 
# 247     835  0.158   
# 248     836  0.0413  
# 249     837  0.0232  
# 250     838  0.0399  
# 251     839  0.141   
# 252     840  0.231   
# 253     841  0.00239 
# 254     843  0.103   
# 255     844 -0.0221  
# 256     846  0.106   
# 257     847  0.0196  
# 258     848  0.0901  
# 259     849  0.110   
# 260     850  0.0354  
# 261     851  0.154   
# 262     852  0.112   
# 263     853  0       
# 264     854  0.266   
# 265     855  0.0484  
# 266     950  0       
# 267     951  0       
# 268     955  0       
# 269     956  0       
# 270     959  0       
# 271     960  0       
# 272     962 -0.136   
# 273     969 -0.136   
# 274     971 -0.0124  
# 275     972  0.0768  
# 276     976 -0.0124  
# 277     978  0.00271 
# 278     979  0.396   
# 279     982 -0.0916  
# 280     985  0       
# 281     988 -0.0692  
# 282     990 -0.0692  
# 283     991  0       
# 284     992 -0.104   
# 285     993  0       
# 286     995 -0.0638  
# 287     997 -0.0176  
# 288     999 -0.0638  
# 289    1001 -0.125   
# 290    1002  0       
# 291    1003  0       
# 292    1004 -0.0386  
# 293    1005 -0.168   
# 294    1006 -0.0651  
# 295    1007 -0.0763  
# 296    1008 -0.0979  
# 297    1009 -0.179   
# 298    1010 -0.127   
# 299    1011  0       
# 300    1012 -0.0176  
# 301    1013  0       
# 302    1015 -0.0236  
# 303    1016 -0.108   
# 304    1017 -0.0375  
# 305    1018 -0.0391  
# 306    1019 -0.0968 

# exact match (not lenient):
#   Ecds_ID mean_diff
# 1     375   0.0599 
# 2     377   0      
# 3     378   0      
# 4     379   0.201  
# 5     380   0.0670 
# 6     381   0      
# 7     382   0      
# 8     383   0      
# 9     391   0      
# 10     400   0      
# 11     401   0      
# 12     402   0      
# 13     408   0      
# 14     409   0      
# 15     410   0      
# 16     411   0      
# 17     412   0      
# 18     413   0      
# 19     541   0      
# 20     542   0      
# 21     543   0      
# 22     544   0      
# 23     545   0      
# 24     546   0      
# 25     547   0      
# 26     549   0      
# 27     550   0      
# 28     551   0      
# 29     552   0      
# 30     553   0      
# 31     554   0      
# 32     555   0      
# 33     556   0      
# 34     557   0      
# 35     558   0      
# 36     560   0      
# 37     561   0      
# 38     562   0      
# 39     564   0      
# 40     565   0      
# 41     566   0      
# 42     567   0      
# 43     568   0      
# 44     569   0      
# 45     570   0      
# 46     571   0      
# 47     572   0      
# 48     586   0.133  
# 49     587   0.361  
# 50     590   0      
# 51     591   0      
# 52     592   0      
# 53     593   0      
# 54     594   0.361  
# 55     595   0.361  
# 56     596   0.246  
# 57     597   0.361  
# 58     599   0.180  
# 59     600   0.180  
# 60     610   0      
# 61     611   0      
# 62     615   0      
# 63     616   0      
# 64     617   0      
# 65     618   0.183  
# 66     619  -0.395  
# 67     622  -0.198  
# 68     623   0      
# 69     628  -0.198  
# 70     629   0      
# 71     630   0      
# 72     631   0.0366 
# 73     650   0      
# 74     660   0      
# 75     669  -0.0392 
# 76     672   0      
# 77     677  -0.0392 
# 78     678   0      
# 79     679   0.0210 
# 80     680   0      
# 81     681  -0.0401 
# 82     682  -0.0601 
# 83     683  -0.0241 
# 84     684   0      
# 85     685  -0.0401 
# 86     686   0      
# 87     687   0.0101 
# 88     688  -0.0980 
# 89     689  -0.0401 
# 90     690   0      
# 91     692   0      
# 92     695   0      
# 93     696  -0.0401 
# 94     698  -0.0401 
# 95     701  -0.0241 
# 96     703   0      
# 97     706  -0.0401 
# 98     708   0.0157 
# 99     709  -0.0304 
# 100     710   0      
# 101     711   0.0550 
# 102     714   0.0611 
# 103     715   0.0925 
# 104     716   0.116  
# 105     717   0.210  
# 106     718   0.00717
# 107     720  -0.0377 
# 108     723   0.00343
# 109     724   0.0561 
# 110     726   0.00992
# 111     727  -0.0152 
# 112     728   0.0197 
# 113     729  -0.00167
# 114     730  -0.00785
# 115     731   0.133  
# 116     732   0.208  
# 117     733  -0.0294 
# 118     734  -0.0183 
# 119     735   0.145  
# 120     736   0.0479 
# 121     737   0.124  
# 122     738   0.0437 
# 123     739   0.0698 
# 124     740   0.0372 
# 125     741   0.179  
# 126     742  -0.119  
# 127     743  -0.0111 
# 128     744  -0.0298 
# 129     745   0.0221 
# 130     746   0.164  
# 131     747   0.175  
# 132     748   0.165  
# 133     749   0.170  
# 134     750   0.118  
# 135     751   0.155  
# 136     752  -0.0975 
# 137     753   0.0800 
# 138     754   0.0627 
# 139     755   0.0917 
# 140     756  -0.0503 
# 141     757   0.137  
# 142     758   0.114  
# 143     759   0.218  
# 144     760  -0.0143 
# 145     761   0.0818 
# 146     762   0.142  
# 147     763   0.00542
# 148     764   0.254  
# 149     765   0.0331 
# 150     766  -0.0324 
# 151     767  -0.0400 
# 152     768   0.496  
# 153     769   0.121  
# 154     770  -0.0125 
# 155     771   0.212  
# 156     772   0.125  
# 157     773  -0.0219 
# 158     774   0.287  
# 159     775   0.0534 
# 160     776   0.0709 
# 161     777   0.189  
# 162     778  -0.108  
# 163     779   0.0579 
# 164     780  -0.0425 
# 165     781  -0.0524 
# 166     782  -0.00543
# 167     783  -0.0335 
# 168     784  -0.0340 
# 169     785  -0.0845 
# 170     786  -0.00887
# 171     787   0.0850 
# 172     788   0.408  
# 173     789   0.0509 
# 174     790   0.427  
# 175     791   0.161  
# 176     792   0.0188 
# 177     793   0.123  
# 178     794   0.203  
# 179     795   0.0658 
# 180     796   0.141  
# 181     797  -0.0102 
# 182     798   0.0749 
# 183     799   0.0799 
# 184     800   0.00246
# 185     801   0.0125 
# 186     802  -0.0268 
# 187     803  -0.00618
# 188     804   0.371  
# 189     805   0.388  
# 190     806   0.315  
# 191     807  -0.00608
# 192     808   0.0198 
# 193     809   0.116  
# 194     810  -0.0153 
# 195     811   0.0194 
# 196     812   0.319  
# 197     813   0.0417 
# 198     814   0.291  
# 199     815   0.0988 
# 200     816  -0.0296 
# 201     817   0.0551 
# 202     818  -0.0294 
# 203     819   0.0417 
# 204     820   0.0367 
# 205     821   0.0375 
# 206     822  -0.00456
# 207     823   0.00105
# 208     824   0.0314 
# 209     825  -0.0217 
# 210     826  -0.00132
# 211     827   0.0338 
# 212     828   0.119  
# 213     829   0.0802 
# 214     830   0.0254 
# 215     831   0.00650
# 216     832   0.0574 
# 217     833   0.131  
# 218     834   0.0367 
# 219     835   0.0554 
# 220     836   0.0674 
# 221     837   0.0639 
# 222     838   0.0845 
# 223     839  -0.0183 
# 224     840   0.161  
# 225     841   0.00442
# 226     843   0.00330
# 227     844  -0.0621 
# 228     846   0.0697 
# 229     847  -0.0249 
# 230     848   0.108  
# 231     849  -0.00557
# 232     850   0.00397
# 233     851   0.172  
# 234     852   0.0517 
# 235     853   0      
# 236     854   0.309  
# 237     855   0.0780 
# 238     955   0      
# 239     956   0      
# 240     960   0      
# 241     972   0.245  
# 242     978   0.0815 
# 243     979   0.115  
# 244     982   0      
# 245     985   0      
# 246     988   0      
# 247     990   0      
# 248     991   0      
# 249     993   0      
# 250     995   0      
# 251     997   0      
# 252     999   0      
# 253    1001  -0.0523 
# 254    1002   0      
# 255    1003   0      
# 256    1004   0      
# 257    1005  -0.158  
# 258    1006  -0.00955
# 259    1007  -0.0535 
# 260    1008  -0.0595 
# 261    1009  -0.0816 
# 262    1010   0.0219 
# 263    1011   0      
# 264    1012   0      
# 265    1013   0      
# 266    1015   0      
# 267    1016   0      
# 268    1017   0      
# 269    1018   0      
# 270    1019   0 

# Relationship to sample size by ecodistrict 
accuracy_with_fallback_grouped_cons2 %>%
  group_by(Ecds_ID) %>%
  summarise(mean_acc = mean(native_metric, na.rm = TRUE),
              mean_n = mean(N_evaluated)) %>%
  print(n = Inf)

ecodistrict_summary <- accuracy_with_fallback_grouped_cons2 %>%
  group_by(Ecds_ID) %>%
  summarise(
    mean_acc = mean(native_metric, na.rm = TRUE),
    mean_n = mean(N_evaluated),
    .groups = "drop"
  )

ggplot(
  ecodistrict_summary,
  aes(mean_n, mean_acc)
) +
  geom_point() +
  scale_x_log10() +
  theme_bw() +
  labs(
    x = "Mean N evaluated (log10 scale)",
    y = "Mean native metric",
    title = "Native metric versus sample size by ecodistrict"
  )

# Which ecodistricts drive native omission?
accuracy_with_fallback_grouped_cons2 %>%
  group_by(Ecds_ID) %>%
  summarise(native_precision = mean(Precision_2, na.rm=TRUE),
               native_recall = mean(Recall_2, na.rm=TRUE)) %>%
  mutate(omission_gap = native_precision - native_recall) %>%
  arrange(omission_gap) %>%
  print(n = Inf)

# with ACI2025:
#   Ecds_ID native_precision native_recall omission_gap
# 1     972           0.328        0.921      -0.593   
# 2     971           0.401        0.988      -0.587   
# 3     976           0.401        0.988      -0.587   
# 4     978           0.516        0.956      -0.440   
# 5    1010           0.492        0.753      -0.261   
# 6     795           0.641        0.863      -0.222   
# 7     834           0.721        0.922      -0.201   
# 8     669           0.200        0.400      -0.200   
# 9     805           0.598        0.786      -0.189   
# 10     979           0.802        0.963      -0.161   
# 11     835           0.749        0.902      -0.153   
# 12     826           0.684        0.808      -0.125   
# 13     997           0.135        0.254      -0.119   
# 14     833           0.728        0.846      -0.118   
# 15     803           0.594        0.711      -0.117   
# 16     814           0.724        0.832      -0.109   
# 17     720           0.237        0.340      -0.103   
# 18     677           0.237        0.340      -0.103   
# 19     626           0.0405       0.143      -0.102   
# 20     821           0.692        0.782      -0.0901  
# 21     809           0.731        0.810      -0.0789  
# 22     797           0.678        0.753      -0.0750  
# 23     832           0.790        0.864      -0.0742  
# 24     811           0.731        0.802      -0.0715  
# 25     806           0.722        0.789      -0.0670  
# 26     823           0.660        0.727      -0.0668  
# 27     848           0.330        0.390      -0.0603  
# 28     824           0.686        0.743      -0.0572  
# 29     807           0.692        0.744      -0.0526  
# 30     827           0.752        0.799      -0.0468  
# 31     831           0.679        0.723      -0.0441  
# 32     992           0.0833       0.125      -0.0417  
# 33     836           0.767        0.804      -0.0369  
# 34     794           0.601        0.633      -0.0327  
# 35     804           0.713        0.732      -0.0191  
# 36     847           0.308        0.326      -0.0179  
# 37     829           0.694        0.712      -0.0174  
# 38     837           0.760        0.777      -0.0167  
# 39    1011           0.333        0.348      -0.0145  
# 40     830           0.788        0.800      -0.0126  
# 41     672           0.274        0.280      -0.00526 
# 42     787           0.723        0.723      -0.000150
# 43     586           0.25         0.25        0       
# 44     590           0            0           0       
# 45     617           0            0           0       
# 46     622           0            0           0       
# 47     623           0            0           0       
# 48     629           0            0           0       
# 49     630           0            0           0       
# 50     647           0            0           0       
# 51     650           0            0           0       
# 52     652           0            0           0       
# 53     657           0            0           0       
# 54     678           0            0           0       
# 55     710           0            0           0       
# 56     962           0            0           0       
# 57     969           0            0           0       
# 58     982           0            0           0       
# 59     817           0.743        0.742       0.000348
# 60     812           0.707        0.704       0.00316 
# 61    1001           0.481        0.468       0.0126  
# 62     825           0.652        0.640       0.0127  
# 63     838           0.746        0.727       0.0195  
# 64     796           0.681        0.661       0.0197  
# 65     818           0.691        0.671       0.0198  
# 66     810           0.704        0.681       0.0229  
# 67     802           0.696        0.671       0.0257  
# 68     820           0.726        0.699       0.0268  
# 69     822           0.742        0.713       0.0298  
# 70     813           0.731        0.696       0.0353  
# 71     815           0.740        0.696       0.0444  
# 72     724           0.285        0.240       0.0456  
# 73     816           0.713        0.655       0.0579  
# 74     595           0.314        0.254       0.0598  
# 75     597           0.314        0.254       0.0598  
# 76     790           0.558        0.492       0.0662  
# 77     379           0.425        0.358       0.0668  
# 78     828           0.721        0.653       0.0677  
# 79     593           0.333        0.263       0.0702  
# 80    1015           0.344        0.273       0.0708  
# 81     819           0.767        0.696       0.0714  
# 82     762           0.647        0.566       0.0810  
# 83     618           0.5          0.417       0.0833  
# 84     844           0.344        0.256       0.0883  
# 85     792           0.642        0.553       0.0890  
# 86     793           0.628        0.530       0.0984  
# 87     723           0.307        0.207       0.0998  
# 88    1005           0.515        0.406       0.108   
# 89    1012           0.414        0.305       0.109   
# 90     791           0.696        0.585       0.111   
# 91     841           0.371        0.258       0.112   
# 92     999           0.354        0.236       0.118   
# 93    1006           0.392        0.256       0.136   
# 94     855           0.409        0.267       0.142   
# 95     799           0.663        0.519       0.143   
# 96     789           0.721        0.571       0.150   
# 97     777           0.723        0.572       0.151   
# 98     709           0.517        0.366       0.152   
# 99    1009           0.469        0.308       0.161   
# 100     628           0.167        0.00490     0.162   
# 101     800           0.666        0.496       0.170   
# 102     784           0.712        0.537       0.174   
# 103     780           0.683        0.501       0.182   
# 104     808           0.702        0.519       0.183   
# 105     783           0.651        0.467       0.184   
# 106     779           0.730        0.545       0.184   
# 107     380           0.425        0.239       0.186   
# 108    1007           0.25         0.0635      0.187   
# 109     771           0.748        0.561       0.187   
# 110     764           0.541        0.353       0.188   
# 111     591           0.333        0.143       0.190   
# 112     754           0.490        0.299       0.191   
# 113     801           0.687        0.496       0.191   
# 114    1017           0.410        0.209       0.201   
# 115     786           0.740        0.535       0.205   
# 116     768           0.766        0.557       0.209   
# 117     587           0.471        0.254       0.217   
# 118     594           0.471        0.254       0.217   
# 119    1019           0.505        0.286       0.218   
# 120     995           0.509        0.275       0.234   
# 121     798           0.659        0.424       0.236   
# 122     699           0.324        0.0850      0.239   
# 123     775           0.734        0.492       0.242   
# 124     843           0.282        0.0390      0.242   
# 125     714           0.25         0.00575     0.244   
# 126    1018           0.571        0.325       0.246   
# 127     761           0.584        0.330       0.255   
# 128     782           0.674        0.417       0.257   
# 129     846           0.532        0.276       0.257   
# 130     375           0.451        0.190       0.261   
# 131     748           0.475        0.213       0.262   
# 132     734           0.671        0.407       0.264   
# 133     631           0.411        0.146       0.265   
# 134     770           0.719        0.454       0.265   
# 135     776           0.769        0.495       0.274   
# 136     736           0.643        0.361       0.281   
# 137     839           0.333        0.05        0.283   
# 138     760           0.522        0.238       0.284   
# 139     988           0.389        0.104       0.285   
# 140     990           0.389        0.104       0.285   
# 141     767           0.688        0.395       0.292   
# 142     596           0.5          0.204       0.296   
# 143     763           0.578        0.276       0.302   
# 144     718           0.431        0.122       0.309   
# 145     739           0.785        0.473       0.312   
# 146     726           0.493        0.179       0.314   
# 147     851           0.758        0.444       0.315   
# 148     769           0.756        0.436       0.320   
# 149     592           0.333        0.00862     0.325   
# 150     660           0.333        0.00862     0.325   
# 151     680           0.333        0.00687     0.326   
# 152     690           0.333        0.00687     0.326   
# 153     703           0.333        0.00687     0.326   
# 154     850           0.451        0.124       0.327   
# 155     749           0.579        0.247       0.331   
# 156     757           0.760        0.411       0.349   
# 157     600           0.485        0.136       0.350   
# 158     772           0.714        0.361       0.353   
# 159     733           0.580        0.220       0.360   
# 160     716           0.395        0.0330      0.362   
# 161     759           0.572        0.209       0.363   
# 162     774           0.811        0.442       0.369   
# 163     849           0.499        0.128       0.371   
# 164     753           0.6          0.228       0.372   
# 165    1008           0.702        0.329       0.373   
# 166    1016           0.493        0.120       0.373   
# 167     740           0.739        0.353       0.386   
# 168     743           0.718        0.329       0.389   
# 169     788           0.712        0.319       0.394   
# 170     781           0.669        0.275       0.394   
# 171     785           0.702        0.308       0.394   
# 172     735           0.624        0.228       0.395   
# 173     778           0.704        0.305       0.399   
# 174     747           0.564        0.164       0.399   
# 175     742           0.635        0.235       0.400   
# 176     691           0.486        0.0850      0.401   
# 177     741           0.657        0.251       0.406   
# 178     728           0.488        0.0759      0.412   
# 179     758           0.663        0.241       0.422   
# 180     765           0.616        0.194       0.422   
# 181     766           0.5          0.0767      0.423   
# 182     755           0.701        0.272       0.429   
# 183     750           0.694        0.260       0.434   
# 184     711           0.488        0.0409      0.447   
# 185     751           0.698        0.251       0.447   
# 186     773           0.729        0.276       0.454   
# 187     752           0.595        0.130       0.465   
# 188     686           0.5          0.0315      0.469   
# 189     737           0.570        0.0997      0.470   
# 190     682           0.5          0.0279      0.472   
# 191    1004           0.5          0.0256      0.474   
# 192     679           0.5          0.0255      0.475   
# 193     696           0.5          0.0255      0.475   
# 194     706           0.5          0.0255      0.475   
# 195     756           0.645        0.169       0.476   
# 196     684           0.5          0.0211      0.479   
# 197     991           0.5          0.0208      0.479   
# 198     708           0.5          0.0191      0.481   
# 199     610           0.5          0.0172      0.483   
# 200     611           0.5          0.0172      0.483   
# 201     745           0.675        0.190       0.485   
# 202     619           0.5          0.00980     0.490   
# 203     615           0.5          0.00862     0.491   
# 204     616           0.5          0.00862     0.491   
# 205     744           0.679        0.188       0.492   
# 206     731           0.679        0.187       0.492   
# 207     692           0.5          0.00515     0.495   
# 208     852           0.661        0.159       0.502   
# 209     729           0.665        0.154       0.511   
# 210     599           0.647        0.136       0.511   
# 211     695           0.6          0.0869      0.513   
# 212     746           0.788        0.243       0.545   
# 213     688           0.6          0.0330      0.567   
# 214     701           0.6          0.0271      0.573   
# 215     738           0.714        0.133       0.581   
# 216     687           0.667        0.0470      0.620   
# 217     730           0.760        0.137       0.623   
# 218    1002           1            0.366       0.634   
# 219    1003           1            0.366       0.634   
# 220     732           0.834        0.195       0.639   
# 221     681           0.667        0.0255      0.641   
# 222     685           0.667        0.0255      0.641   
# 223     689           0.667        0.0255      0.641   
# 224     698           0.667        0.0255      0.641   
# 225     683           0.75         0.0363      0.714   
# 226     727           0.916        0.111       0.805   
# 227     993           0.974        0.154       0.820   
# 228    1013           0.974        0.154       0.820   
# 229     985           1            0.0625      0.938   
# 230     374           0          NaN         NaN       
# 231     377           0          NaN         NaN       
# 232     378           0          NaN         NaN       
# 233     381           0          NaN         NaN       
# 234     382           0          NaN         NaN       
# 235     383           0          NaN         NaN       
# 236     391           0          NaN         NaN       
# 237     400           0          NaN         NaN       
# 238     401           0          NaN         NaN       
# 239     402           0          NaN         NaN       
# 240     408           0          NaN         NaN       
# 241     409         NaN          NaN         NaN       
# 242     410           0          NaN         NaN       
# 243     411           0          NaN         NaN       
# 244     412         NaN          NaN         NaN       
# 245     413         NaN          NaN         NaN       
# 246     541         NaN          NaN         NaN       
# 247     542         NaN          NaN         NaN       
# 248     543         NaN          NaN         NaN       
# 249     544         NaN          NaN         NaN       
# 250     545         NaN          NaN         NaN       
# 251     546         NaN          NaN         NaN       
# 252     547         NaN          NaN         NaN       
# 253     549           0          NaN         NaN       
# 254     550           0          NaN         NaN       
# 255     551         NaN          NaN         NaN       
# 256     552         NaN          NaN         NaN       
# 257     553         NaN          NaN         NaN       
# 258     554         NaN          NaN         NaN       
# 259     555         NaN          NaN         NaN       
# 260     556         NaN          NaN         NaN       
# 261     557         NaN          NaN         NaN       
# 262     558         NaN          NaN         NaN       
# 263     559           0          NaN         NaN       
# 264     560         NaN          NaN         NaN       
# 265     561         NaN          NaN         NaN       
# 266     562         NaN          NaN         NaN       
# 267     564         NaN          NaN         NaN       
# 268     565         NaN          NaN         NaN       
# 269     566         NaN          NaN         NaN       
# 270     567         NaN          NaN         NaN       
# 271     568         NaN          NaN         NaN       
# 272     569         NaN          NaN         NaN       
# 273     570         NaN          NaN         NaN       
# 274     571         NaN          NaN         NaN       
# 275     572         NaN          NaN         NaN       
# 276     583           0          NaN         NaN       
# 277     585           0          NaN         NaN       
# 278     588           0          NaN         NaN       
# 279     589           0          NaN         NaN       
# 280     598           0          NaN         NaN       
# 281     601           0          NaN         NaN       
# 282     607           0          NaN         NaN       
# 283     612           0          NaN         NaN       
# 284     658           0          NaN         NaN       
# 285     659           0          NaN         NaN       
# 286     661           0          NaN         NaN       
# 287     662           0          NaN         NaN       
# 288     693           0          NaN         NaN       
# 289     694           0          NaN         NaN       
# 290     697           0          NaN         NaN       
# 291     700           0          NaN         NaN       
# 292     702         NaN          NaN         NaN       
# 293     704           0          NaN         NaN       
# 294     705           0          NaN         NaN       
# 295     707           0          NaN         NaN       
# 296     715           0          NaN         NaN       
# 297     717           0          NaN         NaN       
# 298     840           0          NaN         NaN       
# 299     853           0          NaN         NaN       
# 300     854           0          NaN         NaN       
# 301     950         NaN          NaN         NaN       
# 302     951         NaN          NaN         NaN       
# 303     955         NaN          NaN         NaN       
# 304     956         NaN          NaN         NaN       
# 305     959         NaN          NaN         NaN       
# 306     960         NaN            0         NaN

# exact match (not lenient):
#   Ecds_ID native_precision native_recall omission_gap
# 1    1010            0.647       0.700      -0.0525  
# 2     669            0.400       0.400       0       
# 3     677            0.400       0.400       0       
# 4     720            0.400       0.400       0.000904
# 5     978            0.994       0.966       0.0287  
# 6     979            0.993       0.963       0.0304  
# 7     795            0.912       0.863       0.0496  
# 8     972            0.983       0.921       0.0624  
# 9     834            0.956       0.886       0.0698  
# 10     835            0.981       0.902       0.0790  
# 11     833            0.954       0.860       0.0935  
# 12     832            0.985       0.864       0.121   
# 13     848            0.499       0.377       0.122   
# 14     821            0.915       0.780       0.135   
# 15     809            0.948       0.810       0.138   
# 16     805            0.929       0.781       0.148   
# 17     806            0.953       0.805       0.148   
# 18     794            0.784       0.633       0.151   
# 19     829            0.941       0.783       0.158   
# 20     827            0.961       0.778       0.183   
# 21     830            0.984       0.800       0.183   
# 22     811            0.987       0.802       0.185   
# 23     824            0.934       0.745       0.190   
# 24     797            0.944       0.753       0.191   
# 25    1009            0.619       0.426       0.193   
# 26     815            0.957       0.760       0.197   
# 27     814            0.988       0.788       0.200   
# 28     803            0.918       0.711       0.207   
# 29     826            0.975       0.769       0.207   
# 30     817            0.956       0.742       0.213   
# 31     823            0.95        0.727       0.223   
# 32     828            0.985       0.760       0.225   
# 33     807            0.972       0.744       0.228   
# 34     796            0.891       0.661       0.230   
# 35     804            0.970       0.732       0.238   
# 36     836            1           0.760       0.240   
# 37     831            0.929       0.687       0.242   
# 38     820            0.948       0.699       0.248   
# 39    1006            0.622       0.373       0.250   
# 40     586            0.5         0.25        0.25    
# 41     800            0.908       0.654       0.254   
# 42     838            0.943       0.689       0.255   
# 43     787            0.979       0.723       0.256   
# 44     825            0.914       0.654       0.260   
# 45    1007            0.534       0.271       0.263   
# 46     762            0.834       0.566       0.267   
# 47     812            0.973       0.704       0.268   
# 48     808            0.954       0.683       0.271   
# 49     837            0.985       0.713       0.272   
# 50     802            0.947       0.671       0.276   
# 51     813            0.981       0.696       0.285   
# 52     810            0.971       0.681       0.291   
# 53     818            0.965       0.671       0.293   
# 54    1005            0.707       0.406       0.301   
# 55     819            0.998       0.696       0.302   
# 56     822            0.990       0.673       0.317   
# 57     792            0.873       0.553       0.320   
# 58     816            0.990       0.660       0.331   
# 59     724            0.570       0.240       0.331   
# 60    1001            0.848       0.513       0.335   
# 61     855            0.672       0.333       0.339   
# 62     791            0.934       0.585       0.349   
# 63     379            0.7         0.341       0.359   
# 64     844            0.632       0.272       0.359   
# 65     789            0.893       0.532       0.360   
# 66     790            0.875       0.514       0.361   
# 67     793            0.903       0.530       0.373   
# 68     847            0.683       0.310       0.374   
# 69     776            0.879       0.495       0.384   
# 70     782            0.910       0.524       0.386   
# 71     798            0.898       0.508       0.390   
# 72     841            0.707       0.317       0.390   
# 73     777            0.977       0.572       0.405   
# 74     784            0.944       0.537       0.407   
# 75     771            0.969       0.561       0.408   
# 76     764            0.739       0.319       0.420   
# 77     799            0.942       0.519       0.423   
# 78     780            0.922       0.497       0.425   
# 79     709            0.792       0.366       0.426   
# 80     754            0.734       0.299       0.435   
# 81     723            0.648       0.207       0.441   
# 82     783            0.909       0.467       0.441   
# 83     761            0.779       0.330       0.449   
# 84     757            0.896       0.438       0.457   
# 85     760            0.696       0.238       0.458   
# 86     786            0.994       0.535       0.459   
# 87     775            0.951       0.492       0.460   
# 88     768            0.988       0.523       0.465   
# 89     843            0.594       0.122       0.472   
# 90    1008            0.784       0.311       0.474   
# 91     770            0.931       0.454       0.476   
# 92     774            0.927       0.442       0.485   
# 93     375            0.676       0.190       0.486   
# 94     849            0.699       0.211       0.488   
# 95     619            0.5         0.00980     0.490   
# 96     714            0.5         0.00575     0.494   
# 97     734            0.901       0.407       0.494   
# 98     622            0.5         0.00490     0.495   
# 99     628            0.5         0.00490     0.495   
# 100     736            0.860       0.361       0.499   
# 101     801            0.995       0.496       0.500   
# 102     767            0.901       0.395       0.506   
# 103     763            0.783       0.276       0.506   
# 104     788            0.988       0.474       0.514   
# 105     749            0.826       0.303       0.524   
# 106     781            0.900       0.367       0.532   
# 107     769            0.975       0.436       0.539   
# 108     726            0.718       0.179       0.539   
# 109     718            0.718       0.179       0.539   
# 110     748            0.760       0.213       0.547   
# 111     733            0.773       0.220       0.553   
# 112     759            0.763       0.209       0.554   
# 113     741            0.907       0.352       0.555   
# 114     846            0.836       0.276       0.560   
# 115     779            0.981       0.416       0.565   
# 116     773            0.953       0.382       0.571   
# 117     735            0.816       0.245       0.571   
# 118     851            0.948       0.369       0.579   
# 119     766            0.682       0.103       0.579   
# 120     618            1           0.417       0.583   
# 121     850            0.712       0.126       0.586   
# 122     772            0.943       0.356       0.587   
# 123     785            0.903       0.308       0.595   
# 124     732            0.846       0.250       0.596   
# 125     765            0.844       0.243       0.601   
# 126     743            0.934       0.329       0.605   
# 127     753            0.902       0.294       0.608   
# 128     380            0.85        0.239       0.611   
# 129     742            0.846       0.235       0.612   
# 130     758            0.929       0.313       0.615   
# 131     778            0.929       0.305       0.623   
# 132     747            0.789       0.164       0.625   
# 133     716            0.658       0.0330      0.625   
# 134     737            0.797       0.166       0.631   
# 135    1002            1           0.366       0.634   
# 136    1003            1           0.366       0.634   
# 137     752            0.765       0.130       0.635   
# 138     739            0.982       0.331       0.651   
# 139    1011            1           0.348       0.652   
# 140     756            0.829       0.169       0.660   
# 141     745            0.917       0.255       0.662   
# 142     755            0.935       0.272       0.663   
# 143     751            0.917       0.251       0.666   
# 144     750            0.988       0.301       0.686   
# 145    1015            0.991       0.305       0.686   
# 146     587            0.941       0.254       0.687   
# 147     594            0.941       0.254       0.687   
# 148     595            0.941       0.254       0.687   
# 149     597            0.941       0.254       0.687   
# 150     631            0.833       0.146       0.687   
# 151    1018            0.994       0.302       0.692   
# 152    1019            1           0.286       0.714   
# 153     672            1           0.280       0.720   
# 154     990            1           0.275       0.725   
# 155     991            1           0.275       0.725   
# 156     729            0.881       0.154       0.727   
# 157     852            0.909       0.176       0.734   
# 158     593            1           0.263       0.737   
# 159     746            0.985       0.243       0.742   
# 160     997            1           0.254       0.746   
# 161     744            0.951       0.188       0.763   
# 162     999            1           0.236       0.764   
# 163     727            0.944       0.175       0.768   
# 164     731            0.958       0.187       0.771   
# 165     728            0.853       0.0759      0.777   
# 166    1004            1           0.218       0.782   
# 167    1017            1           0.209       0.791   
# 168     730            0.994       0.201       0.792   
# 169     596            1           0.204       0.796   
# 170     740            0.994       0.195       0.798   
# 171     995            0.987       0.168       0.818   
# 172     993            0.974       0.154       0.820   
# 173    1012            0.974       0.154       0.820   
# 174    1013            0.974       0.154       0.820   
# 175     599            0.971       0.136       0.835   
# 176     600            0.971       0.136       0.835   
# 177     738            0.986       0.133       0.853   
# 178     591            1           0.143       0.857   
# 179    1016            0.986       0.120       0.867   
# 180     988            1           0.119       0.881   
# 181     695            1           0.0869      0.913   
# 182     682            1           0.0704      0.930   
# 183     711            0.975       0.0409      0.934   
# 184     985            1           0.0625      0.938   
# 185     681            1           0.0538      0.946   
# 186     685            1           0.0538      0.946   
# 187     689            1           0.0538      0.946   
# 188     696            1           0.0538      0.946   
# 189     698            1           0.0538      0.946   
# 190     706            1           0.0538      0.946   
# 191     683            1           0.0533      0.947   
# 192     839            1           0.05        0.95    
# 193     687            1           0.0470      0.953   
# 194     688            1           0.0330      0.967   
# 195     701            1           0.0323      0.968   
# 196     686            1           0.0315      0.969   
# 197     679            1           0.0255      0.975   
# 198     684            1           0.0211      0.979   
# 199     708            1           0.0191      0.981   
# 200     610            1           0.0172      0.983   
# 201     611            1           0.0172      0.983   
# 202     592            1           0.00862     0.991   
# 203     615            1           0.00862     0.991   
# 204     616            1           0.00862     0.991   
# 205     660            1           0.00862     0.991   
# 206     680            1           0.00687     0.993   
# 207     690            1           0.00687     0.993   
# 208     703            1           0.00687     0.993   
# 209     692            1           0.00515     0.995   
# 210     377            0         NaN         NaN       
# 211     378            0         NaN         NaN       
# 212     381            0         NaN         NaN       
# 213     382            0         NaN         NaN       
# 214     383            0         NaN         NaN       
# 215     391            0         NaN         NaN       
# 216     400            0         NaN         NaN       
# 217     401            0         NaN         NaN       
# 218     402            0         NaN         NaN       
# 219     408            0         NaN         NaN       
# 220     409            0         NaN         NaN       
# 221     410            0         NaN         NaN       
# 222     411          NaN         NaN         NaN       
# 223     412          NaN         NaN         NaN       
# 224     413          NaN         NaN         NaN       
# 225     541          NaN         NaN         NaN       
# 226     542          NaN         NaN         NaN       
# 227     543          NaN         NaN         NaN       
# 228     544          NaN         NaN         NaN       
# 229     545          NaN         NaN         NaN       
# 230     546          NaN         NaN         NaN       
# 231     547          NaN         NaN         NaN       
# 232     549          NaN         NaN         NaN       
# 233     550            0         NaN         NaN       
# 234     551          NaN         NaN         NaN       
# 235     552          NaN         NaN         NaN       
# 236     553          NaN         NaN         NaN       
# 237     554          NaN         NaN         NaN       
# 238     555          NaN         NaN         NaN       
# 239     556          NaN         NaN         NaN       
# 240     557          NaN         NaN         NaN       
# 241     558          NaN         NaN         NaN       
# 242     560          NaN         NaN         NaN       
# 243     561          NaN         NaN         NaN       
# 244     562          NaN         NaN         NaN       
# 245     564          NaN         NaN         NaN       
# 246     565          NaN         NaN         NaN       
# 247     566          NaN         NaN         NaN       
# 248     567          NaN         NaN         NaN       
# 249     568          NaN         NaN         NaN       
# 250     569          NaN         NaN         NaN       
# 251     570          NaN         NaN         NaN       
# 252     571          NaN         NaN         NaN       
# 253     572          NaN         NaN         NaN       
# 254     590          NaN           0         NaN       
# 255     617          NaN           0         NaN       
# 256     623          NaN           0         NaN       
# 257     629          NaN           0         NaN       
# 258     630          NaN           0         NaN       
# 259     650          NaN           0         NaN       
# 260     678          NaN           0         NaN       
# 261     710          NaN           0         NaN       
# 262     715            0         NaN         NaN       
# 263     717            0         NaN         NaN       
# 264     840            0         NaN         NaN       
# 265     853            0         NaN         NaN       
# 266     854            0         NaN         NaN       
# 267     955          NaN         NaN         NaN       
# 268     956          NaN         NaN         NaN       
# 269     960          NaN           0         NaN       
# 270     982          NaN           0         NaN 

# Spatial autocorrelation
# Join accuracy table to ecodistrict boundaries
# Map:
# native_metric
# tame_metric
# not_metric
# native–tame difference
# fallback frequency

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