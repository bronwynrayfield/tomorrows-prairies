
#######################################################################################
##              PRELIMINARY ACCURACY ASSESSMENT: NATIVE AND TAME CLASSES             ##
#######################################################################################

##----------------------- Alison Long - June 2026 -------------------------

# step 8 (multiclass): assess accuracy of datasets mapping grasslands
# against ground-truthed point data (2015 to present), with 3 classes:
#   0 = not grass, 1 = tame grass, 2 = native grass
#
# THIS VERSION COMPARES: native vs non-native (tame, crop, shrub, wetland) points to 
#      native vs non-native (tame, not grass) raster values and tame vs non-tame (native, 
#      crop, shrub, wetland) points to tame vs non-tame (native, not grass) raster values
#
# Key changes from binary version:
#   - Uses LC_class3 (0/1/2) instead of LC_class3 (0/1)
#   - Rasters expected as 0/1/2
#   - Time-aligned rules:
#         native (2): latest allowed_year <= Yer_Cll
#         tame   (1): allowed_year == Yer_Cll
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
library(mltools)

# Small utility: "x %||% y" returns y if x is NULL
`%||%` <- function(x, y) if (is.null(x)) y else x

# ----------------------------
# 2) File paths 
# ----------------------------
points_excel_path <- "R:/CanadianGrasslands/WG_Classification/processed/groundtruthed_data_points/groundtruth_points_with_grass_raster_values_natvtam_pli_sep_ge2015_23Jun2026_clean_withnongrass.xlsx"
metadata_csv_path <- "R:/CanadianGrasslands/WG_Classification/github/tomorrows-prairies/processeddatasets/Grass/metadata_full_extent_nattam_pli_sep.csv"

out_dir <- "R:/CanadianGrasslands/WG_Classification/results/accuracy_assessment/native_v_tame/compiled_derived"
#dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ----------------------------
# 3) Configuration
# ----------------------------
# IMPORTANT: Handle PLI as four separate rasters in the points file and metadata.
pli_variants <- c("pli_ap", "pli_cu", "pli_mg", "pli_mm")

# Expected raster columns (3-class values 0/1/2)
raster_cols <- c("aci","ngi","cgr","eccc",
                 "cnvc","pudl",
                 pli_variants,            # four atomic PLI rasters
                 "gpw","plvi","mgi","pglc")

# Grouping helper: collapse pli_* -> "pli"; all others remain themselves
raster_to_group <- function(x) ifelse(x %in% pli_variants, "pli", x)

# ----------------------------
# 4) Load points check for required columns
# ----------------------------
#points <- readxl::read_excel(points_excel_path)
points <- readxl::read_excel(
  points_excel_path,
  col_types = c("numeric", "text", "text", "numeric", rep("guess", 10), rep("numeric", 19))
)

names(points)
str(points)
unique(points$Training)
table(points$Training, useNA = "ifany")
unique(points$LC_class3)
table(points$LC_class3, useNA = "ifany")

# Remove HFI
points <- points %>%
  select(-hfi)

# Add point ID column
points <- points %>%
  mutate(Pt_ID = row_number()) %>%
  relocate(Pt_ID)

# Drop rows where ecodistrict ID is 0
points <- points %>%
  filter(Ecds_ID != 0)

# Define required columns
required_cols <- c("Pt_ID", "LC_class3", "Training","Sorc_ID", "Ecds_ID","Yer_Cll", "LC_Cl_N")
missing_req <- setdiff(required_cols, names(points))
if (length(missing_req) > 0) {
  stop("Missing required columns in points file: ", paste(missing_req, collapse = ", "))
}

#Check LC types present
table(points$LC_Cl_N, useNA = "ifany")

# Rename all non-grass points as OTHER
points <- points %>%
  mutate(
    LC_Cl_N = if_else(
      LC_Cl_N %in% c("ANNUAL CROPLAND", "WOODY", "WETLAND", "OPEN WATER"),
      "OTHER",
      LC_Cl_N))

# Rename all MIXED points as TAME
points <- points %>%
  mutate(
    LC_Cl_N = if_else(
      LC_Cl_N == "MIXED",
      "TAME",
      LC_Cl_N))

#Check LC types present
table(points$LC_Cl_N, useNA = "ifany")

# How many OTHER points for each year?
pts.oth <- points %>%
  filter(LC_Cl_N == "OTHER")

pts.oth %>%
  group_by(Yer_Cll) %>%
  summarise(total_count = n(), .groups = "drop")

# How many TAME points for each year?
pts.tam <- points %>%
  filter(LC_Cl_N == "TAME")

pts.tam %>%
  group_by(Yer_Cll) %>%
  summarise(total_count = n(), .groups = "drop")

# How many NATIVE points for each year?
pts.nat <- points %>%
  filter(LC_Cl_N == "NATIVE")

pts.nat %>%
  group_by(Yer_Cll) %>%
  summarise(total_count = n(), .groups = "drop")

#Check LC types present
table(points$LC_class3, useNA = "ifany")

# Identify which raster columns actually exist in the points file
present_rasters <- intersect(raster_cols, names(points))
if (length(present_rasters) == 0) {
  stop("None of the expected raster columns were found in the points file.")
}

# ----------------------------
# 5) Load metadata and allow MULTIPLE years per raster
# ----------------------------
meta_raw <- readr::read_csv(metadata_csv_path, show_col_types = FALSE)

parse_year_vector <- function(y) {
  # Accept numbers or delimited strings "2020|2021" or "2020, 2021"
  if (is.null(y) || (length(y) == 1 && is.na(y))) return(integer(0))
  if (is.numeric(y)) return(as.integer(stats::na.omit(y)))
  y_chr <- as.character(y)
  y_chr <- stringr::str_trim(y_chr)
  y_chr[y_chr %in% c("", "NA", "NaN", "null", "NULL", "missing")] <- NA_character_
  y_chr <- y_chr[!is.na(y_chr)]
  if (length(y_chr) == 0) return(integer(0))
  parts <- unlist(strsplit(y_chr, "[|,]"))
  parts <- stringr::str_trim(parts)
  suppressWarnings(as.integer(stats::na.omit(as.numeric(parts))))
}

if (all(c("raster", "year") %in% names(meta_raw))) {
  # Long format
  meta_norm <- meta_raw %>%
    transmute(
      raster = str_trim(tolower(as.character(raster))),
      year_raw = year
    ) %>%
    group_by(raster) %>%
    summarise(years = list(parse_year_vector(year_raw)), .groups = "drop")
} else if (nrow(meta_raw) == 1 && any(intersect(raster_cols, names(meta_raw)) %in% names(meta_raw))) {
  # Wide single-row
  meta_norm <- tibble(
    raster = raster_cols,
    years  = map(raster_cols, ~ parse_year_vector(meta_raw[[.x]]))
  ) %>%
    mutate(raster = tolower(raster))
} else {
  stop(
    "Metadata format not recognized. Provide either:\n",
    "  (A) Long CSV with columns: raster, year (multiple rows allowed per raster), OR\n",
    "  (B) Wide 1-row CSV with raster columns; cells may be single years or delimited lists (e.g., '2020|2021')."
  )
}

# Build named list: raster -> unique sorted integer vector of years
raster_years_map <- meta_norm %>%
  mutate(years = map(years, ~ unique(sort(as.integer(.x))))) %>%
  deframe()

# Warn if any present raster has no metadata years
rasters_missing_years <- setdiff(present_rasters, names(raster_years_map))
if (length(rasters_missing_years) > 0) {
  warning("These rasters are present in points but have no years in metadata (they will be excluded): ",
          paste(rasters_missing_years, collapse = ", "))
}

# ----------------------------
# 6) Reshape points to long and attach allowed years (list-column)
# ----------------------------
# Convert raster columns to integer 0/1/2/NA. Anything not 0/1/2 becomes NA.
points_norm <- points %>%
  mutate(
    across(
      all_of(present_rasters),
      ~ {
        x <- as.character(.)
        x <- str_trim(x)
        x[x %in% c("", "NA", "NaN", "null", "NULL")] <- NA_character_
        ifelse(x %in% c("0","1","2"), as.integer(x), NA_integer_)
      }
    )
  )

# Long format: one row per point × raster
# Also attach:
#   - raster_group (pli_* -> pli)
#   - allowed_years (list of integers from metadata)
long <- points_norm %>%
  select(all_of(required_cols), all_of(present_rasters)) %>%
  pivot_longer(
    cols = all_of(present_rasters),
    names_to = "raster",
    values_to = "pred"
  ) %>%
  mutate(
    raster = tolower(raster),
    raster_group = raster_to_group(raster),
    allowed_years = map(raster, ~ raster_years_map[[.x]] %||% integer(0))
  )

# Collapse lists to strings and export
long_out <- long %>%
  mutate(allowed_years = map_chr(allowed_years, ~ paste(.x, collapse = ";")))

write.csv(long_out, file.path(out_dir, "long_prefiltering_26Jun2026_natvtam_withnongrass_compiled.csv"))

# ----------------------------
# 7) Filter rows via class-based time rules using MULTIPLE allowed years
# ----------------------------
# Set up class-dependent year matching rules:
#   Native grass (lc3 == 2L): accept any year ≤ reference year (choose most recent such year)
#   Tame grass (lc3 == 1L): accept only if year matches exactly
#   Non-grass (lc3 == 0L): can configure rule to be like native or tame rule 
#compute_matched_year_class3 <- function(yrs, yer, lc3, notgrass_rule = c("exact","lte")) {
#  notgrass_rule <- match.arg(notgrass_rule)
#  yrs <- as.integer(yrs)
#  if (length(yrs) == 0L || is.na(yer) || is.na(lc3)) return(NA_integer_)
#  if (lc3 == 2L) {                # native
#    ok <- yrs[yrs <= yer]
#    if (length(ok) == 0L) return(NA_integer_)
#    return(max(ok))
#  } else if (lc3 == 1L) {         # tame
#    if (yer %in% yrs) return(yer) else return(NA_integer_)
#  } else if (lc3 == 0L) {         # not grass
#    if (notgrass_rule == "exact") {
#      if (yer %in% yrs) return(yer) else return(NA_integer_)
#    } else {"lte"
#      ok <- yrs[yrs <= yer]
#      if (length(ok) == 0L) return(NA_integer_)
#      return(max(ok))
#    }
#  } else {
#    return(NA_integer_)
#  }
#}

compute_matched_year_class3 <- function(yrs, yer, lc3) {
  yrs <- as.integer(unlist(yrs))
  if (length(yrs) == 0L || is.na(yer) || is.na(lc3)) return(NA_integer_)
  if (lc3 == 2L) { # native
    ok <- yrs[yrs <= yer]
    if (length(ok) == 0L) return(NA_integer_)
    return(max(ok))
  } else if (lc3 == 1L) { # tame
    if (yer %in% yrs) return(yer) else return(NA_integer_)
  } else if (lc3 == 0L) { # nongrass
    if (yer %in% yrs) return(yer) else return(NA_integer_)
  } else {
    return(NA_integer_)
  }
}

# Apply filtering:
#   - require LC_class3 ∈ {0,1,2}, Yer_Cll present, pred present
#   - DO NOT filter out LC_class3 == 0
#   - DO NOT exclude LC_Cl_N == "SHRUB"
#   - require allowed_years matching via matched_year
#   - Exclude training sources (updated filtering)
eval_data <- long %>%
  filter(
    !is.na(LC_class3),
    LC_class3 %in% c(0L, 1L, 2L),  # keep all 3 classes
    !is.na(Yer_Cll)) %>%
  filter(!is.na(pred)) %>%
  rowwise() %>%
  mutate(matched_year = compute_matched_year_class3(allowed_years, Yer_Cll, LC_class3)) %>%
  ungroup() %>%
  filter(!is.na(matched_year)) %>%
  filter(is.na(Training) | (Training != raster & Training != raster_group)) # updated

# Collapse lists to strings and export
eval_data_out <- eval_data %>%
  mutate(allowed_years = map_chr(allowed_years, ~ paste(.x, collapse = ";")))

write.csv(eval_data_out, file.path(out_dir, "eval_data_updatedtraindatfilt_23Jun2026_natvtam_withnongrass_compiled.csv"))

# Isolate datasets that need to be handled differently
#   Datasets with training data need to be filtered separately
#   so that data is not lost

# ACI
#eval_data_aci <- eval_data %>%
#  filter(raster_group == "aci")

# NGI
#eval_data_ngi <- eval_data %>%
#  filter(raster_group == "ngi")

# PLI
#eval_data_pli <- eval_data %>%
#  filter(raster_group == "pli")

# Filter out training data
#   Remove points used to train the model

# ACI
#eval_data_aci_filt <- eval_data_aci %>%
#  filter(is.na(Source) | Source != "aci")

# HFI
#eval_data_hfi_filt <- eval_data_hfi %>%
#  filter(is.na(Source) | Source != "hfi")

# NGI
#eval_data_ngi_filt <- eval_data_ngi %>%
#  filter(is.na(Source) | Source != "ngi")

# PLI
#eval_data_pli_filt <- eval_data_pli %>%
#  filter(is.na(Source) | Source != "pli")

# Filter ACI, HFI, NGI, and PLI out of eval_data
#eval_data_reg <- eval_data %>%
#  filter(!(raster_group %in% c("aci", "pli", "ngi", "hfi")))

# Combine back together
#eval_data_combined <- bind_rows(
#  eval_data_reg,
#  eval_data_aci_filt,
#  eval_data_hfi_filt,
#  eval_data_ngi_filt,
#  eval_data_pli_filt)

# Check number of rows - 60,049 total
#sapply(
#  list(
#    table1 = eval_data_reg,
#    table2 = eval_data_aci_filt,
#    table3 = eval_data_hfi_filt,
#    table4 = eval_data_ngi_filt,
#    table5 = eval_data_pli_filt),
#  nrow)

# Collapse lists to strings and export
#eval_data_combined_out <- eval_data_combined %>%
#  mutate(allowed_years = map_chr(allowed_years, ~ paste(.x, collapse = ";")))

#write.csv(eval_data_combined_out, file.path(out_dir, "eval_data_29Apr2026_withnongrass_good.csv"))

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
  
  n <- sum(cm)
  if (n == 0) return(NA_real_)
  
  r <- rowSums(cm)
  c <- colSums(cm)
  
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
  
  n     <- sum(tab)
  diagv <- diag(tab)
  rows  <- rowSums(tab)   # support (truth)
  cols  <- colSums(tab)   # predicted totals
  
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
    N_Evaluated        = as.integer(n),
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
  eval_data,
  key_cols = "raster_group"
)

confusions_by_dataset <- confusion_block_multiclass(
  eval_data,
  key_cols = "raster_group"
)

metrics_by_dataset

# ----------------------------
# Export output tables
# ----------------------------
out_xlsx <- file.path(out_dir, "nat_v_tame_v_nongrass_accuracy_assessment_results_23Jun2026.xlsx")

write_xlsx(
  list(
    "Metrics_by_Dataset"   = metrics_by_dataset,
    "Confusion_Long_Form"  = confusions_by_dataset
  ),
  path = out_xlsx
)

