
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
points_excel_path <- "R:/CanadianGrasslands/WG_Classification/processed/groundtruthed_data_points/groundtruth_points_with_grass_raster_values_natvtam_pli_sep_ge2015_18Jul2026_clean_withnongrass.xlsx"
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
  col_types = c("numeric", "text", "text", "numeric", rep("guess", 10), rep("numeric", 20))
)

names(points)
str(points)
unique(points$Training)
table(points$Training, useNA = "ifany")
unique(points$LC_class3)
table(points$LC_class3, useNA = "ifany")

# Remove HFI (and ACI 2025 for now)
points <- points %>%
  select(-hfi) %>%
  select(-aci2025)

# Rename aci2024 column to aci (temporarily)
points <- points %>%
  dplyr::rename(aci = aci2024)

# Add point ID column
points <- points %>%
  mutate(Pt_ID = row_number()) %>%
  relocate(Pt_ID)

# Drop rows where ecodistrict ID is 0
points <- points %>%
  filter(Ecds_ID != 0)

unique(points$Ecrgn_n)

# Create column for ecoprovince
points <- points %>%
  mutate(
    Ecprv = case_when(
      Ecrgn_n %in% c("Abitibi Plains", 
                     "Lake Nipigon",
                     "Rupert River Plateau",
                     "Lac Seul Upland") ~ "6.2",
      Ecrgn_n %in% c("Thunder Bay-Quetico", 
                     "Lake of the Woods",
                     "Rainy River") ~ "6.5",
      Ecrgn_n %in% c("Algonquin-Lake Nipissing", 
                     "Southern Laurentians",
                     "Lac Temiscamingue Lowland") ~ "6.6",
      Ecrgn_n %in% c("St-Laurent Lowlands", 
                     "St. Lawrence Lowlands", 
                     "Manitoulin-Lake Simcoe",
                     "Frontenac Axis") ~ "8.1",
      Ecrgn_n == "Lake Erie Lowland" ~ "8.2",
      Ecrgn_n %in% c("Western Alberta Upland", 
                     "Clear Hills Upland") ~ "9.1",
      Ecrgn_n %in% c("Wabasca Lowland", 
                     "Boreal Transition",
                     "Western Boreal",
                     "Peace Lowland",
                     "Mid-Boreal Uplands") ~ "9.2",
      Ecrgn_n %in% c("Interlake Plain", 
                     "Mid-Boreal Lowland") ~ "9.3",
      Ecrgn_n == "Lake Manitoba Plain" ~ "10.1",
      Ecrgn_n %in% c("Aspen Parkland", 
                     "Southwest Manitoba Uplands",
                     "Southwest Manotoba Uplands") ~ "10.2",
      Ecrgn_n %in% c("Mixed Grassland", 
                     "Moist Mixed Grassland",
                     "Fescue Grassland",
                     "Cypress Upland") ~ "10.3",
      Ecrgn_n %in% c("Lower Mainland", 
                     "Eastern Vancouver Island") ~ "13.1",
      Ecrgn_n %in% c("Pacific Ranges", 
                     "Cascade Ranges") ~ "13.2",
      Ecrgn_n %in% c("Central Canadian Rocky Mountains", 
                     "Skeena Mountains",
                     "Fraser Basin") ~ "14.1",
      Ecrgn_n %in% c("Fraser Plateau",
                     "Bulkley Ranges") ~ "14.2",
      Ecrgn_n %in% c("Thompson-Okanagan Plateau", 
                     "Okanogan Highland",
                     "Okanogan Range",
                     "Okanagan Highland",
                     "Okanagan Range",
                     "Interior Transition Ranges") ~ "14.3",
      Ecrgn_n %in% c("Columbia Mountains and Highlands", 
                     "Southern Rocky Mountain Trench",
                     "Northern Continental Divide",
                     "Eastern Continental Ranges",
                     "Western Continental Ranges",
                     "Selkirk-Bitterroot Foothills") ~ "14.4"))

unique(points$Ecprv)

# Define required columns
required_cols <- c("Pt_ID", "LC_class3", "Training","Sorc_ID", "Ecds_ID", "Ecrg_ID", "Ecrgn_n", "Ecprv", "Yer_Cll", "LC_Cl_N")
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

# Rename aci2024 column in metadata to aci (only when excluding aci2025)
meta_norm <- meta_norm %>%
  mutate(raster = if_else(raster == "aci2024", "aci", raster))

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

write.csv(long_out, file.path(out_dir, "long_prefiltering_18Jul2026_natvtam_withnongrass_compiled.csv"))

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

write.csv(eval_data_out, file.path(out_dir, "eval_data_18Jul2026_natvtam_withnongrass_compiled.csv"))

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
out_xlsx <- file.path(out_dir, "nat_v_tame_v_nongrass_accuracy_assessment_results_18Jul2026.xlsx")

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

if (!"Ecds_ID" %in% names(eval_data)) {
  stop("Ecds_ID column not found in eval_data.")
}

# ---- ATOMIC: raster × ecodistrict ----
accuracy_by_ecodistrict_atomic <- summarize_block_multiclass(
  eval_data,
  key_cols = c("raster", "Ecds_ID")
) %>%
  arrange(raster, Ecds_ID)

# ---- GROUPED: collapse pli_* -> pli ----
accuracy_by_ecodistrict_grouped <- eval_data %>%
  mutate(raster = raster_group) %>%
  summarize_block_multiclass(
    key_cols = c("raster", "Ecds_ID")
  ) %>%
  arrange(raster, Ecds_ID)

# ---- OPTIONAL: atomic with year ----
accuracy_by_ecodistrict_year_atomic <- summarize_block_multiclass(
  eval_data,
  key_cols = c("raster", "matched_year", "Ecds_ID")
) %>%
  rename(raster_year = matched_year) %>%
  arrange(raster, raster_year, Ecds_ID)

# ---- Save ----
readr::write_csv(accuracy_by_ecodistrict_grouped,
                 file.path(out_dir, "accuracy_by_ecodistrict_grouped_compiled.csv"))

readr::write_csv(accuracy_by_ecodistrict_atomic,
                 file.path(out_dir, "accuracy_by_ecodistrict_atomic_compiled.csv"))

readr::write_csv(accuracy_by_ecodistrict_year_atomic,
                 file.path(out_dir, "accuracy_by_ecodistrict_year_atomic_compiled.csv"))

# ----------------------------
# 14) Multiclass fallback hierarchy
# ecodistrict -> ecoregion -> ecoprovince
# ----------------------------

min_n <- 20

# ---- Ecoregion metrics ----
accuracy_by_ecoregion <- summarize_block_multiclass(
  eval_data,
  key_cols = c("raster", "Ecrgn_n")
)

# ---- Ecoprovince metrics ----
accuracy_by_ecoprovince <- summarize_block_multiclass(
  eval_data,
  key_cols = c("raster", "Ecprv")
) %>%
  rename_with(
    ~ paste0(.x, "_ecoP"),
    -c(raster, Ecprv)
  )

# ---- Mapping table ----
ecodist_to_regions <- eval_data %>%
  distinct(Ecds_ID, Ecrgn_n, Ecprv)

# --------------------------------------------------
# ATOMIC
# --------------------------------------------------

accuracy_by_ecodistrict_atomic2 <- accuracy_by_ecodistrict_atomic %>%
  left_join(ecodist_to_regions, by = "Ecds_ID")

fallback_joined_atomic <- accuracy_by_ecodistrict_atomic2 %>%
  left_join(
    accuracy_by_ecoregion,
    by = c("raster", "Ecrgn_n"),
    suffix = c("_ecoD", "_ecoR")
  ) %>%
  left_join(
    accuracy_by_ecoprovince,
    by = c("raster", "Ecprv")
  )

# Make sure suffixes are present 
names(fallback_joined_atomic)

metric_cols <- setdiff(
  names(accuracy_by_ecodistrict_atomic),
  c("raster", "Ecds_ID")
)

all(paste0(metric_cols, "_ecoD") %in% names(fallback_joined_atomic))
all(paste0(metric_cols, "_ecoR") %in% names(fallback_joined_atomic))
all(paste0(metric_cols, "_ecoP") %in% names(fallback_joined_atomic))

accuracy_with_fallback_atomic <- fallback_joined_atomic %>%
  mutate(
    source_level = case_when(
      N_Evaluated_ecoD >= min_n ~ "ecodistrict",
      N_Evaluated_ecoR >= min_n ~ "ecoregion",
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

#accuracy_with_fallback_atomic <- fallback_joined_atomic %>%
#  mutate(
#    source_level = case_when(
#      N_Evaluated_ecoD >= min_n ~ "ecodistrict",
#      N_Evaluated_ecoR >= min_n ~ "ecoregion",
#      TRUE ~ "ecoprovince"
#    )
#  ) %>%
#  mutate(
#    across(
#      all_of(metric_cols),
#      ~ case_when(
#        source_level == "ecodistrict" ~ .data[[paste0(cur_column(), "_ecoD")]],
#        source_level == "ecoregion" ~ .data[[paste0(cur_column(), "_ecoR")]],
#        TRUE ~ .data[[paste0(cur_column(), "_ecoP")]]
#      ),
#      .names = "{.col}"
#    )
#  )

# --------------------------------------------------
# GROUPED
# --------------------------------------------------

accuracy_by_ecodistrict_grouped2 <- accuracy_by_ecodistrict_grouped %>%
  left_join(ecodist_to_regions, by = "Ecds_ID")

fallback_joined_grouped <- accuracy_by_ecodistrict_grouped2 %>%
  left_join(
    accuracy_by_ecoregion,
    by = c("raster", "Ecrgn_n"),
    suffix = c("_ecoD", "_ecoR")
  ) %>%
  left_join(
    accuracy_by_ecoprovince,
    by = c("raster", "Ecprv")
  )

# Make sure suffixes are present 
names(fallback_joined_grouped)

#accuracy_with_fallback_grouped <- fallback_joined_grouped %>%
#  mutate(
#    source_level = case_when(
#      N_Evaluated_ecoD >= min_n ~ "ecodistrict",
#      N_Evaluated_ecoR >= min_n ~ "ecoregion",
#      TRUE ~ "ecoprovince"
#    )
#  ) %>%
#  mutate(
#    across(
#      all_of(metric_cols),
#      ~ case_when(
#        source_level == "ecodistrict" ~ .data[[paste0(cur_column(), "_ecoD")]],
#        source_level == "ecoregion" ~ .data[[paste0(cur_column(), "_ecoR")]],
#        TRUE ~ .data[[paste0(cur_column(), "_ecoP")]]
#      ),
#      .names = "{.col}"
#    )
#  )

accuracy_with_fallback_grouped <- fallback_joined_grouped %>%
  mutate(
    source_level = case_when(
      N_Evaluated_ecoD >= min_n ~ "ecodistrict",
      N_Evaluated_ecoR >= min_n ~ "ecoregion",
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
    "accuracy_by_ecodistrict_atomic_nattam_ecoregprovfallback_compiled_18Jul2026.csv"
  )
)

readr::write_csv(
  accuracy_with_fallback_grouped,
  file.path(
    out_dir,
    "accuracy_by_ecodistrict_grouped_nattam_ecoregprovfallback_compiled_18Jul2026.csv"
  )
)

# ----------------------------
# 15) Create rasters of accuracy by ecodistrict
# ----------------------------

summary_by_ecodistrict <- summarize_block_multiclass(
  eval_data,
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

eco_vect <- vect("R:/CanadianGrasslands/processeddatasets/boundaries/national-ecological-framework-for-canada/Ecodistricts/ecodistricts_clean.shp")
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
  eval_data,
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
