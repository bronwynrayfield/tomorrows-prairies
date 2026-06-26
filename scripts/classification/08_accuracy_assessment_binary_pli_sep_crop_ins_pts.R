
####################################################################################
##           PRELIMINARY ACCURACY ASSESSMENT: BINARY (GRASS VS NON-GRASS)         ##
####################################################################################

##----------------------- Alison Long - June 2026 -------------------------

# step 8: assess accuracy of datasets mapping grasslands by comparing them
#         to ground-truthed point data for 2015 to present

## THIS VERSION USES CROP INSURANCE DERIVED DATA FROM EMILY LINDSAY 

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
##   - class-based time rule using allowed years for that raster:
##       NATIVE: any( allowed_year <= Yer_Cll ); matched_year = max(allowed_year ≤ Yer_Cll)
##       TAME:   any( allowed_year == Yer_Cll ); matched_year = Yer_Cll (if allowed)
##       OTHER:  any( allowed_year == Yer_Cll ); matched_year = Yer_Cll (if allowed)
##
## Exclusions:
##   - Source used for training: exclude if Source_norm matches either atomic raster code or grouped alias
##
## Outputs:
##   - Grouped (legacy filenames): PLI variants combined under "pli"
##   - Atomic (new filenames with _atomic): each PLI variant separate (pli_ap, pli_cu, pli_mg, pli_mm)

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
#points_excel_path <- "R:/CanadianGrasslands/WG_Classification/processed/groundtruthed_data_points/ci_derived/point_raster_extract_grass_v_nongrass.xlsx"
points_excel_path <- "R:/CanadianGrasslands/WG_Classification/processed/groundtruthed_data_points/ci_derived/point_raster_extract_v2_grass_v_nongrass.xlsx"
metadata_csv_path <- "R:/CanadianGrasslands/WG_Classification/github/tomorrows-prairies/processeddatasets/Grass/metadata_full_extent_binary_pli_sep_noplow.csv"

out_dir <- "R:/CanadianGrasslands/WG_Classification/results/accuracy_assessment/binary/crop_ins_derived"

# ----------------------------
# 3) Configuration
# ----------------------------
# Handle PLI as four separate rasters in the points file and metadata.
pli_variants <- c("pli_ap", "pli_cu", "pli_mg", "pli_mm")

# Expected raster columns 
raster_cols <- c(
  "aci2024","luts10","luts30","cec","ngi","cgr","eccc", "esa","glc",
  "cnvc","lcc","pudl",
  pli_variants,            # four atomic PLI rasters
  "lcr","gpw","plvi","gcc","dw","esri","mgi","pglc")

# Grouping helper: collapse pli_* -> "pli"; all others remain themselves
raster_to_group <- function(x) ifelse(x %in% pli_variants, "pli", x)

# ----------------------------
# 4) Load points and normalize basic fields
# ----------------------------
#points <- readxl::read_excel(points_excel_path)
#points <- readxl::read_excel(
#  points_excel_path,
#  col_types = c("guess", "text", rep("guess", 39))
#)

points <- readxl::read_excel(
  points_excel_path,
  col_types = c("numeric", "text", "numeric", rep("guess", 4), rep("numeric", 26))
)

names(points)
str(points)
unique(points$LC_Class_N)
#unique(points$Source)
#table(points$Source, useNA = "ifany")

# Rename columns
points <- points %>%
  rename(
    Pt_ID_O    = Point_ID,
    Sorc_ID  = Source_ID,
    Yer_Cll  = Year_Colle,
    LC_Cl_N  = LC_Class_N,
    Ecds_ID = ecodist,
    Provinc  = PROV,
    Ecrgn_n  = ecoreg,
    aci2024  = aafc_aci_merged_reproj_algn_binary_noshrub,
    luts10   = luts10_reproj_algn_binary_masked,
    luts30   = luts30_reproj_algn_binary_masked,
    cec      = nalcms_reproj_algn_binary_noshrub_masked,
    ngi      = cfga_ngi_merged_reproj_algn_binary,
    cgr      = cgr_reproj_algn_binary_masked, 
    eccc     = eccc_reclass_reproj_algn_binary_masked,
    esa      = esa_reproj_algn_binary_noshrub_masked,
    glc      = glc_reproj_algn_binary_noshrub_masked,
    cnvc     = cnvc_algn_binary_masked,
    lcc      = lcc_reproj_algn_binary_noshrub_masked,
    pudl     = pudl_reclass_reproj_algn_binary_noshrub, 
    pli_ap   = sk_pli_ap_reproj_algn_binary_noshrub,
    pli_cu   = sk_pli_cu_reproj_algn_binary_noshrub,
    pli_mg   = sk_pli_mg_reproj_algn_binary_noshrub,
    pli_mm   = sk_pli_mm_reproj_algn_binary_noshrub,
    lcr      = lcr_reproj_algn_binary_masked,
    gpw      = gpw_reproj_algn_binary_noshrub_masked,
    plvi     = plvi_reproj_algn_binary_gte2015,
    hfi      = hfi_reproj_algn_binary_gte2015,
    gcc      = gcc_reproj_algn_binary,
    dw       = dw_reproj_algn_masked_binary,
    esri     = esri_reproj_algn_masked_binary,
    mgi      = mb_mgi_reproj_algn_binary_noshrub,
    pglc     = pglc_reproj_algn_binary,
    plow     = plow_reproj_algn_intact_masked)

str(points)

# Remove Plowprint and HFI
points <- points %>%
  select(-plow) %>%
  select(-hfi)

# Create column for training data exclusion
points <- points %>%
  mutate(
    Training = case_when(
      Sorc_ID == "ACI_field_2024" ~ "aci2024", # if using 2024 ACI data
      TRUE ~ ""))    # keep unmatched cases empty
    
# Create column with numerical class for ground-truthed cover type
points <- points %>%
  mutate(
    LC_grass = case_when(
      LC_Cl_N %in% c("TAME", "NATIVE") ~ 1,
      LC_Cl_N == "OTHER" ~ 0,
      TRUE ~ NA_real_))

unique(points$LC_grass)
table(points$LC_grass, useNA = "ifany")

# Check for NAs in Pt_ID
points %>%
  summarise(n_na = sum(is.na(Pt_ID_O)))

# Create new Pt_ID column with no NAs
points <- points %>%
  mutate(Pt_ID = row_number()) %>%
  relocate(Pt_ID)

# Drop rows where ecodistrict ID is 0
points <- points %>%
  filter(Ecds_ID != 0)

# Define rquired columns
required_cols <- c("Pt_ID", "Ecds_ID", "Ecrgn_n", "Training", "Sorc_ID", "LC_grass", "Yer_Cll", "LC_Cl_N")
missing_req <- setdiff(required_cols, names(points))
if (length(missing_req) > 0) {
  stop("Missing required columns in points file: ", paste(missing_req, collapse = ", "))
}

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

# Normalize LC_grass to integer 0/1, Yer_Cll to year int, and LC_Cl_N to uppercased labels
#points <- points %>%
#  mutate(
#    LC_grass = suppressWarnings(as.integer(as.numeric(LC_grass))),
#    Yer_Cll  = suppressWarnings(as.integer(as.numeric(Yer_Cll))),
#    LC_Cl_N  = {
#      z <- as.character(LC_Cl_N)
#      z <- stringr::str_trim(toupper(z))
#      z[z %in% c("", "NA", "NAN", "NULL", "MISSING")] <- NA_character_
#      z
#    }
#  )

# Normalize Source (if present) for training exclusion
#has_source <- "Source" %in% names(points)
#if (has_source) {
#  points <- points %>%
#    mutate(
#      Source_norm = {
#        s <- as.character(Source)
#        s <- stringr::str_trim(tolower(s))
#        s[s %in% c("", "na", "nan", "null", "missing")] <- NA_character_
#        s
#      }
#    )
#} else {
#  points <- points %>% mutate(Source_norm = NA_character_)
#}

# Identify which raster columns actually exist in the points file
present_rasters <- intersect(raster_cols, names(points))
if (length(present_rasters) == 0) {
  stop("None of the expected raster columns were found in the points file.")
}

# ----------------------------
# 5) Load metadata and allow MULTIPLE years per raster
# ----------------------------
# The metadata CSV can be:
#   (A) Long: columns [raster, year], with multiple rows allowed per raster
#   (B) Wide 1-row: columns named by rasters; cells may contain single year or delimited list like "2020|2021"
#
# Parse into a named LIST: raster -> integer vector of allowed years.

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
  # Long format: may have multiple rows per raster
  meta_norm <- meta_raw %>%
    transmute(
      raster = str_trim(tolower(as.character(raster))),
      year_raw = year
    ) %>%
    group_by(raster) %>%
    summarise(years = list(parse_year_vector(year_raw)), .groups = "drop")
} else if (nrow(meta_raw) == 1 && any(intersect(raster_cols, names(meta_raw)) %in% names(meta_raw))) {
  # Wide single-row: one column per raster, possibly with "2020|2021"
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
# Convert raster columns to integer 0/1/NA. Anything not 0 or 1 becomes NA.
points_norm <- points %>%
  mutate(
    across(
      all_of(present_rasters),
      ~ {
        x <- as.character(.)
        x <- str_trim(x)
        x[x %in% c("", "NA", "NaN", "null", "NULL")] <- NA_character_
        ifelse(x %in% c("0", "1"), as.integer(x), NA_integer_)
      }
    )
  )

# Long format: one row per point × raster
# Also attach:
#   - raster_group (pli_* -> pli)
#   - allowed_years (list of integers from metadata)
long <- points_norm %>%
  select(all_of(required_cols), all_of(present_rasters)) %>% # Source_norm, 
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

# ----------------------------
# 7) Filter rows via class-based time rules using MULTIPLE allowed years
# ----------------------------
# Compute a per-record matched_year for per-year reporting:
#   - NATIVE: latest allowed_year <= Yer_Cll (if any; else NA -> exclude)
#   - TAME: Yer_Cll if Yer_Cll is in allowed_years (else NA -> exclude)
#   - OTHER: Yer_Cll if Yer_Cll is in allowed_years (else NA -> exclude)
compute_matched_year <- function(yrs, yer, class_label) {
  yrs <- as.integer(unlist(yrs))
  if (length(yrs) == 0L || is.na(yer) || is.na(class_label)) return(NA_integer_)
  if (class_label == "NATIVE") {
    ok <- yrs[yrs <= yer]
    if (length(ok) == 0L) return(NA_integer_)
    return(max(ok))
  } else if (class_label == "TAME") {
    if (yer %in% yrs) return(yer) else return(NA_integer_)
  } else if (class_label == "OTHER") {
    if (yer %in% yrs) return(yer) else return(NA_integer_)
  } else {
    return(NA_integer_)
  }
}

# Apply full filtering:
#   - require LC_grass, Yer_Cll, LC_Cl_N
#   - require pred present and allowed_years matching via matched_year
#   - Exclude training sources (updated filtering)
eval_data <- long %>%
  filter(!is.na(LC_grass), !is.na(Yer_Cll), !is.na(LC_Cl_N)) %>%
  filter(!is.na(pred)) %>%
  rowwise() %>%
  mutate(matched_year = compute_matched_year(allowed_years, Yer_Cll, LC_Cl_N)) %>%
  ungroup() %>%
  filter(!is.na(matched_year)) %>%
  filter(is.na(Training) | (Training != raster & Training != raster_group)) # updated

# Add column for confusion matrix quadrants
eval_data <- eval_data %>%
  mutate(quad = case_when(
    LC_grass == 1 & pred == 1 ~ "tp",
    LC_grass == 0 & pred == 0 ~ "tn",
    LC_grass == 0 & pred == 1 ~ "fp",
    LC_grass == 1 & pred == 0 ~ "fn"
  ))

# Collapse lists to strings and export
eval_data_out <- eval_data %>%
  mutate(allowed_years = map_chr(allowed_years, ~ paste(.x, collapse = ";")))

readr::write_csv(eval_data_out, file.path(out_dir, "accuracy_assessment_data_binary_pli_sep_cropinsV2_26Jun2026.csv"))

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
  eval_data,
  key_cols = c("raster", "matched_year")
) %>%
  rename(raster_year = matched_year) %>%
  arrange(desc(Overall_Accuracy))

# GROUPED: collapse pli_* -> "pli"; combine counts across all matched_year within group
summary_by_raster_grouped <- eval_data %>%
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
  eval_data,
  key_cols = c("raster", "matched_year", "Yer_Cll")
) %>%
  rename(raster_year = matched_year) %>%
  arrange(raster, Yer_Cll)

# GROUPED: per grouped raster × Yer_Cll (mix of multiple matched years inside)
per_year_summary_grouped <- summarize_block(
  eval_data %>% mutate(raster = raster_group),
  key_cols = c("raster", "Yer_Cll")
) %>%
  mutate(raster_year = NA_integer_) %>%
  relocate(raster, raster_year) %>%
  arrange(raster, Yer_Cll)

# ----------------------------
# 12) Save outputs (CSV)
# ----------------------------
# GROUPED keep legacy filenames
readr::write_csv(summary_by_raster_grouped, file.path(out_dir, "raster_accuracy_summary.csv"))
readr::write_csv(confusion_counts_grouped,  file.path(out_dir, "raster_confusion_matrices.csv"))
readr::write_csv(per_year_summary_grouped,  file.path(out_dir, "raster_per_year_accuracy.csv"))

# ATOMIC with _atomic suffix
readr::write_csv(summary_by_raster_atomic, file.path(out_dir, "raster_accuracy_summary_atomic.csv"))
readr::write_csv(confusion_counts_atomic,  file.path(out_dir, "raster_confusion_matrices_atomic.csv"))
readr::write_csv(per_year_summary_atomic,  file.path(out_dir, "raster_per_year_accuracy_atomic.csv"))

# ----------------------------
# 13) Accuracy by Ecodistrict
# ----------------------------

# Ensure ecodistrict column exists
if (!"Ecds_ID" %in% names(eval_data)) {
  stop("Ecds_ID column not found in eval_data.")
}

# ---- ATOMIC: raster × ecodistrict ----
accuracy_by_ecodistrict_atomic <- summarize_block(
  eval_data,
  key_cols = c("raster", "Ecds_ID")
) %>%
  arrange(raster, Ecds_ID)

# ---- GROUPED: collapse pli_* -> pli ----
accuracy_by_ecodistrict_grouped <- eval_data %>%
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
  eval_data,
  key_cols = c("raster", "matched_year", "Ecds_ID")
) %>%
  rename(raster_year = matched_year) %>%
  arrange(raster, raster_year, Ecds_ID)

#------ Save CSVs -----

# GROUPED
readr::write_csv(
  accuracy_by_ecodistrict_grouped,
  file.path(out_dir, "accuracy_by_ecodistrict_grouped.csv")
)

# ATOMIC
readr::write_csv(
  accuracy_by_ecodistrict_atomic,
  file.path(out_dir, "accuracy_by_ecodistrict_atomic.csv")
)

# ATOMIC + YEAR (detailed output)
readr::write_csv(
  accuracy_by_ecodistrict_year_atomic,
  file.path(out_dir, "accuracy_by_ecodistrict_year_atomic.csv")
)

# ----------------------------
# 14) Create ecoregion-level fallbacks for ecodistricts with small N
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

#---- Join ecodistrict and ecoregion results ----

# Get mapping (one row per ecodistrict)
ecodist_to_ecoreg <- eval_data %>%
  distinct(Ecds_ID, Ecrgn_n)

accuracy_by_ecodistrict_atomic2 <- accuracy_by_ecodistrict_atomic %>%
  left_join(ecodist_to_ecoreg, by = "Ecds_ID")

# Join
fallback_joined <- accuracy_by_ecodistrict_atomic2 %>%
  left_join(
    accuracy_by_ecoregion,
    by = c("raster", "Ecrgn_n"),
    suffix = c("_ecoD", "_ecoR")
  )

#---- Apply fallback logic ----
# Replace metrics when N_evaluated is too small
accuracy_with_fallback <- fallback_joined %>%
  mutate(
    use_fallback = N_evaluated_ecoD < min_n,
    
    TP = ifelse(use_fallback, TP_ecoR, TP_ecoD),
    TN = ifelse(use_fallback, TN_ecoR, TN_ecoD),
    FP = ifelse(use_fallback, FP_ecoR, FP_ecoD),
    FN = ifelse(use_fallback, FN_ecoR, FN_ecoD),
    N_evaluated = ifelse(use_fallback, N_evaluated_ecoR, N_evaluated_ecoD),
    
    # recompute metrics from chosen counts (important!)
    metrics = pmap(
      list(TP, TN, FP, FN),
      ~ compute_metrics(..1, ..2, ..3, ..4)
    )
  ) %>%
  tidyr::unnest_wider(metrics) %>%
  select(raster, Ecds_ID, Ecrgn_n, use_fallback, everything())

#---- Keep provenance flags ----
accuracy_with_fallback <- accuracy_with_fallback %>%
  mutate(
    source_level = ifelse(use_fallback, "ecoregion", "ecodistrict")
  )

# ----------------------------
# 15) Create rasters of accuracy by ecodistrict
# ----------------------------

summary_by_ecodistrict <- summarize_block(
  eval_data,
  key_cols = c("raster", "Ecds_ID"))

summary_by_ecodistrict %>%
  group_by(raster) %>%
  summarise(min_n = min(N_evaluated),
            mean_n = mean(N_evaluated),
            max_n = max(N_evaluated))

# global MCC fallback
global_mcc <- summary_by_raster_grouped %>%
  select(raster, MCC_global = MCC)

# clean + stabilize
weights_ecodist_mcc <- summary_by_ecodistrict %>%
  select(raster, Ecds_ID, MCC, N_evaluated) %>%
  
  # remove unreliable small samples
  mutate(
    MCC = ifelse(N_evaluated < 30, NA, MCC)
  ) %>%
  
  # join global fallback
  left_join(global_mcc, by = "raster") %>%
  
  # shrink local MCC toward global MCC
  # Otherwise, small ecodistricts will have noisy weights
  mutate(
    MCC_final = case_when(
      is.na(MCC) ~ MCC_global,
      TRUE ~ 0.7 * MCC + 0.3 * MCC_global
    )
  ) %>%
  
  # scale MCC from [-1,1] to [0,1] so values are usable
  mutate(
    MCC_final = (MCC_final + 1) / 2
  ) %>%
  
  # clamp to [0,1] just in case
  mutate(
    MCC_final = pmin(pmax(MCC_final, 0), 1)
  ) %>%

  select(raster, Ecds_ID, MCC_final)

# Validate ecodistrict coverage
weights_ecodist_mcc %>%
  group_by(raster) %>%
  summarise(n_Ecds = n_distinct(Ecds_ID))

# Build raster stack
library(terra)

# Inputs
eco_vect <- vect("R:/CanadianGrasslands/processeddatasets/boundaries/national-ecological-framework-for-canada/Ecodistricts/ecodistricts_clean.shp")
eco_vect$Ecds_ID <- as.character(eco_vect$ECODISTRIC)

weights_ecodist_mcc$Ecds_ID <- as.character(weights_ecodist_mcc$Ecds_ID)

template <- rast("R:/CanadianGrasslands/WG_Classification/processed/aafc_aci/aligned/aafc_aci_merged_reproj_algn_binary_noshrub.tif")  # ensures perfect alignment

weights_ecodist_mcc <- weights_ecodist_mcc %>%
  mutate(raster = tolower(raster))

raster_names <- unique(weights_ecodist_mcc$raster)

# output container
mcc_rasters <- vector("list", length(raster_names))
names(mcc_rasters) <- raster_names

# Loop
for (r_name in raster_names) {
  
  # extract weights for this raster
  w_df <- weights_ecodist_mcc %>%
    filter(raster == r_name)
  
  # Check for empty weights
  if (nrow(w_df) == 0) {
    message("WARNING: No weights found for ", r_name)
  }
  
  # join to polygons
  eco_joined <- merge(
    eco_vect,
    w_df,
    by = "Ecds_ID",
    all.x = TRUE
  )
  
  # fallback value (global MCC)
  fallback_val <- global_mcc$MCC_global[
    global_mcc$raster == r_name
  ]
  
  # Scale fallback like MCC
  fallback_val <- (fallback_val + 1) / 2
  
  eco_joined$MCC_final[is.na(eco_joined$MCC_final)] <- fallback_val
  
  # rasterize
  r_mcc <- rasterize(
    eco_joined,
    template,
    field = "MCC_final",
    touches = TRUE
  )
  
  names(r_mcc) <- r_name
  
  mcc_rasters[[r_name]] <- r_mcc
}

# Stack
mcc_stack <- rast(mcc_rasters)

# enforce exact layer order match with bin_stack
#mcc_stack <- mcc_stack[[names(bin_stack)]]

# Save stack
mcc_file <- file.path(out_dir, "mcc_weight_stack_bin.tif")

writeRaster(
  mcc_stack,
  filename = mcc_file,
  overwrite = TRUE,
  wopt = list(
    datatype = "FLT4S",
    gdal = c("COMPRESS=ZSTD", "TILED=YES", "BIGTIFF=YES")
  )
)

# ----------------------------
# 16) Plots: Map by ecodistrict
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
# 17) Plots: Overall accuracy (GROUPED and ATOMIC)
# ----------------------------
# GROUPED (legacy) plot
p_grouped <- ggplot(summary_by_raster_grouped,
                    aes(x = reorder(raster, Overall_Accuracy), y = Overall_Accuracy)) +
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

ggsave(file.path(out_dir, "raster_accuracy_plot_binary_cropsins.png"),
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

ggsave(file.path(out_dir, "raster_accuracy_plot_atomic_binary_cropins.png"),
       p_atomic, width = 7, height = 5, dpi = 200)

# ----------------------------
# 18) Console previews
# ----------------------------
cat("\n=== Summary by Raster (GROUPED) ===\n")
print(summary_by_raster_grouped, n = Inf)

cat("\n=== Summary by Raster (ATOMIC) ===\n")
print(summary_by_raster_atomic, n = Inf)

cat("\n=== Confusion Counts (GROUPED) ===\n")
print(confusion_counts_grouped, n = Inf)

cat("\n=== Confusion Counts (ATOMIC) ===\n")
print(confusion_counts_atomic, n = Inf)

#################################################################
#################################################################

#------ Precision vs recall plot -------------
# Evaluate if a dataset is missing grass (low recall) or overpredicting grass (low precision)
# Interpretation:
#   Top-right = best (high precision + recall)
#   Below diagonal → overpredicting grass
#   Above diagonal → conservative (missing grass)
p_pr <- ggplot(summary_by_raster_grouped,
       aes(x = Recall_grass, y = Precision_grass, label = raster)) +
  geom_point(size = 3, color = "#2d7f5e") +
  geom_text(nudge_y = 0.02, size = 3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray60") +
  xlim(0, 1) + ylim(0, 1) +
  labs(title = "Precision vs Recall (Grass)",
       x = "Recall (Producer’s Accuracy)",
       y = "Precision (User’s Accuracy)") +
  theme_minimal()

ggsave(file.path(out_dir, "precision_vs_recall_plot_binary_cropins.png"),
       p_pr, width = 7, height = 5, dpi = 300)

#------- False pos vs false neg rate (bias plot) ---------
# Show error type bias
# FP = mapping non-grass as grass
# FN = missing grass (often more critical in conservation)
p_bias <- ggplot(summary_by_raster_grouped,
       aes(x = False_Positive_Rate, y = False_Negative_Rate, label = raster)) +
  geom_point(size = 3, color = "#a6611a") +
  geom_text(nudge_y = 0.02, size = 3) +
  xlim(0, 1) + ylim(0, 1) +
  labs(title = "Error Tradeoff: FP vs FN",
       x = "False Positive Rate (overpredict grass)",
       y = "False Negative Rate (miss grass)") +
  theme_minimal()

ggsave(file.path(out_dir, "error_bias_plot_binary_cropins.png"),
       p_bias, width = 7, height = 5, dpi = 300)

#-------- MCC vs prevalence -------------
# Robust performance check since MCC shows overall true skill and
# prevalnce shows how imbalanced the evaluation data is
p_mp <- ggplot(summary_by_raster_grouped,
       aes(x = Prevalence_grass, y = MCC, label = raster)) +
  geom_point(size = 3, color = "#018571") +
  geom_text(nudge_y = 0.02, size = 3) +
  xlim(0, 1) +
  labs(title = "MCC vs Grass Prevalence",
       x = "Grass Prevalence in Samples",
       y = "Matthews Correlation Coefficient") +
  theme_minimal()

ggsave(file.path(out_dir, "mcc_vs_prevalence_plot_binary_cropins.png"),
       p_mp, width = 7, height = 5, dpi = 300)

#--------- Metric comparison panel -----------
summary_long <- summary_by_raster_grouped %>%
  select(raster, Overall_Accuracy, Balanced_Accuracy, F1_grass, MCC) %>%
  pivot_longer(-raster, names_to = "metric", values_to = "value")

p_compar <- ggplot(summary_long,
       aes(x = reorder(raster, value), y = value, fill = metric)) +
  geom_col(position = "dodge") +
  coord_flip() +
  labs(title = "Comparison of Key Metrics",
       x = "Raster", y = "Score") +
  theme_minimal()

ggsave(file.path(out_dir, "metric_comparison_plot_binary_cropins.png"),
       p_compar, width = 7, height = 5, dpi = 300)

#----------- Temporal performance stability
p_yr <- ggplot(per_year_summary_atomic,
       aes(x = Yer_Cll, y = Balanced_Accuracy, color = raster)) +
  geom_line() +
  geom_point() +
  labs(title = "Accuracy Over Time",
       x = "Ground Truth Year",
       y = "Balanced Accuracy") +
  theme_minimal()

ggsave(file.path(out_dir, "accuracy_by_year_plot_binary_cropins.png"),
       p_yr, width = 7, height = 5, dpi = 300)

#---- Sample size vs performance
# Reliability check 
p_n <- ggplot(summary_by_raster_atomic,
       aes(x = N_evaluated, y = MCC, label = raster)) +
  geom_point(size = 3, color = "#5e3c99") +
  geom_text(nudge_y = 0.02, size = 3) +
  scale_x_log10() +
  labs(title = "Performance vs Sample Size",
       x = "Sample Size (log scale)",
       y = "MCC") +
  theme_minimal()

ggsave(file.path(out_dir, "mcc_vs_N_plot_binary_cropins.png"),
       p_n, width = 7, height = 5, dpi = 300)
