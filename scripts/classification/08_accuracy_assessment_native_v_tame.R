
#######################################################################################
## PRELIMINARY ACCURACY ASSESSMENT: NATIVE AND TAME CLASSES (NATIVE AND TAME POINTS) ##
#######################################################################################

##----------------------- Alison Long - March 2026 -------------------------

# step 6 (multiclass): assess accuracy of datasets mapping grasslands
# against ground-truthed point data (2015 to present), with 3 classes:
#   0 = not grass, 1 = tame grass, 2 = native grass
#   (filter out groundtruthed points that are classified as 0)
#
# Key changes from binary version:
#   - Uses LC_class3 (0/1/2) instead of LC_grass (0/1)
#   - Rasters expected as 0/1/2
#   - Time-aligned rules:
#         native (2): latest allowed_year <= Yer_Cll
#         tame   (1): allowed_year == Yer_Cll
#   - 3x3 confusion matrices and per-class metrics

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
points_excel_path <- "R:/CanadianGrasslands/WG_Classification/processed/groundtruthed_data_points/groundtruth_points_with_grass_raster_values_natvtam_pli_sep_ge2015_Mar2026_clean_v2.xlsx"
metadata_csv_path <- "R:/CanadianGrasslands/WG_Classification/github/tomorrows-prairies/processeddatasets/Grass/metadata_full_extent_nattam_pli_sep.csv"

out_dir <- "R:/CanadianGrasslands/WG_Classification/results/accuracy_assessment/native_v_tame"
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
                 "gpw","plvi","hfi")

# Grouping helper: collapse pli_* -> "pli"; all others remain themselves
raster_to_group <- function(x) ifelse(x %in% pli_variants, "pli", x)

# ----------------------------
# 4) Load points and normalize basic fields
# ----------------------------
points <- readxl::read_excel(points_excel_path)

required_cols <- c("LC_class3", "Yer_Cll")
missing_req <- setdiff(required_cols, names(points))
if (length(missing_req) > 0) {
  stop("Missing required columns in points file: ", paste(missing_req, collapse = ", "))
}

# Normalize LC_class3 to integer 0/1/2, Yer_Cll to year int
points <- points %>%
  mutate(
    LC_class3 = {
      x <- as.character(LC_class3)
      x <- str_trim(tolower(x))
      x[x %in% c("", "na", "nan", "null", "missing")] <- NA_character_
      # Accept words if present
      x <- dplyr::case_match(
        x,
        "0" ~ "0", "1" ~ "1", "2" ~ "2",
        "not grass" ~ "0", "nongrass" ~ "0", "non-grass" ~ "0",
        "tame" ~ "1", "tame grass" ~ "1",
        "native" ~ "2", "native grass" ~ "2",
        .default = x
      )
      ifelse(x %in% c("0","1","2"), as.integer(x), NA_integer_)
    },
    Yer_Cll  = suppressWarnings(as.integer(as.numeric(Yer_Cll)))
  )

# Normalize LC_Cl_N (used only for SHRUB exclusion if present)
has_lccln <- "LC_Cl_N" %in% names(points)
if (has_lccln) {
  points <- points %>%
    mutate(
      LC_Cl_N = {
        z <- as.character(LC_Cl_N)
        z <- stringr::str_trim(toupper(z))
        z[z %in% c("", "NA", "NAN", "NULL", "MISSING")] <- NA_character_
        z
      }
    )
} else {
  points <- points %>% mutate(LC_Cl_N = NA_character_)
}

# Normalize Source (if present) for training exclusion
has_source <- "Source" %in% names(points)
if (has_source) {
  points <- points %>%
    mutate(
      Source_norm = {
        s <- as.character(Source)
        s <- stringr::str_trim(tolower(s))
        s[s %in% c("", "na", "nan", "null", "missing")] <- NA_character_
        s
      }
    )
} else {
  points <- points %>% mutate(Source_norm = NA_character_)
}

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
  select(all_of(required_cols), LC_Cl_N, Source_norm, all_of(present_rasters)) %>%
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
compute_matched_year_class3 <- function(yrs, yer, lc3, notgrass_rule = c("exact","lte")) {
  notgrass_rule <- match.arg(notgrass_rule)
  yrs <- as.integer(yrs)
  if (length(yrs) == 0L || is.na(yer) || is.na(lc3)) return(NA_integer_)
  if (lc3 == 2L) {                # native
    ok <- yrs[yrs <= yer]
    if (length(ok) == 0L) return(NA_integer_)
    return(max(ok))
  } else if (lc3 == 1L) {         # tame
    if (yer %in% yrs) return(yer) else return(NA_integer_)
  } else if (lc3 == 0L) {         # not grass
    if (notgrass_rule == "exact") {
      if (yer %in% yrs) return(yer) else return(NA_integer_)
    } else { # "lte"
      ok <- yrs[yrs <= yer]
      if (length(ok) == 0L) return(NA_integer_)
      return(max(ok))
    }
  } else {
    return(NA_integer_)
  }
}

# Apply full filtering:
#   - require LC_class3 ∈ {0,1,2}, Yer_Cll present, pred present
#   - filter out data points where LC_class3 is 0 (not grass) - only 2 here
#   - optional exclusion: LC_Cl_N == "SHRUB"
#   - require allowed_years matching via matched_year
#   - exclude training sources that match atomic raster or grouped alias
eval_data <- long %>%
  filter(!is.na(LC_class3), LC_class3 %in% c(1L,2L), !is.na(Yer_Cll)) %>%  # <— keep only tame/native
  { if ("LC_Cl_N" %in% names(.)) filter(., is.na(LC_Cl_N) | LC_Cl_N != "SHRUB") else . } %>%
  filter(!is.na(pred)) %>%
  rowwise() %>%
  mutate(matched_year = compute_matched_year_class3(allowed_years, Yer_Cll, LC_class3)) %>%
  ungroup() %>%
  filter(!is.na(matched_year)) %>%
  filter(is.na(Source_norm) | !(Source_norm %in% c(raster, raster_group)))


# ----------------------------
# 8) Multiclass metrics helpers
# ----------------------------
safe_div <- function(n, d) ifelse(d == 0 | is.na(d), NA_real_, n / d)

compute_confusion_3 <- function(df, truth_col = "LC_class3", pred_col = "pred") {
  lvls <- c(0L,1L,2L)
  truth <- factor(df[[truth_col]], levels = lvls)
  pred  <- factor(df[[pred_col]],  levels = lvls)
  tab <- as.matrix(table(truth, pred))
  # Ensure full 3x3 even if missing levels
  if (!all(dimnames(tab)$truth == as.character(lvls))) {
    dimnames(tab)$truth <- as.character(lvls)
  }
  if (!all(dimnames(tab)$pred == as.character(lvls))) {
    dimnames(tab)$pred <- as.character(lvls)
  }
  n <- sum(tab)
  counts_long <- as_tibble(as.data.frame(as.table(tab))) %>%
    transmute(
      Truth = as.integer(as.character(truth)),
      Pred  = as.integer(as.character(pred)),
      N     = as.integer(Freq)
    )
  list(tab = tab, n = n, counts_long = counts_long)
}

compute_metrics_3 <- function(tab) {
  n <- sum(tab)
  diagv <- diag(tab)
  rows <- rowSums(tab)  # support by truth
  cols <- colSums(tab)  # predicted totals
  
  overall <- safe_div(sum(diagv), n)
  
  prec <- mapply(function(tp, colsum) safe_div(tp, colsum), diagv, cols)
  rec  <- mapply(function(tp, rowsum) safe_div(tp, rowsum), diagv, rows)
  f1   <- ifelse(is.na(prec) | is.na(rec) | (prec + rec) == 0, NA_real_, 2 * prec * rec / (prec + rec))
  
  macro_prec <- mean(prec, na.rm = TRUE)
  macro_rec  <- mean(rec, na.rm  = TRUE)
  macro_f1   <- mean(f1, na.rm   = TRUE)
  bal_acc    <- macro_rec
  
  weights <- safe_div(rows, sum(rows))
  weighted_f1 <- if (all(is.na(weights))) NA_real_ else sum((weights %||% 0) * (f1 %||% 0), na.rm = TRUE)
  
  # Cohen's kappa
  pe <- safe_div(sum(rows * cols), n * n)
  kappa <- ifelse(is.na(pe) | (1 - pe) == 0, NA_real_, (overall - pe) / (1 - pe))
  
  tibble(
    N_evaluated       = as.integer(n),
    Overall_Accuracy  = overall,
    Kappa             = kappa,
    Macro_Precision   = macro_prec,
    Macro_Recall      = macro_rec,
    Macro_F1          = macro_f1,
    Weighted_F1       = weighted_f1,
    Balanced_Accuracy = bal_acc,
    # Per-class metrics (suffix _0, _1, _2)
    Precision_0 = prec[1], Recall_0 = rec[1], F1_0 = f1[1], Support_0 = rows[1],
    Precision_1 = prec[2], Recall_1 = rec[2], F1_1 = f1[2], Support_1 = rows[2],
    Precision_2 = prec[3], Recall_2 = rec[3], F1_2 = f1[3], Support_2 = rows[3]
  )
}

summarize_block_multiclass <- function(df, key_cols) {
  # Group by key_cols; compute a 3x3 confusion + metrics per group
  df %>%
    group_by(across(all_of(key_cols))) %>%
    group_modify(~ {
      c3 <- compute_confusion_3(.x)
      met <- compute_metrics_3(c3$tab)
      # Return counts long in separate table; here only metrics + N
      met
    }) %>%
    ungroup()
}

confusion_block_multiclass <- function(df, key_cols) {
  # Return long-form confusion counts per group
  df %>%
    group_by(across(all_of(key_cols))) %>%
    group_modify(~ {
      compute_confusion_3(.x)$counts_long
    }) %>%
    ungroup()
}

# ----------------------------
# 9) Per-raster summaries (ATOMIC and GROUPED)
# ----------------------------

# ATOMIC: by (raster, matched_year)
summary_by_raster_atomic <- summarize_block_multiclass(
  eval_data,
  key_cols = c("raster", "matched_year")
) %>%
  rename(raster_year = matched_year) %>%
  arrange(desc(Overall_Accuracy))

# GROUPED: collapse pli_* -> "pli"; combine rows, recompute confusion+metrics
eval_grouped <- eval_data %>% mutate(raster = raster_group)

summary_by_raster_grouped <- summarize_block_multiclass(
  eval_grouped,
  key_cols = c("raster")
) %>%
  mutate(raster_year = NA_integer_) %>%
  relocate(raster, raster_year) %>%
  arrange(desc(Overall_Accuracy))

# ----------------------------
# 10) Confusion counts tables (3×3, long form)
# ----------------------------
confusion_counts_atomic <- confusion_block_multiclass(
  eval_data,
  key_cols = c("raster", "matched_year")
) %>%
  rename(raster_year = matched_year) %>%
  arrange(raster, raster_year, Truth, Pred)

confusion_counts_grouped <- confusion_block_multiclass(
  eval_grouped,
  key_cols = c("raster")
) %>%
  mutate(raster_year = NA_integer_) %>%
  relocate(raster, raster_year) %>%
  arrange(raster, Truth, Pred)

# ----------------------------
# 11) Per ground-truth year breakdowns (ATOMIC & GROUPED)
# ----------------------------
per_year_summary_atomic <- summarize_block_multiclass(
  eval_data,
  key_cols = c("raster", "matched_year", "Yer_Cll")
) %>%
  rename(raster_year = matched_year) %>%
  arrange(raster, Yer_Cll)

per_year_summary_grouped <- summarize_block_multiclass(
  eval_grouped,
  key_cols = c("raster", "Yer_Cll")
) %>%
  mutate(raster_year = NA_integer_) %>%
  relocate(raster, raster_year) %>%
  arrange(raster, Yer_Cll)

# ----------------------------
# 12) Save outputs (CSV)
# ----------------------------
readr::write_csv(summary_by_raster_grouped, file.path(out_dir, "raster_accuracy_summary_multiclass.csv"))
readr::write_csv(confusion_counts_grouped,  file.path(out_dir, "raster_confusion_counts_multiclass.csv"))
readr::write_csv(per_year_summary_grouped,  file.path(out_dir, "raster_per_year_accuracy_multiclass.csv"))

readr::write_csv(summary_by_raster_atomic, file.path(out_dir, "raster_accuracy_summary_multiclass_atomic.csv"))
readr::write_csv(confusion_counts_atomic,  file.path(out_dir, "raster_confusion_counts_multiclass_atomic.csv"))
readr::write_csv(per_year_summary_atomic,  file.path(out_dir, "raster_per_year_accuracy_multiclass_atomic.csv"))

# step 6 (multiclass): assess accuracy of datasets mapping grasslands
# against ground-truthed point data (2015 to present), with 3 classes:
#   0 = not grass, 1 = tame grass, 2 = native grass
#   (filter out groundtruthed points that are classified as 0)
#
# Key changes from binary version:
#   - Uses LC_class3 (0/1/2) instead of LC_grass (0/1)
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
points_excel_path <- "R:/CanadianGrasslands/WG_Classification/processed/groundtruthed_data_points/groundtruth_points_with_grass_raster_values_natvtam_pli_sep_ge2015_clean.xlsx"
metadata_csv_path <- "R:/CanadianGrasslands/WG_Classification/github/tomorrows-prairies/processeddatasets/Grass/metadata_full_extent_nattam_pli_sep.csv"

out_dir <- "R:/CanadianGrasslands/WG_Classification/results/accuracy_assessment/native_v_tame"
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
                 "gpw","plvi","hfi")

# Grouping helper: collapse pli_* -> "pli"; all others remain themselves
raster_to_group <- function(x) ifelse(x %in% pli_variants, "pli", x)

# ----------------------------
# 4) Load points and normalize basic fields
# ----------------------------
points <- readxl::read_excel(points_excel_path)

required_cols <- c("LC_class3", "Yer_Cll")
missing_req <- setdiff(required_cols, names(points))
if (length(missing_req) > 0) {
  stop("Missing required columns in points file: ", paste(missing_req, collapse = ", "))
}

# Normalize LC_class3 to integer 0/1/2, Yer_Cll to year int
points <- points %>%
  mutate(
    LC_class3 = {
      x <- as.character(LC_class3)
      x <- str_trim(tolower(x))
      x[x %in% c("", "na", "nan", "null", "missing")] <- NA_character_
      # Accept words if present
      x <- dplyr::case_match(
        x,
        "0" ~ "0", "1" ~ "1", "2" ~ "2",
        "not grass" ~ "0", "nongrass" ~ "0", "non-grass" ~ "0",
        "tame" ~ "1", "tame grass" ~ "1",
        "native" ~ "2", "native grass" ~ "2",
        .default = x
      )
      ifelse(x %in% c("0","1","2"), as.integer(x), NA_integer_)
    },
    Yer_Cll  = suppressWarnings(as.integer(as.numeric(Yer_Cll)))
  )

# Optional: Normalize LC_Cl_N (used only for SHRUB exclusion if present)
has_lccln <- "LC_Cl_N" %in% names(points)
if (has_lccln) {
  points <- points %>%
    mutate(
      LC_Cl_N = {
        z <- as.character(LC_Cl_N)
        z <- stringr::str_trim(toupper(z))
        z[z %in% c("", "NA", "NAN", "NULL", "MISSING")] <- NA_character_
        z
      }
    )
} else {
  points <- points %>% mutate(LC_Cl_N = NA_character_)
}

# Normalize Source (if present) for training exclusion
has_source <- "Source" %in% names(points)
if (has_source) {
  points <- points %>%
    mutate(
      Source_norm = {
        s <- as.character(Source)
        s <- stringr::str_trim(tolower(s))
        s[s %in% c("", "na", "nan", "null", "missing")] <- NA_character_
        s
      }
    )
} else {
  points <- points %>% mutate(Source_norm = NA_character_)
}

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
  select(all_of(required_cols), LC_Cl_N, Source_norm, all_of(present_rasters)) %>%
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
compute_matched_year_class3 <- function(yrs, yer, lc3, notgrass_rule = c("exact","lte")) {
  notgrass_rule <- match.arg(notgrass_rule)
  yrs <- as.integer(yrs)
  if (length(yrs) == 0L || is.na(yer) || is.na(lc3)) return(NA_integer_)
  if (lc3 == 2L) {                # native
    ok <- yrs[yrs <= yer]
    if (length(ok) == 0L) return(NA_integer_)
    return(max(ok))
  } else if (lc3 == 1L) {         # tame
    if (yer %in% yrs) return(yer) else return(NA_integer_)
  } else if (lc3 == 0L) {         # not grass
    if (notgrass_rule == "exact") {
      if (yer %in% yrs) return(yer) else return(NA_integer_)
    } else { # "lte"
      ok <- yrs[yrs <= yer]
      if (length(ok) == 0L) return(NA_integer_)
      return(max(ok))
    }
  } else {
    return(NA_integer_)
  }
}

# Apply full filtering:
#   - require LC_class3 ∈ {0,1,2}, Yer_Cll present, pred present
#   - filter out data points where LC_class3 is 0 (not grass) - only 2 here
#   - optional exclusion: LC_Cl_N == "SHRUB"
#   - require allowed_years matching via matched_year
#   - exclude training sources that match atomic raster or grouped alias
eval_data <- long %>%
  filter(!is.na(LC_class3), LC_class3 %in% c(1L,2L), !is.na(Yer_Cll)) %>%  # <— keep only tame/native
  { if ("LC_Cl_N" %in% names(.)) filter(., is.na(LC_Cl_N) | LC_Cl_N != "SHRUB") else . } %>%
  filter(!is.na(pred)) %>%
  rowwise() %>%
  mutate(matched_year = compute_matched_year_class3(allowed_years, Yer_Cll, LC_class3)) %>%
  ungroup() %>%
  filter(!is.na(matched_year)) %>%
  filter(is.na(Source_norm) | !(Source_norm %in% c(raster, raster_group)))


# ----------------------------
# 8) Multiclass metrics helpers
# ----------------------------
safe_div <- function(n, d) ifelse(d == 0 | is.na(d), NA_real_, n / d)

compute_confusion_3 <- function(df, truth_col = "LC_class3", pred_col = "pred") {
  lvls <- c(0L,1L,2L)
  truth <- factor(df[[truth_col]], levels = lvls)
  pred  <- factor(df[[pred_col]],  levels = lvls)
  tab <- as.matrix(table(truth, pred))
  # Ensure full 3x3 even if missing levels
  if (!all(dimnames(tab)$truth == as.character(lvls))) {
    dimnames(tab)$truth <- as.character(lvls)
  }
  if (!all(dimnames(tab)$pred == as.character(lvls))) {
    dimnames(tab)$pred <- as.character(lvls)
  }
  n <- sum(tab)
  counts_long <- as_tibble(as.data.frame(as.table(tab))) %>%
    transmute(
      Truth = as.integer(as.character(truth)),
      Pred  = as.integer(as.character(pred)),
      N     = as.integer(Freq)
    )
  list(tab = tab, n = n, counts_long = counts_long)
}

compute_metrics_3 <- function(tab) {
  n <- sum(tab)
  diagv <- diag(tab)
  rows <- rowSums(tab)  # support by truth
  cols <- colSums(tab)  # predicted totals
  
  overall <- safe_div(sum(diagv), n)
  
  prec <- mapply(function(tp, colsum) safe_div(tp, colsum), diagv, cols)
  rec  <- mapply(function(tp, rowsum) safe_div(tp, rowsum), diagv, rows)
  f1   <- ifelse(is.na(prec) | is.na(rec) | (prec + rec) == 0, NA_real_, 2 * prec * rec / (prec + rec))
  
  macro_prec <- mean(prec, na.rm = TRUE)
  macro_rec  <- mean(rec, na.rm  = TRUE)
  macro_f1   <- mean(f1, na.rm   = TRUE)
  bal_acc    <- macro_rec
  
  weights <- safe_div(rows, sum(rows))
  weighted_f1 <- if (all(is.na(weights))) NA_real_ else sum((weights %||% 0) * (f1 %||% 0), na.rm = TRUE)
  
  # Cohen's kappa
  pe <- safe_div(sum(rows * cols), n * n)
  kappa <- ifelse(is.na(pe) | (1 - pe) == 0, NA_real_, (overall - pe) / (1 - pe))
  
  tibble(
    N_evaluated       = as.integer(n),
    Overall_Accuracy  = overall,
    Kappa             = kappa,
    Macro_Precision   = macro_prec,
    Macro_Recall      = macro_rec,
    Macro_F1          = macro_f1,
    Weighted_F1       = weighted_f1,
    Balanced_Accuracy = bal_acc,
    # Per-class metrics (suffix _0, _1, _2)
    Precision_0 = prec[1], Recall_0 = rec[1], F1_0 = f1[1], Support_0 = rows[1],
    Precision_1 = prec[2], Recall_1 = rec[2], F1_1 = f1[2], Support_1 = rows[2],
    Precision_2 = prec[3], Recall_2 = rec[3], F1_2 = f1[3], Support_2 = rows[3]
  )
}

summarize_block_multiclass <- function(df, key_cols) {
  # Group by key_cols; compute a 3x3 confusion + metrics per group
  df %>%
    group_by(across(all_of(key_cols))) %>%
    group_modify(~ {
      c3 <- compute_confusion_3(.x)
      met <- compute_metrics_3(c3$tab)
      # Return counts long in separate table; here only metrics + N
      met
    }) %>%
    ungroup()
}

confusion_block_multiclass <- function(df, key_cols) {
  # Return long-form confusion counts per group
  df %>%
    group_by(across(all_of(key_cols))) %>%
    group_modify(~ {
      compute_confusion_3(.x)$counts_long
    }) %>%
    ungroup()
}

# ----------------------------
# 9) Per-raster summaries (ATOMIC and GROUPED)
# ----------------------------

# ATOMIC: by (raster, matched_year)
summary_by_raster_atomic <- summarize_block_multiclass(
  eval_data,
  key_cols = c("raster", "matched_year")
) %>%
  rename(raster_year = matched_year) %>%
  arrange(desc(Overall_Accuracy))

# GROUPED: collapse pli_* -> "pli"; combine rows, recompute confusion+metrics
eval_grouped <- eval_data %>% mutate(raster = raster_group)

summary_by_raster_grouped <- summarize_block_multiclass(
  eval_grouped,
  key_cols = c("raster")
) %>%
  mutate(raster_year = NA_integer_) %>%
  relocate(raster, raster_year) %>%
  arrange(desc(Overall_Accuracy))

# ----------------------------
# 10) Confusion counts tables (3×3, long form)
# ----------------------------
confusion_counts_atomic <- confusion_block_multiclass(
  eval_data,
  key_cols = c("raster", "matched_year")
) %>%
  rename(raster_year = matched_year) %>%
  arrange(raster, raster_year, Truth, Pred)

confusion_counts_grouped <- confusion_block_multiclass(
  eval_grouped,
  key_cols = c("raster")
) %>%
  mutate(raster_year = NA_integer_) %>%
  relocate(raster, raster_year) %>%
  arrange(raster, Truth, Pred)

# ----------------------------
# 11) Per ground-truth year breakdowns (ATOMIC & GROUPED)
# ----------------------------
per_year_summary_atomic <- summarize_block_multiclass(
  eval_data,
  key_cols = c("raster", "matched_year", "Yer_Cll")
) %>%
  rename(raster_year = matched_year) %>%
  arrange(raster, Yer_Cll)

per_year_summary_grouped <- summarize_block_multiclass(
  eval_grouped,
  key_cols = c("raster", "Yer_Cll")
) %>%
  mutate(raster_year = NA_integer_) %>%
  relocate(raster, raster_year) %>%
  arrange(raster, Yer_Cll)

# ----------------------------
# 12) Save outputs (CSV)
# ----------------------------
readr::write_csv(summary_by_raster_grouped, file.path(out_dir, "raster_accuracy_summary_multiclass.csv"))
readr::write_csv(confusion_counts_grouped,  file.path(out_dir, "raster_confusion_counts_multiclass.csv"))
readr::write_csv(per_year_summary_grouped,  file.path(out_dir, "raster_per_year_accuracy_multiclass.csv"))

readr::write_csv(summary_by_raster_atomic, file.path(out_dir, "raster_accuracy_summary_multiclass_atomic.csv"))
readr::write_csv(confusion_counts_atomic,  file.path(out_dir, "raster_confusion_counts_multiclass_atomic.csv"))
readr::write_csv(per_year_summary_atomic,  file.path(out_dir, "raster_per_year_accuracy_multiclass_atomic.csv"))

# ----------------------------
# 13) Plots: Overall accuracy (GROUPED and ATOMIC)
# ----------------------------
p_grouped <- ggplot(summary_by_raster_grouped,
                    aes(x = reorder(raster, Overall_Accuracy), y = Overall_Accuracy)) +
  geom_col(fill = "#2d7f5e") +
  coord_flip() +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.02))) +
  labs(title = "Overall Accuracy by Raster (Grouped: pli combines variants) — Multiclass (0/1/2)",
       x = "Raster", y = "Overall Accuracy") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(out_dir, "raster_accuracy_plot_multiclass.png"),
       p_grouped, width = 7, height = 5, dpi = 200)

p_atomic <- ggplot(summary_by_raster_atomic,
                   aes(x = reorder(raster, Overall_Accuracy), y = Overall_Accuracy)) +
  geom_col(fill = "#2d7f5e") +
  coord_flip() +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.02))) +
  labs(title = "Overall Accuracy by Raster (Atomic: pli_ap, pli_cu, pli_mg, pli_mm) — Multiclass (0/1/2)",
       x = "Raster", y = "Overall Accuracy") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(out_dir, "raster_accuracy_plot_multiclass_atomic.png"),
       p_atomic, width = 7, height = 5, dpi = 200)

# ----------------------------
# 14) Console previews
# ----------------------------
cat("\n=== Summary by Raster (GROUPED) — Multiclass ===\n")
print(summary_by_raster_grouped, n = Inf)

cat("\n=== Summary by Raster (ATOMIC) — Multiclass ===\n")
print(summary_by_raster_atomic, n = Inf)

cat("\n=== Confusion Counts (GROUPED) — Multiclass (first 60 rows) ===\n")
print(confusion_counts_grouped %>% slice_head(n = 60))

cat("\n=== Confusion Counts (ATOMIC) — Multiclass (first 60 rows) ===\n")
print(confusion_counts_atomic %>% slice_head(n = 60))

# ----------------------------
# 15) Tame vs Native (1 vs 2) — 2×2 confusion & per-class metrics
# ----------------------------

# Work only with ground-truth in {1,2} — already enforced in eval_data,
# but keep this filter here to be explicit.
eval_grass <- eval_data %>% filter(LC_class3 %in% c(1L, 2L))

safe_div <- function(n, d) ifelse(d == 0 | is.na(d), NA_real_, n / d)

# Build a true 2×2 table between classes {1,2}; predictions outside {1,2}
# (e.g., pred==0) are counted separately (Pred_outside_12) and excluded from the 2×2.
compute_confusion_2_12 <- function(df) {
  truth <- factor(df$LC_class3, levels = c(1L, 2L))
  pred_restricted <- factor(dplyr::if_else(df$pred %in% c(1L, 2L), df$pred, NA_integer_),
                            levels = c(1L, 2L))
  tab <- as.matrix(table(truth, pred_restricted, useNA = "no"))
  # Ensure 2×2 shape with named dims
  all_rows <- c("1","2"); all_cols <- c("1","2")
  if (nrow(tab) == 0 || ncol(tab) == 0) {
    tab <- matrix(0L, nrow = 2, ncol = 2, dimnames = list(all_rows, all_cols))
  } else {
    # align order
    if (!identical(rownames(tab), all_rows)) tab <- tab[match(all_rows, rownames(tab)), , drop = FALSE]
    if (!identical(colnames(tab), all_cols)) tab <- tab[, match(all_cols, colnames(tab)), drop = FALSE]
    tab[is.na(tab)] <- 0L
  }
  pred_other <- sum(!(df$pred %in% c(1L, 2L)))
  list(tab = tab, pred_other = pred_other)
}

summarize_2x2_block <- function(df, key_cols) {
  df %>%
    group_by(across(all_of(key_cols))) %>%
    group_modify(~{
      c2  <- compute_confusion_2_12(.x)
      tab <- c2$tab
      
      # 2x2 cells (truth rows: 1,2; pred cols: 1,2)
      a <- as.numeric(tab["1","1"])  # true 1 predicted 1
      b <- as.numeric(tab["1","2"])  # true 1 predicted 2
      c <- as.numeric(tab["2","1"])  # true 2 predicted 1
      d <- as.numeric(tab["2","2"])  # true 2 predicted 2
      
      n    <- a + b + c + d
      row1 <- a + b
      row2 <- c + d
      col1 <- a + c
      col2 <- b + d
      
      # Core accuracy
      oa  <- safe_div(a + d, n)
      r1  <- safe_div(a, row1)   # Recall/Sensitivity for class 1 (tame)
      r2  <- safe_div(d, row2)   # Recall/Sensitivity for class 2 (native)
      
      # New per-class metrics
      p1  <- safe_div(a, col1)   # Precision/PPV for class 1
      p2  <- safe_div(d, col2)   # Precision/PPV for class 2
      
      sp1 <- safe_div(d, d + c)  # Specificity for class 1 (== recall of class 2)
      sp2 <- safe_div(a, a + b)  # Specificity for class 2 (== recall of class 1)
      
      fpr1 <- ifelse(is.na(sp1), NA_real_, 1 - sp1)
      fpr2 <- ifelse(is.na(sp2), NA_real_, 1 - sp2)
      fnr1 <- ifelse(is.na(r1),  NA_real_, 1 - r1)
      fnr2 <- ifelse(is.na(r2),  NA_real_, 1 - r2)
      
      f1_1 <- ifelse(is.na(p1) | is.na(r1) | (p1 + r1) == 0, NA_real_, 2 * p1 * r1 / (p1 + r1))
      f1_2 <- ifelse(is.na(p2) | is.na(r2) | (p2 + r2) == 0, NA_real_, 2 * p2 * r2 / (p2 + r2))
      
      bal_acc <- ifelse(all(is.na(c(r1, r2))), NA_real_, mean(c(r1, r2), na.rm = TRUE))
      
      # MCC (symmetric for either class as "positive")
      mcc_num <- (a * d) - (b * c)
      mcc_den <- sqrt((a + c) * (a + b) * (d + c) * (d + b))
      mcc     <- safe_div(mcc_num, mcc_den)
      
      # Kappa (binary)
      pe     <- safe_div(row1 * col1 + row2 * col2, n * n)
      kappa2 <- ifelse(is.na(pe) | (1 - pe) == 0, NA_real_, (oa - pe) / (1 - pe))
      
      # Prevalence (truth proportions)
      prev1 <- safe_div(row1, n)
      prev2 <- safe_div(row2, n)
      
      # Negative Predictive Value (per class perspective)
      npv1 <- safe_div(d, d + b)  # P(true not-1 | predicted not-1)
      npv2 <- safe_div(a, a + c)  # P(true not-2 | predicted not-2)
      
      # Likelihood Ratios
      lrp1 <- safe_div(r1, fpr1)   # LR+ for class 1
      lrm1 <- safe_div(fnr1, sp1)  # LR- for class 1
      lrp2 <- safe_div(r2, fpr2)   # LR+ for class 2
      lrm2 <- safe_div(fnr2, sp2)  # LR- for class 2
      
      # Macro/weighted F1 within 2x2 (optional but handy)
      macro_f1_2x2    <- ifelse(all(is.na(c(f1_1, f1_2))), NA_real_, mean(c(f1_1, f1_2), na.rm = TRUE))
      weighted_f1_2x2 <- ifelse(is.na(prev1) & is.na(prev2), NA_real_, sum(c(prev1, prev2) * c(f1_1, f1_2), na.rm = TRUE))
      
      tibble(
        # Counts / coverage
        N_in_2x2        = as.integer(n),
        Pred_outside_12 = as.integer(c2$pred_other),
        
        # Overall
        OA_2x2_1v2            = oa,
        Balanced_Accuracy_2x2 = bal_acc,
        MCC_1v2               = mcc,
        Kappa_2x2             = kappa2,
        Macro_F1_2x2          = macro_f1_2x2,
        Weighted_F1_2x2       = weighted_f1_2x2,
        Prevalence_1          = prev1,
        Prevalence_2          = prev2,
        
        # Class 1 (tame)
        Recall_1_2x2      = r1,
        Precision_1_2x2   = p1,
        Specificity_1_2x2 = sp1,
        FPR_1_2x2         = fpr1,
        FNR_1_2x2         = fnr1,
        F1_1_2x2          = f1_1,
        NPV_1_2x2         = npv1,
        LRplus_1          = lrp1,
        LRminus_1         = lrm1,
        
        # Class 2 (native)
        Recall_2_2x2      = r2,
        Precision_2_2x2   = p2,
        Specificity_2_2x2 = sp2,
        FPR_2_2x2         = fpr2,
        FNR_2_2x2         = fnr2,
        F1_2_2x2          = f1_2,
        NPV_2_2x2         = npv2,
        LRplus_2          = lrp2,
        LRminus_2         = lrm2,
        
        # Directional confusions (for readability)
        Confuse_1_as_2 = as.integer(b),
        Confuse_2_as_1 = as.integer(c)
      )
    }) %>%
    ungroup()
}

# ATOMIC (by raster × matched_year)
confusion_2x2_atomic <- summarize_2x2_block(
  eval_grass, c("raster", "matched_year")
) %>%
  rename(raster_year = matched_year) %>%
  arrange(desc(OA_2x2_1v2))

# GROUPED (pli_* collapsed)
confusion_2x2_grouped <- summarize_2x2_block(
  eval_grass %>% mutate(raster = raster_group), c("raster")
) %>%
  mutate(raster_year = NA_integer_) %>%
  relocate(raster, raster_year) %>%
  arrange(desc(OA_2x2_1v2))

# Save
readr::write_csv(confusion_2x2_grouped, file.path(out_dir, "raster_confusion_2x2_tame_vs_native.csv"))
readr::write_csv(confusion_2x2_atomic,  file.path(out_dir, "raster_confusion_2x2_tame_vs_native_atomic.csv"))

cat("\n=== 2x2 Tame vs Native (GROUPED) ===\n")
print(confusion_2x2_grouped, n = Inf)
cat("\n=== 2x2 Tame vs Native (ATOMIC) ===\n")
print(confusion_2x2_atomic, n = Inf)

# ----------------------------
# 16) Plots: Per-class Recall (tame=1, native=2)
# ----------------------------
make_recall_bar <- function(df, title, outfile) {
  plot_df <- df %>%
    select(raster, Recall_1, Recall_2) %>%
    tidyr::pivot_longer(
      cols = starts_with("Recall_"),
      names_to = "Class",
      values_to = "Recall"
    ) %>%
    mutate(
      Class = dplyr::recode(Class,
                            "Recall_1" = "Tame (1)",
                            "Recall_2" = "Native (2)")
    )
  p <- ggplot(plot_df,
              aes(x = reorder(raster, Recall, FUN = median, na.rm = TRUE),
                  y = Recall, fill = Class)) +
    geom_col(position = position_dodge(width = 0.72), width = 0.7) +
    coord_flip() +
    scale_fill_manual(values = c("Tame (1)" = "#fdb863", "Native (2)" = "#5ab4ac")) +
    scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.02))) +
    labs(title = title, x = "Raster", y = "Recall (Sensitivity)") +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(), legend.position = "right")
  ggsave(file.path(out_dir, outfile), p, width = 8, height = 5, dpi = 200)
}

# Grouped and Atomic per-class recall plots
make_recall_bar(summary_by_raster_grouped,
                "Recall by Class (Grouped: pli combines variants)",
                "raster_recall_by_class_grouped.png")

make_recall_bar(summary_by_raster_atomic,
                "Recall by Class (Atomic: pli_ap, pli_cu, pli_mg, pli_mm)",
                "raster_recall_by_class_atomic.png")

# ----------------------------
# 17) Plots: Compact 2x2 metrics (tame vs native)
# ----------------------------

# This function expects the extended 2x2 columns from your updated summarize_2x2_block().
# It builds a compact, two-facet plot with key metrics:
#   - For each class: Precision, Recall, F1
#   - Overall (non-faceted): OA_2x2_1v2 (as a subtle reference band)
#
# It supports both GROUPED (raster-level) and ATOMIC (raster + year) tables.
make_2x2_compact_plot <- function(df, title, outfile,
                                  order_by = c("F1", "OA"),
                                  top_n = Inf) {
  
  order_by <- match.arg(order_by)
  
  # ---- 1) Robust raster_label creation (handles grouped or atomic) ----
  df <- df %>%
    mutate(
      raster_label =
        if ("raster_year" %in% names(.)) {
          # vectorized row-wise choice
          dplyr::if_else(
            is.na(raster_year),
            raster,
            paste0(raster, " (", raster_year, ")")
          )
        } else {
          raster
        }
    )
  
  # ---- 2) Required columns ----
  needed_cols <- c(
    "raster_label", "OA_2x2_1v2",
    "Precision_1_2x2","Recall_1_2x2","F1_1_2x2",
    "Precision_2_2x2","Recall_2_2x2","F1_2_2x2"
  )
  
  missing_needed <- setdiff(needed_cols, names(df))
  if (length(missing_needed) > 0) {
    stop("Missing required columns for plotting 2x2 compact metrics: ",
         paste(missing_needed, collapse = ", "))
  }
  
  # ---- 3) Long-format for metrics ----
  plot_df <- df %>%
    select(all_of(needed_cols)) %>%
    tidyr::pivot_longer(
      cols = -c(raster_label, OA_2x2_1v2),
      names_to = "Metric",
      values_to = "Score"
    ) %>%
    mutate(
      Class = dplyr::case_when(
        stringr::str_detect(Metric, "_1_") ~ "Tame (1)",
        stringr::str_detect(Metric, "_2_") ~ "Native (2)",
        TRUE ~ "Overall"
      ),
      Metric_simple = dplyr::recode(Metric,
                                    "Precision_1_2x2" = "Precision",
                                    "Recall_1_2x2"    = "Recall",
                                    "F1_1_2x2"        = "F1",
                                    "Precision_2_2x2" = "Precision",
                                    "Recall_2_2x2"    = "Recall",
                                    "F1_2_2x2"        = "F1"
      )
    ) %>%
    filter(Class %in% c("Tame (1)", "Native (2)"))
  
  # ---- 4) Compute ordering metric (F1 or OA), aggregated per raster ----
  order_df <- plot_df %>%
    group_by(raster_label, Class) %>%
    summarise(F1_avg = mean(Score[Metric_simple == "F1"], na.rm = TRUE),
              .groups = "drop") %>%
    group_by(raster_label) %>%
    summarise(F1_overall = mean(F1_avg, na.rm = TRUE),
              .groups = "drop") %>%
    left_join(df %>% select(raster_label, OA_2x2_1v2) %>% distinct(),
              by = "raster_label")
  
  # ---- FIX #1: Use if() outside mutate (not if_else) ----
  order_df <- order_df %>%
    mutate(order_score = if (order_by == "F1") F1_overall else OA_2x2_1v2) %>%
    arrange(desc(order_score)) %>%
    slice_head(n = top_n)
  
  # ---- 5) Filter plot_df to top_n and apply ordering ----
  plot_df <- plot_df %>%
    semi_join(order_df, by = "raster_label") %>%
    mutate(raster_label = factor(raster_label,
                                 levels = rev(order_df$raster_label)))
  
  overall_band <- df %>%
    semi_join(order_df, by = "raster_label") %>%
    select(raster_label, OA_2x2_1v2) %>%
    distinct() %>%
    mutate(raster_label = factor(raster_label,
                                 levels = levels(plot_df$raster_label)))
  
  # ---- 6) Plot ----
  metric_cols <- c(
    "Precision" = "#80b1d3",
    "Recall"    = "#8dd3c7",
    "F1"        = "#fdb462"
  )
  
  p <- ggplot() +
    geom_rect(
      data = overall_band,
      aes(ymin = 0, ymax = OA_2x2_1v2,
          xmin = as.numeric(raster_label) - 0.45,
          xmax = as.numeric(raster_label) + 0.45),
      fill = "grey92", color = NA
    ) +
    geom_col(
      data = plot_df,
      aes(x = raster_label, y = Score, fill = Metric_simple),
      position = position_dodge(width = 0.72),
      width = 0.7
    ) +
    coord_flip() +
    facet_wrap(~ Class, ncol = 2, scales = "free_y") +
    scale_fill_manual(values = metric_cols) +
    
    scale_y_continuous(
      limits = c(0, 1),
      expand = c(0, 0.02)
    ) +
    labs(
      title    = title,
      subtitle = "Bars: Precision, Recall, F1 per class; light band = Overall Accuracy (OA)",
      x        = "Raster",
      y        = "Score (0–1)",
      fill     = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    )
  
  ggsave(filename = file.path(out_dir, outfile),
         plot = p, width = 10, height = 6, dpi = 200)
}

# ---- Create compact plots for both GROUPED and ATOMIC ----
make_2x2_compact_plot(
  confusion_2x2_grouped,
  title   = "Tame vs Native - Compact Metrics (Grouped: pli combines variants)",
  outfile = "tame_native_compact_metrics_grouped.png",
  order_by = "F1",
  top_n = Inf  # or set e.g., 30 to show only the top 30 by F1
)

make_2x2_compact_plot(
  confusion_2x2_atomic,
  title   = "Tame vs Native - Compact Metrics (Atomic: includes pli_ap, pli_cu, pli_mg, pli_mm)",
  outfile = "tame_native_compact_metrics_atomic.png",
  order_by = "F1",
  top_n = Inf
)

# ----------------------------
# 18) Plots: MCC (GROUPED and ATOMIC)
# ----------------------------

# MCC accounts for all confusion-matrix entries and remains informative under class imbalance; 
# 0 indicates random performance, negative values indicate systematic disagreement

# Raster performance for the ecologically relevant tame–native distinction was summarized 
# using the Matthews Correlation Coefficient (MCC), which incorporates all four cells of the 
# 2×2 confusion matrix and is robust to class imbalance.

# When reading the plot:
#  MCC ≈ 1 → excellent discrimination
#  MCC ≈ 0 → no skill (random / prevalence‑driven)
#  MCC < 0 → systematic disagreement (very informative!)

# HFI dropped since it doesn't contain "native" class

# reusable helper function
make_mcc_2x2_bar <- function(df, title, outfile, top_n = Inf) {
  
  stopifnot("MCC_1v2" %in% names(df))
  
  df_plot <- df %>%
    filter(!is.na(MCC_1v2))   # <-- DROP hfi here
  
  # ---- raster label ----
  df_plot <- df_plot %>%
    mutate(
      raster_label = if ("raster_year" %in% names(.)) {
        dplyr::if_else(
          is.na(raster_year),
          raster,
          paste0(raster, " (", raster_year, ")")
        )
      } else {
        raster
      }
    )
  
  # ---- order ----
  df_plot <- df_plot %>%
    arrange(desc(MCC_1v2)) %>%
    slice_head(n = top_n) %>%
    mutate(raster_label = factor(raster_label, levels = rev(raster_label)))
  
  # ---- plot ----
  p <- ggplot(df_plot, aes(x = raster_label, y = MCC_1v2)) +
    geom_col(fill = "#6a5acd") +
    coord_flip() +
    scale_y_continuous(
      limits = c(-1, 1),
      breaks = seq(-1, 1, by = 0.25),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(
      title = title,
      subtitle = "2×2 Tame (1) vs Native (2); rasters with undefined MCC omitted",
      x = "Raster",
      y = "Matthews Correlation Coefficient (MCC)"
    ) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank())
  
  ggsave(
    file.path(out_dir, outfile),
    p,
    width = 8, height = 5, dpi = 200
  )
}

# grouped plot (pli combined)
make_mcc_2x2_bar(
  confusion_2x2_grouped,
  title   = "Tame vs Native (1 vs 2) — MCC by Raster (Grouped)",
  outfile = "raster_mcc_2x2_grouped.png"
)

# atomic plot (raster x year)
make_mcc_2x2_bar(
  confusion_2x2_atomic,
  title   = "Tame vs Native (1 vs 2) — MCC by Raster-Year (Atomic)",
  outfile = "raster_mcc_2x2_atomic.png",
  top_n   = 40   # optional, helps readability
)

print_mcc_2x2 <- function(df, label = "") {
  
  stopifnot("MCC_1v2" %in% names(df))
  
  out <- df %>%
    mutate(
      raster_label = if ("raster_year" %in% names(.)) {
        dplyr::if_else(
          is.na(raster_year),
          raster,
          paste0(raster, " (", raster_year, ")")
        )
      } else {
        raster
      }
    ) %>%
    select(raster_label, MCC_1v2, N_in_2x2, Pred_outside_12) %>%
    arrange(desc(MCC_1v2))
  
  cat("\n=== MCC (2×2) by Raster", label, "===\n")
  print(out, n = Inf)
  
  invisible(out)
}

print_mcc_2x2(confusion_2x2_grouped, "(GROUPED)")
print_mcc_2x2(confusion_2x2_atomic,  "(ATOMIC)")

# ----------------------------
# 19) QC: Ground-truth rows filtered out (LC_class3 == 0)
# ----------------------------
qc_notgrass_points <- points %>%
  filter(!is.na(LC_class3), LC_class3 == 0L)

# Save full rows for auditing (keep columns present in your input)
readr::write_csv(qc_notgrass_points,
                 file.path(out_dir, "qc_filtered_groundtruth_LC_class3_eq0_points.csv"))

# Counts by year and overall
qc_counts_by_year <- qc_notgrass_points %>%
  group_by(Yer_Cll) %>%
  summarise(N = n(), .groups = "drop") %>%
  arrange(Yer_Cll)

qc_counts_overall <- tibble(Total_filtered = nrow(qc_notgrass_points))

readr::write_csv(qc_counts_by_year,
                 file.path(out_dir, "qc_filtered_LC_class3_eq0_by_year.csv"))

readr::write_csv(qc_counts_overall,
                 file.path(out_dir, "qc_filtered_LC_class3_eq0_summary.csv"))

cat("\nQC filtered LC_class3==0 rows: ",
    nrow(qc_notgrass_points),
    " (see CSVs in ", out_dir, ")\n", sep = "")

####################################################################################
## Grass-only evaluation & Ecodistrict (optional, adapted for multiclass)
##
## Grass-only = LC_class3 in {1,2}; "hit" if pred in {1,2}.
####################################################################################

compute_grass_only_counts <- function(df) {
  df2 <- df %>% filter(LC_class3 %in% c(1L,2L))
  tp <- sum(df2$pred %in% c(1L,2L), na.rm = TRUE)
  fn <- sum(df2$pred %in% c(0L), na.rm = TRUE)
  n  <- tp + fn
  tibble(
    N_grass_evaluated = n,
    TP = tp,
    FN = fn,
    Recall_TPR = ifelse(n == 0, NA_real_, tp / n),
    OmissionRate_FNR = ifelse(n == 0, NA_real_, fn / n)
  )
}

# Grass-only summaries (ATOMIC)
summary_grass_only_atomic <- eval_data %>%
  filter(LC_class3 %in% c(1L,2L)) %>%
  group_by(raster, matched_year) %>%
  group_modify(~ compute_grass_only_counts(.x)) %>%
  ungroup() %>%
  rename(raster_year = matched_year) %>%
  arrange(desc(Recall_TPR))

per_year_grass_only_atomic <- eval_data %>%
  filter(LC_class3 %in% c(1L,2L)) %>%
  group_by(raster, matched_year, Yer_Cll) %>%
  group_modify(~ compute_grass_only_counts(.x)) %>%
  ungroup() %>%
  rename(raster_year = matched_year) %>%
  arrange(raster, Yer_Cll)

# Grass-only summaries (GROUPED)
summary_grass_only_grouped <- eval_grouped %>%
  filter(LC_class3 %in% c(1L,2L)) %>%
  group_by(raster) %>%
  group_modify(~ compute_grass_only_counts(.x)) %>%
  ungroup() %>%
  mutate(raster_year = NA_integer_) %>%
  arrange(desc(Recall_TPR))

per_year_grass_only_grouped <- eval_grouped %>%
  filter(LC_class3 %in% c(1L,2L)) %>%
  group_by(raster, Yer_Cll) %>%
  group_modify(~ compute_grass_only_counts(.x)) %>%
  ungroup() %>%
  mutate(raster_year = NA_integer_) %>%
  arrange(raster, Yer_Cll)

# Save grass-only outputs (optional)
readr::write_csv(summary_grass_only_grouped, file.path(out_dir, "raster_grassonly_summary.csv"))
readr::write_csv(per_year_grass_only_grouped, file.path(out_dir, "raster_grassonly_per_year.csv"))
readr::write_csv(summary_grass_only_atomic, file.path(out_dir, "raster_grassonly_summary_atomic.csv"))
readr::write_csv(per_year_grass_only_atomic, file.path(out_dir, "raster_grassonly_per_year_atomic.csv"))

# ----------------------------
# Optional: Maps of Grass-only Recall by Ecodistrict 
# ----------------------------

# Default CSVs for maps
metrics_csv_grouped <- file.path(out_dir, "raster_grassonly_per_year.csv") # can also map the non-year version
metrics_csv_atomic  <- file.path(out_dir, "raster_grassonly_summary_atomic.csv") # illustrative

# Ecodistrict boundary layer (must have a name column matching ECODISTRIC in points)
ecodistrict_path <- "R:/CanadianGrasslands/processeddatasets/boundaries/national-ecological-framework-for-canada/Ecodistricts/ecodistricts_clean.shp"
eco_name_col     <- "ECODISTRIC"   # in the boundary layer
map_crs <- 5070  # e.g., 5070 (NAD83 / Conus Albers)

# For grouped map: recompute metrics by raster × ECODISTRIC (grass-only)
points_raw <- readxl::read_excel(points_excel_path)
required_cols_eco <- c("LC_class3", "Yer_Cll", "ECODISTRIC")
missing_req_eco <- setdiff(required_cols_eco, names(points_raw))
if (length(missing_req_eco) > 0) {
  warning("Skipping ecodistrict maps: missing required columns: ", paste(missing_req_eco, collapse = ", "))
} else {
  present_rasters_eco <- intersect(raster_cols, names(points_raw))
  points2 <- points_raw %>%
    mutate(
      LC_class3 = {
        x <- as.character(LC_class3)
        x <- str_trim(tolower(x))
        x[x %in% c("", "na", "nan", "null", "missing")] <- NA_character_
        x <- dplyr::case_match(
          x,
          "0" ~ "0", "1" ~ "1", "2" ~ "2",
          "not grass" ~ "0", "nongrass" ~ "0", "non-grass" ~ "0",
          "tame" ~ "1", "tame grass" ~ "1",
          "native" ~ "2", "native grass" ~ "2",
          .default = x
        )
        ifelse(x %in% c("0","1","2"), as.integer(x), NA_integer_)
      },
      Yer_Cll = suppressWarnings(as.integer(as.numeric(Yer_Cll))),
      ECODISTRIC = {
        e <- as.character(ECODISTRIC)
        e <- str_trim(e)
        e[e %in% c("", "NA", "NaN", "null", "NULL", "missing")] <- NA_character_
        e
      }
    )
  has_source2 <- "Source" %in% names(points2)
  if (has_source2) {
    points2 <- points2 %>%
      mutate(
        Source_norm = {
          s <- as.character(Source)
          s <- stringr::str_trim(tolower(s))
          s[s %in% c("", "na", "nan", "null", "missing")] <- NA_character_
          s
        }
      )
  } else {
    points2 <- points2 %>% mutate(Source_norm = NA_character_)
  }
  
  points2_norm <- points2 %>%
    mutate(
      across(
        all_of(present_rasters_eco),
        ~ {
          x <- as.character(.)
          x <- str_trim(tolower(x))
          x[x %in% c("", "na", "nan", "null", "missing")] <- NA_character_
          ifelse(x %in% c("0","1","2"), as.integer(x), NA_integer_)
        }
      )
    )
  
  long2 <- points2_norm %>%
    select(all_of(required_cols_eco), Source_norm, all_of(present_rasters_eco)) %>%
    pivot_longer(
      cols = all_of(present_rasters_eco),
      names_to = "raster",
      values_to = "pred"
    ) %>%
    mutate(
      raster = tolower(raster),
      raster_group = raster_to_group(raster),
      allowed_years = map(raster, ~ raster_years_map[[.x]] %||% integer(0))
    )
  
  eval_data2 <- long2 %>%
    filter(
      !is.na(LC_class3), LC_class3 %in% c(1L,2L),  # <— drop 0s (not grass)
      !is.na(Yer_Cll),
      !is.na(pred)
    ) %>%
    
    rowwise() %>%
    mutate(matched_year = compute_matched_year_class3(allowed_years, Yer_Cll, LC_class3)) %>%
    ungroup() %>%
    filter(!is.na(matched_year)) %>%
    filter(is.na(Source_norm) | !(Source_norm %in% c(raster, raster_group))) %>%
    # restrict to ecodistricts with grass-only ground truth
    filter(!is.na(ECODISTRIC))
  
  # GROUPED grass-only metrics per ecodistrict
  by_raster_ecodistrict_grouped <- eval_data2 %>%
    filter(LC_class3 %in% c(1L,2L)) %>%
    mutate(raster = raster_group) %>%
    group_by(raster, ECODISTRIC) %>%
    summarise(
      TP = sum(pred %in% c(1L,2L), na.rm = TRUE),
      FN = sum(pred %in% c(0L), na.rm = TRUE),
      N_grass_evaluated = TP + FN,
      Recall_TPR = ifelse(N_grass_evaluated > 0, TP / N_grass_evaluated, NA_real_),
      OmissionRate_FNR = ifelse(N_grass_evaluated > 0, FN / N_grass_evaluated, NA_real_),
      .groups = "drop"
    ) %>%
    mutate(raster_year = NA_integer_) %>%
    filter(N_grass_evaluated > 0) %>%
    arrange(raster, ECODISTRIC)
  
  readr::write_csv(
    by_raster_ecodistrict_grouped,
    file.path(out_dir, "raster_grassonly_by_ecodistrict.csv")
  )
  
  # Read polygons
  eco_sf <- sf::st_read(ecodistrict_path, quiet = TRUE)
  
  # Ensure both sides use the same join type (character) and are trimmed
  eco_sf <- eco_sf %>%
    mutate(!!eco_name_col := as.character(.data[[eco_name_col]]) %>% stringr::str_trim())
  by_raster_ecodistrict_grouped <- by_raster_ecodistrict_grouped %>%
    mutate(ECODISTRIC = as.character(ECODISTRIC) %>% stringr::str_trim())
  
  # Reproject polygons if requested
  if (!is.na(map_crs)) eco_sf <- sf::st_transform(eco_sf, map_crs)
  
  # Color scale
  recall_scale <- scale_fill_gradientn(
    colours = c("#c6dbef", "#6baed6", "#2171b5", "#08306b", "#030f22"),
    limits = c(0, 1),
    oob = squish,
    na.value = "grey90",
    name = "Recall (grass)"
  )
  
  # Helper: build a map df for a single raster, restricted to districts with data
  make_map_data_filtered <- function(raster_code, metrics_df, eco_sf, eco_name_col) {
    # Metrics for this raster
    m <- metrics_df %>% filter(raster == raster_code, N_grass_evaluated > 0)
    
    # If no metrics for this raster, return empty result
    if (nrow(m) == 0) {
      return(list(bg_subset = eco_sf[0, ], map_df = eco_sf[0, ], bbox = NULL))
    }
    
    # Subset polygons to only districts that appear in metrics for this raster
    # (this ensures polygons without points are not shown)
    active_keys <- unique(m$ECODISTRIC)
    bg_subset <- eco_sf %>% filter(.data[[eco_name_col]] %in% active_keys)
    
    # Guard: remove empty geometries if any slipped in
    bg_subset <- bg_subset %>% filter(!sf::st_is_empty(geometry))
    
    # Join metrics onto the subset (left join on name col)
    join_by <- setNames("ECODISTRIC", eco_name_col)
    map_df <- bg_subset %>%
      left_join(m, by = join_by) %>%
      mutate(has_data = !is.na(N_grass_evaluated) & N_grass_evaluated > 0)
    
    # Compute a tight bbox for cropping
    if (nrow(map_df) > 0) {
      bbox <- sf::st_bbox(map_df)
    } else {
      bbox <- NULL
    }
    
    list(bg_subset = bg_subset, map_df = map_df, bbox = bbox)
  }
  
  # Loop over rasters and write cropped maps
  rasters_to_map_grouped <- sort(unique(by_raster_ecodistrict_grouped$raster))
  for (r in rasters_to_map_grouped) {
    res <- make_map_data_filtered(r, by_raster_ecodistrict_grouped, eco_sf, eco_name_col)
    
    # If nothing to plot (no active ecodistricts), skip
    if (is.null(res$bbox) || nrow(res$map_df) == 0) {
      message("Skipping map for '", r, "' — no ecodistricts with evaluated points.")
      next
    }
    
    # Build plot — background and filled polygons both restricted to active ecodistricts
    p <- ggplot() +
      geom_sf(data = res$bg_subset, fill = "grey95", color = "white", size = 0.2) +
      geom_sf(
        data = dplyr::filter(res$map_df, has_data),
        aes(fill = Recall_TPR),
        color = "white",
        size = 0.2
      ) +
      recall_scale +
      labs(
        title = paste0("Grass Recall by Ecodistrict — ", toupper(r), " (Grouped, multiclass)"),
        subtitle = "Only ecodistricts with evaluated grass points are shown.\nGrass = {tame, native}. Time rules: native ≤ Yer_Cll; tame = Yer_Cll",
        caption = NULL
      ) +
      theme_minimal(base_size = 12) +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position = "right")
    
    # Crop the visible extent to the active ecodistricts' bbox
    p <- p + coord_sf(
      xlim = c(res$bbox["xmin"], res$bbox["xmax"]),
      ylim = c(res$bbox["ymin"], res$bbox["ymax"]),
      expand = FALSE
    )
    
    ggsave(
      filename = file.path(out_dir, paste0("map_grass_recall_by_ecodistrict_", r, ".png")),
      plot = p, width = 9, height = 7, dpi = 200
    )
  }
  
  cat("\nSaved grouped ecodistrict maps to: ", out_dir, "\n")
  
}  

# ----------------------------
# FOR CONSENSUS MAPPING: Combine across years per raster × ecodistrict (micro-average)
# ----------------------------
# ============================================================
#      Per-class accuracy (tame=1, native=2) by ecodistrict,
#      combined across years, using grouped rasters (pli_* → pli).
#      Outputs required fields: raster_name, ecodistrict_id,
#      recall_native, recall_tame
# ============================================================

#by_raster_ecodistrict_recall_grouped <- eval_data2 %>%
#  mutate(raster = raster_group) %>%            # collapse pli_* → "pli"
#  group_by(raster, ECODISTRIC) %>%
#  group_modify(~{
#    c3  <- compute_confusion_3(.x)
#    met <- compute_metrics_3(c3$tab)
    # only keep the per-class recalls for 1 (tame) & 2 (native)
#    tibble(
#      recall_tame   = met$Recall_1,
#      recall_native = met$Recall_2
#    )
#  }) %>%
#  ungroup() %>%
#  transmute(
#    raster_name    = raster,
#    ecodistrict_id = ECODISTRIC,
#    recall_native,
#    recall_tame
#  ) %>%
#  arrange(raster_name, ecodistrict_id)

# Write the grouped output
#readr::write_csv(
#  by_raster_ecodistrict_recall_grouped,
#  file.path(out_dir, "raster_ecodistrict_recall_by_class_grouped_for_consensus.csv")
#)

by_raster_ecodistrict_metrics_grouped <- eval_data2 %>%
  mutate(raster = raster_group) %>%
  group_by(raster, ECODISTRIC) %>%
  group_modify(~{
    # Build 2×2 confusion (tame=1, native=2) ignoring pred outside {1,2}
    c2 <- compute_confusion_2_12(.x)
    tab <- c2$tab
    
    # Extract all metrics using the extended metrics built earlier
    # (MCC, kappa, balanced accuracy, class-specific precision/recall/F1, etc.)
    # Reuse the same logic by wrapping tab in a tibble that matches summarize_2x2_block.
    
    # Compute everything:
    a <- as.numeric(tab["1","1"])
    b <- as.numeric(tab["1","2"])
    c <- as.numeric(tab["2","1"])
    d <- as.numeric(tab["2","2"])
    
    n    <- a + b + c + d
    row1 <- a + b
    row2 <- c + d
    col1 <- a + c
    col2 <- b + d
    
    safe <- safe_div
    
    # Basic rates
    oa  <- safe(a + d, n)
    r1  <- safe(a, row1)
    r2  <- safe(d, row2)
    p1  <- safe(a, col1)
    p2  <- safe(d, col2)
    
    sp1 <- safe(d, d + c)
    sp2 <- safe(a, a + b)
    
    fpr1 <- 1 - sp1
    fpr2 <- 1 - sp2
    fnr1 <- 1 - r1
    fnr2 <- 1 - r2
    
    f1_1 <- ifelse((p1 + r1) == 0, NA_real_, 2*p1*r1/(p1 + r1))
    f1_2 <- ifelse((p2 + r2) == 0, NA_real_, 2*p2*r2/(p2 + r2))
    
    bal_acc <- mean(c(r1, r2), na.rm = TRUE)
    
    # MCC
    mcc_num <- (a*d) - (b*c)
    mcc_den <- sqrt((a+b)*(a+c)*(d+b)*(d+c))
    mcc     <- safe(mcc_num, mcc_den)
    
    # Kappa
    pe <- safe(row1*col1 + row2*col2, n*n)
    kappa2 <- safe(oa - pe, 1 - pe)
    
    # Prevalence
    prev1 <- safe(row1, n)
    prev2 <- safe(row2, n)
    
    # Likelihood ratios
    lrpos1 <- safe(r1, fpr1)
    lrneg1 <- safe(fnr1, sp1)
    lrpos2 <- safe(r2, fpr2)
    lrneg2 <- safe(fnr2, sp2)
    
    tibble(
      N_in_2x2        = n,
      Pred_outside_12 = c2$pred_other,
      
      OA_2x2_1v2            = oa,
      Balanced_Accuracy_2x2 = bal_acc,
      MCC_1v2               = mcc,
      Kappa_2x2             = kappa2,
      Prevalence_1          = prev1,
      Prevalence_2          = prev2,
      
      # Class 1: tame
      Recall_1_2x2      = r1,
      Precision_1_2x2   = p1,
      Specificity_1_2x2 = sp1,
      FPR_1_2x2         = fpr1,
      FNR_1_2x2         = fnr1,
      F1_1_2x2          = f1_1,
      LRplus_1          = lrpos1,
      LRminus_1         = lrneg1,
      
      # Class 2: native
      Recall_2_2x2      = r2,
      Precision_2_2x2   = p2,
      Specificity_2_2x2 = sp2,
      FPR_2_2x2         = fpr2,
      FNR_2_2x2         = fnr2,
      F1_2_2x2          = f1_2,
      LRplus_2          = lrpos2,
      LRminus_2         = lrneg2,
      
      Confuse_1_as_2 = b,
      Confuse_2_as_1 = c
    )
  }) %>%
  ungroup() %>%
  transmute(
    raster_name    = raster,
    ecodistrict_id = ECODISTRIC,
    pick(everything())
  ) %>%
  arrange(raster_name, ecodistrict_id)

readr::write_csv(
  by_raster_ecodistrict_metrics_grouped,
  file.path(out_dir, "raster_ecodistrict_metrics_grouped_for_consensus.csv")
)

# ------------------------------------------------------------
# (Optional) Atomic version: treat each raster separately,
# including each PLI variant individually; still combined across years.
# Produces the same four columns.
# ------------------------------------------------------------
by_raster_ecodistrict_recall_atomic <- eval_data2 %>%
  group_by(raster, ECODISTRIC) %>%
  group_modify(~{
    c3  <- compute_confusion_3(.x)
    met <- compute_metrics_3(c3$tab)
    tibble(
      recall_tame   = met$Recall_1,
      recall_native = met$Recall_2
    )
  }) %>%
  ungroup() %>%
  transmute(
    raster_name    = raster,
    ecodistrict_id = ECODISTRIC,
    recall_native,
    recall_tame
  ) %>%
  arrange(raster_name, ecodistrict_id)

readr::write_csv(
  by_raster_ecodistrict_recall_atomic,
  file.path(out_dir, "raster_ecodistrict_recall_by_class_atomic.csv")
)  

# ============================================================
# VISUALIZATION: MCC vs F1 imbalance per ecodistrict
# ============================================================

plot_df <- by_raster_ecodistrict_metrics_grouped %>%  # use updated composite table
  mutate(
    F1_gap = F1_2_2x2 - F1_1_2x2,           # positive = favoring native
    F1_ratio = F1_2_2x2 / (F1_1_2x2 + 1e-6) # >1 = native bias
  )

# Scatter: MCC vs F1_native / F1_tame balance
p1 <- ggplot(plot_df, aes(x = MCC_1v2, y = F1_gap, color = raster_name)) +
  geom_point(alpha=0.7) +
  geom_hline(yintercept=0, linetype="dashed") +
  labs(
    title = "Imbalance Between Native and Tame F1 vs MCC",
    x = "MCC (global reliability)",
    y = "F1_native – F1_tame (positive = native-favoring)"
  ) +
  theme_minimal()

ggsave(file.path(out_dir, "viz_mcc_vs_f1_gap.png"), p1, width=7, height=5, dpi=200)

#Warning message:
#Removed 800 rows containing missing values or values outside the scale range (`geom_point()`). 

# Ratio plot: native‑favoring vs tame‑favoring bias
p2 <- ggplot(plot_df, aes(x = MCC_1v2, y = F1_ratio, color = raster_name)) +
  geom_point(alpha=0.7) +
  geom_hline(yintercept=1, linetype="dashed") +
  scale_y_log10() +
  labs(
    title = "F1 Ratio (Native/Tame) vs MCC",
    x = "MCC",
    y = "F1_native / F1_tame (log scale)"
  ) +
  theme_minimal()

ggsave(file.path(out_dir, "viz_mcc_vs_f1_ratio.png"), p2, width=7, height=5, dpi=200)

#Warning message:
#Removed 800 rows containing missing values or values outside the scale range (`geom_point()`). 

# Map of imbalance per ecodistrict (optional)
eco_plot_df <- eco_sf %>%
  left_join(plot_df, by = c("ECODISTRIC" = "ecodistrict_id"))

p3 <- ggplot(eco_plot_df) +
  geom_sf(aes(fill = F1_gap), color = NA) +
  scale_fill_gradient2(low="#ca0020", mid="white", high="#0571b0", midpoint=0) +
  labs(title = "Spatial Pattern of Native–Tame F1 Imbalance") +
  theme_minimal()

ggsave(file.path(out_dir, "map_f1_imbalance.png"), p3, width=8, height=6, dpi=200)
