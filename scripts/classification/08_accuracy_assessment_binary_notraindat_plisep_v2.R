
####################################################################################
##         PRELIMINARY ACCURACY ASSESSMENT: BINARY (NATIVE AND TAME POINTS)       ##
####################################################################################

##----------------------- Alison Long - February 2026 -------------------------

# step 8: assess accuracy of datasets mapping grasslands by comparing them
#         to ground-truthed point data for 2015 to present

## v2 is updated with the version of the script from the "check training" file

## CURRENTLY NO MCC?

## NEED TO REMOVE MORE TRAINING DATA FOR HFI AND ACI??

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
points_excel_path <- "R:/CanadianGrasslands/WG_Classification/processed/groundtruthed_data_points/.xlsx"
metadata_csv_path <- "R:/CanadianGrasslands/WG_Classification/github/tomorrows-prairies/processeddatasets/Grass/metadata_full_extent_binary_pli_sep.csv"

out_dir <- "R:/CanadianGrasslands/WG_Classification/results/accuracy_assessment/binary"

# ----------------------------
# 3) Configuration
# ----------------------------
# IMPORTANT: Handle PLI as four separate rasters in the points file and metadata.
pli_variants <- c("pli_ap", "pli_cu", "pli_mg", "pli_mm")

# Expected raster columns 
raster_cols <- c(
  "aci","luts10","luts30","cec","ngi","cgr","eccc", "esa","glc",
  "cnvc","lcc","pudl",
  pli_variants,            # four atomic PLI rasters
  "lcr","gpw","plow","plvi","hfi","gcc","dw","esri","mgi", "pglc")

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
  col_types = c("numeric", "text", "text", "numeric", rep("guess", 10), rep("numeric", 27))
)

names(points)
str(points)
unique(points$Source)
table(points$Source, useNA = "ifany")

required_cols <- c("LC_grass", "Yer_Cll", "LC_Cl_N")
missing_req <- setdiff(required_cols, names(points))
if (length(missing_req) > 0) {
  stop("Missing required columns in points file: ", paste(missing_req, collapse = ", "))
}

# Normalize LC_grass to integer 0/1, Yer_Cll to year int, and LC_Cl_N to uppercased labels
points <- points %>%
  mutate(
    LC_grass = suppressWarnings(as.integer(as.numeric(LC_grass))),
    Yer_Cll  = suppressWarnings(as.integer(as.numeric(Yer_Cll))),
    LC_Cl_N  = {
      z <- as.character(LC_Cl_N)
      z <- stringr::str_trim(toupper(z))
      z[z %in% c("", "NA", "NAN", "NULL", "MISSING")] <- NA_character_
      z
    }
  )

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
  select(all_of(required_cols), Source_norm, all_of(present_rasters)) %>%
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
#   - exclude training sources that match atomic raster or grouped alias
eval_data <- long %>%
  filter(!is.na(LC_grass), !is.na(Yer_Cll), !is.na(LC_Cl_N)) %>%
  #filter(LC_Cl_N != "SHRUB") %>%
  filter(!is.na(pred)) %>%
  rowwise() %>%
  mutate(matched_year = compute_matched_year(allowed_years, Yer_Cll, LC_Cl_N)) %>%
  ungroup() %>%
  filter(!is.na(matched_year)) %>%
  filter(is.na(Source_norm) | !(Source_norm %in% c(raster, raster_group)))

# Check filtering
table(eval_data$Source_norm, eval_data$raster, useNA = "ifany")

eval_data %>%
  filter(!is.na(Source_norm)) %>%
  distinct(Source_norm, raster)

eval_data %>%
  count(raster_group) %>%
  arrange(n)

# Collapse lists to strings and export
eval_data_out <- eval_data %>%
  mutate(allowed_years = map_chr(allowed_years, ~ paste(.x, collapse = ";")))

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
  prec <- safe_div(tp, tp + fp)           # PPV for grass class
  rec  <- safe_div(tp, tp + fn)           # Sensitivity/TPR for grass
  spec <- safe_div(tn, tn + fp)           # Specificity/TNR for not-grass
  f1   <- ifelse(is.na(prec) | is.na(rec) | (prec + rec) == 0, NA_real_, 2 * prec * rec / (prec + rec))
  bal  <- mean(c(rec, spec), na.rm = TRUE)
  tibble(
    Overall_Accuracy     = acc,
    Precision_grass      = prec,
    Recall_grass         = rec,
    Specificity_notgrass = spec,
    F1_grass             = f1,
    Balanced_Accuracy    = bal
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
# 13) Plots: Overall accuracy (GROUPED and ATOMIC)
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
  labs(title = "Overall Accuracy by raster (training points per EL)",
       x = "Raster", y = "Overall Accuracy") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(out_dir, "raster_accuracy_plot_EL_pts_excld2.png"),
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
  labs(title = "Overall Accuracy by Raster (Atomic: pli_ap, pli_cu, pli_mg, pli_mm)",
       x = "Raster", y = "Overall Accuracy") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(out_dir, "raster_accuracy_plot_atomic_EL_pts_excld.png"),
       p_atomic, width = 7, height = 5, dpi = 200)

# Warning message:
#Removed 1 row containing missing values or values outside the scale range (`geom_col()`).

# ----------------------------
# 14) Console previews
# ----------------------------
cat("\n=== Summary by Raster (GROUPED) ===\n")
print(summary_by_raster_grouped, n = Inf)

cat("\n=== Summary by Raster (ATOMIC) ===\n")
print(summary_by_raster_atomic, n = Inf)

cat("\n=== Confusion Counts (GROUPED) ===\n")
print(confusion_counts_grouped, n = Inf)

cat("\n=== Confusion Counts (ATOMIC) ===\n")
print(confusion_counts_atomic, n = Inf)

####################################################################################

## Positive-only evaluation (LC_grass = 1)
##
## Purpose:
##   Some rasters may only have positive cases at evaluated points. For such
##   rasters, omission-type metrics (Recall/TPR, Omission Rate/FNR) are useful.
##
## Approach:
##   - Reuse eval_data (already time-filtered and exclusions applied)
##   - Filter to LC_grass == 1
##   - Compute metrics by raster (atomic and grouped) and by Yer_Cll (optional)
####################################################################################

compute_pos_only_counts <- function(df) {
  df <- df %>% filter(LC_grass == 1)
  tp <- sum(df$pred == 1, na.rm = TRUE)
  fn <- sum(df$pred == 0, na.rm = TRUE)
  n  <- tp + fn
  tibble(
    N_pos_evaluated = n,
    TP = tp,
    FN = fn,
    Recall_TPR = ifelse(n == 0, NA_real_, tp / n),
    OmissionRate_FNR = ifelse(n == 0, NA_real_, fn / n)
  )
}

# ----------------------------
# Positive-only summaries (ATOMIC)
# ----------------------------
summary_pos_only_atomic <- eval_data %>%
  filter(LC_grass == 1) %>%
  group_by(raster, matched_year) %>%
  group_modify(~ compute_pos_only_counts(.x)) %>%
  ungroup() %>%
  rename(raster_year = matched_year) %>%
  arrange(desc(Recall_TPR))

per_year_pos_only_atomic <- eval_data %>%
  filter(LC_grass == 1) %>%
  group_by(raster, matched_year, Yer_Cll) %>%
  group_modify(~ compute_pos_only_counts(.x)) %>%
  ungroup() %>%
  rename(raster_year = matched_year) %>%
  arrange(raster, Yer_Cll)

# ----------------------------
# Positive-only summaries (GROUPED legacy)
# ----------------------------
summary_pos_only_grouped <- eval_data %>%
  filter(LC_grass == 1) %>%
  mutate(raster = raster_group) %>%
  group_by(raster) %>%
  group_modify(~ compute_pos_only_counts(.x)) %>%
  ungroup() %>%
  mutate(raster_year = NA_integer_) %>%
  arrange(desc(Recall_TPR))

per_year_pos_only_grouped <- eval_data %>%
  filter(LC_grass == 1) %>%
  mutate(raster = raster_group) %>%
  group_by(raster, Yer_Cll) %>%
  group_modify(~ compute_pos_only_counts(.x)) %>%
  ungroup() %>%
  mutate(raster_year = NA_integer_) %>%
  arrange(raster, Yer_Cll)

# ----------------------------
# Save positive-only outputs
# ----------------------------
# GROUPED keep legacy filenames
readr::write_csv(summary_pos_only_grouped, file.path(out_dir, "raster_posonly_summary.csv"))
readr::write_csv(per_year_pos_only_grouped, file.path(out_dir, "raster_posonly_per_year.csv"))

# ATOMIC with _atomic suffix
readr::write_csv(summary_pos_only_atomic, file.path(out_dir, "raster_posonly_summary_atomic.csv"))
readr::write_csv(per_year_pos_only_atomic, file.path(out_dir, "raster_posonly_per_year_atomic.csv"))

# ----------------------------
# Plot: Positive-only Recall (TPR) (GROUPED and ATOMIC)
# ----------------------------
# GROUPED legacy
p_pos_grouped <- ggplot(summary_pos_only_grouped,
                        aes(x = reorder(raster, Recall_TPR), y = Recall_TPR)) +
  geom_col(fill = "#2d7f5e") +
  coord_flip() +
  scale_y_continuous(limits = c(0,1), expand = expansion(mult = c(0, 0.02))) +
  labs(title = "Sensitivity (Recall) on Positive-Only Ground Truth (Grouped)",
       x = "Raster", y = "Recall (TPR)") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(out_dir, "raster_posonly_recall_plot.png"),
       p_pos_grouped, width = 7, height = 5, dpi = 200)

# ATOMIC
p_pos_atomic <- ggplot(summary_pos_only_atomic,
                       aes(x = reorder(raster, Recall_TPR), y = Recall_TPR)) +
  geom_col(fill = "#2d7f5e") +
  coord_flip() +
  scale_y_continuous(limits = c(0,1), expand = expansion(mult = c(0, 0.02))) +
  labs(title = "Sensitivity (Recall) on Positive-Only Ground Truth (Atomic: pli variants shown separately)",
       x = "Raster", y = "Recall (TPR)") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(out_dir, "raster_posonly_recall_plot_atomic.png"),
       p_pos_atomic, width = 7, height = 5, dpi = 200)

#Warning message:
#Removed 1 row containing missing values or values outside the scale range (`geom_col()`).

####################################################################################

## Positive-only Accuracy by Ecodistrict (Per Raster)
##
## Description:
##   Computes omission-focused metrics per raster and ecodistrict using
##   positive-only ground truth (LC_grass = 1) AND the same time rules,
##   plus SHRUB exclusion and training-source exclusion.
##
## Outputs (both grouped and atomic):
##   - metrics CSV: TP, FN, N_pos_evaluated, Recall_TPR, OmissionRate_FNR
##   - counts-only CSV
##   - maps for grouped and atomic rasters

####################################################################################

# ----------------------------
# Reload points to ensure ecodistrict field present and clean
# ----------------------------
points_raw <- readxl::read_excel(points_excel_path)

# Update with column name that matches points file
points_raw <- points_raw |>
  dplyr::mutate(
    ecodistrict_id = `Ecds_ID`   
  )

required_cols_eco <- c("LC_grass", "Yer_Cll", "ecodistrict_id", "LC_Cl_N")
missing_req_eco <- setdiff(required_cols_eco, names(points_raw))
if (length(missing_req_eco) > 0) {
  stop("Missing required columns in points file: ", paste(missing_req_eco, collapse = ", "))
}

present_rasters_eco <- intersect(raster_cols, names(points_raw))
if (length(present_rasters_eco) == 0) {
  stop("None of the expected raster columns were found in the points file.")
}

# Normalize key fields consistently with earlier steps
points2 <- points_raw %>%
  mutate(
    LC_grass = {
      x <- as.character(LC_grass)
      x <- str_trim(tolower(x))
      x[x %in% c("", "na", "nan", "null", "missing")] <- NA_character_
      x <- dplyr::case_match(
        x,
        "1" ~ "1", "0" ~ "0",
        "yes" ~ "1", "no" ~ "0",
        "true" ~ "1", "false" ~ "0",
        "grass" ~ "1", "not grass" ~ "0",
        .default = x
      )
      ifelse(x %in% c("0", "1"), as.integer(x), NA_integer_)
    },
    Yer_Cll = {
      y <- as.character(Yer_Cll)
      y <- str_trim(y)
      y[y %in% c("", "na", "nan", "null", "missing")] <- NA_character_
      suppressWarnings(as.integer(as.numeric(y)))
    },
    ecodistrict_id = {
      e <- as.character(ecodistrict_id)
      e <- str_trim(e)
      e[e %in% c("", "NA", "NaN", "null", "NULL", "missing")] <- NA_character_
      e
    },
    LC_Cl_N = {
      z <- as.character(LC_Cl_N)
      z <- stringr::str_trim(toupper(z))
      z[z %in% c("", "NA", "NAN", "NULL", "MISSING")] <- NA_character_
      z
    }
  )

# Source normalization
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

# Normalize rasters and pivot to long with allowed_years attached
points2_norm <- points2 %>%
  mutate(
    across(
      all_of(present_rasters_eco),
      ~ {
        x <- as.character(.)
        x <- str_trim(tolower(x))
        x[x %in% c("", "na", "nan", "null", "missing")] <- NA_character_
        x <- dplyr::case_match(
          x,
          "1" ~ "1", "0" ~ "0",
          "yes" ~ "1", "no" ~ "0",
          "true" ~ "1", "false" ~ "0",
          "grass" ~ "1", "not grass" ~ "0",
          .default = x
        )
        ifelse(x %in% c("0", "1"), as.integer(x), NA_integer_)
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

# Apply time rules to LC_grass == 1 only; compute matched_year
eval_data2 <- long2 %>%
  filter(
    LC_grass == 1,
    !is.na(Yer_Cll),
    !is.na(LC_Cl_N),
    LC_Cl_N != "SHRUB",
    !is.na(pred)
  ) %>%
  rowwise() %>%
  mutate(matched_year = compute_matched_year(allowed_years, Yer_Cll, LC_Cl_N)) %>%
  ungroup() %>%
  filter(!is.na(matched_year)) %>%
  filter(is.na(Source_norm) | !(Source_norm %in% c(raster, raster_group)))

# Keep rows having ecodistrict info
eval_by_eco <- eval_data2 %>%
  filter(!is.na(ecodistrict_id))

# ----------------------------
# Summaries by raster × ecodistrict — ATOMIC and GROUPED
# ----------------------------

# ATOMIC ecodistrict metrics (positive-only)
by_raster_ecodistrict_atomic <- eval_by_eco %>%
  group_by(raster, matched_year, ecodistrict_id) %>%
  summarise(
    TP = sum(pred == 1, na.rm = TRUE),
    FN = sum(pred == 0, na.rm = TRUE),
    N_pos_evaluated = TP + FN,
    Recall_TPR = ifelse(N_pos_evaluated > 0, TP / N_pos_evaluated, NA_real_),
    OmissionRate_FNR = ifelse(N_pos_evaluated > 0, FN / N_pos_evaluated, NA_real_),
    .groups = "drop"
  ) %>%
  rename(raster_year = matched_year) %>%
  filter(N_pos_evaluated > 0) %>%
  arrange(raster, ecodistrict_id)

# GROUPED ecodistrict metrics (pli_* -> pli)
by_raster_ecodistrict_grouped <- eval_by_eco %>%
  mutate(raster = raster_group) %>%
  group_by(raster, ecodistrict_id) %>%
  summarise(
    TP = sum(pred == 1, na.rm = TRUE),
    FN = sum(pred == 0, na.rm = TRUE),
    N_pos_evaluated = TP + FN,
    Recall_TPR = ifelse(N_pos_evaluated > 0, TP / N_pos_evaluated, NA_real_),
    OmissionRate_FNR = ifelse(N_pos_evaluated > 0, FN / N_pos_evaluated, NA_real_),
    .groups = "drop"
  ) %>%
  mutate(raster_year = NA_integer_) %>%
  filter(N_pos_evaluated > 0) %>%
  arrange(raster, ecodistrict_id)

# Counts-only views
counts_by_raster_ecodistrict_atomic <- by_raster_ecodistrict_atomic %>%
  select(raster, raster_year, ecodistrict_id, N_pos_evaluated) %>%
  arrange(raster, ecodistrict_id)

counts_by_raster_ecodistrict_grouped <- by_raster_ecodistrict_grouped %>%
  select(raster, raster_year, ecodistrict_id, N_pos_evaluated) %>%
  arrange(raster, ecodistrict_id)

# ----------------------------
# Save ecodistrict outputs
# ----------------------------
# GROUPED (legacy filenames)
readr::write_csv(
  by_raster_ecodistrict_grouped,
  file.path(out_dir, "raster_posonly_by_ecodistrict.csv")
)

readr::write_csv(
  counts_by_raster_ecodistrict_grouped,
  file.path(out_dir, "raster_posonly_counts_by_ecodistrict.csv")
)

# ATOMIC (_atomic suffix)
readr::write_csv(
  by_raster_ecodistrict_atomic,
  file.path(out_dir, "raster_posonly_by_ecodistrict_atomic.csv")
)

readr::write_csv(
  counts_by_raster_ecodistrict_atomic,
  file.path(out_dir, "raster_posonly_counts_by_ecodistrict_atomic.csv")
)

# Console previews
cat("\n=== Positive-only metrics by raster × ecodistrict (GROUPED, first 30 rows) ===\n")
print(by_raster_ecodistrict_grouped %>% slice_head(n = 30))
cat("\n=== Counts per raster × ecodistrict (GROUPED, first 30 rows) ===\n")
print(counts_by_raster_ecodistrict_grouped %>% slice_head(n = 30))

cat("\n=== Positive-only metrics by raster × ecodistrict (ATOMIC, first 30 rows) ===\n")
print(by_raster_ecodistrict_atomic %>% slice_head(n = 30))
cat("\n=== Counts per raster × ecodistrict (ATOMIC, first 30 rows) ===\n")
print(counts_by_raster_ecodistrict_atomic %>% slice_head(n = 30))

####################################################################################

# ----------------------------
# FOR CONSENSUS MAPPING: Combine across years per raster × ecodistrict (micro-average)
# ----------------------------
# Pool TP and FN across years first, then compute recall on pooled counts.
# This avoids weighting each year equally (macro-averaging), which can bias
# results when some years have very few points.

combined_by_raster_ecodistrict <- eval_by_eco %>%
  # eval_by_eco already has the correct filters and includes ALL years
  group_by(raster, ecodistrict_id) %>%
  summarise(
    TP = sum(pred == 1, na.rm = TRUE),
    FN = sum(pred == 0, na.rm = TRUE),
    N_pos_evaluated = TP + FN,
    recall = ifelse(N_pos_evaluated > 0, TP / N_pos_evaluated, NA_real_),
    .groups = "drop"
  ) %>%
  # Keep only combos with at least 1 evaluable positive point
  filter(N_pos_evaluated > 0) %>%
  arrange(raster, ecodistrict_id)

# Prepare the minimalist CSV the user requested
combined_ecodistrict_csv <- combined_by_raster_ecodistrict %>%
  transmute(
    raster_name    = raster,
    ecodistrict_id = ecodistrict_id,
    recall         = recall
  )

# Optional: round for readability (remove if you need full precision)
combined_ecodistrict_csv <- combined_ecodistrict_csv %>%
  mutate(recall = round(recall, 3))

# Save
readr::write_csv(
  combined_ecodistrict_csv,
  file.path(out_dir, "raster_ecodistrict_recall_combined_years_for_consensus_map.csv")
)

# Console preview
cat("\n=== Recall by raster × ecodistrict (combined across years) ===\n")
print(combined_ecodistrict_csv %>% slice_head(n = 30))

#####################################################################################

# ----------------------------
# Maps: Recall by Ecodistrict — GROUPED and ATOMIC
# ----------------------------

# Default CSVs for maps
metrics_csv_grouped <- file.path(out_dir, "raster_posonly_by_ecodistrict.csv")
metrics_csv_atomic  <- file.path(out_dir, "raster_posonly_by_ecodistrict_atomic.csv")

# Ecodistrict boundary layer (must have a name column matching ecodistrict in metrics)
ecodistrict_path <- "R:/CanadianGrasslands/processeddatasets/boundaries/national-ecological-framework-for-canada/Ecodistricts/ecodistricts_clean.shp"
eco_bounds <- sf::st_read(ecodistrict_path, quiet = TRUE)

eco_name_col     <- "ECODISTRIC"   # in the boundary layer

eco_bounds <- eco_bounds |>
  dplyr::mutate(
    ecodistrict_id = .data[[eco_name_col]]
  )

# Optional map CRS (use native CRS if NA)
map_crs <- 5070  # e.g., 5070 (NAD83 / Conus Albers)

# Load GROUPED metrics
by_raster_ecodistrict <- readr::read_csv(metrics_csv_grouped, show_col_types = FALSE)
req_cols <- c("raster","ecodistrict_id","N_pos_evaluated","Recall_TPR")
missing <- setdiff(req_cols, names(by_raster_ecodistrict))
if (length(missing) > 0) {
  stop("Missing required columns in grouped ecodistrict CSV: ", paste(missing, collapse = ", "))
}
by_raster_ecodistrict <- by_raster_ecodistrict %>%
  filter(!is.na(ecodistrict_id), N_pos_evaluated > 0)

# Load boundaries (sf)
eco_sf <- sf::st_read(ecodistrict_path, quiet = TRUE)
if (!eco_name_col %in% names(eco_sf)) {
  stop("Ecodistrict layer missing expected name column: ", eco_name_col)
}
if (!is.na(map_crs)) {
  eco_sf <- sf::st_transform(eco_sf, map_crs)
}
bg_sf <- eco_sf

# Helper to join metrics onto polygons (flexible to grouped or atomic)
make_map_data <- function(raster_code, metrics_df) {
  m <- metrics_df %>% filter(raster == raster_code)
  # Left-join polygons (x) with metrics (y) by ecodistrict name
  # by = c(x_col = y_col) -> use a named character vector:
  join_by <- setNames("ecodistrict_id", eco_name_col)
  map_df <- eco_sf %>%
    left_join(m, by = join_by) %>%
    mutate(has_data = !is.na(N_pos_evaluated) & N_pos_evaluated > 0)
  map_df
}

# Shared fill scale
recall_scale <- scale_fill_gradientn(
  colours = c("#f7fbff", "#c6dbef", "#6baed6", "#2171b5", "#08306b"),
  limits = c(0, 1),
  oob = squish,
  na.value = "grey90",
  name = "Recall (TPR)"
)

# --- GROUPED maps ---
rasters_to_map_grouped <- sort(unique(by_raster_ecodistrict$raster))
for (r in rasters_to_map_grouped) {
  map_df <- make_map_data(r, by_raster_ecodistrict)
  p <- ggplot() +
    geom_sf(data = bg_sf, fill = "grey95", color = "white", size = 0.2) +
    geom_sf(
      data = dplyr::filter(map_df, has_data),
      aes(fill = Recall_TPR),
      color = "white",
      size = 0.2
    ) +
    recall_scale +
    labs(
      title = paste0("Sensitivity (Recall) by Ecodistrict — ", toupper(r), " (Grouped)"),
      subtitle = "Shown where this raster had evaluable positive points (LC_grass=1).\nTime rule: NATIVE ≤ Yer_Cll; TAME = Yer_Cll; NA rasters excluded; SHRUB excluded.",
      caption = "Light gray = ecodistricts with no evaluated points for this raster"
    ) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "right")
  ggsave(
    filename = file.path(out_dir, paste0("map_recall_by_ecodistrict_", r, ".png")),
    plot = p, width = 9, height = 7, dpi = 200
  )
}

# --- ATOMIC maps ---
by_raster_ecodistrict_atomic_csv <- readr::read_csv(metrics_csv_atomic, show_col_types = FALSE) %>%
  filter(!is.na(ecodistrict_id), N_pos_evaluated > 0)

rasters_to_map_atomic <- sort(unique(by_raster_ecodistrict_atomic_csv$raster))
for (r in rasters_to_map_atomic) {
  map_df <- make_map_data(r, by_raster_ecodistrict_atomic_csv)
  p <- ggplot() +
    geom_sf(data = bg_sf, fill = "grey95", color = "white", size = 0.2) +
    geom_sf(
      data = dplyr::filter(map_df, has_data),
      aes(fill = Recall_TPR),
      color = "white",
      size = 0.2
    ) +
    recall_scale +
    labs(
      title = paste0("Sensitivity (Recall) by Ecodistrict — ", toupper(r), " (Atomic)"),
      subtitle = "Shown where this raster had evaluable positive points (LC_grass=1).\nTime rule: NATIVE ≤ Yer_Cll; TAME = Yer_Cll; NA rasters excluded; SHRUB excluded.",
      caption = "Light gray = ecodistricts with no evaluated points for this raster"
    ) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "right")
  ggsave(
    filename = file.path(out_dir, paste0("map_recall_by_ecodistrict_atomic_", r, ".png")),
    plot = p, width = 9, height = 7, dpi = 200
  )
}

cat("\nSaved grouped and atomic ecodistrict maps to: ", out_dir, "\n")