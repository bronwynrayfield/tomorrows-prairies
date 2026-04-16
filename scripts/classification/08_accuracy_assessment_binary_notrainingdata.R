
####################################################################################
##         PRELIMINARY ACCURACY ASSESSMENT: BINARY (NATIVE AND TAME POINTS)       ##
####################################################################################

##----------------------- Alison Long - February 2026 -------------------------

# step 6: assess accuracy of datasets mapping grasslands by comparing them
#         to ground-truthed point data for 2015 to present

# ============================
# Accuracy Assessment by Raster
# ============================

# Description:
#   Compares raster classifications against ground-truth labels with time-aligned filtering.
#   For each raster, we only evaluate points where:
#       - raster value is present (not NA), and
#       - raster_year <= point's ground-truth year (Yer_Cll) if LC_Cl_N == NATIVE, and
#       - raster_year == point's ground-truth year (Yer_Cll) if LC_Cl_N == TAME, and
#       - source was not used to train the dataset
#   Outputs summary metrics, confusion counts, optional per-year breakdown, and a plot.

# ----------------------------
# 1) Setup
# ----------------------------

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

# ----------------------------
# 2) File paths 
# ----------------------------
points_excel_path <- "R:/CanadianGrasslands/WG_Classification/processed/groundtruthed_data_points/groundtruth_points_with_grass_raster_values_ecodist_ge2015_clean.xlsx"
metadata_csv_path <- "R:/CanadianGrasslands/WG_Classification/github/tomorrows-prairies/processeddatasets/Grass/metadata_full_extent_binary.csv"

out_dir <- "R:/CanadianGrasslands/WG_Classification/results/accuracy_assessment/binary"

# ----------------------------
# 3) Configuration
# ----------------------------
raster_cols <- c(
  "aci","luts10","luts30","cec","ngi","cgr","eccc", "esa","glc",
  "cnvc","lcc","pudl","pli","lcr","gpw","plow","plvi","hfi","gcc","dw","esri"
  )

# ----------------------------
# 4) Load data
# ----------------------------
points <- readxl::read_excel(points_excel_path)

# Basic checks
required_cols <- c("LC_grass", "Yer_Cll", "LC_Cl_N")  # <-- include LC_Cl_N
missing_req <- setdiff(required_cols, names(points))
if (length(missing_req) > 0) {
  stop("Missing required columns in points file: ", paste(missing_req, collapse = ", "))
}

# Coerce ground-truth fields and LC_Cl_N
points <- points %>%
  mutate(
    LC_grass = suppressWarnings(as.integer(as.numeric(LC_grass))),
    Yer_Cll  = suppressWarnings(as.integer(as.numeric(Yer_Cll))),
    # Normalize LC_Cl_N to uppercase, trimmed, blanks -> NA
    LC_Cl_N = {
      z <- as.character(LC_Cl_N)
      z <- stringr::str_trim(toupper(z))
      z[z %in% c("", "NA", "NAN", "NULL", "MISSING")] <- NA_character_
      z
    }
  )

# Normalize Source if present (lowercase, trimmed, blanks -> NA)
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

# Identify which raster columns are present in the points file
present_rasters <- intersect(raster_cols, names(points))
if (length(present_rasters) == 0) {
  stop("None of the expected raster columns were found in the points file.")
}

# ----------------------------
# 5) Load metadata and normalize to (raster, year)
# ----------------------------
meta_raw <- readr::read_csv(metadata_csv_path, show_col_types = FALSE)

if (all(c("raster", "year") %in% names(meta_raw))) {
  meta <- meta_raw %>%
    transmute(
      raster = str_trim(tolower(as.character(raster))),
      year   = suppressWarnings(as.integer(as.numeric(year)))
    )
} else if (nrow(meta_raw) == 1 && all(raster_cols %in% names(meta_raw))) {
  # Wide format: one row, columns named by rasters with year values in the single row
  meta <- tibble(
    raster = raster_cols,
    year   = map_int(raster_cols, ~ suppressWarnings(as.integer(as.numeric(meta_raw[[.x]][1]))))
  ) %>%
    mutate(raster = tolower(raster))
} else {
  stop(
    "Could not detect metadata columns. Provide a CSV with columns: 'raster','year', ",
    "OR a single-row wide table with raster columns named: ",
    paste(raster_cols, collapse = ", ")
  )
}

# Build raster -> year map
raster_year_map <- meta %>%
  filter(!is.na(year)) %>%
  distinct(raster, year) %>%
  deframe()

# ----------------------------
# 6) Reshape points to long and join years
# ----------------------------

# Coerce all raster columns to 0/1/NA_integer_ BEFORE pivot_longer
points_norm <- points %>%
  mutate(
    across(
      all_of(present_rasters),
      ~ {
        x <- as.character(.)
        x <- str_trim(x)
        x[x %in% c("", "NA", "NaN", "null", "NULL")] <- NA_character_
        # Keep only 0/1; anything else becomes NA
        x_num <- ifelse(x %in% c("0", "1"), as.integer(x), NA_integer_)
        x_num
      }
    )
  )

# Long format since all raster cols share integer type; include LC_Cl_N
long <- points_norm %>%
  select(all_of(required_cols), Source_norm, all_of(present_rasters)) %>%
  pivot_longer(
    cols = all_of(present_rasters),
    names_to = "raster",
    values_to = "pred"
  ) %>%
  mutate(
    raster = tolower(raster),
    raster_year = unname(as.integer(raster_year_map[raster])),
    raster_year = ifelse(is.na(raster_year), NA_integer_, raster_year)
  )

# ----------------------------
# 7) Filter rows according to class-based time rules
# ----------------------------
# Rules:
#   - raster value present (pred not NA)
#   - if LC_Cl_N == NATIVE: raster_year <= Yer_Cll
#   - if LC_Cl_N == TAME:   raster_year == Yer_Cll
#   - exclude LC_Cl_N == SHRUB
#   - exclude records where Source was used to train the dataset (Source_norm == raster)
eval_data <- long %>%
  # Must have ground-truth and collection year and class
  filter(!is.na(LC_grass), !is.na(Yer_Cll), !is.na(LC_Cl_N)) %>%
  # Filter out SHRUB rows entirely
  filter(LC_Cl_N != "SHRUB") %>%
  # Raster prediction must be present (exclude NA rasters at that point)
  filter(!is.na(pred)) %>%
  # Must know raster year
  filter(!is.na(raster_year)) %>%
  # Apply time-aligned logic by class
  mutate(
    keep_time = dplyr::case_when(
      LC_Cl_N == "NATIVE" ~ raster_year <= Yer_Cll,
      LC_Cl_N == "TAME"   ~ raster_year == Yer_Cll,
      TRUE                ~ FALSE  # any other class excluded
    )
  ) %>%
  filter(keep_time) %>%
  # Exclude points whose source matches the raster being evaluated
  filter(is.na(Source_norm) | Source_norm != raster) %>%
  select(-keep_time)

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

# ----------------------------
# 9) Per-raster summary
# ----------------------------
summary_by_raster <- eval_data %>%
  group_by(raster, raster_year) %>%
  group_modify(~ {
    counts <- compute_counts(.x)
    metrics <- compute_metrics(counts$TP, counts$TN, counts$FP, counts$FN)
    bind_cols(counts, metrics)
  }) %>%
  ungroup() %>%
  arrange(desc(Overall_Accuracy))

# ----------------------------
# 10) Confusion counts (explicit table)
# ----------------------------
confusion_counts <- summary_by_raster %>%
  select(raster, raster_year, TP, TN, FP, FN, N_evaluated)

# ----------------------------
# 11) Per ground-truth year breakdown
# ----------------------------
per_year_summary <- eval_data %>%
  group_by(raster, raster_year, Yer_Cll) %>%
  group_modify(~ {
    counts <- compute_counts(.x)
    metrics <- compute_metrics(counts$TP, counts$TN, counts$FP, counts$FN)
    bind_cols(counts, metrics)
  }) %>%
  ungroup() %>%
  arrange(raster, Yer_Cll)

# ----------------------------
# 12) Save outputs
# ----------------------------
readr::write_csv(summary_by_raster, file.path(out_dir, "raster_accuracy_summary.csv"))
readr::write_csv(confusion_counts,  file.path(out_dir, "raster_confusion_matrices.csv"))
readr::write_csv(per_year_summary,  file.path(out_dir, "raster_per_year_accuracy.csv"))

# ----------------------------
# 13) Simple plot of overall accuracy per raster
# ----------------------------
p <- ggplot(summary_by_raster, aes(x = reorder(raster, Overall_Accuracy), y = Overall_Accuracy)) +
  geom_col(fill = "#2d7f5e") +
  coord_flip() +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.02))) +
  labs(title = "Overall Accuracy by Raster",
       x = "Raster",
       y = "Overall Accuracy") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

ggsave(filename = file.path(out_dir, "raster_accuracy_plot.png"), 
       p, width = 7, height = 5, dpi = 200)

# ----------------------------
# 14) Console previews
# ----------------------------
cat("\n=== Summary by Raster ===\n")
print(summary_by_raster, n = Inf)

cat("\n=== Confusion Counts ===\n")
print(confusion_counts, n = Inf)

cat("\n=== Per-Year Summary (first 30 rows) ===\n")
print(per_year_summary %>% slice_head(n = 30))

########################################################################
########################################################################

# Dataset may contain only positives (LC_grass = 1) — omission-type metrics are still useful.

# ============================
# Positive-only evaluation
# ============================

# Helper: counts and positive-only metrics
compute_pos_only_counts <- function(df) {
  # df has LC_grass == 1 by design; we still guard
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

# Summary per raster (positive-only metrics) — uses the already filtered eval_data
summary_pos_only <- eval_data %>%
  # Ensure we only consider LC_grass == 1
  filter(LC_grass == 1) %>%
  group_by(raster, raster_year) %>%
  group_modify(~ compute_pos_only_counts(.x)) %>%
  ungroup() %>%
  arrange(desc(Recall_TPR))

# Optional per ground-truth year breakdown
per_year_pos_only <- eval_data %>%
  filter(LC_grass == 1) %>%
  group_by(raster, raster_year, Yer_Cll) %>%
  group_modify(~ compute_pos_only_counts(.x)) %>%
  ungroup() %>%
  arrange(raster, Yer_Cll)

# Save to out_dir 
readr::write_csv(summary_pos_only, file.path(out_dir, "raster_posonly_summary.csv"))
readr::write_csv(per_year_pos_only, file.path(out_dir, "raster_posonly_per_year.csv"))

# Simple plot: Recall (TPR) by raster
p_pos <- ggplot(summary_pos_only,
                aes(x = reorder(raster, Recall_TPR), y = Recall_TPR)) +
  geom_col(fill = "#2d7f5e") +
  coord_flip() +
  scale_y_continuous(limits = c(0,1), expand = expansion(mult = c(0, 0.02))) +
  labs(title = "Sensitivity (Recall) on Positive-Only Ground Truth",
       x = "Raster", y = "Recall (TPR)") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(out_dir, "raster_posonly_recall_plot.png"),
       p_pos, width = 7, height = 5, dpi = 200)

####################################################################################
####################################################################################

# ============================
# Positive-only Accuracy by Ecodistrict (Per Raster)
# ============================

# Description:
#   Computes omission-focused metrics per raster and ecodistrict
#   using positive-only ground truth (LC_grass = 1) AND the same
#   class-based time rules (NATIVE <= Yer_Cll; TAME == Yer_Cll)
#   and exclusion of LC_Cl_N == SHRUB.
#
#   Filters:
#     - LC_grass == 1
#     - LC_Cl_N != "SHRUB"
#     - raster value present (not NA)
#     - raster_year known
#     - Time rule by class (as above)
#     - ECODISTRIC not missing for ecodistrict-level summaries
#
# Outputs:
#   - CSV with TP, FN, Recall (TPR), Omission Rate (FNR), and counts per raster × ecodistrict
#   - CSV with counts-only per raster × ecodistrict
#
# ----------------------------
# Reload points and rebuild long data for ecodistrict metrics (independent pass)
# ----------------------------
points_raw <- readxl::read_excel(points_excel_path)

required_cols_eco <- c("LC_grass", "Yer_Cll", "ECODISTRIC", "LC_Cl_N")
missing_req_eco <- setdiff(required_cols_eco, names(points_raw))
if (length(missing_req_eco) > 0) {
  stop("Missing required columns in points file: ", paste(missing_req_eco, collapse = ", "))
}

present_rasters_eco <- intersect(raster_cols, names(points_raw))
if (length(present_rasters_eco) == 0) {
  stop("None of the expected raster columns were found in the points file.")
}

# Normalize LC_grass, Yer_Cll, ECODISTRIC, LC_Cl_N
points2 <- points_raw %>%
  mutate(
    # LC_grass to 0/1/NA
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
    # Yer_Cll to integer year
    Yer_Cll = {
      y <- as.character(Yer_Cll)
      y <- str_trim(y)
      y[y %in% c("", "na", "nan", "null", "missing")] <- NA_character_
      suppressWarnings(as.integer(as.numeric(y)))
    },
    # Normalize ecodistrict labels (trim/standardize case)
    ECODISTRIC = {
      e <- as.character(ECODISTRIC)
      e <- str_trim(e)
      e[e %in% c("", "NA", "NaN", "null", "NULL", "missing")] <- NA_character_
      e
    },
    # Normalize LC_Cl_N to uppercase and trimmed
    LC_Cl_N = {
      z <- as.character(LC_Cl_N)
      z <- stringr::str_trim(toupper(z))
      z[z %in% c("", "NA", "NAN", "NULL", "MISSING")] <- NA_character_
      z
    }
  )

# Normalize Source if present; else create placeholder
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

# Build raster -> year map already built (raster_year_map), reuse

# Normalize raster columns and pivot to long
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
    raster_year = unname(as.integer(raster_year_map[raster])),
    ECODISTRIC = ECODISTRIC
  )

# Apply evaluation filters (positive-only + class-based time rules + no SHRUB)
eval_data2 <- long2 %>%
  filter(
    LC_grass == 1,
    !is.na(Yer_Cll),
    !is.na(LC_Cl_N),
    LC_Cl_N != "SHRUB",
    !is.na(pred),
    !is.na(raster_year)
  ) %>%
  mutate(
    keep_time = dplyr::case_when(
      LC_Cl_N == "NATIVE" ~ raster_year <= Yer_Cll,
      LC_Cl_N == "TAME"   ~ raster_year == Yer_Cll,
      TRUE                ~ FALSE
    )
  ) %>%
  filter(keep_time) %>%
  # exclude points whose Source matches the raster being evaluated
  filter(is.na(Source_norm) | Source_norm != raster) %>%
  select(-keep_time)

# ----------------------------
# Summaries by raster × ecodistrict (positive-only metrics)
# ----------------------------
# Drop rows with missing ecodistrict for ecodistrict-level outputs
eval_by_eco <- eval_data2 %>%
  filter(!is.na(ECODISTRIC))

# Compute counts and metrics
by_raster_ecodistrict <- eval_by_eco %>%
  group_by(raster, raster_year, ECODISTRIC) %>%
  summarise(
    TP = sum(pred == 1, na.rm = TRUE),
    FN = sum(pred == 0, na.rm = TRUE),
    N_pos_evaluated = TP + FN,
    Recall_TPR = ifelse(N_pos_evaluated > 0, TP / N_pos_evaluated, NA_real_),
    OmissionRate_FNR = ifelse(N_pos_evaluated > 0, FN / N_pos_evaluated, NA_real_),
    .groups = "drop"
  ) %>%
  # Remove groups with zero evaluable points (filter out NA metrics)
  filter(N_pos_evaluated > 0) %>%
  arrange(raster, ECODISTRIC)

# A counts-only table
counts_by_raster_ecodistrict <- by_raster_ecodistrict %>%
  select(raster, raster_year, ECODISTRIC, N_pos_evaluated) %>%
  arrange(raster, ECODISTRIC)

# Save outputs
readr::write_csv(
  by_raster_ecodistrict,
  file.path(out_dir, "raster_posonly_by_ecodistrict.csv")
)

readr::write_csv(
  counts_by_raster_ecodistrict,
  file.path(out_dir, "raster_posonly_counts_by_ecodistrict.csv")
)

# Console previews
cat("\n=== Positive-only metrics by raster × ecodistrict (first 30 rows) ===\n")
print(by_raster_ecodistrict %>% slice_head(n = 30))

cat("\n=== Counts per raster × ecodistrict ===\n")
print(counts_by_raster_ecodistrict %>% slice_head(n = 30))

# ----------------------------
# FOR CONSENSUS MAPPING: Combine across years per raster × ecodistrict (micro-average)
# ----------------------------
# Pool TP and FN across years first, then compute recall on pooled counts.
# This avoids weighting each year equally (macro-averaging), which can bias
# results when some years have very few points.

combined_by_raster_ecodistrict <- eval_by_eco %>%
  # eval_by_eco already has the correct filters and includes ALL years
  group_by(raster, ECODISTRIC) %>%
  summarise(
    TP = sum(pred == 1, na.rm = TRUE),
    FN = sum(pred == 0, na.rm = TRUE),
    N_pos_evaluated = TP + FN,
    recall = ifelse(N_pos_evaluated > 0, TP / N_pos_evaluated, NA_real_),
    .groups = "drop"
  ) %>%
  # Keep only combos with at least 1 evaluable positive point
  filter(N_pos_evaluated > 0) %>%
  arrange(raster, ECODISTRIC)

# Prepare the minimalist CSV the user requested
combined_ecodistrict_csv <- combined_by_raster_ecodistrict %>%
  transmute(
    raster_name    = raster,
    ecodistrict_id = ECODISTRIC,
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

##############################
# ============================
# Plot the ecodistrict-level accuracy for each raster
# ============================

# Paths 
metrics_csv <- file.path(out_dir, "raster_posonly_by_ecodistrict.csv")

# Ecodistrict boundary layer
ecodistrict_path <- "R:/CanadianGrasslands/processeddatasets/boundaries/national-ecological-framework-for-canada/Ecodistricts/ecodistricts_clean.shp"
eco_name_col     <- "ECODISTRIC"          # field in the boundary layer that matches ECODISTRIC in the CSV

# Optional: set a map CRS (or keep NA to use the layer's native CRS)
map_crs <- 5070  # e.g., 5070 (NAD83 / Conus Albers) or 3857 for Web Mercator

# Load the per-raster × ecodistrict metrics
by_raster_ecodistrict <- readr::read_csv(metrics_csv, show_col_types = FALSE)

# Expecting columns:
# raster, raster_year, ECODISTRIC, TP, FN, N_pos_evaluated, Recall_TPR, OmissionRate_FNR
req_cols <- c("raster","ECODISTRIC","N_pos_evaluated","Recall_TPR")
missing <- setdiff(req_cols, names(by_raster_ecodistrict))
if (length(missing) > 0) {
  stop("Missing required columns in by_raster_ecodistrict CSV: ", paste(missing, collapse = ", "))
}

# Keep only ecodistricts with >= 1 evaluable point
by_raster_ecodistrict <- by_raster_ecodistrict %>%
  filter(!is.na(ECODISTRIC), N_pos_evaluated > 0)

# Load ecodistrict boundaries (sf)
eco_sf <- sf::st_read(ecodistrict_path, quiet = TRUE)

if (!eco_name_col %in% names(eco_sf)) {
  stop("The ecodistrict layer does not contain the expected name column: ", eco_name_col)
}

# Optional reprojection
if (!is.na(map_crs)) {
  eco_sf <- sf::st_transform(eco_sf, map_crs)
}

# Light-gray background (all ecodistricts)
bg_sf <- eco_sf

# Helper: build joined sf for a single raster
make_map_data <- function(raster_code) {
  m <- by_raster_ecodistrict %>%
    filter(raster == raster_code)
  
  # Left-join metrics onto polygons:
  # eco_sf[[eco_name_col]] <-> m$ECODISTRIC
  map_df <- eco_sf %>%
    left_join(m, by = setNames("ECODISTRIC", eco_name_col)) %>%
    mutate(has_data = !is.na(N_pos_evaluated) & N_pos_evaluated > 0)
  
  map_df
}

# Fill scale for Recall
recall_scale <- scale_fill_gradientn(
  colours = c("#f7fbff", "#c6dbef", "#6baed6", "#2171b5", "#08306b"),
  limits = c(0, 1),
  oob = squish,
  na.value = "grey90",
  name = "Recall (TPR)"
)

# The list of rasters to map (only those in the CSV)
rasters_to_map <- sort(unique(by_raster_ecodistrict$raster))

# Create and save one map per raster
for (r in rasters_to_map) {
  map_df <- make_map_data(r)
  
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
      title = paste0("Sensitivity (Recall) by Ecodistrict — ", toupper(r)),
      subtitle = "Ecodistricts shown only where this raster had evaluable positive points (LC_grass=1)\nTime rule: NATIVE ≤ Yer_Cll; TAME = Yer_Cll; NA rasters excluded; SHRUB excluded",
      caption = "Light gray = ecodistricts with no evaluated points for this raster"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "right"
    )
  
  ggsave(
    filename = file.path(out_dir, paste0("map_recall_by_ecodistrict_", r, ".png")),
    plot = p,
    width = 9, height = 7, dpi = 200
  )
}

cat("\nSaved one map per raster to: ", out_dir, "\n")
