
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
acc_cm_file <- "R:/CanadianGrasslands/WG_Classification/results/accuracy_assessment/native_v_tame/compiled_derived/accuracy_assessment_data_binary_pli_sep_compiled_23Jun2026.csv"
acc_ci_file <- "R:/CanadianGrasslands/WG_Classification/results/accuracy_assessment/native_v_tame/crop_ins_derived/accuracy_assessment_data_binary_pli_sep_cropinsV2_lenient_23Jun2026.csv"
metadata_file <- "R:/CanadianGrasslands/WG_Classification/github/tomorrows-prairies/processeddatasets/Grass/metadata_full_extent_nattam_pli_sep.csv"

out_dir <- "R:/CanadianGrasslands/WG_Classification/results/accuracy_assessment/native_v_tame/combined_CI_and_compiled"

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
acc_ci <- acc_ci %>%
  mutate(Training = "")

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
readr::write_csv(combined_out, file.path(out_dir, "accuracy_assessment_data_natvtam_pli_sep_CIandCM_CIlenient_23Jun2026.csv"))

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
    N_evaluated        = as.integer(n),
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
out_xlsx <- file.path(out_dir, "nat_v_tame_v_nongrass_accuracy_assessment_results_23Jun2026_lenient.xlsx")

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
                 file.path(out_dir, "accuracy_by_ecodistrict_grouped_lenient.csv"))

readr::write_csv(accuracy_by_ecodistrict_atomic,
                 file.path(out_dir, "accuracy_by_ecodistrict_atomic_lenient.csv"))

readr::write_csv(accuracy_by_ecodistrict_year_atomic,
                 file.path(out_dir, "accuracy_by_ecodistrict_year_atomic_lenient.csv"))

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
  key_cols = c("raster", "Ecoprv_n")
) %>%
  rename_with(
    ~ paste0(.x, "_ecoP"),
    -c(raster, Ecoprv_n)
  )

# ---- Mapping table ----
ecodist_to_regions <- eval_data %>%
  distinct(Ecds_ID, Ecrgn_n, Ecoprv_n)

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
    by = c("raster", "Ecoprv_n")
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
  ) %>%
  mutate(
    across(
      all_of(metric_cols),
      ~ case_when(
        source_level == "ecodistrict" ~ .data[[paste0(cur_column(), "_ecoD")]],
        source_level == "ecoregion" ~ .data[[paste0(cur_column(), "_ecoR")]],
        TRUE ~ .data[[paste0(cur_column(), "_ecoP")]]
      ),
      .names = "{.col}"
    )
  )

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
    by = c("raster", "Ecoprv_n")
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
  ) %>%
  mutate(
    across(
      all_of(metric_cols),
      ~ case_when(
        source_level == "ecodistrict" ~ .data[[paste0(cur_column(), "_ecoD")]],
        source_level == "ecoregion" ~ .data[[paste0(cur_column(), "_ecoR")]],
        TRUE ~ .data[[paste0(cur_column(), "_ecoP")]]
      ),
      .names = "{.col}"
    )
  )

# ----------------------------
# Export multiclass fallback results
# ----------------------------

readr::write_csv(
  accuracy_with_fallback_atomic,
  file.path(
    out_dir,
    "accuracy_by_ecodistrict_atomic_multiclass_fallback.csv"
  )
)

readr::write_csv(
  accuracy_with_fallback_grouped,
  file.path(
    out_dir,
    "accuracy_by_ecodistrict_grouped_multiclass_fallback.csv"
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