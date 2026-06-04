
####################################################################################
##         PRELIMINARY ACCURACY ASSESSMENT: BINARY (NATIVE AND TAME POINTS)       ##
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
##
## Exclusions:
##   - LC_Cl_N == "SHRUB"
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
points_excel_path <- "R:/CanadianGrasslands/WG_Classification/processed/groundtruthed_data_points/ci_derived/point_raster_extract_grass_v_nongrass.xlsx"
metadata_csv_path <- "R:/CanadianGrasslands/WG_Classification/github/tomorrows-prairies/processeddatasets/Grass/metadata_full_extent_binary_pli_sep_noplow.csv"

out_dir <- "R:/CanadianGrasslands/WG_Classification/results/accuracy_assessment/binary/ci_derived"

# ----------------------------
# 3) Configuration
# ----------------------------
# Handle PLI as four separate rasters in the points file and metadata.
pli_variants <- c("pli_ap", "pli_cu", "pli_mg", "pli_mm")

# Expected raster columns 
raster_cols <- c(
  "aci","luts10","luts30","cec","ngi","cgr","eccc", "esa","glc",
  "cnvc","lcc","pudl",
  pli_variants,            # four atomic PLI rasters
  "lcr","gpw","plvi","hfi","gcc","dw","esri","mgi","pglc")

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

# Remove Plowprint
points <- points %>%
  select(-plow_reproj_algn_intact_masked)

# Rename columns
points <- points %>%
  rename(
    Pt_ID    = Point_ID,
    Sorc_ID  = Source_ID,
    Yer_Cll  = Year_Colle,
    LC_Cl_N  = LC_Class_N,
    Ecods_ID = ecodist,
    Provinc  = PROV,
    Ecrgn_n  = ecoreg,
    aci      = aafc_aci_merged_reproj_algn_binary_noshrub,
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
    pglc     = pglc_reproj_algn_binary)

# Create column with numerical class for ground-truthed cover type
points <- points %>%
  mutate(
    LC_grass = case_when(
      LC_Cl_N %in% c("TAME", "NATIVE") ~ 1,
      LC_Cl_N == "OTHER" ~ 0,
      TRUE ~ NA_real_))

unique(points$LC_grass)
table(points$LC_grass, useNA = "ifany")

required_cols <- c("Pt_ID", "LC_grass", "Yer_Cll", "LC_Cl_N")
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

# Create a canonical evaluation key
#long <- long %>%
#  mutate(
#    raster_eval = raster_group   # <-- canonical evaluation identity
#  )

# ----------------------------
# 7) Filter rows via class-based time rules using MULTIPLE allowed years
# ----------------------------
# Compute a per-record matched_year for per-year reporting:
#   - NATIVE: latest allowed_year <= Yer_Cll (if any; else NA -> exclude)
#   - TAME: Yer_Cll if Yer_Cll is in allowed_years (else NA -> exclude)
#   - OTHER: Yer_Cll if Yer_Cll is in allowed_years (else NA -> exclude)
compute_matched_year <- function(yrs, yer, class_label) {
  yrs <- as.integer(yrs)
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
eval_data <- long %>%
  filter(!is.na(LC_grass), !is.na(Yer_Cll), !is.na(LC_Cl_N)) %>%
  filter(!is.na(pred)) %>%
  rowwise() %>%
  mutate(matched_year = compute_matched_year(allowed_years, Yer_Cll, LC_Cl_N)) %>%
  ungroup() %>%
  filter(!is.na(matched_year)) #%>%
  #filter(is.na(Source_norm) | !(Source_norm %in% c(raster, raster_group)))
  #filter(is.na(Source_norm) | Source_norm != raster_eval) # filter using canonical evaluation identity

# Check filtering
#table(eval_data$Source_norm, eval_data$raster, useNA = "ifany")

#eval_data %>%
#  filter(!is.na(Source_norm)) %>%
#  distinct(Source_norm, raster)

#eval_data %>%
#  count(raster_eval) %>%
#  arrange(n)

# Add column for confusion matrix quadrants
eval_data <- eval_data %>%
  mutate(quad = case_when(
    LC_grass == 1 & pred == 1 ~ "tp",
    LC_grass == 0 & pred == 0 ~ "tn",
    LC_grass == 0 & pred == 1 ~ "fp",
    LC_grass == 1 & pred == 0 ~ "fn"
  ))

# Export
readr::write_csv(eval_data, file.path(out_dir, "accuracy_assessment_data_binary_pli_sep_cropins_3Jun2026.csv"))

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
  
  bal <- mean(c(rec, spec), na.rm = TRUE)
  
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
readr::write_csv(summary_by_raster_grouped, file.path(out_dir, "raster_accuracy_summary2.csv"))
readr::write_csv(confusion_counts_grouped,  file.path(out_dir, "raster_confusion_matrices2.csv"))
readr::write_csv(per_year_summary_grouped,  file.path(out_dir, "raster_per_year_accuracy2.csv"))

# ATOMIC with _atomic suffix
readr::write_csv(summary_by_raster_atomic, file.path(out_dir, "raster_accuracy_summary_atomic2.csv"))
readr::write_csv(confusion_counts_atomic,  file.path(out_dir, "raster_confusion_matrices_atomic2.csv"))
readr::write_csv(per_year_summary_atomic,  file.path(out_dir, "raster_per_year_accuracy_atomic2.csv"))

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
