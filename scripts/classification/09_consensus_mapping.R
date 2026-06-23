
###################################################################################
##                     GRASSLAND CONSENSUS MAPPING ANALYSIS                      ##
###################################################################################

##----------------------- Alison Long - February 2026 -------------------------

# step 9: use data fusion techniques to create a consensus grassland map 
#         to serve as baseline for modelling

# NEED TO ADD ECOREGION-LEVEL ACCURACY FALLBACK FOR ECODISTRICTS W/O ACCURACY DATA

# load packages
library(terra)
library(tidyverse)

# Set the temporary directory to a drive with sufficient space
terraOptions(tempdir = "R:/CanadianGrasslands/WG_classification/scratch")

terraOptions(
  memfrac = 0.6,
  #threads = max(1, parallel::detectCores() - 1),
  tempdir = "R:/CanadianGrasslands/WG_classification/scratch"
)

# GDAL write options (apply to all writeRaster/app calls)
wopt_prob   <- list(datatype = "FLT4S", gdal = c("COMPRESS=LZW", "TILED=YES", "BIGTIFF=YES"))
wopt_byte   <- list(datatype = "INT1U", gdal = c("COMPRESS=LZW", "TILED=YES", "BIGTIFF=YES"))
wopt_count  <- list(datatype = "INT2U", gdal = c("COMPRESS=LZW", "TILED=YES", "BIGTIFF=YES"))
wopt_float  <- list(datatype = "FLT4S", gdal = c("COMPRESS=LZW", "TILED=YES", "BIGTIFF=YES"))

####--- DIRECTORIES ---####

# pilot area
#pa_dir <- "R:/CanadianGrasslands/processeddatasets/pilotarea"

# study area
sa_dir <- "R:/CanadianGrasslands/processeddatasets/studyarea"

# accuracy data
acc_dir <- "R:/CanadianGrasslands/WG_Classification/results/accuracy_assessment"

# grasslands
proc_data_dir <- "R:/CanadianGrasslands/WG_Classification/processed"

# grass raster directories
aci_aligned_dir   <- file.path(proc_data_dir, "aafc_aci/aligned")
luts30_masked_dir <- file.path(proc_data_dir, "aafc_luts_2020_30m_v5/masked")
luts10_masked_dir <- file.path(proc_data_dir, "aafc_luts_2020_10m_v1/masked")
cec_masked_dir    <- file.path(proc_data_dir, "cec_nalcms/masked")
ngi_aligned_dir   <- file.path(proc_data_dir, "cfga_ngi/aligned")
cgr_masked_dir    <- file.path(proc_data_dir, "cgr_gam_adv/masked")
eccc_masked_dir   <- file.path(proc_data_dir, "eccc_grass/masked")
esa_masked_dir    <- file.path(proc_data_dir, "esa_worldcover/masked")
glc_masked_dir    <- file.path(proc_data_dir, "glc_fcs30d/masked")
cnvc_masked_dir   <- file.path(proc_data_dir, "natureserve_cnvc/masked")
lcc_masked_dir    <- file.path(proc_data_dir, "nrcan_lcc/masked")
pudl_aligned_dir  <- file.path(proc_data_dir, "phjv_pudl/aligned")
pli_aligned_dir   <- file.path(proc_data_dir, "sk_pli/aligned")
lcr_masked_dir    <- file.path(proc_data_dir, "statcan_lcr/masked")
gpw_masked_dir    <- file.path(proc_data_dir, "wri_gpw_ann/masked")
plow_masked_dir   <- file.path(proc_data_dir, "wwf_plowprint/masked")
plvi_aligned_dir  <- file.path(proc_data_dir, "ab_plvi/aligned")
hfi_aligned_dir   <- file.path(proc_data_dir, "abmi_hfi/aligned")
gcc_aligned_dir   <- file.path(proc_data_dir, "bc_gcc/aligned")
dw_masked_dir     <- file.path(proc_data_dir, "dw_2025_grass/masked")
esri_masked_dir   <- file.path(proc_data_dir, "esri_rangelands/masked")

# output 
out_dir <- "R:/CanadianGrasslands/WG_Classification/results/consensus_mapping"

####--- FILES ---####

# Study area
sa_file <- file.path(sa_dir, "classification_study_area_aligned_cnvc_Apr2026.tif")

# ecodistricts
ecodist_file   <- "R:/CanadianGrasslands/processeddatasets/boundaries/national-ecological-framework-for-canada/Ecodistricts/ecodistricts_clean.shp"  

# accuracy data
acc_bin_file <- file.path(acc_dir, "binary/raster_ecodistrict_recall_combined_years_for_consensus_map.csv")
acc_nt_file  <- file.path(acc_dir, "native_v_tame/raster_ecodistrict_metrics_grouped_for_consensus.csv")

# binary grass raster files
aci_bin_file    <- file.path(aci_aligned_dir, "aafc_aci_merged_reproj_algn_binary_noshrub.tif")
luts30_bin_file <- file.path(luts30_masked_dir, "luts30_reproj_algn_binary_masked.tif")
luts10_bin_file <- file.path(luts10_masked_dir, "luts10_reproj_algn_binary_masked.tif")
cec_bin_file    <- file.path(cec_masked_dir, "nalcms_reproj_algn_binary_noshrub_masked.tif")
ngi_bin_file    <- file.path(ngi_aligned_dir, "cfga_ngi_merged_reproj_algn_binary.tif")
cgr_bin_file    <- file.path(cgr_masked_dir, "cgr_reproj_algn_binary_masked.tif")
eccc_bin_file   <- file.path(eccc_masked_dir, "eccc_reclass_reproj_algn_binary_masked.tif")
esa_bin_file    <- file.path(esa_masked_dir, "esa_reproj_algn_binary_noshrub_masked.tif")
glc_bin_file    <- file.path(glc_masked_dir, "glc_reproj_algn_binary_noshrub_masked.tif")
cnvc_bin_file   <- file.path(cnvc_masked_dir, "cnvc_algn_binary_masked.tif")
lcc_bin_file    <- file.path(lcc_masked_dir, "lcc_reproj_algn_binary_noshrub_masked.tif")
pudl_bin_file   <- file.path(pudl_aligned_dir, "pudl_reclass_reproj_algn_binary_noshrub.tif")
pli_ap_bin_file <- file.path(pli_aligned_dir, "sk_pli_ap_reproj_algn_binary_noshrub.tif")
pli_cu_bin_file <- file.path(pli_aligned_dir, "sk_pli_cu_reproj_algn_binary_noshrub.tif")
pli_mg_bin_file <- file.path(pli_aligned_dir, "sk_pli_mg_reproj_algn_binary_noshrub.tif")
pli_mm_bin_file <- file.path(pli_aligned_dir, "sk_pli_mm_reproj_algn_binary_noshrub.tif")
lcr_bin_file    <- file.path(lcr_masked_dir, "lcr_reproj_algn_binary_masked.tif")
gpw_bin_file    <- file.path(gpw_masked_dir, "gpw_reproj_algn_binary_noshrub_masked.tif")
plow_bin_file   <- file.path(plow_masked_dir, "plow_reproj_algn_intact_masked.tif")
plvi_bin_file   <- file.path(plvi_aligned_dir, "plvi_reproj_algn_binary_gte2015.tif")
hfi_bin_file    <- file.path(hfi_aligned_dir, "hfi_reproj_algn_binary_gte2015.tif")
gcc_bin_file    <- file.path(gcc_aligned_dir, "gcc_reproj_algn_binary.tif")
dw_bin_file     <- file.path(dw_masked_dir, "dw_reproj_algn_masked_binary.tif")
esri_bin_file   <- file.path(esri_masked_dir, "esri_reproj_algn_masked_binary.tif")

# multi-class grass raster files
aci_nt_file    <- file.path(aci_aligned_dir, "aafc_aci_merged_reproj_algn_natvtam_noshrub.tif")
ngi_nt_file    <- file.path(ngi_aligned_dir, "cfga_ngi_merged_reproj_algn_natvtam_good.tif")
cgr_nt_file    <- file.path(cgr_masked_dir, "cgr_reproj_algn_natvtam_masked.tif")
eccc_nt_file   <- file.path(eccc_masked_dir, "eccc_reclass_reproj_algn_natvtam_masked.tif")
cnvc_nt_file   <- file.path(cnvc_masked_dir, "cnvc_algn_natvtam_masked.tif")
pudl_nt_file   <- file.path(pudl_aligned_dir, "pudl_reclass_reproj_algn_natvtam_noshrub.tif")
pli_ap_nt_file <- file.path(pli_aligned_dir, "sk_pli_ap_reproj_algn_natvtam_noshrub.tif")
pli_cu_nt_file <- file.path(pli_aligned_dir, "sk_pli_cu_reproj_algn_natvtam_noshrub.tif")
pli_mg_nt_file <- file.path(pli_aligned_dir, "sk_pli_mg_reproj_algn_natvtam_noshrub.tif")
pli_mm_nt_file <- file.path(pli_aligned_dir, "sk_pli_mm_reproj_algn_natvtam_noshrub.tif")
gpw_nt_file    <- file.path(gpw_masked_dir, "gpw_reproj_algn_natvtam_noshrub_masked.tif")
plvi_nt_file   <- file.path(plvi_aligned_dir, "plvi_reproj_algn_natvtam_gte2015.tif")
hfi_nt_file    <- file.path(hfi_aligned_dir, "hfi_reproj_algn_binary_gte2015.tif") # tame only - same as binary

####--- READ INTO TERRA ---####

# classification study area (modelling extent + Peace River region)
sa <- rast(sa_file)

# ecodistricts
ecodist <- vect(ecodist_file)

#-------------------------------
# 1) Load binary rasters 
#-------------------------------

bin_files <- c(
  aci_bin_file, luts30_bin_file, luts10_bin_file, cec_bin_file, ngi_bin_file, cgr_bin_file, 
  eccc_bin_file, esa_bin_file, glc_bin_file, cnvc_bin_file, lcc_bin_file, pudl_bin_file, 
  pli_ap_bin_file, pli_cu_bin_file, pli_mg_bin_file, pli_mm_bin_file, lcr_bin_file, 
  gpw_bin_file, plow_bin_file, plvi_bin_file, hfi_bin_file, gcc_bin_file, dw_bin_file,
  esri_bin_file
)

# Load as a SpatRaster stack
r_stack <- rast(bin_files)

# Crop and mask to study area (modelling extent + Peace River region)
r_crop <- crop(r_stack, sa)
r_mask <- mask(r_crop, sa)

# --- Derive short names from file paths ---
# Rule: take the part before the first underscore in the basename.
# Fallback: if there is no underscore, use the basename without extension.
make_short_name <- function(path) {
  bn   <- basename(path)
  stem <- tools::file_path_sans_ext(bn)
  
  # Special case: Saskatchewan PLI layers
  # sk_pli_ap_*, sk_pli_cu_*, etc. → pli_ap, pli_cu
  if (grepl("^sk_pli_", stem)) {
    sub("^sk_(pli_[^_]+).*", "\\1", stem)
    
    # General case: take text before first underscore
  } else if (grepl("_", stem)) {
    sub("^([^_]+)_.*$", "\\1", stem)
    
    # Fallback
  } else {
    stem
  }
}

short_names <- tolower(vapply(bin_files, make_short_name, FUN.VALUE = character(1)))

# Ensure names are unique (in case two files have same prefix)
if (any(duplicated(short_names))) {
  warning("Duplicate short names detected. Making unique with numeric suffixes.")
  short_names <- make.unique(short_names, sep = "_")
}

# Apply to the raster stack
names(r_mask) <- short_names

# Preview
cat("Assigned r_stack names:\n")
print(names(r_mask))


# Optional: harmonize names if accuracy table references basenames
# names(r_stack) <- tools::file_path_sans_ext(basename(raster_files))

#-------------------------------
# 2) Load ecodistrict polygons and accuracy table
#-------------------------------

ecodist <- project(ecodist, crs(r_mask))

# ---- Resolve ecodistrict id field and create a uniform alias ----
id_candidates <- c("ECODISTRIC", "ecodistrict_id")
match_field <- id_candidates[id_candidates %in% names(ecodist)][1]

if (is.na(match_field)) {
  stop("Ecodistrict layer missing an id field. Expected one of: ",
       paste(id_candidates, collapse = ", "))
}

# Create a consistent alias in the polygons
ecodist$ecodistrict_id <- ecodist[[match_field]]

# accuracy_df columns expected: ecodistrict_id, raster_name, recall
accuracy_bin_df <- readr::read_csv(acc_bin_file, show_col_types = FALSE)

# basic checks
req_cols <- c("raster_name", "ecodistrict_id", "recall")
missing_cols <- setdiff(req_cols, names(accuracy_bin_df))
if (length(missing_cols) > 0) {
  stop("CSV is missing required column(s): ", paste(missing_cols, collapse = ", "))
}

# Print unique raster names
unique_rasters <- unique(accuracy_bin_df$raster_name)
print(unique_rasters)

# Check against r_mask names
print(names(r_mask))

# Update accuracy_bin_df names to match
accuracy_bin_df <- accuracy_bin_df %>%
  mutate(
    raster_name = case_when(
      raster_name == "aci"   ~ "aafc",
      raster_name == "cec" ~ "nalcms",
      raster_name == "ngi" ~ "cfga",
      TRUE                    ~ raster_name
    )
  )

# Make sure accuracy_bin_df$raster_name matches names(r_stack)
# accuracy_bin_df$raster_name <- tools::file_path_sans_ext(basename(accuracy_bin_df$raster_name))
stopifnot(all(accuracy_bin_df$raster_name %in% names(r_mask)))

# Pivot accuracy values wide: one row per ecodistrict
acc_bin_wide <- accuracy_bin_df %>%
  pivot_wider(names_from = raster_name, values_from = recall)

# ----------------------
# JOIN USING ECODISTRIC
# ----------------------

if (!"ecodistrict_id" %in% names(ecodist)) {
  stop("ecodistrict_id field is missing from ecodistricts.")
}


if (!"ECODISTRIC" %in% names(ecodist)) {
  stop("The ecodistrict polygon layer is missing field 'ECODISTRIC'.")
}

#ecodist_bin <- merge(ecodist, acc_bin_wide, by = "ecodistrict_id", duplicateGeoms = TRUE)

# Merge: ecodist uses ECODISTRIC; accuracy_bin_df uses ecodistrict_id
ecodist_bin <- merge(
  ecodist,
  acc_bin_wide,
  by.x = "ECODISTRIC",     # polygon field
  by.y = "ecodistrict_id", # CSV field
  duplicateGeoms = TRUE
)

#-------------------------------
# 3) Build per‑raster weight rasters from polygon accuracies
#-------------------------------
weight_list <- vector("list", length = nlyr(r_mask))
names(weight_list) <- names(r_mask)

for (rname in names(r_mask)) {
  # Rasterize the accuracy field for this raster’s name
  w <- rasterize(ecodist_bin, r_mask[[1]], field = rname, touches = TRUE)
  weight_list[[rname]] <- w
}

weight_stack <- rast(weight_list)
# Ensure order matches r_mask
weight_stack <- weight_stack[[names(r_mask)]]

#-------------------------------
# 4) Weighted average that ignores NA but uses 0 as real "not grassland"
#-------------------------------
# present = where the raster has data (0 or 1), i.e., NOT NA
present <- !is.na(r_mask)

# Zero-out weights where the raster has NA (no data) at that pixel
present_weights <- weight_stack * present   # NA or 0 becomes 0 weight here

# Numerator: sum over (value * weight). 0 values contribute 0, 1 contributes weight, NAs contribute nothing.
weighted_sum <- sum(r_mask * present_weights, na.rm = TRUE)

# Denominator: sum of weights from rasters that had data (0 or 1)
weight_total_present <- sum(present_weights, na.rm = TRUE)

# Probability-style weighted mean in [0,1]
weighted_probability <- weighted_sum / weight_total_present

# Where no rasters had data, result should be NA (division produced NaN/Inf if denom=0)
weighted_probability[weight_total_present == 0] <- NA
weighted_probability <- clamp(weighted_probability, 0, 1)

#-------------------------------
# 5) Optional: minimum contributing rasters requirement
#-------------------------------
num_sources <- sum(present, na.rm = TRUE)
min_sources <- 2  # change as needed
weighted_probability_minN <- ifel(num_sources >= min_sources, weighted_probability, NA)

#-------------------------------
# 6) Binary output and exports
#-------------------------------
thr <- 0.5  # adjust or tune via ROC if have validation data
binary_map <- weighted_probability >= thr

writeRaster(
  weighted_probability,
  file.path(out_dir, "weighted_grassland_probability.tif"),
  overwrite = TRUE
)

writeRaster(
  weighted_probability_minN,
  file.path(out_dir, "weighted_grassland_probability_minSources.tif"),
  overwrite = TRUE
)

writeRaster(
  binary_map,
  file.path(out_dir, "weighted_grassland_binary_0p5.tif"),
  overwrite = TRUE
)

# Optional diagnostics
writeRaster(
  num_sources,
  file.path(out_dir, "num_contributing_rasters.tif"),
  overwrite = TRUE
)

writeRaster(
  weight_total_present,
  file.path(out_dir, "sum_weights_present.tif"),
  overwrite = TRUE
)

## MORE IDEAS
# 1. If some accuracies are very low (e.g., <0.5), could downweight them further or drop them:
#    weight_floor <- 0.5
#    present_weights <- ifel(weight_stack >= weight_floor, weight_stack, 0) * present
# 2. Track uncertainty: Export num_sources and sum_weights_present to visualize where 
#    probabilities are supported by strong/weak evidence.
# 3. Use minimum total weight instead of minimum count:
#    min_total_weight <- 1.5
#    weighted_probability_minW <- ifel(weight_total_present >= min_total_weight, 
#                                 weighted_probability, NA)


#####################################################################
#####################################################################

## ADVANCED: two versions that use recall/sensitivity/TPR and give uncertainty

# -----------------------------
# 1) Inputs
# -----------------------------

# Accuracy table with recall (sensitivity) per raster per ecodistrict:
# Columns required: ecodistrict_id, raster_name, recall
acc_bin <- read.csv("path/to/recall_table.csv", stringsAsFactors = FALSE)

# Ensure raster_name matches layer names
# acc$raster_name <- tools::file_path_sans_ext(basename(acc$raster_name))
stopifnot(all(acc_bin$raster_name %in% names(r_mask)))

# -----------------------------
# 2) Prepare per-raster recall rasters
# -----------------------------
acc_bin_wide <- acc_bin %>%
  pivot_wider(names_from = raster_name, values_from = recall)

if (!"ecodistrict_id" %in% names(ecodist)) {
  stop("Field 'ecodistrict_id' is missing in ecodistricts.")
}

ecodist_bin <- merge(ecodist, acc_bin_wide, by = "ecodistrict_id", duplicateGeoms = TRUE)

# Rasterize recall per raster layer
recall_list <- vector("list", nlyr(r_mask))
names(recall_list) <- names(r_mask)
for (rname in names(r_mask)) {
  recall_list[[rname]] <- rasterize(ecodist_bin, r_mask[[1]], field = rname, touches = TRUE)
}
recall_stack <- rast(recall_list)
recall_stack <- recall_stack[[names(r_mask)]]  # ensure same order

# Presence (not NA in data rasters)
present <- !is.na(r_mask)

# -----------------------------
# 3A) Strategy A: recall-weighted positive evidence
# -----------------------------
# Use recall to weight only the “1” (grassland) predictions.
# Ignore NA completely.
# Let “0” (known not grassland) count as a vote against, but unweighted (or with a neutral weight of 1).
# Interpretation: a probability‑like index that reflects strength of positive evidence for grassland 
# based on how well each raster detects grassland (recall), without claiming negative reliability we don’t know.
# Pros: No extra assumptions.
# Cons: Under-penalizes false positives if some layers have low specificity.

# Weight only positive predictions (1) by recall; zeros are negative but unweighted.
# Numerator = sum( (r==1) * recall )
# Denominator = sum( recall ) over rasters that had data (0 or 1)
pos_mask <- r_mask == 1
neg_mask <- r_mask == 0

# Zero-out recall where raster has NA
recall_present <- recall_stack * present

numerator_A <- sum(ifel(pos_mask, recall_present, 0), na.rm = TRUE)
denominator_A <- sum(recall_present, na.rm = TRUE)
prob_A <- numerator_A / denominator_A
prob_A[denominator_A == 0] <- NA
prob_A <- clamp(prob_A, 0, 1)

# Diagnostics for A
num_contributing <- sum(present, na.rm = TRUE)
num_positive <- sum(pos_mask, na.rm = TRUE)
num_negative <- sum(neg_mask, na.rm = TRUE)
sum_recall_present <- denominator_A

# Interpretation:
# Interpreted as recall‑weighted positive evidence (probability‑like index).
# High when many positives with high recall and few/no negatives, regardless of specificity.
# Use when avoiding assumptions about negative reliability.

# A simple uncertainty metric: high when prob~0.5 or when little recall is present
# uncertainty_probOnly is highest near 0.5.
#The composite adds a data sufficiency factor (sum_recall_present), so areas with few/weak 
# inputs are flagged as more uncertain even if the current probability isn’t near 0.5.
uncertainty_prob_A <- prob_A * (1 - prob_A)
uncertainty_data_A <- 1 / (1 + sum_recall_present)  # decreases as recall evidence accumulates
uncertainty_composite_A <- uncertainty_prob_A * (1 + uncertainty_data_A)

# -----------------------------
# 3B) Strategy B: Bayesian fusion with assumed specificity
# -----------------------------
# Treat each raster as a detector with known sensitivity = recall and assumed specificity,
# either a single global constant (e.g., 0.85), or a per-ecodistrict constant if known.
# Use likelihood ratios:
#    If raster predicts 1: multiply odds by LR+ = sens / (1 − spec)
#    If raster predicts 0: multiply odds by LR− = (1 − sens) / spec
# Combine across rasters in log‑odds to avoid numeric issues.
# Yields a proper posterior probability of grassland given inputs and assumed specificity.
# Pros: Proper probability and negatives matter; interpretable uncertainty via information content.
# Cons: Requires choosing specificity (can make it a parameter and sensitivity‑test it).

# Choose an assumed specificity (global). Can sensitivity-test this.
spec_assumed <- 0.85  # change as needed (0.5-0.95 reasonable if unknown)

# Clip recall and specificity to avoid log(0)
eps <- 1e-6
sens <- clamp(recall_present, eps, 1 - eps)   # sensitivity = recall
spec <- spec_assumed
spec <- max(min(spec, 1 - eps), eps)

# Likelihood Ratios
# For Xi = 1: LR+ = sens / (1 - spec)
# For Xi = 0: LR- = (1 - sens) / spec
LRpos <- sens / (1 - spec)
LRneg <- (1 - sens) / spec

# Log-likelihood contributions per raster/pixel
logLR_pos <- log(LRpos)
logLR_neg <- log(LRneg)

# Build per-layer logLR chosen by observed value (1 or 0); NA contributes 0
logLR_stack <- ifel(pos_mask, logLR_pos,
                    ifel(neg_mask, logLR_neg, 0))

# Sum log LR across rasters
logLR_sum <- sum(logLR_stack, na.rm = TRUE)

# Prior: choose a prior grassland probability (global or ecodistrict)
# Use a neutral prior if preferred: 0.5 => logit = 0
prior_p <- 0.5
logit <- function(p) log(p / (1 - p))
inv_logit <- function(x) 1 / (1 + exp(-x))

log_odds_prior <- logit(prior_p)
log_odds_post <- log_odds_prior + logLR_sum
prob_B <- inv_logit(log_odds_post)

# or replace with rasterized ecodistrict field
#ecodist$prior_grass <- <numeric values 0–1>
#prior_p_r <- rasterize(ecodist, r_stack[[1]], field="prior_grass")
#log_odds_prior <- log(prior_p_r / (1 - prior_p_r))
#log_odds_post <- log_odds_prior + logLR_sum
#prob_B <- inv_logit(log_odds_post)

# Interpretation
# A Bayesian posterior probability under the assumed specificity and recall from validation data.
# Positives push probability up by LR+; negatives push down by LR−.
# Tune spec_assumed and prior_p and re‑export to see sensitivity.
# information_sum_absLogLR shows where evidence is strong (large) vs. weak (small).

# Uncertainty diagnostics for B
# Information magnitude (larger |logLR| => more decisive information)
# uncertainty_probOnly is highest near 0.5.
#The composite adds a data sufficiency factor (info_sum), so areas with few/weak inputs are flagged 
# as more uncertain even if the current probability isn’t near 0.5.
info_sum <- sum(abs(logLR_stack), na.rm = TRUE)
uncertainty_prob_B <- prob_B * (1 - prob_B)
uncertainty_info_B <- 1 / (1 + info_sum)  # decreases as information grows
uncertainty_composite_B <- uncertainty_prob_B * (1 + uncertainty_info_B)

# -----------------------------
# 4) Optional masks (minimum support)
# -----------------------------
min_sources <- 2
prob_A_minN <- ifel(num_contributing >= min_sources, prob_A, NA)
prob_B_minN <- ifel(num_contributing >= min_sources, prob_B, NA)

# -----------------------------
# 5) Write outputs
# -----------------------------
writeRaster(prob_A,          "prob_grassland_strategyA_recallWeighted.tif", overwrite = TRUE)
writeRaster(prob_A_minN,     "prob_grassland_strategyA_minN.tif", overwrite = TRUE)
writeRaster(uncertainty_composite_A, "uncertainty_strategyA_composite.tif", overwrite = TRUE)
writeRaster(uncertainty_prob_A,      "uncertainty_strategyA_probOnly.tif", overwrite = TRUE)
writeRaster(sum_recall_present,      "sum_recall_present.tif", overwrite = TRUE)

writeRaster(prob_B,          "prob_grassland_strategyB_bayesian.tif", overwrite = TRUE)
writeRaster(prob_B_minN,     "prob_grassland_strategyB_minN.tif", overwrite = TRUE)
writeRaster(uncertainty_composite_B, "uncertainty_strategyB_composite.tif", overwrite = TRUE)
writeRaster(uncertainty_prob_B,      "uncertainty_strategyB_probOnly.tif", overwrite = TRUE)
writeRaster(info_sum,               "information_sum_absLogLR.tif", overwrite = TRUE)

# Useful diagnostics to map
writeRaster(num_contributing, "num_contributing_rasters.tif", overwrite = TRUE)
writeRaster(num_positive,     "num_positive_votes.tif", overwrite = TRUE)
writeRaster(num_negative,     "num_negative_votes.tif", overwrite = TRUE)

## Ways to improve
# Per‑ecodistrict priors: Get grassland prevalence per ecodistrict and rasterize it to a prior and swap prior_p.
# Per‑raster specificity later: Get TN/FP and extend the table to include specificity; can rasterize it exactly like recall and use per‑layer spec.
# Calibrate a threshold: Use an external validation set and get ROC to pick an optimal threshold for the Bayesian map.

# ============================
# 7) TRI-CLASS SYNTHESIS (0=not, 1=tame, 2=native)
# ============================

# -----------------------------
# 7.1 Identify or load 3-class rasters (0=not, 1=tame, 2=native)
# -----------------------------
# Explicitly list the 3-class rasters
nt_files <- c(
  aci_nt_file, ngi_nt_file, cgr_nt_file, eccc_nt_file,
  cnvc_nt_file, pudl_nt_file, pli_ap_nt_file, pli_cu_nt_file,
  pli_mg_nt_file, pli_mm_nt_file, gpw_nt_file
)

nt_stack <- NULL
if (length(nt_files) > 0) {
  nt_stack <- rast(nt_files)
  # Names: make consistent with accuracy_bin_df$raster_name (lowercase + prefix rule)
  make_short_name <- function(path) {
    bn <- basename(path)
    stem <- tools::file_path_sans_ext(bn)
    if (grepl("_", stem)) sub("^([^_]+)_.*$", "\\1", stem) else stem
  }
  tri_names <- tolower(vapply(nt_files, make_short_name, FUN.VALUE = character(1)))
  tri_names <- make.unique(tri_names, sep = "_")
  names(nt_stack) <- tri_names
  
  # Align to r_stack geometry if needed
  if (!compareGeom(r_mask[[1]], nt_stack[[1]], stopOnError = FALSE)) {
    nt_stack <- project(nt_stack, r_mask[[1]])  # or resample() if same CRS but different grid
  }
} else {
  message("No tri-class rasters provided. The code will fall back to a prior for native/tame share.")
}

# -----------------------------------------------------------
# 7.2 GENERALIZED WEIGHT BUILDER FOR ANY ACCURACY METRIC
# -----------------------------------------------------------
# The CSV may have either:
#  - 'recall' only (grass detection)
#  - OR 'recall_native' and 'recall_tame' per ecodistrict (preferred)
# Will detect and use what's available.

# Choose accuracy metrics here (defaults to recall)
metric_native <- "recall_native" # could be recall, precision, F1, specificity, etc.
metric_tame   <- "recall_tame"
metric_single <- "recall"        # used if only one metric exists - could be MCC, OA, etc.

accuracy_nt_df <- readr::read_csv(acc_nt_file, show_col_types = FALSE)
accuracy_nt_df <- accuracy_nt_df %>%
  mutate(raster_name = tolower(stringr::str_trim(raster_name)))

# Basic checks
if (!all(c("raster_name", "ecodistrict_id") %in% names(accuracy_nt_df))) {
  stop("Tri-class accuracy CSV must contain 'raster_name' and 'ecodistrict_id'.")
}

# Detect metric availability
has_native_tame <- all(c(metric_native, metric_tame) %in% names(accuracy_nt_df))
has_single_metric <- metric_single %in% names(accuracy_nt_df)

# If no tri-class rasters, set flags false
if (is.null(nt_stack)) {
  has_native_tame <- FALSE
  has_single_metric <- FALSE
}

metric_native_stack <- NULL
metric_tame_stack   <- NULL
metric_single_stack <- NULL

if (!is.null(nt_stack)) {
  
  # Filter accuracy rows matching tri-class rasters
  acc_tri <- accuracy_nt_df %>%
    dplyr::filter(raster_name %in% names(nt_stack))
  
  if (nrow(acc_tri) == 0) {
    warning("No accuracy rows match tri-class rasters - using weight = 1.")
    metric_single_stack <- rast(lapply(1:nlyr(nt_stack), function(i) nt_stack[[i]]*0 + 1))
    names(metric_single_stack) <- names(nt_stack)
    
  } else if (has_native_tame) {
    
    # --- WIDE tables for two metrics ---
    acc_native_wide <- acc_tri %>%
      select(ecodistrict_id, raster_name, !!metric_native) %>%
      pivot_wider(names_from = raster_name, values_from = !!metric_native)
    
    acc_tame_wide <- acc_tri %>%
      select(ecodistrict_id, raster_name, !!metric_tame) %>%
      pivot_wider(names_from = raster_name, values_from = !!metric_tame)
    
    # Merge into polygons
    eco_native <- merge(ecodist, acc_native_wide, by = "ecodistrict_id", duplicateGeoms = TRUE)
    eco_tame   <- merge(ecodist, acc_tame_wide,   by = "ecodistrict_id", duplicateGeoms = TRUE)
    
    # Rasterize per raster layer
    metric_native_stack <- rast(lapply(names(nt_stack), function(rn) {
      rasterize(eco_native, r_mask[[1]], field = rn, touches = TRUE)
    }))
    names(metric_native_stack) <- names(nt_stack)
    
    metric_tame_stack <- rast(lapply(names(nt_stack), function(rn) {
      rasterize(eco_tame, r_mask[[1]], field = rn, touches = TRUE)
    }))
    names(metric_tame_stack) <- names(nt_stack)
    
  } else if (has_single_metric) {
    
    # --- Single metric provided for all positive predictions ---
    acc_single_wide <- acc_tri %>%
      select(ecodistrict_id, raster_name, !!metric_single) %>%
      pivot_wider(names_from = raster_name, values_from = !!metric_single)
    
    eco_single <- merge(ecodist, acc_single_wide, by = "ecodistrict_id", duplicateGeoms = TRUE)
    
    metric_single_stack <- rast(lapply(names(nt_stack), function(rn) {
      rasterize(eco_single, r_mask[[1]], field = rn, touches = TRUE)
    }))
    names(metric_single_stack) <- names(nt_stack)
    
  } else {
    
    warning("No usable metric columns found - using weight = 1.")
    metric_single_stack <- rast(lapply(1:nlyr(nt_stack), function(i) nt_stack[[i]]*0 + 1))
    names(metric_single_stack) <- names(nt_stack)
  }
}

# ============================================================
# ALTERNATE 7.2 GENERALIZED MCC + CLASS-SPECIFIC METRIC WEIGHT BUILDER
# ============================================================

# Choose metric names here:
metric_mcc      <- "MCC_1v2"
metric_f1_native <- "F1_2_2x2"
metric_f1_tame   <- "F1_1_2x2"

# Weight exponents (tune if needed)
alpha_mcc <- 1.0 # (increase to 1.5 to prioritize MCC more)
beta_f1   <- 1.0

accuracy_nt_df <- readr::read_csv(acc_nt_file, show_col_types = FALSE) %>%
  mutate(raster_name = tolower(stringr::str_trim(raster_name)))

required_cols <- c("raster_name", "ecodistrict_id", 
                   metric_mcc, metric_f1_native, metric_f1_tame)

if (!all(required_cols %in% names(accuracy_nt_df))) {
  stop("Accuracy CSV must contain: ", paste(required_cols, collapse = ", "))
}

# ---- Build composite weights per raster x ecodistrict ----
accuracy_nt_df <- accuracy_nt_df %>%
  mutate(
    W_native = ( !!sym(metric_mcc)      )^alpha_mcc *
      ( !!sym(metric_f1_native) )^beta_f1,
    
    W_tame   = ( !!sym(metric_mcc)      )^alpha_mcc *
      ( !!sym(metric_f1_tame)   )^beta_f1
  )

# Replace NA or negative values
accuracy_nt_df$W_native[!is.finite(accuracy_nt_df$W_native)] <- NA
accuracy_nt_df$W_tame[!is.finite(accuracy_nt_df$W_tame)]     <- NA

accuracy_nt_df$W_native <- pmax(accuracy_nt_df$W_native, 0)
accuracy_nt_df$W_tame   <- pmax(accuracy_nt_df$W_tame, 0)

# ---- Rasterize weights by ecodistrict ----
acc_native_wide <- accuracy_nt_df %>%
  select(ecodistrict_id, raster_name, W_native) %>%
  tidyr::pivot_wider(names_from = raster_name, values_from = W_native)

acc_tame_wide <- accuracy_nt_df %>%
  select(ecodistrict_id, raster_name, W_tame) %>%
  tidyr::pivot_wider(names_from = raster_name, values_from = W_tame)

eco_native <- merge(ecodist, acc_native_wide, by = "ecodistrict_id", duplicateGeoms = TRUE)
eco_tame   <- merge(ecodist, acc_tame_wide,   by = "ecodistrict_id", duplicateGeoms = TRUE)

metric_native_stack <- rast(lapply(names(nt_stack), function(rn) {
  rasterize(eco_native, r_mask[[1]], field = rn, touches = TRUE)
}))
metric_tame_stack <- rast(lapply(names(nt_stack), function(rn) {
  rasterize(eco_tame, r_mask[[1]], field = rn, touches = TRUE)
}))

names(metric_native_stack) <- names(nt_stack)
names(metric_tame_stack)   <- names(nt_stack)

# -----------------------------
# 7.3 Compute share_native / share_tame from tri-class evidence
# -----------------------------
# Masks on tri-class predictions
share_native <- NULL
share_tame   <- NULL
min_tri_sources <- 1

if (!is.null(nt_stack)) {
  tri_present <- !is.na(nt_stack)
  native_mask <- nt_stack == 2
  tame_mask   <- nt_stack == 1
  
  if (has_native_tame) {
    metric_native_present <- metric_native_stack * tri_present
    metric_tame_present   <- metric_tame_stack   * tri_present
    
    native_votes <- sum(ifel(native_mask, metric_native_present, 0), na.rm = TRUE)
    tame_votes   <- sum(ifel(tame_mask,   metric_tame_present,   0), na.rm = TRUE)
    
  } else {
    metric_present <- metric_single_stack * tri_present
    native_votes <- sum(ifel(native_mask, metric_present, 0), na.rm = TRUE)
    tame_votes   <- sum(ifel(tame_mask,   metric_present, 0), na.rm = TRUE)
  }
  
  vote_total <- native_votes + tame_votes
  tri_num_sources <- sum(tri_present, na.rm = TRUE)
  
  share_native_raw <- native_votes / vote_total
  share_native_raw[vote_total == 0] <- NA
  share_native_raw <- clamp(share_native_raw, 0, 1)
  
  share_native <- ifel(tri_num_sources >= min_tri_sources, share_native_raw, NA)
  share_tame   <- 1 - share_native
}

# ============================================================
# ALTERNATE 7.3 Use MCC+F1 composite weights to derive native/tame shares
# ============================================================

tri_present <- !is.na(nt_stack)
native_mask <- nt_stack == 2
tame_mask   <- nt_stack == 1

metric_native_present <- metric_native_stack * tri_present
metric_tame_present   <- metric_tame_stack   * tri_present

native_votes <- sum(ifel(native_mask, metric_native_present, 0), na.rm=TRUE)
tame_votes   <- sum(ifel(tame_mask,   metric_tame_present,   0), na.rm=TRUE)

vote_total <- native_votes + tame_votes
share_native_raw <- native_votes / vote_total
share_native_raw[vote_total == 0] <- NA
share_native_raw <- clamp(share_native_raw, 0, 1)

tri_num_sources <- sum(tri_present, na.rm=TRUE)
min_tri_sources <- 1

share_native <- ifel(tri_num_sources >= min_tri_sources, share_native_raw, NA)
share_tame   <- 1 - share_native

# -----------------------------
# 7.4 Define a fallback prior for native/tame split where tri evidence is missing
# -----------------------------
# Use an ecodistrict prior for native fraction (e.g., from CNVC/NatureServe) by
# putting it in the polygons as field 'prior_native' and script will rasterize it.
# Otherwise, use a global neutral prior (0.5) or a chosen constant per-ecodistrict.

prior_native_global <- 0.5  # adjust as appropriate

if (!"prior_native" %in% names(ecodist)) {
  ecodist$prior_native <- prior_native_global
}

prior_native_r <- rasterize(ecodist, r_mask[[1]], field = "prior_native", touches = TRUE)
prior_native_r <- clamp(prior_native_r, 0, 1)

# Fill gaps in share_native/share_tame with the prior
if (is.null(share_native)) {
  share_native <- prior_native_r
  share_tame   <- 1 - share_native
} else {
  share_native <- ifel(is.na(share_native), prior_native_r, share_native)
  share_tame   <- 1 - share_native
}

# -----------------------------
# 7.5 Compose final class probabilities
# -----------------------------
# Use chosen P(grass): prob_B (Bayesian) or prob_A
p_grass <- prob_A  # or prob_A
p_grass <- clamp(p_grass, 0, 1)

p_native   <- p_grass * share_native
p_tame     <- p_grass * share_tame
p_nongrass <- 1 - p_grass

# -----------------------------
# 7.6 Hard 3-class map (0=not, 1=tame, 2=native)
# -----------------------------
p_stack <- c(p_nongrass, p_tame, p_native)
names(p_stack) <- c("p_not", "p_tame", "p_native")

idx_max <- which.max(p_stack)  # 1,2,3
# Map to 0,1,2
rcl <- matrix(c(1,0,  2,1,  3,2), ncol = 2, byrow = TRUE)
class_3cat <- classify(idx_max, rcl = rcl)
names(class_3cat) <- "grass_3class"  # 0/1/2

# -----------------------------
# 7.7 Uncertainty diagnostics
# -----------------------------
max_p <- max(p_stack)  # p_stack has layers: p_not, p_tame, p_native
uncertainty_1minus <- 1 - max_p

# Normalized entropy (0=certain, 1=max uncertainty)
eps <- 1e-9
H <- - (p_nongrass * log(p_nongrass + eps) +
          p_tame     * log(p_tame     + eps) +
          p_native   * log(p_native   + eps)) / log(3)
uncertainty_entropy <- clamp(H, 0, 1)

# -----------------------------
# 7.8 Optional masks by support
# -----------------------------
# Already have num_sources from earlier (number of rasters contributing to grass/not)
# For tri share support:
if (!is.null(nt_stack)) {
  tri_num_sources <- sum(!is.na(nt_stack), na.rm = TRUE)
  writeRaster(tri_num_sources, file.path(out_dir, "tri_num_contributing.tif"), overwrite = TRUE)
}

# -----------------------------
# 7.9 Write outputs
# -----------------------------
writeRaster(p_nongrass, file.path(out_dir, "p_not_grass.tif"),    overwrite = TRUE)
writeRaster(p_tame,     file.path(out_dir, "p_tame.tif"),         overwrite = TRUE)
writeRaster(p_native,   file.path(out_dir, "p_native.tif"),       overwrite = TRUE)
writeRaster(class_3cat, file.path(out_dir, "grass_3class_A.tif"),   overwrite = TRUE)

writeRaster(uncertainty_1minus, file.path(out_dir, "uncertainty_1minus_maxP.tif"), overwrite = TRUE)
writeRaster(uncertainty_entropy, file.path(out_dir, "uncertainty_entropy_norm.tif"), overwrite = TRUE)

# ============================
# 8) DIAGNOSTICS
# ============================

# Basic summaries
# If p_grass is mostly near 0, the issue is binary fusion 
# If share_tame is mostly near 0 (and share_native near 1), the issue is tri-class weighting / coding
print(global(p_grass,     fun=function(x,...) c(min=min(x,na.rm=TRUE), mean=mean(x,na.rm=TRUE), max=max(x,na.rm=TRUE)), na.rm=TRUE))
print(global(share_tame,  fun=function(x,...) c(min=min(x,na.rm=TRUE), mean=mean(x,na.rm=TRUE), max=max(x,na.rm=TRUE)), na.rm=TRUE))
print(global(share_native,fun=function(x,...) c(min=min(x,na.rm=TRUE), mean=mean(x,na.rm=TRUE), max=max(x,na.rm=TRUE)), na.rm=TRUE))

# Are equality masks working? Count how many tame/native cells exist in tri-class predictions
# If either is ~0, likely still have fractional classes. Re-do project(..., method="near") and rounding.
# 1) Equality masks (make sure you reprojected tri-class with method="near" and rounded)
tame_mask   <- nt_stack == 1
native_mask <- nt_stack == 2

# 2) Count cells where masks are TRUE (numbers should be > 0 if equality works)
tame_count   <- as.numeric(global(tame_mask,   fun = "sum", na.rm = TRUE)[1,1])
native_count <- as.numeric(global(native_mask, fun = "sum", na.rm = TRUE)[1,1])

cat("tame_mask TRUE cells:   ", tame_count,   "\n")
cat("native_mask TRUE cells: ", native_count, "\n")

# 3) Also check total tri coverage
present_tri <- !is.na(nt_stack)
present_count <- as.numeric(global(present_tri, fun = "sum", na.rm = TRUE)[1,1])
cat("tri-class present cells:", present_count, "\n")

# After constructing recall_*_present and masks
# Note: adapt to branch (has_native_tame vs single recall)
# Weighted vote totals:
# native_votes and tame_votes should already be in env; if not, recompute quickly:
tri_present <- !is.na(nt_stack)
native_mask <- nt_stack == 2
tame_mask   <- nt_stack == 1

if (exists("metric_native_stack") && exists("metric_tame_stack") &&
    !is.null(metric_native_stack) && !is.null(metric_tame_stack)) {
  metric_native_present <- metric_native_stack * tri_present
  metric_tame_present   <- metric_tame_stack   * tri_present
  native_votes <- sum(ifel(native_mask, metric_native_present, 0), na.rm = TRUE)
  tame_votes   <- sum(ifel(tame_mask,   metric_tame_present,   0), na.rm = TRUE)
} else {
  metric_present_tri <- if (exists("metric_nt_stack") && !is.null(metric_nt_stack)) (metric_nt_stack * tri_present) else (tri_present*0 + 1)
  native_votes <- sum(ifel(native_mask, metric_present_tri, 0), na.rm = TRUE)
  tame_votes   <- sum(ifel(tame_mask,   metric_present_tri, 0), na.rm = TRUE)
}

# Where is tame getting outvoted?
vote_total <- native_votes + tame_votes
share_native_raw <- native_votes / vote_total; share_native_raw[vote_total == 0] <- NA
share_native_raw <- clamp(share_native_raw, 0, 1)
plot(share_native_raw, main = "Share native") # quick visual: should reveal native dominance spatially

# Check prior influence: If tri evidence is sparse, share_native fills from prior_native_r.
# Map the fraction of cells using the prior (i.e., where tri votes were NA)
prior_used <- is.na(share_native_raw)  # TRUE where we fell back to prior
plot(prior_used, main = "Fell back to prior?")                       # white/TRUE means driven by the prior

# Validate coding per layer: Pick a layer and examine frequencies
i <- 1 # or which(names(nt_stack)=="aci") etc.
# For a single layer i
terra::freq(nt_stack[[i]])  # returns a matrix of unique values and counts

# Check whether weights (recalls) are skewed
# Using the tri-class accuracy table you loaded: accuracy_nt_df
# (Make sure names are normalized to lower case like earlier)
if (all(c("metric_native", "metric_tame") %in% names(accuracy_nt_df))) {
  # overall summary
  print(accuracy_nt_df %>% 
          summarise(
            mean_native = mean(metric_native, na.rm=TRUE),
            mean_tame   = mean(metric_tame,   na.rm=TRUE),
            q25_native  = quantile(metric_native, 0.25, na.rm=TRUE),
            q50_native  = quantile(metric_native, 0.50, na.rm=TRUE),
            q75_native  = quantile(metric_native, 0.75, na.rm=TRUE),
            q25_tame    = quantile(metric_tame,   0.25, na.rm=TRUE),
            q50_tame    = quantile(metric_tame,   0.50, na.rm=TRUE),
            q75_tame    = quantile(metric_tame,   0.75, na.rm=TRUE)
          ))
}

### Current tri-class stack looks like it has much more native evidence than tame, 
### so share_tame = small, and therefore p_tame = p_grass * share_tame = small.

# Unweighted per-pixel counts across tri-class layers
tri_present <- !is.na(nt_stack)
n_native <- sum(nt_stack == 2, na.rm = TRUE)
n_tame   <- sum(nt_stack == 1, na.rm = TRUE)

# Quick summaries
print(global(n_native, fun=function(x,...) c(min=min(x,na.rm=TRUE), mean=mean(x,na.rm=TRUE), max=max(x,na.rm=TRUE)), na.rm=TRUE))
print(global(n_tame,   fun=function(x,...) c(min=min(x,na.rm=TRUE), mean=mean(x,na.rm=TRUE), max=max(x,na.rm=TRUE)), na.rm=TRUE))

# (Optional) visualize where native votes exceed tame votes
vote_margin <- n_native - n_tame
plot(vote_margin, main = "Vote margin (native exceeds tame)")  # bigger positive => stronger native majority

# Confirm it is vote-driven (not weighting)
# Mask to pilot area first:
vote_margin_in <- mask(vote_margin, pa)
print(global(vote_margin_in, fun=function(x,...) c(min=min(x,na.rm=TRUE),
                                                   mean=mean(x,na.rm=TRUE),
                                                   max=max(x,na.rm=TRUE)), na.rm=TRUE))

# Check per-layer distributions
# If most layers have far more 2s than 1s, and per-pixel n_native >= n_tame almost everywhere, 
# then the stack composition is simply native-heavy. Plausible if several inputs represent natural 
# vegetation products while only a subset map tame (agri/pasture) specifically.
layer_counts <- function(s) {
  tibble::tibble(
    layer = names(s),
    n0 = sapply(1:nlyr(s), function(i) as.numeric(global(s[[i]] == 0, fun="sum", na.rm=TRUE)[1,1])),
    n1 = sapply(1:nlyr(s), function(i) as.numeric(global(s[[i]] == 1, fun="sum", na.rm=TRUE)[1,1])),
    n2 = sapply(1:nlyr(s), function(i) as.numeric(global(s[[i]] == 2, fun="sum", na.rm=TRUE)[1,1]))
  )
}
print(layer_counts(nt_stack))


