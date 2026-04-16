
###################################################################################
##                     GRASSLAND CONSENSUS MAPPING ANALYSIS                      ##
###################################################################################

##----------------------- Alison Long - April 2026 -------------------------

# step 8: use data fusion techniques to create a consensus grassland map 
#         to serve as baseline for modelling

# load packages
library(terra)
library(tidyverse)

# Set the global terra options
terraOptions(
  memfrac = 0.6,
  threads = max(1, parallel::detectCores() - 1),
  tempdir = "R:/CanadianGrasslands/WG_classification/scratch"
)

# write options
wopt_prob  <- list(datatype="FLT4S", gdal=c("COMPRESS=LZW","TILED=YES","BIGTIFF=YES"))
wopt_byte  <- list(datatype="INT1U", gdal=c("COMPRESS=LZW","TILED=YES","BIGTIFF=YES"))
wopt_float <- list(datatype="FLT4S", gdal=c("COMPRESS=LZW","TILED=YES","BIGTIFF=YES"))
wopt_count <- list(datatype="INT2U", gdal=c("COMPRESS=LZW","TILED=YES","BIGTIFF=YES"))

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

bin_files <- c(
  aci_bin_file, luts30_bin_file, luts10_bin_file, cec_bin_file, ngi_bin_file, cgr_bin_file, 
  eccc_bin_file, esa_bin_file, glc_bin_file, cnvc_bin_file, lcc_bin_file, pudl_bin_file, 
  pli_ap_bin_file, pli_cu_bin_file, pli_mg_bin_file, pli_mm_bin_file, lcr_bin_file, 
  gpw_bin_file, plow_bin_file, plvi_bin_file, hfi_bin_file, gcc_bin_file, dw_bin_file,
  esri_bin_file
)

#-------------------------------
# Mask binary rasters (streamed to disk)
#-------------------------------

mask_file <- file.path(out_dir, "binary_stack_masked.tif")

if (!file.exists(mask_file)) {
  
  sa <- rast(sa_file)
  
  r_mask <- mask(
    crop(rast(bin_files), sa),
    sa,
    filename  = mask_file,
    overwrite = TRUE,
    wopt      = wopt_byte
  )
  
  #rm(sa, r_mask)
  gc()
}

r_mask <- rast(mask_file)

# Assign short layer names
make_short_name <- function(path) {
  stem <- tools::file_path_sans_ext(basename(path))
  
  if (grepl("^sk_pli_", stem)) {
    sub("^sk_(pli_[^_]+).*", "\\1", stem)
  } else {
    sub("^([^_]+)_.*$", "\\1", stem)
  }
}

names(r_mask) <- tolower(vapply(bin_files, make_short_name, character(1)))

# read polygons and accuracy table
ecodist <- vect(ecodist_file)
ecodist <- project(ecodist, crs(r_mask))
ecodist$ecodistrict_id <- ecodist$ECODISTRIC
ecodist$ecoregion_id <- ecodist$ECOREGION_

# Minimum sample sizes for pooling
MIN_N_ECODIST  <- 5     # per ecodistrict, per raster
MIN_N_ECOREG   <- 20    # pooled across ecoregion, per raster

accuracy_bin_df <- read_csv(acc_bin_file, show_col_types = FALSE)

# Harmonize raster names by updating accuracy_bin_df names
accuracy_bin_df <- accuracy_bin_df %>%
  mutate(
    raster_name = case_when(
      raster_name == "aci"   ~ "aafc",
      raster_name == "cec" ~ "nalcms",
      raster_name == "ngi" ~ "cfga",
      TRUE                    ~ raster_name
    )
  )

# Preview
cat("Assigned r_stack names:\n")
print(names(r_mask))

# Print unique raster names
unique_rasters <- unique(accuracy_bin_df$raster_name)
print(unique_rasters)

stopifnot(all(accuracy_bin_df$raster_name %in% names(r_mask)))

# Join ecodistrict-level accuracy values to polygons
acc_bin_ecodist <- accuracy_bin_df %>%
  select(ecodistrict_id, raster_name, recall) %>%
  pivot_wider(names_from = raster_name, values_from = recall)

# Fallback 1: Use ecoregion-level accuracy for ecodistricts without accuracy values
acc_bin_ecoreg <- accuracy_bin_df %>%
  left_join(
    ecodist %>%
      as.data.frame() %>%
      select(ecodistrict_id, ecoregion_id),
    by = "ecodistrict_id"
  ) %>%
  group_by(ecoregion_id, raster_name) %>%
  summarise(
    recall_ecoregion = mean(recall, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from  = raster_name,
    values_from = recall_ecoregion
  )

# Fallback 2: Use global accuracy for ecodistricts without ecodist or ecoreg-level accuracy
acc_bin_global <- accuracy_bin_df %>%
  group_by(raster_name) %>%
  summarise(
    recall_global = mean(recall, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from  = raster_name,
    values_from = recall_global
  )

# Join recall to polygons
#ecodist <- merge(
#  ecodist, acc_bin_ecodist,
#  by.x = "ecodistrict_id",
#  by.y = "ecodistrict_id",
#  duplicateGeoms = TRUE
#)

# Build final hierarchical recall table per ecodistrict
# Base table: one row per ecodistrict
ecodist_weights <- ecodist %>%
  as.data.frame() %>%
  select(ecodistrict_id, ecoregion_id)

# Join ecodistrict recall
ecodist_weights <- ecodist_weights %>%
  left_join(acc_bin_ecodist, by = "ecodistrict_id")

# Join ecoregion recall (suffix _eco)
ecodist_weights <- ecodist_weights %>%
  left_join(
    acc_bin_ecoreg,
    by = "ecoregion_id",
    suffix = c("", "_eco")
  )

# Add global recall (single-row replicated)
ecodist_weights <- bind_cols(
  ecodist_weights,
  acc_bin_global[rep(1, nrow(ecodist_weights)), ]
)

# Apply fallback logic raster‑by‑raster
raster_names <- names(r_mask)

for (nm in raster_names) {
  
  ecodist_col <- nm
  eco_col     <- paste0(nm, "_eco")
  global_col  <- nm
  
  ecodist_weights[[nm]] <- dplyr::coalesce(
    ecodist_weights[[ecodist_col]],  # ecodistrict
    ecodist_weights[[eco_col]],      # ecoregion
    ecodist_weights[[global_col]]    # global
  )
}

# Convert back to SpatVector and rasterize weights
# Convert to SpatVector (preserve geometry)
ecodist_weights_vect <- vect(
  ecodist_weights,
  geom = geom(ecodist),
  crs  = crs(ecodist)
)

# Rasterize hierarchical weights
weight_dir <- file.path(out_dir, "weights_hierarchical")
dir.create(weight_dir, showWarnings = FALSE)

weight_files <- setNames(character(length(raster_names)), raster_names)
template <- r_mask[[1]]

for (nm in raster_names) {
  
  out_w <- file.path(weight_dir, paste0("weight_", nm, ".tif"))
  
  if (!file.exists(out_w)) {
    message("Rasterizing hierarchical weight for ", nm)
    
    rasterize(
      ecodist_weights_vect,
      template,
      field     = nm,
      touches   = TRUE,
      filename  = out_w,
      overwrite = TRUE,
      wopt      = wopt_float
    )
  }
  
  weight_files[nm] <- out_w
}

# Reload weight stack
weight_stack <- rast(weight_files)[[raster_names]]

# Rasterize weights directly to disk
#weight_dir <- file.path(out_dir, "weights")
#dir.create(weight_dir, showWarnings = FALSE)

#weight_files <- setNames(character(nlyr(r_mask)), names(r_mask))

#template <- r_mask[[1]]

#for (nm in names(r_mask)) {
  
#  out_w <- file.path(weight_dir, paste0("weight_", nm, ".tif"))
  
#  if (!file.exists(out_w)) {
#    message("Rasterizing weights for ", nm)
    
#    rasterize(
#      ecodist,
#      template,
#      field     = nm,
#      touches   = TRUE,
#      filename  = out_w,
#      overwrite = TRUE,
#      wopt      = wopt_float
#    )
#  }
  
#  weight_files[nm] <- out_w
#}

#rm(ecodist)
#gc()

# Re-open weight stack
#weight_stack <- rast(weight_files)[[names(r_mask)]]

# Check for silent misalignment
stopifnot(nlyr(r_mask) == nlyr(weight_stack))
stopifnot(all(names(r_mask) == names(weight_stack)))

# Weighted probability
prob_file <- file.path(out_dir, "weighted_grassland_probability.tif")

app(
  c(r_mask, weight_stack),
  fun = function(...) {
    
    v <- as.matrix(...)     # dims: cells × layers
    nl <- ncol(v) / 2
    
    x <- v[, 1:nl, drop = FALSE]
    w <- v[, (nl + 1):(2 * nl), drop = FALSE]
    
    # Zero weights where raster is NA
    w[is.na(x)] <- 0
    
    num <- rowSums(x * w, na.rm = TRUE)
    den <- rowSums(w, na.rm = TRUE)
    
    out <- num / den
    out[den == 0] <- NA
    
    pmin(pmax(out, 0), 1)
  },
  filename  = prob_file,
  overwrite = TRUE,
  wopt      = wopt_float
)

# Minimum source constraint
count_file <- file.path(out_dir, "num_contributing_rasters.tif")

num_sources <- app(
  r_mask,
  fun = function(...) rowSums(!is.na(cbind(...))),
  filename  = count_file,
  overwrite = TRUE,
  wopt      = wopt_count
)

minN_file <- file.path(out_dir, "weighted_grassland_probability_minSources.tif")

weighted_probability_minN <- app(
  c(rast(prob_file), num_sources),
  fun = function(p, n) ifelse(n >= 2, p, NA),
  filename  = minN_file,
  overwrite = TRUE,
  wopt      = wopt_float
)

# Binary output and diagnostics (probability 0.5 or higher)
binary_file <- file.path(out_dir, "weighted_grassland_binary_0p5.tif")

binary_map <- app(
  rast(prob_file),
  fun = function(x) x >= 0.5,
  filename  = binary_file,
  overwrite = TRUE,
  wopt      = wopt_byte
)

weight_sum_file <- file.path(out_dir, "sum_weights_present.tif")

app(
  weight_stack,
  fun = function(...) rowSums(cbind(...), na.rm = TRUE),
  filename  = weight_sum_file,
  overwrite = TRUE,
  wopt      = wopt_float
)

#-----------------
# Hierarchical/pooled thresholds with recall-targeted thresholds (to make binary map)
#----------------

# 1. Compute recall‑targeted thresholds
#    Threshold = minimum probability achieving a target recall (e.g. 85%)
# 2. Compute thresholds at multiple spatial levels using hierarchical fallback logic:
#    If ecodistrict has ≥ 25 validation points → ecodistrict threshold
#    Else if ecoregion has ≥ 50 points         → ecoregion threshold
#    Else                                    → global threshold
# 3. Rasterize the final threshold surface
# 4. Apply pixel‑wise adaptive thresholding

# Set parameters

# Minimum sample sizes
MIN_N_ECODIST <- 25
MIN_N_ECOREGION <- 50

# Recall target for grassland
TARGET_RECALL <- 0.85

# Safety bounds on thresholds
THR_MIN <- 0.40
THR_MAX <- 0.60

# Threshold offsets based on num_sources

# Confidence bands (based on contributing rasters)
CONF_BREAKS <- c(2, 4, 7, 11, Inf)

# Threshold offsets for each band
CONF_OFFSETS <- c(
  0.10,   # 2 sources → more strict (require stronger evidence (higher threshold) where fewer datasets contribute)
  0.05,   # 3–4
  0.00,   # 5–7
  -0.025,  # 8–11
  -0.05    # ≥ 12 → more lenient (allow weaker evidence (lower threshold) where many datasets agree)
)

#-------- 1) Load validation points and extract probabilities ------------
# Load validation points
val_df <- read_csv(val_points_file, show_col_types = FALSE)

stopifnot(all(c("x", "y", "reference_grass", "ecodistrict_id") %in% names(val_df)))

# Convert to SpatVector
val_pts <- vect(
  val_df,
  geom = c("x", "y"),
  crs  = crs(weighted_probability_minN)
)

# Extract weighted probability
val_df$weighted_prob <- terra::extract(
  weighted_probability_minN,
  val_pts
)[, 2]

# Drop NA extractions
val_df <- val_df %>% filter(!is.na(weighted_prob))

#-------- 2) Build confidence-offset raster from num_sources -------------
offset_file <- file.path(out_dir, "threshold_offsets_by_num_sources.tif")

offset_raster <- app(
  num_sources,
  fun = function(n) {
    out <- rep(NA_real_, length(n))
    
    out[n == 2]                 <- CONF_OFFSETS[1]
    out[n >= 3 & n <= 4]        <- CONF_OFFSETS[2]
    out[n >= 5 & n <= 7]        <- CONF_OFFSETS[3]
    out[n >= 8 & n <= 11]       <- CONF_OFFSETS[4]
    out[n >= 12]                <- CONF_OFFSETS[5]
    
    out
  },
  filename  = offset_file,
  overwrite = TRUE,
  wopt      = wopt_float
)

#---------- 3) Attach eco-level grouping --------------
# Attach ecoregion from polygons
eco_lookup <- ecodist |>
  as.data.frame() |>
  select(ecodistrict_id, !!ECOREGION_FIELD)

val_df <- val_df |>
  left_join(
    eco_lookup,
    by = "ecodistrict_id"
  ) |>
  rename(ecoregion_id = !!ECOREGION_FIELD)

#------------ 4) Build recall-targeted threshold function -----------
recall_target_threshold <- function(prob, truth, target = 0.85) {
  grass_prob <- prob[truth == 1 & !is.na(prob)]
  if (length(grass_prob) == 0) return(NA_real_)
  as.numeric(quantile(grass_prob, probs = 1 - target, type = 7))
}

#--------- 5) Compute thresholds at each level ----------------

# Ecodistrict level
thr_ecodist <- val_df |>
  group_by(ecodistrict_id) |>
  summarise(
    n_grass = sum(reference_grass == 1),
    thr_ecodist = if (n_grass >= MIN_N_ECODIST)
      recall_target_threshold(weighted_prob, reference_grass, TARGET_RECALL)
    else NA_real_,
    .groups = "drop"
  )

# Ecoregion level
thr_ecoregion <- val_df |>
  group_by(ecoregion_id) |>
  summarise(
    n_grass = sum(reference_grass == 1),
    thr_ecoregion = if (n_grass >= MIN_N_ECOREGION)
      recall_target_threshold(weighted_prob, reference_grass, TARGET_RECALL)
    else NA_real_,
    .groups = "drop"
  )

# Global fallback
thr_global <- recall_target_threshold(
  val_df$weighted_prob,
  val_df$reference_grass,
  TARGET_RECALL
)

#------------ 6) Build final ecodistrict threshold table -----------
thr_table <- ecodist |>
  as.data.frame() |>
  select(ecodistrict_id, !!ECOREGION_FIELD) |>
  rename(ecoregion_id = !!ECOREGION_FIELD) |>
  distinct() |>
  left_join(thr_ecodist, by = "ecodistrict_id") |>
  left_join(thr_ecoregion, by = "ecoregion_id") |>
  mutate(
    threshold_final = case_when(
      !is.na(thr_ecodist) ~ thr_ecodist,
      !is.na(thr_ecoregion) ~ thr_ecoregion,
      TRUE                ~ thr_global
    ),
    threshold_final = pmin(pmax(threshold_final, THR_MIN), THR_MAX)
  ) |>
  select(ecodistrict_id, threshold_final)

#------------ 7) Rasterize ecodistrict-specific thresholds -------------
ecodist_thr <- merge(
  ecodist,
  thr_table,
  by = "ecodistrict_id",
  duplicateGeoms = TRUE
)

thr_raster_file <- file.path(out_dir, "ecodistrict_thresholds.tif")

thr_raster <- rasterize(
  ecodist_thr,
  weighted_probability_minN,
  field     = "threshold_final",
  touches   = TRUE,
  filename  = thr_raster_file,
  overwrite = TRUE,
  wopt      = wopt_float
)

#--------- 8) Combine ecodistrict thresholds + confidence offsets ----------
thr_adj_file <- file.path(
  out_dir,
  "ecodistrict_thresholds_confidence_adjusted.tif"
)

thr_adjusted <- app(
  c(thr_raster, offset_raster),
  fun = function(t, off) {
    pmin(pmax(t + off, THR_MIN), THR_MAX)
  },
  filename  = thr_adj_file,
  overwrite = TRUE,
  wopt      = wopt_float
)

#---------- 9) Apply adaptive thresholding ------------------
binary_conf_file <- file.path(
  out_dir,
  "weighted_grassland_binary_ecodistrict_confidence_adjusted.tif"
)

binary_conf <- app(
  c(weighted_probability_minN, thr_adjusted),
  fun = function(p, t) {
    ifelse(is.na(t), NA, p >= t)
  },
  filename  = binary_conf_file,
  overwrite = TRUE,
  wopt      = wopt_byte
)

#------------------ 10) Diagnostics ------------------

# Source of thresholds - Where were thresholds actually estimated?

# Identify threshold source for each ecodistrict
thr_source_df <- thr_table |>
  left_join(
    thr_ecodist |> select(ecodistrict_id, thr_ecodist),
    by = "ecodistrict_id"
  ) |>
  left_join(
    thr_ecozone |> select(ecozone_id, thr_ecozone),
    by = "ecozone_id"
  ) |>
  mutate(
    threshold_source = case_when(
      !is.na(thr_ecodist) ~ "ecodistrict",
      !is.na(thr_ecozone) ~ "ecozone",
      TRUE                ~ "global"
    )
  )

# Tabulate
table(thr_source_df$threshold_source)

# Threshold vs sample size plot - Are ecodistrict thresholds stable?

# Plot ecodistrict threshold vs grass sample size
plot(
  thr_ecodist$n_grass,
  thr_ecodist$thr_ecodist,
  pch = 16,
  cex = 0.8,
  xlab = "Grass validation points (ecodistrict)",
  ylab = "Recall-targeted threshold",
  main = "Ecodistrict thresholds vs sample size"
)

# Global threshold reference
abline(h = thr_global, col = "red", lwd = 2, lty = 2)

# Stability band (optional, useful)
abline(h = c(THR_MIN, THR_MAX), col = "grey70", lty = 3)

# Sample pixels for diagnostics only
set.seed(1)
diag_df <- spatSample(
  c(thr_raster, offset_raster, thr_adjusted, num_sources),
  size = 50000,
  as.data.frame = TRUE,
  na.rm = TRUE
)

names(diag_df) <- c(
  "thr_base",
  "thr_offset",
  "thr_effective",
  "num_sources"
)

# Effective threshold vs number of sources - Does confidence weighting behave logically?
boxplot(
  thr_effective ~ num_sources,
  data = diag_df,
  outline = FALSE,
  xlab = "Number of contributing rasters",
  ylab = "Effective classification threshold",
  main = "Confidence-weighted thresholds by data support"
)

abline(h = thr_global, col = "red", lwd = 2, lty = 2)

# Distribution of confidence-based threshold adjustments - How strong are the confidence effects?
hist(
  diag_df$thr_effective - diag_df$thr_base,
  breaks = 40,
  col = "grey80",
  border = "white",
  xlab = "Threshold adjustment (effective − base)",
  main = "Distribution of confidence-based threshold adjustments"
)

abline(v = 0, lwd = 2)


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
# 4. Add confidence‑weighted thresholds using num_sources
# 5. Tune recall targets by ecoregion


#####################################################################
#####################################################################

## ADVANCED: two versions that use recall/sensitivity/TPR and give uncertainty

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
prob_A_file <- file.path(out_dir, "prob_grassland_strategyA_recallWeighted.tif")

nl <- nlyr(r_mask)

prob_A <- app(
  c(r_mask, recall_stack),
  fun = function(v) {
    
    x <- v[1:nl]                 # binary predictions
    r <- v[(nl + 1):(2 * nl)]    # recall
    
    r[is.na(x)] <- 0             # drop NA rasters
    pos <- x == 1
    
    den <- sum(r, na.rm = TRUE)
    if (den == 0) return(NA_real_)
    
    num <- sum(r[pos], na.rm = TRUE)
    max(0, min(1, num / den))
  },
  filename  = prob_A_file,
  overwrite = TRUE,
  wopt      = wopt_float
)

# Interpretation:
# Interpreted as recall‑weighted positive evidence (probability‑like index).
# High when many positives with high recall and few/no negatives, regardless of specificity.
# Use when avoiding assumptions about negative reliability.

# A simple uncertainty metric: high when prob~0.5 or when little recall is present
# uncertainty_probOnly is highest near 0.5.
#The composite adds a data sufficiency factor (sum_recall_present), so areas with few/weak 
# inputs are flagged as more uncertain even if the current probability isn’t near 0.5.
uncertainty_A_prob_file <- file.path(out_dir, "uncertainty_strategyA_probOnly.tif")

uncertainty_prob_A <- app(
  prob_A,
  fun = function(p) p * (1 - p),
  filename  = uncertainty_A_prob_file,
  overwrite = TRUE,
  wopt      = wopt_float
)

sum_recall_file <- file.path(out_dir, "sum_recall_present.tif")

sum_recall_present <- app(
  recall_stack,
  fun = function(...) rowSums(cbind(...), na.rm = TRUE),
  filename  = sum_recall_file,
  overwrite = TRUE,
  wopt      = wopt_float
)

uncertainty_A_comp_file <- file.path(out_dir, "uncertainty_strategyA_composite.tif")

app(
  c(uncertainty_prob_A, sum_recall_present),
  fun = function(u, s) u * (1 + 1 / (1 + s)),
  filename  = uncertainty_A_comp_file,
  overwrite = TRUE,
  wopt      = wopt_float
)

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

# Prior: choose a prior grassland probability (global or ecodistrict)
# Use a neutral prior if preferred: 0.5 => logit = 0
prior_p <- 0.5
eps <- 1e-6

logit     <- function(p) log(p / (1 - p))
inv_logit <- function(x) 1 / (1 + exp(-x))

prob_B_file <- file.path(out_dir, "prob_grassland_strategyB_bayesian.tif")

prob_B <- app(
  c(r_mask, recall_stack),
  fun = function(v) {
    
    nl <- length(v) / 2
    x  <- v[1:nl]
    sens <- v[(nl + 1):(2 * nl)]
    
    sens <- pmin(pmax(sens, eps), 1 - eps)
    spec <- pmin(pmax(spec_assumed, eps), 1 - eps)
    
    logLR <- ifelse(
      x == 1, log(sens / (1 - spec)),
      ifelse(x == 0, log((1 - sens) / spec), 0)
    )
    
    inv_logit(logit(prior_p) + sum(logLR, na.rm = TRUE))
  },
  filename  = prob_B_file,
  overwrite = TRUE,
  wopt      = wopt_float
)

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
info_sum_file <- file.path(out_dir, "information_sum_absLogLR.tif")

info_sum <- app(
  c(r_mask, recall_stack),
  fun = function(v) {
    nl <- length(v) / 2
    x <- v[1:nl]
    r <- v[(nl + 1):(2 * nl)]
    r <- pmin(pmax(r, eps), 1 - eps)
    abs(ifelse(x == 1, log(r / (1 - spec_assumed)),
               ifelse(x == 0, log((1 - r) / spec_assumed), 0))) |> sum(na.rm = TRUE)
  },
  filename  = info_sum_file,
  overwrite = TRUE,
  wopt      = wopt_float
)

uncertainty_B_file <- file.path(out_dir, "uncertainty_strategyB_probOnly.tif")

uncertainty_prob_B <- app(
  prob_B,
  fun = function(p) p * (1 - p),
  filename  = uncertainty_B_file,
  overwrite = TRUE,
  wopt      = wopt_float
)

uncertainty_B_comp_file <- file.path(out_dir, "uncertainty_strategyB_composite.tif")

app(
  c(uncertainty_prob_B, info_sum),
  fun = function(u, i) u * (1 + 1 / (1 + i)),
  filename  = uncertainty_B_comp_file,
  overwrite = TRUE,
  wopt      = wopt_float
)

# -----------------------------
# 4) Optional masks (minimum support)
# -----------------------------
min_sources <- 2

prob_A_minN <- app(
  c(rast(prob_A_file), num_contributing),
  fun = function(p, n) ifelse(n >= min_sources, p, NA),
  filename  = file.path(out_dir, "prob_grassland_strategyA_minN.tif"),
  overwrite = TRUE,
  wopt      = wopt_float
)

prob_B_minN <- app(
  c(rast(prob_B_file), num_contributing),
  fun = function(p, n) ifelse(n >= min_sources, p, NA),
  filename  = file.path(out_dir, "prob_grassland_strategyB_minN.tif"),
  overwrite = TRUE,
  wopt      = wopt_float
)

## Ways to improve
# Per‑ecodistrict priors: Get grassland prevalence per ecodistrict and rasterize it to a prior and swap prior_p.
# Per‑raster specificity later: Get TN/FP and extend the table to include specificity; can rasterize it exactly like recall and use per‑layer spec.
# Calibrate a threshold: Use an external validation set and get ROC to pick an optimal threshold for the Bayesian map.

# ============================
# 7) TRI-CLASS SYNTHESIS (0=not, 1=tame, 2=native)
# ============================

#------------- 7.1 Load and align tri-class rasters ------------------------
nt_files <- c(
  aci_nt_file, ngi_nt_file, cgr_nt_file, eccc_nt_file,
  cnvc_nt_file, pudl_nt_file, pli_ap_nt_file, pli_cu_nt_file,
  pli_mg_nt_file, pli_mm_nt_file, gpw_nt_file
)

nt_mask_file <- file.path(out_dir, "tri_class_stack_masked.tif")

if (!file.exists(nt_mask_file) && length(nt_files) > 0) {
  
  nt_stack <- rast(nt_files)
  
  make_short_name <- function(path) {
    stem <- tools::file_path_sans_ext(basename(path))
    if (grepl("^sk_pli_", stem)) {
      sub("^sk_(pli_[^_]+).*", "\\1", stem)
    } else {
      sub("^([^_]+)_.*$", "\\1", stem)
    }
  }
  
  names(nt_stack) <- tolower(vapply(nt_files, make_short_name, character(1)))
  
  # Align to grass mask geometry
  if (!compareGeom(r_mask[[1]], nt_stack[[1]], stopOnError = FALSE)) {
    nt_stack <- project(nt_stack, r_mask[[1]])
  }
  
  # Mask to grass study area (disk-backed)
  nt_stack <- mask(nt_stack, r_mask[[1]],
                   filename  = nt_mask_file,
                   overwrite = TRUE,
                   wopt      = wopt_byte)
  
  rm(nt_stack)
  gc()
}

nt_stack <- if (file.exists(nt_mask_file)) rast(nt_mask_file) else NULL

#-------------- 7.2 Define tri-class accuracy metrics --------------------

# Base accuracy metric used to weight tri-class votes
tri_weight_base_metric <- "f1"
# Options: "recall", "precision", "f1"
#   "recall"  = recall_native / recall_tame or single recall
#   "precision"  = precision_native / precision_tame or single precision
#   "f1" = 

# Weighting strategy
tri_weight_strategy <- "mcc_f1"
# Options:
#   "single" = use chosen base metric only
#   "mcc_f1" = MCC × F1 composite weighting

# Define a suffix that will be used to label outputs:
metric_tag <- paste0(tri_weight_approach, "_", tri_weight_base_metric, "Wt")

#----------- 7.3 Rasterize tri-class accuracy metrics ---------------------

# Choose metrics
# Metric column names
metric_map <- list(
  recall = list(
    native = "recall_native",
    tame   = "recall_tame",
    single = "recall"
  ),
  precision = list(
    native = "precision_native",
    tame   = "precision_tame",
    single = "precision"
  ),
  f1 = list(
    native = "F1_native",
    tame   = "F1_tame",
    single = "F1"
  )
)

if (!tri_weight_base_metric %in% names(metric_map)) {
  stop("Unknown tri_weight_base_metric: ", tri_weight_base_metric)
}

metric_native <- metric_map[[tri_weight_base_metric]]$native
metric_tame   <- metric_map[[tri_weight_base_metric]]$tame
metric_single <- metric_map[[tri_weight_base_metric]]$single

# MCC + F1 composite components
metric_mcc       <- "MCC_1v2"
metric_f1_native <- "F1_2_2x2"
metric_f1_tame   <- "F1_1_2x2"

# Exponents (used only for MCC+F1)
alpha_mcc <- 1.0
beta_f1  <- 1.0

# Read accuracy table
accuracy_nt_df <- readr::read_csv(acc_nt_file, show_col_types = FALSE) %>%
  mutate(raster_name = tolower(stringr::str_trim(raster_name)))

# Detect availability
#has_native_tame  <- all(c(metric_native, metric_tame) %in% names(accuracy_nt_df))
#has_single_metric <- metric_single %in% names(accuracy_nt_df)

#if (is.null(nt_stack)) {
#  has_native_tame  <- FALSE
#  has_single_metric <- FALSE
#}

# Rasterize metrics per layer to disk (branch by accuracy metric chosen)
metric_dir <- file.path(out_dir, paste0("tri_metrics_", metric_tag))
dir.create(metric_dir, showWarnings = FALSE)

template <- r_mask[[1]]

# ============================================================
# Function to rasterize a per‑ecodistrict accuracy metric
# ============================================================
#
# df     : data frame with columns ecodistrict_id, raster_name, and a metric field
# field  : name of the metric column to rasterize (e.g., "recall_native", "W_native")
# prefix : filename prefix used to label output rasters
#
# Returns:
#   A disk‑backed SpatRaster stack with one layer per tri‑class raster,
#   aligned exactly to nt_stack and r_mask geometry.
#
make_metric_stack <- function(df, field, prefix) {
  
  # Pre‑allocate a named character vector to store output file paths.
  # Names correspond to raster names (e.g., aci, cnvc, pli_ap, etc.).
  files <- setNames(character(nlyr(nt_stack)), names(nt_stack))
  
  # Loop over each tri‑class raster layer by name
  for (nm in names(nt_stack)) {
    
    # Construct the on‑disk filename for this raster's metric layer
    f <- file.path(metric_dir, paste0(prefix, "_", nm, ".tif"))
    
    # Only rasterize if the file does not already exist
    # (this makes the pipeline restart‑safe)
    if (!file.exists(f)) {
      
      # Join per‑ecodistrict accuracy values to the polygon layer
      # for this specific raster (nm)
      eco <- merge(
        ecodist,
        df |>
          dplyr::filter(raster_name == nm) |>          # keep rows for this raster
          dplyr::select(ecodistrict_id, value = !!sym(field)),
        by = "ecodistrict_id",
        duplicateGeoms = TRUE
      )
      
      # Rasterize the metric from polygons to the template grid
      # - Each pixel gets the metric value of its ecodistrict
      # - touches = TRUE ensures full ecodistrict coverage
      # - filename makes this disk‑backed (no RAM accumulation)
      rasterize(
        eco,
        template,
        field     = "value",
        touches   = TRUE,
        filename  = f,
        overwrite = TRUE,
        wopt      = wopt_float
      )
    }
    
    # Store the filename (existing or newly created)
    files[nm] <- f
  }
  
  # Load all generated files as a single SpatRaster stack
  # Reorder explicitly to match nt_stack layer order
  rast(files)[[names(nt_stack)]]
}

# Initialize metric stacks (only one or two will be used, depending on strategy)
metric_native_stack <- NULL
metric_tame_stack   <- NULL
metric_single_stack <- NULL

# Proceed only if tri‑class rasters exist
if (!is.null(nt_stack)) {
  
  # Keep only accuracy rows corresponding to available tri‑class rasters
  acc_tri <- accuracy_nt_df |>
    dplyr::filter(raster_name %in% names(nt_stack))
  
  # ------------------------------------------------------------
  # Approach: use a single base metric (recall / precision / F1)
  # ------------------------------------------------------------
  if (tri_weight_approach == "single") {
    
    # Case 1: class‑specific metrics exist (native and tame separately)
    if (all(c(metric_native, metric_tame) %in% names(acc_tri))) {
      
      # Rasterize native‑class metric
      metric_native_stack <- make_metric_stack(
        acc_tri, metric_native, "baseNat"
      )
      
      # Rasterize tame‑class metric
      metric_tame_stack <- make_metric_stack(
        acc_tri, metric_tame, "baseTam"
      )
      
      # Case 2: only a single metric applies to all positive predictions
    } else if (metric_single %in% names(acc_tri)) {
      
      metric_single_stack <- make_metric_stack(
        acc_tri, metric_single, "baseAll"
      )
      
      # Case 3: no usable metric found → stop with clear error
    } else {
      stop("Requested base metric not found in tri-class accuracy table.")
    }
    
    # ------------------------------------------------------------
    # Approach: MCC × F1 composite weighting
    # ------------------------------------------------------------
  } else if (tri_weight_aproach == "mcc_f1") {
    
    # Ensure required columns are present
    required <- c(metric_mcc, metric_f1_native, metric_f1_tame)
    if (!all(required %in% names(acc_tri))) {
      stop(
        "MCC+F1 strategy requires columns: ",
        paste(required, collapse = ", ")
      )
    }
    
    # Compute composite weights per ecodistrict × raster
    acc_tri <- acc_tri |>
      dplyr::mutate(
        # Composite weight for native predictions
        W_native = ( !!sym(metric_mcc) )^alpha_mcc *
          ( !!sym(metric_f1_native) )^beta_f1,
        
        # Composite weight for tame predictions
        W_tame   = ( !!sym(metric_mcc) )^alpha_mcc *
          ( !!sym(metric_f1_tame) )^beta_f1,
        
        # Remove non‑finite or negative values
        W_native = ifelse(is.finite(W_native), pmax(W_native, 0), NA),
        W_tame   = ifelse(is.finite(W_tame),   pmax(W_tame,   0), NA)
      )
    
    # Rasterize composite weights for native and tame votes
    metric_native_stack <- make_metric_stack(
      acc_tri, "W_native", "mccF1Nat"
    )
    
    metric_tame_stack <- make_metric_stack(
      acc_tri, "W_tame", "mccF1Tam"
    )
  }
}

#------------- 7.4 Compute native/tame shares ---------------------
share_native_file <- file.path(
  out_dir,
  paste0("share_native_", metric_tag, ".tif")
)

share_native <- if (!is.null(nt_stack)) {
  
  nl <- nlyr(nt_stack)
  
  app(
    c(nt_stack,
      metric_native_stack %||% metric_single_stack,
      metric_tame_stack   %||% metric_single_stack),
    fun = function(v) {
      
      x  <- v[1:nl]                 # tri-class values
      w1 <- v[(nl + 1):(2 * nl)]    # native weights
      w2 <- v[(2 * nl + 1):(3 * nl)]# tame weights
      
      present <- !is.na(x)
      native  <- x == 2
      tame    <- x == 1
      
      w1[!present] <- 0
      w2[!present] <- 0
      
      nv <- sum(w1[native], na.rm = TRUE)
      tv <- sum(w2[tame],   na.rm = TRUE)
      
      tot <- nv + tv
      if (tot == 0) return(NA_real_)
      
      p <- nv / tot
      max(0, min(1, p))
    },
    filename  = share_native_file,
    overwrite = TRUE,
    wopt      = wopt_float
  )
  
} else NULL

share_tame <- if (!is.null(share_native)) 1 - share_native else NULL

#------------- 7.5 Fallback prior for missing tri evidence ---------------
prior_native_global <- 0.5

if (!"prior_native" %in% names(ecodist)) {
  ecodist$prior_native <- prior_native_global
}

prior_native_r <- rasterize(
  ecodist, r_mask[[1]],
  field = "prior_native",
  touches = TRUE,
  filename = file.path(out_dir, paste0("prior_native_", metric_tag, ".tif")),
  overwrite = TRUE,
  wopt = wopt_float
)

prior_native_r <- clamp(prior_native_r, 0, 1)

share_native <- if (is.null(share_native)) {
  prior_native_r
} else {
  ifel(is.na(share_native), prior_native_r, share_native)
}

share_tame <- 1 - share_native

#-------------- 7.6 Compose final class probabilities -----------------
p_grass <- rast(prob_A_file)  # or prob_B_file
p_grass <- clamp(p_grass, 0, 1)

p_native   <- p_grass * share_native
p_tame     <- p_grass * share_tame
p_nongrass <- 1 - p_grass

#------------- 7.7 Build hard 3-class map ------------------
p_stack <- c(p_nongrass, p_tame, p_native)
names(p_stack) <- c("p_not", "p_tame", "p_native")

idx_max <- which.max(p_stack)

class_3cat <- classify(
  idx_max,
  rcl = matrix(c(1,0, 2,1, 3,2), ncol = 2, byrow = TRUE),
  filename  = file.path(out_dir, paste0("grass_3class_StratA_", metric_tag, ".tif")), # or StratB
  overwrite = TRUE,
  wopt      = wopt_byte
)

#------------ 7.8 Calculate uncertainty metrics ----------------
uncertainty_1minus <- app(
  p_stack,
  fun = function(...) 1 - max(cbind(...)),
  filename  = file.path(out_dir, paste0("uncertainty_1minus_maxP_", metric_tag, ".tif")),
  overwrite = TRUE,
  wopt      = wopt_float
)

eps <- 1e-9
uncertainty_entropy <- app(
  p_stack,
  fun = function(pn, pt, pnative) {
    H <- -(
      pn      * log(pn      + eps) +
        pt      * log(pt      + eps) +
        pnative * log(pnative + eps)
    ) / log(3)
    min(max(H, 0), 1)
  },
  filename  = file.path(out_dir, paste0("uncertainty_entropy_norm_", metric_tag, ".tif")),
  overwrite = TRUE,
  wopt      = wopt_float
)
