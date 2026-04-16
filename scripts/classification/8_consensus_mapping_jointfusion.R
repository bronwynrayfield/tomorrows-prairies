
###################################################################################
##                     GRASSLAND CONSENSUS MAPPING ANALYSIS                      ##
###################################################################################

##----------------------- Alison Long - April 2026 -------------------------

# step 8: use data fusion techniques to create a consensus grassland map 
#         to serve as baseline for modelling

################################################################################
# STRATEGY A: JOINT 3‑CLASS FUSION WITH EMPIRICAL AGREEMENT PENALTIES
#
# Purpose:
#   - Produce joint probabilities for {not grass, tame grass, native grass}
#   - Use binary rasters (grass / not-grass) AND tri‑class rasters together
#   - Handle raster dependence explicitly using precomputed agreement metrics
#   - Non‑Bayesian: no priors, no likelihoods, no log‑odds
#
# Key Ideas:
#   1) Binary rasters inform grass presence only
#   2) Tri‑class rasters inform native vs tame directly
#   3) All votes are weighted by:
#        accuracy_weight × redundancy_penalty
#   4) Redundant rasters are downweighted but never removed
################################################################################

library(terra)
library(dplyr)

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

# Study area (modelling extent + Peace River region)
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

nt_files <- c(
  aci_nt_file, ngi_nt_file, cgr_nt_file, eccc_nt_file,
  cnvc_nt_file, pudl_nt_file, pli_ap_nt_file, pli_cu_nt_file,
  pli_mg_nt_file, pli_mm_nt_file, gpw_nt_file
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

# Binary rasters: 0 = not grass, 1 = grass, NA
bin_stack <- r_mask

#------------- 7.1 Load and align tri-class rasters ------------------------

# Tri‑class rasters: 0 = not, 1 = tame, 2 = native, NA

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

#===============================================================================
# 0. REQUIRED INPUTS 
#===============================================================================

# Accuracy rasters (already rasterized by ecodistrict)
recall_bin_stack    <- recall_stack_bin
recall_tame_stack   <- recall_stack_tame
recall_native_stack <- recall_stack_native

# Empirical agreement summaries (already computed by you)
# Higher values = more redundant / more correlated
agreement_bin <- tibble(
  raster_name = names(bin_stack),
  agree_score = AGREEMENT_VALUES_BIN     # numeric vector you already have
)

agreement_nt <- tibble(
  raster_name = names(nt_stack),
  agree_score = AGREEMENT_VALUES_NT      # numeric vector you already have
)

#===============================================================================
# 1. CONVERT AGREEMENT → REDUNDANCY PENALTIES
#===============================================================================

# Tuning parameters
alpha <- 1.0          # strength of correlation penalty
PENALTY_FLOOR <- 0.25 # never allow a raster to contribute less than this share

# Penalize binary rasters
penalty_bin_vec <- agreement_bin |>
  mutate(
    penalty = 1 / (1 + alpha * agree_score),
    penalty = pmax(penalty, PENALTY_FLOOR)
  ) |>
  deframe()

# Penalize tri‑class rasters
penalty_nt_vec <- agreement_nt |>
  mutate(
    penalty = 1 / (1 + alpha * agree_score),
    penalty = pmax(penalty, PENALTY_FLOOR)
  ) |>
  deframe()

#===============================================================================
# 2. RASTERIZE PENALTY FACTORS (CONSTANT PER LAYER)
#===============================================================================

penalty_bin_stack <- rast(
  lapply(names(bin_stack), function(nm) {
    setValues(rast(bin_stack[[1]]), penalty_bin_vec[nm])
  })
)
names(penalty_bin_stack) <- names(bin_stack)

penalty_nt_stack <- if (!is.null(nt_stack)) {
  r <- rast(
    lapply(names(nt_stack), function(nm) {
      setValues(rast(nt_stack[[1]]), penalty_nt_vec[nm])
    })
  )
  names(r) <- names(nt_stack)
  r
} else NULL

#===============================================================================
# 3. JOINT 3‑CLASS FUSION — STRATEGY A (NON‑BAYESIAN)
#===============================================================================

# Numerical safeguards
eps <- 1e-6
NEG_PENALTY <- 1.0   # strength of penalty to grass subclasses when predicting not-grass

out_prob_file <- file.path(
  out_dir,
  "prob_joint_3class_StrategyA_empiricalPenalty.tif"
)

joint_prob_A <- app(
  c(
    bin_stack,
    recall_bin_stack,
    penalty_bin_stack,
    nt_stack,
    recall_tame_stack,
    recall_native_stack,
    penalty_nt_stack
  ),
  fun = function(v) {
    
    #----------------------------
    # Unpack inputs
    #----------------------------
    nb <- nlyr(bin_stack)
    nt <- if (!is.null(nt_stack)) nlyr(nt_stack) else 0
    idx <- 1
    
    x_bin <- if (nb > 0) v[idx:(idx + nb - 1)] else numeric(0)
    idx <- idx + nb
    
    w_bin <- if (nb > 0) v[idx:(idx + nb - 1)] else numeric(0)
    idx <- idx + nb
    
    p_bin <- if (nb > 0) v[idx:(idx + nb - 1)] else numeric(0)
    idx <- idx + nb
    
    x_nt <- if (nt > 0) v[idx:(idx + nt - 1)] else numeric(0)
    idx <- idx + nt
    
    w_tame <- if (nt > 0) v[idx:(idx + nt - 1)] else numeric(0)
    idx <- idx + nt
    
    w_nat <- if (nt > 0) v[idx:(idx + nt - 1)] else numeric(0)
    idx <- idx + nt
    
    p_nt <- if (nt > 0) v[idx:(idx + nt - 1)] else numeric(0)
    
    #----------------------------
    # Initialize evidence
    #----------------------------
    E_not <- 0
    E_tame <- 0
    E_native <- 0
    
    #----------------------------
    # Binary rasters
    #----------------------------
    for (i in seq_len(nb)) {
      if (is.na(x_bin[i])) next
      
      w_eff <- w_bin[i] * p_bin[i]
      if (w_eff <= 0) next
      
      if (x_bin[i] == 1) {
        # Grass → split vote between tame and native
        E_tame   <- E_tame   + w_eff / 2
        E_native <- E_native + w_eff / 2
      } else {
        # Not grass → vote for not, penalize sub‑classes
        E_not    <- E_not    + w_eff
        E_tame   <- E_tame   - NEG_PENALTY * w_eff / 2
        E_native <- E_native - NEG_PENALTY * w_eff / 2
      }
    }
    
    #----------------------------
    # Tri‑class rasters
    #----------------------------
    for (i in seq_len(nt)) {
      if (is.na(x_nt[i])) next
      
      if (x_nt[i] == 0) {
        E_not <- E_not + max(w_tame[i], w_nat[i]) * p_nt[i]
      }
      if (x_nt[i] == 1) {
        E_tame <- E_tame + w_tame[i] * p_nt[i]
      }
      if (x_nt[i] == 2) {
        E_native <- E_native + w_nat[i] * p_nt[i]
      }
    }
    
    #----------------------------
    # Normalize to probabilities
    #----------------------------
    E_not    <- max(E_not, 0)
    E_tame   <- max(E_tame, 0)
    E_native <- max(E_native, 0)
    
    S <- E_not + E_tame + E_native
    if (S <= eps) return(c(NA_real_, NA_real_, NA_real_))
    
    c(E_not / S, E_tame / S, E_native / S)
  },
  filename  = out_prob_file,
  overwrite = TRUE,
  wopt      = wopt_float
)

names(joint_prob_A) <- c("p_not", "p_tame", "p_native")

#===============================================================================
# 4. HARD 3‑CLASS MAP (ARGMAX)
#===============================================================================

class_3_file <- file.path(
  out_dir,
  "grass_3class_joint_StrategyA_empiricalPenalty.tif"
)

class_3_A <- app(
  joint_prob_A,
  fun = function(pn, pt, pnat) which.max(c(pn, pt, pnat)) - 1,
  filename  = class_3_file,
  overwrite = TRUE,
  wopt      = wopt_byte
)

# Class codes:
#   0 = not grass
#   1 = tame grass
#   2 = native grass

#===============================================================================
# 5. UNCERTAINTY DIAGNOSTICS
#===============================================================================

# Normalized entropy
entropy_file <- file.path(
  out_dir,
  "uncertainty_entropy_joint3class_StrategyA_empiricalPenalty.tif"
)

app(
  joint_prob_A,
  fun = function(pn, pt, pnat) {
    H <- -(
      pn * log(pn + eps) +
        pt * log(pt + eps) +
        pnat * log(pnat + eps)
    ) / log(3)
    min(max(H, 0), 1)
  },
  filename  = entropy_file,
  overwrite = TRUE,
  wopt      = wopt_float
)

# Margin uncertainty (1 − max P)
margin_file <- file.path(
  out_dir,
  "uncertainty_1minus_maxP_joint3class_StrategyA_empiricalPenalty.tif"
)

app(
  joint_prob_A,
  fun = function(...) 1 - max(cbind(...)),
  filename  = margin_file,
  overwrite = TRUE,
  wopt      = wopt_float
)