
###################################################################################
##                     GRASSLAND CONSENSUS MAPPING ANALYSIS                      ##
###################################################################################

##----------------------- Alison Long - July 2026 -------------------------

# step 9: use data fusion techniques to create a consensus grassland map 
#         to serve as baseline for modelling

# THIS VERSION IS OPTIMIZED FOR SAGEMAKER
# Next optimization to try: build stacks directly from S3 using GDAL VRT (zero downloads)

# Strengths
#   Accuracy-weighted
#   Temporally weighted
#   Ecodistrict-specific performance weighting
#   Explicit uncertainty products
#   Explicit support metric
#   Conservative masking
#   Computationally scalable

# Potential concerns
#   Softmax may overly flatten probabilities 
#      - Could make confidence lower than underlying evidence suggest
#      - May want to compare to normalized values (ex: p_native <- E_native / E_sum)
#   Binary evidence may dominate native/tame separation
#   Support-score constants are uncalibrated
#   ESS no longer represents independent evidence
#   High-confidence layer may become very sparse

# THIS VERSION (compared to v7):
# Implements new approach

# THIS VERSION (compared to v6):
# Replaces MCC as the accuracy metric with the following:
# For binary: F1_grass (Overall_Accuracy as fallback)
# For triclass: Macro_F1 (Balanced_Accuracy as fallback)

# THIS VERSION (compared to v5):
# Removes use of penalty rasters based on redundancy (from spatial agreement analysis)
# and instead excludes grassland maps that are not independent of each other.
# Alternative approach: group maps into families (independent evidence sources) and 
# allocate a fixed weight to each family

# THIS VERSION (compared to v4):
# Adds a pixel-wise effective sample size (ESS) output and combines it with entropy/confidence
# to represent how much independent weighted evidence (binary) supports each pixel. Also
# combines ESS with total weight to calculate a support score

# THIS VERSION (compared to v3):
# Crops/masks to pilot area
# Updated to fix errors encountered while running

# THIS VERSION (compared to v2):
# Has higher interpretability and stability
# Uses normalized weight scaling - ensures total influence stays bounded, no single raster predominates
# Uses signed contributions instead of votes - symmetric, easier to interpret, numerically stable
# Has no need for negative penalties
# Uses softmax probability mapping instead of linear - always positive, stable when evidence small 
# Probabilities now mean “relative weighted support across all rasters after accounting for accuracy and redundancy” instead of "scaled vote totals"

# Possible upgrades:
# Dataset-level weighting (avoid over-represented sources)
# Entropy-calibrated masking
# Pseudo-likelihood interpretation (bridge to Bayesian)
# Hybrid version that is mathematically equivalent to a Bayesian model but keeps current structure and speed

################################################################################
# STRATEGY A: JOINT 3‑CLASS FUSION 
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
#   4) Redundant rasters are excluded from analysis
################################################################################

## ---- Setup ---- ##
# Temp dirs
tmp_root <- file.path(Sys.getenv("HOME"), "r_tmp")
#dir.create(tmp_root, recursive = TRUE, showWarnings = FALSE)

# Load packages
library(terra)
library(readr)
library(tidyverse)
library(aws.s3)

# Set terra options
terraOptions(
  tempdir = file.path(tmp_root, "terra"),
  memfrac = 0.6, # could go down to 0.4
  progress = 1
)

# increase GDAL block cache size to 8192 MB (8 GB)
Sys.setenv(GDAL_CACHEMAX = 8192)

# Force a fresh run instead of silently re-using old scratch outputs
FORCE_REBUILD <- TRUE

## ---- Settings ---- ##

out_dir    <- "s3://sagemaker-gst-stage.sharing/alison-long/tomorrows-prairies/results/consensus_map"
log_dir    <- "s3://sagemaker-gst-stage.sharing/alison-long/tomorrows-prairies/logs"

terra_dir    <- file.path(tmp_root, "terra")
local_dir    <- file.path(tmp_root, "data")
stack_dir    <- file.path(tmp_root, "stacks")
diag_dir     <- file.path(tmp_root, "diagnostics")
dir.create(diag_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(local_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(stack_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(terra_dir, recursive = TRUE, showWarnings = FALSE)

meta_file    <- "s3://sagemaker-gst-stage.sharing/alison-long/tomorrows-prairies/processed/metadata_full_extent_binary_nattam_with_path_withaci2024_nohfi_cec_luts10_aci2025_pglc_nodups.csv"
ecodist_file <- "s3://sagemaker-gst-stage.sharing/alison-long/tomorrows-prairies/processed/boundaries/national-ecological-framework-for-canada/Ecodistricts/ecodistricts_clean_rast.tif"
pilot_file   <- "s3://sagemaker-gst-stage.sharing/alison-long/tomorrows-prairies/processed/pilot-area/pilot_area_aligned_cnvc_Jan2026.tif"
sa_file      <- "s3://sagemaker-gst-stage.sharing/alison-long/tomorrows-prairies/processed/study-area/modelling_extent_aligned_cnvc_Jan2026.tif"
acc_bin_file <- "s3://sagemaker-gst-stage.sharing/alison-long/tomorrows-prairies/processed/accuracy_metrics/binary_accuracy_ecoregprovfallback_atomic_combinedCIandCM_28Jul2026.csv"
acc_nt_file  <- "s3://sagemaker-gst-stage.sharing/alison-long/tomorrows-prairies/processed/accuracy_metrics/nattam_accuracy_atomic_combinedCIandCM_forconsensus_29Jul2026_final.csv"

## ---- Message log ---- ##
run_ts <- format(Sys.time(), "%Y%m%d_%H%M%S")

log_file <- file.path(log_dir, paste0("consensus_mapping_log_", run_ts, ".log"))
log_file_local <- file.path(tmp_root, paste0("consensus_mapping_", run_ts, ".log"))

log_msg <- function(..., level = "INFO") {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  msg <- paste0("[", level, "] ", ts, " : ", paste(..., collapse = ""))
  message(msg)
  cat(msg, "\n", file = log_file_local, append = TRUE)
}

on.exit({
  if (file.exists(log_file_local)) {
    tryCatch(
      aws.s3::put_object(file = log_file_local, object = log_file),
      error = function(e) message("Could not upload log file: ", conditionMessage(e))
    )
  }
}, add = TRUE)

log_msg("---- SESSION INFO ----")

script_name <- tryCatch({
  cmd_args <- commandArgs(trailingOnly = FALSE)
  sub("^--file=", "", cmd_args[grep("^--file=", cmd_args)])
}, error = function(e) NA_character_)

log_msg("Script name: ", script_name)

session_lines <- capture.output(sessionInfo())

for (ln in session_lines) {
  log_msg(ln)
}

log_msg("R version: ", R.version.string)
log_msg("Platform: ", R.version$platform)
log_msg("Working directory: ", getwd())
log_msg("terra version: ", as.character(packageVersion("terra")))
log_msg("tidyverse version: ", as.character(packageVersion("tidyverse")))
log_msg("aws.s3 version: ", packageVersion("aws.s3"))
log_msg("readr version: ", packageVersion("readr"))
log_msg("Temp root: ", tmp_root)
log_msg("GDAL_CACHEMAX: ", Sys.getenv("GDAL_CACHEMAX"))
log_msg("FORCE_REBUILD: ", FORCE_REBUILD)
log_msg("Sys.info: ",
        paste(names(Sys.info()), Sys.info(), collapse = "; "))
log_msg("terra options:")
capture.output(print(terraOptions())) |>
  lapply(log_msg)
gc_out <- capture.output(gc())
log_msg(paste(gc_out, collapse = " | "))

log_msg("---- END SESSION INFO ----")

log_msg("Starting consensus mapping analysis")

## ---- Load metadata ---- ##
log_msg("Reading metadata from: ", meta_file)

meta <- aws.s3::s3read_using(read_csv, object = meta_file, show_col_types = FALSE) %>%
  rename(file = file_s3)

## ---- Validate metadata ---- ##
log_msg("Metadata rows loaded: ", nrow(meta))

valid_types <- c("binary", "nt")

if (!all(meta$stack_type %in% valid_types)) {
  bad <- unique(meta$stack_type[!meta$stack_type %in% valid_types])
  stop("Invalid stack_type values found: ", paste(bad, collapse = ", "))
}

if (any(duplicated(meta$raster_id))) {
  dupes <- meta$raster_id[duplicated(meta$raster_id)]
  stop("Duplicate raster_id values found: ", paste(unique(dupes), collapse=", "))
}

if (any(is.na(meta$file)) || any(meta$file == "")) {
  stop("Missing file path(s) in metadata for raster_id(s): ",
       paste(meta$raster_id[is.na(meta$file) | meta$file == ""], collapse = ", "))
}

## ---- Temporal settings ---- ##
target_year <- 2025 # use newest dataset or policy-relevant year
lambda <- 6  # temporal decay scale in years; larger values retain influence from older datasets longer

# time decay function
time_decay <- function(year, target_year, lambda) {
  exp(-abs(year - target_year) / lambda)
}

# Assign weights
meta <- meta %>%
  mutate(
    time_weight = time_decay(
      year,
      target_year,
      lambda))

## ---- Split by stack type ---- ##
bin_meta <- meta %>% filter(stack_type == "binary")
nt_meta  <- meta %>% filter(stack_type == "nt")

## ---- Named file vectors ---- ##
#rfiles_s3 <- setNames(meta$file,     meta$raster_id)
bin_files <- setNames(bin_meta$file, bin_meta$raster_id)
nt_files  <- setNames(nt_meta$file,  nt_meta$raster_id)

log_msg("Binary rasters: ", nrow(bin_meta))
log_msg("Tri-class rasters: ", nrow(nt_meta))

## ---- Load pilot area ---- ##
log_msg("Loading pilot area")

pilot_local <- file.path(local_dir, "pilot_area.tif")

aws.s3::save_object(
  object = pilot_file,
  file = pilot_local
)

pilot <- rast(pilot_local)

log_msg("Pilot area loaded")

## ---- Load study area ---- ##
log_msg("Loading study area")

sa_local <- file.path(local_dir, "study_area.tif")

aws.s3::save_object(
  object = sa_file,
  file = sa_local
)

sa <- rast(sa_local)

log_msg("Study area loaded")

## ---- Download rasters locally ---- ##
message("Downloading rasters locally...")
log_msg("Downloading rasters locally")

#download_files <- function(files_named) {
#  sapply(names(files_named), function(nm) {
#    dest <- file.path(local_dir, paste0(nm, ".tif"))
#    if (!file.exists(dest)) {
#      aws.s3::save_object(files_named[[nm]], file = dest)
#    }
#    dest
#  }, USE.NAMES = TRUE)
#}

download_files <- function(files_named) {
  
  if (length(files_named) == 0)
    return(setNames(character(0), character(0)))
  
  sapply(names(files_named), function(nm) {
    
    raw_file <- file.path(local_dir, paste0(nm, "_raw.tif"))
    crop_file <- file.path(local_dir, paste0(nm, ".tif"))
    
    if (!file.exists(crop_file) || FORCE_REBUILD) {
      
      # download original
      aws.s3::save_object(files_named[[nm]], file = raw_file)
      
      # crop + mask to pilot
      r <- rast(raw_file)
      
      r <- crop(r, pilot)
      r <- mask(r, pilot)
      
      # crop + mask to study area
      # r <- rast(raw_file)
      # 
      # r <- crop(r, sa)
      # r <- mask(r, sa)
      
      writeRaster(
        r,
        crop_file,
        overwrite = TRUE,
        wopt = list(
          gdal = c(
            "COMPRESS=ZSTD",
            "TILED=YES",
            "BIGTIFF=YES"
          )
        )
      )
      
      rm(r)
      gc()
    }
    
    crop_file
    
  }, USE.NAMES = TRUE)
  
}

bin_local <- download_files(bin_files)
nt_local  <- download_files(nt_files)

log_msg("All rasters downloaded locally")

## ---- Confirm raster alignment ---- ##
message("Checking alignment...")
log_msg("Checking raster alignment")

check_alignment <- function(files) {
  template <- rast(files[[1]])
  
  all(vapply(files, function(f) {
    compareGeom(template, rast(f), stopOnError = FALSE)
  }, logical(1)))
}

if (!check_alignment(bin_local)) stop("Binary rasters not aligned")
if (length(nt_local) > 0 && !check_alignment(nt_local)) {
  stop("Tri-class rasters not aligned")
}

message("Rasters confirmed aligned — skipping resampling")
log_msg("Rasters confirmed aligned — skipping resampling")

## ---- Align rasters to a common grid (if needed) ---- ##
#message("Aligning rasters...")

#template <- rast(all_local[[1]])

#aligned_files <- lapply(names(all_local), function(nm) {
#  f <- all_local[[nm]]
#  r <- rast(f)

#  if (!compareGeom(r, template, stopOnError = FALSE)) {
#    r <- resample(r, template, method = "near")
#  }

#  out <- file.path(local_dir, paste0(nm, "_aligned.tif"))
#  writeRaster(r, out, overwrite = TRUE)
#  out
#})

#names(aligned_files) <- names(all_local)

## ---- Build raster stack ---- ##
log_msg("Building optimized raster stacks")

wopt_stack <- list(
  datatype = "INT1U",
  gdal = c("COMPRESS=ZSTD", "TILED=YES", "BIGTIFF=YES")
)

build_stack <- function(local_files, out_path, force) {
  if (!file.exists(out_path) || force) {
    writeRaster(rast(local_files), out_path, overwrite = TRUE, wopt = wopt_stack)
  }
  r <- rast(out_path)
  names(r) <- names(local_files)
  r
}

bin_stack_file <- file.path(stack_dir, "bin_stack.tif")
bin_stack <- build_stack(bin_local, bin_stack_file, FORCE_REBUILD)

nt_stack_file <- file.path(stack_dir, "nt_stack.tif")
nt_stack <- if (length(nt_local) > 0) {
  build_stack(nt_local, nt_stack_file, FORCE_REBUILD)
} else {
  NULL
}
has_nt <- !is.null(nt_stack)

if (!has_nt) {
  stop("No tri-class rasters available")
}

## ---- Create temporal-weight stacks ----##

#bin_time <- bin_meta$time_weight
#names(bin_time) <- names(bin_local)

#stopifnot(
#  identical(
#    names(bin_time),
#    names(bin_stack)
#  )
#)

#time_bin_stack <- rast(
#  lapply(bin_local, function(f) {
#    r <- rast(f)
#    values(r) <- bin_time[which(bin_local == f)]
#    r
#  })
#)

#nt_time <- nt_meta$time_weight
#names(nt_time) <- names(nt_local)

#stopifnot(
#  identical(
#    names(nt_time),
#    names(nt_stack)
#  )
#)

#time_nt_stack <- rast(
#  lapply(nt_local, function(f) {
#    r <- rast(f)
#    values(r) <- nt_time[which(nt_local == f)]
#    r
#  })
#)

# =========================
# GLOBAL OPTIONS
# =========================
#terraOptions(
#  memfrac = 0.8,
#  tempdir = "/mnt/tmp", # use fast scratch disk
#  progress = 1
#)

# =========================
# INTEGRATE ACCURACY DATA
# =========================

## ---- Load binary accuracy table ---- ##
log_msg("Loading binary accuracy table")

acc_bin <- aws.s3::s3read_using(
  read_csv,
  object = acc_bin_file,
  show_col_types = FALSE) %>%
  mutate(
    raster_name = paste0(
      tolower(trimws(raster)),
      "_bin"))

## ---- Load tri-class accuracy table ---- ##
log_msg("Loading tri-class accuracy table")

acc_nt <- aws.s3::s3read_using(
  read_csv,
  object = acc_nt_file,
  show_col_types = FALSE) %>%
  mutate(
    raster_name = paste0(
      tolower(trimws(raster)),
      "_nt"))

## ---- Validate binary accuracy table ---- ##
log_msg("Binary accuracy rows loaded: ", nrow(acc_bin))

req_cols <- c(
  "raster_name",
  "Ecds_ID",
  "acc"
)

if (!all(req_cols %in% names(acc_bin))) {
  stop(
    "Missing required columns in binary accuracy table: ",
    paste(setdiff(req_cols, names(acc_bin)), collapse = ", ")
  )
}

acc_bin <- acc_bin %>%
  mutate(
    Ecds_ID = as.integer(round(Ecds_ID))
  )

dup_check <- acc_bin %>%
  count(raster_name, Ecds_ID) %>%
  filter(n > 1)

if (nrow(dup_check) > 0) {
  stop(
    "Duplicate raster/ecodistrict combinations found in acc_bin"
  )
}

missing_bin_names <- setdiff(
  unique(acc_bin$raster_name),
  meta$raster_id
)

if(length(missing_bin_names) > 0)
  stop("Raster names do not match")

## ---- Validate tri-class accuracy table ---- ##
log_msg("Tri-class accuracy rows loaded: ", nrow(acc_nt))

req_cols_nt <- c(
  "raster_name",
  "Ecds_ID",
  "acc",
  "native_metric",
  "tame_metric",
  "not_metric"
)

if (!all(req_cols_nt %in% names(acc_nt))) {
  
  stop(
    "Missing required columns in tri-class accuracy table: ",
    paste(
      setdiff(req_cols_nt, names(acc_nt)),
      collapse = ", "
    )
  )
}

acc_nt <- acc_nt %>%
  mutate(
    Ecds_ID = as.integer(round(Ecds_ID))
  )

dup_check_nt <- acc_nt %>%
  count(raster_name, Ecds_ID) %>%
  filter(n > 1)

if (nrow(dup_check_nt) > 0) {
  stop(
    "Duplicate raster/ecodistrict combinations found in acc_nt"
  )
}

missing_nt_names <- setdiff(
  unique(acc_nt$raster_name),
  meta$raster_id
)

if(length(missing_nt_names) > 0)
  stop("Raster names do not match")

## ---- Add temporal weighting ---- ##
log_msg("Applying accuracy + temporal weighting")

# acc_to_weight <- function(x) {
#   w <- (x + 1) / 2
#   w[!is.finite(w)] <- NA
#   pmax(w, 1e-4)
# }

if (any(acc_bin$acc < 0 | acc_bin$acc > 1, na.rm = TRUE)) {
  stop("acc_bin contains values outside [0,1]")
}

if (any(acc_nt$acc < 0 | acc_nt$acc > 1, na.rm = TRUE)) {
  stop("acc_nt contains values outside [0,1]")
}

acc_to_weight <- function(x) {
  w <- x
  w[!is.finite(w)] <- NA
  pmax(w, 1e-4)
}

#time_decay <- function(year, target_year, lambda) {
#  exp(-lambda * abs(year - target_year))
#}

acc_bin <- acc_bin %>%
  left_join(
    meta %>%
      select(
        raster_id,
        year
      ),
    by = c(
      "raster_name" = "raster_id"
    )
  )

if (any(is.na(acc_bin$year))) {
  stop("Missing year information in acc_bin")
}

acc_bin <- acc_bin %>%
  mutate(
    w_acc  = acc_to_weight(acc),
    w_time = time_decay(
      year,
      target_year,
      lambda
    ),
    weight = w_acc * w_time
  )

## ---- Add temporal weighting ---- ##
log_msg("Applying tri-class metric weighting")

acc_nt <- acc_nt %>%
  left_join(
    meta %>%
      select(
        raster_id,
        year
      ),
    by = c(
      "raster_name" = "raster_id"
    )
  )

if (any(is.na(acc_nt$year))) {
  stop("Missing year information in acc_nt")
}

# make sure F1 scores are between 0 and 1
# stopifnot(
#   all(acc_nt$F1_1 >= 0 & acc_nt$F1_1 <= 1),
#   all(acc_nt$F1_2 >= 0 & acc_nt$F1_2 <= 1)
# )

# acc_nt <- acc_nt %>%
#   mutate(
#     w_acc = acc_to_weight(acc),
#     
#     w_time = time_decay(
#       year,
#       target_year,
#       lambda
#     ),
#     
#     W_tame   = w_acc * w_time,
#     W_native = w_acc * w_time
#   )

# validate class-specific F1s
if (any(acc_nt$native_metric < 0 | acc_nt$native_metric > 1, na.rm = TRUE)) {
  stop("acc_nt$native_metric contains values outside [0,1]")
}

if (any(acc_nt$tame_metric < 0 | acc_nt$tame_metric > 1, na.rm = TRUE)) {
  stop("acc_nt$tame_metric contains values outside [0,1]")
}

if (any(acc_nt$not_metric < 0 | acc_nt$not_metric > 1, na.rm = TRUE)) {
  stop("acc_nt$not_metric contains values outside [0,1]")
}

if(any(is.na(acc_nt$native_metric)))
  stop("native_metric still contains NAs")

if(any(is.na(acc_nt$tame_metric)))
  stop("tame_metric still contains NAs")

if(any(is.na(acc_nt$not_metric)))
  stop("not_metric still contains NAs")

acc_nt <- acc_nt %>%
  mutate(
    
    # native_metric =
    #   coalesce(
    #     F1_2,
    #     Macro_F1,
    #     Balanced_Accuracy
    #   ),
    # 
    # tame_metric =
    #   coalesce(
    #     F1_1,
    #     Macro_F1,
    #     Balanced_Accuracy
    #   ),
    # 
    # not_metric =
    #   coalesce(
    #     F1_0,
    #     Macro_F1,
    #     Balanced_Accuracy
    #   ),
    
    W_native = native_metric * w_time,
    W_tame   = tame_metric * w_time,
    W_not    = not_metric   * w_time
  )

# acc_nt <- acc_nt %>%
#   mutate(
#     
#     w_time = time_decay(
#       year,
#       target_year,
#       lambda
#     ),
#     
#     # class-specific weights
#     W_tame   = pmax(F1_1, 1e-4) * w_time,
#     W_native = pmax(F1_2, 1e-4) * w_time,
#     W_not    = pmax(Macro_F1, 1e-4) * w_time # must decide whether to use Macro_F1 or F1_0
#   )                                   # Macro_F1 = balanced class influence, prevents nongrass from dominating
#                                       # F1_0 = statistically pure, nongrass will receive larger weights

## ---- Validate tri-class raster coverage ---- ##

missing_nt <- setdiff(
  names(nt_stack),
  acc_nt$raster_name
)

if (length(missing_nt) > 0) {
  stop(
    "Missing tri-class accuracy values for: ",
    paste(
      missing_nt,
      collapse = ", "
    )
  )
}

## ---- Load ecodistrict raster ---- ##
log_msg("Loading pilot-clipped ecodistrict raster")

ecodist_r <- rast(ecodist_file)

ecodist_r <- crop(ecodist_r, pilot)
ecodist_r <- mask(ecodist_r, pilot)

# ecodist_r <- crop(ecodist_r, sa)
# ecodist_r <- mask(ecodist_r, sa)

ecodist_r <- as.int(ecodist_r)

if (!compareGeom(ecodist_r, bin_stack, stopOnError=FALSE)) {
  stop("ecodist_r does not match bin_stack geometry")
}

## ---- Validate ecodistrict coverage ---- ##
log_msg("Validating ecodistrict coverage")

ids_r <- unique(freq(ecodist_r)[, "value"])

na_count <- sum(is.na(ids_r))

log_msg("NA values in ecodistrict raster IDs: ", na_count)

if (na_count > 0) {
  warning("Ecodistrict raster contains ", na_count, " NA values")
}

# if there are NAs, use:
ids_r <- ids_r[!is.na(ids_r)]

ids_tbl_bin <- unique(acc_bin$Ecds_ID)
missing_ids_bin <- setdiff(ids_r, ids_tbl_bin)

ids_tbl_nt <- unique(acc_nt$Ecds_ID)
missing_ids_nt <- setdiff(ids_r, ids_tbl_nt)

if (length(missing_ids_bin) > 0) {
  stop("Ecodistrict IDs missing from binary accuracy table: ",
       paste(
         utils::head(missing_ids_bin),
         collapse = ", "
       ))}

if (length(missing_ids_nt) > 0) {
  stop("Ecodistrict IDs missing from triclass accuracy table: ",
       paste(
         utils::head(missing_ids_nt),
         collapse = ", "
       ))}

## ---- Build accuracy weight stack helper ---- ##

build_weight_stack <- function(
    accuracy_table,
    stack_names,
    value_col,
    ecodist_r,
    out_prefix
) {
  
  weight_files <- lapply(
    stack_names,
    function(nm) {
      
      log_msg(
        "Building weight raster: ",
        nm
      )
      
      lookup <- accuracy_table %>%
        filter(
          raster_name == nm
        ) %>%
        select(
          Ecds_ID,
          all_of(value_col)
        ) %>%
        as.matrix()
      
      if (nrow(lookup) == 0) {
        stop(
          "No weights found for raster: ",
          nm
        )
      }
      
      out <- file.path(
        terra_dir,
        paste0(
          out_prefix,
          "_",
          nm,
          ".tif"
        )
      )
      
      if (!file.exists(out) ||
          FORCE_REBUILD) {
        
        #w <- classify(
        #  ecodist_r,
        #  rcl = lookup
        #)
        
        w <- subst(
          ecodist_r,
          from = lookup[,1],
          to   = lookup[,2]
        )
        
        writeRaster(
          w,
          out,
          overwrite = TRUE,
          wopt = list(
            datatype = "FLT4S"
          )
        )
      }
      
      out
      
    }
  )
  
  s <- rast(
    unlist(weight_files)
  )
  
  names(s) <- stack_names
  
  s
}

## ---- Build binary accuracy stack ---- ##
log_msg("Building binary accuracy weight stack")

acc_bin_stack <- build_weight_stack(
  accuracy_table = acc_bin,
  stack_names    = names(bin_stack),
  value_col      = "weight",
  ecodist_r      = ecodist_r,
  out_prefix     = "bin_acc"
)

stopifnot(
  identical(
    names(acc_bin_stack),
    names(bin_stack)
  )
)

## ---- Build tame weight stack ---- ##
log_msg("Building tame weight stack")

metric_tame_stack <- build_weight_stack(
  accuracy_table = acc_nt,
  stack_names    = names(nt_stack),
  value_col      = "W_tame",
  ecodist_r      = ecodist_r,
  out_prefix     = "tame"
)

stopifnot(
  identical(
    names(metric_tame_stack),
    names(nt_stack)
  )
)

## ---- Build native weight stack ---- ##
log_msg("Building native weight stack")

metric_native_stack <- build_weight_stack(
  accuracy_table = acc_nt,
  stack_names    = names(nt_stack),
  value_col      = "W_native",
  ecodist_r      = ecodist_r,
  out_prefix     = "native"
)

stopifnot(
  identical(
    names(metric_native_stack),
    names(nt_stack)
  )
)

log_msg("Building nongrass weight stack")

metric_not_stack <- build_weight_stack(
  accuracy_table = acc_nt,
  stack_names    = names(nt_stack),
  value_col      = "W_not",
  ecodist_r      = ecodist_r,
  out_prefix     = "not"
)

stopifnot(
  identical(
    names(metric_not_stack),
    names(nt_stack)
  )
)


#mcc_bin_stack    <- rast("mcc_bin.tif")
#mcc_tame_stack   <- rast("mcc_tame.tif")
#mcc_native_stack <- rast("mcc_native.tif")

#penalty_bin_stack <- rast("penalty_bin.tif")
#penalty_nt_stack  <- rast("penalty_nt.tif")

#out_dir <- "tiles_out"
#dir.create(out_dir, showWarnings = FALSE)

# Confirm weights are within expected bounds
tmp_chk1 <- global(acc_bin_stack, c("min","max"), na.rm=TRUE)
log_msg(capture.output(tmp_chk1))

tmp_chk2 <- global(metric_tame_stack, c("min","max"), na.rm = TRUE)
log_msg(capture.output(tmp_chk2))

tmp_chk3 <- global(metric_native_stack, c("min","max"), na.rm = TRUE)
log_msg(capture.output(tmp_chk3))

tmp_chk4 <- global(metric_not_stack, c("min","max"), na.rm = TRUE)
log_msg(capture.output(tmp_chk4))

# =========================
# TILE GRID
# =========================
# ~2000–5000 tiles is ideal, 500-1500 for single node
tiles <- makeTiles(pilot, n = 10, overlap = 0) # for sa, make n = 1000

log_msg("Tile count: ", length(tiles))

# Save tiles for reproducibility
writeVector(vect(tiles), "tile_index.gpkg", overwrite=TRUE)

# =========================
# PROCESS FUNCTION
# =========================
tile_dir <- "/mnt/tmp/tiles"
dir.create(tile_dir, showWarnings=FALSE)

# Define maximum entropy
max_entropy <- log(3)

process_tile <- function(i) {
  
  tile <- tiles[[i]]
  
  out_file <- file.path(tile_dir, sprintf("tile_%05d.tif", i))
  if (file.exists(out_file)) return(out_file)
  
  # crop all rasters
  b   <- crop(bin_stack, tile)
  
  # skip empty tiles
  if (is.na(global(b[[1]], "max")[[1]])) {
    log_msg(sprintf("Tile %d skipped (no data)", i))
    return(NULL)
  }
  
  nt  <- crop(nt_stack, tile)
  
  #time_bin <- crop(time_bin_stack, tile)
  #time_nt  <- crop(time_nt_stack, tile)
  
  rbin <- crop(acc_bin_stack, tile)
  rt   <- crop(metric_tame_stack, tile)
  rn   <- crop(metric_native_stack, tile)
  rnot <- crop(metric_not_stack, tile)
  
  #pbin <- crop(penalty_bin_stack, tile)
  #pnt  <- crop(penalty_nt_stack, tile)
  
  # effective weight = accuracy (no more redundancy penalty)
  #bin_w_eff   <- rbin * pbin * time_bin
  #nt_w_tame   <- rt * pnt * time_nt
  #nt_w_native <- rn * pnt * time_nt
  
  bin_present <- !is.na(b)
  
  bin_w_eff <- rbin # * pbin
  
  bin_w_eff_present <- bin_w_eff * bin_present
  
  #for (j in seq_len(nlyr(bin_w_eff))) {
  #  bin_w_eff[[j]] <- bin_w_eff[[j]] * bin_time[j]
  #}
  
  nt_present <- !is.na(nt)
  
  nt_w_dataset <- (rt + rn + rnot) / 3
  
  all_weights <- c(bin_w_eff_present, nt_w_dataset * nt_present)
  
  nt_w_tame <- rt # * pnt
  nt_w_native <- rn # * pnt
  
  nt_w_tame_present   <- nt_w_tame   * nt_present
  nt_w_native_present <- nt_w_native * nt_present
  
  nt_not_present <- rnot * nt_present
  
  nt_dataset_weight <- (rt + rn + rnot)/3
  
  # ----------------------
  # GLOBAL NORMALIZATION
  # ----------------------
  # create global denominators - counts nt 3 times
  #global_weight_sum <- app(c(bin_w_eff_present, nt_w_tame_present,
  #    nt_w_native_present, nt_not_present), sum, na.rm = TRUE)
  
  # create global denominators
  global_weight_sum <- app(
    c(
      bin_w_eff_present,
      nt_dataset_weight * nt_present
    ),
    sum,
    na.rm = TRUE
  )
  
  # Normalize all evidence against the same denominator
  bin_w_norm <- ifel(
    global_weight_sum > 0,
    bin_w_eff_present / global_weight_sum,
    0
  )
  
  nt_w_tame_norm <- ifel(
    global_weight_sum > 0,
    nt_w_tame_present / global_weight_sum,
    0
  )
  
  nt_w_native_norm <- ifel(
    global_weight_sum > 0,
    nt_w_native_present / global_weight_sum,
    0
  )
  
  nt_not_norm <- ifel(
    global_weight_sum > 0,
    nt_not_present / global_weight_sum,
    0
  )
  
  # ----------------------
  # TRICLASS WEIGHT RATIOS
  # ----------------------
  
  ratio_not_dataset <- ifel(
    nt_dataset_weight > 0,
    rnot / nt_dataset_weight,
    NA
  )
  
  ratio_native_dataset <- ifel(
    nt_dataset_weight > 0,
    rn / nt_dataset_weight,
    NA
  )
  
  ratio_tame_dataset <- ifel(
    nt_dataset_weight > 0,
    rt / nt_dataset_weight,
    NA
  )
  
  names(ratio_not_dataset)    <- "ratio_not_dataset"
  names(ratio_native_dataset) <- "ratio_native_dataset"
  names(ratio_tame_dataset)   <- "ratio_tame_dataset"
  
  max_class_ratio <- app(
    c(
      ratio_not_dataset,
      ratio_native_dataset,
      ratio_tame_dataset
    ),
    max,
    na.rm = TRUE
  )
  
  names(max_class_ratio) <- "max_class_ratio"
  
  # ----------------------
  # DATASET COUNTS
  # ----------------------
  num_bin_sources <- app(bin_present, sum, na.rm = TRUE)
  names(num_bin_sources) <- "num_bin_sources"
  
  num_nt_sources <- app(nt_present, sum, na.rm = TRUE)
  names(num_nt_sources) <- "num_nt_sources"
  
  # ----------------------
  # EFFECTIVE SAMPLE SIZE (average class-specific reliability + temporal weighting)
  # ----------------------
  
  # Sum of weights - binary only
  #w_sum <- app(bin_w_eff_present, sum, na.rm = TRUE)
  
  # Sum of weights
  w_sum <- app(all_weights, sum, na.rm = TRUE)
  
  # total effective weight contributing to the pixel (for export)
  sum_weight_present <- w_sum
  names(sum_weight_present) <- "sum_weight_present"
  
  # Sum of squared weights
  w_sum_sq <- app(all_weights^2, sum, na.rm = TRUE)
  
  # Effective sample size - amount of dataset-level support (not class evidence)
  ess <- w_sum^2 / w_sum_sq
  
  # Handle divide-by-zero
  ess[w_sum_sq == 0] <- NA
  
  names(ess) <- "ess"
  
  # support score - uses ESS and total weight
  support_score <- (
    ess / (ess + 3)
  ) * (
    w_sum / (w_sum + 1)
  )
  
  names(support_score) <- "support_score"
  
  # Convert ESS into uncertainty metric
  # Higher ESS = lower uncertainty
  ess_uncertainty <- 1 / (1 + ess)
  
  names(ess_uncertainty) <- "ess_uncertainty"
  
  # Convert support into uncertainty metric
  support_uncertainty <- 1 - support_score
  
  names(support_uncertainty) <- "support_uncertainty"
  
  # ----------------------
  # SUPPORT PROVENANCE
  # ----------------------
  
  binary_support_total <- app(
    bin_w_eff_present,
    sum,
    na.rm = TRUE
  )
  
  triclass_support_total <- app(
    nt_dataset_weight * nt_present,
    sum,
    na.rm = TRUE
  )
  
  fraction_binary_support <- ifel(
    global_weight_sum > 0,
    binary_support_total / global_weight_sum,
    NA
  )
  
  fraction_triclass_support <- ifel(
    global_weight_sum > 0,
    triclass_support_total / global_weight_sum,
    NA
  )
  
  names(fraction_binary_support)   <- "fraction_binary_support"
  names(fraction_triclass_support) <- "fraction_triclass_support"
  
  # ----------------------
  # EVIDENCE COMPUTATION
  # ----------------------
  
  # ----------------------
  # NORMALIZED WEIGHTS
  # ----------------------
  
  # binary weights
  #bin_w_sum <- app(bin_w_eff_present, sum, na.rm=TRUE)
  #bin_w_norm <- ifel(bin_w_sum > 0, bin_w_eff_present / bin_w_sum, 0)
  
  # tri-class weights
  #nt_w_sum <- app(nt_w_tame_present + nt_w_native_present, sum, na.rm=TRUE)
  #nt_w_tame_norm   <- ifel(nt_w_sum > 0, nt_w_tame_present / nt_w_sum, 0)
  #nt_w_native_norm <- ifel(nt_w_sum > 0, nt_w_native_present / nt_w_sum, 0)
  
  # ----------------------
  # HIERARCHICAL EVIDENCE
  # ----------------------
  
  # Binary evidence
  E_bin_grass <- app(bin_w_norm * (b == 1), sum, na.rm = TRUE)
  E_bin_not   <- app(bin_w_norm * (b == 0), sum, na.rm = TRUE)
  
  # ----------------------
  # Triclass class-specific evidence
  # ----------------------
  
  #nt_not_sum  <- app(nt_not_present, sum, na.rm = TRUE)
  #nt_not_norm <- ifel(nt_not_sum > 0, nt_not_present / nt_not_sum, 0)
  
  E_not_tri    <- app(nt_not_norm * (nt == 0), sum, na.rm = TRUE)
  E_tame_tri   <- app(nt_w_tame_norm * (nt == 1), sum, na.rm = TRUE)
  E_native_tri <- app(nt_w_native_norm * (nt == 2), sum, na.rm = TRUE)
  
  # ----------------------
  # Grass vs nongrass
  # ----------------------
  
  Grass_support <- E_bin_grass + E_tame_tri + E_native_tri
  
  Not_support <- E_bin_not + E_not_tri
  
  grass_total <- Grass_support + Not_support
  
  p_grass <- Grass_support / grass_total
  p_grass <- ifel(grass_total > 0, p_grass, NA)
  
  p_not <- Not_support / grass_total
  p_not <- ifel(grass_total > 0, p_not, NA)
  
  # ----------------------
  # Native vs tame within grass
  # ----------------------
  
  grass_type_total <- E_tame_tri + E_native_tri
  
  # Pixels with grass evidence but no native/tame evidence
  #native_tame_missing <- ifel(grass_type_total == 0, 1, 0)
  native_tame_missing <- ifel((grass_type_total == 0) & (p_grass > 0.5),
                              1, 0)
  names(native_tame_missing) <- "native_tame_missing"
  
  # Shares
  share_native <- E_native_tri / grass_type_total
  share_tame   <- E_tame_tri / grass_type_total
  
  # If no native/tame evidence exists
  share_native <- ifel(grass_type_total > 0, share_native, NA)
  share_tame   <- ifel(grass_type_total > 0, share_tame, NA)
  
  # Optional neutral fallback
  neutral <- p_grass * 0 + 0.5
  
  share_native <- cover(share_native, neutral)
  
  share_tame <- cover(share_tame, neutral)
  
  # ----------------------
  # Final class probabilities
  # ----------------------
  
  p_native <- p_grass * share_native
  p_tame   <- p_grass * share_tame
  
  # numerical cleanup
  prob_sum <- p_not + p_tame + p_native
  
  #p_not    <- p_not    / prob_sum
  #p_tame   <- p_tame   / prob_sum
  #p_native <- p_native / prob_sum
  
  p_not <- ifel(prob_sum > 0, p_not / prob_sum, NA)
  p_tame <- ifel(prob_sum > 0, p_tame / prob_sum, NA)
  p_native <- ifel(prob_sum > 0, p_native / prob_sum, NA)
  
  # ----------------------
  # BINARY VOTE COUNTS
  # ----------------------
  
  num_positive <- app(bin_present * (b == 1), sum, na.rm = TRUE)
  names(num_positive) <- "num_positive"
  
  num_negative <- app(bin_present * (b == 0), sum, na.rm = TRUE)
  names(num_negative) <- "num_negative"
  
  vote_total <- num_positive + num_negative
  
  vote_margin <- ifel(
    vote_total > 0,
    abs(num_positive - num_negative) / vote_total,
    NA
  )
  
  names(vote_margin) <- "vote_margin"
  
  weighted_vote_total <- Grass_support + Not_support
  
  weighted_vote_margin <- ifel(
    weighted_vote_total > 0,
    abs(Grass_support - Not_support) /
      weighted_vote_total,
    NA
  )
  
  names(weighted_vote_margin) <- "weighted_vote_margin"
  
  # ----------------------
  # TRICLASS VOTE COUNTS
  # ----------------------
  
  num_not_votes <- app(
    nt_present * (nt == 0),
    sum,
    na.rm = TRUE
  )
  names(num_not_votes) <- "num_not_votes"
  
  num_tame_votes <- app(
    nt_present * (nt == 1),
    sum,
    na.rm = TRUE
  )
  names(num_tame_votes) <- "num_tame_votes"
  
  num_native_votes <- app(
    nt_present * (nt == 2),
    sum,
    na.rm = TRUE
  )
  names(num_native_votes) <- "num_native_votes"
  
  # votes
  vote_stack <- c(
    num_not_votes,
    num_tame_votes,
    num_native_votes
  )
  
  max_votes <- app(
    vote_stack,
    max
  )
  
  total_votes <-
    num_not_votes +
    num_tame_votes +
    num_native_votes
  
  triclass_agreement <- ifel(
    total_votes > 0,
    max_votes / total_votes,
    NA
  )
  
  names(triclass_agreement) <- "triclass_agreement"
  
  vote_margin_tri <- app(
    c(
      num_not_votes,
      num_tame_votes,
      num_native_votes
    ),
    fun = function(x) {
      
      x <- sort(x, decreasing = TRUE)
      
      total <- sum(x)
      
      if (total == 0)
        return(NA)
      
      (x[1] - x[2]) / total
    }
  )
  
  names(vote_margin_tri) <- "vote_margin_tri"
  
  tri_support <- c(
    E_not_tri,
    E_tame_tri,
    E_native_tri
  )
  
  max_support <- app(
    tri_support,
    max
  )
  
  total_support <-
    E_not_tri +
    E_tame_tri +
    E_native_tri
  
  weighted_agreement_tri <- ifel(
    total_support > 0,
    max_support / total_support,
    NA
  )
  
  names(weighted_agreement_tri) <- "weighted_agreement_tri"
  
  # ----------------------
  # Entropy
  # ----------------------
  
  eps <- 1e-6
  
  entropy <- -(
    p_not * log(p_not + eps) +
      p_tame * log(p_tame + eps) +
      p_native * log(p_native + eps)
  )
  
  # Normalize entropy (0 to 1)
  entropy_norm <- entropy / max_entropy
  entropy_norm <- clamp(entropy_norm, 0, 1)
  
  names(entropy_norm) <- "entropy_norm"
  
  # stack outputs
  out <- c(p_not, p_tame, p_native, entropy, entropy_norm, ess, ess_uncertainty, support_score,
           support_uncertainty, native_tame_missing, num_bin_sources, num_nt_sources,
           sum_weight_present, num_positive, num_negative, vote_margin, weighted_vote_margin,
           num_not_votes, num_tame_votes, num_native_votes, triclass_agreement, vote_margin_tri, 
           weighted_agreement_tri, ratio_not_dataset, ratio_native_dataset, ratio_tame_dataset,
           max_class_ratio, fraction_binary_support, fraction_triclass_support)
  
  names(out) <- c("p_not", "p_tame", "p_native", "entropy", "entropy_norm",
                  "ess", "ess_uncertainty", "support_score", "support_uncertainty",
                  "native_tame_missing", "num_bin_sources", "num_nt_sources",
                  "sum_weight_present", "num_positive", "num_negative", "vote_margin", 
                  "weighted_vote_margin", "num_not_votes", "num_tame_votes", "num_native_votes", 
                  "triclass_agreement", "vote_margin_tri", "weighted_agreement_tri",
                  "ratio_not_dataset", "ratio_native_dataset", "ratio_tame_dataset",
                  "max_class_ratio", "fraction_binary_support", "fraction_triclass_support")
  
  writeRaster(out, out_file, overwrite=TRUE)
  
  return(out_file)
}

# =========================
# MOSAIC AND EXPORT
# =========================
results <- lapply(seq_along(tiles), process_tile)

print("All tiles completed.")

# Get tile outputs
files <- Filter(Negate(is.null), results)

if (length(files) == 0) {
  files <- list.files("/mnt/tmp/tiles", full.names=TRUE, pattern="\\.tif$")
}

# Ensure files exist
files <- files[file.exists(files)]

if (length(files) == 0) {
  stop("No tile files found for mosaicking")
}

cat("Mosaicking", length(files), "tiles...\n")

# Check layer consistency - make sure correct number of layers
#stopifnot(all(sapply(files, function(f) nlyr(rast(f)) == 4)))
stopifnot(all(vapply(files, function(f) nlyr(rast(f)) == 29, logical(1))))

# Mosaic all tiles together
mosaic_out <- mosaic(rast(files), fun="first")
# or can try
#mosaic_out <- do.call(mosaic, c(lapply(files, rast), list(fun="mean")))

# Set names
names(mosaic_out) <- c("p_not", "p_tame", "p_native", "entropy", "entropy_norm",
                       "ess", "ess_uncertainty", "support_score", "support_uncertainty",
                       "native_tame_missing", "num_bin_sources", "num_nt_sources",
                       "sum_weight_present", "num_positive", "num_negative", "vote_margin", 
                       "weighted_vote_margin", "num_not_votes", "num_tame_votes", "num_native_votes", 
                       "triclass_agreement", "vote_margin_tri", "weighted_agreement_tri",
                       "ratio_not_dataset", "ratio_native_dataset", "ratio_tame_dataset",
                       "max_class_ratio", "fraction_binary_support", "fraction_triclass_support")

#------ Derived layers ------#
p_not    <- mosaic_out[["p_not"]]
p_tame   <- mosaic_out[["p_tame"]]
p_native <- mosaic_out[["p_native"]]

# Probability stack
prob_stack <- mosaic_out[[c("p_not", "p_tame", "p_native")]]
entropy <- mosaic_out[["entropy"]]

# support score
support_score <- mosaic_out[["support_score"]]
support_uncertainty <- mosaic_out[["support_uncertainty"]]

ess <- mosaic_out[["ess"]]
ess_uncertainty <- mosaic_out[["ess_uncertainty"]]

entropy_norm <- mosaic_out[["entropy_norm"]]

# Missing
native_tame_missing <- mosaic_out[["native_tame_missing"]]

# number of sources and votes
num_bin_sources <- mosaic_out[["num_bin_sources"]]
num_nt_sources  <- mosaic_out[["num_nt_sources"]]

sum_weight_present <- mosaic_out[["sum_weight_present"]]

num_positive <- mosaic_out[["num_positive"]]
num_negative <- mosaic_out[["num_negative"]]

vote_margin <- mosaic_out[["vote_margin"]]
weighted_vote_margin <- mosaic_out[["weighted_vote_margin"]]

num_not_votes <- mosaic_out[["num_not_votes"]]
num_tame_votes <- mosaic_out[["num_tame_votes"]]
num_native_votes <- mosaic_out[["num_native_votes"]] 

triclass_agreement <- mosaic_out[["triclass_agreement"]]

vote_margin_tri <- mosaic_out[["vote_margin_tri"]]
weighted_agreement_tri <- mosaic_out[["weighted_agreement_tri"]]

ratio_not_dataset <- mosaic_out[["ratio_not_dataset"]]
ratio_native_dataset <- mosaic_out[["ratio_native_dataset"]]
ratio_tame_dataset <- mosaic_out[["ratio_tame_dataset"]]
max_class_ratio <- mosaic_out[["max_class_ratio"]]

fraction_binary_support <- mosaic_out[["fraction_binary_support"]]
fraction_triclass_support <- mosaic_out[["fraction_triclass_support"]]

# Check
#hist(values(support_score))
vals_support <- spatSample(
  support_score,
  size = 1e6,
  na.rm = TRUE,
  values = TRUE
)

quantile(vals_support)

png(
  file.path(diag_dir, "support_score_hist.png"),
  width = 2400,
  height = 1600,
  res = 300
)

hist(vals_support,
  main = "Support Score Distribution",
  xlab = "Support Score"
)

dev.off()
#hist(spatSample(support_score, size = 1e6, na.rm = TRUE, values = TRUE))
#quantile(values(spatSample(support_score, 1e6, na.rm = TRUE)))

# Confidence = max probability
confidence <- app(prob_stack, max)
names(confidence) <- "confidence"

# Check probabilities sum to ~1
check <- app(prob_stack, sum)
stats <- global(check, c("min", "max", "mean"), na.rm=TRUE)
print(stats)
log_msg(capture.output(stats))

if (stats[1, "min"] < 0.99 ||
    stats[1, "max"] > 1.01) {
  stop("Probability normalization failed")
}

# check drivers of masking
tmp_chk5 <- global(entropy, c("min","max","mean"), na.rm = TRUE)
log_msg(capture.output(tmp_chk5))

tmp_chk6 <- global(support_score, c("min","max","mean"), na.rm = TRUE)
log_msg(capture.output(tmp_chk6))

tmp_chk7 <- global(p_not, c("min","max","mean"), na.rm = TRUE)
log_msg(capture.output(tmp_chk7))

tmp_chk8 <- global(p_tame, c("min","max","mean"), na.rm = TRUE)
log_msg(capture.output(tmp_chk8))

tmp_chk9 <- global(p_native, c("min","max","mean"), na.rm = TRUE)
log_msg(capture.output(tmp_chk9))

# Convert confidence to uncertainty
confidence_uncertainty <- 1 - confidence

names(confidence_uncertainty) <- "confidence_uncertainty"

# adjusted confidence
#confidence <- app(prob_stack, max)

decision_score <- confidence * support_score
names(decision_score) <- "decision_score"

# Build a combined uncertainty layer
# Equal weighting: classifcation confidence, class ambiguity, amount of supporting eveidence
uncertainty_total <- (
  confidence_uncertainty +
    entropy_norm +
    support_uncertainty
) / 3

names(uncertainty_total) <- "uncertainty_total"

# Weight entropy and confidence more than support:
uncertainty_total_w <- (
  0.35 * confidence_uncertainty +
    0.35 * entropy_norm +
    0.30 * support_uncertainty
)

names(uncertainty_total_w) <- "uncertainty_total_weighted"

#------ Hard classification ------#
hard_class <- which.max(prob_stack)

# Check to finetune confidence threshold
vals_conf <- spatSample(
  decision_score,
  size = 1e6,
  na.rm = TRUE,
  values = TRUE)

png(
  file.path(diag_dir, "decision_score_hist.png"),
  width = 2400,
  height = 1600,
  res = 300
)

hist(
  vals_conf,
  main = "Decision Score Distribution",
  xlab = "Decision Score"
)

dev.off()

# hard classification (1=not, 2=tame, 3=native)
# mask low-confidence or high-uncertainty pixels
# max_prob: require >0.6 confidence (fairly strict) - avoid overconfident misclassification
#           May be too high for softmax outputs - check using histogram above
# entropy: uncertainty masking - remove unreliable areas
# max entropy for 3 classes ≈ 1.10, 0.8=moderately conservative
# could auto-calibrate thresholds based on data distribution instead of fixed values
# goal: keep coverage high, penalize low-support areas smoothly instead of abruptly, 
# and prevent a large number of tiny-weight rasters from producing deceptively good ESS values
hard_class <- ifel(
  (decision_score < 0.4) | (entropy > 0.8),
  NA,
  hard_class
)

# could add rule for how to deal with ties
# currently slightly biased toward “non_grass” in rare tie case

names(hard_class) <- "class"

# Optional labeling (may not persist in GeoTIFF)
levels(hard_class) <- data.frame(
  value = 1:3,
  class = c("non_grass", "tame", "native")
)

#------ High-confidence-only classification ------#
# Conservative/strict: confidence ≥ 0.8, entropy ≤ 0.5
# Balanced: confidence ≥ 0.75, entropy ≤ 0.6
# Inclusive: confidence ≥ 0.7, entropy ≤ 0.7
high_conf_class <- ifel(
  (decision_score >= 0.75) & (entropy <= 0.6),
  hard_class,
  NA)

names(high_conf_class) <- "high_conf_class"

# optional labeling
levels(high_conf_class) <- data.frame(
  value = 1:3,
  class = c("non_grass", "tame", "native"))

###################################################

#------ Tiled temporal workflow (grass vs nongrass) ------#

# Check for uneven year representation
year_counts <- table(bin_meta$year)

log_msg("Rasters per year:")
print(year_counts)

if (any(year_counts < 2)) {
  warning(
    "Some years are represented by fewer than 2 rasters."
  )
}

# Aggregate if some years only have 1-2 rasters
bin_meta_pd <- bin_meta %>%
  mutate(
    period = case_when(
      year <= 2017 ~ "2015_2017",
      year <= 2020 ~ "2018_2020",
      year <= 2023 ~ "2021_2023",
      TRUE         ~ "2024_2025"
    )
  )

# Build annual stack
years <- sort(unique(bin_meta$year))

year_rasters <- lapply(years, function(y) {
  
  idx <- which(bin_meta$year == y)
  
  if (length(idx) == 1) {
    rast(bin_local[idx])
  } else {
    app(
      rast(bin_local[idx]),
      mean,
      na.rm = TRUE
    )
  }
})

year_stack <- rast(year_rasters)

names(year_stack) <- paste0("y", years)

# OR aggregated periods stack
periods <- sort(unique(bin_meta_pd$period))

period_rasters <- lapply(periods, function(p) {
  
  idx <- which(bin_meta_pd$period == p)
  
  app(
    rast(bin_local[idx]),
    mean,
    na.rm = TRUE
  )
})

year_stack_pds <- rast(period_rasters)

names(year_stack_pds) <- periods

# CHOOSE WHETHER TO AGGREGATE YEARS OR NOT
USE_PERIODS <- TRUE

temporal_stack <- if (USE_PERIODS) {
  log_msg("Using period-aggregated temporal stack")
  year_stack_pds
} else {
  log_msg("Using annual temporal stack")
  year_stack
}

# Create temporary tile directory
temp_tile_dir <- "/mnt/tmp/temporal_tiles"
dir.create(temp_tile_dir, showWarnings = FALSE)

# Tile processor
process_temporal_tile <- function(i, temporal_stack) {
  
  tile <- tiles[[i]]
  
  out_file <- file.path(
    temp_tile_dir,
    sprintf("temporal_%05d.tif", i)
  )
  
  if (file.exists(out_file)) {
    return(out_file)
  }
  
  ys <- crop(temporal_stack, tile)
  
  if (is.na(global(ys[[1]], "max")[[1]])) {
    return(NULL)
  }
  
  # -----------------------
  # Temporal variance
  # -----------------------
  
  temporal_var <- app(
    ys,
    var,
    na.rm = TRUE
  )
  
  temporal_sd <- sqrt(temporal_var)
  names(temporal_sd) <- "temporal_sd"
  
  # -----------------------
  # Temporal entropy
  # -----------------------
  
  grass_freq <- app(
    ys,
    mean,
    na.rm = TRUE
  )
  
  eps <- 1e-6
  
  temporal_entropy <- -(
    grass_freq * log(grass_freq + eps) +
      (1 - grass_freq) *
      log(1 - grass_freq + eps)
  )
  
  names(temporal_entropy) <- "temporal_entropy"
  
  # -----------------------
  # Temporal range
  # -----------------------
  
  temporal_range <- app(
    ys,
    function(x) {
      max(x, na.rm = TRUE) -
        min(x, na.rm = TRUE)
    }
  )
  
  names(temporal_range) <- "temporal_range"
  
  # -----------------------
  # Temporal trend - will be slowest operation! few points if years aggregated! cut if needed
  # -----------------------
  
  #yrs <- years
  
  # temporal_slope <- app(
  #   ys,
  #   function(x) {
  
  #     if (all(is.na(x))) {
  #       return(NA)
  #     }
  
  #     if (sum(!is.na(x)) < 3) {
  #       return(NA)
  #     }
  
  #     ok <- !is.na(x)
  
  #     if (sum(ok) < 3) {
  #       return(NA)
  #     }
  
  # could be very slow
  #coef(
  #  lm(
  #    x[ok] ~ yrs[ok]
  #  )
  #)[2]
  
  # faster to compute slope analytically
  #     cov(x[ok], yrs[ok])/var(yrs[ok])
  #   }
  # )
  
  #names(temporal_slope) <- "temporal_slope"
  
  out <- c(
    temporal_sd,
    temporal_entropy,
    temporal_range #,
    #temporal_slope
  )
  
  writeRaster(
    out,
    out_file,
    overwrite = TRUE,
    wopt = list(
      datatype = "FLT4S",
      gdal = c(
        "COMPRESS=ZSTD",
        "TILED=YES",
        "BIGTIFF=YES"
      )
    )
  )
  
  return(out_file)
}

# Run tiled temporal processing
temporal_files <- lapply(
  seq_along(tiles),
  process_temporal_tile,
  temporal_stack = temporal_stack
)

temporal_files <- Filter(
  Negate(is.null),
  temporal_files
)

temporal_files <- temporal_files[
  file.exists(temporal_files)
]

# Mosaic temporal outputs
cat(
  "Mosaicking",
  length(temporal_files),
  "temporal tiles...\n"
)

temporal_out <- mosaic(
  rast(temporal_files),
  fun = "mean"
)

names(temporal_out) <- c(
  "temporal_sd",
  "temporal_entropy",
  "temporal_range"#,
  #"temporal_slope"
)

# Extract layers
temporal_sd <- temporal_out[["temporal_sd"]]

temporal_entropy <- temporal_out[["temporal_entropy"]]

temporal_range <- temporal_out[["temporal_range"]]

#temporal_slope <- temporal_out[["temporal_slope"]]

#------ Tiled temporal workflow (native vs tame vs nongrass) ------#

# Build annual/period triclass stacks
nt_years <- sort(unique(nt_meta$year))

# nt_year_rasters <- lapply(nt_years, function(y) {
#   
#   idx <- which(nt_meta$year == y)
#   
#   s <- rast(nt_local[idx])
#   
#   if (length(idx) == 1) {
#     
#     r <- s
#     
#   } else {
#     
#     # average class code NOT appropriate
#     # instead build class proportions
#     
#     stop("Need proportion workflow below")
#   }
# })

# Create yearly class frequencies
nt_year_probs <- lapply(nt_years, function(y) {
  
  idx <- which(nt_meta$year == y)
  
  s <- rast(nt_local[idx])
  
  p_not <- app(s == 0, mean, na.rm = TRUE)
  p_tame <- app(s == 1, mean, na.rm = TRUE)
  p_native <- app(s == 2, mean, na.rm = TRUE)
  
  names(p_not) <- paste0("not_", y)
  names(p_tame) <- paste0("tame_", y)
  names(p_native) <- paste0("native_", y)
  
  list(
    p_not = p_not,
    p_tame = p_tame,
    p_native = p_native
  )
})

# Build time series stacks for each class
not_stack <- rast(
  lapply(nt_year_probs, function(x) x$p_not)
)

tame_stack <- rast(
  lapply(nt_year_probs, function(x) x$p_tame)
)

native_stack <- rast(
  lapply(nt_year_probs, function(x) x$p_native)
)

names(not_stack) <- paste0("y", nt_years)
names(tame_stack) <- paste0("y", nt_years)
names(native_stack) <- paste0("y", nt_years)

# Tiled triclass temporal processor
process_triclass_temporal_tile <- function(i) {
  
  tile <- tiles[[i]]
  
  out_file <- file.path(
    temp_tile_dir,
    sprintf("triclass_temporal_%05d.tif", i)
  )
  
  if (file.exists(out_file)) {
    return(out_file)
  }
  
  not_y <- crop(not_stack, tile)
  tame_y <- crop(tame_stack, tile)
  native_y <- crop(native_stack, tile)
  
  if (is.na(global(not_y[[1]], "max")[[1]])) {
    return(NULL)
  }
  
  # -----------------------------------
  # Temporal SD by class
  # -----------------------------------
  
  sd_not <- app(not_y, sd, na.rm = TRUE)
  sd_tame <- app(tame_y, sd, na.rm = TRUE)
  sd_native <- app(native_y, sd, na.rm = TRUE)
  
  names(sd_not) <- "sd_not"
  names(sd_tame) <- "sd_tame"
  names(sd_native) <- "sd_native"
  
  # -----------------------------------
  # Maximum temporal instability
  # -----------------------------------
  
  max_sd <- app(
    c(sd_not, sd_tame, sd_native),
    max,
    na.rm = TRUE
  )
  
  names(max_sd) <- "class_instability"
  
  # -----------------------------------
  # Dominant class over whole time period
  # -----------------------------------
  
  mean_not <- app(not_y, mean, na.rm = TRUE)
  mean_tame <- app(tame_y, mean, na.rm = TRUE)
  mean_native <- app(native_y, mean, na.rm = TRUE)
  
  dominant_class <- as.int(
    c(mean_not, mean_tame, mean_native),
    which.max
  )
  
  names(dominant_class) <- "dominant_class"
  
  # -----------------------------------
  # Temporal class entropy
  # -----------------------------------
  
  eps <- 1e-6
  
  class_entropy <- -(
    mean_not * log(mean_not + eps) +
      mean_tame * log(mean_tame + eps) +
      mean_native * log(mean_native + eps)
  )
  
  class_entropy <- class_entropy / log(3)
  
  names(class_entropy) <- "class_entropy"
  
  out <- c(
    sd_not,
    sd_tame,
    sd_native,
    max_sd,
    dominant_class,
    class_entropy
  )
  
  writeRaster(
    out,
    out_file,
    overwrite = TRUE,
    wopt = list(
      datatype = "FLT4S",
      gdal = c(
        "COMPRESS=ZSTD",
        "TILED=YES",
        "BIGTIFF=YES"
      )
    )
  )
  
  out_file
}

# Run
triclass_temporal_files <- lapply(
  seq_along(tiles),
  process_triclass_temporal_tile
)

triclass_temporal_files <- Filter(
  Negate(is.null),
  triclass_temporal_files
)

triclass_temporal_files <- triclass_temporal_files[
  file.exists(triclass_temporal_files)
]

triclass_temporal_out <- mosaic(
  rast(triclass_temporal_files),
  fun = "first"
)

names(triclass_temporal_out) <- c(
  "sd_not",
  "sd_tame",
  "sd_native",
  "class_instability",
  "dominant_class",
  "class_entropy"
)

# Extract layers
sd_not <- triclass_temporal_out[["sd_not"]]
sd_tame <- triclass_temporal_out[["sd_tame"]]
sd_native <- triclass_temporal_out[["sd_native"]]

class_instability <- triclass_temporal_out[["class_instability"]]

dominant_class <- triclass_temporal_out[["dominant_class"]]
dominant_class <- as.int(dominant_class)

class_entropy <- triclass_temporal_out[["class_entropy"]]

#------ Diagnostics ------#

# summary statistics
# use to refine thresholds especially for support and conf_adj
global(
  c(
    support_score,
    confidence,
    decision_score
  ),
  c("min","max","mean"),
  na.rm = TRUE
)

# calculate coverage after masking
coverage <- global(
  !is.na(hard_class),
  "mean",
  na.rm = TRUE
)

log_msg(
  "Classified coverage: ",
  round(100 * coverage[1,1], 1),
  "%"
)

# check confidence to inform calibration
global(decision_score,
       c("min","max","mean"),
       na.rm = TRUE)

# Visualize confidence vs entropy
# to understand whether masking rules are too strict/loose

# Sample ~1e6 pixels max (adjust if needed)
set.seed(42)

vals <- spatSample(
  c(decision_score, entropy),
  size = 1e6,
  method = "random",
  na.rm = TRUE,
  as.df = TRUE
)

colnames(vals) <- c("decision_score", "entropy")

# Confidence vs entropy plot
# bottom-right = high confidence, low entropy (best pixels)
# top-left = low confidence, high entropy (worst pixels)
# middle cloud = uncertain transitions
# Thresholds appear as: blue line: confidence = 0.6, red line: entropy = 0.8
library(ggplot2)

p1 <- ggplot(vals, aes(x = decision_score, y = entropy)) +
  stat_binhex(bins = 60) +
  scale_fill_viridis_c(trans = "log10") +
  geom_vline(xintercept = 0.4, color = "blue", linetype = "dashed") +
  geom_hline(yintercept = 0.8, color = "red", linetype = "dashed") +
  labs(
    title = "Confidence vs Entropy",
    x = "Confidence (max probability)",
    y = "Entropy",
    fill = "Pixel count"
  ) +
  theme_minimal()

ggsave(
  file.path(diag_dir, "confidence_vs_entropy.png"),
  p1,
  width = 8,
  height = 6,
  dpi = 300
)

# Density plot (smooth view)
p2 <- ggplot(vals, aes(x = decision_score, y = entropy)) +
  geom_density_2d_filled() +
  geom_vline(xintercept = 0.4, color = "blue", linetype = "dashed") +
  geom_hline(yintercept = 0.8, color = "red", linetype = "dashed") +
  labs(
    title = "Density of Confidence vs Entropy"
  ) +
  theme_minimal()

ggsave(
  file.path(diag_dir, "confidence_vs_entropy_density.png"),
  p2,
  width = 8,
  height = 6,
  dpi = 300
)

# Confidence distribution
# If most points cluster around confidence > 0.5 then threshold is appropriate
p3 <- ggplot(vals, aes(x = decision_score)) +
  geom_histogram(bins = 50, fill = "steelblue") +
  geom_vline(xintercept = 0.4, color = "red", linetype = "dashed") +
  labs(title = "Confidence distribution") +
  theme_minimal()

ggsave(
  file.path(diag_dir, "confidence_dist.png"),
  p3,
  width = 8,
  height = 6,
  dpi = 300
)

# Entropy distribution
# If most points cluster around entropy < 0.6 then threshold is appropriate
p4 <- ggplot(vals, aes(x = entropy)) +
  geom_histogram(bins = 50, fill = "darkgreen") +
  geom_vline(xintercept = 0.8, color = "red", linetype = "dashed") +
  labs(title = "Entropy distribution") +
  theme_minimal()

ggsave(
  file.path(diag_dir, "entropy_dist.png"),
  p4,
  width = 8,
  height = 6,
  dpi = 300
)

log_msg("---- FINAL DIAGNOSTICS ----")

log_msg("Probability sum check:")
log_msg(capture.output(
  global(check, c("min","max","mean"), na.rm=TRUE)
))

log_msg("ESS:")
log_msg(capture.output(
  global(ess, c("min","max","mean"), na.rm=TRUE)
))

log_msg("Support score:")
log_msg(capture.output(
  global(support_score, c("min","max","mean"), na.rm=TRUE)
))

log_msg("Confidence adjusted:")
log_msg(capture.output(
  global(decision_score, c("min","max","mean"), na.rm=TRUE)
))

log_msg("Entropy:")
log_msg(capture.output(
  global(entropy, c("min","max","mean"), na.rm=TRUE)
))

log_msg("Coverage:")
log_msg(round(100 * coverage[1,1], 1), "%")

log_msg("Year counts:")
log_msg(capture.output(print(year_counts)))

for (r in list(
  num_bin_sources,
  num_nt_sources,
  sum_weight_present,
  vote_margin,
  weighted_vote_margin
)) {
  log_msg(capture.output(
    global(r, c("min","max","mean"), na.rm=TRUE)
  ))
}

#------ Mapping results ------#

# Map of entropy
png(
  file.path(diag_dir, "entropy_map.png"),
  width = 2400,
  height = 1600,
  res = 300
)

plot(
  entropy,
  main = "Entropy (Uncertainty)",
  col = rev(terrain.colors(100))
)

dev.off()

# Map of confidence
png(
  file.path(diag_dir, "confidence_map.png"),
  width = 2400,
  height = 1600,
  res = 300
)

plot(
  decision_score,
  main = "Confidence (Max Probability)",
  col = viridis::viridis(100)
)

dev.off()

# Binary "uncertain areas" map
uncertain_mask <- ifel(
  (decision_score < 0.4) | (entropy > 0.8),
  1,
  0
)

names(uncertain_mask) <- "uncertain"

png(
  file.path(diag_dir, "binary_uncertain_areas_map.png"),
  width = 2400,
  height = 1600,
  res = 300
)

plot(
  uncertain_mask,
  main = "Uncertain Areas",
  col = c("transparent", "red"),
  legend = FALSE
)

dev.off()

# Overlay uncertainty on classification
png(
  file.path(diag_dir, "binary_uncertain_areas_hardclasss_map.png"),
  width = 2400,
  height = 1600,
  res = 300
)

plot(hard_class, main="Classification with Uncertainty Overlay")

plot(
  uncertain_mask,
  col = c(NA, adjustcolor("red", alpha.f = 0.5)),
  add = TRUE
)

dev.off()

# Continuous uncertainty "hotspot" map
png(
  file.path(diag_dir, "continuous_uncertainty_hotspots_map.png"),
  width = 2400,
  height = 1600,
  res = 300
)

plot(
  uncertainty_total,
  main = "Uncertainty Hotspots",
  col = hcl.colors(100, "Reds")
)

dev.off()

# Mask only uncertain areas
uncertain_only <- mask(entropy, uncertain_mask, maskvalues = 0)

png(
  file.path(diag_dir, "entropy_uncertain_areas_map.png"),
  width = 2400,
  height = 1600,
  res = 300
)

plot(
  uncertain_only,
  main = "Entropy (Only Uncertain Areas)",
  col = rev(heat.colors(100))
)

dev.off()

#------- Write outputs to S3 ----------#
log_msg("Writing outputs")

aws.s3::put_object(
  file = file.path(diag_dir, "support_score_hist.png"),
  object = paste0(out_dir, "/diagnostics/support_score_hist.png"))

aws.s3::put_object(
  file = file.path(diag_dir, "decision_score_hist.png"),
  object = paste0(out_dir, "/diagnostics/decision_score_hist.png"))

aws.s3::put_object(
  file = file.path(diag_dir, "confidence_vs_entropy.png"),
  object = paste0(out_dir, "/diagnostics/confidence_vs_entropy.png"))

aws.s3::put_object(
  file = file.path(diag_dir, "confidence_vs_entropy_density.png"),
  object = paste0(out_dir, "/diagnostics/confidence_vs_entropy_density.png"))

aws.s3::put_object(
  file = file.path(diag_dir, "confidence_dist.png"),
  object = paste0(out_dir, "/diagnostics/confidence_dist.png"))

aws.s3::put_object(
  file = file.path(diag_dir, "entropy_dist.png"),
  object = paste0(out_dir, "/diagnostics/entropy_dist.png"))

aws.s3::put_object(
  file = file.path(diag_dir, "entropy_map.png"),
  object = paste0(out_dir, "/diagnostics/entropy_map.png"))

aws.s3::put_object(
  file = file.path(diag_dir, "confidence_map.png"),
  object = paste0(out_dir, "/diagnostics/confidence_map.png"))

aws.s3::put_object(
  file = file.path(diag_dir, "binary_uncertain_areas_map.png"),
  object = paste0(out_dir, "/diagnostics/binary_uncertain_areas_map.png"))

aws.s3::put_object(
  file = file.path(diag_dir, "binary_uncertain_areas_hardclasss_map.png"),
  object = paste0(out_dir, "/diagnostics/binary_uncertain_areas_hardclasss_map.png"))

aws.s3::put_object(
  file = file.path(diag_dir, "continuous_uncertainty_hotspots_map.png"),
  object = paste0(out_dir, "/diagnostics/continuous_uncertainty_hotspots_map.png"))

aws.s3::put_object(
  file = file.path(diag_dir, "entropy_uncertain_areas_map.png"),
  object = paste0(out_dir, "/diagnostics/entropy_uncertain_areas_map.png"))

aws.s3::put_object(
  file = "tile_index.gpkg",
  object = paste0(out_dir, "/tile_index.gpkg")
)

# Define wopt
wopt       <- list(datatype="FLT4S", gdal=c("COMPRESS=ZSTD","TILED=YES","BIGTIFF=YES"))
wopt_int   <- list(datatype="INT2U", gdal=c("COMPRESS=ZSTD","TILED=YES","BIGTIFF=YES"))
wopt_class <- list(datatype="INT1U", gdal=c("COMPRESS=ZSTD","TILED=YES","BIGTIFF=YES"))

# Write locally first, then save to S3
outputs <- list(
  list(raster = p_native, local = "p_native.tif", s3 = "p_native.tif", wopt = wopt),
  list(raster = p_tame, local = "p_tame.tif", s3 = "p_tame.tif", wopt = wopt),
  list(raster = p_not, local = "p_not_grass.tif", s3 = "p_not_grass.tif", wopt = wopt),
  list(raster = entropy, local = "entropy.tif", s3 = "entropy.tif", wopt = wopt),
  list(raster = entropy_norm, local = "entropy_norm.tif", s3 = "entropy_norm.tif", wopt = wopt),
  list(raster = ess, local = "ess.tif", s3 = "ess.tif", wopt = wopt),
  list(raster = ess_uncertainty, local = "ess_uncertainty.tif", s3 = "ess_uncertainty.tif", wopt = wopt),
  list(raster = support_score, local = "support_score.tif", s3 = "support_score.tif", wopt = wopt),
  list(raster = support_uncertainty, local = "support_uncertainty.tif", s3 = "support_uncertainty.tif", wopt = wopt),
  list(raster = native_tame_missing, local = "native_tame_missing.tif", s3 = "native_tame_missing.tif", wopt = wopt),
  list(raster = confidence, local = "confidence.tif", s3 = "confidence.tif", wopt = wopt),
  list(raster = decision_score, local = "decision_score.tif", s3 = "decision_score.tif", wopt = wopt),
  list(raster = hard_class, local = "hard_class.tif", s3 = "hard_class.tif", wopt = wopt_class),
  list(raster = high_conf_class, local = "high_conf_class.tif", s3 = "high_conf_class.tif", wopt = wopt_class),
  list(raster = uncertainty_total, local = "uncertainty_total.tif", s3 = "uncertainty_total.tif", wopt = wopt),
  list(raster = uncertainty_total_w, local = "uncertainty_total_w.tif", s3 = "uncertainty_total_w.tif", wopt = wopt),
  list(raster = temporal_sd, local = "temporal_sd.tif", s3 = "temporal_sd_bin.tif", wopt = wopt),
  list(raster = temporal_entropy, local = "temporal_entropy.tif", s3 = "temporal_entropy_bin.tif", wopt = wopt),
  list(raster = temporal_range, local = "temporal_range.tif", s3 = "temporal_range_bin.tif", wopt = wopt),
  list(raster = num_bin_sources, local = "num_bin_sources.tif", s3 = "num_bin_sources.tif", wopt = wopt_int),
  list(raster = num_nt_sources, local = "num_nt_sources.tif", s3 = "num_nt_sources.tif", wopt = wopt_int),
  list(raster = num_positive, local = "num_positive.tif", s3 = "num_positive.tif", wopt = wopt_int),
  list(raster = num_negative, local = "num_negative.tif", s3 = "num_negative.tif", wopt = wopt_int),
  list(raster = sum_weight_present, local = "sum_weight_present.tif", s3 = "sum_weight_present.tif", wopt = wopt),
  list(raster = vote_margin, local = "vote_margin.tif", s3 = "vote_margin.tif", wopt = wopt),
  list(raster = weighted_vote_margin, local = "weighted_vote_margin.tif", s3 = "weighted_vote_margin.tif", wopt = wopt),
  list(raster = num_not_votes, local = "num_not_votes.tif", s3 = "num_not_votes.tif", wopt = wopt_int),
  list(raster = num_tame_votes, local = "num_tame_votes.tif", s3 = "num_tame_votes.tif", wopt = wopt_int),
  list(raster = num_native_votes, local = "num_native_votes.tif", s3 = "num_native_votes.tif", wopt = wopt_int),
  list(raster = triclass_agreement, local = "triclass_agreement.tif", s3 = "triclass_agreement.tif", wopt = wopt),
  list(raster = vote_margin_tri, local = "vote_margin_tri.tif", s3 = "vote_margin_tri.tif", wopt = wopt),
  list(raster = weighted_agreement_tri, local = "weighted_agreement_tri.tif", s3 = "weighted_agreement_tri.tif", wopt = wopt),
  list(raster = ratio_not_dataset, local = "ratio_not_dataset.tif", s3 = "ratio_not_dataset.tif", wopt = wopt),
  list(raster = ratio_native_dataset, local = "ratio_native_dataset.tif", s3 = "ratio_native_dataset.tif", wopt = wopt),
  list(raster = ratio_tame_dataset, local = "ratio_tame_dataset.tif", s3 = "ratio_tame_dataset.tif", wopt = wopt),
  list(raster = max_class_ratio, local = "max_class_ratio.tif", s3 = "max_class_ratio.tif", wopt = wopt),
  list(raster = fraction_binary_support, local = "fraction_binary_support.tif", s3 = "fraction_binary_support.tif", wopt = wopt),
  list(raster = fraction_triclass_support, local = "fraction_triclass_support.tif", s3 = "fraction_triclass_support.tif", wopt = wopt),
  list(raster = sd_not, local = "temporal_sd_not.tif", s3 = "temporal_sd_not.tif", wopt = wopt),
  list(raster = sd_tame, local = "temporal_sd_tame.tif", s3 = "temporal_sd_tame.tif", wopt = wopt),
  list(raster = sd_native, local = "temporal_sd_native.tif", s3 = "temporal_sd_native.tif", wopt = wopt),
  list(raster = class_instability, local = "class_instability.tif", s3 = "class_instability_triclass.tif", wopt = wopt),
  list(raster = dominant_class, local = "dominant_temporal_class.tif", s3 = "dominant_temporal_class_triclass.tif", wopt = wopt_class),
  list(raster = class_entropy, local = "class_entropy.tif", s3 = "class_entropy_triclass.tif", wopt = wopt)
)

write_and_collect <- function(spec) {
  local_path <- file.path(tmp_root, spec$local)
  writeRaster(spec$raster, local_path, overwrite = TRUE, wopt = spec$wopt)
  spec$local_path <- local_path
  spec
}

outputs <- lapply(outputs, write_and_collect)
log_msg("All output rasters written locally; beginning S3 upload")

rm(
  bin_stack,
  nt_stack
)

gc()

upload_results <- lapply(outputs, function(spec) {
  s3_object <- paste0(out_dir, "/", spec$s3)
  ok <- tryCatch({
    aws.s3::put_object(file = spec$local_path, object = s3_object)
    TRUE
  }, error = function(e) {
    log_msg("FAILED to upload ", spec$s3, ": ", conditionMessage(e), level = "ERROR")
    FALSE
  })
  if (ok) log_msg("Uploaded: ", s3_object)
  data.frame(file = spec$s3, uploaded = ok)
})

upload_summary <- do.call(rbind, upload_results)
n_failed <- sum(!upload_summary$uploaded)

if (n_failed > 0) {
  log_msg("WARNING: ", n_failed, " of ", nrow(upload_summary),
          " output uploads failed: ",
          paste(upload_summary$file[!upload_summary$uploaded], collapse = ", "),
          level = "ERROR")
} else {
  log_msg("All ", nrow(upload_summary), " output rasters uploaded successfully")
}

log_msg("Analysis complete")

# Upload log to S3
aws.s3::put_object(
  file = log_file_local,
  object = paste0(
    "alison-long/tomorrows-prairies/logs/",
    basename(log_file_local)
  ),
  bucket = "sagemaker-gst-stage.sharing"
)
