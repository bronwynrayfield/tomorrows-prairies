
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
dir.create(local_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(stack_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(terra_dir, recursive = TRUE, showWarnings = FALSE)

meta_file  <- "s3://sagemaker-gst-stage.sharing/alison-long/tomorrows-prairies/processed/metadata_full_extent_binary_nattam_with_path_nohfi.csv"
ecodist_file <- "s3://sagemaker-gst-stage.sharing/alison-long/tomorrows-prairies/processed/boundaries/national-ecological-framework-for-canada/Ecodistricts/ecodistricts_clean_rast_Oct2026.tif"
pilot_file   <- "s3://sagemaker-gst-stage.sharing/alison-long/tomorrows-prairies/processed/pilot-area/pilot_area_aligned_cnvc_Jan2026.tif"
sa_file      <- "s3://sagemaker-gst-stage.sharing/alison-long/tomorrows-prairies/processed/study-area/modelling_extent_aligned_cnvc_Jan2026.tif"
acc_bin_file <- "s3://sagemaker-gst-stage.sharing/alison-long/tomorrows-prairies/processed/accuracy_metrics/binary_accuracy_by_ecodistrict_atomic_combinedCIandCM_23Jun2026.csv"
acc_nt_file  <- "s3://sagemaker-gst-stage.sharing/alison-long/tomorrows-prairies/processed/accuracy_metrics/nattam_accuracy_by_ecodistrict_atomic_lenient_23Jun2026.csv"

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
global(
  acc_bin_stack,
  c("min","max"),
  na.rm = TRUE
)

global(
  metric_tame_stack,
  c("min","max"),
  na.rm = TRUE
)

global(
  metric_native_stack,
  c("min","max"),
  na.rm = TRUE
)

# =========================
# TILE GRID
# =========================
# ~2000–5000 tiles is ideal, 500-1500 for single node
tiles <- makeTiles(sa, n = 1000, overlap = 0)

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
  
  present <- !is.na(b)
  
  bin_w_eff <- rbin # * pbin
  
  bin_w_eff_present <- bin_w_eff * present
  
  #for (j in seq_len(nlyr(bin_w_eff))) {
  #  bin_w_eff[[j]] <- bin_w_eff[[j]] * bin_time[j]
  #}
  
  nt_present <- !is.na(nt)
  
  nt_w_tame <- rt # * pnt
  nt_w_native <- rn # * pnt
  
  nt_w_tame_present   <- nt_w_tame   * nt_present
  nt_w_native_present <- nt_w_native * nt_present
  
  nt_not_present <- rnot * nt_present
  
  #for (j in seq_len(nlyr(nt_w_tame))) {
  #  nt_w_tame[[j]] <- nt_w_tame[[j]] * nt_time[j]
  #  nt_w_native[[j]] <- nt_w_native[[j]] * nt_time[j]
  #}
  
  #NEG_HALF <- 0.5
  
  # ----------------------
  # EFFECTIVE SAMPLE SIZE (effective weighted support - does not use redundancy adjustment)
  # ----------------------
  
  # Sum of weights - use only binary rasters so those with native vs tame don't get double counted
  w_sum <- app(
    bin_w_eff_present,
    sum,
    na.rm = TRUE
  )
  
  # Sum of squared weights
  w_sum_sq <- app(
    bin_w_eff_present^2,
    sum,
    na.rm = TRUE
  )
  
  # Effective sample size
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
  # EVIDENCE COMPUTATION
  # ----------------------
  
  # ----------------------
  # NORMALIZED WEIGHTS
  # ----------------------
  
  # binary weights
  #bin_w_eff <- rbin * pbin
  bin_w_sum <- app(bin_w_eff_present, sum, na.rm=TRUE)
  bin_w_norm <- ifel(bin_w_sum > 0, bin_w_eff_present / bin_w_sum, 0)
  
  # tri-class weights
  #nt_w_tame   <- rt * pnt
  #nt_w_native <- rn * pnt
  
  nt_w_sum <- app(nt_w_tame_present + nt_w_native_present, sum, na.rm=TRUE)
  nt_w_tame_norm   <- ifel(nt_w_sum > 0, nt_w_tame_present / nt_w_sum, 0)
  nt_w_native_norm <- ifel(nt_w_sum > 0, nt_w_native_present / nt_w_sum, 0)
  
  # ----------------------
  # HIERARCHICAL EVIDENCE
  # ----------------------
  
  # Binary evidence
  E_bin_grass <- app(
    bin_w_norm * (b == 1),
    sum,
    na.rm = TRUE
  )
  
  E_bin_not <- app(
    bin_w_norm * (b == 0),
    sum,
    na.rm = TRUE
  )
  
  # ----------------------
  # Triclass class-specific evidence
  # ----------------------
  
  #nt_not_present <- metric_not_stack * nt_present
  
  nt_not_sum <- app(
    nt_not_present,
    sum,
    na.rm = TRUE
  )
  
  nt_not_norm <- ifel(
    nt_not_sum > 0,
    nt_not_present / nt_not_sum,
    0
  )
  
  E_not_tri <- app(
    nt_not_norm * (nt == 0),
    sum,
    na.rm = TRUE
  )
  
  E_tame_tri <- app(
    nt_w_tame_norm * (nt == 1),
    sum,
    na.rm = TRUE
  )
  
  E_native_tri <- app(
    nt_w_native_norm * (nt == 2),
    sum,
    na.rm = TRUE
  )
  
  # ----------------------
  # Grass vs nongrass
  # ----------------------
  
  Grass_support <-
    E_bin_grass +
    E_tame_tri +
    E_native_tri
  
  Not_support <-
    E_bin_not +
    E_not_tri
  
  grass_total <- Grass_support + Not_support
  
  p_grass <- Grass_support / grass_total
  p_grass <- ifel(grass_total > 0, p_grass, NA)
  
  p_not <- Not_support / grass_total
  p_not <- ifel(grass_total > 0, p_not, NA)
  
  # ----------------------
  # Native vs tame within grass
  # ----------------------
  
  grass_type_total <-
    E_tame_tri +
    E_native_tri
  
  # Pixels with grass evidence but no native/tame evidence
  #native_tame_missing <- ifel(grass_type_total == 0, 1, 0)
  native_tame_missing <- ifel((grass_type_total == 0) & (p_grass > 0.5),
    1, 0)
  names(native_tame_missing) <- "native_tame_missing"
  
  # Shares
  share_native <-
    E_native_tri / grass_type_total
  
  share_tame <-
    E_tame_tri / grass_type_total
  
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
  
  p_not    <- p_not    / prob_sum
  p_tame   <- p_tame   / prob_sum
  p_native <- p_native / prob_sum
  
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
           support_uncertainty, native_tame_missing)
  names(out) <- c("p_not", "p_tame", "p_native", "entropy", "entropy_norm",
                  "ess", "ess_uncertainty", "support_score", "support_uncertainty",
                  "native_tame_missing")
  
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
stopifnot(all(vapply(files, function(f) nlyr(rast(f)) == 10, logical(1))))

# Mosaic all tiles together
mosaic_out <- mosaic(rast(files), fun="first")
# or can try
#mosaic_out <- do.call(mosaic, c(lapply(files, rast), list(fun="mean")))

# Set names
names(mosaic_out) <- c("p_not", "p_tame", "p_native", "entropy", "entropy_norm",
                       "ess", "ess_uncertainty", "support_score", "support_uncertainty",
                       "native_tame_missing")

#------ Derived layers ------#
# Probability stack
prob_stack <- mosaic_out[[c("p_not", "p_tame", "p_native")]]
entropy <- mosaic_out[["entropy"]]

# support score
support_score <- mosaic_out[["support_score"]]
support_uncertainty <- mosaic_out[["support_uncertainty"]]

# Missing
native_tame_missing <- mosaic_out[["native_tame_missing"]]

# Check
hist(values(support_score))

quantile(values(spatSample(support_score, 1e6, na.rm = TRUE)))

# Confidence = max probability
confidence <- app(prob_stack, max)
names(confidence) <- "confidence"

# Check probabilities sum to ~1
check <- app(prob_stack, sum)
stats <- global(check, c("min", "max", "mean"), na.rm=TRUE)
print(stats)
log_msg(capture.output(stats))

# check drivers of masking
global(
  entropy,
  c("min","max","mean"),
  na.rm = TRUE)

global(
  support_score,
  c("min","max","mean"),
  na.rm = TRUE)

# Convert confidence to uncertainty
confidence_uncertainty <- 1 - confidence

names(confidence_uncertainty) <- "confidence_uncertainty"

# adjusted confidence
confidence <- app(prob_stack, max)

confidence_adj <- confidence * support_score
names(confidence_adj) <- "confidence_adjusted"

entropy_norm <- mosaic_out[["entropy_norm"]]
ess_uncertainty <- mosaic_out[["ess_uncertainty"]]

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
  confidence_adj,
  size = 1e6,
  na.rm = TRUE,
  values = TRUE)

hist(vals_conf)

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
  (confidence_adj < 0.4) | (entropy > 0.8),
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
  (confidence_adj >= 0.75) & (entropy <= 0.6),
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
bin_meta <- bin_meta %>%
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
periods <- unique(bin_meta$period)

period_rasters <- lapply(periods, function(p) {
  
  idx <- which(bin_meta$period == p)
  
  app(
    rast(bin_local[idx]),
    mean,
    na.rm = TRUE
  )
})

year_stack_pds <- rast(period_rasters)

# Create temporary tile directory
temp_tile_dir <- "/mnt/tmp/temporal_tiles"
dir.create(temp_tile_dir, showWarnings = FALSE)

# Tile processor
process_temporal_tile <- function(i) {
  
  tile <- tiles[[i]]
  
  out_file <- file.path(
    temp_tile_dir,
    sprintf("temporal_%05d.tif", i)
  )
  
  if (file.exists(out_file)) {
    return(out_file)
  }
  
  ys <- crop(year_stack, tile)
  
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
  process_temporal_tile
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

#------ Temporal analysis (grass vs nongrass) ------#

# Mean grassland probability by year
#years <- sort(unique(bin_meta$year))

#year_rasters <- lapply(years, function(y) {

#  idx <- which(bin_meta$year == y)

#  if (length(idx) == 1) {
#    rast(bin_local[idx])
#  } else {
#    app(
#      rast(bin_local[idx]),
#      mean,
#      na.rm = TRUE
#    )
#  }
#})

#year_stack <- rast(year_rasters)

#names(year_stack) <- paste0("y", years)

# Temporal variance
#temporal_var <- app(
#  year_stack,
#  var,
#  na.rm = TRUE
#)

#names(temporal_var) <- "temporal_variance"

# Temporal SD (easier to interpret)
# ~0: Stable through time
# 0.1–0.3: Some change
# >0.4: Strong temporal disagreement/change
# ~0.5: Frequent switching
#temporal_sd <- sqrt(temporal_var)

#names(temporal_sd) <- "temporal_sd"

# highly dependent on uneven sampling by year
# may want to require a minimum number of rasters per year
# or aggregate into broader periods
# check:
#table(bin_meta$year)

#---- Temporal entropy ----#

# Frequency of grass occurrence
#grass_freq <- app(
#  year_stack,
#  mean,
#  na.rm = TRUE
#)

# Temporal entropy
# 0: Always grass or always non-grass
# Moderate: Occasional change
# High: Frequently changing through time
#eps <- 1e-6

#temporal_entropy <- -(
#  grass_freq * log(grass_freq + eps) +
#    (1 - grass_freq) * log(1 - grass_freq + eps)
#)

#names(temporal_entropy) <- "temporal_entropy"

#---- Temporal range ----#
# Measure difference between the maximum and minimum grassland probability over time
# 0 = no temporal change
# 1 = switched completely between grass and non-grass
#temporal_range <- app(
#  year_stack,
#  function(x) max(x, na.rm=TRUE) - min(x, na.rm=TRUE)
#)

#names(temporal_range) <- "temporal_range"

#---- Temporal trend layer ----#
# positive = becoming more likely grassland
# negative = becoming less likely grassland
# near zero = stable
#yrs <- years

#temporal_slope <- app(
#  year_stack,
#  function(x) {
#    if (all(is.na(x))) return(NA)
#    if(sum(!is.na(x)) < 3) return(NA) # minimum to increase stability 
#    coef(lm(x ~ yrs))[2]
#  }
#)

#names(temporal_slope) <- "temporal_slope"

#------ Diagnostics ------#

# summary statistics
# use to refine thresholds especially for support and conf_adj
global(
  c(
    support_score,
    confidence,
    confidence_adj
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
global(confidence_adj,
       c("min","max","mean"),
       na.rm = TRUE)

# Visualize confidence vs entropy
# to understand whether masking rules are too strict/loose

# Sample ~1e6 pixels max (adjust if needed)
set.seed(42)

vals <- spatSample(
  c(confidence_adj, entropy),
  size = 1e6,
  method = "random",
  na.rm = TRUE,
  as.df = TRUE
)

colnames(vals) <- c("confidence_adj", "entropy")

# Confidence vs entropy plot
# bottom-right = high confidence, low entropy (best pixels)
# top-left = low confidence, high entropy (worst pixels)
# middle cloud = uncertain transitions
# Thresholds appear as: blue line: confidence = 0.6, red line: entropy = 0.8
library(ggplot2)

ggplot(vals, aes(x = confidence_adj, y = entropy)) +
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

# Density plot (smooth view)
ggplot(vals, aes(x = confidence_adj, y = entropy)) +
  geom_density_2d_filled() +
  geom_vline(xintercept = 0.4, color = "blue", linetype = "dashed") +
  geom_hline(yintercept = 0.8, color = "red", linetype = "dashed") +
  labs(
    title = "Density of Confidence vs Entropy"
  ) +
  theme_minimal()

# Confidence distribution
# If most points cluster around confidence > 0.5 then threshold is appropriate
ggplot(vals, aes(x = confidence_adj)) +
  geom_histogram(bins = 50, fill = "steelblue") +
  geom_vline(xintercept = 0.4, color = "red", linetype = "dashed") +
  labs(title = "Confidence distribution") +
  theme_minimal()

# Entropy distribution
# If most points cluster around entropy < 0.6 then threshold is appropriate
ggplot(vals, aes(x = entropy)) +
  geom_histogram(bins = 50, fill = "darkgreen") +
  geom_vline(xintercept = 0.8, color = "red", linetype = "dashed") +
  labs(title = "Entropy distribution") +
  theme_minimal()

#------ Mapping results ------#

# Map of entropy
plot(
  entropy,
  main = "Entropy (Uncertainty)",
  col = rev(terrain.colors(100))
)


# Map of confidence
plot(
  confidence_adj,
  main = "Confidence (Max Probability)",
  col = viridis::viridis(100)
)

# Binary "uncertain areas" map
uncertain_mask <- ifel(
  (confidence_adj < 0.4) | (entropy > 0.8),
  1,
  0
)

names(uncertain_mask) <- "uncertain"

plot(
  uncertain_mask,
  main = "Uncertain Areas",
  col = c("transparent", "red"),
  legend = FALSE
)

# Overlay uncertainty on classification
plot(hard_class, main="Classification with Uncertainty Overlay")

plot(
  uncertain_mask,
  col = c(NA, adjustcolor("red", alpha.f = 0.5)),
  add = TRUE
)

# Continuous uncertainty "hotspot" map
plot(
  uncertainty_total,
  main = "Uncertainty Hotspots",
  col = hcl.colors(100, "Reds")
)

# Mask only uncertain areas
uncertain_only <- mask(entropy, uncertain_mask, maskvalues = 0)

plot(
  uncertain_only,
  main = "Entropy (Only Uncertain Areas)",
  col = rev(heat.colors(100))
)


#------ Write outputs ------#

ess <- mosaic_out[["ess"]]
ess_uncertainty <- mosaic_out[["ess_uncertainty"]]

# Combined raster 
final_stack <- c(mosaic_out, 
                 confidence, confidence_adj, hard_class, high_conf_class,
                 uncertainty_total, uncertainty_total_w,
                 temporal_sd, temporal_entropy, temporal_range) #, temporal_slope)

names(final_stack) <- c(
  "p_not",
  "p_tame",
  "p_native",
  "entropy",
  "entropy_norm",
  "ess",
  "ess_uncertainty",
  "support_score",
  "support_uncertainty",
  "native_tame_missing",
  "confidence",
  "confidence_adj",
  "class",
  "high_conf_class",
  "uncertainty_total",
  "uncertainty_total_w",
  "temporal_sd",
  "temporal_entropy",
  "temporal_range"#,
  #  "temporal_slope"
)

stopifnot(nlyr(final_stack) == length(names(final_stack)))

final_out <- "/mnt/tmp/joint_full_stack.tif"

writeRaster(
  final_stack,
  final_out,
  overwrite=TRUE,
  wopt=list(
    datatype="FLT4S",   # classification will still store fine
    gdal=c("COMPRESS=ZSTD","TILED=YES","BIGTIFF=YES")
  )
)

aws.s3::put_object(
  file   = final_out,
  object = paste0(out_dir, "/joint_full_stack.tif")
)

# upload log to S3
aws.s3::put_object(
  file = log_file_local,
  object = log_file
)

# -----------------------------------
# Optional: separate outputs
# -----------------------------------

# Probability raster
prob_out <- "/mnt/tmp/joint_prob_full.tif"

writeRaster(
  mosaic_out,
  prob_out,
  overwrite=TRUE,
  wopt=list(
    datatype="FLT4S",
    gdal=c("COMPRESS=ZSTD","TILED=YES","BIGTIFF=YES")
  )
)

aws.s3::put_object(
  file   = prob_out,
  object = paste0(out_dir, "/joint_prob_full.tif")
)

# Confidence raster
conf_out <- "/mnt/tmp/joint_confidence.tif"

writeRaster(
  confidence,
  conf_out,
  overwrite=TRUE,
  wopt=list(
    datatype="FLT4S",
    gdal=c("COMPRESS=ZSTD","TILED=YES","BIGTIFF=YES")
  )
)

aws.s3::put_object(
  file   = conf_out,
  object = paste0(out_dir, "/joint_confidence.tif")
)

# Classification raster
hard_out <- "/mnt/tmp/joint_hard_class.tif"

writeRaster(
  hard_class,
  hard_out,
  overwrite=TRUE,
  wopt=list(
    datatype="INT1U",
    gdal=c("COMPRESS=ZSTD","TILED=YES","BIGTIFF=YES")
  )
)

aws.s3::put_object(
  file   = hard_out,
  object = paste0(out_dir, "/joint_hard_class.tif")
)

# High-confidence raster
high_conf_out <- "/mnt/tmp/joint_high_conf_class.tif"

writeRaster(
  high_conf_class,
  high_conf_out,
  overwrite = TRUE,
  wopt = list(
    datatype = "INT1U",
    gdal = c("COMPRESS=ZSTD","TILED=YES","BIGTIFF=YES")
  )
)

aws.s3::put_object(
  file   = high_conf_out,
  object = paste0(out_dir, "/joint_high_conf_class.tif")
)

# Temporal SD
writeRaster(
  temporal_sd,
  "/mnt/tmp/temporal_sd.tif",
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

aws.s3::put_object(
  file = "/mnt/tmp/temporal_sd.tif",
  object = paste0(out_dir, "/temporal_sd.tif")
)

# Temporal entropy
writeRaster(
  temporal_entropy,
  "/mnt/tmp/temporal_entropy.tif",
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

aws.s3::put_object(
  file = "/mnt/tmp/temporal_entropy.tif",
  object = paste0(out_dir, "/temporal_entropy.tif")
)

# Temporal range
writeRaster(
  temporal_range,
  "/mnt/tmp/temporal_range.tif",
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

aws.s3::put_object(
  file = "/mnt/tmp/temporal_range.tif",
  object = paste0(out_dir, "/temporal_range.tif")
)

#--------- Possible next steps (examples) ---------------

# Thresholding - high confidence native grass
p_native > 0.8 & entropy < 0.3