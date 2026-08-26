
###################################################################################
##                     GRASSLAND CONSENSUS MAPPING ANALYSIS                      ##
###################################################################################

##----------------------- Alison Long - August 2026 -------------------------

# step 9: use data fusion techniques to create a consensus grassland map 
#         to serve as baseline for modelling

# THIS VERSION IS OPTIMIZED FOR SAGEMAKER
# Next optimization to try: build stacks directly from S3 using GDAL VRT (zero downloads)

# RAN PILOT USING ml.m5d.4xlarge
# trying sa with ml.r5.16xlarge  (64 vCPUs, 512 GiB RAM) - hung after creating tiles?
# trying sa with ml.r5.24xlarge
# could also try ml.m5d.16xlarge (64 vCPUs, 256 GiB RAM, Local NVMe SSD)

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
remove.packages("terra")
unlink(
  "~/R/x86_64-pc-linux-gnu-library/4.4/terra",
  recursive = TRUE,
  force = TRUE
)
file.exists("~/R/x86_64-pc-linux-gnu-library/4.4/terra")
unlink(tempdir(), recursive = TRUE)
remotes::install_version(
  "terra",
  version = "1.9-1",
  upgrade = "never"
)
library(terra)

remove.packages("readr")
unlink(
  "~/R/x86_64-pc-linux-gnu-library/4.4/readr",
  recursive = TRUE,
  force = TRUE
)
file.exists("~/R/x86_64-pc-linux-gnu-library/4.4/readr")
remotes::install_version(
  "readr",
  version = "2.2.0",
  upgrade = "never"
)
library(readr)

system("ldd --version")

install.packages(
  "readr",
  type = "source"
)

remove.packages("tidyverse")
unlink(
  "~/R/x86_64-pc-linux-gnu-library/4.4/tidyverse",
  recursive = TRUE,
  force = TRUE
)
file.exists("~/R/x86_64-pc-linux-gnu-library/4.4/tidyverse")

install.packages(
  "tidyverse",
  type = "source"
)
library(tidyverse)

remove.packages("ggplot2")
remove.packages("scales")

install.packages(
  "ggplot2",
  type = "source"
)
install.packages(
  "scales",
  type = "source"
)

getOption("repos")
options("pkgType")

.libPaths()

find.package("farver")
packageDescription("farver")[c("Package","Version","Built")]

sessionInfo()

system("g++ --version")
system("strings /lib/x86_64-linux-gnu/libstdc++.so.6 | grep GLIBCXX | tail")

dir.exists("~/R/x86_64-pc-linux-gnu-library/4.4")
file.rename(
  "~/R/x86_64-pc-linux-gnu-library/4.4",
  "~/R/x86_64-pc-linux-gnu-library/4.4_backup"
)

list.files("~/R/x86_64-pc-linux-gnu-library")

dir.create(
  "~/R/x86_64-pc-linux-gnu-library/4.4",
  recursive = TRUE,
  showWarnings = FALSE
)

installed.packages()
.libPaths()

list.files("~/R/x86_64-pc-linux-gnu-library/4.4")

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
acc_bin_file <- "s3://sagemaker-gst-stage.sharing/alison-long/tomorrows-prairies/processed/accuracy_metrics/binary_accuracy_atomic_combinedCIandCM_forconsensus_28Jul2026.csv"
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
target_year <- 2025 # year represented by modelling baseline

# sensitivity analysis for lambda
# temporal decay scale in years; larger values retain influence from older datasets longer
# with lambda = 6, weight by years away from 2025:
# 0 years = 1, 1 year = 0.846, 2 years = 0.717, 3 years = 0.607, etc.
lambda_values <- c(5, 6, 7, 10)
yrs_l <- 2015:2025

png(file.path(diag_dir, "lambda_sensitivity.png"),
    width = 2000,
    height = 1500,
    res = 300)

plot(yrs_l, exp(-(2025-yrs_l)/5),
     type = "l",
     ylim = c(0,1),
     xlab = "Dataset year",
     ylab = "Temporal weight")

lines(yrs_l, exp(-(2025-yrs_l)/6), col = 2)
lines(yrs_l, exp(-(2025-yrs_l)/7), col = 3)
lines(yrs_l, exp(-(2025-yrs_l)/10), col = 4)

abline(h = 0.5, lty = 2)

legend("bottomright",
       legend = c("lambda=5","lambda=6","lambda=7","lambda=10"),
       col = 1:4,
       lwd = 2)

dev.off()

# also plot years before baseline
age <- 0:10

png(file.path(diag_dir, "lambda_sensitivity2.png"),
    width = 2000,
    height = 1500,
    res = 300)

plot(age, exp(-age/5),
     type="l",
     ylim=c(0,1),
     xlab="Years before 2025",
     ylab="Weight")

lines(age, exp(-age/6), col=2)
lines(age, exp(-age/7), col=3)
lines(age, exp(-age/10), col=4)

dev.off()

# export
aws.s3::put_object(
  file = file.path(diag_dir, "lambda_sensitivity.png"),
  object = paste0(out_dir, "/diagnostics/lambda_sensitivity.png"))

aws.s3::put_object(
  file = file.path(diag_dir, "lambda_sensitivity2.png"),
  object = paste0(out_dir, "/diagnostics/lambda_sensitivity2.png"))

# choose lambda
lambda <- 6  

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
      # r <- rast(raw_file)
      # 
      # r <- crop(r, pilot)
      # r <- mask(r, pilot)
      
      # crop + mask to study area
      r <- rast(raw_file)

      r <- crop(r, sa)
      r <- mask(r, sa)
      
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

# ----------------------------------
# Remove rasters with no pilot area/study area data
# ----------------------------------

has_data <- function(f) {
  r <- rast(f)
  global(!is.na(r), "sum", na.rm = TRUE)[1,1] > 0
}

bin_keep <- vapply(bin_local, has_data, logical(1))
nt_keep  <- vapply(nt_local, has_data, logical(1))

if (any(!bin_keep)) {
  log_msg(
    "Removing empty binary rasters: ",
    paste(names(bin_local)[!bin_keep], collapse = ", ")
  )
}

if (any(!nt_keep)) {
  log_msg(
    "Removing empty tri-class rasters: ",
    paste(names(nt_local)[!nt_keep], collapse = ", ")
  )
}

bin_local <- bin_local[bin_keep]
nt_local  <- nt_local[nt_keep]

bin_meta <- bin_meta %>%
  filter(raster_id %in% names(bin_local))

nt_meta <- nt_meta %>%
  filter(raster_id %in% names(nt_local))

log_msg("Binary rasters retained: ", length(bin_local))
log_msg("Tri-class rasters retained: ", length(nt_local))

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

req_cols_bin <- c(
  "raster_name",
  "Ecds_ID",
  "acc"
)

if (!all(req_cols_bin %in% names(acc_bin))) {
  stop(
    "Missing required columns in binary accuracy table: ",
    paste(setdiff(req_cols_bin, names(acc_bin)), collapse = ", ")
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

extra_bin_names <- setdiff(
  unique(acc_bin$raster_name),
  meta$raster_id
)

if (length(extra_bin_names) > 0) {
  
  warning(
    "Accuracy table contains unused rasters: ",
    paste(extra_bin_names, collapse = ", ")
  )
  
  acc_bin <- acc_bin %>%
    filter(!raster_name %in% extra_bin_names)
  
  log_msg(
    "Removed ",
    length(extra_bin_names),
    " unused raster(s) from acc_bin"
  )
  
}

## ---- Validate tri-class accuracy table ---- ##
log_msg("Tri-class accuracy rows loaded: ", nrow(acc_nt))

req_cols_nt <- c(
  "raster_name",
  "Ecds_ID",
  #"acc",
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

extra_nt_names <- setdiff(
  unique(acc_nt$raster_name),
  meta$raster_id
)

if (length(extra_nt_names) > 0) {
  
  warning(
    "Accuracy table contains unused rasters: ",
    paste(extra_nt_names, collapse = ", ")
  )
  
  acc_nt <- acc_nt %>%
    filter(!raster_name %in% extra_nt_names)
  
  log_msg(
    "Removed ",
    length(extra_nt_names),
    " unused raster(s) from acc_nt"
  )
  
}

## ---- Add temporal weighting ---- ##
log_msg("Applying accuracy + temporal weighting")

if (any(acc_bin$acc < 0 | acc_bin$acc > 1, na.rm = TRUE)) {
  stop("acc_bin contains values outside [0,1]")
}

for (col in c("native_metric", "tame_metric", "not_metric")) {
  
  if (any(acc_nt[[col]] < 0 | acc_nt[[col]] > 1, na.rm = TRUE)) {
    
    stop(
      sprintf("acc_nt$%s contains values outside [0,1]", col)
    )
    
  }
  
}

acc_to_weight <- function(x) {
  w <- x
  w[!is.finite(w)] <- NA
  pmax(w, 1e-4)
}

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
    w_time = time_decay(
      year,
      target_year,
      lambda
    ),
    
    W_native = native_metric * w_time,
    W_tame   = tame_metric * w_time,
    W_not    = not_metric   * w_time
  )

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
# log_msg("Loading pilot-clipped ecodistrict raster")
# 
# ecodist_r <- rast(ecodist_file)
# 
# ecodist_r <- crop(ecodist_r, pilot)
# ecodist_r <- mask(ecodist_r, pilot)

log_msg("Loading sa-clipped ecodistrict raster")

ecodist_r <- rast(ecodist_file)

ecodist_r <- crop(ecodist_r, sa)
ecodist_r <- mask(ecodist_r, sa)

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

# How much area is represented by missing ecodistricts?
ecodist_freq <- freq(ecodist_r)

ecodist_freq %>%
  filter(value %in% missing_ids_bin) %>%
  arrange(desc(count))

acc_bin %>%
  filter(Ecds_ID %in% missing_ids_bin) %>%
  count(Ecds_ID)

acc_bin %>%
  group_by(Ecds_ID) %>%
  summarise(n_rasters = n_distinct(raster_name))

# fill with raster-level global
acc_global_bin <- acc_bin %>%
  group_by(raster_name) %>%
  summarise(
    global_acc = mean(acc, na.rm = TRUE),
    .groups = "drop"
  )

acc_global_nt <- acc_nt %>%
  group_by(raster_name) %>%
  summarise(
    native_metric = mean(native_metric, na.rm = TRUE),
    tame_metric   = mean(tame_metric, na.rm = TRUE),
    not_metric    = mean(not_metric, na.rm = TRUE),
    .groups = "drop"
  )

# Fill in missing ecodistrict-level accuracy (binary)
if (length(missing_ids_bin) > 0) {
  
  log_msg(
    "Adding missing ecodistricts using raster-level mean accuracy for: ",
    paste(missing_ids_bin, collapse = ", ")
  )
  
  # Create a row for every raster × missing ecodistrict combination
  missing_rows_bin <- tidyr::expand_grid(
    raster_name = unique(acc_bin$raster_name),
    Ecds_ID = missing_ids_bin
  ) %>%
    left_join(acc_global_bin, by = "raster_name") %>%
    left_join(
      meta %>% select(raster_id, year),
      by = c("raster_name" = "raster_id")
    ) %>%
    mutate(
      acc = global_acc,
      accuracy_imputed = TRUE,
      w_acc = acc_to_weight(acc),
      w_time = time_decay(year, target_year, lambda),
      weight = w_acc * w_time
    ) %>%
    select(-global_acc)
  
  # Append to accuracy table
  acc_bin <- acc_bin %>%
    mutate(accuracy_imputed = FALSE) %>%
    bind_rows(missing_rows_bin)
}

log_msg("Added ", nrow(missing_rows_bin), " imputed binary accuracy records")

remaining_missing_bin <- setdiff(ids_r, unique(acc_bin$Ecds_ID))

if (length(remaining_missing_bin) > 0) {
  stop(
    "Still missing ecodistrict IDs after filling: ",
    paste(remaining_missing_bin, collapse = ", ")
  )
}

#---------- Fill in missing ecodistrict-level accuracy (triclass) --------#
if (length(missing_ids_nt) > 0) {
  
  log_msg(
    "Adding missing ecodistricts using raster-level mean accuracy for: ",
    paste(missing_ids_nt, collapse = ", ")
  )
  
  # Create a row for every raster × missing ecodistrict combination
  missing_rows_nt <- tidyr::expand_grid(
    raster_name = unique(acc_nt$raster_name),
    Ecds_ID = missing_ids_nt
  ) %>%
    left_join(acc_global_nt, by = "raster_name") %>%
    left_join(
      meta %>% select(raster_id, year),
      by = c("raster_name" = "raster_id")
    ) %>%
    mutate(
      accuracy_imputed = TRUE,
      w_time = time_decay(year, target_year, lambda),
      
      W_native = native_metric * w_time,
      W_tame   = tame_metric   * w_time,
      W_not    = not_metric    * w_time
    )
  
  # Append to accuracy table
  acc_nt <- acc_nt %>%
    mutate(accuracy_imputed = FALSE) %>%
    bind_rows(missing_rows_nt)
}

log_msg("Added ", nrow(missing_rows_nt), " imputed triclass accuracy records")

remaining_missing_nt <- setdiff(ids_r, unique(acc_nt$Ecds_ID))

if (length(remaining_missing_nt) > 0) {
  stop(
    "Still missing triclass ecodistrict IDs after filling: ",
    paste(remaining_missing_nt, collapse = ", ")
  )
}

# Check that every raster has an accuracy value for every ecodistrict
length(unique(acc_nt$raster_name))
acc_nt %>%
  count(Ecds_ID) %>%
  filter(n < length(unique(acc_nt$raster_name)))

## ---- Build accuracy weight stack helper ---- ##

build_weight_stack <- function(
    accuracy_table,
    stack_names,
    value_col,
    ecodist_r,
    out_prefix
) {
  
  all_ids <- unique(freq(ecodist_r)[, "value"])
  all_ids <- all_ids[!is.na(all_ids)]
  
  weight_files <- lapply(stack_names, function(nm) {
    
    log_msg("Building weight raster: ", nm)
    
    lookup <- accuracy_table %>%
      filter(raster_name == nm) %>%
      select(Ecds_ID, all_of(value_col)) %>%
      as.matrix()
    
    if (nrow(lookup) == 0)
      stop("No weights found for raster: ", nm)
    
    out <- file.path(terra_dir, paste0(out_prefix, "_", nm, ".tif"))
    
    if (!file.exists(out) || FORCE_REBUILD) {
      
      w <- ecodist_r
      
      missing_ids <- setdiff(all_ids, lookup[, 1])
      
      if (length(missing_ids) > 0) {
        log_msg(
          nm, ": ",
          length(missing_ids),
          " ecodistrict(s) have no weight and will be set to NA:",
          paste(sort(missing_ids), collapse = ", ")
        )
        
        w <- subst(w,from = missing_ids, to = rep(NA, length(missing_ids)))
      }
      
      w <- subst(w, from = lookup[, 1], to = lookup[, 2])
      
      # QA
      mx <- global(w, "max", na.rm = TRUE)[1,1]
      
      if (mx > max(lookup[,2], na.rm = TRUE))
        stop("Unexpected values remain in weight raster for ", nm)
      
      writeRaster(w, out,
                  overwrite = TRUE,
                  wopt = list(datatype = "FLT4S"))
    }
    
    out
  })
  
  s <- rast(unlist(weight_files))
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

# Confirm weights are within expected bounds
tmp_chk1 <- global(acc_bin_stack, c("min","max"), na.rm=TRUE)
log_msg(capture.output(tmp_chk1))

tmp_chk2 <- global(metric_tame_stack, c("min","max"), na.rm = TRUE)
log_msg(capture.output(tmp_chk2))

tmp_chk3 <- global(metric_native_stack, c("min","max"), na.rm = TRUE)
log_msg(capture.output(tmp_chk3))

tmp_chk4 <- global(metric_not_stack, c("min","max"), na.rm = TRUE)
log_msg(capture.output(tmp_chk4))

# --------------------------------------------
# GLOBAL NATIVE/TAME FALLBACK RATIOS
# --------------------------------------------
# When evidence is missing, allocate probability according to observed native/tame prevalence
log_msg("Computing global native/tame fallback ratios")

global_native_freq <- global(
  app(nt_stack == 2, mean, na.rm = TRUE),
  "mean",
  na.rm = TRUE
)[1,1]

global_tame_freq <- global(
  app(nt_stack == 1, mean, na.rm = TRUE),
  "mean",
  na.rm = TRUE
)[1,1]

fallback_native_global <-
  global_native_freq /
  (global_native_freq + global_tame_freq)

fallback_tame_global <-
  global_tame_freq /
  (global_native_freq + global_tame_freq)

# global_native_freq <- global(
#   app(nt_stack == 2, mean, na.rm = TRUE),
#   "mean",
#   na.rm = TRUE
# )[1,1]
# 
# global_tame_freq <- global(
#   app(nt_stack == 1, mean, na.rm = TRUE),
#   "mean",
#   na.rm = TRUE
# )[1,1]
# 
# fallback_native_global <-
#   global_native_freq /
#   (global_native_freq + global_tame_freq)
# 
# fallback_tame_global <-
#   global_tame_freq /
#   (global_native_freq + global_tame_freq)

if (!is.finite(fallback_native_global) ||
    !is.finite(fallback_tame_global)) {
  stop("Global fallback ratios could not be computed")
}

log_msg(
  sprintf(
    "Global fallback ratios: native=%.4f tame=%.4f",
    fallback_native_global,
    fallback_tame_global
  )
)

# should always be exactly 1:
log_msg(
  sprintf(
    "Fallback sum = %.6f",
    fallback_native_global +
      fallback_tame_global
  )
)

if (
  abs(
    fallback_native_global +
    fallback_tame_global - 1
  ) > 1e-6
) {
  stop("Fallback ratios do not sum to 1")
}

# =========================
# TILE GRID
# =========================
# ~2000–5000 tiles is ideal, 500-1500 for single node
#tiles <- makeTiles(pilot, n = 10) # for sa, make n = 1000

#tiles <- getTileExtents(pilot, c(500, 500))
tiles <- getTileExtents(sa, c(1000, 1000))

log_msg("Tile count: ", nrow(tiles))

# Save tiles for reproducibility
writeVector(vect(tiles), "tile_index.gpkg", overwrite=TRUE)

# =========================
# PROCESS FUNCTION
# =========================
tile_dir <- file.path(tmp_root, "tiles")
dir.create(tile_dir, recursive = TRUE, showWarnings = FALSE)

# Define maximum entropy
max_entropy <- log(3)

process_tile <- function(i) {
  
  tile <- ext(tiles[i, ])
  
  out_file <- file.path(tile_dir, sprintf("tile_%05d.tif", i))
  if (file.exists(out_file)) return(out_file)
  
  # crop all rasters
  b   <- crop(bin_stack, tile)
  
  # skip empty tiles
  has_cells <- global(
    !is.na(b[[1]]),
    "sum",
    na.rm = TRUE
  )[1,1]
  
  if (has_cells == 0) {
    log_msg(sprintf("Tile %d skipped (no data)", i))
    return(NULL)
  }
  
  nt  <- crop(nt_stack, tile)
  
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
  # create global denominators (counts nt 3 times)
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
  
  # take mean or will have one layer per dataset
  ratio_not_dataset <- app(
    ifel(nt_dataset_weight > 0,
         rnot / nt_dataset_weight,
         NA),
    mean,
    na.rm = TRUE
  )
  
  ratio_native_dataset <- app(
    ifel(nt_dataset_weight > 0,
         rn / nt_dataset_weight,
         NA),
    mean,
    na.rm = TRUE
  )
  
  ratio_tame_dataset <- app(
    ifel(nt_dataset_weight > 0,
         rt / nt_dataset_weight,
         NA),
    mean,
    na.rm = TRUE
  )
  
  # ratio_not_dataset <- ifel(
  #   nt_dataset_weight > 0,
  #   rnot / nt_dataset_weight,
  #   NA
  # )
  # 
  # ratio_native_dataset <- ifel(
  #   nt_dataset_weight > 0,
  #   rn / nt_dataset_weight,
  #   NA
  # )
  # 
  # ratio_tame_dataset <- ifel(
  #   nt_dataset_weight > 0,
  #   rt / nt_dataset_weight,
  #   NA
  # )
  
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
  K_ESS <- 3
  K_WEIGHT <- 1
  support_score <- (ess / (ess + K_ESS)) * (w_sum / (w_sum + K_WEIGHT))
  
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
  
  # neutral fallback
  # neutral <- p_grass * 0 + 0.5
  # share_native <- cover(share_native, neutral)
  # share_tame <- cover(share_tame, neutral)
  
  # use global fallbacks calculated before tiling 
  fallback_native <- fallback_native_global
  fallback_tame   <- fallback_tame_global
  
  # use fallback that reflects observed landscape instead of neutral value
  share_native <- ifel(is.na(share_native), fallback_native, share_native)
  share_tame   <- ifel(is.na(share_tame), fallback_tame, share_tame)
  
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
results <- lapply(seq_len(nrow(tiles)), process_tile)

print("All tiles completed.")

# to use tiles in r_tmp (must put old tiles in separate folder)
tile_dir <- file.path(tmp_root, "tiles")
files <- list.files(tile_dir, full.names = TRUE, pattern = "\\.tif$")

# Get tile outputs
#files <- unlist(results)
files <- unlist(Filter(Negate(is.null), results))

if (length(files) == 0) {
  files <- list.files(tile_dir, full.names=TRUE, pattern="\\.tif$")
}

# Ensure files exist
files <- files[file.exists(files)]

if (length(files) == 0) {
  stop("No tile files found for mosaicking")
}

# plot the tiles
tile_vect <- vect("tile_index.gpkg")
tile_vect <- vect(tiles)

tile_vect$tile_id <- seq_len(nrow(tile_vect))

plot(tile_vect,
     border = "grey")

text(
  centroids(tile_vect),
  labels = tile_vect$tile_id,
  cex = 0.4
)

# QA
cat("Tile files found:", length(files), "\n")
print(basename(files))
sapply(files, function(f) nlyr(rast(f)))

cat("Mosaicking", length(files), "tiles...\n")

# Check layer consistency - make sure correct number of layers
#stopifnot(all(sapply(files, function(f) nlyr(rast(f)) == 4)))
stopifnot(all(vapply(files, function(f) nlyr(rast(f)) == 29, logical(1))))

# Mosaic all tiles together
mosaic_out <- do.call(merge, lapply(files, rast))

#mosaic_out <- mosaic(rast(files), fun="first")
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

# write locally now
checkpoint_dir <- file.path(tmp_root, "mosaic_checkpoint")
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)

log_msg("Writing mosaic_out checkpoint")

writeRaster(
  mosaic_out,
  file.path(checkpoint_dir, "mosaic_out.tif"),
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

mosaic_out <- rast(file.path(tmp_root, "mosaic_checkpoint", "mosaic_out.tif"))

checkpoint_dir <- file.path(tmp_root, "mosaic_layers")
dir.create(checkpoint_dir,
           recursive = TRUE,
           showWarnings = FALSE)

for (nm in names(mosaic_out)) {
  
  log_msg("Checkpointing ", nm)
  
  writeRaster(
    mosaic_out[[nm]],
    file.path(
      checkpoint_dir,
      paste0(nm, ".tif")
    ),
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
}

log_msg("Finished writing mosaic_out checkpoint")

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

plot(num_bin_sources)
plot(num_nt_sources)

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
tmp_fraction <- fraction_binary_support +
  fraction_triclass_support

tmp_stats <- global(
  tmp_fraction,
  c("min","max","mean"),
  na.rm = TRUE
)

print(tmp_stats) # should be ~1 everywhere

#---------- calibrate support_score ----------------#
# Check distribution
#hist(values(support_score))
vals_support <- spatSample(
  support_score,
  size = 1000,
  na.rm = TRUE,
  values = TRUE
)

quantile(vals_support$support_score)

# 0%       25%       50%       75%      100% 
# 0.3201789 0.5831215 0.6392761 0.6885397 0.7182373 

png(
  file.path(diag_dir, "support_score_hist.png"),
  width = 2400,
  height = 1600,
  res = 300
)

hist(vals_support$support_score,
     main = "Support Score Distribution",
     xlab = "Support Score"
)

dev.off()
#hist(spatSample(support_score, size = 1e6, na.rm = TRUE, values = TRUE))
#quantile(values(spatSample(support_score, 1e6, na.rm = TRUE)))

aws.s3::put_object(
  file = file.path(diag_dir, "support_score_hist.png"),
  object = paste0(out_dir, "/diagnostics/support_score_hist.png"))

calc_support <- function(ess, w_sum, k_ess, k_w) {
  
  (ess / (ess + k_ess)) *
    (w_sum / (w_sum + k_w))
  
}

support_k1 <- calc_support(ess, sum_weight_present, k_ess = 1, k_w = 1)
support_k3 <- calc_support(ess, sum_weight_present, k_ess = 3, k_w = 1)
support_k5 <- calc_support(ess, sum_weight_present, k_ess = 5, k_w = 1)

# Compare distributions
# Increasing the ESS regularization constant reduces support scores across the landscape, 
# compressing the distribution toward lower values while preserving relative spatial variation. 
# This reflects stronger penalization of low effective sample sizes and results in more 
# conservative confidence estimates.
stats <- global(
  c(support_k1, support_k3, support_k5),
  c("min","max","mean"),
  na.rm = TRUE)

log_msg(
  "Support score sensitivity analysis (K_ESS = 1, 3, 5): min, max, and mean values. ",
  "Distribution summaries:\n",
  capture.output(print(stats)))

# Examine quantiles
# Distributions are very compressed - suggests support metric is being driven by a relatively 
# narrow range of ESS values and/or w_sum values - support is not separating pixels very strongly
global(c(support_k1, support_k3, support_k5), 
       quantile, probs=c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm=TRUE)

#        X5.      X25.      X50.      X75.      X95.
# ess   0.6500709 0.7062533 0.7498243 0.7931196 0.8078078
# ess.1 0.5122457 0.5853020 0.6392762 0.6889475 0.7071444
# ess.2 0.4202574 0.5005994 0.5571364 0.6103342 0.6290974

# ESS distribution
# effective sample size is high almost everywhere (median pixel has ESS = 10.6)
# support distribution is narrow because ESS itself is narrow (5% = 6.1, 95% = 13.1)
# Most pixels have similar effective support, meaning most pixels are informed by many datasets,
# extreme low-support areas are rare, and the support metric is behaving stably
global(
  ess,
  quantile,
  probs = c(0.05, 0.25, 0.5, 0.75, 0.95),
  na.rm = TRUE
)

#      X5.     X25.     X50.     X75.     X95.
# ess 6.121585 8.715552 10.55651 12.59069 13.10684

# testing larger K values more aligned with actual ESS distribution
# K=10 pushes median below 0.5, K=15 makes even the best-supported pixels barely reach 0.4.
# Given that median ESS is 10.6, a pixel with a fairly typical level of support ends up being 
# classified as only "moderately supported" under K=10 and "weakly supported" under K=15.
support_k8  <- calc_support(ess, sum_weight_present, 8, 1)
support_k10 <- calc_support(ess, sum_weight_present, 10, 1)
support_k15 <- calc_support(ess, sum_weight_present, 15, 1)

global(
  c(support_k8, support_k10, support_k15),
  quantile,
  probs = c(.05, 0.25, 0.5, 0.75, 0.95),
  na.rm = TRUE)

#             X5.      X25.      X50.      X75.      X95.
# ess   0.4202574 0.5005994 0.5571364 0.6103342 0.6290974
# ess.1 0.2896904 0.3667801 0.4216827 0.4752489 0.4926966
# ess.2 0.2211476 0.2905429 0.3392119 0.3891240 0.4049051

 #             X5.      X25.      X50.      X75.      X95.
 # ess   0.4202574 0.5005994 0.5571364 0.6103342 0.6290974
 # ess.1 0.2896904 0.3667801 0.4216827 0.4752489 0.4926966
 # ess.2 0.2211476 0.2905429 0.3392119 0.3891240 0.4049051

# Evaluate behavior versus source counts
# Support scores should increase when evidence increases
n_sources <- num_bin_sources + num_nt_sources
names(n_sources) <- "n_sources"

# mean support by source count
# support function should show more sources = higher support with diminishing returns
df <- spatSample(c(n_sources, support_score), size = 1000, na.rm = TRUE, as.df = TRUE)

# crashed memory?
# df %>%
#   mutate(n_sources = num_bin_sources + num_nt_sources) %>%
#   group_by(n_sources) %>%
#   summarise(mean_support = mean(support_score))

# plot support vs source count
png(
  file.path(diag_dir, "support_score_vs_source_count.png"),
  width = 2400,
  height = 1600,
  res = 300
)

plot(df$n_sources, df$support_score, pch = 16, cex = 0.4)
lines(aggregate(support_score ~ n_sources, data = df, mean), col = "red", lwd = 3)

dev.off()

aws.s3::put_object(
  file = file.path(diag_dir, "support_score_vs_source_count.png"),
  object = paste0(out_dir, "/diagnostics/support_score_vs_source_count.png"))

##### MAKE DECISION RE CALIBRATION ######
# Sensitivity analysis of the support metric showed that increasing the ESS regularization 
# constant progressively reduced support scores (median support: 0.75 for K=1, 0.64 for K=3, 
# 0.56 for K=5, 0.42 for K=10). Given a median effective sample size of 10.6 datasets, K=3 
# provided a balance between optimistic (K=1) and overly conservative (K≥10) support estimates 
# and was retained for subsequent analyses.

#--------------- Confidence --------------------#
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
#           min      max      mean  entropy 
# -9.999995e-07 1.098609 0.3041583

tmp_chk6 <- global(support_score, c("min","max","mean"), na.rm = TRUE)
log_msg(capture.output(tmp_chk6))
#       min       max      mean support_score 
# 0.2425772 0.7182373 0.6300239

tmp_chk7 <- global(p_not, c("min","max","mean"), na.rm = TRUE)
log_msg(capture.output(tmp_chk7))
# min max      mean      p_not   
# 0   1        0.7536173

tmp_chk8 <- global(p_tame, c("min","max","mean"), na.rm = TRUE)
log_msg(capture.output(tmp_chk8))
# min max       mean    p_tame   
# 0   1         0.09693468

tmp_chk9 <- global(p_native, c("min","max","mean"), na.rm = TRUE)
log_msg(capture.output(tmp_chk9))
# min max     mean     p_native   
# 0   1       0.149448

# check ESS more
global(ess, c("min","max","mean"), na.rm = TRUE)
#          min      max     mean
# ess 2.607269 13.48475 10.31373

vals_ess <- spatSample(
  ess,
  size = 1000,
  na.rm = TRUE,
  values = TRUE
)

quantile(vals_ess[[1]])
# 0%       25%       50%       75%      100% 
# 4.754179  8.653579 10.556506 12.456938 13.484747 

# check weight support
global(sum_weight_present,
       c("min","max","mean"),
       na.rm = TRUE)
#                        min      max   mean
# sum_weight_present 1.09071 7.236026 4.8036

quantile(
  spatSample(
    sum_weight_present,
    size = 1000,
    na.rm = TRUE,
    values = TRUE
  )[[1]]
)
# 0%      25%      50%      75%     100% 
# 2.115514 3.795466 4.595938 6.065835 7.236026 

global(num_bin_sources, c("min","max")) # NaN for both
global(num_nt_sources, c("min","max"))  # NaN for both

global(num_bin_sources, c("min","max","mean"), na.rm=TRUE)
#                 min max     mean
# num_bin_sources   0   8 6.327684

global(num_nt_sources, c("min","max","mean"), na.rm=TRUE)
#                min max     mean
# num_nt_sources   0   9 6.109177

freq(num_bin_sources)
#   layer value     count
# 1     1     0  72526859
# 2     1     3        12
# 3     1     4    460536
# 4     1     5     53628
# 5     1     6     43818
# 6     1     7 684249299
# 7     1     8     31688

freq(num_nt_sources)
#    layer value     count
# 1      1     0  72529091
# 2      1     1    274574
# 3      1     2  21720348
# 4      1     3  43970777
# 5      1     4  42343664
# 6      1     5  43414436
# 7      1     6  17743230
# 8      1     7 164619630
# 9      1     8 350740369
# 10     1     9      9721

# Convert confidence to uncertainty
confidence_uncertainty <- 1 - confidence

names(confidence_uncertainty) <- "confidence_uncertainty"

# DECISION SCORE
# A pixel should only receive a high final score if both the class prediction is 
# confident AND there is strong support behind that prediction
# NOTE: decision score can never exceed the smaller term (scale can be compressed)

# adjusted confidence
#confidence <- app(prob_stack, max)

decision_score <- confidence * support_score
names(decision_score) <- "decision_score"

# Build a combined uncertainty layer

# Equal weighting: classifcation confidence, class ambiguity, amount of supporting eveidence
uncertainty_total <- (confidence_uncertainty + entropy_norm + support_uncertainty) / 3

# data-driven weighting:
# standardize inputs
z_conf <- scale(confidence_uncertainty)
z_ent  <- scale(entropy_norm)
z_sup  <- scale(support_uncertainty)

# examine correlations
# if two layers are highly redundant (r > 0.8), don't give equal weight
tmp <- spatSample(
  c(
    confidence_uncertainty,
    entropy_norm,
    support_uncertainty
  ),
  size = 10000,
  as.df = TRUE,
  na.rm = TRUE
)

cor(tmp)

# with size = 1,000:
#                        confidence_uncertainty entropy_norm support_uncertainty
# confidence_uncertainty              1.0000000    0.9694997          -0.1564002
# entropy_norm                        0.9694997    1.0000000          -0.1804322
# support_uncertainty                -0.1564002   -0.1804322           1.0000000

# with size = 10,000:
#                        confidence_uncertainty entropy_norm support_uncertainty
# confidence_uncertainty              1.0000000    0.9657242          -0.1355265
# entropy_norm                        0.9657242    1.0000000          -0.1787427
# support_uncertainty                -0.1355265   -0.1787427           1.0000000

# confidence_uncertainty and entropy_norm are highly correlated (0.97) - only use one
uncertainty_total <- (entropy_norm + support_uncertainty) / 2
names(uncertainty_total) <- "uncertainty_total"

# Weight entropy and confidence more than support:
uncertainty_total_w <- (
  0.6 * entropy_norm +
  0.4 * support_uncertainty)

names(uncertainty_total_w) <- "uncertainty_total_weighted"

#------ Hard classification ------#
hard_class <- which.max(prob_stack)

names(hard_class) <- "hard_class"

# Optional labeling (may not persist in GeoTIFF)
levels(hard_class) <- data.frame(
  value = 1:3,
  class = c("non_grass", "tame", "native")
)

writeRaster(hard_class, file.path(checkpoint_dir, "hard_class.tif"),
  overwrite = TRUE, wopt = list(datatype = "INT1U",
    gdal = c("COMPRESS=ZSTD", "TILED=YES", "BIGTIFF=YES")))

#------ High-confidence-only classification ------#
# Check to calibrate confidence threshold
vals_conf <- spatSample(
  decision_score,
  size = 1000000,
  na.rm = TRUE,
  values = TRUE)

png(
  file.path(diag_dir, "decision_score_hist4.png"),
  width = 2400,
  height = 1600,
  res = 300
)

hist(
  vals_conf$decision_score,
  main = "Decision Score Distribution",
  xlab = "Decision Score"
)

dev.off()

aws.s3::put_object(
  file = file.path(diag_dir, "decision_score_hist4.png"),
  object = paste0(out_dir, "/diagnostics/decision_score_hist4.png"))

summary(vals_conf$decision_score)
# with size = 1,000:
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.2231  0.4714  0.5798  0.5478  0.6394  0.7182 

# with size = 10,000:
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.1903  0.4596  0.5746  0.5428  0.6384  0.7182 

# with size = 100,000:
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.1702  0.4639  0.5759  0.5446  0.6384  0.7182 

# with size = 1,000,000:
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.1644  0.4630  0.5759  0.5445  0.6384  0.7182 

decision_q <- quantile(
  vals_conf$decision_score,
  probs = c( 0, 0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99, 1))

log_msg("Quantiles for decision_score:")
print(decision_q)
# with size = 1,000:
# 0%        1%        5%        10%       25%       50%       75%       90%       95%       99%       100% 
# 0.2230858 0.2574506 0.3087435 0.3620848 0.4713789 0.5797557 0.6393719 0.6857520 0.6977146 0.7071473 0.7182373 

# with size = 1,000,000:
# 0%        1%        5%        10%       25%       50%       75%       90%       95%       99%       100% 
# 0.1644407 0.2534290 0.3073677 0.3502129 0.4630198 0.5759259 0.6383783 0.6868626 0.6979046 0.7074378 0.7182373 

# 100k and 1M samples nearly identical - decision score distribution is well characterized

# coverage should be going down not up:
# 10% = 0.35 = 35% coverage = inclusive map, keep nearly everything (threshold at 10% retains 90% of valid pixels)
# 25% = 0.46 = 46% coverage = moderate filtering
# 50% = 0.58 = 58% coverage = retain only above-median certainty
# 75% = 0.64 = 64% coverage = high-confidence map
# 90% = 0.69 = 69% coverage = very high-confidence map (threshold at 90% retains 10% of valid pixels)

# check
mean(vals_conf$decision_score >= decision_q["10%"]) # = 0.9
mean(vals_conf$decision_score >= decision_q["25%"]) # = 0.750083
mean(vals_conf$decision_score >= decision_q["50%"]) # = 0.503716

# generate several candidate maps
# allows choice of coverage rather than arbitrary score -
# use instead of "decision_score < 0.4" in hard_class
thr_candidates <- c(
  decision_q[["10%"]], # 90% coverage
  decision_q[["25%"]], # 75% coverage
  decision_q[["50%"]]  # 50% coverage
)

coverage_tbl <- data.frame()

valid_n <- global(!is.na(decision_score), "sum", na.rm = TRUE)[1,1]

for (thr in thr_candidates) {
  retained_n <- global(ifel(decision_score >= thr, 1, 0), "sum", na.rm = TRUE)[1,1]
  coverage_tbl <- rbind(coverage_tbl, 
                        data.frame(threshold = thr, coverage = retained_n / valid_n))
}

coverage_tbl
#   threshold  coverage
# 1 0.3502977 0.9004411
# 2 0.4632833 0.7501316
# 3 0.5759259 0.5041226

# Threshold  Coverage  Interpretation
# 0.35       90%       Very inclusive
# 0.46       75%       Moderate filtering
# 0.58       50%       High-confidence half of the landscape
# 0.64       25%       Very high confidence
# 0.69       10%       Exceptional confidence only

log_msg(
  "Decision-score threshold calibration complete. ",
  "Thresholds corresponding to approximately 90%, 75%, and 50% coverage are ",
  sprintf("%.3f, %.3f, %.3f",
          decision_q[["10%"]],
          decision_q[["25%"]],
          decision_q[["50%"]])
)


# do same for entropy
vals_entropy <- spatSample(entropy, size = 10000, na.rm = TRUE, values = TRUE)[[1]] # hung with size = 100,000

entropy_q <- quantile(vals_entropy, probs = c(0.50, 0.75, 0.90, 0.95))

entropy_q
# with size = 1,000:
#   50%       75%       90%       95% 
#   0.1033858 0.5905285 0.9000026 1.0083001 

# with size = 10,000:
#   50%       75%       90%       95% 
#   0.1304808 0.5766452 0.8724542 1.0041316

# instead of spatSample:
global(entropy, quantile, probs = c(0.5, 0.75, 0.9, 0.95), na.rm = TRUE)
#              X50.      X75.     X90.     X95.
# entropy 0.1347117 0.5817126 0.892818 1.014815

entropy_q <- global(
  entropy,
  quantile,
  probs = c(0.5, 0.75, 0.9, 0.95),
  na.rm = TRUE
)

entropy_q

# generate several candidate maps
# allows choice of coverage rather than arbitrary threshold -
# use instead of "entropy > 0.8" in hard_class
# goal is retain low entropy, remove high entropy
thr_candidates_e <- as.numeric(entropy_q[1, ])

coverage_tbl_e <- data.frame()

valid_n_e <- global(!is.na(entropy), "sum", na.rm = TRUE)[1,1]

for (thr in thr_candidates_e) {
  retained_n_e <- global(ifel(entropy <= thr, 1, 0), "sum", na.rm = TRUE)[1,1]
  coverage_tbl_e <- rbind(coverage_tbl_e, 
                          data.frame(threshold = thr, coverage = retained_n_e / valid_n_e))
}

coverage_tbl_e

#   threshold  coverage
# 1 0.1347117 0.5000050
# 2 0.5817126 0.7500000
# 3 0.8928180 0.9000001
# 4 1.0148147 0.9500000

# Entropy quantile  Threshold  Coverage retained  Pixels removed
# 50th              0.135      50%                highest 50% entropy
# 75th              0.582      75%                highest 25% entropy
# 90th              0.893      90%                highest 10% entropy
# 95th              1.015      95%                highest 5% entropy

log_msg(
  "Entropy threshold calibration complete. ",
  "Retaining the lowest-entropy 50%, 75%, 90%, and 95% of valid pixels yields thresholds of ",
  sprintf("%.3f, %.3f, %.3f, %.3f",
          entropy_q[1, "X50."],
          entropy_q[1, "X75."],
          entropy_q[1, "X90."],
          entropy_q[1, "X95."])
)

# hard classification (1=not, 2=tame, 3=native)
# mask low-confidence or high-uncertainty pixels
# entropy: uncertainty masking - remove unreliable areas
# goal: keep coverage high, penalize low-support areas smoothly instead of abruptly, 
# and prevent a large number of tiny-weight rasters from producing deceptively good ESS values
# decision_threshold <- decision_q[["25%"]]    # ~0.463
# entropy_threshold  <- entropy_q[1, "X90."]   # ~0.893
# 
# high_conf_class <- ifel(
#   (decision_score < decision_threshold) |
#     (entropy > entropy_threshold),
#   NA,
#   high_conf_class
# )

#hard_class <- rast(file.path(tmp_root, "mosaic_layers", "hard_class.tif"))

# decision score already combines classification confidence + support/evidence and
# entropy is highly correlated with confidence_uncertainty so adds relatively little 
# independent information beyond what is already captured in the class probabilities.
# can just use decision score:

# high confidence map (50% coverage)
high_conf_class <- ifel(
  decision_score < decision_q[["50%"]],
  NA,
  hard_class)

# moderate map (retain 75% coverage)
mod_class <- ifel(
  decision_score < decision_q[["25%"]],
  NA,
  hard_class)

# more inclusive map (retain 90% coverage)
inclusive_class <- ifel(
  decision_score < decision_q[["10%"]],
  NA,
  hard_class)

plot(hard_class)
plot(inclusive_class)
plot(mod_class)
plot(high_conf_class)

checkpoint_dir <- file.path(tmp_root, "mosaic_layers")

writeRaster(high_conf_class, file.path(checkpoint_dir, "high_conf_class.tif"),
            overwrite = TRUE, wopt = list(datatype = "INT1U",
            gdal = c("COMPRESS=ZSTD", "TILED=YES", "BIGTIFF=YES")))

writeRaster(mod_class, file.path(checkpoint_dir, "mod_class.tif"),
            overwrite = TRUE, wopt = list(datatype = "INT1U",
                                          gdal = c("COMPRESS=ZSTD", "TILED=YES", "BIGTIFF=YES")))

writeRaster(inclusive_class, file.path(checkpoint_dir, "inclusive_class.tif"),
            overwrite = TRUE, wopt = list(datatype = "INT1U",
                                          gdal = c("COMPRESS=ZSTD", "TILED=YES", "BIGTIFF=YES")))

# mosaic_dir <- file.path(tmp_root, "mosaic_layers")
# 
# files_to_upload <- c(
#   "inclusive_class.tif",
#   "mod_class.tif",
#   "high_conf_class.tif"
# )
# 
# for (f in files_to_upload) {
#   
#   local_file <- file.path(mosaic_dir, f)
#   
#   log_msg("Uploading ", f)
#   
#   aws.s3::put_object(
#     file = local_file,
#     object = paste0(
#       out_dir,
#       "/mosaic_layers/",
#       f
#     ),
#     multipart = TRUE
#   )
# }
# 
# log_msg("Finished uploading filtered classification maps")

# save as PNG
png(
  file.path(diag_dir, "hard_class_map.png"),
  width = 2400,
  height = 1600,
  res = 300
)

plot(
  hard_class,
  main = "Final classification (unfiltered)",
  col = viridis::viridis(100)
)

dev.off()

aws.s3::put_object(
  file = file.path(diag_dir, "hard_class_map.png"),
  object = paste0(out_dir, "/diagnostics/hard_class_map.png"))

# save as PNG
png(
  file.path(diag_dir, "inclusive_filt_class_map.png"),
  width = 2400,
  height = 1600,
  res = 300
)

plot(
  inclusive_class,
  main = "Final classification (filtered, drop lowest 10% by decision score)",
  col = viridis::viridis(100)
)

dev.off()

aws.s3::put_object(
  file = file.path(diag_dir, "inclusive_filt_class_map.png"),
  object = paste0(out_dir, "/diagnostics/inclusive_filt_class_map.png"))

# save as PNG
png(
  file.path(diag_dir, "mod_filt_class_map.png"),
  width = 2400,
  height = 1600,
  res = 300
)

plot(
  mod_class,
  main = "Final classification (filtered, drop lowest 25% by decision score)",
  col = viridis::viridis(100)
)

dev.off()

aws.s3::put_object(
  file = file.path(diag_dir, "mod_filt_class_map.png"),
  object = paste0(out_dir, "/diagnostics/mod_filt_class_map.png"))

# save as PNG
png(
  file.path(diag_dir, "high_conf_filt_class_map.png"),
  width = 2400,
  height = 1600,
  res = 300
)

plot(
  high_conf_class,
  main = "Final classification (filtered, drop lowest 50% by decision score)",
  col = viridis::viridis(100)
)

dev.off()

aws.s3::put_object(
  file = file.path(diag_dir, "high_conf_filt_class_map.png"),
  object = paste0(out_dir, "/diagnostics/high_conf_filt_class_map.png"))

# is tame underrepresented?
global(p_not, "mean", na.rm = TRUE)    # mean = 0.7536173
global(p_tame, "mean", na.rm = TRUE)   # mean = 0.09693468
global(p_native, "mean", na.rm = TRUE) # mean = 0.149448
freq(hard_class)

#   layer     value     count
# 1     1 non_grass 533382526 = 78% (+3%)
# 2     1      tame  53304622 = 8%  (-2%)
# 3     1    native  98148288 = 14% (-1%)

# tame probability distribution
# if tame rarely exceeds 0.4 then it almost never wins argmax competition
tame_vals <- spatSample(p_tame, size = 10000, values = TRUE)
quantile(tame_vals[[1]], na.rm = TRUE)

# 0%        25%       50%       75%       100% 
# 0.0000000 0.0000000 0.0000000 0.1161091 0.9962528 

# the scarcity of tame in the hard map originates in the probability surfaces themselves

# distribution only where tame wins
# if tame pixels require much higher probabilities to win than native pixels, 
# then there may still be a structural bias
tame_mask <- hard_class == 2
global(mask(p_tame, tame_mask), c("mean", "min", "max"), na.rm = TRUE)
#              mean min max                 
# p_tame 0.09693468   0   1

# compare to native
native_mask <- hard_class == 3
global(mask(p_native, native_mask), c("mean", "min", "max"), na.rm = TRUE)
#              mean min max                 
# p_native 0.149448   0   1

# look at class-specific weights
summary(acc_nt$tame_metric)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0000  0.2376  0.5082  0.4907  0.7465  1.0000 
summary(acc_nt$native_metric)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0000  0.3471  0.6116  0.5522  0.7863  1.0000 
summary(acc_nt$not_metric)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0000  0.3596  0.7254  0.6280  0.9193  1.0000 

# are models disagreeing more on tame?
tame_tmp <- spatSample(
  c(p_tame,
    triclass_agreement,
    entropy,
    support_score),
  size = 1000,
  as.df = TRUE,
  na.rm = TRUE
)

tame_tmp

# correlations
cor(
  tame_tmp[, c(
    "p_tame",
    "triclass_agreement",
    "entropy",
    "support_score"
  )],
  use = "complete.obs"
)

#                        p_tame triclass_agreement    entropy support_score
# p_tame              1.0000000         -0.6967515  0.6875212     0.0285825
# triclass_agreement -0.6967515          1.0000000 -0.8535302    -0.1599136
# entropy             0.6875212         -0.8535302  1.0000000     0.1843801
# support_score       0.0285825         -0.1599136  0.1843801     1.0000000

# compare class proportions before and after masking
get_class_props <- function(r) {
  
  f <- freq(r)
  
  out <- data.frame(
    class = as.character(f$value),
    count = f$count
  )
  
  out$prop <- out$count / sum(out$count)
  
  out
}

full_props <- get_class_props(hard_class)
full_props$map <- "full"

high_props <- get_class_props(high_conf_class)
high_props$map <- "high"

mod_props <- get_class_props(mod_class)
mod_props$map <- "moderate"

incl_props <- get_class_props(inclusive_class)
incl_props$map <- "inclusive"

class_comp <- rbind(
  full_props,
  incl_props,
  mod_props,
  high_props
)

class_comp
#        class     count         prop       map
# 1  non_grass 533382526 0.7788477318      full
# 2       tame  53304622 0.0778356656      full
# 3     native  98148288 0.1433166026      full
# 4  non_grass 505345774 0.8194965010 inclusive
# 5       tame  28872362 0.0468210102 inclusive
# 6     native  82435838 0.1336824888 inclusive
# 7  non_grass 444730051 0.8657106918  moderate
# 8       tame   5893526 0.0114723268  moderate
# 9     native  63093136 0.1228169814  moderate
# 10 non_grass 302562933 0.8763816999      high
# 11      tame    107155 0.0003103773      high
# 12    native  42570956 0.1233079228      high

class_comp %>%
  select(map, class, prop) %>%
  pivot_wider(
    names_from = class,
    values_from = prop
  )
#   map       non_grass     tame native
# 1 full          0.779 0.0778    0.143
# 2 inclusive     0.819 0.0468    0.134
# 3 moderate      0.866 0.0115    0.123
# 4 high          0.876 0.000310  0.123

baseline <- full_props %>%
  select(class, full_prop = prop)

comparison <- class_comp %>%
  filter(map != "full") %>%
  left_join(baseline, by = "class") %>%
  mutate(
    pct_change = 100 * (prop - full_prop) / full_prop
  )

comparison
#       class     count         prop       map  full_prop pct_change
# 1 non_grass 505345774 0.8194965010 inclusive 0.77884773   5.219091
# 2      tame  28872362 0.0468210102 inclusive 0.07783567 -39.846329
# 3    native  82435838 0.1336824888 inclusive 0.14331660  -6.722259
# 4 non_grass 444730051 0.8657106918  moderate 0.77884773  11.152753
# 5      tame   5893526 0.0114723268  moderate 0.07783567 -85.260836
# 6    native  63093136 0.1228169814  moderate 0.14331660 -14.303731
# 7 non_grass 302562933 0.8763816999      high 0.77884773  12.522854
# 8      tame    107155 0.0003103773      high 0.07783567 -99.601240
# 9    native  42570956 0.1233079228      high 0.14331660 -13.961174

# absolute retention by class
retention_by_class <- function(full_map, filtered_map) {
  
  full_freq <- freq(full_map)
  filt_freq <- freq(filtered_map)
  
  merge(
    data.frame(class = full_freq$value,
               full_n = full_freq$count),
    data.frame(class = filt_freq$value,
               kept_n = filt_freq$count),
    by = "class",
    all.x = TRUE
  ) %>%
    mutate(
      kept_n = replace_na(kept_n, 0),
      retention = kept_n / full_n
    )
}

retention_by_class(hard_class, inclusive_class)
#       class    full_n    kept_n retention
# 1    native  98148288  82435838 0.8399111
# 2 non_grass 533382526 505345774 0.9474359
# 3      tame  53304622  28872362 0.5416484

retention_by_class(hard_class, mod_class)
#       class    full_n    kept_n retention
# 1    native  98148288  63093136 0.6428348
# 2 non_grass 533382526 444730051 0.8337919
# 3      tame  53304622   5893526 0.1105631

retention_by_class(hard_class, high_conf_class)
#       class    full_n    kept_n   retention
# 1    native  98148288  42570956 0.433741198
# 2 non_grass 533382526 302562933 0.567253178
# 3      tame  53304622    107155 0.002010238

# map pixels where tame is plausible but not dominant
tame_possible <- p_tame >= 0.3
tame_dominant <- hard_class == 2
missed_tame <- tame_possible & !tame_dominant

png(
  file.path(diag_dir, "missed_tame.png"),
  width = 2400,
  height = 1600,
  res = 300
)

plot(missed_tame, main = "Non-tame pixels where tame possible (p >= 0.3)")

dev.off()

aws.s3::put_object(
  file = file.path(diag_dir, "missed_tame.png"),
  object = paste0(out_dir, "/diagnostics/missed_tame.png"))


#------------- Confidence zones -------------#
# use quantile values to define:
# Zone 1 = Very low confidence (cutoff is 10% value)
# Zone 2 = Low confidence (cutoffs are 10 & 25% values)
# Zone 3 = Moderate confidence (cutoffs are 25 & 50% values)
# Zone 4 = High confidence (cutoff is 50% but max value reached at 90%)
# confidence_zone <- classify(
#   decision_score,
#   matrix(
#     c(
#       -Inf, 0.36, 1,  # very low
#       0.36, 0.50, 2,  # low
#       0.50, 0.62, 3,  # moderate
#       0.62, Inf,  4   # high
#     ),
#     ncol = 3,
#     byrow = TRUE
#   )
# )

q10 <- decision_q[["10%"]]
q25 <- decision_q[["25%"]]
q50 <- decision_q[["50%"]]

confidence_zone <- classify(
  decision_score,
  matrix(c(
    -Inf,q10,1,
    q10,q25,2,
    q25,q50,3,
    q50,Inf,4
  ), ncol = 3,
  byrow = TRUE)
)

names(confidence_zone) <- "confidence_zone"
plot(confidence_zone)

# combine class and decision score
# 11 = Not grass, very low confidence
# 14 = Not grass, high confidence
# 
# 21 = Tame, very low confidence
# 24 = Tame, high confidence
# 
# 31 = Native, very low confidence
# 34 = Native, high confidence
combined <- hard_class * 10 + confidence_zone
names(combined) <- "combined"

freq(combined)
plot(combined)

###################################################

#------ Tiled temporal workflow (grass vs nongrass) ------#

# if starting here
#bin_stack <- rast(file.path(stack_dir, "bin_stack.tif"))
#nt_stack <- rast(file.path(stack_dir, "nt_stack.tif"))

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

log_msg("Rasters per period:")
print(table(bin_meta_pd$period))

all(bin_meta$raster_id %in% names(bin_stack))
all(nt_meta$raster_id %in% names(nt_stack))

# Build annual stack
years <- sort(unique(bin_meta$year))

year_rasters <- lapply(years, function(y) {
  
  rasters_this_year <- bin_meta$raster_id[bin_meta$year == y]
  
  s <- bin_stack[[rasters_this_year]]
  
  if (nlyr(s) == 1) {
    s
  } else {
    app(s, mean, na.rm = TRUE)
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
temporal_tile_dir <- file.path(tmp_root, "temporal_tiles")
dir.create(temporal_tile_dir, recursive = TRUE, showWarnings = FALSE)

# Tile processor
process_temporal_tile <- function(i, temporal_stack) {
  
  tile <- ext(tiles[i, ])
  
  out_file <- file.path(
    temporal_tile_dir,
    sprintf("temporal_%05d.tif", i)
  )
  
  if (file.exists(out_file)) {
    return(out_file)
  }
  
  ys <- crop(temporal_stack, tile)
  
  has_cells <- global(!is.na(ys[[1]]), "sum", na.rm = TRUE)[1,1]
  
  if (has_cells == 0) {
    return(NULL)
  }
  
  # -----------------------
  # Temporal variance
  # -----------------------
  
  temporal_var <- app(
    ys,
    function(x) {
      
      x <- x[!is.na(x)]
      
      if (length(x) < 2) {
        return(NA)
      }
      
      var(x)
    }
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
      
      if (all(is.na(x))) {
        return(NA)
      }
      
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
temporal_files <- lapply(seq_len(nrow(tiles)), process_temporal_tile, temporal_stack = temporal_stack)
temporal_files <- unlist(Filter(Negate(is.null), temporal_files))
temporal_files <- temporal_files[file.exists(temporal_files)]

length(temporal_files)

rast(temporal_files[1])

tmp <- rast(temporal_files[1])

names(tmp)

temporal_range <- tmp[["temporal_range"]]

global(
  temporal_range,
  c("min", "max", "mean"),
  na.rm = TRUE
)

global(temporal_range,
       c("min","max","mean"),
       na.rm = TRUE)

freq(is.infinite(temporal_range))

sapply(
  temporal_files,
  function(f) nlyr(rast(f))
)


# Mosaic temporal outputs
cat("Mosaicking", length(temporal_files), "temporal tiles...\n")

temporal_bin_out <- do.call(merge,lapply(temporal_files, rast))
#temporal_out <- mosaic(rast(temporal_files), fun = "mean")

names(temporal_bin_out) <- c(
  "temporal_sd",
  "temporal_entropy",
  "temporal_range"#,
  #"temporal_slope"
)

# Extract layers
temporal_sd <- temporal_bin_out[["temporal_sd"]]
temporal_entropy <- temporal_bin_out[["temporal_entropy"]]
temporal_range <- temporal_bin_out[["temporal_range"]]

#temporal_slope <- temporal_out[["temporal_slope"]]

#------ Tiled temporal workflow (native vs tame vs nongrass) ------#

year_counts_nt <- table(nt_meta$year)

log_msg("Rasters per year:")
print(year_counts_nt)

# Aggregate if some years only have 1-2 rasters
nt_meta_pd <- nt_meta %>%
  mutate(
    period = case_when(
      year <= 2017 ~ "2015_2017",
      year <= 2020 ~ "2018_2020",
      year <= 2023 ~ "2021_2023",
      TRUE         ~ "2024_2025"
    )
  )

log_msg("Rasters per period:")
print(table(nt_meta_pd$period))

# Build annual stack
years_nt <- sort(unique(nt_meta$year))

year_rasters_nt <- lapply(years_nt, function(y) {
  
  rasters_this_year <- nt_meta$raster_id[nt_meta$year == y]
  
  s <- nt_stack[[rasters_this_year]]
  
  if (nlyr(s) == 1) {
    s
  } else {
    app(
      s,
      mean,
      na.rm = TRUE
    )
  }
})

year_stack_nt <- rast(year_rasters_nt)

names(year_stack_nt) <- paste0("y", years_nt)

# OR aggregated periods stack
periods_nt <- sort(unique(nt_meta_pd$period))

period_rasters_nt <- lapply(periods_nt, function(p) {
  
  idx <- which(nt_meta_pd$period == p)
  
  app(
    rast(nt_local[idx]),
    mean,
    na.rm = TRUE
  )
})

year_stack_pds_nt <- rast(period_rasters_nt)

names(year_stack_pds_nt) <- periods_nt

# CHOOSE WHETHER TO AGGREGATE YEARS OR NOT
USE_PERIODS_nt <- TRUE

temporal_stack_nt <- if (USE_PERIODS_nt) {
  log_msg("Using period-aggregated temporal stack (triclass)")
  year_stack_pds_nt
} else {
  log_msg("Using annual temporal stack (triclass)")
  year_stack_nt
}

# Create yearly class frequencies
# nt_year_probs <- lapply(nt_years, function(y) {
#   
#   idx <- which(nt_meta$year == y)
#   
#   s <- rast(nt_local[idx])
#   
#   p_not <- app(s == 0, mean, na.rm = TRUE)
#   p_tame <- app(s == 1, mean, na.rm = TRUE)
#   p_native <- app(s == 2, mean, na.rm = TRUE)
#   
#   names(p_not) <- paste0("not_", y)
#   names(p_tame) <- paste0("tame_", y)
#   names(p_native) <- paste0("native_", y)
#   
#   list(
#     p_not = p_not,
#     p_tame = p_tame,
#     p_native = p_native
#   )
# })

# Build time series stacks for each class
# not_stack <- rast(
#   lapply(nt_year_probs, function(x) x$p_not)
# )
# 
# tame_stack <- rast(
#   lapply(nt_year_probs, function(x) x$p_tame)
# )
# 
# native_stack <- rast(
#   lapply(nt_year_probs, function(x) x$p_native)
# )
# 
# names(not_stack) <- paste0("y", nt_years)
# names(tame_stack) <- paste0("y", nt_years)
# names(native_stack) <- paste0("y", nt_years)

# Create period class frequencies
nt_period_probs <- lapply(periods_nt, function(p) {
  
  idx <- which(nt_meta_pd$period == p)
  
  s <- rast(nt_local[idx])
  
  p_not <- app(s == 0, mean, na.rm = TRUE)
  p_tame <- app(s == 1, mean, na.rm = TRUE)
  p_native <- app(s == 2, mean, na.rm = TRUE)
  
  names(p_not) <- paste0("not_", p)
  names(p_tame) <- paste0("tame_", p)
  names(p_native) <- paste0("native_", p)
  
  list(
    p_not = p_not,
    p_tame = p_tame,
    p_native = p_native
  )
})

# Build period stacks
not_stack <- rast(
  lapply(nt_period_probs, function(x) x$p_not)
)

tame_stack <- rast(
  lapply(nt_period_probs, function(x) x$p_tame)
)

native_stack <- rast(
  lapply(nt_period_probs, function(x) x$p_native)
)

names(not_stack) <- periods_nt
names(tame_stack) <- periods_nt
names(native_stack) <- periods_nt

log_msg("Periods used in triclass temporal analysis:")
log_msg(paste(periods_nt, collapse = ", "))

log_msg("Layers in not_stack:")
log_msg(nlyr(not_stack))

log_msg("Layers in tame_stack:")
log_msg(nlyr(tame_stack))

log_msg("Layers in native_stack:")
log_msg(nlyr(native_stack))

# Tiled triclass temporal processor
process_triclass_temporal_tile <- function(i) {
  
  tile <- ext(tiles[i, ])
  
  out_file <- file.path(
    temporal_tile_dir,
    sprintf("triclass_temporal_%05d.tif", i)
  )
  
  if (file.exists(out_file)) {
    return(out_file)
  }
  
  not_y <- crop(not_stack, tile)
  tame_y <- crop(tame_stack, tile)
  native_y <- crop(native_stack, tile)
  
  has_cells <- global(
    !is.na(not_y[[1]]),
    "sum",
    na.rm = TRUE
  )[1,1]
  
  if (has_cells == 0) {
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
  
  dominant_class <- app(
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
triclass_temporal_files <- lapply(seq_len(nrow(tiles)), process_triclass_temporal_tile)
triclass_temporal_files <- unlist(Filter(Negate(is.null), triclass_temporal_files))
triclass_temporal_files <- triclass_temporal_files[file.exists(triclass_temporal_files)]

# Mosaic temporal outputs
cat("Mosaicking", length(triclass_temporal_files), "temporal tiles...\n")

triclass_temporal_out <- do.call(merge,lapply(triclass_temporal_files, rast))

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

plot(class_entropy)

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

tmp_chk10 <- global(
  c(
    support_score,
    confidence,
    decision_score
  ),
  c("min", "max", "mean"),
  na.rm = TRUE
)

for (i in seq_len(nrow(tmp_chk10))) {
  
  log_msg(
    sprintf(
      "%s: min=%.4f max=%.4f mean=%.4f",
      rownames(tmp_chk10)[i],
      tmp_chk10[i, "min"],
      tmp_chk10[i, "max"],
      tmp_chk10[i, "mean"]
    )
  )
}

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

# Visualize confidence vs entropy
# to understand whether masking rules are too strict/loose

# Sample ~1e6 pixels max (adjust if needed)
set.seed(42)

vals <- spatSample(
  c(decision_score, entropy),
  size = 1000,
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
decision_threshold <- decision_q[["10%"]]   

p1 <- ggplot(vals, aes(x = decision_score, y = entropy)) +
  stat_binhex(bins = 60) +
  scale_fill_viridis_c(trans = "log10") +
  geom_vline(xintercept = decision_threshold, color = "blue", linetype = "dashed") +
  #geom_hline(yintercept = entropy_threshold, color = "red", linetype = "dashed") +
  labs(
    title = "Decision score vs Entropy",
    x = "Decision score",
    y = "Entropy",
    fill = "Pixel count"
  ) +
  theme_minimal()

p1

ggsave(
  file.path(diag_dir, "decisionscore_vs_entropy.png"),
  p1,
  width = 8,
  height = 6,
  dpi = 300
)

# Density plot (smooth view)
p2 <- ggplot(vals, aes(x = decision_score, y = entropy)) +
  geom_density_2d_filled() +
  geom_vline(xintercept = decision_threshold, color = "blue", linetype = "dashed") +
  #geom_hline(yintercept = entropy_threshold, color = "red", linetype = "dashed") +
  labs(
    title = "Density of decision score vs Entropy"
  ) +
  theme_minimal()

p2

ggsave(
  file.path(diag_dir, "decisionscore_vs_entropy_density.png"),
  p2,
  width = 8,
  height = 6,
  dpi = 300
)

# Confidence distribution
# If most points cluster around confidence > 0.5 then threshold is appropriate
p3 <- ggplot(vals, aes(x = decision_score)) +
  geom_histogram(bins = 50, fill = "steelblue") +
  geom_vline(xintercept = decision_threshold, color = "red", linetype = "dashed") +
  labs(title = "Decision score distribution") +
  theme_minimal()

p3

ggsave(
  file.path(diag_dir, "decisionscore_dist.png"),
  p3,
  width = 8,
  height = 6,
  dpi = 300
)

# Entropy distribution
# If most points cluster around entropy < 0.6 then threshold is appropriate
p4 <- ggplot(vals, aes(x = entropy)) +
  geom_histogram(bins = 50, fill = "darkgreen") +
  geom_vline(xintercept = entropy_threshold, color = "red", linetype = "dashed") +
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

log_msg("Decision score:")
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

aws.s3::put_object(
  file = file.path(diag_dir, "entropy_map.png"),
  object = paste0(out_dir, "/diagnostics/entropy_map.png"))

# Map of confidence
png(
  file.path(diag_dir, "decision_score_map.png"),
  width = 2400,
  height = 1600,
  res = 300
)

plot(
  decision_score,
  main = "Decision score",
  col = viridis::viridis(100)
)

dev.off()

aws.s3::put_object(
  file = file.path(diag_dir, "decision_score_map.png"),
  object = paste0(out_dir, "/diagnostics/decision_score_map.png"))

# Binary "uncertain areas" map
uncertain_mask <- ifel(
  (decision_score < decision_threshold) | (entropy > entropy_threshold),
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

#------- Write outputs to S3 (if they are already in r_tmp) ----------#

# PNG files in r_tmp/diagnostics
diag_files <- list.files(
  diag_dir,
  recursive = TRUE,
  full.names = TRUE
)

upload_results <- sapply(diag_files, function(f) {
  
  rel_path <- gsub(
    paste0("^", normalizePath(diag_dir), "/?"),
    "",
    normalizePath(f)
  )
  
  aws.s3::put_object(
    file = f,
    object = file.path(out_dir, "diagnostics", rel_path)
  )
  
}, simplify = TRUE)

log_msg(
  sum(upload_results),
  " of ",
  length(upload_results),
  " diagnostic files uploaded successfully"
)

# TIFFs in r_tmp/mosaic_layers
tif_files <- list.files(
  file.path(tmp_root, "mosaic_layers"),
  pattern = "\\.tif$",
  full.names = TRUE
)

for (f in tif_files) {
  
  aws.s3::put_object(
    file = f,
    object = paste0(
      out_dir,
      "/mosaic_layers/",
      basename(f)
    )
  )
}

log_msg(
  "Uploaded ",
  length(tif_files),
  " TIFF files from mosaic_layers"
)

# log files
log_files <- list.files(
  tmp_root,
  pattern = "^consensus_mapping_2026082.*\\.log$",
  full.names = TRUE
)

log_files

for (f in log_files) {
  
  log_msg("Uploading ", basename(f))
  
  aws.s3::put_object(
    file = f,
    object = paste0(
      out_dir,
      "/logs/",
      basename(f)
    )
  )
}

log_msg("Uploaded ", length(log_files), " log file(s)")

#------- Write outputs to S3 (if they are not yet in r_tmp) ----------#
log_msg("Writing outputs")

aws.s3::put_object(
  file = file.path(diag_dir, "lambda_sensitivity.png"),
  object = paste0(out_dir, "/diagnostics/lambda_sensitivity.png"))

aws.s3::put_object(
  file = file.path(diag_dir, "lambda_sensitivity2.png"),
  object = paste0(out_dir, "/diagnostics/lambda_sensitivity2.png"))

aws.s3::put_object(
  file = file.path(diag_dir, "support_score_hist.png"),
  object = paste0(out_dir, "/diagnostics/support_score_hist.png"))

aws.s3::put_object(
  file = file.path(diag_dir, "decision_score_hist.png"),
  object = paste0(out_dir, "/diagnostics/decision_score_hist.png"))

aws.s3::put_object(
  file = file.path(diag_dir, "decisionscore_vs_entropy.png"),
  object = paste0(out_dir, "/diagnostics/decisionscore_vs_entropy.png"))

aws.s3::put_object(
  file = file.path(diag_dir, "decisionscore_vs_entropy_density.png"),
  object = paste0(out_dir, "/diagnostics/decisionscore_vs_entropy_density.png"))

aws.s3::put_object(
  file = file.path(diag_dir, "decisionscore_dist.png"),
  object = paste0(out_dir, "/diagnostics/decisionscore_dist.png"))

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
  list(raster = confidence_zone, local = "confidence_zone.tif", s3 = "confidence_zone.tif", wopt = wopt_class),
  list(raster = combined, local = "combined.tif", s3 = "combined.tif", wopt = wopt_class),
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
