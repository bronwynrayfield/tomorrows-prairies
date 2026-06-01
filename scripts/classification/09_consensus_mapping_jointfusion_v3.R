
###################################################################################
##                     GRASSLAND CONSENSUS MAPPING ANALYSIS                      ##
###################################################################################

##----------------------- Alison Long - May 2026 -------------------------

# step 8: use data fusion techniques to create a consensus grassland map 
#         to serve as baseline for modelling

# THIS VERSION IS OPTIMIZED FOR SAGEMAKER
# Next optimization to try: build stacks directly from S3 using GDAL VRT (zero downloads)

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

## ---- Setup ---- ##
# Temp dirs
#tmp_root <- file.path(Sys.getenv("HOME"), "r_tmp")
#dir.create(tmp_root, recursive = TRUE, showWarnings = FALSE)

scratch_dir <- "/mnt/tmp"

terra_dir <- file.path(scratch_dir, "terra")
data_dir  <- file.path(scratch_dir, "data")
tile_dir  <- file.path(scratch_dir, "tiles")
log_dir_local <- "/mnt/tmp/logs"

dir.create(terra_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(tile_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir_local, showWarnings = FALSE, recursive = TRUE)

terraOptions(
  tempdir = scratch_dir,
  memfrac = 0.8,
  progress = 1
)

# Load packages
library(terra)
library(future.apply)
library(readr)
library(dplyr)
library(aws.s3)

# Set terra options
#terraOptions(
#  tempdir = file.path(tmp_root, "terra"),
#  memfrac = 0.6, # could go down to 0.4
#  progress = 1
#)

## ---- Use limited parallel workers ---- ##
n_cores <- min(4, parallel::detectCores(logical = FALSE) - 1)

plan(multisession, workers = n_cores)

log_msg("Using parallel workers: ", n_cores)

## ---- Settings ---- ##
meta_file  <- "s3://sagemaker-gst-stage.sharing/alison-long/tomorrows-prairies/processed/metadata_full_extent_binary_nattam_with_path.csv"
out_dir    <- "s3://sagemaker-gst-stage.sharing/alison-long/tomorrows-prairies/results/consensus_map"
log_dir    <- "s3://sagemaker-gst-stage.sharing/alison-long/tomorrows-prairies/logs"

## ---- Message log ---- ##
log_file <- file.path(log_dir, paste0(
  "consensus_map_log_",
  format(Sys.time(), "%Y%m%d_%H%M%S"),
  ".log"
))

# Ensure local temp file first 
log_file_local <- file.path(
  log_dir_local,
  paste0("consensus_map_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log")
)

log_msg <- function(..., level = "INFO") {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  msg <- paste0("[", level, "] ", ts, " : ", paste(..., collapse = ""))
  
  # Console
  message(msg)
  
  # Append to log file
  cat(msg, "\n", file = log_file_local, append = TRUE)
}

log_msg("Starting data fusion analysis")
log_msg("Parallel workers detected: ", n_cores)

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

## ---- Split by stack type ---- ##
bin_meta <- meta %>% filter(stack_type == "binary")
nt_meta  <- meta %>% filter(stack_type == "nt")

## ---- Named file vectors ---- ##
#rfiles_s3 <- setNames(meta$file,     meta$raster_id)
bin_files <- setNames(bin_meta$file, bin_meta$raster_id)
nt_files  <- setNames(nt_meta$file,  nt_meta$raster_id)

log_msg("Binary rasters: ", nrow(bin_meta))
log_msg("Tri-class rasters: ", nrow(nt_meta))

## ---- Download rasters locally ---- ##
message("Downloading rasters locally...")
log_msg("Downloading rasters locally")

download_files <- function(files_named) {
  sapply(names(files_named), function(nm) {
    dest <- file.path(local_dir, paste0(nm, ".tif"))
    if (!file.exists(dest)) {
      aws.s3::save_object(files_named[[nm]], file = dest)
    }
    dest
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

# Grass vs nongrass
bin_stack_file <- "/mnt/tmp/bin_stack.tif"

if (!file.exists(bin_stack_file)) {
  writeRaster(
    rast(bin_local),
    bin_stack_file,
    overwrite = TRUE,
    wopt = list(
      datatype="INT1U",
      gdal=c("COMPRESS=ZSTD","TILED=YES","BIGTIFF=YES")
    )
  )
}

bin_stack <- rast(bin_stack_file)
names(bin_stack) <- names(bin_local)

# Native vs tame vs nongrass
nt_stack_file <- "/mnt/tmp/nt_stack.tif"

if (!file.exists(nt_stack_file)) {
  writeRaster(
    rast(nt_local),
    nt_stack_file,
    overwrite = TRUE,
    wopt = list(
      datatype="INT1U",
      gdal=c("COMPRESS=ZSTD","TILED=YES","BIGTIFF=YES")
    )
  )
}

nt_stack <- rast(nt_stack_file)
names(nt_stack) <- names(nt_local)

# =========================
# GLOBAL OPTIONS
# =========================
terraOptions(
  memfrac = 0.8,
  tempdir = "/mnt/tmp", # use fast scratch disk
  progress = 1
)

# increase GDAL block cache size to 8192 MB (8 GB)
Sys.setenv(GDAL_CACHEMAX = 8192)

# =========================
# INPUTS
# =========================
sa <- rast("classification_study_area_aligned.tif")

#bin_stack   <- rast("binary_stack_masked.tif")
#nt_stack    <- rast("tri_class_stack_masked.tif")

recall_bin_stack    <- rast("recall_bin.tif")
recall_tame_stack   <- rast("recall_tame.tif")
recall_native_stack <- rast("recall_native.tif")

penalty_bin_stack <- rast("penalty_bin.tif")
penalty_nt_stack  <- rast("penalty_nt.tif")

out_dir <- "tiles_out"
dir.create(out_dir, showWarnings = FALSE)

# =========================
# TILE GRID
# =========================
# ~2000–5000 tiles is ideal, 500-1500 for single node
tiles <- makeTiles(sa, n = 2000, overlap = 0)

# Save tiles for reproducibility
writeVector(vect(tiles), "tile_index.gpkg", overwrite=TRUE)

# =========================
# PROCESS FUNCTION
# =========================
process_tile <- function(i) {
  
  tile <- tiles[[i]]
  
  out_file <- file.path(out_dir, sprintf("tile_%05d.tif", i))
  if (file.exists(out_file)) return(out_file)
  
  # crop all rasters
  b   <- crop(bin_stack, tile)
  
  # skip empty tiles
  if (is.na(global(b[[1]], "max")[[1]])) {
    log_msg(sprintf("Tile %d skipped (no data)", i))
    return(NULL)
  }
  
  nt  <- crop(nt_stack, tile)
  
  rbin <- crop(recall_bin_stack, tile)
  rt   <- crop(recall_tame_stack, tile)
  rn   <- crop(recall_native_stack, tile)
  
  pbin <- crop(penalty_bin_stack, tile)
  pnt  <- crop(penalty_nt_stack, tile)
  
  # effective weight = accuracy x redundancy penalty
  bin_w_eff   <- rbin * pbin
  nt_w_tame   <- rt * pnt
  nt_w_native <- rn * pnt
  
  NEG_HALF <- 0.5
  
  # ----------------------
  # EVIDENCE COMPUTATION
  # ----------------------
  
  # ----------------------
  # NORMALIZED WEIGHTS
  # ----------------------
  
  # binary weights
  bin_w_eff <- rbin * pbin
  bin_w_sum <- sum(bin_w_eff, na.rm=TRUE)
  bin_w_norm <- ifel(bin_w_sum > 0, bin_w_eff / bin_w_sum, 0)
  
  # tri-class weights
  nt_w_tame   <- rt * pnt
  nt_w_native <- rn * pnt
  
  nt_w_sum <- sum(nt_w_tame + nt_w_native, na.rm=TRUE)
  nt_w_tame_norm   <- ifel(nt_w_sum > 0, nt_w_tame / nt_w_sum, 0)
  nt_w_native_norm <- ifel(nt_w_sum > 0, nt_w_native / nt_w_sum, 0)
  
  # ----------------------
  # EVIDENCE (POSITIVE ONLY)
  # ----------------------
  
  E_not <-
    sum(bin_w_norm * (b == 0), na.rm=TRUE) +
    sum(pmax(nt_w_tame_norm, nt_w_native_norm) * (nt == 0), na.rm=TRUE)
  
  E_tame <-
    sum(bin_w_norm * (b == 1), na.rm=TRUE) * 0.5 +
    sum(nt_w_tame_norm * (nt == 1), na.rm=TRUE)
  
  E_native <-
    sum(bin_w_norm * (b == 1), na.rm=TRUE) * 0.5 +
    sum(nt_w_native_norm * (nt == 2), na.rm=TRUE)
  
  # ----------------------
  # SOFTMAX NORMALIZATION
  # ----------------------
  
  eps <- 1e-6
  
  E_stack <- c(E_not, E_tame, E_native)
  
  # stabilize before exp
  E_max <- max(E_stack)
  expE <- exp(E_stack - E_max)
  
  # can add temperature parameter (tau) to control sharpness:
  # 1   = sharper classification
  # 2-3 = smoother probabilities
  # 3+  = very conservative
  #tau <- 2
  #expE <- exp((E_stack - E_max) / tau)
  
  S <- sum(expE)
  
  p_not    <- expE[1] / S
  p_tame   <- expE[2] / S
  p_native <- expE[3] / S

# =========================
# RUN PARALLEL TILES
# =========================
results <- future_lapply(seq_along(tiles), process_tile)

print("All tiles completed.")

# =========================
# MOSAIC AND EXPORT
# =========================

# Get tile outputs
files <- Filter(Negate(is.null), results)

if (length(files) == 0) {
  files <- list.files("tiles_out", full.names=TRUE, pattern="\\.tif$")
}

# Ensure files exist
files <- files[file.exists(files)]

if (length(files) == 0) {
  stop("No tile files found for mosaicking")
}

cat("Mosaicking", length(files), "tiles...\n")

# Check layer consistency
#stopifnot(all(sapply(files, function(f) nlyr(rast(f)) == 4)))
stopifnot(all(vapply(files, function(f) nlyr(rast(f)) == 4, logical(1))))

# Mosaic all tiles together
mosaic_out <- mosaic(rast(files), fun="mean")

#rast_list <- lapply(files, rast)
#mosaic_out <- do.call(mosaic, c(rast_list, list(
#  fun="mean"
#)))

# Assign layer names
names(mosaic_out) <- c("p_not", "p_tame", "p_native", "entropy")

# Export raster to TIF
out_file <- "/mnt/tmp/joint_prob_full.tif"

writeRaster(
  mosaic_out,
  out_file,
  overwrite=TRUE,
  wopt=list(
    datatype="FLT4S",
    gdal=c("COMPRESS=ZSTD","TILED=YES","BIGTIFF=YES")
  )
)

aws.s3::put_object(
  file   = out_file,
  object = paste0(out_dir, "/joint_prob_full.tif")
)

# Local
#writeRaster(
#  mosaic_out,
#  "joint_prob_full.tif",
#  overwrite=TRUE,
#  wopt=list(
#    datatype="FLT4S",
#    gdal=c("COMPRESS=ZSTD","TILED=YES","BIGTIFF=YES")))

# upload log to S3
aws.s3::put_object(
  file = log_file_local,
  object = log_file
)

#--------- Possible next steps (examples) ---------------

# Hard classification
class <- which.max(joint_prob) - 1

# Thresholding - high confidence native grass
p_native > 0.8 & entropy < 0.3

# Uncertainty masking - remove unreliable areas
entropy > 0.7 

