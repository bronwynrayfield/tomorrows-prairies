
###################################################################################
##                     GRASSLAND CONSENSUS MAPPING ANALYSIS                      ##
###################################################################################

##----------------------- Alison Long - July 2026 -------------------------

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
tmp_root <- file.path(Sys.getenv("HOME"), "r_tmp")
#dir.create(tmp_root, recursive = TRUE, showWarnings = FALSE)

# Load packages
library(terra)
library(readr)
library(dplyr)
library(aws.s3)

# Set terra options
terraOptions(
  tempdir = file.path(tmp_root, "terra"),
  memfrac = 0.6, # could go down to 0.4
  progress = 1
)

## ---- Settings ---- ##
meta_file  <- "s3://sagemaker-gst-stage.sharing/alison-long/tomorrows-prairies/processed/metadata_full_extent_binary_nattam_with_path.csv"
out_dir    <- "s3://sagemaker-gst-stage.sharing/alison-long/tomorrows-prairies/results/consensus_map"
log_dir    <- "s3://sagemaker-gst-stage.sharing/alison-long/tomorrows-prairies/logs"
local_dir  <- file.path(tmp_root, "rasters")
dir.create(local_dir, recursive = TRUE, showWarnings = FALSE)

## ---- Message log ---- ##
log_file <- file.path(log_dir, paste0(
  "consensus_map_log_",
  format(Sys.time(), "%Y%m%d_%H%M%S"),
  ".log"
))

# Ensure local temp file first 
log_dir_local <- file.path(tmp_root, "logs")
dir.create(log_dir_local, recursive = TRUE, showWarnings = FALSE)

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

## ---- Temporal settings ---- ##
target_year <- 2025 # use newest dataset or policy-relevant year
lambda <- 3  # decay rate (0.1–0.4 typical, 0.1 = slow decay (older data still matters))
#  0.3 = faster decay (focus on recent)

# Assign weights
meta <- meta %>%
  mutate(
    time_weight = exp(
      -abs(year - target_year) / lambda
    )
  )

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

## ---- Create temporal-weight stacks ----##

bin_time <- bin_meta$time_weight
names(bin_time) <- names(bin_local)

stopifnot(
  identical(
    names(bin_time),
    names(bin_stack)
  )
)

#time_bin_stack <- rast(
#  lapply(bin_local, function(f) {
#    r <- rast(f)
#    values(r) <- bin_time[which(bin_local == f)]
#    r
#  })
#)

nt_time <- nt_meta$time_weight
names(nt_time) <- names(nt_local)

stopifnot(
  identical(
    names(nt_time),
    names(nt_stack)
  )
)

#time_nt_stack <- rast(
#  lapply(nt_local, function(f) {
#    r <- rast(f)
#    values(r) <- nt_time[which(nt_local == f)]
#    r
#  })
#)

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

if (length(nt_local) > 0) {
  if (!file.exists(nt_stack_file)) {
    writeRaster(
      rast(nt_local),
      nt_stack_file,
      overwrite = TRUE,
      wopt = list(datatype="INT1U",
                  gdal=c("COMPRESS=ZSTD","TILED=YES","BIGTIFF=YES")))
  }
  
  nt_stack <- rast(nt_stack_file)
  names(nt_stack) <- names(nt_local)
}

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

mcc_bin_stack    <- rast("mcc_bin.tif")
mcc_tame_stack   <- rast("mcc_tame.tif")
mcc_native_stack <- rast("mcc_native.tif")

penalty_bin_stack <- rast("penalty_bin.tif")
penalty_nt_stack  <- rast("penalty_nt.tif")

#out_dir <- "tiles_out"
#dir.create(out_dir, showWarnings = FALSE)

# =========================
# TILE GRID
# =========================
# ~2000–5000 tiles is ideal, 500-1500 for single node
tiles <- makeTiles(sa, n = 1000, overlap = 0)

# Save tiles for reproducibility
writeVector(vect(tiles), "tile_index.gpkg", overwrite=TRUE)

# =========================
# PROCESS FUNCTION
# =========================
tile_dir <- "/mnt/tmp/tiles"
dir.create(tile_dir, showWarnings=FALSE)

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
  
  rbin <- crop(mcc_bin_stack, tile)
  rt   <- crop(mcc_tame_stack, tile)
  rn   <- crop(mcc_native_stack, tile)
  
  pbin <- crop(penalty_bin_stack, tile)
  pnt  <- crop(penalty_nt_stack, tile)
  
  # effective weight = accuracy x redundancy penalty
  #bin_w_eff   <- rbin * pbin * time_bin
  #nt_w_tame   <- rt * pnt * time_nt
  #nt_w_native <- rn * pnt * time_nt
  
  bin_w_eff <- rbin * pbin
  
  for (j in seq_len(nlyr(bin_w_eff))) {
    bin_w_eff[[j]] <- bin_w_eff[[j]] * bin_time[j]
  }
  
  nt_w_tame <- rt * pnt
  nt_w_native <- rn * pnt
  
  for (j in seq_len(nlyr(nt_w_tame))) {
    nt_w_tame[[j]] <- nt_w_tame[[j]] * nt_time[j]
    nt_w_native[[j]] <- nt_w_native[[j]] * nt_time[j]
  }
  
  #NEG_HALF <- 0.5
  
  # ----------------------
  # EVIDENCE COMPUTATION
  # ----------------------
  
  # ----------------------
  # NORMALIZED WEIGHTS
  # ----------------------
  
  # binary weights
  #bin_w_eff <- rbin * pbin
  bin_w_sum <- app(bin_w_eff, sum, na.rm=TRUE)
  bin_w_norm <- ifel(bin_w_sum > 0, bin_w_eff / bin_w_sum, 0)
  
  # tri-class weights
  #nt_w_tame   <- rt * pnt
  #nt_w_native <- rn * pnt
  
  nt_w_sum <- app(nt_w_tame + nt_w_native, sum, na.rm=TRUE)
  nt_w_tame_norm   <- ifel(nt_w_sum > 0, nt_w_tame / nt_w_sum, 0)
  nt_w_native_norm <- ifel(nt_w_sum > 0, nt_w_native / nt_w_sum, 0)
  
  # ----------------------
  # EVIDENCE (POSITIVE ONLY)
  # ----------------------
  
  nt_max <- app(c(nt_w_tame_norm, nt_w_native_norm), max, na.rm=TRUE)
  E_not <-
    app(bin_w_norm * (b == 0), sum, na.rm=TRUE) +
    app(nt_max * (nt == 0), sum, na.rm=TRUE)
  
  E_tame <-
    app(bin_w_norm * (b == 1), sum, na.rm=TRUE) * 0.5 +
    app(nt_w_tame_norm * (nt == 1), sum, na.rm=TRUE)
  
  E_native <-
    app(bin_w_norm * (b == 1), sum, na.rm=TRUE) * 0.5 +
    app(nt_w_native_norm * (nt == 2), sum, na.rm=TRUE)
  
  # ----------------------
  # SOFTMAX NORMALIZATION
  # ----------------------
  
  eps <- 1e-6
  
  E_stack <- c(E_not, E_tame, E_native)
  
  # stabilize before exp
  E_max <- app(E_stack, max)
  expE <- exp(E_stack - E_max)
  
  
  # can add temperature parameter (tau) to control sharpness:
  # 1   = sharper classification
  # 2-3 = smoother probabilities
  # 3+  = very conservative
  #tau <- 2
  #expE <- exp((E_stack - E_max) / tau)
  
  S <- app(expE, sum)
  
  p_not    <- expE[[1]] / S
  p_tame   <- expE[[2]] / S
  p_native <- expE[[3]] / S
  
  # stack outputs
  out <- c(p_not, p_tame, p_native)
  
  entropy <- -(
    p_not * log(p_not + eps) +
      p_tame * log(p_tame + eps) +
      p_native * log(p_native + eps)
  )
  
  out <- c(out, entropy)
  names(out) <- c("p_not", "p_tame", "p_native", "entropy")
  
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

# Check layer consistency
#stopifnot(all(sapply(files, function(f) nlyr(rast(f)) == 4)))
stopifnot(all(vapply(files, function(f) nlyr(rast(f)) == 4, logical(1))))

# Mosaic all tiles together
mosaic_out <- mosaic(rast(files), fun="mean")
# or can try
#mosaic_out <- do.call(mosaic, c(lapply(files, rast), list(fun="mean")))

# Set names
names(mosaic_out) <- c("p_not", "p_tame", "p_native", "entropy")

#------ Derived layers ------#
# Probability stack
prob_stack <- mosaic_out[[c("p_not", "p_tame", "p_native")]]
entropy <- mosaic_out[["entropy"]]

# Confidence = max probability
confidence <- app(prob_stack, max)
names(confidence) <- "confidence"

# Check probabilities sum to ~1
check <- app(prob_stack, sum)
stats <- global(check, c("min", "max", "mean"), na.rm=TRUE)
print(stats)
log_msg(capture.output(stats))

#rast_list <- lapply(files, rast)
#mosaic_out <- do.call(mosaic, c(rast_list, list(
#  fun="mean"
#)))

#------ Hard classification ------#
hard_class <- which.max(prob_stack)

# hard classification (1=not, 2=tame, 3=native)
# mask low-confidence or high-uncertainty pixels
# max_prob: require >0.6 confidence (fairly strict) - avoid overconfident misclassification
# entropy: uncertainty masking - remove unreliable areas
# max entropy for 3 classes ≈ 1.10, 0.8=moderately conservative
# could auto-calibrate thresholds based on data distribution instead of fixed values
hard_class <- ifel(
  (confidence < 0.6) | (entropy > 0.8),
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
  (confidence >= 0.75) & (entropy <= 0.6),
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

year_stack <- rast(period_rasters)

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
  
  yrs <- years
  
  temporal_slope <- app(
    ys,
    function(x) {
      
      if (all(is.na(x))) {
        return(NA)
      }
      
      if (sum(!is.na(x)) < 3) {
        return(NA)
      }
      
      ok <- !is.na(x)
      
      if (sum(ok) < 3) {
        return(NA)
      }
      
      coef(
        lm(
          x[ok] ~ yrs[ok]
        )
      )[2]
    }
  )
  
  names(temporal_slope) <- "temporal_slope"
  
  out <- c(
    temporal_sd,
    temporal_entropy,
    temporal_range,
    temporal_slope
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
  "temporal_range",
  "temporal_slope"
)

# Extract layers
temporal_sd <- temporal_out[["temporal_sd"]]

temporal_entropy <- temporal_out[["temporal_entropy"]]

temporal_range <- temporal_out[["temporal_range"]]

temporal_slope <- temporal_out[["temporal_slope"]]

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

# Visualize confidence vs entropy
# to understand whether masking rules are too strict/loose

# Sample ~1e6 pixels max (adjust if needed)
set.seed(42)

vals <- spatSample(
  c(confidence, entropy),
  size = 1e6,
  method = "random",
  na.rm = TRUE,
  as.df = TRUE
)

colnames(vals) <- c("confidence", "entropy")

# Confidence vs entropy plot
# bottom-right = high confidence, low entropy (best pixels)
# top-left = low confidence, high entropy (worst pixels)
# middle cloud = uncertain transitions
# Thresholds appear as: blue line: confidence = 0.6, red line: entropy = 0.8
library(ggplot2)

ggplot(vals, aes(x = confidence, y = entropy)) +
  stat_binhex(bins = 60) +
  scale_fill_viridis_c(trans = "log10") +
  geom_vline(xintercept = 0.6, color = "blue", linetype = "dashed") +
  geom_hline(yintercept = 0.8, color = "red", linetype = "dashed") +
  labs(
    title = "Confidence vs Entropy",
    x = "Confidence (max probability)",
    y = "Entropy",
    fill = "Pixel count"
  ) +
  theme_minimal()

# Density plot (smooth view)
ggplot(vals, aes(x = confidence, y = entropy)) +
  geom_density_2d_filled() +
  geom_vline(xintercept = 0.6, color = "blue", linetype = "dashed") +
  geom_hline(yintercept = 0.8, color = "red", linetype = "dashed") +
  labs(
    title = "Density of Confidence vs Entropy"
  ) +
  theme_minimal()

# Confidence distribution
# If most points cluster around confidence > 0.7 then threshold is appropriate
ggplot(vals, aes(x = confidence)) +
  geom_histogram(bins = 50, fill = "steelblue") +
  geom_vline(xintercept = 0.6, color = "red", linetype = "dashed") +
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
  confidence,
  main = "Confidence (Max Probability)",
  col = viridis::viridis(100)
)

# Binary "uncertain areas" map
uncertain_mask <- ifel(
  (confidence < 0.6) | (entropy > 0.8),
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
uncertainty_score <- (1 - confidence) + (entropy / 1.1)

plot(
  uncertainty_score,
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

# Combined raster 
final_stack <- c(mosaic_out, confidence, hard_class, high_conf_class,
                  temporal_sd, temporal_entropy, temporal_range, temporal_slope)

names(final_stack) <- c(
  "p_not", "p_tame", "p_native", "entropy",
  "confidence",
  "class", "high_conf_class",
  "temporal_sd",
  "temporal_entropy",
  "temporal_range",
  "temporal_slope")

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

# Temporal slope
writeRaster(
  temporal_slope,
  "/mnt/tmp/temporal_slope.tif",
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
  file = "/mnt/tmp/temporal_slope.tif",
  object = paste0(out_dir, "/temporal_slope.tif")
)

#--------- Possible next steps (examples) ---------------

# Thresholding - high confidence native grass
p_native > 0.8 & entropy < 0.3