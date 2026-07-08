
###################################################################################
##                     GRASSLAND CONSENSUS MAPPING ANALYSIS                      ##
###################################################################################

##----------------------- Alison Long - June 2026 -------------------------

# step 9: use data fusion techniques to create a consensus grassland map 
#         to serve as baseline for modelling

# THIS VERSION IS OPTIMIZED FOR SAGEMAKER
# using ml.m5d.12xlarge instance (48 vCPU, 192 GB RAM, NVMe SSD)

# THIS VERSION IS IMPROVED VERSION OF V2

# This version addresses temporality by:
#   Adding time decay weights to MCC weights - prioritize recent data
#   Adding temporal variability layer

# Temporal versions of script need to have PLI atomic (not grouped)

# Temporally-informed outputs:
#   prob_grass: current best estimate (recent data prioritized)
#   temporal_sd: low unstable this pixel is across time
#       High = dynamic landscape or conflicting years
#       Low = consistent classification
#   uncertainty_total: combined model + data + temporal uncertainty

## ---- Setup ---- ##
# Temp dirs
tmp_root <- file.path(Sys.getenv("HOME"), "r_tmp")
dir.create(tmp_root, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(tmp_root, "terra"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(tmp_root, "data"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(tmp_root, "stacks"), recursive = TRUE, showWarnings = FALSE)

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

# GDAL cache
Sys.setenv(GDAL_CACHEMAX = 4096)

# Force a fresh run instead of silently re-using old scratch outputs
FORCE_REBUILD <- TRUE

## ---- Settings ---- ##

out_dir      <- "s3://sagemaker-gst-stage.sharing/alison-long/tomorrows-prairies/results/consensus_map"
log_dir      <- "s3://sagemaker-gst-stage.sharing/alison-long/tomorrows-prairies/logs"
terra_dir    <- file.path(tmp_root, "terra")
local_dir    <- file.path(tmp_root, "data")
stack_dir    <- file.path(tmp_root, "stacks")

meta_file    <- "s3://sagemaker-gst-stage.sharing/alison-long/tomorrows-prairies/processed/metadata_full_extent_binary_nattam_with_path_nohfi.csv"
ecodist_file <- "s3://sagemaker-gst-stage.sharing/alison-long/tomorrows-prairies/processed/boundaries/national-ecological-framework-for-canada/Ecodistricts/ecodistricts_clean_rast.tif"
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
  stop("Duplicate raster_id values found: ", paste(unique(dupes), collapse = ", "))
}

if (any(is.na(meta$file)) || any(meta$file == "")) {
  stop("Missing file path(s) in metadata for raster_id(s): ",
       paste(meta$raster_id[is.na(meta$file) | meta$file == ""], collapse = ", "))
}

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

#download_files <- function(files_named) {
#  if (length(files_named) == 0) return(setNames(character(0), character(0)))
#  sapply(names(files_named), function(nm) {
#    dest <- file.path(local_dir, paste0(nm, ".tif"))
#    if (!file.exists(dest) || FORCE_REBUILD) {
#      ok <- tryCatch({
#        aws.s3::save_object(files_named[[nm]], file = dest)
#        TRUE
#      }, error = function(e) {
#        log_msg("FAILED to download ", nm, " from ", files_named[[nm]], ": ",
#                conditionMessage(e), level = "ERROR")
#        FALSE
#      })
#      if (!ok) stop("Aborting: raster download failed for ", nm)
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
      
      # crop + mask to pilot/study area
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

## ---- Confirm raster alignment ---- ##
message("Checking alignment...")
log_msg("Checking raster alignment")

check_alignment <- function(files) {
  if (length(files) == 0) return(TRUE)
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

## ---- Temporal settings ---- ##
target_year <- 2025 # use newest dataset or policy-relevant year
lambda <- 0.25  # decay rate (0.1–0.4 typical, 0.1 = slow decay (older data still matters))
                                            #  0.3 = faster decay (focus on recent)

#--------- Build MCC weight rasters (binary stack) -------#
log_msg("Building MCC weight stack (binary)")

# Load accuracy (MCC)
acc_bin <- aws.s3::s3read_using(read_csv, object = acc_bin_file, show_col_types = FALSE) %>%
  mutate(raster_name = paste0(tolower(trimws(raster)), "_bin"))

stopifnot(all(c("raster_name", "Ecds_ID", "MCC") %in% names(acc_bin)))

acc_bin <- acc_bin %>% mutate(Ecds_ID = as.integer(round(Ecds_ID)))

dup_check <- acc_bin %>% count(raster_name, Ecds_ID) %>% filter(n > 1)
if (nrow(dup_check) > 0) {
  stop("Duplicate (raster_name, Ecds_ID) rows in accuracy table acc_bin: ",
       paste(utils::head(paste(dup_check$raster_name, dup_check$Ecds_ID), 10), collapse = "; "))
}

# Convert MCC → weight [0,1]
mcc_to_weight <- function(x) {
  w <- (x + 1) / 2
  w[!is.finite(w)] <- NA
  pmax(w, 1e-4)
}

# Join year information into accuracy table
acc_bin <- acc_bin %>%
  left_join(
    meta %>% select(raster_id, year),
    by = c("raster_name" = "raster_id")
  )

# For now, assign pli_bin to year 2022
#acc_bin <- acc_bin %>%
#  mutate(year = if_else(raster_name == "pli_bin", 2022L, year))

if (any(is.na(acc_bin$year))) {
  stop("Missing year info for some rasters in acc_bin")
}

# Temporal decay function
time_decay <- function(year, target_year, lambda) {
  exp(-lambda * abs(year - target_year))
}

#acc_bin <- acc_bin %>% mutate(weight = mcc_to_weight(MCC))

# Apply MCC + temporal weighting
acc_bin <- acc_bin %>%
  mutate(
    w_mcc = mcc_to_weight(MCC),
    w_time = time_decay(year, target_year, lambda),
    weight = w_mcc * w_time
  )

# check for missing columns or naming mismatch
missing <- setdiff(names(bin_stack), acc_bin$raster_name)
if (length(missing) > 0) {
  stop("Missing MCC values for: ", paste(missing, collapse = ", "))
}

# Load rasterized ecodistricts 
#ecodist_r <- rast(ecodist_file)
#ecodist_r <- as.int(ecodist_r) 

log_msg("Loading clipped ecodistrict raster")

ecodist_r <- rast(ecodist_file)

ecodist_r <- crop(ecodist_r, sa)
ecodist_r <- mask(ecodist_r, sa)

ecodist_r <- as.int(ecodist_r)

# check ecodistrict raster values vs IDs
ecodist_freq <- freq(ecodist_r)
log_msg("Ecodistrict raster value/cell counts (head): ",
        paste(utils::capture.output(print(utils::head(ecodist_freq))), collapse = " | "))

ids_r <- unique(ecodist_freq[, "value"])
ids_tbl <- unique(acc_bin$Ecds_ID)

length(missing_ids)

acc_bin %>%
  filter(Ecds_ID %in% missing_ids)

ecodist_freq %>%
  filter(value %in% missing_ids)

missing_ids <- setdiff(ids_r, ids_tbl)
if (length(missing_ids) > 0) {
  stop("Ecodistrict IDs missing from accuracy table: ",
       paste(utils::head(missing_ids), collapse = ", "))
}

# For now, fill MCC for missing ecodistricts
if (length(missing_ids) > 0) {
  
  log_msg(
    "Adding missing ecodistricts with neutral MCC=0: ",
    paste(missing_ids, collapse = ", ")
  )
  
  # Create a row for every raster × missing ecodistrict combination
  missing_rows <- expand.grid(
    raster_name = unique(acc_bin$raster_name),
    Ecds_ID = missing_ids,
    stringsAsFactors = FALSE
  ) %>%
    as_tibble() %>%
    left_join(
      meta %>%
        select(raster_id, year),
      by = c("raster_name" = "raster_id")
    ) %>%
    mutate(
      MCC = 0
    )
  
  # Append to accuracy table
  acc_bin <- bind_rows(acc_bin, missing_rows)
}

# Build weight stack layer-by-layer
weight_files <- lapply(names(bin_stack), function(nm) {
  log_msg("Rasterizing weight for: ", nm)
  
  lookup_bin <- acc_bin %>%
    filter(raster_name == nm) %>%
    select(Ecds_ID, weight) %>%
    as.matrix()
  
  if (nrow(lookup_bin) == 0) stop("No weights found for raster: ", nm)
  
  out <- file.path(terra_dir, paste0("w_", nm, ".tif"))
  if (!file.exists(out) || FORCE_REBUILD) {
    w <- classify(ecodist_r, rcl = lookup_bin)
    writeRaster(w, out, overwrite = TRUE, wopt = list(datatype = "FLT4S"))
  }
  out
})

weight_stack <- rast(unlist(weight_files))
names(weight_stack) <- names(bin_stack)

stopifnot(identical(names(weight_stack), names(bin_stack)))

#--------- MCC-weighted positive evidence -------#
log_msg("Running Strategy A (MCC-weighted fusion)")

present <- !is.na(bin_stack)
pos_mask <- bin_stack == 1
neg_mask <- bin_stack == 0

# Zero-out weights where no data
weight_present <- weight_stack * present

# Numerator
numerator_A <- app(
  pos_mask * weight_present,
  fun = "sum", na.rm = TRUE,
  filename = file.path(terra_dir, "numerator_A.tif"),
  overwrite = TRUE
)

# Denominator
denominator_A <- app(
  weight_present,
  fun = "sum", na.rm = TRUE,
  filename = file.path(terra_dir, "denominator_A.tif"),
  overwrite = TRUE
)

# Probability
prob_A <- numerator_A / denominator_A
prob_A[denominator_A == 0] <- NA

prob_A <- clamp(
  prob_A, 0, 1,
  filename = file.path(terra_dir, "prob_A.tif"),
  overwrite = TRUE
)

#--------- Temporal variability ---------#

log_msg("Computing temporal variability")

# Weighted mean (already computed as prob_A)
mu <- prob_A

# Weighted disagreement across sources (not just temporal)
temporal_var_bin <- app(
  c(bin_stack, weight_present, mu),
  fun = function(x) {
    n <- (length(x) - 1) / 2
    y <- x[1:n]
    w <- x[(n+1):(2*n)]
    m <- x[length(x)]
    
    valid <- !is.na(y) & !is.na(w) & w > 0
    
    if (!any(valid)) return(NA)
    
    sum(w[valid] * (y[valid] - m)^2) / sum(w[valid])
  },
  filename = file.path(terra_dir, "temporal_variance_bin.tif"),
  overwrite = TRUE
)

# Optional: sqrt for interpretability (std dev)
temporal_sd <- sqrt(temporal_var_bin)

writeRaster(
  temporal_sd,
  file.path(terra_dir, "temporal_sd.tif"),
  overwrite = TRUE
)

#---- Diagnostics ----
num_contributing <- app(present, fun = "sum", na.rm = TRUE,
                        filename = file.path(terra_dir, "num_contributing.tif"), overwrite = TRUE)
num_positive <- app(pos_mask, fun = "sum", na.rm = TRUE,
                    filename = file.path(terra_dir, "num_positive.tif"), overwrite = TRUE)
num_negative <- app(neg_mask, fun = "sum", na.rm = TRUE,
                    filename = file.path(terra_dir, "num_negative.tif"), overwrite = TRUE)

sum_weight_present <- denominator_A

#---- Uncertainty ----
uncertainty_prob_A <- prob_A * (1 - prob_A)

uncertainty_data_A <- app(
  sum_weight_present,
  fun = function(x) 1 / (1 + x),
  filename = file.path(terra_dir, "uncertainty_data.tif"),
  overwrite = TRUE
)

#uncertainty_composite_A <- uncertainty_prob_A * (1 + uncertainty_data_A)

#writeRaster(uncertainty_composite_A, file.path(terra_dir, "uncertainty_composite_A.tif"),
#            overwrite = TRUE, wopt = list(datatype = "FLT4S"))

# Normalize temporal SD (optional but recommended)
# Can flatten signal/compress values heavily - meaningful differences above 1 are dampened
temporal_sd_norm <- temporal_sd / (temporal_sd + 1)

uncertainty_composite_A <- (
  uncertainty_prob_A +
    uncertainty_data_A +
    temporal_sd_norm
) / 3

writeRaster(
  uncertainty_composite_A,
  file.path(terra_dir, "uncertainty_total.tif"),
  overwrite = TRUE
)

# Minimum support
min_sources <- 2

prob_A_minN <- ifel(
  num_contributing >= min_sources,
  prob_A, NA,
  filename = file.path(terra_dir, "prob_A_minN.tif"),
  overwrite = TRUE
)

#--------- Tri-class (MCC + F1 weights) ---------#

if (has_nt) {
  
  log_msg("Running tri-class fusion")
  
  acc_nt <- aws.s3::s3read_using(read_csv, object = acc_nt_file, show_col_types = FALSE) %>%
    mutate(raster_name = paste0(tolower(trimws(raster)), "_nt"))
  
  acc_nt <- acc_nt %>%
    left_join(
      meta %>%
        select(raster_id, year),
      by = c("raster_name" = "raster_id")
    )
  
  if (any(is.na(acc_nt$year))) {
    stop("Missing year info for some rasters in acc_nt")
  }
  
  req <- c("raster_name", "Ecds_ID", "MCC_Multiclass", "F1_1", "F1_2")
  stopifnot(all(req %in% names(acc_nt)))
  
  acc_nt <- acc_nt %>% mutate(Ecds_ID = as.integer(round(Ecds_ID)))
  
  dup_check_nt <- acc_nt %>% count(raster_name, Ecds_ID) %>% filter(n > 1)
  if (nrow(dup_check_nt) > 0) {
    stop("Duplicate (raster_name, Ecds_ID) rows in accuracy table acc_nt: ",
         paste(utils::head(paste(dup_check_nt$raster_name, dup_check_nt$Ecds_ID), 10), collapse = "; "))
  }
  
  acc_nt <- acc_nt %>%
    mutate(
      w_time = time_decay(year, target_year, lambda),
      
      W_native = pmax((MCC_Multiclass + 1) / 2 * F1_2 * w_time, 0),
      W_tame   = pmax((MCC_Multiclass + 1) / 2 * F1_1 * w_time, 0)
    )
  
  missing_nt <- setdiff(names(nt_stack), acc_nt$raster_name)
  if (length(missing_nt) > 0) {
    stop("Missing tri-class accuracy values for: ", paste(missing_nt, collapse = ", "))
  }
  
  native_files <- lapply(names(nt_stack), function(nm) {
    lookup_native <- acc_nt %>% filter(raster_name == nm) %>% select(Ecds_ID, W_native) %>% as.matrix()
    if (nrow(lookup_native) == 0) stop("No native weights for raster: ", nm)
    out <- file.path(terra_dir, paste0("nt_native_", nm, ".tif"))
    if (!file.exists(out) || FORCE_REBUILD) {
      w <- classify(ecodist_r, rcl = lookup_native)
      writeRaster(w, out, overwrite = TRUE, wopt = list(datatype = "FLT4S"))
    }
    out
  })
  
  tame_files <- lapply(names(nt_stack), function(nm) {
    lookup_tame <- acc_nt %>% filter(raster_name == nm) %>% select(Ecds_ID, W_tame) %>% as.matrix()
    if (nrow(lookup_tame) == 0) stop("No tame weights for raster: ", nm)
    out <- file.path(terra_dir, paste0("nt_tame_", nm, ".tif"))
    if (!file.exists(out) || FORCE_REBUILD) {
      w <- classify(ecodist_r, rcl = lookup_tame)
      writeRaster(w, out, overwrite = TRUE, wopt = list(datatype = "FLT4S"))
    }
    out
  })
  
  metric_native_stack <- rast(unlist(native_files))
  metric_tame_stack   <- rast(unlist(tame_files))
  
  stopifnot(nlyr(metric_native_stack) == nlyr(nt_stack))
  stopifnot(nlyr(metric_tame_stack) == nlyr(nt_stack))
  
  names(metric_native_stack) <- names(nt_stack)
  names(metric_tame_stack)   <- names(nt_stack)
  
  stopifnot(identical(names(metric_native_stack), names(nt_stack)))
  stopifnot(identical(names(metric_tame_stack), names(nt_stack)))
  
  native_mask <- nt_stack == 2
  tame_mask <- nt_stack == 1
  
  native_votes <- app(
    native_mask * metric_native_stack,
    fun = "sum", na.rm = TRUE,
    filename = file.path(terra_dir, "native_votes.tif"), overwrite = TRUE
  )
  
  tame_votes <- app(
    tame_mask * metric_tame_stack,
    fun = "sum", na.rm = TRUE,
    filename = file.path(terra_dir, "tame_votes.tif"), overwrite = TRUE
  )
  
  vote_total <- native_votes + tame_votes
  
  share_native <- native_votes / vote_total
  share_native[vote_total == 0] <- NA
  share_native <- clamp(share_native, 0, 1)
  share_tame <- 1 - share_native
  
} else {
  
  log_msg("No tri-class rasters present this run — native/tame shares set to NA", level = "WARN")
  
  share_native <- prob_A; values(share_native) <- NA
  share_tame <- prob_A; values(share_tame) <- NA
  
}

# Compute temporal variability
log_msg("Computing temporal variability (tri-class)")

tri_present <- !is.na(nt_stack)

# Total weights (native + tame)
metric_total_stack <- metric_native_stack + metric_tame_stack
weight_present_nt <- metric_total_stack * tri_present

# Convert class to numeric signal:
# 2 = native, 1 = tame, NA otherwise
#signal_nt <- classify(nt_stack,
#                      rcl = matrix(c(
#                        1,1,
#                        2,2
#                      ), ncol=2, byrow=TRUE))

signal_nt <- nt_stack
signal_nt[!(nt_stack %in% c(1,2))] <- NA

# Weighted mean signal
mean_signal <- app(
  c(signal_nt, weight_present_nt),
  fun = function(x) {
    n <- length(x) / 2
    y <- x[1:n]
    w <- x[(n+1):(2*n)]
    
    valid <- !is.na(y) & !is.na(w) & w > 0
    
    if (!any(valid)) return(NA)
    
    sum(y[valid] * w[valid]) / sum(w[valid])
  },
  filename=file.path(terra_dir,"tri_mean_signal.tif"),
  overwrite=TRUE
)

# Weighted variance
temporal_var_nt <- app(
  c(signal_nt, weight_present_nt, mean_signal),
  fun = function(x) {
    n <- (length(x) - 1) / 2
    y <- x[1:n]
    w <- x[(n+1):(2*n)]
    m <- x[length(x)]
    
    valid <- !is.na(y) & !is.na(w) & w > 0
    
    if (!any(valid) || is.na(m)) return(NA)
    
    sum(w[valid] * (y[valid] - m)^2) / sum(w[valid])
  },
  filename=file.path(terra_dir,"temporal_var_nt.tif"),
  overwrite=TRUE
)

temporal_sd_nt <- sqrt(temporal_var_nt)

writeRaster(
  temporal_sd_nt,
  file.path(terra_dir,"temporal_sd_nt.tif"),
  overwrite=TRUE
)

#--------- Final probabilities -----------#
log_msg("Combining binary + tri-class outputs")

# Probabilities
p_grass <- prob_A

p_native <- p_grass * share_native # assumes "is it grass" evidence (bin_stack) and "is it native vs tame" evidence (nt_stack) are independent and always consistent -
p_tame   <- p_grass * share_tame   # potential disagreement is not surfaced/resolved - produces NA or silently-inconsistent output
p_not    <- 1 - p_grass

log_msg("Probabilities  calculated")

# Hard classification with tiebreak rule
log_msg("Running hard classification")
classify_with_tiebreak <- function(p_not, p_tame, p_native) {
  eps <- 1e-9
  p_native_biased <- p_native + 2 * eps
  p_tame_biased <- p_tame + eps
  p_not_biased <- p_not
  
  p_stack_biased <- c(p_not_biased, p_tame_biased, p_native_biased)
  which.max(p_stack_biased)
}

class_map <- classify_with_tiebreak(p_not, p_tame, p_native)
# class_map: 1 = not, 2 = tame, 3 = native (which.max layer index)

rcl <- matrix(c(1, 0, 2, 1, 3, 2), ncol = 2, byrow = TRUE)

# For max-probability / uncertainty reporting, use the unbiased stack
log_msg("Stacking unbiased rasters")
p_stack <- c(p_not, p_tame, p_native)
max_p <- app(p_stack, fun = "max", na.rm = TRUE,
             filename = file.path(terra_dir, "maxp.tif"), overwrite = TRUE)

#----------- Uncertainty ----------#
log_msg("Calculating uncertainty")
uncertainty <- 1 - max_p

# Combine binary + tri-class temporal uncertainty

# Normalize tri-class SD
log_msg("Normalizing tri-class SD")
temporal_sd_nt_norm <- temporal_sd_nt / (temporal_sd_nt + 1)

# Combine all uncertainty sources
log_msg("Combining all uncertainty sources")
uncertainty_total <- (
  uncertainty +
    temporal_sd_nt_norm
) / 2

writeRaster(
  uncertainty_total,
  file.path(terra_dir,"uncertainty_total_combined.tif"), # need to export!
  overwrite=TRUE
)

#------- Write outputs to S3 ----------#
log_msg("Writing outputs")

# Define wopt
wopt <- list(
  datatype="FLT4S",
  gdal=c("COMPRESS=ZSTD","TILED=YES","BIGTIFF=YES"))

wopt_int <- list(
  datatype="INT2U",
  gdal=c("COMPRESS=ZSTD","TILED=YES","BIGTIFF=YES"))

wopt_class <- list(
  datatype="INT1U",
  gdal=c("COMPRESS=ZSTD","TILED=YES","BIGTIFF=YES"))

class_3cat <- classify(class_map, rcl)

# Write locally first, then save to S3
outputs <- list(
  list(raster = prob_A, local = "prob_A.tif", s3 = "prob_grass_MCC.tif", wopt = wopt),
  list(raster = prob_A_minN, local = "prob_grass_minN.tif", s3 = "prob_grass_minN.tif", wopt = wopt),
  list(raster = p_native, local = "p_native.tif", s3 = "p_native.tif", wopt = wopt),
  list(raster = p_tame, local = "p_tame.tif", s3 = "p_tame.tif", wopt = wopt),
  list(raster = p_not, local = "p_not_grass.tif", s3 = "p_not_grass.tif", wopt = wopt),
  list(raster = uncertainty, local = "uncertainty.tif", s3 = "uncertainty.tif", wopt = wopt),
  list(raster = uncertainty_composite_A,local = "uncertainty_composite.tif", s3 = "uncertainty_composite.tif", wopt = wopt),
  list(raster = num_contributing, local = "num_sources.tif", s3 = "num_sources.tif", wopt = wopt_int),
  list(raster = num_positive, local = "num_positive.tif", s3 = "num_positive.tif", wopt = wopt_int),
  list(raster = num_negative, local = "num_negative.tif", s3 = "num_negative.tif", wopt = wopt_int),
  list(raster = class_3cat, local = "grass_3class.tif", s3 = "grass_3class.tif", wopt = wopt_class)
)

write_and_collect <- function(spec) {
  local_path <- file.path(tmp_root, spec$local)
  writeRaster(spec$raster, local_path, overwrite = TRUE, wopt = spec$wopt)
  spec$local_path <- local_path
  spec
}

outputs <- lapply(outputs, write_and_collect)
log_msg("All output rasters written locally; beginning S3 upload")

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
