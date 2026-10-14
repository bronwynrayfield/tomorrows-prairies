
###################################################################################
##          ANALYSIS: SPATIAL AGREEMENT (NATIVE VS TAME VS NON-GRASS)            ##
###################################################################################

##----------------------- Alison Long - May 2026 -------------------------

# step 6: calculate pairwise spatial agreement metrics to assess similarities and 
#         differences among existing grassland datasets. This version handles all 
#         grassland data products regardless of extent and calculates metrics at
#         the ecodistrict level (in addition to global).

# NEED TO UPDATE TO PLI SEP SO METADATA IS CORRECT

# CROP TO A STUDY AREA?

# For instance type, best option is r6i.4xlarge (budget option r6i.2xlarge) 
# ran using ml.m5d.8xlarge
# Prefer instances with NVMe local storage

# if unable to load terra after R updates:
system("ldd ~/R/x86_64-pc-linux-gnu-library/4.4/terra/libs/terra.so")
system("find / -name 'libproj.so*' 2>/dev/null")
system("find / -name 'libgdal.so*' 2>/dev/null")
system("find / -name 'libgeos_c.so*' 2>/dev/null")

Sys.getenv("CONDA_PREFIX")

system("ls $CONDA_PREFIX/lib | grep -E 'proj|gdal|geos'")

system("which gdalinfo")
system("gdalinfo --version")
system("which proj")
system("proj")

# then in terminal:

# sudo apt-get update

# then:
# sudo apt-get install -y \
# gdal-bin \
# libgdal-dev \
# proj-bin \
# libproj-dev \
# libgeos-dev

# then reinstall terra

# then:

# Sys.setenv(
#   AWS_ACCESS_KEY_ID     = creds$AccessKeyId,
#   AWS_SECRET_ACCESS_KEY = creds$SecretAccessKey,
#   AWS_SESSION_TOKEN     = creds$Token
# )

# aws.signature::locate_credentials()

## ---- Setup ---- ##
# Temp dirs
tmp_root <- file.path(Sys.getenv("HOME"), "r_tmp")
#dir.create(tmp_root, recursive = TRUE, showWarnings = FALSE)

# Load packages
library(terra)
#library(future.apply)  # running sequentially to avoid disk contention
library(readr)
library(dplyr)
library(aws.s3)

# Set terra options
terraOptions(
  tempdir = file.path(tmp_root, "terra"),
  memfrac = 0.6, # could go down to 0.4
  progress = 1
)

# Force a fresh run instead of silently re-using old scratch outputs
FORCE_REBUILD <- TRUE

## ---- Parallel (light use only) ---- ##
#n_cores <- min(6, parallel::detectCores(logical = FALSE) - 2)
#future::plan(sequential)  # Avoids disk contention

## ---- Settings ---- ##
#meta_file  <- "s3://sagemaker-gst-stage.sharing/alison-long/tomorrows-prairies/processed/metadata_full_extent_nattam_with_path.csv"
meta_file  <- "s3://sagemaker-gst-stage.sharing/alison-long/tomorrows-prairies/processed/metadata_full_extent_nattam_with_path_plicombinedandseparated.csv"
ecodist_file <- "s3://sagemaker-gst-stage.sharing/alison-long/tomorrows-prairies/processed/boundaries/national-ecological-framework-for-canada/Ecodistricts/ecodistricts_clean_rast_Oct2026.tif"
ccs_file     <- "s3://sagemaker-gst-stage.sharing/alison-long/tomorrows-prairies/processed/boundaries/census-consolidated-subdivisions/ccs_aligned_class_sa_rast.tif"
pilot_file   <- "s3://sagemaker-gst-stage.sharing/alison-long/tomorrows-prairies/processed/pilot-area/pilot_area_aligned_cnvc_Jan2026.tif"
#sa_file      <- "s3://sagemaker-gst-stage.sharing/alison-long/tomorrows-prairies/processed/study-area/modelling_extent_aligned_cnvc_Jan2026.tif"
sa_file      <- "s3://sagemaker-gst-stage.sharing/alison-long/tomorrows-prairies/processed/study-area/temperate_scope_aligned_cnvc_Apr2026.tif"
out_dir    <- "s3://sagemaker-gst-stage.sharing/alison-long/tomorrows-prairies/results/spatial_agreement_analysis"
log_dir    <- "s3://sagemaker-gst-stage.sharing/alison-long/tomorrows-prairies/logs"
local_dir <- file.path(tmp_root, "data")
dir.create(local_dir, showWarnings = FALSE, recursive = TRUE)

## ---- Message log ---- ##
log_file <- file.path(log_dir, paste0(
  "spatial_agreement_log_",
  format(Sys.time(), "%Y%m%d_%H%M%S"),
  ".log"
))

# Ensure local temp file first (upload to S3 at end)
log_file_local <- file.path(
  Sys.getenv("HOME"),
  paste0("spatial_agreement_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log")
)

log_con <- file(log_file_local, open = "wt")

sink(log_con, split = TRUE)          # captures console output

on.exit({
  sink()
  close(log_con)
}, add = TRUE)

log_msg <- function(..., level = "INFO") {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  msg <- paste0(
    "[", level, "] ",
    ts,
    " : ",
    paste(..., collapse = "")
  )
  
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
#log_msg("tidyverse version: ", as.character(packageVersion("tidyverse")))
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

log_msg("Starting spatial agreement analysis")
#log_msg("Parallel workers detected: ", n_cores)

## ---- Load metadata ---- ##
log_msg("Reading metadata from: ", meta_file)

meta <- aws.s3::s3read_using(read_csv, object = meta_file, show_col_types = FALSE) %>%
  rename(file = file_s3)

rfiles_s3 <- setNames(meta$file, meta$raster)

stopifnot(!anyDuplicated(meta$raster))

log_msg("Metadata rows loaded: ", nrow(meta))

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

#rfiles_local <- sapply(names(rfiles_s3), function(nm) {
#  dest <- file.path(local_dir, paste0(nm, ".tif"))
#  if (!file.exists(dest)) {
#    aws.s3::save_object(rfiles_s3[[nm]], file = dest)
#  }
#  dest
#}, USE.NAMES = TRUE)

# rfiles_local <- sapply(names(rfiles_s3), function(nm) {
#   
#   s3_path <- rfiles_s3[[nm]]
#   
#   s3_clean <- sub("^s3://", "", s3_path)
#   parts <- strsplit(s3_clean, "/", fixed = TRUE)[[1]]
#   
#   bucket <- parts[1]
#   key <- paste(parts[-1], collapse = "/")
#   
#   dest <- file.path(local_dir, paste0(nm, ".tif"))
#   
#   aws.s3::save_object(object = key, bucket = bucket, file = dest)
#   
#   dest
# }, USE.NAMES = TRUE)
# 
# log_msg("All rasters downloaded locally")

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

rfiles_local <- download_files(rfiles_s3)

# Check file size to ensure correct download
bad_files <- names(rfiles_local)[file.size(rfiles_local) < 10000]
if (length(bad_files) > 0) {
  stop(paste("These rasters failed to download correctly:",
             paste(bad_files, collapse = ", ")))
}

#file.info(rfiles_local)
#rfiles_s3[c("aci", "ngi", "cgr")]

#unlink(local_dir, recursive = TRUE)
#dir.create(local_dir, recursive = TRUE)

## ---- Confirm raster alignment ---- ##
message("Checking alignment...")
log_msg("Checking raster alignment")

template <- rast(rfiles_local[[1]])

aligned_ok <- all(sapply(rfiles_local, function(f) {
  r <- rast(f)
  compareGeom(r, template, stopOnError = FALSE)
}))

if (!aligned_ok) {
  stop("Rasters are NOT perfectly aligned. You must align them.")
}

message("Rasters confirmed aligned — skipping resampling")
log_msg("Rasters confirmed aligned — skipping resampling")

## ---- Align rasters to a common grid (if needed) ---- ##
#message("Aligning rasters...")

#template <- rast(rfiles_local[[1]])

#aligned_files <- lapply(names(rfiles_local), function(nm) {
#  f <- rfiles_local[[nm]]
#  r <- rast(f)

#  if (!compareGeom(r, template, stopOnError = FALSE)) {
#    r <- resample(r, template, method = "near")
#  }

#  out <- file.path(local_dir, paste0(nm, "_aligned.tif"))
#  writeRaster(r, out, overwrite = TRUE)
#  out
#})

#names(aligned_files) <- names(rfiles_local)

## ---- Load ecodistrict raster ---- ##
log_msg("Loading pilot-clipped ecodistrict raster")

ecodist_r <- rast(ecodist_file)

if (!compareGeom(template, pilot, stopOnError = FALSE)) {
  stop("Pilot raster is not aligned with template raster")
}

ecodist_r <- crop(ecodist_r, pilot)
ecodist_r <- mask(ecodist_r, pilot)

ecodist_rast <- as.int(ecodist_r)

# Ensure alignment
if (!compareGeom(ecodist_rast, template, stopOnError = FALSE)) {
  stop("Ecodistrict raster is not aligned with input rasters")
}

log_msg("Loaded ecodistrict raster")

## ---- Load census consolidated subdivision raster ---- ##
log_msg("Loading pilot-clipped census consolidated subdivision raster")

ccs_r <- rast(ccs_file)

if (!compareGeom(template, pilot, stopOnError = FALSE)) {
  stop("Pilot raster is not aligned with template raster")
}

ccs_r <- crop(ccs_r, pilot)
ccs_r <- mask(ccs_r, pilot)

ccs_rast <- as.int(ccs_r)

# Ensure alignment
if (!compareGeom(ccs_rast, template, stopOnError = FALSE)) {
  stop("Census consolidated subdivision raster is not aligned with input rasters")
}

log_msg("Loaded census consolidated subdivision raster")

## ---- Build raster stack ---- ##
log_msg("Building raster stack")
r_stack <- rast(rfiles_local)
n_rasters <- nlyr(r_stack)
names(r_stack) <- names(rfiles_local)
ras_names <- names(r_stack)

## ---- Initialize ecodistricts ---- ##

# Get unique district IDs (could also delay until first block)
#district_vals <- as.integer(sort(unique(na.omit(values(ecodist_rast, mat = FALSE)))))
freq_table <- freq(ecodist_rast)
district_vals <- as.integer(freq_table$value)

if (anyDuplicated(district_vals)) {
  stop("Duplicate district IDs found")
}

# initialize empty set
#district_vals <- integer(0)

n_districts <- length(district_vals)

# Map district IDs → index
district_index <- setNames(seq_along(district_vals), district_vals)

## ---- Initialize CCS ---- ##

# Get unique IDs (could also delay until first block)
#district_vals <- as.integer(sort(unique(na.omit(values(ecodist_rast, mat = FALSE)))))
freq_table <- freq(ccs_rast)
ccs_vals <- as.integer(freq_table$value)

if (anyDuplicated(ccs_vals)) {
  stop("Duplicate CCS IDs found")
}

# initialize empty set
#district_vals <- integer(0)

n_ccs <- length(ccs_vals)

# Map CCS IDs → index
ccs_index <- setNames(seq_along(ccs_vals), ccs_vals)

## ---- Initialize count arrays ---- ##
classes <- c(0,1,2)
#stopifnot(all(values(r_stack) %in% classes | is.na(values(r_stack))))

n_classes <- length(classes)

counts <- array(
  0,
  dim = c(n_rasters, n_rasters, n_classes, n_classes),
  dimnames = list(ras_names, ras_names, classes, classes)
)

counts_ecodist <- array(
  0,
  dim = c(n_rasters, n_rasters, n_classes, n_classes, n_districts)
)

counts_ccs <- array(
  0,
  dim = c(n_rasters, n_rasters, n_classes, n_classes, n_ccs)
)

# Count pixels where ecodistrict or ccs = NA
na_district_pixels <- 0
na_ccs_pixels <- 0

## ---- Block processing ---- ##
message("Starting blockwise processing...")
log_msg("Starting blockwise processing")

# Block size determined by terraOptions(memfrac=...), available RAM, and number of layers
# Can force smaller chunks for lower RAM per iteration (increase # of blocks)
bs <- blocks(r_stack) #, n = 10000) # Try 8-20k
total_blocks <- bs$n

log_msg(
  "Block count: ", bs$n,
  "; min rows=", min(bs$nrows),
  "; max rows=", max(bs$nrows),
  "; mean rows=", round(mean(bs$nrows), 1)
)

# Derive block size from terra's default
#bs_default <- blocks(r_stack)

#mean_rows <- mean(bs_default$nrows)

# shrink blocks moderately (not drastically)
#target_rows <- floor(mean_rows / 3)

#rows <- seq(1, nrow(r_stack), by = target_rows)

#bs <- list(
#  row = rows,
#  nrows = pmin(target_rows, nrow(r_stack) - rows + 1),
#  n = length(rows)
#)

#total_blocks <- bs$n

# OPEN FILE CONNECTIONS
readStart(r_stack)
readStart(ecodist_rast)
readStart(ccs_rast)

on.exit({
  try(readStop(r_stack), silent = TRUE)
  try(readStop(ecodist_rast), silent = TRUE)
  try(readStop(ccs_rast), silent = TRUE)
}, add = TRUE)

for (i in 1:total_blocks) {
  
  message(sprintf("Block %d / %d", i, total_blocks))
  
  vals <- readValues(
    r_stack,
    row = bs$row[i],
    nrows = bs$nrows[i],
    mat = TRUE
  )
  
  districts <- readValues(
    ecodist_rast,
    row = bs$row[i],
    nrows = bs$nrows[i],
    mat = FALSE
  )
  
  ccs <- readValues(
    ccs_rast,
    row = bs$row[i],
    nrows = bs$nrows[i],
    mat = FALSE
  )
  
  # Count NA ecodistrict and ccs pixels
  na_district_pixels <- na_district_pixels + sum(is.na(districts))
  na_ccs_pixels <- na_ccs_pixels + sum(is.na(ccs))
  
  # Safe indexing
  districts <- as.integer(districts)
  ccs <- as.integer(ccs)

  # Skip empty block
  if (all(is.na(vals))) next
  
  # Keep vals 
  #valid_rows <- rowSums(!is.na(vals)) > 0 & !is.na(districts)
  valid_rows_global  <- rowSums(!is.na(vals)) > 0
  valid_rows_ecodist <- valid_rows_global & !is.na(districts)
  valid_rows_ccs <- valid_rows_global & !is.na(ccs)
  
  #vals <- vals[valid_rows, , drop = FALSE]
  vals_global <- vals[valid_rows_global, , drop = FALSE]
  vals_dist   <- vals[valid_rows_ecodist, , drop = FALSE]
  vals_ccs <- vals[valid_rows_ccs, , drop = FALSE]
  districts_block <- districts[valid_rows_ecodist]
  ccs_block <- ccs[valid_rows_ccs]
  
  bad_vals <- unique(vals[!is.na(vals) & !(vals %in% classes)])
  if (length(bad_vals) > 0) {
    stop("Unexpected class values detected: ", paste(bad_vals, collapse = ", "))
  }
  
  # Use integers to avoid double precision overhead and cut RAM
  storage.mode(vals) <- "integer"
  
  if (nrow(vals_global) == 0) {
    rm(vals, vals_global, vals_dist, districts, districts_block, ccs, ccs_block, vals_ccs)
    gc()
    next
  }
  
  # Global matrix
  M_global <- vals_global
  
  # District matrix
  M_dist <- vals_dist
  
  M_c <- NULL
  split_idx_eco <- NULL
  split_idx_ccs <- NULL
  
  #split_idx <- split(seq_len(nrow(vals_dist)), districts_block)
  
  ## ---- Global class-wise counting ---- ##
  # For each pair, only count pixels where both rasters are non-NA
  for (ri in 1:(n_rasters - 1)) {
    for (rj in (ri + 1):n_rasters) {
      
      valid <- !is.na(M_global[,ri]) & !is.na(M_global[,rj])
      if (!any(valid)) next
      
      Mi <- M_global[valid, ri]
      Mj <- M_global[valid, rj]
      
      # Build contingency table
      tab <- table(factor(Mi, levels = classes),
                   factor(Mj, levels = classes))
      
      # Add to counts
      counts[ri, rj, , ] <- counts[ri, rj, , ] + tab
      counts[rj, ri, , ] <- counts[rj, ri, , ] + t(tab)
    }
  }
  
  ## ---- Ecodistrict-level class-wise counting ---- ##
  if (nrow(vals_dist) > 0) {
    
    split_idx_eco <- split(seq_len(nrow(vals_dist)), districts_block)
    
    for (d in names(split_idx_eco)) {
      
      rows <- split_idx_eco[[d]]
      d_id <- district_index[as.character(d)]
      if (is.na(d_id)) next
      d_id <- as.integer(d_id)
      
      M_d <- M_dist[rows, , drop = FALSE]
      
      for (ri in 1:(n_rasters - 1)) {
        for (rj in (ri + 1):n_rasters) {
          
          valid <- !is.na(M_d[,ri]) & !is.na(M_d[,rj])
          if (!any(valid)) next
          
          Mi <- M_d[valid, ri]
          Mj <- M_d[valid, rj]
          
          tab <- table(
            factor(Mi, levels = classes),
            factor(Mj, levels = classes)
          )
          
          counts_ecodist[ri, rj, , , d_id] <-
            counts_ecodist[ri, rj, , , d_id] + tab
          
          counts_ecodist[rj, ri, , , d_id] <-
            counts_ecodist[rj, ri, , , d_id] + t(tab)
        }
      }
    }
  }
  
  ## ---- CCS-level class-wise counting ---- ##
  if (nrow(vals_ccs) > 0) {
    
    split_idx_ccs <- split(seq_len(nrow(vals_ccs)), ccs_block)
    
    for (c in names(split_idx_ccs)) {
      
      rows <- split_idx_ccs[[c]]
      
      c_id <- ccs_index[as.character(c)]
      if (is.na(c_id)) next
      
      c_id <- as.integer(c_id)
      
      M_c <- vals_ccs[rows, , drop = FALSE]
      
      for (ri in 1:(n_rasters - 1)) {
        for (rj in (ri + 1):n_rasters) {
          
          valid <- !is.na(M_c[, ri]) &
            !is.na(M_c[, rj])
          
          if (!any(valid)) next
          
          Mi <- M_c[valid, ri]
          Mj <- M_c[valid, rj]
          
          tab <- table(
            factor(Mi, levels = classes),
            factor(Mj, levels = classes)
          )
          
          counts_ccs[ri, rj, , , c_id] <-
            counts_ccs[ri, rj, , , c_id] + tab
          
          counts_ccs[rj, ri, , , c_id] <-
            counts_ccs[rj, ri, , , c_id] + t(tab)
        }
      }
    }
  }
  
  # Force garbage collection to prevent memory accumulation
  rm(vals, vals_global, vals_dist, vals_ccs, districts, districts_block, 
     ccs, ccs_block, M_global, M_dist)
  #if (exists("M_c", inherits = FALSE)) rm(M_c)
  if (exists("split_idx_eco", inherits = FALSE)) rm(split_idx_eco)
  if (exists("split_idx_ccs", inherits = FALSE)) rm(split_idx_ccs)
  gc()
}

# ✅ CLOSE FILE CONNECTIONS
readStop(r_stack)
readStop(ecodist_rast)
readStop(ccs_rast)

message("Finished block processing")
log_msg("Finished block processing")

message("Total pixels where district = NA: ", na_district_pixels)
log_msg("Total pixels where ecodistrict raster is NA: ", na_district_pixels)

message("Total pixels where CCS = NA: ", na_ccs_pixels)
log_msg("Total pixels where CCS raster is NA: ", na_ccs_pixels)

## ---- Check that all classes occur ---- ##

for (k in seq_along(classes)) {
  
  row_total <- sum(counts[, , k, ], na.rm = TRUE)
  col_total <- sum(counts[, , , k], na.rm = TRUE)
  
  message(
    sprintf(
      "Class %s totals: row=%s col=%s",
      classes[k],
      format(row_total, big.mark = ","),
      format(col_total, big.mark = ",")))
  
  if (row_total == 0 && col_total == 0) {
    warning(
      sprintf(
        "Class %s never occurs in any raster",
        classes[k]))
  }
}

## ---- Compute metrics ---- ##
message("Computing metrics...")
log_msg("Computing metrics")

results_global <- list()

for (i in 1:(n_rasters - 1)) {
  for (j in (i + 1):n_rasters) {
    
    C <- counts[i, j, , ]  # 3x3 matrix
    if (any(C < 0)) stop("Negative counts detected")
    N <- sum(C)
    
    if (is.na(N) || N == 0) {
      kappa <- NA
      accuracy <- NA
      jacc <- rep(NA, n_classes)
      jaccard_macro <- NA
      f1_class <- rep(NA, n_classes)
      f1_macro <- NA
      mcc_multi <- NA
    } else {
      
      # Multi-class Kappa
      po <- sum(diag(C)) / N
      
      row_marg <- rowSums(C) / N
      col_marg <- colSums(C) / N
      
      pe <- sum(row_marg * col_marg)
      
      kappa <- ifelse(pe == 1, NA, (po - pe)/(1 - pe))
      
      # Class-specific Jaccard
      jacc <- numeric(n_classes)
      
      for (k in 1:n_classes) {
        TP <- C[k, k]
        FN <- sum(C[k, ]) - TP
        FP <- sum(C[, k]) - TP
        
        denom <- TP + FP + FN
        jacc[k] <- ifelse(denom == 0, NA, TP / denom)
      }
      
      # Macro Jaccard
      jaccard_macro <- if (all(is.na(jacc))) NA_real_ else mean(jacc, na.rm = TRUE)
      
      # Overall accuracy
      accuracy <- sum(diag(C)) / N
      
      # Per-class F1
      f1_class <- numeric(n_classes)
      
      for (k in 1:n_classes) {
        
        TP <- C[k, k]
        FN <- sum(C[k, ]) - TP
        FP <- sum(C[, k]) - TP
        
        precision <- ifelse((TP + FP) == 0, NA, TP / (TP + FP))
        recall    <- ifelse((TP + FN) == 0, NA, TP / (TP + FN))
        
        f1_class[k] <- ifelse(
          is.na(precision) | is.na(recall) | (precision + recall == 0),
          NA,
          2 * precision * recall / (precision + recall)
        )
      }
      
      # Macro F1
      f1_macro <- if (all(is.na(f1_class))) NA_real_ else mean(f1_class, na.rm = TRUE)
      
      # Multi-class MCC (Gorodkin 2004 formulation)
      t_sum <- sum(C)                     # total N
      sum_row_prod <- sum(rowSums(C) * colSums(C))
      trace_C <- sum(diag(C))
      
      numerator <- (trace_C * t_sum) - sum_row_prod
      
      denom_row <- t_sum^2 - sum(colSums(C)^2)
      denom_col <- t_sum^2 - sum(rowSums(C)^2)
      
      denominator <- sqrt(denom_row * denom_col)
      
      mcc_multi <- ifelse(denominator == 0, NA, numerator / denominator)
      
    }
    
    # Define output
    results_global[[length(results_global)+1]] <- tibble(
      region_type = "global",
      region_id = "ALL",
      raster_i = ras_names[i],
      raster_j = ras_names[j],
      c00 = C[1,1], # 9 cells in confusion matrix
      c01 = C[1,2],
      c02 = C[1,3],
      c10 = C[2,1],
      c11 = C[2,2],
      c12 = C[2,3],
      c20 = C[3,1],
      c21 = C[3,2],
      c22 = C[3,3],
      row0 = sum(C[1,]), # pixels class 0 in raster_i
      row1 = sum(C[2,]), # pixels class 1 in raster_i
      row2 = sum(C[3,]), # pixels class 2 in raster_i
      col0 = sum(C[,1]), # pixels class 0 in raster_j
      col1 = sum(C[,2]), # pixels class 1 in raster_j
      col2 = sum(C[,3]), # pixels class 2 in raster_j
      kappa = kappa,
      accuracy = accuracy,
      jaccard_macro = jaccard_macro,
      jaccard_class0 = jacc[1],
      jaccard_class1 = jacc[2],
      jaccard_class2 = jacc[3],
      f1_macro = f1_macro,
      f1_class0 = f1_class[1],
      f1_class1 = f1_class[2],
      f1_class2 = f1_class[3],
      mcc_multi = mcc_multi,
      N = N
    )
  }
}

results_ecodist <- list()

for (d in seq_len(n_districts)) {
  
  d_name <- district_vals[d]
  
  for (i in 1:(n_rasters - 1)) {
    for (j in (i + 1):n_rasters) {
      
      C <- counts_ecodist[i, j, , , d]
      N <- sum(C)
      
      if (is.na(N) || N == 0) {
        kappa <- NA
        accuracy <- NA
        jacc <- rep(NA, n_classes)
        jaccard_macro <- NA
        f1_class <- rep(NA, n_classes)
        f1_macro <- NA
        mcc_multi <- NA
      } else {
        
        po <- sum(diag(C)) / N
        row_marg <- rowSums(C) / N
        col_marg <- colSums(C) / N
        pe <- sum(row_marg * col_marg)
        
        kappa <- ifelse(pe == 1, NA, (po - pe)/(1 - pe))
        
        jacc <- numeric(n_classes)
        for (k in 1:n_classes) {
          TP <- C[k, k]
          FN <- sum(C[k, ]) - TP
          FP <- sum(C[, k]) - TP
          denom <- TP + FP + FN
          jacc[k] <- ifelse(denom == 0, NA, TP / denom)
        }
        
        jaccard_macro <- if (all(is.na(jacc))) NA_real_ else mean(jacc, na.rm = TRUE)
        accuracy <- sum(diag(C)) / N
        
        f1_class <- numeric(n_classes)
        for (k in 1:n_classes) {
          TP <- C[k, k]
          FN <- sum(C[k, ]) - TP
          FP <- sum(C[, k]) - TP
          
          precision <- ifelse((TP+FP)==0, NA, TP/(TP+FP))
          recall    <- ifelse((TP+FN)==0, NA, TP/(TP+FN))
          
          f1_class[k] <- ifelse(
            is.na(precision) | is.na(recall) | (precision+recall)==0,
            NA,
            2*precision*recall/(precision+recall)
          )
        }
        
        f1_macro <- if (all(is.na(f1_class))) NA_real_ else mean(f1_class, na.rm = TRUE)
        
        t_sum <- sum(C)
        sum_row_prod <- sum(rowSums(C) * colSums(C))
        trace_C <- sum(diag(C))
        
        numerator <- (trace_C * t_sum) - sum_row_prod
        
        denom_row <- t_sum^2 - sum(colSums(C)^2)
        denom_col <- t_sum^2 - sum(rowSums(C)^2)
        
        denominator <- sqrt(denom_row * denom_col)
        
        mcc_multi <- ifelse(denominator == 0, NA, numerator / denominator)
      }
      
      results_ecodist[[length(results_ecodist)+1]] <- tibble(
        region_type = "ecodistrict", 
        region_id = as.character(d_name),
        raster_i = ras_names[i],
        raster_j = ras_names[j],
        c00 = C[1,1], # 9 cells in confusion matrix
        c01 = C[1,2],
        c02 = C[1,3],
        c10 = C[2,1],
        c11 = C[2,2],
        c12 = C[2,3],
        c20 = C[3,1],
        c21 = C[3,2],
        c22 = C[3,3],
        row0 = sum(C[1,]), # pixels class 0 in raster_i
        row1 = sum(C[2,]), # pixels class 1 in raster_i
        row2 = sum(C[3,]), # pixels class 2 in raster_i
        col0 = sum(C[,1]), # pixels class 0 in raster_j
        col1 = sum(C[,2]), # pixels class 1 in raster_j
        col2 = sum(C[,3]), # pixels class 2 in raster_j
        kappa = kappa,
        accuracy = accuracy,
        jaccard_macro = jaccard_macro,
        jaccard_class0 = jacc[1],
        jaccard_class1 = jacc[2],
        jaccard_class2 = jacc[3],
        f1_macro = f1_macro,
        f1_class0 = f1_class[1],
        f1_class1 = f1_class[2],
        f1_class2 = f1_class[3],
        mcc_multi = mcc_multi,
        N = N
      )
    }
  }
}

results_ccs <- list()

for (c in seq_len(n_ccs)) {
  
  c_name <- ccs_vals[c]
  
  for (i in 1:(n_rasters - 1)) {
    for (j in (i + 1):n_rasters) {
      
      C <- counts_ccs[i, j, , , c]
      N <- sum(C)
      
      if (is.na(N) || N == 0) {
        kappa <- NA
        accuracy <- NA
        jacc <- rep(NA, n_classes)
        jaccard_macro <- NA
        f1_class <- rep(NA, n_classes)
        f1_macro <- NA
        mcc_multi <- NA
      } else {
        
        po <- sum(diag(C)) / N
        row_marg <- rowSums(C) / N
        col_marg <- colSums(C) / N
        pe <- sum(row_marg * col_marg)
        
        kappa <- ifelse(pe == 1, NA, (po - pe)/(1 - pe))
        
        jacc <- numeric(n_classes)
        for (k in 1:n_classes) {
          TP <- C[k, k]
          FN <- sum(C[k, ]) - TP
          FP <- sum(C[, k]) - TP
          denom <- TP + FP + FN
          jacc[k] <- ifelse(denom == 0, NA, TP / denom)
        }
        
        jaccard_macro <- if (all(is.na(jacc))) NA_real_ else mean(jacc, na.rm = TRUE)
        accuracy <- sum(diag(C)) / N
        
        f1_class <- numeric(n_classes)
        for (k in 1:n_classes) {
          TP <- C[k, k]
          FN <- sum(C[k, ]) - TP
          FP <- sum(C[, k]) - TP
          
          precision <- ifelse((TP+FP)==0, NA, TP/(TP+FP))
          recall    <- ifelse((TP+FN)==0, NA, TP/(TP+FN))
          
          f1_class[k] <- ifelse(
            is.na(precision) | is.na(recall) | (precision+recall)==0,
            NA,
            2*precision*recall/(precision+recall)
          )
        }
        
        f1_macro <- if (all(is.na(f1_class))) NA_real_ else mean(f1_class, na.rm = TRUE)
        
        t_sum <- sum(C)
        sum_row_prod <- sum(rowSums(C) * colSums(C))
        trace_C <- sum(diag(C))
        
        numerator <- (trace_C * t_sum) - sum_row_prod
        
        denom_row <- t_sum^2 - sum(colSums(C)^2)
        denom_col <- t_sum^2 - sum(rowSums(C)^2)
        
        denominator <- sqrt(denom_row * denom_col)
        
        mcc_multi <- ifelse(denominator == 0, NA, numerator / denominator)
      }
      
      results_ccs[[length(results_ccs)+1]] <- tibble(
        region_type = "ccs",
        region_id = as.character(c_name),
        raster_i = ras_names[i],
        raster_j = ras_names[j],
        c00 = C[1,1], # 9 cells in confusion matrix
        c01 = C[1,2],
        c02 = C[1,3],
        c10 = C[2,1],
        c11 = C[2,2],
        c12 = C[2,3],
        c20 = C[3,1],
        c21 = C[3,2],
        c22 = C[3,3],
        row0 = sum(C[1,]), # pixels class 0 in raster_i
        row1 = sum(C[2,]), # pixels class 1 in raster_i
        row2 = sum(C[3,]), # pixels class 2 in raster_i
        col0 = sum(C[,1]), # pixels class 0 in raster_j
        col1 = sum(C[,2]), # pixels class 1 in raster_j
        col2 = sum(C[,3]), # pixels class 2 in raster_j
        kappa = kappa,
        accuracy = accuracy,
        jaccard_macro = jaccard_macro,
        jaccard_class0 = jacc[1],
        jaccard_class1 = jacc[2],
        jaccard_class2 = jacc[3],
        f1_macro = f1_macro,
        f1_class0 = f1_class[1],
        f1_class1 = f1_class[2],
        f1_class2 = f1_class[3],
        mcc_multi = mcc_multi,
        N = N
      )
    }
  }
}

pairwise_long <- bind_rows(results_global, results_ecodist, results_ccs)

log_msg("Done computing metrics")

##------ Final checks ----##

validate_spatial_agreement_multiclass <- function(
    counts,
    counts_ecodist,
    counts_ccs,
    pairwise_long,
    district_vals,
    ccs_vals,
    ras_names,
    classes,
    na_district_pixels = NA,
    na_ccs_pixels = NA,
    tol = 1e-8
) {
  
  message("---- Running multiclass validation checks ----")
  
  ## ---- 0. Basic dimensions ---- ##
  if (length(dim(counts)) != 4) {
    stop("counts must be a 4D array: [raster_i, raster_j, class_i, class_j]")
  }
  
  if (length(dim(counts_ecodist)) != 5) {
    stop("counts_ecodist must be a 5D array: [raster_i, raster_j, class_i, class_j, district]")
  }
  
  if (length(dim(counts_ccs)) != 5) {
    stop("counts_ccs must be a 5D array: [raster_i, raster_j, class_i, class_j, ccs]")
  }
  
  n_rasters <- length(ras_names)
  n_classes <- length(classes)
  n_districts <- length(district_vals)
  n_ccs <- length(ccs_vals)
  
  expected_counts_dim <- c(n_rasters, n_rasters, n_classes, n_classes)
  expected_ecodist_dim <- c(n_rasters, n_rasters, n_classes, n_classes, n_districts)
  expected_ccs_dim <- c(n_rasters, n_rasters, n_classes, n_classes, n_ccs)
  
  if (!all(dim(counts) == expected_counts_dim)) {
    stop(sprintf(
      "counts has wrong dimensions. Expected %s, got %s",
      paste(expected_counts_dim, collapse = " x "),
      paste(dim(counts), collapse = " x ")
    ))
  }
  
  if (!all(dim(counts_ecodist) == expected_ecodist_dim)) {
    stop(sprintf(
      "counts_ecodist has wrong dimensions. Expected %s, got %s",
      paste(expected_ecodist_dim, collapse = " x "),
      paste(dim(counts_ecodist), collapse = " x ")
    ))
  }
  
  if (!all(dim(counts_ccs) == expected_ccs_dim)) {
    stop(sprintf(
      "counts_ccs has wrong dimensions. Expected %s, got %s",
      paste(expected_ccs_dim, collapse = " x "),
      paste(dim(counts_ccs), collapse = " x ")
    ))
  }
  
  message("✓ Dimension checks passed")
  
  ## ---- 1. Count integrity ---- ##
  if (any(is.na(counts))) {
    stop("NA values detected in counts")
  }
  if (any(is.na(counts_ecodist))) {
    stop("NA values detected in counts_ecodist")
  }
  if (any(is.na(counts_ccs))) {
    stop("NA values detected in counts_ccs")
  }
  
  if (any(counts < 0)) {
    stop("Negative values detected in counts")
  }
  if (any(counts_ecodist < 0)) {
    stop("Negative values detected in counts_ecodist")
  }
  if (any(counts_ccs < 0)) {
    stop("Negative values detected in counts_ccs")
  }
  
  message("✓ Count integrity checks passed")
  
  ## ---- 2. Symmetry ---- ##
  for (i in seq_len(n_rasters)) {
    for (j in seq_len(n_rasters)) {
      
      C_ij <- counts[i, j, , ]
      C_ji <- counts[j, i, , ]
      
      if (!isTRUE(all.equal(C_ij, t(C_ji), tolerance = tol))) {
        stop(sprintf(
          "Global symmetry check failed for raster pair (%s, %s)",
          ras_names[i], ras_names[j]
        ))
      }
    }
  }
  
  for (d in seq_len(n_districts)) {
    for (i in seq_len(n_rasters)) {
      for (j in seq_len(n_rasters)) {
        
        C_ij <- counts_ecodist[i, j, , , d]
        C_ji <- counts_ecodist[j, i, , , d]
        
        if (!isTRUE(all.equal(C_ij, t(C_ji), tolerance = tol))) {
          stop(sprintf(
            "Ecodistrict symmetry check failed for district %s, raster pair (%s, %s)",
            district_vals[d], ras_names[i], ras_names[j]
          ))
        }
      }
    }
  }
  
  for (c in seq_len(n_ccs)) {
    for (i in seq_len(n_rasters)) {
      for (j in seq_len(n_rasters)) {
        
        C_ij <- counts_ccs[i, j, , , c]
        C_ji <- counts_ccs[j, i, , , c]
        
        if (!isTRUE(all.equal(C_ij, t(C_ji), tolerance = tol))) {
          stop(sprintf(
            "CCS symmetry check failed for district %s, raster pair (%s, %s)",
            ccs_vals[c], ras_names[i], ras_names[j]
          ))
        }
      }
    }
  }
  
  message("✓ Symmetry checks passed")
  
  ## ---- 3. Diagonal self-comparison ---- ##
  for (i in seq_len(n_rasters)) {
    
    Cii <- counts[i, i, , ]
    offdiag <- Cii
    diag(offdiag) <- 0

    if (any(offdiag != 0)) {
      stop(sprintf(
        "Global diagonal self-comparison failed for raster %s: off-diagonal entries are nonzero",
        ras_names[i]
      ))
    }
  }
  
  for (d in seq_len(n_districts)) {
    for (i in seq_len(n_rasters)) {
      
      Cii <- counts_ecodist[i, i, , , d]
      offdiag <- Cii
      diag(offdiag) <- 0
      
      if (any(offdiag != 0)) {
        stop(sprintf(
          "Ecodistrict diagonal self-comparison failed for district %s, raster %s",
          district_vals[d], ras_names[i]
        ))
      }
    }
  }
  
  for (c in seq_len(n_ccs)) {
    for (i in seq_len(n_rasters)) {
      
      Cii <- counts_ccs[i, i, , , c]
      offdiag <- Cii
      diag(offdiag) <- 0
      
      if (any(offdiag != 0)) {
        stop(sprintf(
          "CCS diagonal self-comparison failed for district %s, raster %s",
          ccs_vals[c], ras_names[i]
        ))
      }
    }
  }
  
  message("✓ Diagonal self-comparison checks passed")
  
  ## ---- 4. Global vs ecodistrict summed totals ---- ##
  counts_sum <- apply(counts_ecodist, c(1,2,3,4), sum)
  diff <- counts - counts_sum
  
  if (any(is.na(diff))) {
    stop("NA values detected in (counts - summed ecodistrict counts)")
  }
  
  # district totals should never exceed global totals
  if (any(diff < -tol)) {
    bad_idx <- which(diff < -tol, arr.ind = TRUE)[1, , drop = FALSE]
    stop(sprintf(
      paste0(
        "Summed ecodistrict counts exceed global counts at: ",
        "raster_i=%s, raster_j=%s, class_i=%s, class_j=%s"
      ),
      ras_names[bad_idx[1,1]],
      ras_names[bad_idx[1,2]],
      classes[bad_idx[1,3]],
      classes[bad_idx[1,4]]
    ))
  }
  
  max_diff_eco <- max(abs(diff), na.rm = TRUE)
  total_diff_eco <- sum(diff) # sum of pairwise counts located in pixels
                          # where ecodistrict is NA (not number of NA-district pixels)
  
  message("Max difference (global - summed ecodistrict): ", max_diff_eco)
  message("Total difference (global - summed ecodistrict): ", total_diff_eco)
  
  if (!is.na(na_district_pixels)) {
    message("Total pixels where ecodistrict raster is NA: ", na_district_pixels)
  }
  
  message("✓ Global vs ecodistrict sum checks passed")
  
  counts_sum <- apply(counts_ccs, c(1,2,3,4), sum)
  diff <- counts - counts_sum
  
  if (any(is.na(diff))) {
    stop("NA values detected in (counts - summed CCS counts)")
  }
  
  # CCS totals should never exceed global totals
  if (any(diff < -tol)) {
    bad_idx <- which(diff < -tol, arr.ind = TRUE)[1, , drop = FALSE]
    stop(sprintf(
      paste0(
        "Summed CCS counts exceed global counts at: ",
        "raster_i=%s, raster_j=%s, class_i=%s, class_j=%s"
      ),
      ras_names[bad_idx[1,1]],
      ras_names[bad_idx[1,2]],
      classes[bad_idx[1,3]],
      classes[bad_idx[1,4]]
    ))
  }
  
  max_diff_ccs <- max(abs(diff), na.rm = TRUE)
  total_diff_ccs <- sum(diff) # sum of pairwise counts located in pixels
  # where CCS is NA (not number of NA-CCS pixels)
  
  message("Max difference (global - summed CCS): ", max_diff_ccs)
  message("Total difference (global - summed CCS): ", total_diff_ccs)
  
  if (!is.na(na_ccs_pixels)) {
    message("Total pixels where CCS raster is NA: ", na_ccs_pixels)
  }
  
  message("✓ Global vs CCS sum checks passed")
  
  ## ---- 5. Empty ecodistricts or CCS ---- ##
  district_totals <- sapply(seq_len(n_districts), function(d) {
    sum(counts_ecodist[,,,,d])
  })
  
  empty_districts <- district_vals[district_totals == 0]
  
  if (length(empty_districts) > 0) {
    warning(
      "Empty ecodistricts detected: ",
      paste(empty_districts, collapse = ", ")
    )
  } else {
    message("✓ All ecodistricts contain data")
  }
  
  ccs_totals <- sapply(seq_len(n_ccs), function(c) {
    sum(counts_ccs[,,,,c])
  })
  
  empty_ccs <- ccs_vals[ccs_totals == 0]
  
  if (length(empty_ccs) > 0) {
    warning(
      "Empty CCS detected: ",
      paste(empty_ccs, collapse = ", ")
    )
  } else {
    message("✓ All CCS contain data")
  }
  
  ## ---- 6. Per-pair total count consistency ---- ##
  # for each raster pair, total across classes should match mirrored pair
  for (i in seq_len(n_rasters)) {
    for (j in seq_len(n_rasters)) {
      Nij <- sum(counts[i, j, , ])
      Nji <- sum(counts[j, i, , ])
      
      if (!isTRUE(all.equal(Nij, Nji, tolerance = tol))) {
        stop(sprintf(
          "Global total count mismatch for pair (%s, %s) vs (%s, %s)",
          ras_names[i], ras_names[j], ras_names[j], ras_names[i]
        ))
      }
    }
  }
  
  for (d in seq_len(n_districts)) {
    for (i in seq_len(n_rasters)) {
      for (j in seq_len(n_rasters)) {
        Nij <- sum(counts_ecodist[i, j, , , d])
        Nji <- sum(counts_ecodist[j, i, , , d])
        
        if (!isTRUE(all.equal(Nij, Nji, tolerance = tol))) {
          stop(sprintf(
            "Ecodistrict total count mismatch for district %s, pair (%s, %s)",
            district_vals[d], ras_names[i], ras_names[j]
          ))
        }
      }
    }
  }
  
  for (c in seq_len(n_ccs)) {
    for (i in seq_len(n_rasters)) {
      for (j in seq_len(n_rasters)) {
        Nij <- sum(counts_ccs[i, j, , , c])
        Nji <- sum(counts_ccs[j, i, , , c])
        
        if (!isTRUE(all.equal(Nij, Nji, tolerance = tol))) {
          stop(sprintf(
            "CCS total count mismatch for CCS %s, pair (%s, %s)",
            ccs_vals[c], ras_names[i], ras_names[j]
          ))
        }
      }
    }
  }
  
  message("✓ Total count consistency checks passed")
  
  ## ---- 7. Validate pairwise_long structure ---- ##
  required_cols <- c(
    "region_type", "region_id", "raster_i", "raster_j",
    "kappa", "accuracy",
    "jaccard_macro", "jaccard_class0", "jaccard_class1", "jaccard_class2",
    "f1_macro", "f1_class0", "f1_class1", "f1_class2",
    "mcc_multi", "N"
  )
  
  missing_cols <- setdiff(required_cols, names(pairwise_long))
  if (length(missing_cols) > 0) {
    stop("pairwise_long is missing required columns: ",
         paste(missing_cols, collapse = ", "))
  }
  
  message("✓ pairwise_long structure checks passed")
  
  ## ---- 8. Duplicated rows ---- ##
  dup_rows <- duplicated(pairwise_long[, c("region_type", "region_id", "raster_i", "raster_j")])
  if (any(dup_rows)) {
    stop("Duplicate rows detected in pairwise_long for (region_type, region_id, raster_i, raster_j)")
  }
  
  # dup_rows_ccs <- duplicated(pairwise_long[, c("region_type", "region_id", "raster_i", "raster_j")])
  # if (any(dup_rows_ccs)) {
  #   stop("Duplicate rows detected in pairwise_long for (ccs, raster_i, raster_j)")
  # }
  
  message("✓ Duplicate row checks passed")
  
  ## ---- 9. Pair completeness ---- ##
  expected_pairs <- choose(n_rasters, 2)
  
  # Global
  n_global <- sum(pairwise_long$region_type == "global")
  if (n_global != expected_pairs) {
    warning(sprintf(
      "Expected %d global raster pairs, found %d",
      expected_pairs, n_global
    ))
  } else {
    message("✓ Global pair count is complete")
  }
  
  # Per district
  district_counts <- pairwise_long %>%
    filter(region_type == "ecodistrict") %>%
    count(region_id, name = "n_pairs")
  
  incomplete <- district_counts %>%
    dplyr::filter(n_pairs != expected_pairs)
  
  if (nrow(incomplete) > 0) {
    warning(
      "Some districts do not have the expected number of raster pairs: ",
      paste0(incomplete$region_id, " (", incomplete$n_pairs, ")", collapse = ", ")
    )
  } else {
    message("✓ All districts have complete pair counts")
  }
  
  # Per CCS
  ccs_counts <- pairwise_long %>%
    filter(region_type == "ccs") %>%
    count(region_id, name = "n_pairs")
  
  incomplete <- ccs_counts %>%
    dplyr::filter(n_pairs != expected_pairs)
  
  if (nrow(incomplete) > 0) {
    warning(
      "Some CCS do not have the expected number of raster pairs: ",
      paste0(incomplete$region_id, " (", incomplete$n_pairs, ")", collapse = ", ")
    )
  } else {
    message("✓ All CCS have complete pair counts")
  }
  
  ## ---- 10. Metric bounds ---- ##
  check_range <- function(x, minv, maxv, name) {
    bad <- which(!is.na(x) & (x < minv - tol | x > maxv + tol))
    if (length(bad) > 0) {
      stop(sprintf("Metric %s out of bounds", name))
    }
  }
  
  check_range(pairwise_long$accuracy, 0, 1, "accuracy")
  check_range(pairwise_long$jaccard_macro, 0, 1, "jaccard_macro")
  check_range(pairwise_long$jaccard_class0, 0, 1, "jaccard_class0")
  check_range(pairwise_long$jaccard_class1, 0, 1, "jaccard_class1")
  check_range(pairwise_long$jaccard_class2, 0, 1, "jaccard_class2")
  check_range(pairwise_long$f1_macro, 0, 1, "f1_macro")
  check_range(pairwise_long$f1_class0, 0, 1, "f1_class0")
  check_range(pairwise_long$f1_class1, 0, 1, "f1_class1")
  check_range(pairwise_long$f1_class2, 0, 1, "f1_class2")
  check_range(pairwise_long$kappa, -1, 1, "kappa")
  check_range(pairwise_long$mcc_multi, -1, 1, "mcc_multi")
  
  if (any(!is.na(pairwise_long$N) & pairwise_long$N < 0)) {
    stop("Negative N values detected in pairwise_long")
  }
  
  message("✓ Metric bounds checks passed")
  
  ## ---- 11. Cross-check N against count arrays ---- ##
  # global rows
  global_rows <- pairwise_long %>% dplyr::filter(region_type == "global")
  
  for (r in seq_len(nrow(global_rows))) {
    i <- match(global_rows$raster_i[r], ras_names)
    j <- match(global_rows$raster_j[r], ras_names)
    N_expect <- sum(counts[i, j, , ])
    
    if (!isTRUE(all.equal(global_rows$N[r], N_expect, tolerance = tol))) {
      stop(sprintf(
        "Global N mismatch for pair (%s, %s): output=%s, expected=%s",
        global_rows$raster_i[r], global_rows$raster_j[r],
        global_rows$N[r], N_expect
      ))
    }
  }
  
  # district rows
  district_rows <- pairwise_long %>% dplyr::filter(region_type == "ecodistrict")
  
  for (r in seq_len(nrow(district_rows))) {
    i <- match(district_rows$raster_i[r], ras_names)
    j <- match(district_rows$raster_j[r], ras_names)
    d <- match(as.integer(district_rows$region_id[r]), district_vals)
    
    if (is.na(d)) {
      stop(sprintf("District %s in pairwise_long not found in district_vals",
                   district_rows$region_id[r]))
    }
    
    N_expect <- sum(counts_ecodist[i, j, , , d])
    
    if (!isTRUE(all.equal(district_rows$N[r], N_expect, tolerance = tol))) {
      stop(sprintf(
        "District N mismatch for district %s, pair (%s, %s): output=%s, expected=%s",
        district_rows$region_id[r],
        district_rows$raster_i[r],
        district_rows$raster_j[r],
        district_rows$N[r],
        N_expect
      ))
    }
  }
  
  # CCS rows
  ccs_rows <- pairwise_long %>% dplyr::filter(region_type == "ccs")
  
  for (r in seq_len(nrow(ccs_rows))) {
    i <- match(ccs_rows$raster_i[r], ras_names)
    j <- match(ccs_rows$raster_j[r], ras_names)
    c <- match(as.integer(ccs_rows$region_id[r]), ccs_vals)
    
    if (is.na(c)) {
      stop(sprintf("CCS %s in pairwise_long not found in ccs_vals",
                   ccs_rows$region_id[r]))
    }
    
    N_expect <- sum(counts_ccs[i, j, , , c])
    
    if (!isTRUE(all.equal(ccs_rows$N[r], N_expect, tolerance = tol))) {
      stop(sprintf(
        "CCS N mismatch for CCS %s, pair (%s, %s): output=%s, expected=%s",
        ccs_rows$region_id[r],
        ccs_rows$raster_i[r],
        ccs_rows$raster_j[r],
        ccs_rows$N[r],
        N_expect
      ))
    }
  }
  
  message("✓ Output N cross-checks passed")
  
  ## ---- 12. Validate exported confusion matrix counts ---- ##
  
  required_confusion_cols <- c(
    "c00","c01","c02",
    "c10","c11","c12",
    "c20","c21","c22"
  )
  
  missing_confusion_cols <- setdiff(
    required_confusion_cols,
    names(pairwise_long)
  )
  
  if (length(missing_confusion_cols) > 0) {
    stop(
      "Missing confusion matrix columns: ",
      paste(missing_confusion_cols, collapse = ", ")
    )
  }
  
  N_from_confusion <-
    pairwise_long$c00 +
    pairwise_long$c01 +
    pairwise_long$c02 +
    pairwise_long$c10 +
    pairwise_long$c11 +
    pairwise_long$c12 +
    pairwise_long$c20 +
    pairwise_long$c21 +
    pairwise_long$c22
  
  mismatch <- which(
    N_from_confusion != pairwise_long$N
  )
  
  if (length(mismatch) > 0) {
    stop(
      sprintf(
        "Confusion matrix totals do not equal N for %d row(s)",
        length(mismatch)
      )
    )
  }
  
  message("✓ Confusion matrix totals match N")
  
  ## ---- 13. Validate row/column marginal totals ---- ##
  
  if (any(pairwise_long$row0 !=
          pairwise_long$c00 +
          pairwise_long$c01 +
          pairwise_long$c02)) {
    stop("row0 does not match confusion matrix")
  }
  
  if (any(pairwise_long$row1 !=
          pairwise_long$c10 +
          pairwise_long$c11 +
          pairwise_long$c12)) {
    stop("row1 does not match confusion matrix")
  }
  
  if (any(pairwise_long$row2 !=
          pairwise_long$c20 +
          pairwise_long$c21 +
          pairwise_long$c22)) {
    stop("row2 does not match confusion matrix")
  }
  
  if (any(pairwise_long$col0 !=
          pairwise_long$c00 +
          pairwise_long$c10 +
          pairwise_long$c20)) {
    stop("col0 does not match confusion matrix")
  }
  
  if (any(pairwise_long$col1 !=
          pairwise_long$c01 +
          pairwise_long$c11 +
          pairwise_long$c21)) {
    stop("col1 does not match confusion matrix")
  }
  
  if (any(pairwise_long$col2 !=
          pairwise_long$c02 +
          pairwise_long$c12 +
          pairwise_long$c22)) {
    stop("col2 does not match confusion matrix")
  }
  
  if (any(pairwise_long$N !=
          (pairwise_long$row0 +
           pairwise_long$row1 +
           pairwise_long$row2))) {
    stop("N does not match row totals")
  }
  
  if (any(pairwise_long$N !=
          (pairwise_long$col0 +
           pairwise_long$col1 +
           pairwise_long$col2))) {
    stop("N does not match column totals")
  }
  
  message("✓ Row and column totals match confusion matrix")
  
  message("---- All multiclass validation checks completed ----")
  
  invisible(list(
    max_diff_eco = max_diff_eco,
    total_diff_eco = total_diff_eco,
    empty_districts = empty_districts,
    max_diff_ccs = max_diff_ccs,
    total_diff_ccs = total_diff_ccs,
    empty_ccs = empty_ccs
  ))
}

# call the validation function to run checks 
validate_spatial_agreement_multiclass(
  counts,
  counts_ecodist,
  counts_ccs,
  pairwise_long,
  district_vals,
  ccs_vals,
  ras_names,
  classes,
  na_district_pixels,
  na_ccs_pixels
)

## ---- Save outputs ---- ##
out_csv_local <- file.path(tmp_root, "pairwise_results.csv")
write_csv(pairwise_long, out_csv_local)

aws.s3::put_object(
  file = out_csv_local,
  object = file.path(out_dir, "pairwise_nattam_spatial_agreement_ecodist_CCS_14Oct2026.csv")
)

message("Finished pairwise spatial agreement analysis")
log_msg("Spatial agreement analysis completed successfully")

# Upload log to S3
aws.s3::put_object(
  file = log_file_local,
  object = log_file
)

message("Log written to: ", log_file)
