
###################################################################################
##          ANALYSIS: SPATIAL AGREEMENT (BINARY GRASS VS NON-GRASS)              ##
###################################################################################

##----------------------- Alison Long - June 2026 -------------------------

# step 6: calculate pairwise spatial agreement metrics to assess similarities and 
#         differences among existing grassland datasets. This version handles all 
#         grassland data products regardless of extent and calculates metrics at
#         the ecodistrict level.

#### SAME AS V4 BUT ADD F1 and MCC AS METRICS ####

# For instance type, best option is r6i.4xlarge (budget option r6i.2xlarge) 
# trying ml.m5d.8xlarge
# Prefer instances with NVMe local storage

## ---- Setup ---- ##
# Temp dirs
tmp_root <- file.path(Sys.getenv("HOME"), "r_tmp")
#dir.create(tmp_root, recursive = TRUE, showWarnings = FALSE)

# Load packages
library(terra)
library(future.apply)
library(readr)
library(dplyr)
library(aws.s3)

# Set terra options
terraOptions(
  tempdir = file.path(tmp_root, "terra"),
  memfrac = 0.6, # could go down to 0.4
  progress = 1
)

## ---- Parallel (light use only) ---- ##
n_cores <- min(6, parallel::detectCores(logical = FALSE) - 2)
future::plan(sequential)  # Avoids disk contention

## ---- Settings ---- ##
meta_file    <- "s3://sagemaker-gst-stage.sharing/alison-long/tomorrows-prairies/processed/metadata_full_extent_binary_with_path.csv"
ecodist_file <- "s3://sagemaker-gst-stage.sharing/alison-long/tomorrows-prairies/processed/ecodistricts_clean_rast.tif"
out_dir      <- "s3://sagemaker-gst-stage.sharing/alison-long/tomorrows-prairies/results/spatial_agreement_analysis"
log_dir      <- "s3://sagemaker-gst-stage.sharing/alison-long/tomorrows-prairies/logs"
local_dir    <- file.path(tmp_root, "data")
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

log_msg <- function(..., level = "INFO") {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  msg <- paste0("[", level, "] ", ts, " : ", paste(..., collapse = ""))
  
  # Console
  message(msg)
  
  # Append to log file
  cat(msg, "\n", file = log_file_local, append = TRUE)
}

log_msg("Starting spatial agreement analysis")
log_msg("Parallel workers detected: ", n_cores)

## ---- Load metadata ---- ##
log_msg("Reading metadata from: ", meta_file)

meta <- aws.s3::s3read_using(read_csv, object = meta_file, show_col_types = FALSE) %>%
  rename(file = file_s3)

rfiles_s3 <- setNames(meta$file, meta$raster)

log_msg("Metadata rows loaded: ", nrow(meta))

## ---- Download rasters locally ---- ##
message("Downloading rasters locally...")
log_msg("Downloading rasters locally")

rfiles_local <- sapply(names(rfiles_s3), function(nm) {
  dest <- file.path(local_dir, paste0(nm, ".tif"))
  if (!file.exists(dest)) {
    aws.s3::save_object(rfiles_s3[[nm]], file = dest)
  }
  dest
}, USE.NAMES = TRUE)

log_msg("All rasters downloaded locally")

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
ecodist_local <- file.path(local_dir, "ecodistricts.tif")
aws.s3::save_object(ecodist_file, file = ecodist_local)

ecodist_rast <- rast(ecodist_local)

# Ensure alignment
if (!compareGeom(ecodist_rast, template, stopOnError = FALSE)) {
  stop("Ecodistrict raster is not aligned with input rasters")
}

log_msg("Loaded ecodistrict raster")

## ---- Build raster stack ---- ##
log_msg("Building raster stack")
r_stack <- rast(rfiles_local)
n_rasters <- nlyr(r_stack)
names(r_stack) <- names(rfiles_local)
ras_names <- names(r_stack)

## ---- Initialize district-level count arrays ---- ##

# Get unique district IDs (could also delay until first block)
district_vals <- as.integer(sort(unique(na.omit(values(ecodist_rast, mat = FALSE)))))

# initialize empty set
#district_vals <- integer(0)

n_districts <- length(district_vals)

# Map district IDs → index
district_index <- setNames(seq_along(district_vals), district_vals)

# Safe indexing
#district_vals <- as.integer(district_vals)

# 3D arrays: [raster_i, raster_j, district]
n11 <- array(0, dim = c(n_rasters, n_rasters, n_districts))
n10 <- array(0, dim = c(n_rasters, n_rasters, n_districts))
n01 <- array(0, dim = c(n_rasters, n_rasters, n_districts))
n00 <- array(0, dim = c(n_rasters, n_rasters, n_districts))

# Global counts 
n11_global <- matrix(0, n_rasters, n_rasters)
n10_global <- matrix(0, n_rasters, n_rasters)
n01_global <- matrix(0, n_rasters, n_rasters)
n00_global <- matrix(0, n_rasters, n_rasters)

# Count pixels where ecodistrict = NA
na_district_pixels <- 0

## ---- Block processing ---- ##
message("Starting blockwise processing...")
log_msg("Starting blockwise processing")

# Block size determined by terraOptions(memfrac=...), available RAM, and number of layers
# Can force smaller chunks for lower RAM per iteration (increase # of blocks)
bs <- blocks(r_stack) #, n = 10000) # Try 8-20k
str(bs)
total_blocks <- bs$n
summary(bs$nrows)

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

# ✅ OPEN FILE CONNECTIONS
readStart(r_stack)
readStart(ecodist_rast)

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
  
  # Count NA district pixels
  na_district_pixels <- na_district_pixels + sum(is.na(districts))
  
  # Safe indexing
  districts <- as.integer(districts)
  
  #storage.mode(vals) <- "integer"
  
  # Drop fully empty rows
  valid_rows <- rowSums(!is.na(vals)) > 0 & !is.na(districts)
  vals <- vals[valid_rows, , drop = FALSE]
  districts <- districts[valid_rows]
  
  if (nrow(vals) == 0) next
  
  # Binary presence
  M <- vals > 0
  NA_mask <- is.na(vals)
  
  valid_mat <- !NA_mask
  M_int <- M * 1
  valid_int <- valid_mat * 1
  M_int[valid_int == 0] <- 0
  
  storage.mode(valid_int) <- "double"
  storage.mode(M_int) <- "double"
  
  # ---- GLOBAL (all pixels in block) ----
  N_block_all <- t(valid_int) %*% valid_int
  n11_block_all <- t(M_int) %*% M_int
  
  present_i_all <- t(M_int) %*% valid_int
  present_j_all <- t(valid_int) %*% M_int
  
  n10_block_all <- present_i_all - n11_block_all
  n01_block_all <- present_j_all - n11_block_all
  n00_block_all <- N_block_all - (n11_block_all + n10_block_all + n01_block_all)
  
  # accumulate
  n11_global <- n11_global + n11_block_all
  n10_global <- n10_global + n10_block_all
  n01_global <- n01_global + n01_block_all
  n00_global <- n00_global + n00_block_all
  
  ## ---- Split by district ---- ##
  split_idx <- split(seq_len(nrow(vals)), districts)
  
  for (d in names(split_idx)) {
    
    rows <- split_idx[[d]]
    
    d_id <- district_index[as.character(d)]
    if (is.na(d_id)) next
    d_id <- as.integer(d_id)
    
    M_d <- M_int[rows, , drop = FALSE]
    valid_d <- valid_int[rows, , drop = FALSE]
    
    # ---- Matrix algebra ----
    N_block <- t(valid_d) %*% valid_d
    n11_block <- t(M_d) %*% M_d
    
    present_i <- t(M_d) %*% valid_d
    present_j <- t(valid_d) %*% M_d
    
    n10_block <- present_i - n11_block
    n01_block <- present_j - n11_block
    n00_block <- N_block - (n11_block + n10_block + n01_block)
    
    # ---- Accumulate per district ----
    n11[,,d_id] <- n11[,,d_id] + n11_block
    n10[,,d_id] <- n10[,,d_id] + n10_block
    n01[,,d_id] <- n01[,,d_id] + n01_block
    n00[,,d_id] <- n00[,,d_id] + n00_block
  }
  
  rm(vals, districts, M, M_int, valid_int, split_idx)
  if (i %% 10 == 0) gc()
}

readStop(r_stack)
readStop(ecodist_rast)

message("Finished block processing")
log_msg("Finished block processing")

message("Total pixels where district = NA: ", na_district_pixels)
log_msg("Total pixels where district = NA: ", na_district_pixels)

## ---- Compute metrics per district ---- ##
message("Computing metrics...")
log_msg("Computing metrics")

results <- list()

for (d in seq_len(n_districts)) {
  
  d_name <- district_vals[d]
  
  for (i in 1:(n_rasters - 1)) {
    for (j in (i + 1):n_rasters) {
      
      tp <- n11[i,j,d]
      fp <- n10[i,j,d]
      fn <- n01[i,j,d]
      tn <- n00[i,j,d]
      
      N <- tp + fp + fn + tn
      
      if (is.na(N) || N == 0) {
        d_j <- 1; d_s <- 1
        kappa <- NA; f1 <- NA; mcc <- NA
      } else {
        denom_j <- tp + fp + fn
        d_j <- ifelse(denom_j == 0, 0, 1 - tp / denom_j)
        
        denom_s <- 2*tp + fp + fn
        d_s <- ifelse(denom_s == 0, 0, 1 - (2*tp / denom_s))
        
        po <- (tp + tn) / N
        p1_yes <- (tp + fp) / N
        p1_no  <- (fn + tn) / N
        p2_yes <- (tp + fn) / N
        p2_no  <- (fp + tn) / N
        pe <- p1_yes*p2_yes + p1_no*p2_no
        
        kappa <- ifelse(pe == 1, NA, (po - pe)/(1 - pe))
        
        precision <- ifelse((tp+fp)==0, NA, tp/(tp+fp))
        recall    <- ifelse((tp+fn)==0, NA, tp/(tp+fn))
        
        f1 <- ifelse(is.na(precision) | is.na(recall) | (precision+recall==0),
                     NA,
                     2*precision*recall/(precision+recall))
        
        denom_mcc <- sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
        mcc <- ifelse(denom_mcc==0, NA,
                      (tp*tn - fp*fn)/denom_mcc)
      }
      
      results[[length(results) + 1]] <- tibble(
        district = d_name,
        raster_i = ras_names[i],
        raster_j = ras_names[j],
        d_jaccard = d_j,
        d_sorensen = d_s,
        kappa = kappa,
        f1 = f1,
        mcc = mcc,
        N = N
      )
    }
  }
}

global_results <- list()

for (i in 1:(n_rasters - 1)) {
  for (j in (i + 1):n_rasters) {
    
    tp <- n11_global[i,j]
    fp <- n10_global[i,j]
    fn <- n01_global[i,j]
    tn <- n00_global[i,j]
    
    N <- tp + fp + fn + tn
    
    if (is.na(N) || N == 0) {
      d_j <- 1; d_s <- 1
      kappa <- NA; f1 <- NA; mcc <- NA
    } else {
      denom_j <- tp + fp + fn
      d_j <- ifelse(denom_j == 0, 0, 1 - tp / denom_j)
      
      denom_s <- 2*tp + fp + fn
      d_s <- ifelse(denom_s == 0, 0, 1 - (2*tp / denom_s))
      
      po <- (tp + tn) / N
      p1_yes <- (tp + fp) / N
      p1_no  <- (fn + tn) / N
      p2_yes <- (tp + fn) / N
      p2_no  <- (fp + tn) / N
      pe <- p1_yes*p2_yes + p1_no*p2_no
      
      kappa <- ifelse(pe == 1, NA, (po - pe)/(1 - pe))
      
      precision <- ifelse((tp+fp)==0, NA, tp/(tp+fp))
      recall    <- ifelse((tp+fn)==0, NA, tp/(tp+fn))
      
      f1 <- ifelse(is.na(precision) | is.na(recall) | (precision+recall==0),
                   NA,
                   2*precision*recall/(precision+recall))
      
      denom_mcc <- sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
      mcc <- ifelse(denom_mcc==0, NA,
                    (tp*tn - fp*fn)/denom_mcc)
    }
    
    global_results[[length(global_results) + 1]] <- tibble(
      district = "ALL",
      raster_i = ras_names[i],
      raster_j = ras_names[j],
      d_jaccard = d_j,
      d_sorensen = d_s,
      kappa = kappa,
      f1 = f1,
      mcc = mcc,
      N = N
    )
  }
}

pairwise_long <- bind_rows(
  bind_rows(global_results),
  bind_rows(results)
)

log_msg("Done computing metrics")

##------ Final checks ----##
validate_spatial_agreement <- function(
    n11_global, n10_global, n01_global, n00_global,
    n11, n10, n01, n00,
    pairwise_long,
    district_vals,
    tol = 1e-6
) {
  
  message("---- Running validation checks ----")
  
  ## ---- 1. Symmetry ---- ##
  if (!all.equal(n11_global, t(n11_global))) {
    stop("n11_global is not symmetric")
  }
  if (!all.equal(n00_global, t(n00_global))) {
    stop("n00_global is not symmetric")
  }
  if (!all.equal(n10_global, t(n01_global))) {
    stop("n10_global != transpose(n01_global)")
  }
  message("✓ Symmetry checks passed")
  
  ## ---- 2. Count reconstruction ---- ##
  N_global <- n11_global + n10_global + n01_global + n00_global
  
  if (any(N_global < 0)) {
    stop("Negative counts detected in global totals")
  }
  if (any(is.na(N_global))) {
    stop("NA detected in global totals")
  }
  message("✓ Count reconstruction checks passed")
  
  ## ---- 3. Diagonal (self-comparison) ---- ##
  diag_ok <- TRUE
  for (i in seq_len(nrow(n11_global))) {
    if (n10_global[i,i] != 0 || n01_global[i,i] != 0) {
      diag_ok <- FALSE
    }
  }
  
  if (!diag_ok) {
    stop("Diagonal check failed: self-comparison not consistent")
  }
  message("✓ Diagonal checks passed")
  
  ## ---- 4. Global vs district totals ---- ##
  if (length(dim(n11)) == 3) {
    
    n11_sum <- apply(n11, c(1,2), sum)
    
    diff <- n11_global - n11_sum
    max_diff <- max(abs(diff), na.rm = TRUE)
    
    message("Max difference (global - sum of districts): ", max_diff)
    
    if (max_diff > 0) {
      message("  (This should reflect pixels with NA district)")
    }
  }
  
  ## ---- 5. Empty districts ---- ##
  if (length(dim(n11)) == 3) {
    
    district_totals <- sapply(seq_along(district_vals), function(d) {
      sum(n11[,,d] + n10[,,d] + n01[,,d] + n00[,,d])
    })
    
    empty <- which(district_totals == 0)
    
    if (length(empty) > 0) {
      warning("Empty districts detected: ",
              paste(district_vals[empty], collapse = ", "))
    } else {
      message("✓ All districts contain data")
    }
  }
  
  ## ---- 6. Metric bounds ---- ##
  check_range <- function(x, minv, maxv, name) {
    bad <- which(!is.na(x) & (x < minv | x > maxv))
    if (length(bad) > 0) {
      stop(sprintf("Metric %s out of bounds", name))
    }
  }
  
  check_range(pairwise_long$d_jaccard, 0, 1, "Jaccard")
  check_range(pairwise_long$d_sorensen, 0, 1, "Sorensen")
  check_range(pairwise_long$f1, 0, 1, "F1")
  check_range(pairwise_long$mcc, -1, 1, "MCC")
  
  message("✓ Metric bounds checks passed")
  
  message("---- All validation checks completed ----")
}

log_msg("Running validation checks...")

validate_spatial_agreement(
  n11_global, n10_global, n01_global, n00_global,
  n11, n10, n01, n00,
  pairwise_long,
  district_vals
)

## ---- Save outputs ---- ##
out_csv_local <- file.path(tmp_root, "pairwise_results.csv")
write_csv(pairwise_long, out_csv_local)

aws.s3::put_object(
  file = out_csv_local,
  object = file.path(out_dir, "pairwise_binary_spatial_agreement_ecodist_24Jun2026.csv")
)

message("Finished pairwise spatial agreement analysis")
log_msg("Spatial agreement analysis completed successfully")

# Upload log to S3
aws.s3::put_object(
  file = log_file_local,
  object = file.path(log_dir, basename(log_file_local))
)

message("Log written to: ", log_file)
