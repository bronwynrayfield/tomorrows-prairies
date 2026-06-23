
###################################################################################
##          ANALYSIS: SPATIAL AGREEMENT (BINARY GRASS VS NON-GRASS)              ##
###################################################################################

##----------------------- Alison Long - May 2026 -------------------------

# step 6: calculate pairwise spatial agreement metrics to assess similarities and 
#         differences among existing grassland datasets. This version handles all 
#         grassland data products regardless of extent.

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
meta_file  <- "s3://sagemaker-gst-stage.sharing/alison-long/tomorrows-prairies/processed/metadata_full_extent_binary_with_path.csv"
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

## ---- Build raster stack ---- ##
log_msg("Building raster stack")
r_stack <- rast(rfiles_local)
n_rasters <- nlyr(r_stack)
names(r_stack) <- names(rfiles_local)
ras_names <- names(r_stack)

## ---- Initialize count matrices ---- ##
n11 <- matrix(0, n_rasters, n_rasters)
n10 <- matrix(0, n_rasters, n_rasters)
n01 <- matrix(0, n_rasters, n_rasters)
n00 <- matrix(0, n_rasters, n_rasters)

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

for (i in 1:total_blocks) {
  
  message(sprintf("Block %d / %d", i, total_blocks))
  
  vals <- readValues(
    r_stack,
    row = bs$row[i],
    nrows = bs$nrows[i],
    mat = TRUE
  )
  
  # Use integers to avoid double precision overhead and cut RAM
  storage.mode(vals) <- "integer"
  
  # Skip empty block
  if (all(is.na(vals))) next
  
  # Keep vals (not M yet)
  valid_rows <- rowSums(!is.na(vals)) > 0
  vals <- vals[valid_rows, , drop = FALSE]
  
  if (nrow(vals) == 0) next
  
  # Convert to binary (0/1)
  M <- vals > 0
  #storage.mode(M) <- "integer"
  
  # Build NA mask per raster so NA pixels can be excluded from comparison
  # Pixels beyond raster extent should not be included in analysis
  # Converting NAs to 0s would inflate agreement and penalize rasters with smaller extent
  NA_mask <- is.na(vals)
  
  ## ---- Matrix algebra ---- ##
  # Valid co-observation mask
  valid_mat <- !NA_mask
  
  # Convert logical to integer
  M_int <- M * 1
  valid_int <- valid_mat * 1
  
  # Zero out invalid cells
  M_int[valid_int == 0] <- 0
  
  # valid counts for each pair
  N_block_mat <- t(valid_int) %*% valid_int
  
  # n11: both 1 AND both observed
  n11_block <- t(M_int) %*% M_int
  
  # n10 and n01 only where both rasters are valid
  present_i <- t(M_int) %*% valid_int
  present_j <- t(valid_int) %*% M_int
  
  # per-raster presence
  #col_sum <- colSums(M_int)
  
  # n10 and n01 (only where valid)
  #n10_block <- outer(col_sum, col_sum, "-") - n11_block
  #n01_block <- t(n10_block)
  
  n10_block <- present_i - n11_block
  n01_block <- present_j - n11_block
  
  # n00 derived correctly - does not include NAs
  n00_block <- N_block_mat - (n11_block + n10_block + n01_block)
  
  # Validate correct handling of NAs
  if (anyNA(n11_block) || anyNA(n10_block) ||
      anyNA(n01_block) || anyNA(n00_block)) {
    stop("NA detected in block matrices")
  }
  
  ## ---- Accumulate ---- ##
  n11 <- n11 + n11_block
  n10 <- n10 + n10_block
  n01 <- n01 + n01_block
  n00 <- n00 + n00_block
  
  # Force garbage collection per block to prevent memory accumulation
  rm(vals, M, n11_block, n10_block, n01_block, n00_block)
  gc()
}

# ✅ CLOSE FILE CONNECTIONS
readStop(r_stack)

message("Finished block processing")
log_msg("Finished block processing")

## ---- Compute metrics ---- ##
message("Computing metrics...")
log_msg("Computing metrics")

results <- list()

for (i in 1:(n_rasters - 1)) {
  for (j in (i + 1):n_rasters) {
    
    tp <- n11[i,j]
    fp <- n10[i,j]
    fn <- n01[i,j]
    tn <- n00[i,j]
    
    N <- tp + fp + fn + tn
    
    if (is.na(N) || N == 0) {
      d_j <- 1
      d_s <- 1
      kappa <- NA
      f1 <- NA
      mcc <- NA
      
    } else {
      
      # ---- Jaccard ----
      denom_j <- tp + fp + fn
      d_j <- ifelse(denom_j == 0, 0, 1 - tp / denom_j)
      
      # ---- Sørensen ----
      denom_s <- 2*tp + fp + fn
      d_s <- ifelse(denom_s == 0, 0, 1 - (2*tp / denom_s))
      
      # ---- Kappa ----
      po <- (tp + tn) / N
      p1_yes <- (tp + fp) / N
      p1_no  <- (fn + tn) / N
      p2_yes <- (tp + fn) / N
      p2_no  <- (fp + tn) / N
      pe <- p1_yes*p2_yes + p1_no*p2_no
      
      kappa <- ifelse(pe == 1, NA, (po - pe) / (1 - pe))
      
      # ---- F1 score ----
      precision <- ifelse((tp + fp) == 0, NA, tp / (tp + fp))
      recall    <- ifelse((tp + fn) == 0, NA, tp / (tp + fn))
      
      f1 <- ifelse(is.na(precision) | is.na(recall) | (precision + recall == 0),
                   NA,
                   2 * precision * recall / (precision + recall))
      
      # ---- MCC ----
      denom_mcc <- sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
      
      mcc <- ifelse(denom_mcc == 0,
                    NA,
                    (tp * tn - fp * fn) / denom_mcc)
    }
    
    results[[length(results) + 1]] <- tibble(
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

pairwise_long <- bind_rows(results)

log_msg("Done computing metrics")

## ---- Save outputs ---- ##
out_csv_local <- file.path(tmp_root, "pairwise_results.csv")
write_csv(pairwise_long, out_csv_local)

aws.s3::put_object(
  file = out_csv_local,
  object = file.path(out_dir, "pairwise_binary_spatial_agreement_01Jun2026.csv")
)

message("Finished pairwise spatial agreement analysis")
log_msg("Spatial agreement analysis completed successfully")

# Upload log to S3
aws.s3::put_object(
  file = log_file_local,
  object = log_file
)

message("Log written to: ", log_file)
