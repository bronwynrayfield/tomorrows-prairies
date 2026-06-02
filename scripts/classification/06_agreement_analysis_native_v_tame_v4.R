
###################################################################################
##          ANALYSIS: SPATIAL AGREEMENT (NATIVE VS TAME VS NON-GRASS)            ##
###################################################################################

##----------------------- Alison Long - May 2026 -------------------------

# step 6: calculate pairwise spatial agreement metrics to assess similarities and 
#         differences among existing grassland datasets. This version handles all 
#         grassland data products regardless of extent.

# For instance type, best option is r6i.4xlarge (budget option r6i.2xlarge) 
# ran using ml.m5d.8xlarge
# Prefer instances with NVMe local storage

## ---- Setup ---- ##
# Temp dirs
tmp_root <- file.path(Sys.getenv("HOME"), "r_tmp")
#dir.create(tmp_root, recursive = TRUE, showWarnings = FALSE)

# Load packages
library(terra)
library(future.apply)  # running sequentially to avoid disk contention
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
#n_cores <- min(6, parallel::detectCores(logical = FALSE) - 2)
future::plan(sequential)  # Avoids disk contention

## ---- Settings ---- ##
meta_file  <- "s3://sagemaker-gst-stage.sharing/alison-long/tomorrows-prairies/processed/metadata_full_extent_nattam_with_path.csv"
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
#log_msg("Parallel workers detected: ", n_cores)

## ---- Load metadata ---- ##
log_msg("Reading metadata from: ", meta_file)

meta <- aws.s3::s3read_using(read_csv, object = meta_file, show_col_types = FALSE) %>%
  rename(file = file_s3)

rfiles_s3 <- setNames(meta$file, meta$raster)

log_msg("Metadata rows loaded: ", nrow(meta))

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

rfiles_local <- sapply(names(rfiles_s3), function(nm) {
  
  s3_path <- rfiles_s3[[nm]]
  
  s3_clean <- sub("^s3://", "", s3_path)
  parts <- strsplit(s3_clean, "/", fixed = TRUE)[[1]]
  
  bucket <- parts[1]
  key <- paste(parts[-1], collapse = "/")
  
  dest <- file.path(local_dir, paste0(nm, ".tif"))
  
  aws.s3::save_object(object = key, bucket = bucket, file = dest)
  
  dest
}, USE.NAMES = TRUE)

log_msg("All rasters downloaded locally")

# Check file size to ensure correct download
bad_files <- names(rfiles_local)[file.size(rfiles_local) < 10000]
if (length(bad_files) > 0) {
  stop(paste("These rasters failed to download correctly:",
             paste(bad_files, collapse = ", ")))
}

file.info(rfiles_local)
rfiles_s3[c("aci", "ngi", "cgr")]

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

## ---- Build raster stack ---- ##
log_msg("Building raster stack")
r_stack <- rast(rfiles_local)
n_rasters <- nlyr(r_stack)
names(r_stack) <- names(rfiles_local)
ras_names <- names(r_stack)

## ---- Initialize count matrices ---- ##
classes <- c(0,1,2)
#stopifnot(all(values(r_stack) %in% classes | is.na(values(r_stack))))

n_classes <- length(classes)

counts <- array(
  0,
  dim = c(n_rasters, n_rasters, n_classes, n_classes),
  dimnames = list(ras_names, ras_names, classes, classes)
)

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

# OPEN FILE CONNECTIONS
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
  
  # Keep all 3 categorical values
  M <- vals 
  #storage.mode(M) <- "integer"
  
  ## ---- Class-wise counting ---- ##
  # For each pair, only count pixels where both rasters are non-NA
  for (ri in 1:(n_rasters - 1)) {
    for (rj in (ri + 1):n_rasters) {
      
      valid <- !is.na(M[,ri]) & !is.na(M[,rj])
      if (!any(valid)) next
      
      Mi <- M[valid, ri]
      Mj <- M[valid, rj]
      
      # Build contingency table
      tab <- table(factor(Mi, levels = classes),
                   factor(Mj, levels = classes))
      
      # Add to counts
      counts[ri, rj, , ] <- counts[ri, rj, , ] + tab
      counts[rj, ri, , ] <- counts[rj, ri, , ] + t(tab)
    }
  }
  
  # Force garbage collection to prevent memory accumulation
  rm(vals, M)
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
    jaccard_macro <- mean(jacc, na.rm = TRUE)
    
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
    f1_macro <- mean(f1_class, na.rm = TRUE)
    
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
    results[[length(results)+1]] <- tibble(
      raster_i = ras_names[i],
      raster_j = ras_names[j],
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

pairwise_long <- bind_rows(results)

log_msg("Done computing metrics")

## ---- Save outputs ---- ##
out_csv_local <- file.path(tmp_root, "pairwise_results.csv")
write_csv(pairwise_long, out_csv_local)

aws.s3::put_object(
  file = out_csv_local,
  object = file.path(out_dir, "pairwise_nattam_spatial_agreement_01Jun2026.csv")
)

message("Finished pairwise spatial agreement analysis")
log_msg("Spatial agreement analysis completed successfully")

# Upload log to S3
aws.s3::put_object(
  file = log_file_local,
  object = log_file
)

message("Log written to: ", log_file)
