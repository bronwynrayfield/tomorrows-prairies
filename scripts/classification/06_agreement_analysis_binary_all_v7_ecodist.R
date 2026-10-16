
###################################################################################
##          ANALYSIS: SPATIAL AGREEMENT (BINARY GRASS VS NON-GRASS)              ##
###################################################################################

##----------------------- Alison Long - June 2026 -------------------------

# step 6: calculate pairwise spatial agreement metrics to assess similarities and 
#         differences among existing grassland datasets. This version handles all 
#         grassland data products regardless of extent and calculates metrics at
#         the ecodistrict level.

# CROP TO A STUDY OR PILOT AREA - update script depending on desired extent!

#### SAME AS V4 BUT ADD F1 and MCC AS METRICS ####

# For instance type, best option is r6i.4xlarge (budget option r6i.2xlarge) 
# trying ml.m5d.8xlarge
# Prefer instances with NVMe local storage

## ---- Setup ---- ##
# Temp dirs
tmp_root <- file.path(Sys.getenv("HOME"), "r_tmp")
#dir.create(tmp_root, recursive = TRUE, showWarnings = FALSE)

remove.packages("terra")
install.packages("terra")

# Load packages
library(terra)
#library(future.apply)
library(readr)
library(dplyr)
#library(aws.s3)
library(paws)
library(paws.storage)
library(ggplot2)

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
meta_file    <- "s3://sagemaker-gst-stage.sharing/alison-long/tomorrows-prairies/processed/metadata_full_extent_binary_with_path_Oct2026.csv"
#ecodist_file <- "s3://sagemaker-gst-stage.sharing/alison-long/tomorrows-prairies/processed/ecodistricts_clean_rast.tif"
ecodist_file <- "s3://sagemaker-gst-stage.sharing/alison-long/tomorrows-prairies/processed/boundaries/national-ecological-framework-for-canada/Ecodistricts/ecodistricts_clean_rast_Oct2026.tif"
ccs_file     <- "s3://sagemaker-gst-stage.sharing/alison-long/tomorrows-prairies/processed/boundaries/census-consolidated-subdivisions/ccs_aligned_class_sa_rast.tif"
lookup_file  <- "s3://sagemaker-gst-stage.sharing/alison-long/tomorrows-prairies/processed/boundaries/census-consolidated-subdivisions/ccs_aligned_class_sa_rast_lookup.csv"
pilot_file   <- "s3://sagemaker-gst-stage.sharing/alison-long/tomorrows-prairies/processed/pilot-area/pilot_area_aligned_cnvc_Jan2026.tif"
#sa_file      <- "s3://sagemaker-gst-stage.sharing/alison-long/tomorrows-prairies/processed/study-area/modelling_extent_aligned_cnvc_Jan2026.tif"
sa_file      <- "s3://sagemaker-gst-stage.sharing/alison-long/tomorrows-prairies/processed/study-area/temperate_scope_aligned_cnvc_Apr2026.tif"
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

log_con <- file(log_file_local, open = "wt")

sink(log_con, split = TRUE)

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
  
  cat(msg, "\n",
      file = log_file_local,
      append = TRUE)
}

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
#log_msg("aws.s3 version: ", packageVersion("aws.s3"))
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

log_msg("Starting spatial agreement analysis")
#log_msg("Parallel workers detected: ", n_cores)

## ---- Reusable helper function for reading rasters ---- ##
read_s3_raster <- function(s3_path) {
  
  s3_uri <- sub("^s3://", "", s3_path)
  
  bucket <- sub("/.*", "", s3_uri)
  key <- sub("^[^/]*/", "", s3_uri)
  
  tmp_file <- tempfile(fileext = tools::file_ext(key))
  
  s3 <- paws.storage::s3()
  
  obj <- s3$get_object(
    Bucket = bucket,
    Key = key
  )
  
  writeBin(obj$Body, tmp_file)
  
  terra::rast(tmp_file)
}

## ---- Load metadata ---- ##
log_msg("Reading metadata from: ", meta_file)

s3_uri <- sub("^s3://", "", meta_file)

bucket <- sub("/.*", "", s3_uri)
key <- sub("^[^/]*/", "", s3_uri)

s3 <- paws.storage::s3()

obj <- s3$get_object(Bucket = bucket, Key = key)

meta <- readr::read_csv(rawConnection(obj$Body), show_col_types = FALSE) %>%
  rename(file = file_s3)

#meta <- aws.s3::s3read_using(read_csv, object = meta_file, show_col_types = FALSE) %>%
#  rename(file = file_s3)

rfiles_s3 <- setNames(meta$file, meta$raster)

stopifnot(!anyDuplicated(meta$raster))

log_msg("Metadata rows loaded: ", nrow(meta))

## ---- Load pilot area ---- ##
log_msg("Loading pilot area")

# pilot_local <- file.path(local_dir, "pilot_area.tif")
# 
# aws.s3::save_object(
#   object = pilot_file,
#   file = pilot_local
# )

# Read raster
pilot <- read_s3_raster(pilot_file)

log_msg("Pilot area loaded")

## ---- Load study area ---- ##
log_msg("Loading study area")

# sa_local <- file.path(local_dir, "study_area.tif")
# 
# aws.s3::save_object(
#   object = sa_file,
#   file = sa_local
# )
# 
# sa <- rast(sa_local)

# Read raster
sa <- read_s3_raster(sa_file)

log_msg("Study area loaded")

## ---- Download rasters locally ---- ##
message("Downloading rasters locally...")
log_msg("Downloading rasters locally")

# download_files <- function(files_named) {
#   
#   if (length(files_named) == 0)
#     return(setNames(character(0), character(0)))
#   
#   sapply(names(files_named), function(nm) {
#     
#     raw_file <- file.path(local_dir, paste0(nm, "_raw.tif"))
#     crop_file <- file.path(local_dir, paste0(nm, ".tif"))
#     
#     if (!file.exists(crop_file) || FORCE_REBUILD) {
#       
#       # download original
#       aws.s3::save_object(files_named[[nm]], file = raw_file)
#       
#       # crop + mask to pilot
#       r <- rast(raw_file)
#       
#       r <- crop(r, pilot)
#       r <- mask(r, pilot)
#       
#       writeRaster(
#         r,
#         crop_file,
#         overwrite = TRUE,
#         wopt = list(
#           gdal = c(
#             "COMPRESS=ZSTD",
#             "TILED=YES",
#             "BIGTIFF=YES"
#           )
#         )
#       )
#       
#       rm(r)
#       gc()
#     }
#     
#     crop_file
#     
#   }, USE.NAMES = TRUE)
#   
# }

download_files <- function(files_named) {
  
  if (length(files_named) == 0)
    return(setNames(character(0), character(0)))
  
  sapply(names(files_named), function(nm) {
    
    log_msg("Processing raster: ", nm)
    log_msg("S3 path: ", files_named[[nm]])
    
    crop_file <- file.path(local_dir, paste0(nm, ".tif"))
    
    if (!file.exists(crop_file) || FORCE_REBUILD) {
      
      # download and read raster from S3
      r <- read_s3_raster(files_named[[nm]])
      
      # crop + mask to pilot
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
log_msg("Loading pilot-clipped ecodistrict raster")

ecodist_r <- read_s3_raster(ecodist_file)

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

ccs_r <- read_s3_raster(ccs_file)
ccs_r
#plot(ccs_r)
str(levels(ccs_r))

if (!compareGeom(template, pilot, stopOnError = FALSE)) {
  stop("Pilot raster is not aligned with template raster")
}

ccs_r <- crop(ccs_r, pilot)
ccs_r <- mask(ccs_r, pilot)

ccs_rast <- as.int(ccs_r)
# ccs_rast <- (ccs_r)
# ccs_rast

# Ensure alignment
if (!compareGeom(ccs_rast, template, stopOnError = FALSE)) {
  stop("Census consolidated subdivision raster is not aligned with input rasters")
}

# create lookup table for CCSUIDs
#ccs_lookup <- levels(ccs_rast)[[1]]

# load lookup file
s3_uri <- sub("^s3://", "", lookup_file)

bucket <- sub("/.*", "", s3_uri)
key <- sub("^[^/]*/", "", s3_uri)

obj <- s3$get_object(Bucket = bucket, Key = key)

ccs_lookup <- readr::read_csv(rawConnection(obj$Body), show_col_types = FALSE)

stopifnot(all(c("ID", "CCSUID") %in% names(ccs_lookup)))

log_msg("Loaded census consolidated subdivision raster")

## ---- Build raster stack ---- ##
log_msg("Building raster stack")

r_stack <- rast(rfiles_local)
n_rasters <- nlyr(r_stack)
names(r_stack) <- names(rfiles_local)
ras_names <- names(r_stack)

log_msg("Raster stack built")

## ---- Create faceted raster figure ---- ##
log_msg("Creating faceted raster figure")

# aggregate for plotting (adjust factor as needed)
plot_stack <- aggregate(r_stack, fact = 20, fun = "modal", na.rm = TRUE)

# convert to data frame
plot_df <- as.data.frame(plot_stack, xy = TRUE, na.rm = FALSE)
plot_df <- tidyr::pivot_longer(plot_df, cols = all_of(ras_names), 
                               names_to = "raster", values_to = "value")

# convert binary values to factor for cleaner legend
plot_df$value <- factor(plot_df$value, levels = c(0, 1), labels = c("Non-grass", "Grassland"))

# create plot
facet_png <- file.path(tmp_root, "input_rasters_faceted.png")

p <- ggplot(plot_df) +
  geom_raster(aes(x = x, y = y, fill = value)) +
  facet_wrap(~ raster, ncol = 4) +
  scale_fill_manual(values = c("Non-grass" = "grey90", "Grassland" = "darkgreen"), na.value = "white") +
  coord_equal() +
  theme_void() +
  theme(strip.text = element_text(size = 8, face = "bold"), legend.position = "bottom") + 
  labs(fill = NULL, title = "Input binary raster layers")

ggsave(facet_png, p, width = 16, height = ceiling(n_rasters / 4) * 3, dpi = 300, bg = "white")

log_msg("Faceted raster figure saved")

## ---- Initialize district-level count arrays ---- ##
log_msg("Initializing ecodistrict count arrays")

# Get unique district IDs (could also delay until first block)
#district_vals <- as.integer(sort(unique(na.omit(values(ecodist_rast, mat = FALSE)))))
freq_table_eco <- freq(ecodist_rast)
district_vals <- as.integer(freq_table_eco$value)

if (anyDuplicated(district_vals)) {
  stop("Duplicate district IDs found")
}

# initialize empty set
#district_vals <- integer(0)

n_districts <- length(district_vals)

# Map district IDs → index
district_index <- setNames(seq_along(district_vals), district_vals)

# Safe indexing
#district_vals <- as.integer(district_vals)

# 3D arrays: [raster_i, raster_j, district]
n11_eco <- array(0, dim = c(n_rasters, n_rasters, n_districts))
n10_eco <- array(0, dim = c(n_rasters, n_rasters, n_districts))
n01_eco <- array(0, dim = c(n_rasters, n_rasters, n_districts))
n00_eco <- array(0, dim = c(n_rasters, n_rasters, n_districts))

# Count pixels where ecodistrict = NA
na_district_pixels <- 0

log_msg("Ecodistrict count arrays intialized")

## ---- Initialize CCS ---- ##
log_msg("Initializing CCS count arrays")

# Get unique IDs (could also delay until first block)
freq_table_ccs <- freq(ccs_rast)
ccs_vals <- as.integer(freq_table_ccs$value)

if (anyDuplicated(ccs_vals)) {
  stop("Duplicate CCS IDs found")
}

# initialize empty set
#district_vals <- integer(0)

n_ccs <- length(ccs_vals)

# Map CCS IDs → index
ccs_index <- setNames(seq_along(ccs_vals), ccs_vals)

# 3D arrays: [raster_i, raster_j, district]
n11_ccs <- array(0, dim = c(n_rasters, n_rasters, n_ccs))
n10_ccs <- array(0, dim = c(n_rasters, n_rasters, n_ccs))
n01_ccs <- array(0, dim = c(n_rasters, n_rasters, n_ccs))
n00_ccs <- array(0, dim = c(n_rasters, n_rasters, n_ccs))

# Count pixels where ecodistrict = NA
na_ccs_pixels <- 0

log_msg("CCS count arrays intialized")

## ---- Initialize global counts ---- ##
log_msg("Initializing global count arrays")

# Global counts 
n11_global <- matrix(0, n_rasters, n_rasters)
n10_global <- matrix(0, n_rasters, n_rasters)
n01_global <- matrix(0, n_rasters, n_rasters)
n00_global <- matrix(0, n_rasters, n_rasters)

log_msg("Global count arrays intialized")

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

# ✅ OPEN FILE CONNECTIONS
readStart(r_stack)
readStart(ecodist_rast)
readStart(ccs_rast)

on.exit({
  try(readStop(r_stack), silent = TRUE)
  try(readStop(ecodist_rast), silent = TRUE)
  try(readStop(ccs_rast), silent = TRUE)
}, add = TRUE)

make_binary_matrices <- function(x) {
  
  M <- x > 0
  valid <- !is.na(x)
  
  M_int <- M * 1
  valid_int <- valid * 1
  
  M_int[valid_int == 0] <- 0
  
  storage.mode(M_int) <- "double"
  storage.mode(valid_int) <- "double"
  
  list(
    M = M_int,
    valid = valid_int
  )
}

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
  
  # Count NA district & CCS pixels
  na_district_pixels <- na_district_pixels + sum(is.na(districts))
  na_ccs_pixels      <- na_ccs_pixels + sum(is.na(ccs))
  
  # Safe indexing
  districts <- as.integer(districts)
  ccs       <- as.integer(ccs)
  
  #storage.mode(vals) <- "integer"
  
  # Drop fully empty rows
  # valid_rows <- rowSums(!is.na(vals)) > 0 & !is.na(districts)
  # vals <- vals[valid_rows, , drop = FALSE]
  # districts <- districts[valid_rows]
  # 
  # if (nrow(vals) == 0) next
  # 
  # # Binary presence
  # M <- vals > 0
  # NA_mask <- is.na(vals)
  # 
  # valid_mat <- !NA_mask
  # M_int <- M * 1
  # valid_int <- valid_mat * 1
  # M_int[valid_int == 0] <- 0
  
  # Keep global pixels even if district/CCS is NA
  valid_rows_global  <- rowSums(!is.na(vals)) > 0
  valid_rows_ecodist <- valid_rows_global & !is.na(districts)
  valid_rows_ccs     <- valid_rows_global & !is.na(ccs)
  
  vals_global <- vals[valid_rows_global, , drop = FALSE]
  vals_eco    <- vals[valid_rows_ecodist, , drop = FALSE]
  vals_ccs    <- vals[valid_rows_ccs, , drop = FALSE]
  
  districts_block <- districts[valid_rows_ecodist]
  ccs_block       <- ccs[valid_rows_ccs]
  
  if (nrow(vals_global) == 0) {
    rm(vals, vals_global, vals_eco, vals_ccs, districts, districts_block, ccs, ccs_block)
    gc()
    next
  }
  
  global_obj <- make_binary_matrices(vals_global)
  M_global <- global_obj$M
  valid_global <- global_obj$valid
  
  eco_obj <- make_binary_matrices(vals_eco)
  M_eco <- eco_obj$M
  valid_eco <- eco_obj$valid
  
  ccs_obj <- make_binary_matrices(vals_ccs)
  M_ccs <- ccs_obj$M
  valid_ccs <- ccs_obj$valid
  
  #storage.mode(valid_int) <- "double"
  #storage.mode(M_int) <- "double"
  
  # ---- GLOBAL (all pixels in block) ----
  # N_block_all <- t(valid_int) %*% valid_int
  # n11_block_all <- t(M_int) %*% M_int
  # 
  # present_i_all <- t(M_int) %*% valid_int
  # present_j_all <- t(valid_int) %*% M_int
  
  N_block_all <- t(valid_global) %*% valid_global
  n11_block_all <- t(M_global) %*% M_global
  
  present_i_all <- t(M_global) %*% valid_global
  present_j_all <- t(valid_global) %*% M_global
  
  n10_block_all <- present_i_all - n11_block_all
  n01_block_all <- present_j_all - n11_block_all
  n00_block_all <- N_block_all - (n11_block_all + n10_block_all + n01_block_all)
  
  # accumulate
  n11_global <- n11_global + n11_block_all
  n10_global <- n10_global + n10_block_all
  n01_global <- n01_global + n01_block_all
  n00_global <- n00_global + n00_block_all
  
  if (any(n11_block_all < 0, na.rm = TRUE))
    stop("Negative n11 detected")
  
  if (any(n10_block_all < 0, na.rm = TRUE))
    stop("Negative n10 detected")
  
  if (any(n01_block_all < 0, na.rm = TRUE))
    stop("Negative n01 detected")
  
  if (any(n00_block_all < 0, na.rm = TRUE))
    stop("Negative n00 detected")
  
  if (i %% 50 == 0) {
    stopifnot(all(n11_global >= 0))
  }
  
  ## ---- Split by district ---- ##
  split_idx_eco <- split(seq_len(nrow(vals_eco)), districts_block)
  
  for (d in names(split_idx_eco)) {
    
    rows <- split_idx_eco[[d]]
    
    d_id <- district_index[as.character(d)]
    if (is.na(d_id)) {
      warning("Unknown district encountered: ", d)
      next
    }
    
    d_id <- as.integer(d_id)
    
    M_d <- M_eco[rows, , drop = FALSE]
    valid_d <- valid_eco[rows, , drop = FALSE]
    
    # ---- Matrix algebra ----
    N_block <- t(valid_d) %*% valid_d
    n11_block <- t(M_d) %*% M_d
    
    present_i <- t(M_d) %*% valid_d
    present_j <- t(valid_d) %*% M_d
    
    n10_block <- present_i - n11_block
    n01_block <- present_j - n11_block
    n00_block <- N_block - (n11_block + n10_block + n01_block)
    
    if (any(n00_block < 0, na.rm = TRUE)) {
      stop("Negative n00 detected")
    }
    
    # ---- Accumulate per district ----
    n11_eco[,,d_id] <- n11_eco[,,d_id] + n11_block
    n10_eco[,,d_id] <- n10_eco[,,d_id] + n10_block
    n01_eco[,,d_id] <- n01_eco[,,d_id] + n01_block
    n00_eco[,,d_id] <- n00_eco[,,d_id] + n00_block
  }
  
  ## ---- Split by CCS ---- ##
  
  if (nrow(vals_ccs) > 0) {
    
    split_idx_ccs <- split(
      seq_len(nrow(vals_ccs)),
      ccs_block
    )
    
    for (c in names(split_idx_ccs)) {
      
      rows <- split_idx_ccs[[c]]
      
      c_id <- ccs_index[as.character(c)]
      
      if (is.na(c_id)) {
        warning("Unknown CCS encountered: ", c)
        next
      }
      
      c_id <- as.integer(c_id)
      
      M_c <- M_ccs[rows, , drop = FALSE]
      valid_c <- valid_ccs[rows, , drop = FALSE]
      
      N_block <- t(valid_c) %*% valid_c
      n11_block <- t(M_c) %*% M_c
      
      present_i <- t(M_c) %*% valid_c
      present_j <- t(valid_c) %*% M_c
      
      n10_block <- present_i - n11_block
      n01_block <- present_j - n11_block
      n00_block <- N_block - (n11_block + n10_block + n01_block)
      
      n11_ccs[,,c_id] <- n11_ccs[,,c_id] + n11_block
      n10_ccs[,,c_id] <- n10_ccs[,,c_id] + n10_block
      n01_ccs[,,c_id] <- n01_ccs[,,c_id] + n01_block
      n00_ccs[,,c_id] <- n00_ccs[,,c_id] + n00_block
    }
  }
  
  rm(vals, vals_global, vals_eco, vals_ccs, districts, districts_block, 
     ccs, ccs_block, M_global, valid_global, M_eco, valid_eco, M_ccs, valid_ccs)
  
  if (exists("split_idx_eco", inherits = FALSE)) rm(split_idx_eco)
  if (exists("split_idx_ccs", inherits = FALSE)) rm(split_idx_ccs)
  if (i %% 10 == 0)
    gc()
}

readStop(r_stack)
readStop(ecodist_rast)
readStop(ccs_rast)

message("Finished block processing")
log_msg("Finished block processing")

message("Total pixels where district = NA: ", na_district_pixels)
log_msg("Total pixels where district = NA: ", na_district_pixels)

## ---- Compute metrics per district ---- ##
message("Computing metrics...")
log_msg("Computing metrics")

message(
  "First district_vals: ",
  paste(head(district_vals, 10), collapse = ",")
)

# ecodistrict results
ecodist_results <- list()

for (d in seq_len(n_districts)) {
  
  d_name <- district_vals[d]
  
  for (i in 1:(n_rasters - 1)) {
    for (j in (i + 1):n_rasters) {
      
      tp <- n11_eco[i,j,d]
      fp <- n10_eco[i,j,d]
      fn <- n01_eco[i,j,d]
      tn <- n00_eco[i,j,d]
      
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
      
      ecodist_results[[length(ecodist_results) + 1]] <- tibble(
        region_type = "ecodistrict",
        region_id = as.character(d_name),
        raster_i = ras_names[i],
        raster_j = ras_names[j],
        n11 = tp,
        n10 = fp,
        n01 = fn,
        n00 = tn,
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

# CCS results
ccs_results <- list()

for (c in seq_len(n_ccs)) {
  
  c_name <- ccs_vals[c]
  
  ccsuid <- ccs_lookup$CCSUID[match(c_name, ccs_lookup$ID)]
  
  if (is.na(ccsuid)) {
    stop(paste("Missing CCSUID for raster ID", c_name))}
  
  for (i in 1:(n_rasters - 1)) {
    for (j in (i + 1):n_rasters) {
      
      tp <- n11_ccs[i,j,c]
      fp <- n10_ccs[i,j,c]
      fn <- n01_ccs[i,j,c]
      tn <- n00_ccs[i,j,c]
      
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
      
      ccs_results[[length(ccs_results) + 1]] <- tibble(
        region_type = "ccs",
        region_id = as.character(ccsuid),
        raster_i = ras_names[i],
        raster_j = ras_names[j],
        n11 = tp,
        n10 = fp,
        n01 = fn,
        n00 = tn,
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
# global results
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
      region_type = as.character("global"),
      region_id = as.character("ALL"),
      raster_i = ras_names[i],
      raster_j = ras_names[j],
      n11 = tp,
      n10 = fp,
      n01 = fn,
      n00 = tn,
      d_jaccard = d_j,
      d_sorensen = d_s,
      kappa = kappa,
      f1 = f1,
      mcc = mcc,
      N = N
    )
  }
}

pairwise_long <- bind_rows(global_results, ecodist_results, ccs_results)

# pairwise_long <- bind_rows(
#   bind_rows(global_results) %>% mutate(district = as.character(district)),
#   bind_rows(ecodist_results) %>% mutate(district = as.character(district)))

log_msg("Done computing metrics")

##------ Final checks ----##
validate_spatial_agreement <- function(
    n11_global, n10_global, n01_global, n00_global,
    n11_eco, n10_eco, n01_eco, n00_eco,
    n11_ccs, n10_ccs, n01_ccs, n00_ccs,
    pairwise_long,
    district_vals,
    ccs_vals,
    tol = 1e-6
) {
  
  message("---- Running validation checks ----")
  
  ## ---- 1. Symmetry ---- ##
  #if (!all.equal(n11_global, t(n11_global))) {
  if (!isTRUE(all.equal(n11_global, t(n11_global)))) {
    stop("n11_global is not symmetric")
  }
  #if (!all.equal(n00_global, t(n00_global))) {
  if (!isTRUE(all.equal(n00_global, t(n00_global)))) {
    stop("n00_global is not symmetric")
  }
  #if (!all.equal(n10_global, t(n01_global))) {
  if (!isTRUE(all.equal(n10_global, t(n01_global)))) {
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
  if (length(dim(n11_eco)) == 3) {
    
    count_checks <- list(
      n11_eco = list(global = n11_global, district = n11_eco),
      n10_eco = list(global = n10_global, district = n10_eco),
      n01_eco = list(global = n01_global, district = n01_eco),
      n00_eco = list(global = n00_global, district = n00_eco)
    )
    
    for (nm in names(count_checks)) {
      
      district_sum <- apply(
        count_checks[[nm]]$district,
        c(1, 2),
        sum
      )
      
      diff_eco <- count_checks[[nm]]$global - district_sum
      max_diff_eco <- max(abs(diff_eco), na.rm = TRUE)
      
      message(
        sprintf(
          "%s: max difference (global - sum of districts) = %s",
          nm, max_diff_eco
        )
      )
      
      if (max_diff_eco > 0) {
        message("  (This should reflect pixels with NA district)")
      }
      
      if (any(diff_eco < 0)) {
        stop(sprintf("%s district total exceeds global total", nm))
      }
      
      # if (max_diff > na_district_pixels) {
      #   stop(sprintf("%s totals do not reconcile", nm))
      # }
    }
  }
  
  ## ---- 5. Global vs CCS totals ---- ##
  if (length(dim(n11_ccs)) == 3) {
    
    count_checks <- list(
      n11_ccs = list(global = n11_global, ccs = n11_ccs),
      n10_ccs = list(global = n10_global, ccs = n10_ccs),
      n01_ccs = list(global = n01_global, ccs = n01_ccs),
      n00_ccs = list(global = n00_global, ccs = n00_ccs)
    )
    
    for (nm in names(count_checks)) {
      
      ccs_sum <- apply(
        count_checks[[nm]]$ccs,
        c(1, 2),
        sum
      )
      
      diff_ccs <- count_checks[[nm]]$global - ccs_sum
      max_diff_ccs <- max(abs(diff_ccs), na.rm = TRUE)
      
      message(
        sprintf(
          "%s: max difference (global - sum of CCS) = %s",
          nm, max_diff_ccs
        )
      )
      
      if (max_diff_ccs > 0) {
        message("  (This should reflect pixels with NA CCS)")
      }
      
      if (any(diff_ccs < 0)) {
        stop(sprintf("%s CCS total exceeds global total", nm))
      }
      
      # if (max_diff > na_district_pixels) {
      #   stop(sprintf("%s totals do not reconcile", nm))
      # }
    }
  }
  
  ## ---- 6. Empty districts ---- ##
  if (length(dim(n11_eco)) == 3) {
    
    district_totals <- sapply(seq_along(district_vals), function(d) {
      sum(n11_eco[,,d] + n10_eco[,,d] + n01_eco[,,d] + n00_eco[,,d])
      #sum(n11[1,1,d] + n10[1,1,d] + n01[1,1,d] + n00[1,1,d])
    })
    
    empty_eco <- which(district_totals == 0)
    
    if (length(empty_eco) > 0) {
      warning("Empty districts detected: ",
              paste(district_vals[empty_eco], collapse = ", "))
    } else {
      message("✓ All districts contain data")
    }
  }
  
  ## ---- 7. Empty CCS ---- ##
  if (length(dim(n11_ccs)) == 3) {
    
    ccs_totals <- sapply(seq_along(ccs_vals), function(c) {
      sum(n11_ccs[,,c] + n10_ccs[,,c] + n01_ccs[,,c] + n00_ccs[,,c])
      #sum(n11[1,1,c] + n10[1,1,c] + n01[1,1,c] + n00[1,1,c])
    })
    
    empty_ccs <- which(ccs_totals == 0)
    
    if (length(empty_ccs) > 0) {
      warning("Empty CCS detected: ",
              paste(ccs_vals[empty_ccs], collapse = ", "))
    } else {
      message("✓ All CCS contain data")
    }
  }
  
  ## ---- 8. Metric bounds ---- ##
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
  
  ## ---- 9. Lookup integrity ---- ##
  missing_ids <- setdiff(
    ccs_vals,
    ccs_lookup$ID
  )
  
  if (length(missing_ids) > 0) {
    stop(
      paste(
        "CCS IDs missing from lookup:",
        paste(missing_ids, collapse = ", ")
      )
    )
  }
  
  message("✓ Lookup integrity checks passed")
  
  message("---- All validation checks completed ----")
}

log_msg("Running validation checks...")

validate_spatial_agreement(
  n11_global, n10_global, n01_global, n00_global,
  n11_eco, n10_eco, n01_eco, n00_eco,
  n11_ccs, n10_ccs, n01_ccs, n00_ccs,
  pairwise_long,
  district_vals,
  ccs_vals
)

## ---- Save outputs ---- ##
out_csv_local <- file.path(tmp_root, "pairwise_results.csv")
write_csv(pairwise_long, out_csv_local)

# aws.s3::put_object(
#   file = out_csv_local,
#   object = file.path(out_dir, "pairwise_binary_spatial_agreement_ecodist_CCS_16Oct2026.csv")
# )

# reusable function to upload files
upload_file_s3 <- function(local_file, s3_path) {
  
  s3_uri <- sub("^s3://", "", s3_path)
  
  bucket <- sub("/.*", "", s3_uri)
  key    <- sub("^[^/]*/", "", s3_uri)
  
  raw_data <- readBin(
    local_file,
    what = "raw",
    n = file.info(local_file)$size
  )
  
  s3$put_object(
    Bucket = bucket,
    Key = key,
    Body = raw_data
  )
}

upload_file_s3(out_csv_local, file.path(out_dir, "pairwise_binary_spatial_agreement_ecodist_CCS_16Oct2026.csv"))

# upload plot
upload_file_s3(facet_png, file.path(out_dir, "input_binary_rasters_faceted_16Oct2026.png"))

message("Finished pairwise spatial agreement analysis")
log_msg("Spatial agreement analysis completed successfully")

# Upload log to S3
# aws.s3::put_object(
#   file = log_file_local,
#   object = file.path(log_dir, basename(log_file_local))
# )

upload_file_s3(log_file_local, file.path(out_dir, basename(log_file_local)))

message("Log written to: ", log_file)
