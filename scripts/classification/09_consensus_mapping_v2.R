
###################################################################################
##                     GRASSLAND CONSENSUS MAPPING ANALYSIS                      ##
###################################################################################

##----------------------- Alison Long - June 2026 -------------------------

# step 9: use data fusion techniques to create a consensus grassland map 
#         to serve as baseline for modelling

# THIS VERSION IS OPTIMIZED FOR SAGEMAKER

# NEED TO ADD ECOREGION-LEVEL ACCURACY FALLBACK FOR ECODISTRICTS W/O ACCURACY DATA
# (or will this be needed once we have crop ins data?)

## ---- Setup ---- ##
# Temp dirs
tmp_root <- file.path(Sys.getenv("HOME"), "r_tmp")
#dir.create(tmp_root, recursive = TRUE, showWarnings = FALSE)

# Load packages
library(terra)
library(future.apply)
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

## ---- Parallel (light use only) ---- ##
#n_cores <- min(6, parallel::detectCores(logical = FALSE) - 2)
future::plan(sequential)  # Avoids disk contention

## ---- Settings ---- ##

out_dir      <- "s3://sagemaker-gst-stage.sharing/alison-long/tomorrows-prairies/results/consensus_map"
log_dir      <- "s3://sagemaker-gst-stage.sharing/alison-long/tomorrows-prairies/logs"
terra_dir    <- file.path(tmp_root, "terra")
local_dir    <- file.path(tmp_root, "data")

meta_file    <- "s3://sagemaker-gst-stage.sharing/alison-long/tomorrows-prairies/processed/metadata_full_extent_binary_nattam_with_path.csv"
ecodist_file <- "s3://sagemaker-gst-stage.sharing/alison-long/tomorrows-prairies/processed/ecodistricts_clean_rast.tif"
acc_bin_file <- "s3://sagemaker-gst-stage.sharing/alison-long/tomorrows-prairies/processed/accuracy_metrics/binary_accuracy_by_ecodistrict_grouped_combinedCIandCM_23Jun2026.csv"
acc_nt_file  <- "s3://sagemaker-gst-stage.sharing/alison-long/tomorrows-prairies/processed/accuracy_metrics/nattam_accuracy_by_ecodistrict_grouped_lenient_23Jun2026.csv"

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

#--------- Build MCC weight rasters (binary stack) -------#
log_msg("Building MCC weight stack (binary)")

# Load accuracy (MCC)
acc_bin <- aws.s3::s3read_using(read_csv, object = acc_bin_file, show_col_types = FALSE) %>%
  mutate(raster_name = tolower(trimws(raster)))

stopifnot(all(c("raster_name","Ecds_ID","MCC") %in% names(acc_bin)))

# Convert MCC → weight [0,1]
mcc_to_weight <- function(x) {
  w <- (x + 1) / 2
  w[!is.finite(w)] <- NA
  pmax(w, 1e-4)
}

acc_bin <- acc_bin %>%
  mutate(weight = mcc_to_weight(MCC))

acc_bin_wide <- acc_bin %>%
  select(Ecds_ID, raster_name, weight) %>%
  pivot_wider(names_from = raster_name, values_from = weight)

# check for missing columns or naming mismatch
missing <- setdiff(names(bin_stack), acc_bin$raster_name)
if (length(missing) > 0) {
  stop("Missing MCC values for: ", paste(missing, collapse=", "))
}

# Load rasterized ecodistricts 
ecodist_r <- rast(ecodist_file)

# check ecodistrict raster values vs IDs
print(head(freq(ecodist_r)))

ids_r <- unique(freq(ecodist_r)[,1])
ids_tbl <- unique(acc_bin$Ecds_ID)

missing_ids <- setdiff(ids_r, ids_tbl)

if (length(missing_ids) > 0) {
  stop("Ecodistrict IDs missing from accuracy table: ",
       paste(head(missing_ids), collapse=", "))
}

# Build weight stack layer-by-layer
weight_files <- lapply(names(bin_stack), function(nm) {
  
  log_msg("Rasterizing weight for: ", nm)
  
  lookup_bin <- acc_bin %>%
    filter(raster_name == nm) %>%
    select(Ecds_ID, weight) %>%
    as.matrix()
  
  if (nrow(lookup_bin) == 0) {
    stop("No weights found for raster: ", nm)
  }
  
  w <- classify(ecodist_r, rcl = lookup_bin)
  
  out <- file.path(terra_dir, paste0("w_", nm, ".tif"))
  
  writeRaster(w, out, overwrite = TRUE,
              wopt = list(datatype="FLT4S"))
  
  out
})

weight_stack <- rast(weight_files)
names(weight_stack) <- names(bin_stack)

#--------- MCC-weighted positive evidence -------#
log_msg("Running Strategy A (MCC-weighted fusion)")

present <- !is.na(bin_stack)
pos_mask <- bin_stack == 1
neg_mask <- bin_stack == 0

# Zero-out weights where no data
weight_present <- weight_stack * present

# Numerator
numerator_A <- app(
  c(pos_mask, weight_present),
  fun = function(x) {
    n <- length(x) / 2
    pos <- x[1:n]
    w   <- x[(n+1):(2*n)]
    sum(pos * w, na.rm = TRUE)
  },
  filename = file.path(terra_dir, "numerator_A.tif"),
  overwrite = TRUE
)

# Denominator
denominator_A <- app(
  weight_present,
  sum,
  na.rm = TRUE,
  filename = file.path(terra_dir, "denominator_A.tif"),
  overwrite = TRUE
)

# Probability 
prob_A <- numerator_A / denominator_A
prob_A[denominator_A == 0] <- NA # division safety

prob_A <- clamp(
  prob_A, 0, 1,
  filename = file.path(terra_dir, "prob_A.tif"),
  overwrite = TRUE
)

# Diagnostics
num_contributing <- app(present, sum, na.rm=TRUE,
                        filename=file.path(terra_dir, "num_contributing.tif"),
                        overwrite=TRUE
)

num_positive <- app(pos_mask, sum, na.rm=TRUE,
                    filename=file.path(terra_dir, "num_positive.tif"),
                    overwrite=TRUE
)

num_negative <- app(neg_mask, sum, na.rm=TRUE,
                    filename=file.path(terra_dir, "num_negative.tif"),
                    overwrite=TRUE
)

sum_weight_present <- denominator_A

# Uncertainty
uncertainty_prob_A <- prob_A * (1 - prob_A)

uncertainty_data_A <- app(
  sum_weight_present,
  fun = function(x) 1 / (1 + x),
  filename = file.path(terra_dir, "uncertainty_data.tif"),
  overwrite = TRUE
)

uncertainty_composite_A <- uncertainty_prob_A * (1 + uncertainty_data_A)

# Minimum support
min_sources <- 2

prob_A_minN <- ifel(
  num_contributing >= min_sources,
  prob_A, NA,
  filename = file.path(terra_dir, "prob_A_minN.tif"),
  overwrite = TRUE
)

#--------- Tri-class (MCC + F1 weights) ---------#
log_msg("Running tri-class fusion")

acc_nt <- aws.s3::s3read_using(read_csv, object = acc_nt_file, show_col_types = FALSE) %>%
  mutate(raster_name = tolower(trimws(raster)))

acc_nt <- acc_nt %>%
  mutate(
    MCC_1v2 = MCC_Multiclass
    #MCC_1v2 = (MCC_OvR_1 + MCC_OvR_2) / 2, # use multiclass for now, fix later
    F1_1_2x2 = F1_1,
    F1_2_2x2 = F1_2
  )

# Required columns
req <- c("raster_name","Ecds_ID","MCC_1v2","F1_2_2x2","F1_1_2x2")
stopifnot(all(req %in% names(acc_nt)))

# Build weights
acc_nt <- acc_nt %>%
  mutate(
    W_native = pmax((MCC_1v2 + 1)/2 * F1_2_2x2, 0),
    W_tame   = pmax((MCC_1v2 + 1)/2 * F1_1_2x2, 0)
  )

# Pivot
acc_native <- acc_nt %>%
  select(Ecds_ID, raster_name, W_native) %>%
  pivot_wider(names_from=raster_name, values_from=W_native)

acc_tame <- acc_nt %>%
  select(Ecds_ID, raster_name, W_tame) %>%
  pivot_wider(names_from=raster_name, values_from=W_tame)

# Rasterize weights (fast categorical lookup)
native_files <- lapply(names(nt_stack), function(nm) {
  out <- file.path(terra_dir, paste0("nt_native_", nm, ".tif"))
  
  lookup_native <- acc_nt %>%
    filter(raster_name == nm) %>%
    select(Ecds_ID, W_native) %>%
    as.matrix()
  
  if (nrow(lookup_native) == 0) {
    stop("No native weights for raster: ", nm)
  }
  
  w <- classify(ecodist_r, rcl = lookup_native)
  
  writeRaster(w, out, overwrite=TRUE)
  out
})

tame_files <- lapply(names(nt_stack), function(nm) {
  out <- file.path(terra_dir, paste0("nt_tame_", nm, ".tif"))
  
  lookup_tame <- acc_nt %>%
    filter(raster_name == nm) %>%
    select(Ecds_ID, W_tame) %>%
    as.matrix()
  
  if (nrow(lookup_tame) == 0) {
    stop("No tame weights for raster: ", nm)
  }
  
  w <- classify(ecodist_r, rcl = lookup_tame)
  
  writeRaster(w, out, overwrite=TRUE)
  out
})

metric_native_stack <- rast(native_files)
metric_tame_stack   <- rast(tame_files)

# Native/tame voting
tri_present <- !is.na(nt_stack)
native_mask <- nt_stack == 2
tame_mask   <- nt_stack == 1

native_votes <- app(
  c(native_mask, metric_native_stack),
  fun = function(x) {
    n <- length(x)/2
    sum(x[1:n] * x[(n+1):(2*n)], na.rm=TRUE)
  },
  filename=file.path(terra_dir,"native_votes.tif"),
  overwrite=TRUE
)

tame_votes <- app(
  c(tame_mask, metric_tame_stack),
  fun = function(x) {
    n <- length(x)/2
    sum(x[1:n] * x[(n+1):(2*n)], na.rm=TRUE)
  },
  filename=file.path(terra_dir,"tame_votes.tif"),
  overwrite=TRUE
)

# Shares
vote_total <- native_votes + tame_votes

share_native <- native_votes / vote_total
share_native[vote_total == 0] <- NA
share_native <- clamp(share_native, 0, 1)

share_tame <- 1 - share_native

#--------- Final probabilities -----------#
log_msg("Combining binary + tri-class outputs")

# Probabilities
p_grass <- prob_A

p_native <- p_grass * share_native
p_tame   <- p_grass * share_tame
p_not    <- 1 - p_grass

# Hard classification
p_stack <- c(p_not, p_tame, p_native)

class_map <- which.max(p_stack)

rcl <- matrix(c(1,0,2,1,3,2), ncol=2, byrow=TRUE)

#----------- Uncertainty ----------#
max_p <- app(p_stack, max,
             filename=file.path(terra_dir,"maxp.tif"),
             overwrite=TRUE
)

uncertainty <- 1 - max_p

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

# Write locally first, then save to S3
prob_A_local <- file.path(tmp_root, "prob_A.tif")
writeRaster(prob_A, prob_A_local, overwrite=TRUE, wopt = wopt)
aws.s3::put_object(
  file   = prob_A_local,
  object = paste0(out_dir, "/prob_grass_MCC.tif"))

prob_A_minN_local <- file.path(tmp_root, "prob_grass_minN.tif")
writeRaster(prob_A_minN, prob_A_minN_local, overwrite=TRUE, wopt = wopt)
aws.s3::put_object(
  file   = prob_A_minN_local,
  object = paste0(out_dir, "/prob_grass_minN.tif"))

p_native_local <- file.path(tmp_root, "p_native.tif")
writeRaster(p_native, p_native_local, overwrite=TRUE, wopt = wopt)
aws.s3::put_object(
  file   = p_native_local,
  object = paste0(out_dir, "/p_native.tif"))

p_tame_local <- file.path(tmp_root, "p_tame.tif")
writeRaster(p_tame, p_tame_local, overwrite=TRUE, wopt = wopt)
aws.s3::put_object(
  file   = p_tame_local,
  object = paste0(out_dir, "/p_tame.tif"))

p_not_local <- file.path(tmp_root, "p_not_grass.tif")
writeRaster(p_not, p_not_local, overwrite=TRUE, wopt = wopt)
aws.s3::put_object(
  file   = p_not_local,
  object = paste0(out_dir, "/p_not_grass.tif"))

uncertainty_local <- file.path(tmp_root, "uncertainty.tif")
writeRaster(uncertainty, uncertainty_local, overwrite=TRUE, wopt = wopt)
aws.s3::put_object(
  file   = uncertainty_local,
  object = paste0(out_dir, "/uncertainty.tif"))

num_contributing_local <- file.path(tmp_root, "num_sources.tif")
writeRaster(num_contributing, num_contributing_local, overwrite=TRUE, wopt = wopt_int)
aws.s3::put_object(
  file   = num_contributing_local,
  object = paste0(out_dir, "/num_sources.tif"))

num_positive_local <- file.path(tmp_root, "num_positive.tif")
writeRaster(num_positive, num_positive_local, overwrite=TRUE, wopt = wopt_int)
aws.s3::put_object(
  file   = num_positive_local,
  object = paste0(out_dir, "/num_positive.tif"))

num_negative_local <- file.path(tmp_root, "num_negative.tif")
writeRaster(num_negative, num_negative_local, overwrite=TRUE, wopt = wopt_int)
aws.s3::put_object(
  file   = num_negative_local,
  object = paste0(out_dir, "/num_negative.tif"))

class_3cat <- classify(class_map, rcl)

class_3cat_local <- file.path(tmp_root, "grass_3class.tif")
writeRaster(class_3cat, class_3cat_local, overwrite=TRUE, wopt = wopt_class)
aws.s3::put_object(
  file   = class_3cat_local,
  object = paste0(out_dir, "/grass_3class.tif"))

# upload log to S3
aws.s3::put_object(
  file = log_file_local,
  object = log_file
)
