
###################################################################################
##          ANALYSIS: SPATIAL AGREEMENT (BINARY GRASS VS NON-GRASS)              ##
###################################################################################

##----------------------- Alison Long - May 2026 -------------------------

# step 6: calculate pairwise spatial agreement metrics to assess similarities and 
#         differences among existing grassland datasets. This version handles all 
#         grassland data products regardless of extent.

## Optimizations vs original (v3): 
# - Rasters copied from S3 to local NVMe once (avoids repeated S3 reads per pair) 
# - mask() removed — crosstab handles NA exclusion natively 
# - compareGeom/resample removed — input rasters are pre-aligned 
# - ncell() used for overlap size estimate (header-only, no pixel read) 
# - n_over derived from crosstab sum (no extra global() read) 
# - Bitwise integer counts replace table() in sampling path 
# - memfrac lowered to prevent worker memory contention under future.callr 

# ---- Set temp dirs BEFORE loading packages ---- 
tmp_root <- file.path(Sys.getenv("HOME"), "r_tmp") 
dir.create(tmp_root, showWarnings = FALSE, recursive = TRUE) 
dir.create(file.path(tmp_root, "terra"), showWarnings = FALSE, recursive = TRUE) 
dir.create(file.path(tmp_root, "curl"), showWarnings = FALSE, recursive = TRUE) 
Sys.setenv(TMPDIR = tmp_root, TEMP = tmp_root, TMP = tmp_root, CURL_TMPDIR = file.path(tmp_root, "curl")) 

## ---- Packages ---- ## 
library(terra) 
library(future.callr) 
library(future.apply) 
library(readr) 
library(dplyr) 
library(aws.s3) 

## ---- Terra memory tuning ---- ## 
# Each future.callr worker is a separate process with its own terra memory pool. 
# Total potential usage = n_cores x memfrac x RAM, so keep memfrac modest. 
# Adjust based on instance RAM (example assumes ~64 GB; raise for larger instances). 
terraOptions( memfrac = 0.10, # ~6 GB per worker x 6 workers = ~36 GB total 
              progress = 0, # suppress per-tile progress in parallel workers 
              tempdir = file.path(tmp_root, "terra") ) 

## ---- Parallel setup ---- ## 
n_cores <- min(6, parallel::detectCores(logical = FALSE) - 2) 
future::plan(callr, workers = n_cores) 
message("Using ", n_cores, " workers") 

## ---- Settings ---- ## 
meta_file <- "s3://sagemaker-gst-stage.sharing/alison-long/tomorrows-prairies/processed/metadata_full_extent_binary_with_path.csv" 
out_dir <- "s3://sagemaker-gst-stage.sharing/alison-long/tomorrows-prairies/results/spatial_agreement_analysis" 
log_dir <- "s3://sagemaker-gst-stage.sharing/alison-long/tomorrows-prairies/logs" 
local_raster_dir <- "/tmp/rasters" # local NVMe — fast disk, survives the session 

sample_n <- 1e6 # pixels sampled per pair when overlap exceeds max_exact 
max_exact <- 5e7 # overlap cell threshold above which sampling is used 
no_overlap_dist <- 1 # distance assigned when two rasters do not overlap 

set.seed(42) 

## ---- Local log file (uploaded to S3 at end) ---- ## 
log_file_local <- file.path( 
  Sys.getenv("HOME"), 
  paste0("spatial_agreement_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log") ) 

log_file_s3 <- file.path(log_dir, paste0( "spatial_agreement_log_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log" )) 

log_msg <- function(..., level = "INFO") { 
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S") 
  msg <- paste0("[", level, "] ", ts, " : ", paste(..., collapse = "")) 
  message(msg) 
  cat(msg, "\n", file = log_file_local, append = TRUE) } 

log_msg("Starting spatial agreement analysis") 
log_msg("Parallel workers: ", n_cores) 

## ---- Read metadata ---- ## 
log_msg("Reading metadata from: ", meta_file) 

meta <- aws.s3::s3read_using( 
  FUN = readr::read_csv, 
  object = meta_file, 
  show_col_types = FALSE 
  ) %>% 
  dplyr::rename(file = file_s3) 

required_cols <- c("raster", "file") 
if (!all(required_cols %in% names(meta))) { 
  log_msg("Metadata missing required columns", level = "ERROR") 
  stop("Metadata must contain: ", paste(required_cols, collapse = ", ")) 
  } 

if (anyDuplicated(meta$raster)) { 
  log_msg("Duplicate raster names detected in metadata", level = "ERROR") 
  stop("Raster names must be unique") 
  } 

log_msg("Metadata rows loaded: ", nrow(meta)) 

## ---- Copy rasters from S3 to local disk (once) ---- ## 
# Reading each raster from S3 on every pairwise operation is the dominant 
# bottleneck at this scale. Copying all files locally means each of the 
# ~400 pairs reads from fast local NVMe instead of repeatedly hitting S3. 
log_msg("Copying rasters from S3 to local disk: ", local_raster_dir) 

dir.create(local_raster_dir, showWarnings = FALSE, recursive = TRUE) 

rfiles_local <- setNames( 
  vapply(seq_len(nrow(meta)), function(i) { 
    nm <- meta$raster[i] 
    src <- meta$file[i] 
    dest <- file.path(local_raster_dir, paste0(nm, ".tif")) 
    
    if (!file.exists(dest)) { 
      log_msg(" Downloading: ", nm) 
      aws.s3::save_object(object = src, file = dest) 
      } else { 
        log_msg(" Already cached: ", nm) 
        } 
    dest 
    }, character(1)), 
  meta$raster 
  ) 

log_msg("All rasters available locally") 

## ---- Audit raster geometry ---- ## 
log_msg("Auditing raster geometry") 

for (nm in names(rfiles_local)) { 
  r <- terra::rast(rfiles_local[[nm]]) 
  log_msg( 
    nm, " | ", 
    "rows=", nrow(r), 
    " cols=", ncol(r), 
    " res=", paste(res(r), collapse = "x"), 
    " CRS=", ifelse(is.na(crs(r)), "NA", "OK") 
    ) 
  } 

## ---- Helper: confusion counts → metrics ---- ## 
.metrics_from_counts <- function(n00, n01, n10, n11) { 
  N <- n00 + n01 + n10 + n11 
  if (N == 0) { 
    return(list( 
      d_jaccard = no_overlap_dist, 
      d_sorensen = no_overlap_dist, 
      kappa = NA_real_, 
      N = 0L 
      )) 
    } 
  
  # Jaccard distance (presence-only) 
  denom_j <- n11 + n10 + n01 
  d_jac <- if (denom_j == 0) 0 else 1 - (n11 / denom_j) 
  
  # Sørensen–Dice distance 
  denom_s <- 2 * n11 + n10 + n01 
  d_sor <- if (denom_s == 0) 0 else 1 - (2 * n11 / denom_s) 
  
  # Cohen's kappa 
  po <- (n11 + n00) / N 
  p1_yes <- (n11 + n10) / N 
  p1_no <- (n01 + n00) / N 
  p2_yes <- (n11 + n01) / N 
  p2_no <- (n10 + n00) / N 
  pe <- p1_yes * p2_yes + p1_no * p2_no 
  kappa <- if (pe == 1) NA_real_ else (po - pe) / (1 - pe) 
  
  list( 
    d_jaccard = d_jac, 
    d_sorensen = d_sor, 
    kappa = kappa, 
    N = N 
    ) 
  } 

## ---- Helper: fast pairwise confusion counts ---- ## 
.fast_pair_counts <- function(file1, file2) { 
  r1 <- terra::rast(file1) 
  r2 <- terra::rast(file2) 
  
  # Header-only spatial reject — no pixels read 
  if (!terra::relate(ext(r1), ext(r2), "intersects")) { 
    return(list(ok = FALSE)) 
    } 
  
  e <- terra::intersect(ext(r1), ext(r2)) 
  if (is.null(e)) return(list(ok = FALSE)) 
  
  # Crop to intersection (lazy — pixel values not yet loaded) 
  r1o <- terra::crop(r1, e) 
  r2o <- terra::crop(r2, e) 
  
  # Input rasters are pre-aligned: skip compareGeom / resample entirely. 
  s <- c(r1o, r2o) 
  names(s) <- c("a", "b") 
  
  # Use ncell() as a conservative upper bound for the overlap size check. 
  # It reads only the header, not pixel values. 
  n_cells_approx <- terra::ncell(r1o) 
  
  # ---- Sampling path (very large overlaps) ---- 
  if (!is.null(sample_n) && n_cells_approx > max_exact) { 
    df <- terra::spatSample( 
      s, 
      size = sample_n, 
      method = "random", 
      na.rm = TRUE, # drops cells where either layer is NA 
      as.df = TRUE 
      ) 
    
    if (nrow(df) == 0) return(list(ok = FALSE)) 
    
    a <- as.integer(df$a) 
    b <- as.integer(df$b) 
    
    # Bitwise integer comparisons — faster than table() at 1 M rows 
    n11 <- sum(a == 1L & b == 1L) 
    n10 <- sum(a == 1L & b == 0L) 
    n01 <- sum(a == 0L & b == 1L) 
    n00 <- sum(a == 0L & b == 0L) 
    
    # ---- Exact path ---- 
    } else { 
      
      # crosstab excludes NAs natively (useNA = FALSE): 
      # no mask() step needed. 
      tab <- terra::crosstab(s, useNA = FALSE, long = FALSE) 
      
      get_cell <- function(r, c) 
        if (r %in% rownames(tab) && c %in% colnames(tab)) tab[r, c] else 0L 
      
      n00 <- get_cell("0", "0") 
      n01 <- get_cell("0", "1") 
      n10 <- get_cell("1", "0") 
      n11 <- get_cell("1", "1") 
      } 
  
  # n_over derived from counts — no extra global() read required 
  n_over <- n00 + n01 + n10 + n11 
  if (n_over == 0) return(list(ok = FALSE)) 
  list(ok = TRUE, n00 = n00, n01 = n01, n10 = n10, n11 = n11) 
  } 

## ---- Pairwise computation ---- ## 
ras_names <- names(rfiles_local) 
n <- length(rfiles_local) 
pairs <- combn(n, 2, simplify = FALSE) 

log_msg("Raster count: ", n) 
log_msg("Raster pairs to compute: ", length(pairs)) 

results <- future_lapply(pairs, function(p) { 
  i <- p[1]; j <- p[2] 
  cnt <- .fast_pair_counts(rfiles_local[[i]], rfiles_local[[j]]) 
  if (!cnt$ok) { 
    return(tibble::tibble( 
      raster_i = ras_names[i], 
      raster_j = ras_names[j], 
      d_jaccard = no_overlap_dist, 
      d_sorensen = no_overlap_dist, 
      kappa = NA_real_, 
      N = 0L 
      )) 
    } 
  
  m <- .metrics_from_counts(cnt$n00, cnt$n01, cnt$n10, cnt$n11) 
  
  tibble::tibble( 
    raster_i = ras_names[i], 
    raster_j = ras_names[j], 
    d_jaccard = m$d_jaccard, 
    d_sorensen = m$d_sorensen, 
    kappa = m$kappa, 
    N = m$N 
    ) 
  }, future.seed = TRUE) 

pairwise_long <- dplyr::bind_rows(results) 

log_msg("Pairwise computation completed. Rows: ", nrow(pairwise_long)) 

## ---- Save outputs ---- ## 
out_csv <- file.path(out_dir, "pairwise_binary_spatial_agreement.csv") 
out_rds <- file.path(out_dir, "pairwise_binary_spatial_agreement.rds") 

# Write locally first, then push to S3 
local_csv <- file.path(tmp_root, "pairwise_binary_spatial_agreement.csv") 
local_rds <- file.path(tmp_root, "pairwise_binary_spatial_agreement.rds") 

readr::write_csv(pairwise_long, local_csv) 
saveRDS(pairwise_long, local_rds) 

aws.s3::put_object(file = local_csv, object = out_csv) 
aws.s3::put_object(file = local_rds, object = out_rds) 

log_msg("Outputs written:") 
log_msg(" - ", out_csv) 
log_msg(" - ", out_rds) 

message("Finished pairwise spatial agreement analysis.") 
log_msg("Spatial agreement analysis completed successfully") 

## ---- Upload log to S3 ---- ## 
aws.s3::put_object(file = log_file_local, object = log_file_s3) 
message("Log written to: ", log_file_s3)
