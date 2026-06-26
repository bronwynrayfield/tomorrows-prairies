
###################################################################################
##          ANALYSIS: SPATIAL AGREEMENT (BINARY GRASS VS NON-GRASS)              ##
###################################################################################

##----------------------- Alison Long - April 2026 -------------------------

# step 6: calculate pairwise spatial agreement metrics to assess similarities and 
#         differences among existing grassland datasets. This version handles all 
#         grassland data products regardless of extent.

# ---- Set temp dirs BEFORE loading packages ----
tmp_root <- file.path(Sys.getenv("HOME"), "r_tmp")

dir.create(tmp_root, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(tmp_root, "terra"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(tmp_root, "curl"),  showWarnings = FALSE, recursive = TRUE)

Sys.setenv(
  TMPDIR = tmp_root,
  TEMP   = tmp_root,
  TMP    = tmp_root,
  CURL_TMPDIR = file.path(tmp_root, "curl")
)

## ---- Packages ---- ##
library(terra)           # raster engine
library(future.callr)    # safest parallel backend on SageMaker
library(future.apply)    # parallel lapply
library(readr)
library(dplyr)

## ---- Terra memory tuning ---- ##
terraOptions(
  memfrac = 0.6,          # allow terra to use most RAM
  progress = 1,
  tempdir = file.path(tmp_root, "terra") 
)

## ---- Parallel setup (safe on EC2/SageMaker) ---- ##
# Use physical cores only; callr avoids shared-memory crashes
#n_cores <- max(1, parallel::detectCores(logical = FALSE) - 10)
n_cores <- min(6, parallel::detectCores(logical = FALSE) - 2)
future::plan(callr, workers = n_cores)

message("Using ", n_cores, " workers")

## ---- Settings ---- ##
meta_file  <- "s3://sagemaker-gst-stage.sharing/alison-long/tomorrows-prairies/processed/metadata_full_extent_binary_with_path.csv"
out_dir    <- "s3://sagemaker-gst-stage.sharing/alison-long/tomorrows-prairies/results/spatial_agreement_analysis"
log_dir    <- "s3://sagemaker-gst-stage.sharing/alison-long/tomorrows-prairies/logs"

sample_n  <- 1e6          # max pixels sampled per raster pair (NULL = exact only)
max_exact <- 5e7          # max overlapping cells before sampling kicks in
no_overlap_dist <- 1      # if two rasters do not overlap at all, set their distance to this value

set.seed(42)

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

## ---- Read metadata ---- ##
log_msg("Reading metadata from: ", meta_file)

meta <- aws.s3::s3read_using(
  FUN = readr::read_csv,
  object = meta_file,
  show_col_types = FALSE
) %>%
  dplyr::rename(file = file_s3) # column in csv file is called file_s3 not file


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

## ---- Load rasters lazily (no values read yet) ---- ##
log_msg("Reading rasters")

#rlist <- lapply(meta$file, terra::rast)
#names(rlist) <- meta$raster

rfiles <- setNames(meta$file, meta$raster)

# Force binary validity (0/1/NA)
#rlist <- lapply(rlist, function(r) {
#  r <- round(r)
#  r[r < 0] <- 0
#  r[r > 1] <- 1
#  r
#})

#log_msg("Rasters loaded: ", length(rlist))

exists("rlist")

log_msg("Auditing raster geometry")

for (nm in names(rfiles)) {
  r <- terra::rast(rfiles[[nm]])
  log_msg(
    nm, " | ",
    "rows=", nrow(r),
    " cols=", ncol(r),
    " res=", paste(res(r), collapse = "x"),
    " CRS=", ifelse(is.na(crs(r)), "NA", "OK")
  )
}

## ---- Helper: metrics from confusion counts ----
.metrics_from_counts <- function(n00, n01, n10, n11) {
  
  N <- n00 + n01 + n10 + n11
  if (N == 0) {
    return(list(
      d_jaccard   = no_overlap_dist,
      d_sorensen  = no_overlap_dist,
      kappa       = NA_real_,
      N           = 0
    ))
  }
  
  # Jaccard (presence-only)
  denom_j <- n11 + n10 + n01
  d_jac <- if (denom_j == 0) 0 else 1 - (n11 / denom_j)
  
  # Sørensen–Dice
  denom_s <- 2 * n11 + n10 + n01
  d_sor <- if (denom_s == 0) 0 else 1 - (2 * n11 / denom_s)
  
  # Cohen's kappa
  po <- (n11 + n00) / N
  p1_yes <- (n11 + n10) / N
  p1_no  <- (n01 + n00) / N
  p2_yes <- (n11 + n01) / N
  p2_no  <- (n10 + n00) / N
  pe <- p1_yes * p2_yes + p1_no * p2_no
  kappa <- if (pe == 1) NA_real_ else (po - pe) / (1 - pe)
  
  list(
    d_jaccard  = d_jac,
    d_sorensen = d_sor,
    kappa      = kappa,
    N          = N
  )
}

## ---- Helper: fast pairwise counts ----
.fast_pair_counts <- function(file1, file2) {
  
  r1 <- terra::rast(file1)
  r2 <- terra::rast(file2)
  
  # Quick spatial reject
  if (!terra::relate(ext(r1), ext(r2), "intersects")) {
    return(list(ok = FALSE))
  }
  
  e <- terra::intersect(ext(r1), ext(r2))
  if (is.null(e)) return(list(ok = FALSE))
  
  # Crop geometry only (fast, no values read)
  r1o <- crop(r1, e)
  r2o <- crop(r2, e)
  
  # Align grids if needed
  if (!compareGeom(r1o, r2o, stopOnError = FALSE)) {
    r2o <- resample(r2o, r1o, method = "near")
  }
  
  s <- c(r1o, r2o)
  names(s) <- c("a", "b")
  
  # Mask to cells where both rasters have data
  s <- mask(s, !is.na(s$a) & !is.na(s$b), maskvalues = FALSE)
  
  n_over <- global(!is.na(s$a), "sum", na.rm = TRUE)[1,1]
  if (n_over == 0) return(list(ok = FALSE))
  
  # Decide: exact or sampled
  do_exact <- is.null(sample_n) || n_over <= max_exact
  
  if (do_exact) {
    tab <- terra::crosstab(s, useNA = FALSE, long = FALSE)
    
    get <- function(r, c) if (r %in% rownames(tab) && c %in% colnames(tab)) tab[r, c] else 0
    
    return(list(
      ok  = TRUE,
      n00 = get("0", "0"),
      n01 = get("0", "1"),
      n10 = get("1", "0"),
      n11 = get("1", "1")
    ))
  }
  
  # Sampling path (very large overlaps)
  df <- spatSample(
    s, size = sample_n,
    method = "random",
    na.rm = TRUE,
    as.df = TRUE
  )
  
  df$a <- as.integer(df$a > 0)
  df$b <- as.integer(df$b > 0)
  
  tab <- table(df$a, df$b)
  
  list(
    ok  = TRUE,
    n00 = ifelse("0" %in% rownames(tab) && "0" %in% colnames(tab), tab["0","0"], 0),
    n01 = ifelse("0" %in% rownames(tab) && "1" %in% colnames(tab), tab["0","1"], 0),
    n10 = ifelse("1" %in% rownames(tab) && "0" %in% colnames(tab), tab["1","0"], 0),
    n11 = ifelse("1" %in% rownames(tab) && "1" %in% colnames(tab), tab["1","1"], 0)
  )
}

## ---- Pairwise computation ----
ras_names <- names(rfiles)
n <- length(rfiles)
pairs <- combn(n, 2, simplify = FALSE)

log_msg("Raster pairs to compute: ", length(pairs))

results <- future_lapply(pairs, function(p) {
  
  i <- p[1]; j <- p[2]
  cnt <- .fast_pair_counts(rfiles[[i]], rfiles[[j]])
  
  if (!cnt$ok) {
    return(tibble(
      raster_i = ras_names[i],
      raster_j = ras_names[j],
      d_jaccard = no_overlap_dist,
      d_sorensen = no_overlap_dist,
      kappa = NA_real_,
      N = 0
    ))
  }
  
  m <- .metrics_from_counts(cnt$n00, cnt$n01, cnt$n10, cnt$n11)
  
  tibble(
    raster_i   = ras_names[i],
    raster_j   = ras_names[j],
    d_jaccard  = m$d_jaccard,
    d_sorensen = m$d_sorensen,
    kappa      = m$kappa,
    N          = m$N
  )
})

pairwise_long <- bind_rows(results)

log_msg("Pairwise computation completed")

## ---- Save outputs ----
out_csv <- file.path(out_dir, "pairwise_binary_spatial_agreement.csv")
out_rds <- file.path(out_dir, "pairwise_binary_spatial_agreement.rds")

write_csv(pairwise_long, out_csv)
saveRDS(pairwise_long, out_rds)

log_msg("Outputs written:")
log_msg(" - ", out_csv)
log_msg(" - ", out_rds)

message("Finished pairwise spatial agreement analysis")
log_msg("Spatial agreement analysis completed successfully")

# Upload log to S3
aws.s3::put_object(
  file = log_file_local,
  object = log_file
)

message("Log written to: ", log_file)

