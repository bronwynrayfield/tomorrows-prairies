
###################################################################################
##                         ACCURACY ASSESSMENT - BINARY                          ##
###################################################################################

##----------------------- Alison Long - February 2026 -------------------------

# step 6: assess accuracy of datasets mapping grasslands by comparing them
#         to ground-truthed point data

# ============================================================
# grass / not grass raster accuracy assessment (year-matched)
# WITH per-raster exclusion of training points based on source and/or point IDs
#
# Methods included:
#   Method 1: Strict point extraction (pixel at point)
#   Method 2: Buffered extraction + majority vote
#   Method 3: Overall summaries across rasters (macro + weighted)
#
# Outputs:
#   - matched points with predictions (method 1 & 2)
#   - accuracy metrics per raster (method 1 & 2)
#   - optional by-source summaries (assessment sources only)
# ============================================================

# -----------------------------
# 0) PACKAGES
# -----------------------------
pkgs <- c("terra", "sf", "dplyr", "readr", "lubridate", "purrr", "tidyr", "yardstick", "stringr")
to_install <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(to_install) > 0) install.packages(to_install)

library(terra)
library(sf)
library(dplyr)
library(readr)
library(lubridate)
library(purrr)
library(tidyr)
library(yardstick)
library(stringr)

# -----------------------------
# 1) DEFINE SETTINGS AND INPUTS
# -----------------------------
raster_specs_csv <- "path/to/raster_specs.csv"          # raster_id,file,year,exclude_sources,exclude_point_ids
points_shp       <- "path/to/groundtruth_points.shp"    # shapefile

# Column names in shapefile:
truth_col  <- "grass_obs"   # 0/1 ground truth
date_col   <- "obs_date"    # date/datetime/character
source_col <- "source"      # REQUIRED for source filtering
id_col     <- NULL          # e.g. "pt_id" (set NULL if no have stable IDs)

# Buffer distances for Method 2 (meters). Use c(0) to mimic strict extraction,
# or include multiple to test sensitivity.
buffer_m_values <- c(30)

# Output folder:
out_dir <- "accuracy_outputs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# 2) HELPER FUNCTIONS
# -----------------------------

parse_obs_datetime <- function(x) {
  if (inherits(x, "Date") || inherits(x, "POSIXt")) return(x)
  dt <- ymd_hms(x, quiet = TRUE, tz = "UTC")
  dt <- if_else(is.na(dt), ymd(x, quiet = TRUE, tz = "UTC"), dt)
  dt
}

split_semicolon <- function(x) {
  # returns character vector; handles NA/blank gracefully
  if (length(x) == 0 || is.na(x) || trimws(x) == "") return(character(0))
  str_split(x, pattern = ";", simplify = FALSE)[[1]] %>%
    trimws() %>%
    (\(v) v[v != ""])()
}

majority_vote <- function(vals) {
  vals <- vals[!is.na(vals)]
  if (length(vals) == 0) return(NA_integer_)
  as.integer(names(which.max(table(vals))))
}

calc_metrics <- function(df, truth_col = "truth", pred_col = "pred") {
  df <- df %>%
    filter(!is.na(.data[[truth_col]]), !is.na(.data[[pred_col]])) %>%
    mutate(
      truth = factor(.data[[truth_col]], levels = c(0, 1)),
      pred  = factor(.data[[pred_col]],  levels = c(0, 1))
    )
  
  if (nrow(df) == 0) {
    return(tibble(
      n = 0, prevalence = NA_real_,
      TP = NA_integer_, FP = NA_integer_, TN = NA_integer_, FN = NA_integer_,
      accuracy = NA_real_, bal_accuracy = NA_real_, sens = NA_real_, spec = NA_real_,
      precision = NA_real_, recall = NA_real_, f_meas = NA_real_, kap = NA_real_
    ))
  }
  
  cm <- table(df$truth, df$pred)
  for (r in c("0","1")) for (c in c("0","1")) if (is.na(cm[r,c])) cm[r,c] <- 0
  TN <- cm["0","0"]; FP <- cm["0","1"]; FN <- cm["1","0"]; TP <- cm["1","1"]
  
  m <- metric_set(accuracy, bal_accuracy, sens, spec, precision, recall, f_meas, kap)
  met <- m(df, truth = truth, estimate = pred) %>%
    select(.metric, .estimate) %>%
    pivot_wider(names_from = .metric, values_from = .estimate)
  
  met %>%
    mutate(
      n = nrow(df),
      prevalence = mean(as.numeric(as.character(df$truth)) == 1),
      TP = as.integer(TP), FP = as.integer(FP),
      TN = as.integer(TN), FN = as.integer(FN)
    ) %>%
    select(n, prevalence, TP, FP, TN, FN, everything())
}

overall_summaries <- function(acc_tbl) {
  metric_cols <- intersect(
    c("accuracy", "bal_accuracy", "sens", "spec", "precision", "recall", "f_meas", "kap"),
    names(acc_tbl)
  )
  macro <- acc_tbl %>%
    summarise(across(all_of(metric_cols), ~ mean(.x, na.rm = TRUE))) %>%
    mutate(summary_type = "macro_mean") %>%
    select(summary_type, everything())
  
  weighted <- acc_tbl %>%
    summarise(across(all_of(metric_cols), ~ weighted.mean(.x, w = n, na.rm = TRUE))) %>%
    mutate(summary_type = "weighted_by_n") %>%
    select(summary_type, everything())
  
  bind_rows(macro, weighted)
}

# Filter points for a given raster row (year match + exclusion rules)
filter_points_for_raster <- function(pts_sf, raster_year, excl_sources, excl_ids = NULL, id_col = NULL) {
  out <- pts_sf %>%
    filter(year == raster_year)
  
  if (length(excl_sources) > 0) {
    out <- out %>% filter(!(.data[[source_col]] %in% excl_sources))
  }
  
  if (!is.null(excl_ids) && length(excl_ids) > 0 && !is.null(id_col) && id_col %in% names(out)) {
    # coerce to character for safe matching
    out <- out %>% filter(!(as.character(.data[[id_col]]) %in% as.character(excl_ids)))
  }
  
  out
}

# -----------------------------
# 3) READ RASTER SPECS + LOAD RASTERS
# -----------------------------
specs <- read_csv(raster_specs_csv, show_col_types = FALSE)

required_cols <- c("raster_id", "file", "year")
missing_cols <- setdiff(required_cols, names(specs))
if (length(missing_cols) > 0) stop("Raster specs CSV missing columns: ", paste(missing_cols, collapse = ", "))

specs <- specs %>%
  mutate(
    year = as.integer(year),
    exclude_sources = if ("exclude_sources" %in% names(.)) exclude_sources else "",
    exclude_point_ids = if ("exclude_point_ids" %in% names(.)) exclude_point_ids else ""
  )

if (any(!file.exists(specs$file))) {
  bad <- specs$file[!file.exists(specs$file)]
  stop("Some raster files do not exist:\n", paste(bad, collapse = "\n"))
}

# Load rasters into a named list by raster_id (NOT by year)
r_list <- set_names(lapply(specs$file, terra::rast), specs$raster_id)

# Assume same CRS across rasters; take CRS from first raster
ref_r <- r_list[[1]]
r_crs <- terra::crs(ref_r)

# -----------------------------
# 4) READ POINTS + DERIVE YEAR
# -----------------------------
pts <- st_read(points_shp, quiet = TRUE)

if (!truth_col %in% names(pts)) stop("Truth column not found: ", truth_col)
if (!date_col  %in% names(pts)) stop("Date column not found: ", date_col)
if (!source_col %in% names(pts)) stop("Source column not found (required for filtering): ", source_col)

pts <- pts %>%
  mutate(
    .obs_dt = parse_obs_datetime(.data[[date_col]]),
    year    = lubridate::year(.obs_dt),
    truth   = as.integer(.data[[truth_col]]),
    .source = as.character(.data[[source_col]])
  )

if (any(is.na(pts$year))) stop("Some points have unparseable dates in: ", date_col)
if (any(!pts$truth %in% c(0,1), na.rm = TRUE)) stop("Truth column must be 0/1: ", truth_col)

# Project points to raster CRS
pts <- st_transform(pts, crs = r_crs)

# -----------------------------
# 5) METHOD 1: STRICT POINT EXTRACTION (per raster)
# -----------------------------
extract_strict_for_raster <- function(raster_id) {
  r <- r_list[[raster_id]]
  sp <- specs %>% filter(raster_id == !!raster_id) %>% slice(1)
  
  raster_year <- sp$year
  excl_sources <- split_semicolon(sp$exclude_sources)
  excl_ids <- split_semicolon(sp$exclude_point_ids)
  
  pts_use <- filter_points_for_raster(
    pts_sf = pts,
    raster_year = raster_year,
    excl_sources = excl_sources,
    excl_ids = excl_ids,
    id_col = id_col
  )
  
  if (nrow(pts_use) == 0) {
    return(list(
      matched = tibble(),
      metrics = tibble(raster_id = raster_id, year = raster_year, method = "strict_point",
                       n = 0, prevalence = NA_real_, TP = NA_integer_, FP = NA_integer_,
                       TN = NA_integer_, FN = NA_integer_)
    ))
  }
  
  ex <- terra::extract(r, terra::vect(pts_use))
  pred <- as.integer(ex[, 2])
  
  matched <- pts_use %>%
    st_drop_geometry() %>%
    mutate(
      raster_id = raster_id,
      raster_year = raster_year,
      method = "strict_point",
      pred = pred
    ) %>%
    select(method, raster_id, raster_year, truth, pred, everything())
  
  metrics <- matched %>%
    calc_metrics(truth_col = "truth", pred_col = "pred") %>%
    mutate(method = "strict_point", raster_id = raster_id, year = raster_year) %>%
    select(method, raster_id, year, everything())
  
  list(matched = matched, metrics = metrics)
}

res_strict <- map(specs$raster_id, extract_strict_for_raster)
matched_strict <- bind_rows(map(res_strict, "matched"))
metrics_strict <- bind_rows(map(res_strict, "metrics"))

write_csv(matched_strict, file.path(out_dir, "matched_points_strict_by_raster.csv"))
write_csv(metrics_strict, file.path(out_dir, "accuracy_strict_by_raster.csv"))

overall_strict <- overall_summaries(metrics_strict) %>%
  mutate(method = "strict_point") %>%
  select(method, everything())
write_csv(overall_strict, file.path(out_dir, "overall_summary_strict.csv"))

# Optional: assessment accuracy by SOURCE (after training sources removed)
if (!is.null(source_col) && nrow(matched_strict) > 0) {
  metrics_strict_by_source <- matched_strict %>%
    group_by(raster_id, raster_year, .data[[source_col]]) %>%
    group_modify(~ calc_metrics(.x, truth_col = "truth", pred_col = "pred")) %>%
    ungroup() %>%
    mutate(method = "strict_point") %>%
    rename(year = raster_year) %>%
    select(method, raster_id, year, source = .data[[source_col]], everything())
  
  write_csv(metrics_strict_by_source, file.path(out_dir, "accuracy_strict_by_raster_by_source.csv"))
}

# -----------------------------
# 6) METHOD 2: BUFFERED EXTRACTION + MAJORITY VOTE (per raster, per buffer)
# -----------------------------
extract_buffer_for_raster <- function(raster_id, buffer_m) {
  r <- r_list[[raster_id]]
  sp <- specs %>% filter(raster_id == !!raster_id) %>% slice(1)
  
  raster_year <- sp$year
  excl_sources <- split_semicolon(sp$exclude_sources)
  excl_ids <- split_semicolon(sp$exclude_point_ids)
  
  pts_use <- filter_points_for_raster(
    pts_sf = pts,
    raster_year = raster_year,
    excl_sources = excl_sources,
    excl_ids = excl_ids,
    id_col = id_col
  )
  
  if (nrow(pts_use) == 0) {
    return(list(
      matched = tibble(),
      metrics = tibble(raster_id = raster_id, year = raster_year,
                       method = paste0("buffer_majority_", buffer_m, "m"),
                       n = 0, prevalence = NA_real_, TP = NA_integer_, FP = NA_integer_,
                       TN = NA_integer_, FN = NA_integer_)
    ))
  }
  
  # Buffer points (assumes projected CRS in meters; if geographic, reproject to a metric CRS before this)
  buf <- st_buffer(pts_use, dist = buffer_m)
  
  ex <- terra::extract(r, terra::vect(buf))  # multiple rows per polygon: ID, value
  pred_by_id <- tapply(ex[, 2], ex[, 1], majority_vote)
  pred <- as.integer(pred_by_id[as.character(seq_len(nrow(pts_use)))])
  
  matched <- pts_use %>%
    st_drop_geometry() %>%
    mutate(
      raster_id = raster_id,
      raster_year = raster_year,
      method = paste0("buffer_majority_", buffer_m, "m"),
      pred = pred
    ) %>%
    select(method, raster_id, raster_year, truth, pred, everything())
  
  metrics <- matched %>%
    calc_metrics(truth_col = "truth", pred_col = "pred") %>%
    mutate(method = paste0("buffer_majority_", buffer_m, "m"),
           raster_id = raster_id,
           year = raster_year) %>%
    select(method, raster_id, year, everything())
  
  list(matched = matched, metrics = metrics)
}

res_buffers <- list()
for (buffer_m in buffer_m_values) {
  res <- map(specs$raster_id, extract_buffer_for_raster, buffer_m = buffer_m)
  res_buffers[[as.character(buffer_m)]] <- res
  
  matched_buf <- bind_rows(map(res, "matched"))
  metrics_buf <- bind_rows(map(res, "metrics"))
  
  write_csv(matched_buf, file.path(out_dir, paste0("matched_points_buffer_majority_", buffer_m, "m_by_raster.csv")))
  write_csv(metrics_buf, file.path(out_dir, paste0("accuracy_buffer_majority_", buffer_m, "m_by_raster.csv")))
  
  overall_buf <- overall_summaries(metrics_buf) %>%
    mutate(method = paste0("buffer_majority_", buffer_m, "m")) %>%
    select(method, everything())
  
  write_csv(overall_buf, file.path(out_dir, paste0("overall_summary_buffer_majority_", buffer_m, "m.csv")))
  
  # Optional by-source for buffered
  if (!is.null(source_col) && nrow(matched_buf) > 0) {
    metrics_buf_by_source <- matched_buf %>%
      group_by(raster_id, raster_year, .data[[source_col]]) %>%
      group_modify(~ calc_metrics(.x, truth_col = "truth", pred_col = "pred")) %>%
      ungroup() %>%
      mutate(method = paste0("buffer_majority_", buffer_m, "m")) %>%
      rename(year = raster_year) %>%
      select(method, raster_id, year, source = .data[[source_col]], everything())
    
    write_csv(metrics_buf_by_source,
              file.path(out_dir, paste0("accuracy_buffer_majority_", buffer_m, "m_by_raster_by_source.csv")))
  }
}

# -----------------------------
# 7) METHOD 3: COMBINED OUTPUTS ACROSS ALL METHODS
# -----------------------------
metrics_all <- metrics_strict
for (buffer_m in buffer_m_values) {
  mpath <- file.path(out_dir, paste0("accuracy_buffer_majority_", buffer_m, "m_by_raster.csv"))
  metrics_all <- bind_rows(metrics_all, read_csv(mpath, show_col_types = FALSE))
}

write_csv(metrics_all, file.path(out_dir, "accuracy_all_methods_by_raster.csv"))

message("\n=== DONE ===")
message("Outputs written to: ", normalizePath(out_dir))
message("Key outputs:")
message(" - accuracy_strict_by_raster.csv")
message(" - accuracy_buffer_majority_*m_by_raster.csv")
message(" - overall_summary_*.csv")