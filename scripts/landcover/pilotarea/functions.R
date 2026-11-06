# Functions to write spatial data to disk with helpful tags (e.g. cropped)
# Generic tagging
collapse_tags <- function(...) {
  tags <- c(...)
  tags <- tags[nzchar(tags)]
  if (length(tags) == 0) NULL else paste(tags, collapse = "_")
}

tagged_path <- function(root, rel_dir, filename, ...) {
  stem <- tools::file_path_sans_ext(filename)
  ext  <- tools::file_ext(filename)
  tag  <- collapse_tags(...)
  name <- if (is.null(tag)) stem else sprintf("%s_%s", stem, tag)
  if (identical(ext, "")) ext <- "tif"
  file.path(root, rel_dir, sprintf("%s.%s", name, ext))
}
# Write vector data
write_sf_tagged <- function(x, root, rel_dir, filename, ..., quiet = TRUE) {
  stopifnot(inherits(x, "sf"))
  stem <- tools::file_path_sans_ext(filename)
  tag  <- collapse_tags(...)
  layer <- if (is.null(tag)) stem else sprintf("%s_%s", stem, tag)
  dsn <- file.path(root, rel_dir)
  dir.create(dsn, recursive = TRUE, showWarnings = FALSE)
  sf::st_write(
    x,
    dsn = dsn,
    layer = layer,
    driver = "ESRI Shapefile",
    delete_layer = TRUE,
    quiet = quiet
  )
}
# Write raster data
write_rast_tagged <- function(r, root, rel_dir, filename, ..., overwrite = TRUE,
                              gdal = c("COMPRESS=LZW", "TILED=YES")) {
  stopifnot(inherits(r, "SpatRaster"))
  out <- tagged_path(root, rel_dir, filename, ...)
  dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
  terra::writeRaster(r, out, overwrite = overwrite, gdal = gdal)
}




#--- Reusable function to summarize grassland loss area by destination class and year ----
# Assumptions per user request:
#   • class and year rasters are ALREADY reprojected to an equal‑area CRS
#   • rasters are ALREADY clipped to the study area (same CRS, res, extent)
#   • area is uniform per pixel, computed from raster resolution
# Optional: summarize within spatial strata (fast rasterized IDs or exact, area‑weighted)
#
# Inputs
#   class_tif  : path to *_loss_c05_grass_b01_class.tif (destination class codes)
#   year_tif   : path to *_loss_c05_grass_b03_year.tif  (encoded years: year - 1900)
#   strata     : OPTIONAL sf/SpatVector of polygons for strata (same CRS as rasters)
#   strata_id  : OPTIONAL name of ID field in strata (required if strata provided)
#   strata_exact : if TRUE use exactextractr (area‑weighted); else rasterize IDs (fast)
#   out_units : 'ha' (default) or 'km2'
#   class_lookup : OPTIONAL named vector mapping codes→labels; defaults provided below
#   year_min/year_max : OPTIONAL calendar year filters (after decoding)
#   drop_classes : OPTIONAL numeric vector of destination class codes to drop
#   keep_codes : include numeric codes in output (default TRUE)
#
# Output
#   Data frame with columns (no strata): year | to_class | class_to | pixels | area_ha (or area_km2)
#   With strata: stratum_id | year | to_class | class_to | pixels | area_ha (or area_km2)

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
})

.default_class_lookup <- function() {
  c(
    `1`  = "forest",
    `2`  = "shrub",
    `3`  = "crop",
    `4`  = "upast",    # pasture (non‑intensive)
    `5`  = "grass",
    `6`  = "builtup",
    `7`  = "bare",
    `8`  = "wetland",
    `9`  = "water",
    `10` = "ipast",    # pasture (intensive)
    `11` = "forest/shrub",
    `12` = "forest/wetland",
    `13` = "shrub/forest",
    `14` = "wetland/forest",
    `15` = "upast/grass",
    `16` = "grass/upast"
  )
}

summarize_loss <- function(
    r_class,
    r_year,
    strata = NULL,
    strata_id = NULL,
    strata_exact = FALSE,
    out_units = c("ha", "km2"),
    class_lookup = NULL,
    year_min = NULL,
    year_max = NULL,
    drop_classes = NULL,
    keep_codes = TRUE
) {
  out_units <- match.arg(out_units)
  if (is.null(class_lookup)) class_lookup <- .default_class_lookup()
  
  # Geometry check
  if (!compareGeom(r_class, r_year, stopOnError = FALSE)) {
    stop("class and year rasters do not share the same geometry (CRS/res/extent).")
  }
  
  # Pixel area (uniform), assuming equal‑area CRS already
  resm <- res(r_class)
  pix_area_ha <- (resm[1] * resm[2]) / 10000
  
  # Decode years (calendar); drop 0/NA (no change) rows later
  r_year_cal <- r_year + 1900
  
  # Build stack
  s <- c(r_class, r_year_cal)
  names(s) <- c("class_code_raw", "year")
  
  # Helper to split from/to codes (supports either combined or simple codes)
  decode_from_to <- function(x) {
    x <- as.integer(round(x))
    from_code <- ifelse(x >= 100, x %/% 100, NA_integer_)
    to_code   <- ifelse(x >= 100, x %% 100, x)
    data.frame(class_code_raw = x, from_code = from_code, to_code = to_code)
  }
  
  # Optional strata
  has_strata <- !is.null(strata)
  if (has_strata && is.null(strata_id)) stop("Provide strata_id when supplying strata.")
  
  if (!has_strata) {
    df <- as.data.frame(s, na.rm = TRUE)
    names(df) <- c("class_code_raw", "year")
    
    # Drop empty / invalid encodings: year < 1985 (i.e., 1900 from encoded 0), or class code 0
    df <- df |> filter(!is.na(class_code_raw), !is.na(year)) |> filter(year >= 1985)
    
    # Decode from/to
    dec <- decode_from_to(df$class_code_raw)
    df <- bind_cols(df, dec[, c("from_code", "to_code")])
    
    # Remove obviously empty codes (0) and optional drops
    df <- df |> filter(!is.na(to_code), to_code != 0)
    if (!is.null(drop_classes)) df <- df |> filter(!(to_code %in% drop_classes))
    
    # Year filters
    if (!is.null(year_min)) df <- df |> filter(year >= year_min)
    if (!is.null(year_max)) df <- df |> filter(year <= year_max)
    
    # Map labels
    df$from_class <- unname(class_lookup[as.character(df$from_code)])
    df$to_class   <- unname(class_lookup[as.character(df$to_code)])
    # Fallback labels for unexpected/missing codes
    df$from_class[is.na(df$from_class) & !is.na(df$from_code)] <- paste0("code_", df$from_code)
    df$to_class[is.na(df$to_class) & !is.na(df$to_code)]       <- paste0("code_", df$to_code)
    
    # Summarize by year × destination (and retain from_code for QA/optionality)
    sum_df <- df |>
      group_by(year, from_code, from_class, to_code, to_class) |>
      summarise(pixels = n(), area_ha = pixels * pix_area_ha, .groups = "drop")
    
  } else if (has_strata && !strata_exact) {
    # Fast: rasterize strata IDs to grid (assign by cell center)
    v_strata <- if (inherits(strata, "SpatVector")) strata else vect(strata)
    if (!identical(crs(v_strata), crs(r_class))) stop("strata CRS must match raster CRS (already equal‑area).")
    if (!(strata_id %in% names(v_strata))) stop("strata_id not found in strata.")
    
    r_strata <- rasterize(v_strata, r_class, field = strata_id)
    s2 <- c(s, r_strata)
    names(s2) <- c("class_code_raw", "year", "stratum_id")
    
    df <- as.data.frame(s2, na.rm = TRUE)
    df <- df |> filter(!is.na(class_code_raw), !is.na(year), !is.na(stratum_id)) |> filter(year >= 1985)
    
    dec <- decode_from_to(df$class_code_raw)
    df <- bind_cols(df, dec[, c("from_code", "to_code")])
    df <- df |> filter(!is.na(to_code), to_code != 0)
    
    if (!is.null(drop_classes)) df <- df |> filter(!(to_code %in% drop_classes))
    if (!is.null(year_min)) df <- df |> filter(year >= year_min)
    if (!is.null(year_max)) df <- df |> filter(year <= year_max)
    
    df$from_class <- unname(class_lookup[as.character(df$from_code)])
    df$to_class   <- unname(class_lookup[as.character(df$to_code)])
    df$from_class[is.na(df$from_class) & !is.na(df$from_code)] <- paste0("code_", df$from_code)
    df$to_class[is.na(df$to_class) & !is.na(df$to_code)]       <- paste0("code_", df$to_code)
    
    sum_df <- df |>
      group_by(stratum_id, year, from_code, from_class, to_code, to_class) |>
      summarise(pixels = n(), area_ha = pixels * pix_area_ha, .groups = "drop")
    
  } else {
    # Exact: area‑weighted using exactextractr
    if (!requireNamespace("exactextractr", quietly = TRUE)) {
      stop("strata_exact=TRUE requires the 'exactextractr' package.")
    }
    
    v_strata <- if (inherits(strata, "SpatVector")) strata else vect(strata)
    if (!identical(crs(v_strata), crs(r_class))) stop("strata CRS must match raster CRS (already equal‑area).")
    v_strata_df <- terra::as.data.frame(v_strata, geom = FALSE)
    if (!(strata_id %in% names(v_strata_df))) stop("strata_id not found in strata.")
    
    s3 <- s  # class_code_raw, year
    
    sum_list <- list()
    n_poly <- nrow(v_strata_df)
    for (i in seq_len(n_poly)) {
      poly_i <- v_strata[i]
      sid    <- v_strata_df[[strata_id]][i]
      
      ee <- exactextractr::exact_extract(s3, sf::st_as_sf(poly_i))[[1]]
      ee <- ee[!is.na(ee$class_code_raw) & !is.na(ee$year) & !is.na(ee$coverage_fraction), ]
      ee <- ee[ee$year >= 1985, ]
      if (!nrow(ee)) next
      
      # Decode from/to
      dec <- decode_from_to(ee$class_code_raw)
      ee$from_code <- dec$from_code
      ee$to_code   <- dec$to_code
      ee <- ee[!is.na(ee$to_code) & ee$to_code != 0, ]
      
      if (!is.null(drop_classes)) ee <- ee[!(ee$to_code %in% drop_classes), ]
      if (!is.null(year_min)) ee <- ee[ee$year >= year_min, ]
      if (!is.null(year_max)) ee <- ee[ee$year <= year_max, ]
      if (!nrow(ee)) next
      
      # Area per record = coverage_fraction * pixel_area
      ee$area_ha <- ee$coverage_fraction * pix_area_ha
      
      ee$from_class <- unname(class_lookup[as.character(ee$from_code)])
      ee$to_class   <- unname(class_lookup[as.character(ee$to_code)])
      ee$from_class[is.na(ee$from_class) & !is.na(ee$from_code)] <- paste0("code_", ee$from_code)
      ee$to_class[is.na(ee$to_class) & !is.na(ee$to_code)]       <- paste0("code_", ee$to_code)
      
      sum_i <- ee |>
        group_by(year, from_code, from_class, to_code, to_class) |>
        summarise(pixels = n(), area_ha = sum(area_ha, na.rm = TRUE), .groups = "drop") |>
        mutate(stratum_id = sid)
      
      sum_list[[length(sum_list) + 1]] <- sum_i
    }
    
    sum_df <- dplyr::bind_rows(sum_list)
    if (nrow(sum_df)) sum_df <- sum_df |> select(stratum_id, year, from_code, from_class, to_code, to_class, pixels, area_ha)
  }
  
  # Units
  if (out_units == "km2") {
    sum_df <- mutate(sum_df, area_km2 = area_ha / 100) |> select(-area_ha)
  }
  
  # Optionally hide numeric codes
  if (!keep_codes) {
    sum_df <- sum_df |> select(-matches("_code$"))
  }
  
  # Order
  ord_cols <- intersect(c("stratum_id", "year", "from_class", "to_class"), names(sum_df))
  sum_df <- arrange(sum_df, across(all_of(ord_cols)))
  
  as.data.frame(sum_df)
}

# --- Minimal example (commented) ---
# res <- summarize_grass_loss(
#   class_tif = ".../recentchange_2022_loss_c05_grass_b01_class.tif",
#   year_tif  = ".../recentchange_2022_loss_c05_grass_b03_year.tif",
#   out_units = "ha"
# )
# head(res)

# --- Helper: quick QA to list unexpected 'from' codes (should mostly be 5 or 16 for grass) ---
# qa_from_codes <- function(df) dplyr::count(df, from_code, from_class, sort = TRUE) (commented) ---
# res <- summarize_grass_loss(
#   class_tif = ".../recentchange_2022_loss_c05_grass_b01_class.tif",
#   year_tif  = ".../recentchange_2022_loss_c05_grass_b03_year.tif",
#   out_units = "ha"
# )
# head(res)

#--- Reusable function to build a categorical SpatRaster of destination group for a given decade ----
make_decade_layer <- function(r_class, r_year, start, end_exclusive) {
  # year mask for the decade (e.g., 1990–1999 with end_exclusive=2000)
  r_year_cal <- r_year + 1900
  decade_mask <- r_year_cal >= start & r_year_cal < end_exclusive
  
  # destination class code = to_code
  to_code <- r_class %% 100
  
  # map to_code -> destination group (cropland, builtup, past, else NA)
  # 3=crop, 6=builtup, 4=upast, 10=ipast
  to_group <- classify(
    to_code,
    rcl = matrix(
      c(
        3, 1,   # cropland → 1
        6, 2,   # builtup  → 2
        4, 3,   # upast    → 3 (past)
        10, 3   # ipast    → 3 (past)
      ),
      ncol = 2, byrow = TRUE
    ),
    others = NA
  )
  
  # keep only pixels that changed in the decade
  out <- mask(to_group, decade_mask, maskvalues = 0, updatevalue = NA)
  
  # attach labels so legend shows names
  levels(out) <- data.frame(ID = c(1,2,3),
                            to_group = c("cropland","builtup","past"))
  
  out
}
