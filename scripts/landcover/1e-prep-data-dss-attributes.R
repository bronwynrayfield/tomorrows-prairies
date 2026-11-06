# ============================================================
# Purpose: Extract soil covariates from DSS v3 shapefiles
#          and rasterize to the modelling extent template
#
# Inputs:
#   - DSS v3 shapefiles for AB, SK, MB:
#       dss_v3_XX.shp          (polygons / PAT)
#       dss_v3_XX_cmp.dbf      (component table)
#       soil_name_XX_v2.dbf    (soil name table)
#       soil_layer_XX_v2.dbf   (soil layer table)
#   - Modelling extent raster template (me_r)
#
# Outputs (all in processeddatasets/dss_filepath/):
#   - dss_drainage_modelling-extent.tif   INT1U 1-7
#   - dss_slope_modelling-extent.tif      FLT4S percent
#   - dss_watertbl_modelling-extent.tif   INT1U 0-4
#   - dss_kp0_modelling-extent.tif        INT2S percent vol, layer 1
#   - dss_kp33_modelling-extent.tif       INT2S percent vol, layer 1
#   - dss_kp1500_modelling-extent.tif     INT2S percent vol, layer 1
#   - dss_ldepth_modelling-extent.tif     INT2S cm, layer 1 (lower horizon depth)
#   - dss_ldepth_max_modelling-extent.tif INT2S cm, all layers (total soil depth)
#
# Notes:
#   - DRAINAGE ordinal:
#       VR = 1  Very rapidly drained
#       R  = 2  Rapidly drained
#       W  = 3  Well drained
#       MW = 4  Moderately well drained
#       I  = 5  Imperfectly drained
#       P  = 6  Poorly drained
#       VP = 7  Very poorly drained
#   - WATERTBL ordinal:
#       NO = 0  Never present
#       YU = 1  Present during unspecified period
#       YN = 2  Present during non-growing season
#       YG = 3  Present during growing season
#       YB = 4  Always present
#   - SLOPE_P: -9 treated as NA
#   - KP10 excluded: no valid layer 1 values for SK in source data
#   - UDEPTH excluded: layer 1 UDEPTH is always 0;
#     redundant with LDEPTH at all other layers
#   - All variables: percent-weighted mean across CMP components
#     per polygon, weighted by PERCENT; -9 sentinel values -> NA
#   - SLT variables (KP*, LDEPTH layer 1): surface horizon only
#   - Provinces rasterized individually then mosaicked
#
#   - PROVINCIAL SEAM. Provinces digitized the boundary separately,
#     leaving a sliver no polygon covers; the cell-centre rule leaves
#     those cells NA in the mosaic. Filled in two passes:
#       (a) touches = TRUE, limited to cells touched by two provinces:
#           real surveyed attributes, reaches ~1 cell in from each side
#       (b) focal() on the interior remainder: modal (ordinal) or mean,
#           i.e. INTERPOLATED, not surveyed
#     Genuine no-data (water) touches one province only, so stays NA,
#     though its edge erodes 1 cell where it falls in the seam mask.
#   - seam_cells = 2: sets mask width and focal window. I also tried 3 and 
#     checked output in QGIS. 2 was sufficient to fix slivers
#
# ============================================================

source(file.path("scripts/landcover", "0a-setup.R"))

# --- Parameters -------------------------------------------------------

dss_root <- file.path("existingdatasets", "Soil/detailed-soil-survey")
out_dir  <- file.path("processeddatasets", dss_filepath)

provinces <- c("ab", "sk", "mb")

drain_lookup    <- c(VR = 1L, R = 2L, W = 3L, MW = 4L, I = 5L, P = 6L, VP = 7L)
watertbl_lookup <- c(NO = 0L, YU = 1L, YN = 2L, YG = 3L, YB = 4L)

me_r <- rast(file.path(sa_filepath, me_raster_filename))

# --- Helper -----------------------------------------------------------

load_province <- function(prov, field) {
  
  pat <- st_read(
    file.path(dss_root, prov, paste0("dss_v3_", prov, ".shp")),
    quiet = TRUE
  )
  cmp <- read.dbf(
    file.path(dss_root, prov, paste0("dss_v3_", prov, "_cmp.dbf")),
    as.is = TRUE
  )
  names(cmp) <- toupper(names(cmp))
  cmp <- cmp |> mutate(PERCENT = as.numeric(PERCENT))
  
  # SNT-sourced fields
  if (field %in% c("drain_wtmean", "watertbl_wtmean")) {
    snt <- read.dbf(
      file.path(dss_root, prov, paste0("soil_name_", prov, "_v2.dbf")),
      as.is = TRUE
    )
    names(snt) <- toupper(names(snt))
    
    if (field == "drain_wtmean") {
      cmp <- cmp |>
        left_join(snt |> select(SOIL_ID, DRAINAGE), by = "SOIL_ID") |>
        mutate(drain_wtmean = drain_lookup[DRAINAGE])
      
      poly_val <- cmp |>
        filter(!is.na(drain_wtmean), !is.na(PERCENT), PERCENT > 0) |>
        group_by(POLY_ID) |>
        summarise(
          drain_wtmean = round(
            sum(drain_wtmean * PERCENT) / sum(PERCENT)),
          .groups = "drop")
    }
    
    if (field == "watertbl_wtmean") {
      cmp <- cmp |>
        left_join(snt |> select(SOIL_ID, WATERTBL), by = "SOIL_ID") |>
        mutate(watertbl_wtmean = watertbl_lookup[WATERTBL])
      
      poly_val <- cmp |>
        filter(!is.na(watertbl_wtmean), !is.na(PERCENT), PERCENT > 0) |>
        group_by(POLY_ID) |>
        summarise(
          watertbl_wtmean = round(
            sum(watertbl_wtmean * PERCENT) / sum(PERCENT)),
          .groups = "drop")
    }
  }
  
  # CMP-sourced field
  if (field == "slope_wtmean") {
    cmp <- cmp |>
      mutate(slope_wtmean = ifelse(as.numeric(SLOPE_P) == -9, NA_real_,
                                   as.numeric(SLOPE_P)))
    poly_val <- cmp |>
      filter(!is.na(slope_wtmean), !is.na(PERCENT), PERCENT > 0) |>
      group_by(POLY_ID) |>
      summarise(
        slope_wtmean = sum(slope_wtmean * PERCENT) / sum(PERCENT),
        .groups = "drop")
  }
  
  # SLT-sourced fields
  if (field %in% c("kp0_wtmean", "kp10_wtmean", "kp33_wtmean",
                   "kp1500_wtmean", "ldepth_wtmean")) {
    slt <- read.dbf(
      file.path(dss_root, prov, paste0("soil_layer_", prov, "_v2.dbf")),
      as.is = TRUE
    )
    names(slt) <- toupper(names(slt))
    
    src_col <- toupper(sub("_wtmean", "", field))  # e.g. "kp0_wtmean" -> "KP0"
    # LDEPTH: -9 is a valid depth (organic LFH layer above mineral surface) 
    slt1 <- slt |>
      filter(as.integer(LAYER_NO) == 1) |>
      mutate(val = ifelse(as.numeric(.data[[src_col]]) == -9 & src_col != "LDEPTH", 
                          NA_real_,
                          as.numeric(.data[[src_col]]))) |>
      select(SOIL_ID, val)
    
    cmp <- cmp |> left_join(slt1, by = "SOIL_ID")
    
    poly_val <- cmp |>
      filter(!is.na(val), !is.na(PERCENT), PERCENT > 0) |>
      group_by(POLY_ID) |>
      summarise(
        !!field := sum(val * PERCENT) / sum(PERCENT),
        .groups = "drop")
  }
  
  if (field == "ldepth_max") {
    slt <- read.dbf(
      file.path(dss_root, prov, paste0("soil_layer_", prov, "_v2.dbf")),
      as.is = TRUE
    )
    names(slt) <- toupper(names(slt))
    
    slt_depth <- slt |>
      mutate(LDEPTH = ifelse(as.numeric(LDEPTH) == -9, NA_real_, as.numeric(LDEPTH))) |>
      group_by(SOIL_ID) |>
      summarise(ldepth_max = max(LDEPTH, na.rm = TRUE), .groups = "drop") |>
      mutate(ldepth_max = ifelse(is.infinite(ldepth_max), NA_real_, ldepth_max))
    
    cmp <- cmp |> left_join(slt_depth, by = "SOIL_ID")
    
    poly_val <- cmp |>
      filter(!is.na(ldepth_max), !is.na(PERCENT), PERCENT > 0) |>
      group_by(POLY_ID) |>
      summarise(
        ldepth_max = sum(ldepth_max * PERCENT) / sum(PERCENT),
        .groups = "drop")
  }
  
  pat |>
    left_join(poly_val, by = "POLY_ID") |>
    select(POLY_ID, all_of(field)) |>
    st_transform(crs(me_r))
}

# --- Variables --------------------------------------------------------

vars <- list(
   list(field = "drain_wtmean",    file = dss_drainage_filename,   dt = "INT1U"),
   list(field = "slope_wtmean",    file = dss_slope_filename,      dt = "FLT4S"),
   list(field = "watertbl_wtmean", file = dss_watertbl_filename,   dt = "INT1U"),
   list(field = "kp0_wtmean",      file = dss_kp0_filename,        dt = "INT2S"),
   list(field = "kp33_wtmean",     file = dss_kp33_filename,       dt = "INT2S"),
   list(field = "kp1500_wtmean",   file = dss_kp1500_filename,     dt = "INT2S"),
   list(field = "ldepth_wtmean",   file = dss_ldepth_filename,     dt = "INT2S"),
   list(field = "ldepth_max",      file = dss_ldepth_max_filename, dt = "INT2S")
)

# --- Seam mask --------------------------------------------------------
# Built from province geometry only, so it is computed once and is
# identical across every DSS variable.

prov_geom <- lapply(provinces, function(prov)
  vect(st_transform(
    st_read(file.path(dss_root, prov, paste0("dss_v3_", prov, ".shp")), quiet = TRUE),
    crs(me_r))))
names(prov_geom) <- provinces

seam_cells <- 2L
seam_buf   <- seam_cells * res(me_r)[1]

# touches = TRUE reaches only ~1 cell, so buffer the dissolved province
# outlines by the measured gap width to make the mask span it. aggregate()
# first: buffering ~100k individual polygons is prohibitively slow, and the
# mask needs province identity only, not attributes.
prov_touch <- lapply(prov_geom, function(p)
  rasterize(buffer(aggregate(p), width = seam_buf), me_r, field = 1, touches = TRUE))

# Multiply gives 1 where both provinces reach, NA otherwise (NA propagates).
# Kept as separate seams so the focal pass can run cropped to each: the two
# boundaries sit at opposite ends of the study area, so their combined
# extent is effectively the full raster.
seams    <- list(prov_touch$ab * prov_touch$sk,
                 prov_touch$sk * prov_touch$mb)
seam     <- cover(seams[[1]], seams[[2]])
seam_ext <- lapply(seams, function(s) ext(trim(s)))

# Focal window: after the touches pass the residual gap is seam_cells - 2
# wide, and its centre cell sits ceiling(k/2) cells from the nearest valid
# data. Smallest window that reaches it - a larger one would pull values
# from soils well off the boundary.
k <- max(1L, seam_cells - 2L)
w <- 2L * ceiling(k / 2) + 1L

# --- Rasterize --------------------------------------------------------

for (v in vars) {
  message(v$field)
  
  pols <- lapply(provinces, function(prov)
    crop(vect(load_province(prov, v$field)), me_r))
  
  # Primary pass: cell-centre rule
  r <- do.call(mosaic, c(lapply(pols, function(p)
    rasterize(p, me_r, field = v$field, fun = "max")), fun = "first"))
  
  # Pass (a): surveyed attributes into seam cells the polygons overlap.
  # cover() writes into NA only, so nothing already valid is altered.
  patch <- do.call(mosaic, c(lapply(pols, function(p)
    rasterize(p, me_r, field = v$field, fun = "max", touches = TRUE)), fun = "first"))
  r <- cover(r, mask(patch, seam))
  
  # Pass (b): interpolate the interior cells no polygon reaches
  fill_fun <- if (v$dt == "INT1U") "modal" else "mean"  # ordinal classes: no averaging
  for (se in seam_ext) {
    rf <- focal(crop(r, se), w = w, fun = fill_fun, na.policy = "only", na.rm = TRUE)
    r  <- cover(r, mask(extend(rf, r), seam))
  }
  
  r <- mask(r, me_r)
  
  write_rast_tagged(
    r,
    root      = "./processeddatasets",
    rel_dir   = dss_filepath,
    filename = sub("\\.tif$", "_modelling-extent.tif", v$file),
    datatype  = v$dt,
    overwrite = TRUE,
    gdal      = c("COMPRESS=DEFLATE", "PREDICTOR=2")
  )
  
  rm(pols, r, patch)
  gc()
}
