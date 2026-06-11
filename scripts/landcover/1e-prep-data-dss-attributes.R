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
#   - dss_drainage_wtmean.tif   INT1U 1-7
#   - dss_slope_wtmean.tif      FLT4S percent
#   - dss_watertbl_dom.tif      INT1U 0-4
#   - dss_kp0_wtmean.tif        INT2S percent vol, layer 1
#   - dss_kp10_wtmean.tif       INT2S percent vol, layer 1
#   - dss_kp33_wtmean.tif       INT2S percent vol, layer 1
#   - dss_kp1500_wtmean.tif     INT2S percent vol, layer 1
#   - dss_udepth_wtmean.tif     INT2S cm, layer 1
#   - dss_ldepth_wtmean.tif     INT2S cm, layer 1
#
# Notes:
#   - DRAINAGE ordinal: VR=1, R=2, W=3, MW=4, I=5, P=6, VP=7
#   - WATERTBL ordinal: NO=0, YU=1, YN=2, YG=3, YB=4
#   - SLOPE_P: -9 treated as NA
#   - SLT variables (KP*, UDEPTH, LDEPTH): layer 1 only;
#     -9 sentinel values treated as NA
#   - All CMP-level variables: percent-weighted mean per polygon
#     across components, weighted by PERCENT
#   - WATERTBL: dominant component (highest PERCENT) per polygon
#   - Provinces rasterized individually then mosaicked
#
# Dependencies: terra, sf, dplyr, foreign
# ============================================================

source(file.path("scripts/landcover", "0-setup.R"))

library(foreign)

# --- Parameters -------------------------------------------------------

dss_root <- file.path("existingdatasets", "Soil/detailed-soil-survey")
out_dir  <- file.path("processeddatasets", dss_filepath)

provinces <- c("ab", "sk", "mb")

drain_lookup   <- c(VR = 1L, R = 2L, W = 3L, MW = 4L, I = 5L, P = 6L, VP = 7L)
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
                   "kp1500_wtmean", "udepth_wtmean", "ldepth_wtmean")) {
    slt <- read.dbf(
      file.path(dss_root, prov, paste0("soil_layer_", prov, "_v2.dbf")),
      as.is = TRUE
    )
    names(slt) <- toupper(names(slt))
    
    src_col <- toupper(sub("_wtmean", "", field))  # e.g. "kp0_wtmean" -> "KP0"
    
    slt1 <- slt |>
      filter(as.integer(LAYER_NO) == 1) |>
      mutate(val = ifelse(as.numeric(.data[[src_col]]) == -9, NA_real_,
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

# --- Rasterize --------------------------------------------------------

vars <- list(
#  list(field = "drain_wtmean",  file = "dss_drainage_wtmean.tif",  dt = "INT1U"),
#  list(field = "slope_wtmean",  file = "dss_slope_wtmean.tif",     dt = "FLT4S"),
#  list(field = "watertbl_wtmean",  file = "dss_watertbl_wtmean.tif",     dt = "INT1U"),
#  list(field = "kp0_wtmean",    file = "dss_kp0_wtmean.tif",       dt = "INT2S"),
##  list(field = "kp10_wtmean",   file = "dss_kp10_wtmean.tif",      dt = "INT2S"),
#  list(field = "kp33_wtmean",   file = "dss_kp33_wtmean.tif",      dt = "INT2S"),
#  list(field = "kp1500_wtmean", file = "dss_kp1500_wtmean.tif",    dt = "INT2S"),
#  list(field = "udepth_lyr1_wtmean", file = "dss_udepth_wtmean.tif",    dt = "INT2S"),
#  list(field = "ldepth_lyr1_wtmean", file = "dss_ldepth_wtmean.tif",    dt = "INT2S"),
  list(field = "ldepth_max", file = "dss_ldepth_max.tif", dt = "INT2S")
)


for (v in vars) {
  message(v$field)
  
  tiles <- lapply(provinces, function(prov) {
    pol <- load_province(prov, field = v$field)
    vec <- crop(vect(pol), me_r)
    rasterize(vec, me_r, field = v$field, fun = "max")
  })
  
  r <- do.call(mosaic, c(tiles, fun = "first"))
  r <- mask(r, me_r)
  
  writeRaster(
    r,
    file.path(out_dir, v$file),
    datatype  = v$dt,
    overwrite = TRUE,
    gdal      = c("COMPRESS=DEFLATE", "PREDICTOR=2")
  )
  
  rm(tiles, r)
  gc()
}