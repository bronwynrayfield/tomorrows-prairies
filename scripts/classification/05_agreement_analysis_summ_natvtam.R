
###################################################################################
##           ANALYSIS: SPATIAL AGREEMENT SUMMARIES (NATIVE VS TAME)              ##
###################################################################################

##----------------------- Alison Long - January 2026 -------------------------

# step 5: run descriptive data summaries to assess similarities and differences
#        among existing grassland datasets

# load packages
library(terra)
library(sf)
library(tidyverse)

# Set the temporary directory to a drive with sufficient space
terraOptions(tempdir = "R:/CanadianGrasslands/WG_classification/scratch")

####--- DIRECTORIES ---####

# pilot area
pa_dir <- "R:/CanadianGrasslands/processeddatasets/pilotarea"

# study area
sa_dir <- "R:/CanadianGrasslands/processeddatasets/studyarea"

# grasslands
proc_grass_data_dir <- "R:/CanadianGrasslands/WG_Classification/processed"

# grass raster directories
aci_aligned_dir  <- file.path(proc_grass_data_dir, "aafc_aci/aligned")
ngi_aligned_dir  <- file.path(proc_grass_data_dir, "cfga_ngi/aligned")
cgr_masked_dir   <- file.path(proc_grass_data_dir, "cgr_gam_adv/masked")
eccc_masked_dir  <- file.path(proc_grass_data_dir, "eccc_grass/masked")
cnvc_masked_dir  <- file.path(proc_grass_data_dir, "natureserve_cnvc/masked")
pudl_aligned_dir <- file.path(proc_grass_data_dir, "phjv_pudl/aligned")
pli_aligned_dir  <- file.path(proc_grass_data_dir, "sk_pli/aligned")
gpw_masked_dir   <- file.path(proc_grass_data_dir, "wri_gpw_ann/masked")
plow_masked_dir  <- file.path(proc_grass_data_dir, "wwf_plowprint/masked")
plvi_aligned_dir <- file.path(proc_grass_data_dir, "ab_plvi/aligned")
hfi_aligned_dir  <- file.path(proc_grass_data_dir, "abmi_hfi/aligned")

# output 
out_dir <- "R:/CanadianGrasslands/WG_Classification/results/spatial_agreement_analysis"

####--- FILES ---####

# Study area
sa_file <- file.path(sa_dir, "class_study_area_aligned_cnvc_Dec2025.tif")

# ecodistricts
ecodist_file   <- "R:/CanadianGrasslands/processeddatasets/boundaries/national-ecological-framework-for-canada/Ecodistricts/ecodistricts_clean.shp"  

# grass raster files
aci_file  <- file.path(aci_aligned_dir,  "aafc_aci_merged_reproj_algn_natvtam_noshrub.tif")
ngi_file  <- file.path(ngi_aligned_dir,  "cfga_ngi_merged_reproj_algn_natvtam_good.tif")
cgr_file  <- file.path(cgr_masked_dir,   "cgr_reproj_algn_natvtam_masked.tif")
eccc_file <- file.path(eccc_masked_dir,  "eccc_reclass_reproj_algn_natvtam_masked.tif")
cnvc_file <- file.path(cnvc_masked_dir,  "cnvc_algn_natvtam_masked.tif")
pudl_file <- file.path(pudl_aligned_dir, "pudl_reclass_reproj_algn_natvtam_noshrub.tif")
pli_file  <- file.path(pli_aligned_dir,  "sk_pli_merged_reproj_algn_natvtam_noshrub.tif")
gpw_file  <- file.path(gpw_masked_dir,   "gpw_reproj_algn_natvtam_noshrub_masked.tif")
plow_file <- file.path(plow_masked_dir,  "??.tif")
plvi_file <- file.path(plvi_aligned_dir, "plvi_reproj_algn_natvtam.tif")
hfi_file  <- file.path(hfi_aligned_dir, "hfi_reproj_algn_binary.tif") # tame only - same as binary

####--- READ INTO TERRA ---####

# ecodistricts
ecodist <- vect(ecodist_file)

# study area
sa <- rast(sa_file)

# grass rasters
aci  <- rast(aci_file)
ngi  <- rast(ngi_file)
cgr  <- rast(cgr_file)
eccc <- rast(eccc_file)
cnvc <- rast(cnvc_file)
pudl <- rast(pudl_file)
pli  <- rast(pli_file)
gpw  <- rast(gpw_file)
plow <- rast(plow_file)
plvi <- rast(plvi_file)
hfi  <- rast(hfi_file)

files <- c(
  aci_file, ngi_file, cgr_file, eccc_file, cnvc_file,  pudl_file, pli_file,
  gpw_file, plow_file, plvi_file, hfi_file
)

# Load as a SpatRaster stack
rstack <- rast(files)

#-------- PIXELWISE AGREEMENT ------------

# Count how many rasters have value == 1 for each pixel (ignores NA)
count_tame <- terra::app(rstack, fun = function(x) sum(x == 1, na.rm = TRUE))

count_tame
plot(count_tame)

# Save result 
writeRaster(count_tame, file.path(out_dir, "pixel_data_count_tame.tif"), overwrite = TRUE)

# Count how many rasters have value == 2 for each pixel (ignores NA)
count_native <- terra::app(rstack, fun = function(x) sum(x == 2, na.rm = TRUE))

count_native
plot(count_native)

# Save result 
writeRaster(count_native, file.path(out_dir, "pixel_data_count_native.tif"), overwrite = TRUE)

## For each agreement class, ha/ecodistrict

#-------- TOTAL GRASSLAND AREA -----------------
# Sum grassland area (ha) for a single raster r inside study area sa
# r: SpatRaster with values {1=grass, 0=not grass, NA=no data}
# sa: SpatRaster with values {1=inside study area, NA=outside}

# Named list of rasters (so names carry into the output)
r_list <- list(
  aafc_aci = aci,
  ab_plvi = plvi,
  abmi_hfi = hfi,
  cfga_ngi = ngi,
  cgr_gam_adv = cgr,
  eccc_grass = eccc,
  natureserve_cnvc = cnvc,
  phjv_pudl = pudl,
  sk_pli = pli,
  wri_gpw_ann = gpw,
  wwf_plowprint = plow
)

# Precompute cell areas once on the shared grid
cs <- cellSize(r_list[[1]], unit = "m")

native_hectares <- function(r, sa, cs) {
  # Safety: ensure sa is on same grid as r (should be if truly same geometry)
  if (!compareGeom(r, sa, stopOnError = FALSE)) {
    stop("sa is not on the same grid as r; please project/align first.")
  }
  # Keep area only where r==2 and inside study area (sa non-NA)
  cs_native_inside <- ifel((r == 2) & !is.na(sa), cs, NA)
  
  # Sum and convert to hectares
  total_m2 <- global(cs_native_inside, "sum", na.rm = TRUE)[1, 1]
  as.numeric(total_m2) / 10000
}

# Compute total hectares for each raster
native_ha_results <- purrr::map_dbl(r_list, ~ native_hectares(.x, sa, cs))
native_ha_results

# Tidy tibble with results
native_ha_df <- tibble::tibble(
  raster = names(native_ha_results),
  native_ha = unname(native_ha_results)
)

# Print to console
native_ha_df

# write to CSV in your outputs directory
readr::write_csv(native_ha_df, file.path(out_dir, "native_total_hectares_per_raster_5prov.csv"))

tame_hectares <- function(r, sa, cs) {
  # Safety: ensure sa is on same grid as r (should be if truly same geometry)
  if (!compareGeom(r, sa, stopOnError = FALSE)) {
    stop("sa is not on the same grid as r; please project/align first.")
  }
  # Keep area only where r==2 and inside study area (sa non-NA)
  cs_tame_inside <- ifel((r == 2) & !is.na(sa), cs, NA)
  
  # Sum and convert to hectares
  total_m2 <- global(cs_tame_inside, "sum", na.rm = TRUE)[1, 1]
  as.numeric(total_m2) / 10000
}

# Compute total hectares for each raster
tame_ha_results <- purrr::map_dbl(r_list, ~ tame_hectares(.x, sa, cs))
tame_ha_results

# Tidy tibble with results
tame_ha_df <- tibble::tibble(
  raster = names(tame_ha_results),
  tame_ha = unname(tame_ha_results)
)

# Print to console
tame_ha_df

# write to CSV in your outputs directory
readr::write_csv(tame_ha_df, file.path(out_dir, "tame_total_hectares_per_raster_5prov.csv"))

#-------- GRASSLAND AREA PER ECODISTRICT -----------------
# - All rasters in r_list share the SAME grid and CRS as each other and as sa.
# - sa: 1 inside, NA outside
# - each r: 1 = grassland, 0 = not grassland, NA = nodata
# - ecodist: SpatVector with a key column named ECODISTRIC (ecodistrict ID)

# 0) Validate ecodistrict ID column
if (!"ECODISTRIC" %in% names(ecodist)) {
  stop("The ecodistrict SpatVector must contain a column named 'ECODISTRIC'.")
}

# 1) Sanity check grid alignment
template <- r_list[[1]]
ok <- compareGeom(template, sa, stopOnError = FALSE) &&
  all(vapply(r_list, function(x) compareGeom(x, template, stopOnError = FALSE), logical(1)))
if (!ok) stop("Not all rasters and 'sa' share the exact same grid/CRS.")

# 2) Precompute per-cell area ONCE (m^2) and inside-study-area mask
cs     <- cellSize(template, unit = "m")  # m^2 per cell, on the shared grid
inside <- !is.na(sa)                      # TRUE inside study area, FALSE/NA outside

# 3) Create a numeric zone_id for rasterization (robust even if ECODISTRIC is character)
ecodist$zone_id <- seq_len(nrow(ecodist))

# 4) Rasterize ecodistricts ONCE to create a zones raster matching the shared grid
#    touches = FALSE => cell belongs to a zone if its center falls inside the polygon
zones <- rasterize(
  ecodist, template, field = "zone_id",
  touches = FALSE, background = NA
  # filename = "ecodist_zones.tif", overwrite = TRUE  # optionally persist to disk
)

# Keep a small attribute table for joins later
ecodist_df <- as.data.frame(ecodist, geom = FALSE) %>%
  select(zone_id, ECODISTRIC)

# 5) Helper: per-raster zonal sum of grassland area (ha) by ecodistrict
zonal_grass_ha <- function(r) {
  # Keep cell area only where r==1 and inside study area; NA elsewhere
  cs_grass <- ifel((r == 1) & inside, cs, NA)
  
  # Sum area per ecodistrict (m^2); returns only zones with at least one valid cell
  z <- zonal(cs_grass, zones, fun = "sum", na.rm = TRUE)  # columns: zone (zone_id), sum
  
  # Convert to tibble, to hectares
  z <- as_tibble(z) %>%
    transmute(zone_id = .data$zone, grass_ha = .data$sum / 10000)
  
  # Ensure all ecodistricts appear (fill missing with 0)
  full <- tibble(zone_id = ecodist_df$zone_id)
  full %>%
    left_join(z, by = "zone_id") %>%
    mutate(grass_ha = replace_na(grass_ha, 0)) %>%
    left_join(ecodist_df, by = "zone_id") %>%
    relocate(zone_id, ECODISTRIC)
}

# 6) Compute per-ecodistrict grassland hectares for each raster in r_list (named)
by_raster <- imap(r_list, ~{
  out <- zonal_grass_ha(.x)
  out$raster <- .y
  out
})

# 7) Long and wide outputs

# Long tidy table: one row per (raster, ecodistrict)
result_long <- bind_rows(by_raster) %>%
  relocate(raster, zone_id, ECODISTRIC, grass_ha)

# Wide table: one row per ecodistrict, one column per raster (hectares)
result_wide <- result_long %>%
  select(zone_id, ECODISTRIC, raster, grass_ha) %>%
  pivot_wider(names_from = raster, values_from = grass_ha)

# Inspect results
result_long
result_wide

# Write to disk
readr::write_csv(result_long, "grassland_ha_by_ecodistrict_long_5prov.csv")
readr::write_csv(result_wide, "grassland_ha_by_ecodistrict_wide_5prov.csv")

# To save polygon layer with the attribute joined:
#ed_out <- ed
#ed_out$hectares_with_data <- out$hectares_with_data[match(ed_out$ecod_id, out$ecod_id)]
#writeVector(ed_out, "ecodistricts_with_hectares.gpkg", filetype = "GPKG", overwrite = TRUE)

#################################################################################
#################################################################################

# Version optimized for large rasters:

# ================================
# Large-Raster Optimized Workflow
# ================================

# ---- Packages ----
library(terra)
library(sf)
library(tidyverse)

# ---- Terra options: stream to disk, use large temp, show progress ----
terraOptions(
  tempdir  = "R:/CanadianGrasslands/WG_classification/scratch",  # your scratch
  todisk   = TRUE,        # encourage disk-based processing
  memfrac  = 0.8,         # proportion of RAM terra may use (hint)
  progress = 1            # show progress bars
)

# ==================
# DIRECTORIES
# ==================

# pilot area
pa_dir <- "R:/CanadianGrasslands/processeddatasets/pilotarea"

# study area
sa_dir <- "R:/CanadianGrasslands/processeddatasets/studyarea"

# grasslands
proc_grass_data_dir <- "R:/CanadianGrasslands/WG_Classification/processed"

# grass raster directories
aci_aligned_dir  <- file.path(proc_grass_data_dir, "aafc_aci/aligned")
ngi_aligned_dir  <- file.path(proc_grass_data_dir, "cfga_ngi/aligned")
cgr_masked_dir   <- file.path(proc_grass_data_dir, "cgr_gam_adv/masked")
eccc_masked_dir  <- file.path(proc_grass_data_dir, "eccc_grass/masked")
cnvc_masked_dir  <- file.path(proc_grass_data_dir, "natureserve_cnvc/masked")
pudl_aligned_dir <- file.path(proc_grass_data_dir, "phjv_pudl/aligned")
pli_aligned_dir  <- file.path(proc_grass_data_dir, "sk_pli/aligned")
gpw_masked_dir   <- file.path(proc_grass_data_dir, "wri_gpw_ann/masked")
plow_masked_dir  <- file.path(proc_grass_data_dir, "wwf_plowprint/masked")
plvi_aligned_dir <- file.path(proc_grass_data_dir, "ab_plvi/aligned")
hfi_aligned_dir  <- file.path(proc_grass_data_dir, "abmi_hfi/aligned")

# output 
out_dir <- "R:/CanadianGrasslands/WG_Classification/results/spatial_agreement_analysis"
#if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ==================
# FILES
# ==================

# ecodistricts
ecodist_file <- "R:/CanadianGrasslands/processeddatasets/boundaries/national-ecological-framework-for-canada/Ecodistricts/ecodistricts_clean.shp"

# grass raster files (aligned rasters where 0==not grass, 1==tame, 2==native)
aci_file  <- file.path(aci_aligned_dir,  "aafc_aci_merged_reproj_algn_natvtam_noshrub.tif")
ngi_file  <- file.path(ngi_aligned_dir,  "cfga_ngi_merged_reproj_algn_natvtam_good.tif")
cgr_file  <- file.path(cgr_masked_dir,   "cgr_reproj_algn_natvtam_masked.tif")
eccc_file <- file.path(eccc_masked_dir,  "eccc_reclass_reproj_algn_natvtam_masked.tif")
cnvc_file <- file.path(cnvc_masked_dir,  "cnvc_algn_natvtam_masked.tif")
pudl_file <- file.path(pudl_aligned_dir, "pudl_reclass_reproj_algn_natvtam_noshrub.tif")
pli_file  <- file.path(pli_aligned_dir,  "sk_pli_merged_reproj_algn_natvtam_noshrub.tif")
gpw_file  <- file.path(gpw_masked_dir,   "gpw_reproj_algn_natvtam_noshrub_masked.tif")
plow_file <- file.path(plow_masked_dir,  "??.tif")
plvi_file <- file.path(plvi_aligned_dir, "plvi_reproj_algn_natvtam.tif")
hfi_file  <- file.path(hfi_aligned_dir,  "??.tif")

files <- c(
  aci_file, ngi_file, cgr_file, eccc_file, cnvc_file, pudl_file, pli_file,
  gpw_file, plow_file, plvi_file, hfi_file
)

# ==================
# READ DATA (lazy links, not full read)
# ==================

# vector
ecodist <- vect(ecodist_file)

# raster stack (SpatRaster of layers on disk)
rstack <- rast(files)

# OPTIONAL: Quick alignment check (same extent/res/crs/dim)
# This is inexpensive compared to analysis and helps catch a misaligned file.
message("Checking alignment across rasters...")
for (i in 2:nlyr(rstack)) {
  ok <- compareGeom(rstack[[1]], rstack[[i]], crs = TRUE, ext = TRUE, rowcol = TRUE, res = TRUE, stopOnError = FALSE)
  if (!ok) stop(sprintf("Alignment mismatch between layer 1 and layer %d: %s", i, names(rstack)[i]))
}
message("Alignment OK.")

# Name the layers (nice for outputs)
names(rstack) <- c(
  "aafc_aci","cfga_ngi","cgr_gam_adv","eccc_grass","natureserve_cnvc","phjv_pudl","sk_pli",
  "wri_gpw_ann","wwf_plowprint","ab_plvi","abmi_hfi"
)

# Build a named list as in your original code (pointing to the same layers)
#r_list <- lapply(seq_len(nlyr(rstack)), function(i) rstack[[i]])
#names(r_list) <- names(rstack)

# Split lists by role - HFI has tame class only so need to exclude from native list
r_list_native <- as.list(rstack[[setdiff(names(rstack), "abmi_hfi")]])
names(r_list_native) <- setdiff(names(rstack), "abmi_hfi")

r_list_tame <- as.list(rstack)
names(r_list_tame) <- names(rstack)

# ==============================
# 1) TOTAL NATIVE GRASSLAND AREA (ha)
# ==============================

# Compute per-cell area ONCE on the shared grid; persist (reused many times)
message("Computing per-cell area once and reusing...")
cs_m2 <- cellSize(rstack[[1]], unit = "m")

# (Optional) write to disk so it's cached across sessions and reusable
cs_m2 <- writeRaster(
  cs_m2,
  filename = file.path(out_dir, "cellsize_m2.tif"),
  wopt     = list(datatype = "FLT4S", gdal = c("COMPRESS=LZW")),
  overwrite = TRUE
)

# Efficient hectare calculation using ifel() + global() sum
native_hectares <- function(r, cs_precomp) {
  # Area where grass == 2, else 0 (no mask raster created)
  native_area_m2 <- ifel(r == 2, cs_precomp, 0)
  total_m2 <- global(native_area_m2, "sum", na.rm = TRUE)[1, 1]
  as.numeric(total_m2) / 10000  # m2 -> ha
}

message("Summarizing native grassland hectares per raster...")
native_ha_results <- purrr::map_dbl(r_list_native, ~ native_hectares(.x, cs_m2))

native_ha_tbl <- tibble::tibble(
  raster = names(native_ha_results),
  native_hectares = as.numeric(native_ha_results)
) |>
  dplyr::arrange(raster)

print(native_ha_tbl)

readr::write_csv(native_ha_tbl, file.path(out_dir, "native_total_hectares_per_raster.csv"))

# ==============================
# 2) TOTAL TAME GRASSLAND AREA (ha)
# ==============================

# Efficient hectare calculation using ifel() + global() sum
tame_hectares <- function(r, cs_precomp) {
  # Area where grass == 1, else 0 (no mask raster created)
  tame_area_m2 <- ifel(r == 1, cs_precomp, 0)
  total_m2 <- global(tame_area_m2, "sum", na.rm = TRUE)[1, 1]
  as.numeric(total_m2) / 10000  # m2 -> ha
}

message("Summarizing tame grassland hectares per raster...")
tame_ha_results <- purrr::map_dbl(r_list_tame, ~ tame_hectares(.x, cs_m2))

tame_ha_tbl <- tibble::tibble(
  raster = names(tame_ha_results),
  tame_hectares = as.numeric(tame_ha_results)
) |>
  dplyr::arrange(raster)

print(tame_ha_tbl)

readr::write_csv(tame_ha_tbl, file.path(out_dir, "tame_total_hectares_per_raster.csv"))

# ==========================================
# 3) NATIVE GRASSLAND AREA (ha) PER ECODISTRICT
# ==========================================

# Project ecodistricts ONCE to raster CRS/grid
message("Projecting ecodistricts once to raster CRS...")
ecodist_proj <- if (!terra::same.crs(rstack, ecodist)) terra::project(ecodist, rstack) else ecodist

# Prepare attribute lookup (keep only needed fields)
eco_attrs <- as.data.frame(ecodist_proj, geom = FALSE)
keep_fields <- intersect(
  c("ECODISTRIC", "ECOZONE", "ECOPROVIN", "NAME", "ENGLISH_NA"),
  names(eco_attrs)
)
eco_lookup <- tibble::as_tibble(eco_attrs[, keep_fields, drop = FALSE]) |>
  dplyr::distinct(ECODISTRIC, .keep_all = TRUE)

# Function: grass hectares per ecodistrict (reuses projected polygons + cell areas)
native_ha_by_ecodist <- function(r, eco_proj, cs_precomp, id_field = "ECODISTRIC") {
  native_area_m2 <- ifel(r == 2, cs_precomp, 0)
  
  ex <- terra::extract(
    native_area_m2, eco_proj,
    fun = "sum", na.rm = TRUE, weights = TRUE
  )
  
  names(ex)[names(ex) == "sum"] <- "native_m2"
  
  eco_vals <- as.data.frame(eco_proj, geom = FALSE)
  ex[[id_field]] <- eco_vals[[id_field]][ex$ID]
  
  out <- tibble::tibble(
    !!id_field := ex[[id_field]],
    native_hectares = as.numeric(ex$native_m2) / 10000
  )
  
  # If 0 preferred instead of NA for polygons with no overlap, uncomment:
  # out$native_hectares[is.na(out$native_hectares)] <- 0
  
  out
}

message("Extracting native grassland hectares by ecodistrict for each raster (areal-weighted)...")
eco_native_list <- purrr::imap(r_list_native, function(r, nm) {
  res <- native_ha_by_ecodist(r, ecodist_proj, cs_m2, id_field = "ECODISTRIC")
  res$raster <- nm
  res
})

eco_native_long <- dplyr::bind_rows(eco_native_list) |>
  dplyr::left_join(eco_lookup, by = "ECODISTRIC") |>
  dplyr::relocate(raster, ECODISTRIC) |>
  dplyr::arrange(raster, ECODISTRIC)

# Save long table
readr::write_csv(eco_native_long, file.path(out_dir, "native_hectares_by_ecodistrict_long.csv"))

# Also produce a wide table (one column per raster)
eco_native_wide <- eco_native_long |>
  dplyr::select(ECODISTRIC, dplyr::all_of(setdiff(keep_fields, "ECODISTRIC")), raster, native_hectares) |>
  tidyr::pivot_wider(
    names_from  = raster,
    values_from = native_hectares
  ) |>
  dplyr::arrange(ECODISTRIC)

# To fill missing with 0 in the wide table, uncomment:
# eco_native_wide <- eco_native_wide |>
#   dplyr::mutate(across(-c(ECODISTRIC, all_of(setdiff(keep_fields, "ECODISTRIC"))),
#                        ~ tidyr::replace_na(., 0)))

readr::write_csv(eco_native_wide, file.path(out_dir, "native_hectares_by_ecodistrict_wide.csv"))

# Optional: quick peek
# print(dplyr::glimpse(eco_native_long))

# ==========================================
# 4) TAME GRASSLAND AREA (ha) PER ECODISTRICT
# ==========================================

# Project ecodistricts ONCE to raster CRS/grid
message("Projecting ecodistricts once to raster CRS...")
ecodist_proj <- if (!terra::same.crs(rstack, ecodist)) terra::project(ecodist, rstack) else ecodist

# Prepare attribute lookup (keep only needed fields)
eco_attrs <- as.data.frame(ecodist_proj, geom = FALSE)
keep_fields <- intersect(
  c("ECODISTRIC", "ECOZONE", "ECOPROVIN", "NAME", "ENGLISH_NA"),
  names(eco_attrs)
)
eco_lookup <- tibble::as_tibble(eco_attrs[, keep_fields, drop = FALSE]) |>
  dplyr::distinct(ECODISTRIC, .keep_all = TRUE)

# Function: tame grass hectares per ecodistrict (reuses projected polygons + cell areas)
tame_ha_by_ecodist <- function(r, eco_proj, cs_precomp, id_field = "ECODISTRIC") {
  tame_area_m2 <- ifel(r == 1, cs_precomp, 0)
  
  ex <- terra::extract(
    tame_area_m2, eco_proj,
    fun = "sum", na.rm = TRUE, weights = TRUE
  )
  
  names(ex)[names(ex) == "sum"] <- "tame_m2"
  
  eco_vals <- as.data.frame(eco_proj, geom = FALSE)
  ex[[id_field]] <- eco_vals[[id_field]][ex$ID]
  
  out <- tibble::tibble(
    !!id_field := ex[[id_field]],
    tame_hectares = as.numeric(ex$tame_m2) / 10000
  )
  
  # If 0 preferred instead of NA for polygons with no overlap, uncomment:
  # out$tame_hectares[is.na(out$tame_hectares)] <- 0
  
  out
}

message("Extracting tame grassland hectares by ecodistrict for each raster (areal-weighted)...")
eco_tame_list <- purrr::imap(r_list_tame, function(r, nm) {
  res <- tame_ha_by_ecodist(r, ecodist_proj, cs_m2, id_field = "ECODISTRIC")
  res$raster <- nm
  res
})

eco_tame_long <- dplyr::bind_rows(eco_tame_list) |>
  dplyr::left_join(eco_lookup, by = "ECODISTRIC") |>
  dplyr::relocate(raster, ECODISTRIC) |>
  dplyr::arrange(raster, ECODISTRIC)

# Save long table
readr::write_csv(eco_tame_long, file.path(out_dir, "tame_hectares_by_ecodistrict_long.csv"))

# Also produce a wide table (one column per raster)
eco_tame_wide <- eco_tame_long |>
  dplyr::select(ECODISTRIC, dplyr::all_of(setdiff(keep_fields, "ECODISTRIC")), raster, tame_hectares) |>
  tidyr::pivot_wider(
    names_from  = raster,
    values_from = tame_hectares
  ) |>
  dplyr::arrange(ECODISTRIC)

# To fill missing with 0 in the wide table, uncomment:
# eco_tame_wide <- eco_tame_wide |>
#   dplyr::mutate(across(-c(ECODISTRIC, all_of(setdiff(keep_fields, "ECODISTRIC"))),
#                        ~ tidyr::replace_na(., 0)))

readr::write_csv(eco_tame_wide, file.path(out_dir, "tame_hectares_by_ecodistrict_wide.csv"))

# Optional: quick peek
# print(dplyr::glimpse(eco_tame_long))

message("All done.")
# End of script