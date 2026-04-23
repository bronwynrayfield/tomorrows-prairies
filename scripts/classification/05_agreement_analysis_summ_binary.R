
###################################################################################
##               ANALYSIS: SPATIAL AGREEMENT SUMMARIES (BINARY)                  ##
###################################################################################

##----------------------- Alison Long - January 2026 -------------------------

# step 5: run descriptive data summaries to assess similarities and differences
#        among existing grassland datasets

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

# write options
wopt_prob  <- list(datatype="FLT4S", gdal=c("COMPRESS=LZW","TILED=YES","BIGTIFF=YES"))
wopt_byte  <- list(datatype="INT1U", gdal=c("COMPRESS=LZW","TILED=YES","BIGTIFF=YES"))
wopt_float <- list(datatype="FLT4S", gdal=c("COMPRESS=LZW","TILED=YES","BIGTIFF=YES"))
wopt_count <- list(datatype="INT2U", gdal=c("COMPRESS=LZW","TILED=YES","BIGTIFF=YES"))

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
aci_aligned_dir   <- file.path(proc_grass_data_dir, "aafc_aci/aligned")
luts30_masked_dir <- file.path(proc_grass_data_dir, "aafc_luts_2020_30m_v5/masked")
luts10_masked_dir <- file.path(proc_grass_data_dir, "aafc_luts_2020_10m_v1/masked")
cec_masked_dir    <- file.path(proc_grass_data_dir, "cec_nalcms/masked")
ngi_aligned_dir   <- file.path(proc_grass_data_dir, "cfga_ngi/aligned")
cgr_masked_dir    <- file.path(proc_grass_data_dir, "cgr_gam_adv/masked")
eccc_masked_dir   <- file.path(proc_grass_data_dir, "eccc_grass/masked")
esa_masked_dir    <- file.path(proc_grass_data_dir, "esa_worldcover/masked")
glc_masked_dir    <- file.path(proc_grass_data_dir, "glc_fcs30d/masked")
cnvc_masked_dir   <- file.path(proc_grass_data_dir, "natureserve_cnvc/masked")
lcc_masked_dir    <- file.path(proc_grass_data_dir, "nrcan_lcc/masked")
pudl_aligned_dir  <- file.path(proc_grass_data_dir, "phjv_pudl/aligned")
pli_aligned_dir   <- file.path(proc_grass_data_dir, "sk_pli/aligned")
lcr_masked_dir    <- file.path(proc_grass_data_dir, "statcan_lcr/masked")
gpw_masked_dir    <- file.path(proc_grass_data_dir, "wri_gpw_ann/masked")
plow_masked_dir   <- file.path(proc_grass_data_dir, "wwf_plowprint/masked")
plvi_aligned_dir  <- file.path(proc_grass_data_dir, "ab_plvi/aligned")
hfi_aligned_dir   <- file.path(proc_grass_data_dir, "abmi_hfi/aligned")
gcc_aligned_dir   <- file.path(proc_grass_data_dir, "bc_gcc/aligned")
dw_masked_dir     <- file.path(proc_grass_data_dir, "dw_2025_grass/masked")
esri_masked_dir   <- file.path(proc_grass_data_dir, "esri_rangelands/masked")

# output 
out_dir <- "R:/CanadianGrasslands/WG_Classification/results/spatial_agreement_analysis"
#if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ==================
# FILES
# ==================

# ecodistricts
ecodist_file <- "R:/CanadianGrasslands/processeddatasets/boundaries/national-ecological-framework-for-canada/Ecodistricts/ecodistricts_clean.shp"

# Study area (modellng extent + Peace River region)
sa_file <- file.path(sa_dir, "classification_study_area_aligned_cnvc_Apr2026.tif")

# grass raster files (binary 0/1 aligned)
aci_file    <- file.path(aci_aligned_dir, "aafc_aci_merged_reproj_algn_binary_noshrub.tif")
luts30_file <- file.path(luts30_masked_dir, "luts30_reproj_algn_binary_masked.tif")
luts10_file <- file.path(luts10_masked_dir, "luts10_reproj_algn_binary_masked.tif")
cec_file    <- file.path(cec_masked_dir, "nalcms_reproj_algn_binary_noshrub_masked.tif")
ngi_file    <- file.path(ngi_aligned_dir, "cfga_ngi_merged_reproj_algn_binary.tif")
cgr_file    <- file.path(cgr_masked_dir, "cgr_reproj_algn_binary_masked.tif")
eccc_file   <- file.path(eccc_masked_dir, "eccc_reclass_reproj_algn_binary_masked.tif")
esa_file    <- file.path(esa_masked_dir, "esa_reproj_algn_binary_noshrub_masked.tif")
glc_file    <- file.path(glc_masked_dir, "glc_reproj_algn_binary_noshrub_masked.tif")
cnvc_file   <- file.path(cnvc_masked_dir, "cnvc_algn_binary_masked.tif")
lcc_file    <- file.path(lcc_masked_dir, "lcc_reproj_algn_binary_noshrub_masked.tif")
pudl_file   <- file.path(pudl_aligned_dir, "pudl_reclass_reproj_algn_binary_noshrub.tif")
pli_file    <- file.path(pli_aligned_dir, "sk_pli_merged_reproj_algn_binary_noshrub.tif")
lcr_file    <- file.path(lcr_masked_dir, "lcr_reproj_algn_binary_masked.tif")
gpw_file    <- file.path(gpw_masked_dir, "gpw_reproj_algn_binary_noshrub_masked.tif")
plow_file   <- file.path(plow_masked_dir, "plow_reproj_algn_intact_masked.tif")
plvi_file   <- file.path(plvi_aligned_dir, "plvi_reproj_algn_binary_gte2015.tif")
hfi_file    <- file.path(hfi_aligned_dir, "hfi_reproj_algn_binary_gte2015.tif")
gcc_file    <- file.path(gcc_aligned_dir, "gcc_reproj_algn_binary.tif")
dw_file     <- file.path(dw_masked_dir, "dw_reproj_algn_masked_binary.tif")
esri_file   <- file.path(esri_masked_dir, "esri_reproj_algn_masked_binary.tif")

files <- c(
  aci_file, luts30_file, luts10_file, cec_file, ngi_file, cgr_file, eccc_file,
  esa_file, glc_file, cnvc_file, lcc_file, pudl_file, pli_file, lcr_file, 
  gpw_file, plow_file, plvi_file, hfi_file, gcc_file, dw_file, esri_file
)

# ==================
# READ DATA (lazy links, not full read)
# ==================

# ecodistricts
ecodist <- vect(ecodist_file)

# classification study area (modelling extent + Peace River region)
sa <- rast(sa_file)

# raster stack (SpatRaster of layers on disk)
rstack <- rast(files)

# Quick alignment check (same extent/res/crs/dim)
message("Checking alignment across rasters...")
for (i in 2:nlyr(rstack)) {
  ok <- compareGeom(rstack[[1]], rstack[[i]], crs = TRUE, ext = TRUE, rowcol = TRUE, res = TRUE, stopOnError = FALSE)
  if (!ok) stop(sprintf("Alignment mismatch between layer 1 and layer %d: %s", i, names(rstack)[i]))
}
message("Alignment OK.")

# Mask to study area (modelling extent + Peace River region)
mask_file <- file.path(out_dir, "binary_stack_masked.tif")

if (!file.exists(mask_file)) {
  
  sa <- rast(sa_file)
  
  r_mask <- mask(
    crop(rast(files), sa),
    sa,
    filename  = mask_file,
    overwrite = TRUE,
    wopt      = wopt_byte
  )
  
  #rm(sa, r_mask)
  gc()
}

r_mask <- rast(mask_file)

# Name the layers (nice for outputs)
names(r_mask) <- c(
  "aafc_aci","aafc_luts30m","aafc_luts10m","cec_nalcms","cfga_ngi","cgr_gam_adv","eccc_grass",
  "esa_worldcover","glc_fcs30d","natureserve_cnvc","nrcan_lcc","phjv_pudl","sk_pli","statcan_lcr",
  "wri_gpw_ann","wwf_plowprint","ab_plvi","abmi_hfi","bc_gcc","google_dw","esri_lulc"
)

# Build a named list as in your original code (pointing to the same layers)
r_list <- lapply(seq_len(nlyr(r_mask)), function(i) r_mask[[i]])
names(r_list) <- names(r_mask)

# ==============================
# 1) PIXELWISE COUNT OF NON-NA
# ==============================

# Use logical stack + sum; stream to disk; compact datatype + compression
message("Computing pixelwise non-NA count...")
has_data <- !is.na(r_mask)

count_nonNA <- app(
  has_data,
  fun      = sum,
  na.rm    = TRUE,
  filename = file.path(out_dir, "pixel_data_count_binary.tif"),
  wopt     = list(datatype = "INT2U", gdal = c("COMPRESS=LZW")),
  overwrite = TRUE
)

plot(count_nonNA)

# ==============================
# 2) TOTAL GRASSLAND AREA (ha)
# ==============================

# Compute per-cell area ONCE on the shared grid; persist (reused many times)
message("Computing per-cell area once and reusing...")
cs_m2 <- cellSize(r_mask[[1]], unit = "m")

# (Optional) write to disk so it's cached across sessions and reusable
cs_m2 <- writeRaster(
  cs_m2,
  filename = file.path(out_dir, "cellsize_m2.tif"),
  wopt     = list(datatype = "FLT4S", gdal = c("COMPRESS=LZW")),
  overwrite = TRUE
)

# Efficient hectare calculation using ifel() + global() sum
grass_hectares <- function(r, cs_precomp) {
  # Area where grass == 1, else 0 (no mask raster created)
  grass_area_m2 <- ifel(r == 1, cs_precomp, 0)
  total_m2 <- global(grass_area_m2, "sum", na.rm = TRUE)[1, 1]
  as.numeric(total_m2) / 10000  # m2 -> ha
}

message("Summarizing total grassland hectares per raster...")
grass_ha_results <- purrr::map_dbl(r_list, ~ grass_hectares(.x, cs_m2))

grass_ha_tbl <- tibble::tibble(
  raster = names(grass_ha_results),
  grass_hectares = as.numeric(grass_ha_results)
) |>
  dplyr::arrange(raster)

print(grass_ha_tbl)

readr::write_csv(grass_ha_tbl, file.path(out_dir, "grass_total_hectares_per_raster.csv"))

# ==========================================
# 3) GRASSLAND AREA (ha) PER ECODISTRICT
# ==========================================

# Project ecodistricts ONCE to raster CRS/grid
message("Projecting ecodistricts once to raster CRS...")
ecodist_proj <- if (!terra::same.crs(r_mask, ecodist)) terra::project(ecodist, r_mask) else ecodist

# Prepare attribute lookup (keep only needed fields)
eco_attrs <- as.data.frame(ecodist_proj)
keep_fields <- intersect(
  c("ECODISTRIC", "ECOZONE_ID", "ECOREGION_", "ECOPROVIN"),
  names(eco_attrs)
)
eco_lookup <- tibble::as_tibble(eco_attrs[, keep_fields, drop = FALSE]) |>
  dplyr::distinct(ECODISTRIC, .keep_all = TRUE)

# Function: grass hectares per ecodistrict (reuses projected polygons + cell areas)
grass_ha_by_ecodist <- function(r, eco_proj, cs_precomp, id_field = "ECODISTRIC") {
  grass_area_m2 <- ifel(r == 1, cs_precomp, 0)
  
  ex <- terra::extract(
    grass_area_m2, eco_proj,
    fun = "sum", na.rm = TRUE, weights = TRUE
  )
  
  names(ex)[names(ex) == "sum"] <- "grass_m2"
  
  eco_vals <- as.data.frame(eco_proj, geom = FALSE)
  ex[[id_field]] <- eco_vals[[id_field]][ex$ID]
  
  out <- tibble::tibble(
    !!id_field := ex[[id_field]],
    grass_hectares = as.numeric(ex$grass_m2) / 10000
  )
  
  # If 0 preferred instead of NA for polygons with no overlap, uncomment:
  # out$grass_hectares[is.na(out$grass_hectares)] <- 0
  
  out
}

message("Extracting grassland hectares by ecodistrict for each raster (areal-weighted)...")
eco_grass_list <- purrr::imap(r_list, function(r, nm) {
  res <- grass_ha_by_ecodist(r, ecodist_proj, cs_m2, id_field = "ECODISTRIC")
  res$raster <- nm
  res
})

eco_grass_long <- dplyr::bind_rows(eco_grass_list) |>
  dplyr::left_join(eco_lookup, by = "ECODISTRIC") |>
  dplyr::relocate(raster, ECODISTRIC) |>
  dplyr::arrange(raster, ECODISTRIC)

# Save long table
readr::write_csv(eco_grass_long, file.path(out_dir, "grass_hectares_by_ecodistrict_long.csv"))

# Also produce a wide table (one column per raster)
eco_grass_wide <- eco_grass_long |>
  dplyr::select(ECODISTRIC, dplyr::all_of(setdiff(keep_fields, "ECODISTRIC")), raster, grass_hectares) |>
  tidyr::pivot_wider(
    names_from  = raster,
    values_from = grass_hectares
  ) |>
  dplyr::arrange(ECODISTRIC)

# To fill missing with 0 in the wide table, uncomment:
# eco_grass_wide <- eco_grass_wide |>
#   dplyr::mutate(across(-c(ECODISTRIC, all_of(setdiff(keep_fields, "ECODISTRIC"))),
#                        ~ tidyr::replace_na(., 0)))

readr::write_csv(eco_grass_wide, file.path(out_dir, "grass_hectares_by_ecodistrict_wide.csv"))

# Optional: quick peek
# print(dplyr::glimpse(eco_grass_long))

message("All done.")
# End of script

########################################################################################
########################################################################################

# Original version:

# load packages
library(terra)
library(sf)
library(tidyverse)

# Set the temporary directory to a drive with sufficient space
terraOptions(tempdir = "R:/CanadianGrasslands/WG_classification/scratch")

####--- DIRECTORIES ---####

# pilot area
#pa_dir <- "R:/CanadianGrasslands/processeddatasets/pilotarea"

# study area
sa_dir <- "R:/CanadianGrasslands/processeddatasets/studyarea"

# grasslands
proc_grass_data_dir <- "R:/CanadianGrasslands/WG_Classification/processed"

# grass raster directories
aci_aligned_dir   <- file.path(proc_grass_data_dir, "aafc_aci/aligned")
luts30_masked_dir <- file.path(proc_grass_data_dir, "aafc_luts_2020_30m_v5/masked")
luts10_masked_dir <- file.path(proc_grass_data_dir, "aafc_luts_2020_10m_v1/masked")
cec_masked_dir    <- file.path(proc_grass_data_dir, "cec_nalcms/masked")
ngi_aligned_dir   <- file.path(proc_grass_data_dir, "cfga_ngi/aligned")
cgr_masked_dir    <- file.path(proc_grass_data_dir, "cgr_gam_adv/masked")
eccc_masked_dir   <- file.path(proc_grass_data_dir, "eccc_grass/masked")
esa_masked_dir    <- file.path(proc_grass_data_dir, "esa_worldcover/masked")
glc_masked_dir    <- file.path(proc_grass_data_dir, "glc_fcs30d/masked")
cnvc_masked_dir   <- file.path(proc_grass_data_dir, "natureserve_cnvc/masked")
lcc_masked_dir    <- file.path(proc_grass_data_dir, "nrcan_lcc/masked")
pudl_aligned_dir  <- file.path(proc_grass_data_dir, "phjv_pudl/aligned")
pli_aligned_dir   <- file.path(proc_grass_data_dir, "sk_pli/aligned")
lcr_masked_dir    <- file.path(proc_grass_data_dir, "statcan_lcr/masked")
gpw_masked_dir    <- file.path(proc_grass_data_dir, "wri_gpw_ann/masked")
plow_masked_dir   <- file.path(proc_grass_data_dir, "wwf_plowprint/masked")
plvi_aligned_dir  <- file.path(proc_grass_data_dir, "ab_plvi/aligned")
hfi_aligned_dir   <- file.path(proc_grass_data_dir, "abmi_hfi/aligned")
gcc_aligned_dir   <- file.path(proc_grass_data_dir, "bc_gcc/aligned")
dw_masked_dir     <- file.path(proc_grass_data_dir, "dw_2025_grass/masked")
esri_masked_dir   <- file.path(proc_grass_data_dir, "esri_rangelands/masked")

# output 
out_dir <- "R:/CanadianGrasslands/WG_Classification/results/spatial_agreement_analysis"

####--- FILES ---####

# Study area
sa_file <- file.path(sa_dir, "class_study_area_aligned_cnvc_Dec2025.tif")

# ecodistricts
ecodist_file   <- "R:/CanadianGrasslands/processeddatasets/boundaries/national-ecological-framework-for-canada/Ecodistricts/ecodistricts_clean.shp"  

# grass raster files
aci_file    <- file.path(aci_aligned_dir, "aafc_aci_merged_reproj_algn_binary_noshrub.tif")
luts30_file <- file.path(luts30_masked_dir, "luts30_reproj_algn_binary_masked.tif")
luts10_file <- file.path(luts10_masked_dir, "luts10_reproj_algn_binary_masked.tif")
cec_file    <- file.path(cec_masked_dir, "nalcms_reproj_algn_binary_noshrub_masked.tif")
ngi_file    <- file.path(ngi_aligned_dir, "cfga_ngi_merged_reproj_algn_binary.tif")
cgr_file    <- file.path(cgr_masked_dir, "cgr_reproj_algn_binary_masked.tif")
eccc_file   <- file.path(eccc_masked_dir, "eccc_reclass_reproj_algn_binary_masked.tif")
esa_file    <- file.path(esa_masked_dir, "esa_reproj_algn_binary_noshrub_masked.tif")
glc_file    <- file.path(glc_masked_dir, "glc_reproj_algn_binary_noshrub_masked.tif")
cnvc_file   <- file.path(cnvc_masked_dir, "cnvc_algn_binary_masked.tif")
lcc_file    <- file.path(lcc_masked_dir, "lcc_reproj_algn_binary_noshrub_masked.tif")
pudl_file   <- file.path(pudl_aligned_dir, "pudl_reclass_reproj_algn_binary_noshrub.tif")
pli_file <- file.path(pli_aligned_dir, "sk_pli_merged_reproj_algn_binary_noshrub.tif")
lcr_file    <- file.path(lcr_masked_dir, "lcr_reproj_algn_binary_masked.tif")
gpw_file    <- file.path(gpw_masked_dir, "gpw_reproj_algn_binary_noshrub_masked.tif")
plow_file   <- file.path(plow_masked_dir, "plow_reproj_algn_intact_masked.tif")
plvi_file   <- file.path(plvi_aligned_dir, "plvi_reproj_algn_binary_gte2015.tif")
hfi_file    <- file.path(hfi_aligned_dir, "hfi_reproj_algn_binary_gte2015.tif")
gcc_file    <- file.path(gcc_aligned_dir, "gcc_reproj_algn_binary.tif")
dw_file     <- file.path(dw_masked_dir, "dw_reproj_algn_masked_binary.tif")
esri_file   <- file.path(esri_masked_dir, "esri_reproj_algn_masked_binary.tif")

####--- READ INTO TERRA ---####

# ecodistricts
ecodist <- vect(ecodist_file)

# study area
sa <- rast(sa_file)

# grass rasters
aci <- rast(aci_file)
luts <- rast(luts_file)
cec <- rast(cec_file)
ngi <- rast(ngi_file)
cgr <- rast(cgr_file)
eccc <- rast(eccc_file)
esa <- rast(esa_file)
glc <- rast(glc_file)
cnvc <- rast(cnvc_file)
lcc <- rast(lcc_file)
pudl <- rast(pudl_file)
pli <- rast(pli_file)
lcr <- rast(lcr_file)
gpw <- rast(gpw_file)
plow <- rast(plow_file)

files <- c(
  aci_file, luts10_file, luts30_file, cec_file, ngi_file, cgr_file, eccc_file,
  esa_file, glc_file, cnvc_file, lcc_file, pudl_file, pli_file,
  lcr_file, gpw_file, plow_file,
  plvi_file, hfi_file, gcc_file, dw_file, esri_file
)

# Load as a SpatRaster stack
rstack <- rast(files)

#-------- PIXELWISE AGREEMENT ------------

# Count how many rasters have *non-NA* data for each pixel
count_nonNA <- app(rstack, fun = function(x) sum(!is.na(x)))

count_nonNA
plot(count_nonNA)

#class       : SpatRaster 
#size        : 86539, 130223, 1  (nrow, ncol, nlyr)
#resolution  : 30, 30  (x, y)
#extent      : -2321003, 1585687, 270327.4, 2866497  (xmin, xmax, ymin, ymax)
#coord. ref. : North_America_Albers_Equal_Area_Conic 
#source      : spat_344bc5ece498a_214204_o8up7yuDdJ6XHBX.tif 
#name        : lyr.1 
#min value   :     0 
#max value   :    15

# Save result 
writeRaster(count_nonNA, file.path(out_dir, "pixel_data_count_nonNA.tif"), overwrite = TRUE)

# Count how many rasters have value == 1 for each pixel (ignores NA)
count_ones <- terra::app(rstack, fun = function(x) sum(x == 1, na.rm = TRUE))

count_ones
plot(count_ones)

#class       : SpatRaster 
#size        : 86539, 130223, 1  (nrow, ncol, nlyr)
#resolution  : 30, 30  (x, y)
#extent      : -2321003, 1585687, 270327.4, 2866497  (xmin, xmax, ymin, ymax)
#coord. ref. : North_America_Albers_Equal_Area_Conic 
#source      : spat_344bc348c3169_214204_zxzxx9jyyx6UnnI.tif 
#name        : lyr.1 
#min value   :     0 
#max value   :    15

# Save result 
writeRaster(count_ones, file.path(out_dir, "pixel_data_count_binary.tif"), overwrite = TRUE)

## For each agreement class, ha/ecodistrict

#-------- TOTAL GRASSLAND AREA -----------------
# Sum grassland area (ha) for a single raster r inside study area sa
# r: SpatRaster with values {1=grass, 0=not grass, NA=no data}
# sa: SpatRaster with values {1=inside study area, NA=outside}

# Named list of rasters (so names carry into the output)
r_list <- list(
  aafc_aci = aci,
  aafc_luts = luts,
  cec_nalcms = cec,
  cfga_ngi = ngi,
  cgr_gam_adv = cgr,
  eccc_grass = eccc,
  esa_worldcover = esa,
  glc_fcs30d = glc,
  natureserve_cnvc = cnvc,
  nrcan_lcc = lcc,
  phjv_pudl = pudl,
  sk_pli = pli,
  statcan_lcr = lcr,
  wri_gpw_ann = gpw,
  wwf_plowprint = plow
)

# Precompute cell areas once on the shared grid
cs <- cellSize(r_list[[1]], unit = "m")

grass_hectares <- function(r, sa, cs) {
  # Safety: ensure sa is on same grid as r (should be if truly same geometry)
  if (!compareGeom(r, sa, stopOnError = FALSE)) {
    stop("sa is not on the same grid as r; please project/align first.")
  }
  # Keep area only where r==1 and inside study area (sa non-NA)
  cs_grass_inside <- ifel((r == 1) & !is.na(sa), cs, NA)
  
  # Sum and convert to hectares
  total_m2 <- global(cs_grass_inside, "sum", na.rm = TRUE)[1, 1]
  as.numeric(total_m2) / 10000
}

# Compute total hectares for each raster
grass_ha_results <- purrr::map_dbl(r_list, ~ grass_hectares(.x, sa, cs))
grass_ha_results

# Tidy tibble with results
grass_ha_df <- tibble::tibble(
  raster = names(grass_ha_results),
  grass_ha = unname(grass_ha_results)
)

# Print to console
grass_ha_df

# write to CSV in your outputs directory
readr::write_csv(grass_ha_df, file.path(out_dir, "grass_total_hectares_per_raster_5prov.csv"))

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
