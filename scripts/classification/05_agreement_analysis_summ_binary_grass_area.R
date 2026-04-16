
###################################################################################
##               ANALYSIS: SPATIAL AGREEMENT SUMMARIES (BINARY)                  ##
###################################################################################

##----------------------- Alison Long - April 2026 -------------------------

# step 5: run descriptive data summaries to assess similarities and differences
#        among existing grassland datasets

# ==================
# SETUP
# ==================

# Load packages
library(terra)
library(future.apply)
library(readr)
library(aws.s3)

# Configuration
N_WORKERS <- max(1, parallel::detectCores() %/% 4)
             #parallel::detectCores() - 2
MEM_FRAC <- 0.4 # fraction of RAM terra can use per worker

# ==================
# DIRECTORIES
# ==================

# pilot area
#pa_dir <- "R:/CanadianGrasslands/processeddatasets/pilotarea"

# study area
#sa_dir <- "R:/CanadianGrasslands/processeddatasets/studyarea"
sa_dir <- "s3://sagemaker-gst-stage.sharing/alison-long/tomorrows-prairies/processed/study-area"

# grasslands
#proc_grass_data_dir <- "R:/CanadianGrasslands/WG_Classification/processed"
proc_grass_data_dir <- "s3://sagemaker-gst-stage.sharing/alison-long/tomorrows-prairies/processed"

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
mgi_aligned_dir   <- file.path(proc_grass_data_dir, "mb_mgi/aligned")

# output 
#out_dir <- "R:/CanadianGrasslands/WG_Classification/results/spatial_agreement_analysis"
out_dir <- "s3://sagemaker-gst-stage.sharing/alison-long/tomorrows-prairies/results/spatial_agreement_analysis"

# ==================
# FILES
# ==================

# ecodistricts
#ecodist_file <- "R:/CanadianGrasslands/processeddatasets/boundaries/national-ecological-framework-for-canada/Ecodistricts/ecodistricts_clean.shp"

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
mgi_file    <- file.path(mgi_aligned_dir, "mb_mgi_reproj_algn_binary_noshrub.tif")

# ==================
# LOAD STACK 
# ==================
files <- c(
  aci_file, luts30_file, luts10_file, cec_file, ngi_file, cgr_file, eccc_file,
  esa_file, glc_file, cnvc_file, lcc_file, pudl_file, pli_file, lcr_file, 
  gpw_file, plow_file, plvi_file, hfi_file, gcc_file, dw_file, esri_file, mgi_file
)

rast_stack <- rast(files)

sa <- rast(sa_file)

# ==================
# CHECK CRS
# ==================
# 1. CRS must exist
if (is.na(crs(rast_stack))) {
  stop("Raster has no CRS defined. Area calculations would be invalid.")
}

# 2. Must not be geographic (lon/lat)
if (isTRUE(is.lonlat(rast_stack))) {
  stop(
    "Raster is in a geographic (lon/lat) CRS. ",
    "Reproject to an equal-area projected CRS before computing area."
  )
}

# 3. CRS must be projected; units ideally meters
crs_units <- tryCatch(
  terra::units(rast_stack),
  error = function(e) NA
)

if (all(is.na(crs_units))) {
  warning(
    "CRS units could not be determined. ",
    "Proceeding because raster is projected (not lon/lat). ",
    "Ensure CRS is equal-area with meter units."
  )
} else if (any(crs_units == "")) {
  warning(
    "CRS units are unspecified in metadata. ",
    "Assuming projected CRS uses linear units (typically meters). ",
    "Verify projection is equal-area."
  )
} else if (!all(tolower(crs_units) %in% c("m", "meter", "metre", "meters", "metres"))) {
  stop(sprintf(
    "CRS units are '%s', not meters. Reproject to an equal-area CRS in meters.",
    paste(unique(crs_units), collapse = ", ")
  ))
}

cat("CRS check passed:\n")
cat(sprintf("  CRS: %s\n", crs(rast_stack)))
cat(sprintf("  Units: %s\n\n", crs_units))

# ==================
# Pixel area scalar (equal-area projection, computed once)
# ==================
px_area_ha <- prod(res(rast_stack)) / 10000
cat(sprintf("Pixel area: %.6f ha\n", px_area_ha))
cat(sprintf("Processing %d layers across %d workers\n\n", nlyr(rast_stack), N_WORKERS))

# ==================
# Parallelize across layers
# ==================
plan(multisession, workers = N_WORKERS)

area_ha <- future_sapply(1:length(files), function(i) {
  library(terra)
  terraOptions(memfrac = 0.4) # lower per-worker fraction
  
  stack_i <- rast(files[[i]])
  sa_w <- rast(sa_file)
  
  lyr <- crop(stack_i, sa_w, snap = "near")
  lyr <- mask(lyr, sa_w)
  result <- global(lyr, fun = "sum", na.rm = TRUE)[[1]] * px_area_ha
  
  rm(stack_i, sa_w, lyr)
  #gc()
  
  result
}, future.seed = NULL)

plan(sequential)


# ==================
# Results
# ==================
results <- data.frame(
  layer = names(rast_stack),
  file = basename(files),
  area_ha = round(area_ha, 4)
)

#readr::write_csv(results, file.path(out_dir, "grass_total_hectares_per_raster_16Apr2026.csv"))
#write.csv(results, OUTPUT_CSV, row.names = FALSE)

s3write_using(results, 
              FUN = write.csv, 
              object = "grass_total_hectares_per_raster_16Apr2026.csv", 
              bucket = "s3://sagemaker-gst-stage.sharing")

cat("Done. Results:\n")
print(results)

#############################################################################################
#############################################################################################

## For Prairies Ecozone instead of modelling extent + Peace River

# Study area (Prairies Ecozone)
sa_file <- file.path(sa_dir, "modelling_study_area_aligned_cnvc_Jan2026.tif")

# ==================
# LOAD STACK 
# ==================
files <- c(
  aci_file, luts30_file, luts10_file, cec_file, ngi_file, cgr_file, eccc_file,
  esa_file, glc_file, cnvc_file, lcc_file, pudl_file, pli_file, lcr_file, 
  gpw_file, plow_file, plvi_file, hfi_file, gcc_file, dw_file, esri_file, mgi_file
)

rast_stack <- rast(files)

sa <- rast(sa_file)

# ==================
# CHECK CRS
# ==================
# 1. CRS must exist
if (is.na(crs(rast_stack))) {
  stop("Raster has no CRS defined. Area calculations would be invalid.")
}

# 2. Must not be geographic (lon/lat)
if (isTRUE(is.lonlat(rast_stack))) {
  stop(
    "Raster is in a geographic (lon/lat) CRS. ",
    "Reproject to an equal-area projected CRS before computing area."
  )
}

# 3. CRS must be projected; units ideally meters
crs_units <- tryCatch(
  terra::units(rast_stack),
  error = function(e) NA
)

if (all(is.na(crs_units))) {
  warning(
    "CRS units could not be determined. ",
    "Proceeding because raster is projected (not lon/lat). ",
    "Ensure CRS is equal-area with meter units."
  )
} else if (any(crs_units == "")) {
  warning(
    "CRS units are unspecified in metadata. ",
    "Assuming projected CRS uses linear units (typically meters). ",
    "Verify projection is equal-area."
  )
} else if (!all(tolower(crs_units) %in% c("m", "meter", "metre", "meters", "metres"))) {
  stop(sprintf(
    "CRS units are '%s', not meters. Reproject to an equal-area CRS in meters.",
    paste(unique(crs_units), collapse = ", ")
  ))
}

cat("CRS check passed:\n")
cat(sprintf("  CRS: %s\n", crs(rast_stack)))
cat(sprintf("  Units: %s\n\n", crs_units))

# ==================
# Pixel area scalar (equal-area projection, computed once)
# ==================
px_area_ha <- prod(res(rast_stack)) / 10000
cat(sprintf("Pixel area: %.6f ha\n", px_area_ha))
cat(sprintf("Processing %d layers across %d workers\n\n", nlyr(rast_stack), N_WORKERS))

# ==================
# Parallelize across layers
# ==================
plan(multisession, workers = N_WORKERS)

area_ha <- future_sapply(1:length(files), function(i) {
  library(terra)
  terraOptions(memfrac = 0.4) # lower per-worker fraction
  
  stack_i <- rast(files[[i]])
  sa_w <- rast(sa_file)
  
  lyr <- crop(stack_i, sa_w, snap = "near")
  lyr <- mask(lyr, sa_w)
  result <- global(lyr, fun = "sum", na.rm = TRUE)[[1]] * px_area_ha
  
  rm(stack_i, sa_w, lyr)
  #gc()
  
  result
}, future.seed = NULL)

plan(sequential)

# ==================
# Results
# ==================
results <- data.frame(
  layer = names(rast_stack),
  file = basename(files),
  area_ha = round(area_ha, 4)
)

s3write_using(results, 
              FUN = write.csv, 
              object = "grass_total_hectares_per_raster_pe_16Apr2026.csv", 
              bucket = "s3://sagemaker-gst-stage.sharing")

cat("Done. Results:\n")
print(results)
