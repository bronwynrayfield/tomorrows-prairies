
###################################################################################
##           ANALYSIS: SPATIAL AGREEMENT SUMMARIES (NATIVE VS TAME)              ##
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
aci_aligned_dir  <- file.path(proc_grass_data_dir, "aafc_aci/aligned")
ngi_aligned_dir  <- file.path(proc_grass_data_dir, "cfga_ngi/aligned")
cgr_masked_dir   <- file.path(proc_grass_data_dir, "cgr_gam_adv/masked")
eccc_masked_dir  <- file.path(proc_grass_data_dir, "eccc_grass/masked")
cnvc_masked_dir  <- file.path(proc_grass_data_dir, "natureserve_cnvc/masked")
pudl_aligned_dir <- file.path(proc_grass_data_dir, "phjv_pudl/aligned")
pli_aligned_dir  <- file.path(proc_grass_data_dir, "sk_pli/aligned")
gpw_masked_dir   <- file.path(proc_grass_data_dir, "wri_gpw_ann/masked")
#plow_masked_dir  <- file.path(proc_grass_data_dir, "wwf_plowprint/masked")
plvi_aligned_dir <- file.path(proc_grass_data_dir, "ab_plvi/aligned")
hfi_aligned_dir  <- file.path(proc_grass_data_dir, "abmi_hfi/aligned")
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
aci_file  <- file.path(aci_aligned_dir,  "aafc_aci_merged_reproj_algn_natvtam_noshrub.tif")
ngi_file  <- file.path(ngi_aligned_dir,  "cfga_ngi_merged_reproj_algn_natvtam_good.tif")
cgr_file  <- file.path(cgr_masked_dir,   "cgr_reproj_algn_natvtam_masked.tif")
eccc_file <- file.path(eccc_masked_dir,  "eccc_reclass_reproj_algn_natvtam_masked.tif")
cnvc_file <- file.path(cnvc_masked_dir,  "cnvc_algn_natvtam_masked.tif")
pudl_file <- file.path(pudl_aligned_dir, "pudl_reclass_reproj_algn_natvtam_noshrub.tif")
pli_file  <- file.path(pli_aligned_dir,  "sk_pli_merged_reproj_algn_natvtam_noshrub.tif")
gpw_file  <- file.path(gpw_masked_dir,   "gpw_reproj_algn_natvtam_noshrub_masked.tif")
#plow_file <- file.path(plow_masked_dir,  "??.tif")
plvi_file <- file.path(plvi_aligned_dir, "plvi_reproj_algn_natvtam_gte2015.tif")
hfi_file  <- file.path(hfi_aligned_dir, "hfi_reproj_algn_binary_gte2015.tif") # tame only - same as binary
mgi_file    <- file.path(mgi_aligned_dir, "mb_mgi_merged_reproj_algn_natvtam_noshrub.tif")

# ==================
# LOAD STACK 
# ==================
files <- c(
  aci_file, ngi_file, cgr_file, eccc_file, cnvc_file, pudl_file, pli_file,  
  gpw_file, plvi_file, hfi_file, mgi_file
)

rast_stack <- rast(files)

sa <- rast(sa_file)

plvi <- rast(plvi_file)
plvi
plot(plvi)

hfi <- rast(hfi_file)
hfi
plot(hfi)

# crop to study area
crop_ext <- ext(sa)
plvi_sa <- crop(plvi, crop_ext)
plvi_sa
plot (plvi_sa, main = "PLVI cropped (native vs tame)")

# mask to pilot area
plvi_masked <- mask(plvi_sa, sa)
plvi_masked
plot(plvi_masked, main = "PLVI masked (native vs tame)")

# print frequency table
freq_table <- freq(plvi_masked) 
print(freq_table)

#  layer value    count
#1     1     0 18310256
#2     1     1  4056463 x 0.09 = 365,081.67 ha tame
#3     1     2   841174 x 0.09 = 75,705.66 ha native

# crop to study area
crop_ext <- ext(sa)
hfi_sa <- crop(hfi, crop_ext)
hfi_sa
plot (hfi_sa, main = "HFI cropped (native vs tame)")

# mask to pilot area
hfi_masked <- mask(hfi_sa, sa)
hfi_masked
plot(hfi_masked, main = "HFI masked (native vs tame)")

# print frequency table
freq_table <- freq(hfi_masked) 
print(freq_table)

#  layer value   count
#1     1     1 1449098 x 0.09 = 130,418.82 ha tame

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
# Tell the future framework to use separate R sessions (processes),
# with one worker per CPU core specified in N_WORKERS
plan(multisession, workers = N_WORKERS)

# Apply a function in parallel over each raster file index
area_ha <- future_sapply(1:length(files), function(i) {
  
  # Load terra inside each worker process
  # (workers do not inherit loaded packages)
  library(terra)
  
  # Limit how much RAM terra can use in this worker
  # This helps prevent memory exhaustion when processing large raster
  terraOptions(memfrac = 0.4) # lower per-worker fraction
  
  # Load the i‑th raster file as a SpatRaster
  # Each worker reads only the file it needs
  stack_i <- rast(files[[i]])
  
  # Load the study area raster (mask) inside the worker
  sa_w <- rast(sa_file)
  
  # Crop the input raster to the study area extent (snap = "near" aligns pixels to minimize resampling)
  lyr <- crop(stack_i, sa_w, snap = "near")
  
  # Mask the cropped raster to the exact study area shape (pixels outside the study area are set to NA)
  lyr <- mask(lyr, sa_w)
  
  # Area for each class
  # Calculate total area (ha) where raster value == 1 ("tame") or == 2 ("native")
  # 1) Convert to logical raster (TRUE/FALSE)
  # 2) Sum TRUE pixels (treated as 1)
  # 3) Multiply by per‑pixel area (hectares)
  tame_ha   <- global(lyr==1, fun = "sum", na.rm = TRUE)[[1]] * px_area_ha
  native_ha <- global(lyr==2, fun = "sum", na.rm = TRUE)[[1]] * px_area_ha
  
  # Remove large raster objects to free worker memory
  # (gc() deliberately NOT called to avoid terra pointer crashes)
  rm(stack_i, sa_w, lyr)
  #gc()
  
  # Return both class areas as a named numeric vector
  # future_sapply() will combine these into a 2 × N matrix
  c(tame_ha = tame_ha, native_ha = native_ha)
  
}, future.seed = NULL) # Disable RNG seeding (not needed; deterministic)

plan(sequential) # Return to normal (non‑parallel) execution

# ==================
# Results
# ==================
# Check structure
str(area_ha)

# Extract values 
tame_ha   <- area_ha["tame_ha", ]
native_ha <- area_ha["native_ha", ]

#Build results table
results <- data.frame(
  layer = names(rast_stack),
  file = basename(files),
  tame_ha = round(tame_ha, 4),
  native_ha = round(native_ha, 4)
)

#readr::write_csv(results, file.path(out_dir, "grass_total_hectares_per_raster_natvtam_modextpr_16Apr2026.csv"))
#write.csv(results, OUTPUT_CSV, row.names = FALSE)

# Write to csv
s3write_using(results, 
              FUN = write.csv, 
              object = "grass_total_hectares_per_raster_natvtam_modextpr_16Apr2026.csv", 
              bucket = "s3://sagemaker-gst-stage.sharing")

# Print results
cat("Done. Results:\n")
print(results)

######################################################################################
######################################################################################

## For Prairies Ecozone instead of modelling extent + Peace River

# Study area (Prairies Ecozone)
sa_file <- file.path(sa_dir, "modelling_study_area_aligned_cnvc_Jan2026.tif")

# ==================
# LOAD STACK 
# ==================
files <- c(
  aci_file, ngi_file, cgr_file, eccc_file, cnvc_file, pudl_file, pli_file,  
  gpw_file, plvi_file, hfi_file, mgi_file
)

rast_stack <- rast(files)

sa <- rast(sa_file)

plvi <- rast(plvi_file)
plvi
plot(plvi)

hfi <- rast(hfi_file)
hfi
plot(hfi)

# crop to study area
crop_ext <- ext(sa)
plvi_sa <- crop(plvi, crop_ext)
plvi_sa
plot (plvi_sa, main = "PLVI cropped (native vs tame)")

# mask to pilot area
plvi_masked <- mask(plvi_sa, sa)
plvi_masked
plot(plvi_masked, main = "PLVI masked (native vs tame)")

# print frequency table
freq_table <- freq(plvi_masked) 
print(freq_table)

#  layer value    count
#1     1     0 17472181
#2     1     1  3843452 x 0.09 = 345,910.68 ha tame
#3     1     2   773768 x 0.09 = 69,639.12 ha native

# crop to study area
crop_ext <- ext(sa)
hfi_sa <- crop(hfi, crop_ext)
hfi_sa
plot (hfi_sa, main = "HFI cropped (native vs tame)")

# mask to pilot area
hfi_masked <- mask(hfi_sa, sa)
hfi_masked
plot(hfi_masked, main = "HFI masked (native vs tame)")

# print frequency table
freq_table <- freq(hfi_masked) 
print(freq_table)

#  layer value   count
#1     1     1  465519 x 0.09 = 41,896.71 ha tame

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
# Tell the future framework to use separate R sessions (processes),
# with one worker per CPU core specified in N_WORKERS
plan(multisession, workers = N_WORKERS)

# Apply a function in parallel over each raster file index
area_ha <- future_sapply(1:length(files), function(i) {
  
  # Load terra inside each worker process
  # (workers do not inherit loaded packages)
  library(terra)
  
  # Limit how much RAM terra can use in this worker
  # This helps prevent memory exhaustion when processing large raster
  terraOptions(memfrac = 0.4) # lower per-worker fraction
  
  # Load the i‑th raster file as a SpatRaster
  # Each worker reads only the file it needs
  stack_i <- rast(files[[i]])
  
  # Load the study area raster (mask) inside the worker
  sa_w <- rast(sa_file)
  
  # Crop the input raster to the study area extent (snap = "near" aligns pixels to minimize resampling)
  lyr <- crop(stack_i, sa_w, snap = "near")
  
  # Mask the cropped raster to the exact study area shape (pixels outside the study area are set to NA)
  lyr <- mask(lyr, sa_w)
  
  # Area for each class
  # Calculate total area (ha) where raster value == 1 ("tame") or == 2 ("native")
  # 1) Convert to logical raster (TRUE/FALSE)
  # 2) Sum TRUE pixels (treated as 1)
  # 3) Multiply by per‑pixel area (hectares)
  tame_ha   <- global(lyr==1, fun = "sum", na.rm = TRUE)[[1]] * px_area_ha
  native_ha <- global(lyr==2, fun = "sum", na.rm = TRUE)[[1]] * px_area_ha
  
  # Remove large raster objects to free worker memory
  # (gc() deliberately NOT called to avoid terra pointer crashes)
  rm(stack_i, sa_w, lyr)
  #gc()
  
  # Return both class areas as a named numeric vector
  # future_sapply() will combine these into a 2 × N matrix
  c(tame_ha = tame_ha, native_ha = native_ha)
  
}, future.seed = NULL) # Disable RNG seeding (not needed; deterministic)

plan(sequential) # Return to normal (non‑parallel) execution

# ==================
# Results
# ==================
# Check structure
str(area_ha)

# Extract values 
tame_ha   <- area_ha["tame_ha", ]
native_ha <- area_ha["native_ha", ]

#Build results table
results <- data.frame(
  layer = names(rast_stack),
  file = basename(files),
  tame_ha = round(tame_ha, 4),
  native_ha = round(native_ha, 4)
)

# Write to csv
s3write_using(results, 
              FUN = write.csv, 
              object = "grass_total_hectares_per_raster_natvtam_pe_16Apr2026.csv", 
              bucket = "s3://sagemaker-gst-stage.sharing")

# Print results
cat("Done. Results:\n")
print(results)
