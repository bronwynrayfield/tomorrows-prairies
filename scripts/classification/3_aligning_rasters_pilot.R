
###################################################################################
##               PREPROCESSING: ALIGNING RASTERS - PILOT AREA                    ##
###################################################################################

##----------------------- Alison Long - January 2026 -------------------------

# step 3: align all grassland layers to classification study area file from step 1
#         and crop to pilot area to prep for analysis

# load packages
library(terra)
library(sf)
library(tidyverse)

# Set the temporary directory to a drive with sufficient space
terraOptions(tempdir = "R:/CanadianGrasslands")

####--- DIRECTORIES ---####

# pilot area
pa_dir <- "./processeddatasets/pilot_area"

# study area
sa_dir <- "./processeddatasets/Studyarea"

# grasslands
raw_grass_data_dir <- "./existingdatasets/class_inputs"
proc_grass_data_dir <- "./processeddatasets/grass"
out_grass_dir <- "./outputs/grass"

# NGI tiles directory
ngi_tile_root <- file.path(raw_grass_data_dir, "cfga_ngi")
ngi_merged_dir <- file.path(proc_grass_data_dir, "cfga_ngi/merged")
ngi_cropped_dir <- file.path(proc_grass_data_dir, "cfga_ngi/cropped")
ngi_aligned_dir <- file.path(proc_grass_data_dir, "cfga_ngi/aligned")

# ECCC directory
eccc_root <- file.path(raw_grass_data_dir, "eccc_grass_summary")
eccc_cropped_dir <- file.path(proc_grass_data_dir, "eccc_grass_summary/cropped")

# PHJV PUDL directory
pudl_root <- file.path(raw_grass_data_dir, "phjv_pudl")
pudl_cropped_dir <- file.path(proc_grass_data_dir, "phjv_pudl/cropped")

# Central Grasslands Roadmap directory
cgr_root <- file.path(raw_grass_data_dir, "cgr_gam_adv")
cgr_cropped_dir <- file.path(proc_grass_data_dir, "cgr_gam_adv/cropped")

# WRI annual grass directory
wri_root <- file.path(raw_grass_data_dir, "wri_annual_grass")
wri_cropped_dir <- file.path(proc_grass_data_dir, "wri_annual_grass/cropped")


####--- Files ---####
# Pilot Area
pa_file <- file.path(pa_dir, "GNP_pilot_area_crop_CSD.shp")

# Study Area
sa_file <- file.path(sa_dir, "class_study_area_aligned_cnvc_Dec2025.tif")

# CNVC raster
cnvc_file <- file.path(raw_grass_data_dir, "ns_cnvc/NatureServe_CNVC_COMPRESSED_20251023.tif")

# NGI rasters
ngi_ab_file <- file.path(ngi_tile_root, "AB-002.tif")
ngi_mb_file <- file.path(ngi_tile_root, "MB-003.tif")
ngi_sk_file <- file.path(ngi_tile_root, "SK-001.tif")

# AAFC LUTS 

# ECCC grass summary
eccc_file <- file.path(eccc_root, "grass_summary_a_class.tif")

# WRI
wri_file <- file.path(wri_root, "gpw_grassland_rf.savgol.bthr_c_30m_20220101_20221231_go_epsg.4326_v1.tif")

# PUDL
pudl_file <- file.path(pudl_root, "PHJV_PUDLmask.tif")

# Central Grasslands Roadmap - Grassland Assessment Map (Advanced)
cgr_file <- file.path(cgr_root, "CGR_GAM_V2_separated.tif")

##############################################

# prep rasters for alignment

#----- merge NGI rasters ----

# Load rasters
r_ab <- rast(ngi_ab_file)
r_mb <- rast(ngi_mb_file)
r_sk <- rast(ngi_sk_file)

res(r_ab)
res(r_mb)
res(r_sk) # 10 x 10
origin(r_ab)
origin(r_mb)
origin(r_sk) # 0 0 

ext(r_ab)
ext(r_mb)
ext(r_sk)
nrow(r_ab); ncol(r_ab)
nrow(r_mb); ncol(r_mb)
nrow(r_sk); ncol(r_sk)

# Set pixels with value 0 to NA
r_ab[r_ab == 0] <- NA
r_mb[r_mb == 0] <- NA
r_sk[r_sk == 0] <- NA

# Merge using mosaic with max function
merged_ngi <- mosaic(r_ab, r_mb, r_sk, fun = max)

res(merged_ngi)
origin(merged_ngi)

# Optional: plot the result
plot(merged_ngi, main = "Merged NGI Raster")

# Save the raster as 8-bit unsigned integer
ngi_merge_out_path <- file.path(ngi_merged_dir, "cfga_ngi_merged_no0s.tif")
writeRaster(merged_ngi, ngi_merge_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

##############################################

# use study area file as snap raster

# Load study area raster
sa <- rast(sa_file)

# print raster details 
origin(sa)
res(sa)
crs(sa)

# visualize study area raster
plot(sa)

# Load pilot area shapefile
pa_vect <- vect(pa_file)

# Load classification study area shapefile
#sa_vect <- vect(sa_file)

# Reproject pilot area to match study area raster CRS
pa_vect_proj <- project(pa_vect, crs(sa))

# Rasterize pilot area using CNVC raster as template
pa_rast <- rasterize(pa_vect_proj, sa, field=1, touches=TRUE) # field=1 assigns constant value to rasterized cells
# touches=TRUE ensures that any cell touched by the polygon is included

# print raster details 
origin(pa_rast)
res(pa_rast)
crs(pa_rast)

# Get frequency table of raster values
freq_table <- freq(pa_rast)

# Print the table - 1723584 pixels with value 1
print(freq_table)

# visualize the new raster
plot(pa_rast)

# Trim NA-only rows/columns
pa_rast_trimmed <- trim(pa_rast)
plot(pa_rast_trimmed)

# visualize the trimmed raster
plot(pa_rast_trimmed)

# print raster details 
origin(pa_rast_trimmed)
res(pa_rast_trimmed)
crs(pa_rast_trimmed)

# Compare extents and columns/rows
ext(pa_rast)
ext(pa_rast_trimmed)

nrow(pa_rast); ncol(pa_rast)
nrow(pa_rast_trimmed); ncol(pa_rast_trimmed)

# Save the raster as 8-bit unsigned integer
pa_out_path <- file.path(pa_dir, "pilot_area_test_aligned_cnvc_trim_crop_CSD.tif")
writeRaster(pa_rast_trimmed, pa_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# Convert raster to polygons (excluding NA values)
pa_polygons <- as.polygons(pa_rast_trimmed, na.rm = TRUE, dissolve = TRUE)

# Export shapefile
pa_out_shp <- file.path(pa_dir, "pilot_area_test_aligned_cnvc_crop_CSD.shp")
writeVector(pa_polygons, pa_out_shp, overwrite = TRUE)

#------------------ CROP RASTERS TO PILOT AREA -------------------------

# designate file paths
#merged_ngi_path <- file.path(ngi_merged_dir, "cfga_ngi_merged_no0s.tif")
#pa_rast_path <- file.path(pa_dir, "pilot_area_test_aligned_cnvc.tif")

# read rasters
#merged_ngi <- rast(merged_ngi_path)
#pa_rast <- rast(pa_rast_path)

# reproject raster to desired crs
#merged_ngi <- project(merged_ngi, crs(pa_rast))

# crop raster to pilot area 
#cropped_ngi <- crop(merged_ngi, ext(pa_rast))

# align cropped raster to pilot area (resolution, origin, CRS)
#crp_algn_ngi <- resample(cropped_ngi, pa_rast, method = "near") # nearest neighbor for categorical

# check alignment - returns TRUE if rasters are aligned in terms of CRS, resolution, origin, and extent
#compareGeom(pa_rast, crp_algn_ngi, stopOnError = FALSE)

# print raster details 
#origin(crp_algn_ngi)
#res(crp_algn_ngi)
#crs(crp_algn_ngi)

# visualize CNVC raster
#plot(crp_algn_ngi)

# Save the raster as 8-bit unsigned integer
#ngi_crp_algn_out_path <- file.path(ngi_cropped_dir, "cfga_ngi_crp_algn_cnvc.tif")
#writeRaster(crp_algn_ngi, ngi_crp_algn_out_path, datatype = "INT1U", 
#            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
#            NAflag = 255,
#            overwrite = TRUE)

##### use SK raster only since mosaic is producing different resolution

# import pilot area rater created above so you don't have to repeat steps
pa_in_path <- file.path(pa_dir, "pilot_area_test_aligned_cnvc_trim_crop_CSD.tif")
pa_rast <- rast(pa_in_path)
plot(pa_rast)
origin(pa_rast)
res(pa_rast)
crs(pa_rast)

# Load raster
r_sk <- rast(ngi_sk_file)

# 1) Reproject source to the target CRS, keeping the source's native resolution
#    Pass the CRS string and set res=res(source) to avoid adopting target's resolution.
sk_rp <- project(
  x      = r_sk,
  y      = crs(pa_rast),
  res    = res(r_sk),
  method = "near"   # for categorical data
)

origin(sk_rp) # -2.668551 -1.848842
res(sk_rp) # 10 10
crs(sk_rp)

# 2) Align the origin to match the target
origin(sk_rp) <- origin(pa_rast)

origin(sk_rp) # -3.303879 -2.550334
res(sk_rp) # 10 10
crs(sk_rp)

##### another way

# Choose interpolation method: "near" for categorical data
method <- "near"     

# --- Create a template raster ---
# The template defines: target CRS + target extent + source resolution
tmpl <- rast(
  ext = ext(pa_rast),   # crop area will match target's rectangle
  crs = crs(pa_rast),   # CRS to match target
  res = res(r_sk)       # keep native resolution from source
)

# --- Project the source to the template ---
# This simultaneously reprojects, aligns grid, and (effectively) crops to target extent
sk_rp <- project(r_sk, tmpl, method = method)

origin(sk_rp) # -3.303879 -2.550334
res(sk_rp) # 10 10
crs(sk_rp)

origin(sk_rp) <- origin(pa_rast)

############ another way

# Choose interpolation method:
method <- "near"     # use "near" for categorical; "bilinear" for continuous

# --- Build a template from the target's grid ---
# This clones target's CRS, origin, and extent.
tmpl2 <- rast(pa_in_path)

# Replace the template's resolution with the source's native resolution.
# This preserves the target’s origin but increases/decreases ncol/nrow to match the new res.
res(tmpl2) <- res(r_sk)

# --- Project the source onto the template grid ---
# Output now has: target CRS, target origin, source's resolution, target extent (cropped)
out <- project(r_sk, tmpl2, method = method)

origin(out) # -3.303879 -2.550334
res(out) # 10 10
crs(out)

#############################

# Can have native resolution or identical numeric origin values - cannot have both unless 
# the resolutions are the same. CAN ensure (and usually what matters for alignment) 
# is that xmin/ymin and CRS match the target, so the grid is anchored the same way

method <- "near"  # "near" for categorical; "bilinear" for continuous

# Build a template that clones the target's CRS and extent (thus xmin/ymin),
# but uses the source's resolution
tmpl3 <- rast(ext = ext(pa_rast), crs = crs(pa_rast), res = res(r_sk))

# Project the source onto the template's grid (reproject + align + crop)
out3 <- project(r_sk, tmpl3, method = method)

origin(out3) # -3.303879 -2.550334
res(out3) # 10 10
crs(out3)

# Checking true alignment via coordinates, not numeric origin - TRUE!
stopifnot(identical(crs(out3), crs(pa_rast)))
stopifnot(all.equal(res(out3), res(r_sk)))
stopifnot(all.equal(xmin(out3), xmin(pa_rast)))
stopifnot(all.equal(ymin(out3), ymin(pa_rast)))
stopifnot(all.equal(xmax(out3), xmax(pa_rast)))
stopifnot(all.equal(ymax(out3), ymax(pa_rast)))
stopifnot(all.equal(ext(out3), ext(pa_rast)))
#stopifnot(identical(origin(out3), origin(pa_rast))) # not true

plot(out3)

# Save the raster as 8-bit unsigned integer
test3_out_path <- file.path(ngi_aligned_dir, "sk_aligned_cnvc_out3.tif")
writeRaster(out3, test3_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# might be aligned but not cropped? - SAME OUTPUT AS out3
# crop raster to pilot area 
#out3_crop <- crop(out3, ext(pa_rast))
#test3_crop_out_path <- file.path(ngi_aligned_dir, "sk_aligned_cnvc_out3_crop.tif")
#writeRaster(out3_crop, test3_crop_out_path, datatype = "INT1U", 
#            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
#            NAflag = 255,
#            overwrite = TRUE)

# need to mask? must have same number of rows/columns to mask so must resample
out3_resampled <- resample(out3, pa_rast)
out3_resamp_mask <- mask(out3_resampled, pa_rast)
plot(out3_resamp_mask)

# Save the raster as 8-bit unsigned integer
test3_mask_out_path <- file.path(ngi_aligned_dir, "sk_aligned_cnvc_out3_resampled_masked.tif")
writeRaster(out3_resamp_mask, test3_mask_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)





# reproject NGI SK raster to match CRS of pilot area
ngi_sk_proj <- project(r_sk, crs(pa_rast))

# check output
origin(ngi_sk_proj) # diff from study area
res(ngi_sk_proj) # prior to reprojecting was 10x10 but now 9.868212x9.868212
crs(ngi_sk_proj)

# align raster to pilot area (resolution, origin, CRS) BUT DON'T WANT TO CHANGE RESOLUTION
ngi_sk_proj_algn <- resample(ngi_sk_proj, pa_rast, method = "near") # nearest neighbor for categorical

# check output
origin(ngi_sk_proj_algn) 
res(ngi_sk_proj_algn) 
crs(ngi_sk_proj_algn)

# crop raster to pilot area 
ngi_sk_proj_algn_crop <- crop(ngi_sk_proj_algn, ext(pa_rast))

# snap origin of NGI raster to match pilot area while keeping its resolution
#ngi_sk_proj_algn <- align(ngi_sk_proj, pa_rast, snap = "near")

# crop to extent of pilot area
#cropped_ngi_sk <- crop(aligned_ngi_sk, ext(pa_rast))

# print raster details
origin(ngi_sk_proj_crop_algn)
res(ngi_sk_proj_crop_algn)
crs(ngi_sk_proj_crop_algn)

# visualize
plot(ngi_sk_proj_crop_algn)

# Save the raster
ngi_crp_algn_out_path <- file.path(ngi_cropped_dir, "cfga_ngi_crp_algn_native_res.tif")
writeRaster(ngi_sk_proj_crop_algn, ngi_crp_algn_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)



#------------------- ALIGN RASTERS WITH CNVC ---------------------------


