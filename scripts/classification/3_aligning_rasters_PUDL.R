
###################################################################################
##               PREPROCESSING: ALIGNING RASTERS TO REFERENCE                    ##
###################################################################################

##----------------------- Alison Long - January 2026 -------------------------

# step 3: align all grassland layers to classification study area file from step 1
#         to prep for analysis

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
raw_grass_data_dir <- "R:/CanadianGrasslands/WG_Classification/class_inputs"
proc_grass_data_dir <- "R:/CanadianGrasslands/WG_Classification/processed"
out_grass_dir <- "./outputs/grass"

# PUDL directory
pudl_root <- file.path(raw_grass_data_dir, "phjv_pudl")
pudl_aligned_dir <- file.path(proc_grass_data_dir, "phjv_pudl/aligned")
pudl_cropped_dir <- file.path(proc_grass_data_dir, "phjv_pudl/cropped")

####--- Files ---####
# Pilot area
pa_file <- file.path(pa_dir, "pilot_area_aligned_cnvc_Jan2026.tif")

# Study area
sa_file <- file.path(sa_dir, "class_study_area_aligned_cnvc_Dec2025.tif")

# CNVC raster
#cnvc_file <- file.path(raw_grass_data_dir, "ns_cnvc/NatureServe_CNVC_COMPRESSED_20251023.tif")

# PUDL raster
pudl_file <- file.path(pudl_root, "PHJV_PUDLmask.tif")

##### use study area file as snap raster #####

# Load study area raster
sa <- rast(sa_file)

# print raster details 
origin(sa) # 6.696121 -2.550334
res(sa) # 30 30
crs(sa) # EPSG:9001
ext(sa) # SpatExtent : -2321003.30387889, 1585686.69612111, 270327.4496659, 2866497.4496659 (xmin, xmax, ymin, ymax)

####--- ALGIN RASTERS ---#### 

#----- PUDL ----#

# Load rasters
r_pudl <- rast(pudl_file)
r_pudl

#class       : SpatRaster 
#size        : 118947, 169708, 1  (nrow, ncol, nlyr)
#resolution  : 10, 10  (x, y)
#extent      : -1708665, -11585, 2888225, 4077695  (xmin, xmax, ymin, ymax)
#coord. ref. : Albers_Conical_Equal_Area 
#source      : PHJV_PUDLmask.tif 
#name        : PHJV_PUDLmask 
#min value   :             0 
#max value   :             3 

res(r_pudl) # 10 10
origin(r_pudl) # 5 5
crs(r_pudl) # 
ext(r_pudl) # SpatExtent : -1708665, -11585, 2888225, 4077695 (xmin, xmax, ymin, ymax)

plot(r_pudl)

# print frequency table
freq_table <- freq(r_pudl) 
print(freq_table) # 0: 5981498360; 1: 616168759; 2: 193793951; 3: 200728132

# Mark as categorical (factor)
r_pudl <- as.factor(r_pudl)

# reproject before resampling
pudl_reproj <- project(r_pudl, crs(sa)) # or should be project(ngi_reclass, sa) ??

freq_table2 <- freq(pudl_reproj) 
print(freq_table2) # 0: 5981498360; 1: 616168759; 2: 193793951; 3; 200728132

res(pudl_reproj) # 9.576849 9.576849
origin(pudl_reproj) # 0.7826532 0.9658756
crs(pudl_reproj) # EPSG:9001
ext(pudl_reproj) # SpatExtent : -1601851.79178791, -10447.5601285845, 1036206.50267669, 2356652.93202072 (xmin, xmax, ymin, ymax)

# align to reference raster 
pudl_reproj_algn <- resample(pudl_reproj, sa, method = "near")
plot(pudl_reproj_algn)

# Get values and counts
freq_table3 <- freq(pudl_reproj_algn) 
print(freq_table3) # 0: 609557112; 1: 62790661; 2: 19747161; 3: 20456239

res(pudl_reproj_algn) # 30 30
origin(pudl_reproj_algn) # 6.696121 -2.550334
crs(pudl_reproj_algn) # EPSG:9001
ext(pudl_reproj_algn) # SpatExtent : -2321003.30387889, 1585686.69612111, 270327.4496659, 2866497.4496659 (xmin, xmax, ymin, ymax)

# Save the raster as 8-bit unsigned integer
pudl_out_path <- file.path(pudl_aligned_dir, "pudl_reproj_algn.tif")
writeRaster(pudl_reproj_algn, pudl_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# Trim NA-only rows/columns
pudl_reproj_algn_trim <- trim(pudl_reproj_algn)

# visualize the trimmed raster
plot(pudl_reproj_algn_trim)

# print raster details 
origin(pudl_reproj_algn_trim) # 6.696121 -2.550334
res(pudl_reproj_algn_trim) # 30 30 
crs(pudl_reproj_algn_trim) # EPSG:9001
ext(pudl_reproj_algn_trim) # SpatExtent : -1565033.30387889, -10793.3038788857, 1062567.4496659, 2329887.4496659 (xmin, xmax, ymin, ymax)

# Save the raster as 8-bit unsigned integer
pudl_trim_out_path <- file.path(pudl_aligned_dir, "pudl_reproj_algn_trim.tif")
writeRaster(pudl_reproj_algn_trim, pudl_trim_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# reclassify to binary (values 1-10 become 1 for grass)
pudl_algn_binary <- classify(pudl_reproj_algn, 
                             rcl = matrix(c(0, 0, # Potentially Disturbed/Other Cover
                                            1, 1, # Potentially Undisturbed Grass (In PUDL)
                                            2, 1, # Potentially Disturbed Grass (In PUDL)
                                            3, 0), # Shrub (In PUDL)
                                          ncol = 2, byrow = TRUE))

# Save the raster as 8-bit unsigned integer
pudl_binary_out_path <- file.path(pudl_aligned_dir, "pudl_reclass_reproj_algn_binary_noshrub.tif")
writeRaster(pudl_algn_binary, pudl_binary_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# Load study area raster
pa <- rast(pa_file)

# print raster details 
origin(pa) # 6.696121 -2.550334
res(pa) # 30 30
crs(pa)

# crop to pilot area
crop_ext <- ext(pa)
pudl_pa <- crop(pudl_reproj_algn, crop_ext)
pudl_pa
plot (pudl_pa)

# class       : SpatRaster 
# size        : 1367, 1692, 1  (nrow, ncol, nlyr)
# resolution  : 30, 30  (x, y)
# extent      : -828263.3, -777503.3, 1111437, 1152447  (xmin, xmax, ymin, ymax)
# coord. ref. : North_America_Albers_Equal_Area_Conic 
# source(s)   : memory
# varname     : class_study_area_aligned_cnvc_Dec2025 
# categories  : PHJV_PUDLmask 
# name        : PHJV_PUDLmask 
# min value   :             0 
# max value   :             3 

# mask to pilot area
pudl_masked <- mask(pudl_pa, pa)
pudl_masked
plot(pudl_masked)

# class       : SpatRaster 
# size        : 1367, 1692, 1  (nrow, ncol, nlyr)
# resolution  : 30, 30  (x, y)
# extent      : -828263.3, -777503.3, 1111437, 1152447  (xmin, xmax, ymin, ymax)
# coord. ref. : North_America_Albers_Equal_Area_Conic 
# source(s)   : memory
# varname     : class_study_area_aligned_cnvc_Dec2025 
# categories  : PHJV_PUDLmask 
# name        : PHJV_PUDLmask 
# min value   :             0 
# max value   :             3 

# print raster details 
origin(pudl_masked) # 6.696121 -2.550334
res(pudl_masked) # 30 30
crs(pudl_masked)

# Save the raster as 8-bit unsigned integer
pudl_pa_out_path <- file.path(pudl_cropped_dir, "pudl_reproj_algn_crop_pa.tif")
writeRaster(pudl_masked, pudl_pa_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)
