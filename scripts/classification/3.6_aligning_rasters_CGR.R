
###################################################################################
##               PREPROCESSING: ALIGNING RASTERS TO REFERENCE                    ##
###################################################################################

##----------------------- Alison Long - January 2026 -------------------------##

# step 3: align all grassland layers to classification study area file from step 1
#         to prep for analysis

# load packages
library(terra)
library(sf)
library(tidyverse)

# Set the temporary directory to a drive with sufficient space
terraOptions(tempdir = "R:/CanadianGrasslands/WG_classification/github/tomorrows-prairies/scratch")

####--- DIRECTORIES ---####

# pilot area
pa_dir <- "R:/CanadianGrasslands/processeddatasets/pilotarea"

# study area
sa_dir <- "R:/CanadianGrasslands/processeddatasets/studyarea"

# grasslands
raw_grass_data_dir <- "R:/CanadianGrasslands/WG_Classification/class_inputs"
proc_grass_data_dir <- "R:/CanadianGrasslands/WG_Classification/processed"
out_grass_dir <- "./outputs/grass"

# cgr directory
cgr_root <- file.path(raw_grass_data_dir, "cgr_gam_adv")
cgr_aligned_dir <- file.path(proc_grass_data_dir, "cgr_gam_adv/aligned")
cgr_cropped_dir <- file.path(proc_grass_data_dir, "cgr_gam_adv/cropped")
cgr_masked_dir <- file.path(proc_grass_data_dir, "cgr_gam_adv/masked")
cgr_pilot_dir <- file.path(proc_grass_data_dir, "cgr_gam_adv/pilot")

####--- Files ---####
# Pilot area
pa_file <- file.path(pa_dir, "pilot_area_aligned_cnvc_Jan2026.tif")

# Study area
sa_file <- file.path(sa_dir, "class_study_area_aligned_cnvc_Dec2025.tif")

# CNVC raster
#cnvc_file <- file.path(raw_grass_data_dir, "ns_cnvc/NatureServe_CNVC_COMPRESSED_20251023.tif")

# cgr raster
cgr_file <- file.path(cgr_root, "CGR_GAM_V2_separated.tif")

##### use study area file as snap raster #####

# Load study area raster
sa <- rast(sa_file)

# print raster details 
origin(sa) # 6.696121 -2.550334
res(sa) # 30 30
crs(sa) # EPSG:9001
ext(sa) # SpatExtent : -2321003.30387889, 1585686.69612111, 270327.4496659, 2866497.4496659 (xmin, xmax, ymin, ymax)

####--- ALGIN RASTERS ---#### 

#----- Central Grasslands Roadmap ----#

# Load rasters
r_cgr <- rast(cgr_file)
r_cgr

#class       : SpatRaster 
#size        : 126119, 68331, 1  (nrow, ncol, nlyr)
#resolution  : 29.99976, 29.99976  (x, y)
#extent      : -1568228, 481685.7, -229033, 3554507  (xmin, xmax, ymin, ymax)
#coord. ref. : NAD_1983_Albers 
#source      : CGR_GAM_V2_separated.tif 
#categories  : acres, percentage, Category, Count 
#name        :             Category 
#min value   :                      
#max value   : Masks (all combined) 

res(r_cgr) # 29.99976 29.99976
origin(r_cgr) # 9.504995 -14.844033
crs(r_cgr) # NAD_1983_Albers
ext(r_cgr) # SpatExtent : -1568228.07428792, 481685.690035148, -229033.030168777, 3554507.00352985 (xmin, xmax, ymin, ymax)

plot(r_cgr)

# print frequency table
freq_table <- freq(r_cgr) 
print(freq_table) 

#   layer                     value      count
#1     1                           4230908988
#2     1                      Core  835084749
#3     1   Vulnerable - conversion  348066885
#4     1 Vulnerable - encroachment  180846264
#5     1         Vulnerable - both   75696446
#6     1                 Plowprint 1114867146
#7     1                Encroached  775568223
#8     1         Plowed/encroached   67700100
#9     1      Masks (all combined)  989098588

# Mark as categorical (factor)
r_cgr <- as.factor(r_cgr)

# reproject before resampling
cgr_reproj <- project(r_cgr, crs(sa)) # or should be project(ngi_reclass, sa) ??

freq_table2 <- freq(cgr_reproj) 
print(freq_table2) 

#layer value      count
#1     1     0 4058709080
#2     1     5  801099734
#3     1   100  333886524
#4     1   111  173482878
#5     1   200   72617949
#6     1   500 1069507316
#7     1   550  744008209
#8     1   800   64948001
#9     1  1000  948823436

res(cgr_reproj) # 30.62956 30.62956
origin(cgr_reproj) # -8.873951 -12.406117
crs(cgr_reproj) # EPSG:9001
ext(cgr_reproj) # SpatExtent : -1521318.08906894, 466969.468113426, -2209537.30848753, 1781893.14471371 (xmin, xmax, ymin, ymax)

# align to reference raster 
cgr_reproj_algn <- resample(cgr_reproj, sa, method = "near")
plot(cgr_reproj_algn)

# Get values and counts
freq_table3 <- freq(cgr_reproj_algn) 
print(freq_table3) 

#  layer value      count
#1     1     0 1637719493
#2     1     5  352074725
#3     1   100  212473664
#4     1   111   55738531
#5     1   200   34742860
#6     1   500  641427421
#7     1   550  128468432
#8     1   800   20419641
#9     1  1000  121810555

res(cgr_reproj_algn) # 30 30
origin(cgr_reproj_algn) # 6.696121 -2.550334
crs(cgr_reproj_algn) # EPSG:9001
ext(cgr_reproj_algn) # SpatExtent : -2321003.30387889, 1585686.69612111, 270327.4496659, 2866497.4496659 (xmin, xmax, ymin, ymax)

# Save the raster as 16-bit unsigned integer
cgr_out_path <- file.path(cgr_aligned_dir, "cgr_reproj_algn.tif")
writeRaster(cgr_reproj_algn, cgr_out_path, datatype = "INT2U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 65535,
            overwrite = TRUE)

# Trim NA-only rows/columns
cgr_reproj_algn_trim <- trim(cgr_reproj_algn)

# visualize the trimmed raster
plot(cgr_reproj_algn_trim)

# print raster details 
origin(cgr_reproj_algn_trim) # 6.696121 -2.550334
res(cgr_reproj_algn_trim) # 30 30 
crs(cgr_reproj_algn_trim) # EPSG:9001
ext(cgr_reproj_algn_trim) # SpatExtent : -1487093.30387889, 455916.696121114, 270327.4496659, 1781907.4496659 (xmin, xmax, ymin, ymax)

# Save the raster as 16-bit unsigned integer
cgr_trim_out_path <- file.path(cgr_aligned_dir, "cgr_reproj_algn_trim.tif")
writeRaster(cgr_reproj_algn_trim, cgr_trim_out_path, datatype = "INT2U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 65535,
            overwrite = TRUE)

# reclassify to binary (values 1-10 become 1 for grass)
#cgr_algn_binary <- classify(cgr_reproj_algn,
cgr_algn_binary <- classify(cgr_algn,
                             rcl = matrix(c(0, NA, # 
                                            5, 1, # Core grasslands
                                            100, 1, # Vulnerable - conversion
                                            111, 1, # Vulnerable - encroachment
                                            200, 1, # Vulnerable - both
                                            500, 0, # Plowprint
                                            550, 0, # Encroached
                                            800, 0, # Plowed/encroached
                                            1000, 0), # Masks
                                          ncol = 2, byrow = TRUE))

plot(cgr_algn_binary)
cgr_algn_binary

#class       : SpatRaster 
#size        : 86539, 130223, 1  (nrow, ncol, nlyr)
#resolution  : 30, 30  (x, y)
#extent      : -2321003, 1585687, 270327.4, 2866497  (xmin, xmax, ymin, ymax)
#coord. ref. : North_America_Albers_Equal_Area_Conic 
#source      : spat_5ba705fbb3cf3_375408_DvggOkjkTOQeAAI.tif 
#varname     : class_study_area_aligned_cnvc_Dec2025 
#name        : Category 
#min value   :        0 
#max value   :        1

# Save the raster as 8-bit unsigned integer
cgr_binary_out_path <- file.path(cgr_aligned_dir, "cgr_reproj_algn_binary.tif")
writeRaster(cgr_algn_binary, cgr_binary_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# load raster
cgr <- file.path(cgr_aligned_dir, "cgr_reproj_algn.tif")
cgr_algn <- rast(cgr)

# reclassify to native vs. tame (0 = not grass, 1 = tame, 2 = native)
cgr_algn_nattam <- classify(cgr_algn, 
                            rcl = matrix(c(0, NA, # 
                                           5, 2, # Core grasslands
                                           100, 2, # Vulnerable - conversion
                                           111, 1, # Vulnerable - encroachment
                                           200, 1, # Vulnerable - both
                                           500, 0, # Plowprint
                                           550, 0, # Encroached
                                           800, 0, # Plowed/encroached
                                           1000, 0), # Masks
                                         ncol = 2, byrow = TRUE))
plot(cgr_algn_nattam)
cgr_algn_nattam

# Save the raster as 8-bit unsigned integer
cgr_nattam_out_path <- file.path(cgr_aligned_dir, "cgr_reproj_algn_natvtam.tif")
writeRaster(cgr_algn_nattam, cgr_nattam_out_path, datatype = "INT1U", 
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
cgr_pa <- crop(cgr_reproj_algn, crop_ext)
cgr_pa
plot (cgr_pa)

#class       : SpatRaster 
#size        : 1367, 1692, 1  (nrow, ncol, nlyr)
#resolution  : 30, 30  (x, y)
#extent      : -828263.3, -777503.3, 1111437, 1152447  (xmin, xmax, ymin, ymax)
#coord. ref. : North_America_Albers_Equal_Area_Conic 
#source(s)   : memory
#varname     : class_study_area_aligned_cnvc_Dec2025 
#categories  : Category 
#name        : Category 
#min value   :        5 
#max value   :     1000 

# mask to pilot area
cgr_masked <- mask(cgr_pa, pa)
cgr_masked
plot(cgr_masked)

# print raster details 
origin(cgr_masked) # 6.696121 -2.550334
res(cgr_masked) # 30 30
crs(cgr_masked)

# Save the raster as 16-bit unsigned integer
cgr_pa_out_path <- file.path(cgr_pilot_dir, "cgr_reproj_algn_crop_pa.tif")
writeRaster(cgr_masked, cgr_pa_out_path, datatype = "INT2U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 65535,
            overwrite = TRUE)

# load binary raster
cgr_binary_file <- file.path(cgr_aligned_dir, "cgr_reclass_reproj_algn_binary.tif")
cgr_algn_binary <- rast(cgr_binary_file)

# crop binary raster to pilot area
cgr_bin_pa <- crop(cgr_algn_binary, crop_ext)
cgr_bin_pa
plot (cgr_bin_pa)

# mask to pilot area
cgr_bin_masked <- mask(cgr_bin_pa, pa)
cgr_bin_masked
plot(cgr_bin_masked)

# Save the raster as 8-bit unsigned integer
cgr_bin_pa_out_path <- file.path(cgr_pilot_dir, "cgr_reproj_algn_crop_pa_binary.tif")
writeRaster(cgr_bin_masked, cgr_bin_pa_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# crop native vs. tame raster to pilot area
cgr_nt_pa <- crop(cgr_algn_nattam, crop_ext)
cgr_nt_pa
plot (cgr_nt_pa)

# mask to pilot area
cgr_nt_masked <- mask(cgr_nt_pa, pa)
cgr_nt_masked
plot(cgr_nt_masked)

# Save the raster as 8-bit unsigned integer
cgr_nt_pa_out_path <- file.path(cgr_pilot_dir, "cgr_reproj_algn_crop_pa_natvtam.tif")
writeRaster(cgr_nt_masked, cgr_nt_pa_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

################################################

# isolate each class

# load raster
cgr <- file.path(cgr_aligned_dir, "cgr_reproj_algn.tif")
cgr_algn <- rast(cgr)

# reclassify 
cgr_algn_core <- classify(cgr_algn, 
                            rcl = matrix(c(0, NA, # 
                                           5, 1, # Core grasslands
                                           100, 0, # Vulnerable - conversion
                                           111, 0, # Vulnerable - encroachment
                                           200, 0, # Vulnerable - both
                                           500, 0, # Plowprint
                                           550, 0, # Encroached
                                           800, 0, # Plowed/encroached
                                           1000, 0), # Masks
                                         ncol = 2, byrow = TRUE))
plot(cgr_algn_core)
cgr_algn_core

# Save the raster as 8-bit unsigned integer
cgr_core_out_path <- file.path(cgr_aligned_dir, "cgr_reproj_algn_core_grass.tif")
writeRaster(cgr_algn_core, cgr_core_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# reclassify 
cgr_algn_vulconv <- classify(cgr_algn, 
                          rcl = matrix(c(0, NA, # 
                                         5, 0, # Core grasslands
                                         100, 1, # Vulnerable - conversion
                                         111, 0, # Vulnerable - encroachment
                                         200, 0, # Vulnerable - both
                                         500, 0, # Plowprint
                                         550, 0, # Encroached
                                         800, 0, # Plowed/encroached
                                         1000, 0), # Masks
                                       ncol = 2, byrow = TRUE))
plot(cgr_algn_vulconv)
cgr_algn_vulconv

# Save the raster as 8-bit unsigned integer
cgr_vulconv_out_path <- file.path(cgr_aligned_dir, "cgr_reproj_algn_vuln_conv.tif")
writeRaster(cgr_algn_vulconv, cgr_vulconv_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# reclassify 
cgr_algn_vulencr <- classify(cgr_algn, 
                             rcl = matrix(c(0, NA, # 
                                            5, 0, # Core grasslands
                                            100, 0, # Vulnerable - conversion
                                            111, 1, # Vulnerable - encroachment
                                            200, 0, # Vulnerable - both
                                            500, 0, # Plowprint
                                            550, 0, # Encroached
                                            800, 0, # Plowed/encroached
                                            1000, 0), # Masks
                                          ncol = 2, byrow = TRUE))
plot(cgr_algn_vulencr)
cgr_algn_vulencr

# Save the raster as 8-bit unsigned integer
cgr_vulencr_out_path <- file.path(cgr_aligned_dir, "cgr_reproj_algn_vuln_encrch.tif")
writeRaster(cgr_algn_vulencr, cgr_vulencr_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# reclassify 
cgr_algn_vulboth <- classify(cgr_algn, 
                             rcl = matrix(c(0, NA, # 
                                            5, 0, # Core grasslands
                                            100, 0, # Vulnerable - conversion
                                            111, 0, # Vulnerable - encroachment
                                            200, 1, # Vulnerable - both
                                            500, 0, # Plowprint
                                            550, 0, # Encroached
                                            800, 0, # Plowed/encroached
                                            1000, 0), # Masks
                                          ncol = 2, byrow = TRUE))
plot(cgr_algn_vulboth)
cgr_algn_vulboth

# Save the raster as 8-bit unsigned integer
cgr_vulboth_out_path <- file.path(cgr_aligned_dir, "cgr_reproj_algn_vuln_both.tif")
writeRaster(cgr_algn_vulboth, cgr_vulboth_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# reclassify 
cgr_algn_encrch <- classify(cgr_algn, 
                             rcl = matrix(c(0, NA, # 
                                            5, 0, # Core grasslands
                                            100, 0, # Vulnerable - conversion
                                            111, 0, # Vulnerable - encroachment
                                            200, 0, # Vulnerable - both
                                            500, 0, # Plowprint
                                            550, 1, # Encroached
                                            800, 0, # Plowed/encroached
                                            1000, 0), # Masks
                                          ncol = 2, byrow = TRUE))
plot(cgr_algn_encrch)
cgr_algn_encrch

# Save the raster as 8-bit unsigned integer
cgr_encrch_out_path <- file.path(cgr_aligned_dir, "cgr_reproj_algn_encroached.tif")
writeRaster(cgr_algn_encrch, cgr_encrch_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# reclassify 
cgr_algn_plow <- classify(cgr_algn, 
                            rcl = matrix(c(0, NA, # 
                                           5, 0, # Core grasslands
                                           100, 0, # Vulnerable - conversion
                                           111, 0, # Vulnerable - encroachment
                                           200, 0, # Vulnerable - both
                                           500, 1, # Plowprint
                                           550, 0, # Encroached
                                           800, 0, # Plowed/encroached
                                           1000, 0), # Masks
                                         ncol = 2, byrow = TRUE))
plot(cgr_algn_plow)
cgr_algn_plow

# Save the raster as 8-bit unsigned integer
cgr_plow_out_path <- file.path(cgr_aligned_dir, "cgr_reproj_algn_plowed.tif")
writeRaster(cgr_algn_plow, cgr_plow_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# reclassify 
cgr_algn_plwencrch <- classify(cgr_algn, 
                          rcl = matrix(c(0, NA, # 
                                         5, 0, # Core grasslands
                                         100, 0, # Vulnerable - conversion
                                         111, 0, # Vulnerable - encroachment
                                         200, 0, # Vulnerable - both
                                         500, 0, # Plowprint
                                         550, 0, # Encroached
                                         800, 1, # Plowed/encroached
                                         1000, 0), # Masks
                                       ncol = 2, byrow = TRUE))
plot(cgr_algn_plwencrch)
cgr_algn_plwencrch

# Save the raster as 8-bit unsigned integer
cgr_plwencrch_out_path <- file.path(cgr_aligned_dir, "cgr_reproj_algn_plowed_encroached.tif")
writeRaster(cgr_algn_plwencrch, cgr_plwencrch_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

################

######### crop and mask to study area

# load raster
cgr_file <- file.path(cgr_aligned_dir, "cgr_reproj_algn.tif")
cgr <- rast(cgr_file)

crop_ext <- ext(sa)

# crop raster to study area
cgr_sa <- crop(cgr, crop_ext)
cgr_sa
plot (cgr_sa)

# mask to study area
cgr_masked <- mask(cgr_sa, sa)
cgr_masked
plot(cgr_masked)

# Save the raster as 8-bit unsigned integer
cgr_sa_out_path <- file.path(cgr_masked_dir, "cgr_reproj_algn_masked.tif")
writeRaster(cgr_masked, cgr_sa_out_path, datatype = "INT2U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 65535,
            overwrite = TRUE)

######### crop and mask to study area

# load binary raster
cgr_binary_file <- file.path(cgr_aligned_dir, "cgr_reproj_algn_binary.tif")
cgr_binary <- rast(cgr_binary_file)

crop_ext <- ext(sa)

# crop raster to study area
cgr_bin_sa <- crop(cgr_binary, crop_ext)
cgr_bin_sa
plot (cgr_bin_sa)

# mask to study area
cgr_bin_masked <- mask(cgr_bin_sa, sa)
cgr_bin_masked
plot(cgr_bin_masked)

# Save the raster as 8-bit unsigned integer
cgr_bin_sa_out_path <- file.path(cgr_masked_dir, "cgr_reproj_algn_binary_masked.tif")
writeRaster(cgr_bin_masked, cgr_bin_sa_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

######### crop and mask to study area

# load native vs tame raster
cgr_nattam_file <- file.path(cgr_aligned_dir, "cgr_reproj_algn_natvtam.tif")
cgr_nattam <- rast(cgr_nattam_file)

crop_ext <- ext(sa)

# crop raster to study area
cgr_nt_sa <- crop(cgr_nattam, crop_ext)
cgr_nt_sa
plot (cgr_nt_sa)

# mask to study area
cgr_nt_masked <- mask(cgr_nt_sa, sa)
cgr_nt_masked
plot(cgr_nt_masked)

# Save the raster as 8-bit unsigned integer
cgr_nt_sa_out_path <- file.path(cgr_masked_dir, "cgr_reproj_algn_natvtam_masked.tif")
writeRaster(cgr_nt_masked, cgr_nt_sa_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

######### crop and mask to study area

# load class raster
cgr_cg_file <- file.path(cgr_aligned_dir, "cgr_reproj_algn_core_grass.tif")
cgr_cg <- rast(cgr_cg_file)

crop_ext <- ext(sa)

# crop raster to study area
cgr_cg_sa <- crop(cgr_cg, crop_ext)
cgr_cg_sa
plot (cgr_cg_sa)

# mask to study area
cgr_cg_masked <- mask(cgr_cg_sa, sa)
cgr_cg_masked
plot(cgr_cg_masked)

# Save the raster as 8-bit unsigned integer
cgr_cg_sa_out_path <- file.path(cgr_masked_dir, "cgr_reproj_algn_core_grass_masked.tif")
writeRaster(cgr_cg_masked, cgr_cg_sa_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

######### crop and mask to study area

# load class raster
cgr_encr_file <- file.path(cgr_aligned_dir, "cgr_reproj_algn_encroached.tif")
cgr_encr <- rast(cgr_encr_file)

crop_ext <- ext(sa)

# crop raster to study area
cgr_encr_sa <- crop(cgr_encr, crop_ext)
cgr_encr_sa
plot (cgr_encr_sa)

# mask to study area
cgr_encr_masked <- mask(cgr_encr_sa, sa)
cgr_encr_masked
plot(cgr_encr_masked)

# Save the raster as 8-bit unsigned integer
cgr_encr_sa_out_path <- file.path(cgr_masked_dir, "cgr_reproj_algn_encroached_masked.tif")
writeRaster(cgr_encr_masked, cgr_encr_sa_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

######### crop and mask to study area

# load class raster
cgr_pl_file <- file.path(cgr_aligned_dir, "cgr_reproj_algn_plowed.tif")
cgr_pl <- rast(cgr_pl_file)

crop_ext <- ext(sa)

# crop raster to study area
cgr_pl_sa <- crop(cgr_pl, crop_ext)
cgr_pl_sa
plot (cgr_pl_sa)

# mask to study area
cgr_pl_masked <- mask(cgr_pl_sa, sa)
cgr_pl_masked
plot(cgr_pl_masked)

# Save the raster as 8-bit unsigned integer
cgr_pl_sa_out_path <- file.path(cgr_masked_dir, "cgr_reproj_algn_plowed_masked.tif")
writeRaster(cgr_pl_masked, cgr_pl_sa_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

######### crop and mask to study area

# load class raster
cgr_pe_file <- file.path(cgr_aligned_dir, "cgr_reproj_algn_plowed_encroached.tif")
cgr_pe <- rast(cgr_pe_file)

crop_ext <- ext(sa)

# crop raster to study area
cgr_pe_sa <- crop(cgr_pe, crop_ext)
cgr_pe_sa
plot(cgr_pe_sa)

# mask to study area
cgr_pe_masked <- mask(cgr_pe_sa, sa)
cgr_pe_masked
plot(cgr_pe_masked)

# Save the raster as 8-bit unsigned integer
cgr_pe_sa_out_path <- file.path(cgr_masked_dir, "cgr_reproj_algn_plowed_encroached_masked.tif")
writeRaster(cgr_pe_masked, cgr_pe_sa_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

######### crop and mask to study area

# load class raster
cgr_vb_file <- file.path(cgr_aligned_dir, "cgr_reproj_algn_vuln_both.tif")
cgr_vb <- rast(cgr_vb_file)

crop_ext <- ext(sa)

# crop raster to study area
cgr_vb_sa <- crop(cgr_vb, crop_ext)
cgr_vb_sa
plot(cgr_vb_sa)

# mask to study area
cgr_vb_masked <- mask(cgr_vb_sa, sa)
cgr_vb_masked
plot(cgr_vb_masked)

# Save the raster as 8-bit unsigned integer
cgr_vb_sa_out_path <- file.path(cgr_masked_dir, "cgr_reproj_algn_vuln_both_masked.tif")
writeRaster(cgr_vb_masked, cgr_vb_sa_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

######### crop and mask to study area

# load class raster
cgr_vc_file <- file.path(cgr_aligned_dir, "cgr_reproj_algn_vuln_conv.tif")
cgr_vc <- rast(cgr_vc_file)

crop_ext <- ext(sa)

# crop raster to study area
cgr_vc_sa <- crop(cgr_vc, crop_ext)
cgr_vc_sa
plot(cgr_vc_sa)

# mask to study area
cgr_vc_masked <- mask(cgr_vc_sa, sa)
cgr_vc_masked
plot(cgr_vc_masked)

# Save the raster as 8-bit unsigned integer
cgr_vc_sa_out_path <- file.path(cgr_masked_dir, "cgr_reproj_algn_vuln_conv_masked.tif")
writeRaster(cgr_vc_masked, cgr_vc_sa_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

######### crop and mask to study area

# load class raster
cgr_ve_file <- file.path(cgr_aligned_dir, "cgr_reproj_algn_vuln_encrch.tif")
cgr_ve <- rast(cgr_ve_file)

crop_ext <- ext(sa)

# crop raster to study area
cgr_ve_sa <- crop(cgr_ve, crop_ext)
cgr_ve_sa
plot(cgr_ve_sa)

# mask to study area
cgr_ve_masked <- mask(cgr_ve_sa, sa)
cgr_ve_masked
plot(cgr_ve_masked)

# Save the raster as 8-bit unsigned integer
cgr_ve_sa_out_path <- file.path(cgr_masked_dir, "cgr_reproj_algn_vuln_encrch_masked.tif")
writeRaster(cgr_ve_masked, cgr_ve_sa_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)
