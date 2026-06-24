
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
raw_grass_data_dir <- "R:/CanadianGrasslands/existingdatasets"
proc_grass_data_dir <- "R:/CanadianGrasslands/WG_Classification/processed"
out_grass_dir <- "./outputs/grass"

# ACI tiles directory
#aci_tile_root <- file.path(raw_grass_data_dir, "Cropland/AAFC_Annual_Crop_Inventory_2024/v2")
aci_tile_root <- file.path(raw_grass_data_dir, "Cropland/AAFC_Annual_Crop_Inventory_2025")
aci_merged_dir <- file.path(proc_grass_data_dir, "aafc_aci/merged")
aci_cropped_dir <- file.path(proc_grass_data_dir, "aafc_aci/cropped")
aci_aligned_dir <- file.path(proc_grass_data_dir, "aafc_aci/aligned")
aci_pilot_dir <- file.path(proc_grass_data_dir, "aafc_aci/pilot")

####--- Files ---####
# Pilot Area
pa_file <- file.path(pa_dir, "pilot_area_aligned_cnvc_Jan2026.tif")

# Study Area
sa_file <- file.path(sa_dir, "class_study_area_aligned_cnvc_Dec2025.tif")

# CNVC raster
#cnvc_file <- file.path(raw_grass_data_dir, "ns_cnvc/NatureServe_CNVC_COMPRESSED_20251023.tif")

# aci rasters
aci_bc_file <- file.path(aci_tile_root, "aci_2025_bc_v1/aci_2025_bc_v1.tif")
aci_ab_file <- file.path(aci_tile_root, "aci_2025_ab_v1/aci_2025_ab_v1.tif")
aci_sk_file <- file.path(aci_tile_root, "aci_2025_sk_v1/aci_2025_sk_v1.tif")
aci_mb_file <- file.path(aci_tile_root, "aci_2025_mb_v1/aci_2025_mb_v1.tif")
aci_on_file <- file.path(aci_tile_root, "aci_2025_on_v1/aci_2025_on_v1.tif")

##### use study area file as snap raster #####

# Load study area raster
sa <- rast(sa_file)

# print raster details 
origin(sa)
sa

####--- ALGIN RASTERS ---#### 

#----- AAFC Annual Crop Inventory ----#

# Load rasters
r_bc <- rast(aci_bc_file)
r_ab <- rast(aci_ab_file)
r_sk <- rast(aci_sk_file)
r_mb <- rast(aci_mb_file)
r_on <- rast(aci_on_file)

r_ab

# Set pixels with value 0 to NA - 0s are only on periphery
r_bc[r_bc == 0] <- NA
r_ab[r_ab == 0] <- NA
r_sk[r_sk == 0] <- NA
r_mb[r_mb == 0] <- NA
r_on[r_on == 0] <- NA

# Merge 
merged_aci <- mosaic(r_bc, r_ab, r_sk, r_mb, r_on, fun = "max") 
merged_aci_no0s <- merge(r_bc, r_ab, r_sk, r_mb, r_on) 

origin(merged_aci)  
merged_aci
merged_aci_no0s

# plot the result
plot(merged_aci, main = "Merged ACI raster")
plot(merged_aci_no0s, main = "Merged ACI raster - no zeroes")

# print frequency table
freq_table <- freq(merged_aci) 
print(freq_table) 

# Save the raster as 8-bit unsigned integer
aci_merge_out_path <- file.path(aci_merged_dir, "aafc_aci_merged.tif")
writeRaster(merged_aci, aci_merge_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# Save the raster as 8-bit unsigned integer
aci_merge_out_path <- file.path(aci_merged_dir, "aafc_aci_merged_no0s.tif")
writeRaster(merged_aci_no0s, aci_merge_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# Mark as categorical (factor)
merged_aci_no0s <- as.factor(merged_aci_no0s)

# reproject before resampling
aci_merged_reproj <- project(merged_aci_no0s, sa)
aci_merged_reproj

freq_table2 <- freq(aci_merged_reproj) 
print(freq_table2) 

# align to reference raster 
aci_merged_reproj_algn <- resample(aci_merged_reproj, sa, method = "near")
plot(aci_merged_reproj_algn)
aci_merged_reproj_algn

#class       : SpatRaster
#size        : 86539, 130223, 1  (nrow, ncol, nlyr)
#resolution  : 30, 30  (x, y)
#extent      : -2321003, 1585687, 270327.4, 2866497  (xmin, xmax, ymin, ymax)
#coord. ref. : North_America_Albers_Equal_Area_Conic (ESRI:102008)
#source      : spat_b64849645cad_46664_VuY9VpSXXzR4qRb.tif
#color table : 1
#varname     : class_study_area_aligned_cnvc_Dec2025
#categories  : aci_2025_bc_v1
#name        : aci_2025_bc_v1
#min value   :             20
#max value   :            230 

# Get values and counts
freq_table3 <- freq(aci_merged_reproj_algn) 
print(freq_table3) 

# Save the raster as 8-bit unsigned integer
aci_out_path <- file.path(aci_aligned_dir, "aafc_aci_2025_merged_reproj_algn.tif")
writeRaster(aci_merged_reproj_algn, aci_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# reclassify to binary (0 = not grass, 1 = grass)
aci_algn_binary <- classify(aci_merged_reproj_algn, 
                            rcl = matrix(c(
                              10, 0, # cloud
                              20, 0, # water
                              30, 0, # exposed/barren land
                              34, 0, # urban/developed
                              35, 0, # greenhouses
                              50, 0, # shrubland
                              60, 0, # forest fire/burnt area
                              80, 0, # wetland
                              85, 0, # peatland
                              110, 1, # grassland
                              120, 0, # agriculture (undifferentiated)
                              121, 0, # cropland
                              122, 1, # pasture/forages
                              130, 0, # too wet to be seeded
                              131, 0, # fallow
                              132, 0, # cereals
                              133, 0, # barley
                              134, 0, # other grains
                              135, 0, # millet
                              136, 0, # oats
                              137, 0, # rye
                              138, 0, # spelt
                              139, 0, # triticale
                              140, 0, # wheat
                              141, 0, # switchgrass
                              142, 0, # sorghum
                              143, 0, # quinoa
                              145, 0, # winter wheat
                              146, 0, # spring wheat
                              147, 0, # corn for grain
                              148, 0, # tobacco
                              149, 0, # ginseng
                              150, 0, # oilseeds
                              151, 0, # borage
                              152, 0, # camelina
                              153, 0, # canola/rapeseed
                              154, 0, # flaxseed
                              155, 0, # mustard
                              156, 0, # safflower
                              157, 0, # sunflower
                              158, 0, # soybeans
                              159, 0, # other oilseeds
                              160, 0, # pulses
                              161, 0, # other pulses
                              162, 0, # peas
                              163, 0, # chickpeas
                              167, 0, # beans
                              168, 0, # fababeans
                              174, 0, # lentils
                              175, 0, # vegetables
                              176, 0, # tomatoes
                              177, 0, # potatoes
                              178, 0, # sugarbeets
                              179, 0, # other vegetables
                              180, 0, # fruits
                              181, 0, # berries
                              182, 0, # blueberry
                              183, 0, # cranberry
                              185, 0, # other berries
                              188, 0, # orchards
                              189, 0, # other fruits
                              190, 0, # vineyards
                              191, 0, # hops
                              192, 0, # sod
                              193, 0, # herbs
                              194, 0, # nursery
                              195, 0, # buckwheat
                              196, 0, # canaryseed
                              197, 0, # hemp
                              198, 0, # vetch
                              199, 0, # other crops
                              200, 0, # forest (undifferentiated)
                              210, 0, # coniferous
                              220, 0, # broadleaf
                              230, 0), # mixedwood
                              ncol = 2, byrow = TRUE))

aci_algn_binary

#class       : SpatRaster
#size        : 86539, 130223, 1  (nrow, ncol, nlyr)
#resolution  : 30, 30  (x, y)
#extent      : -2321003, 1585687, 270327.4, 2866497  (xmin, xmax, ymin, ymax)
#coord. ref. : North_America_Albers_Equal_Area_Conic (ESRI:102008)
#source      : spat_b6486f0d2eb_46664_6Uxy02myoMrg5LO.tif
#varname     : class_study_area_aligned_cnvc_Dec2025
#name        : aci_2025_bc_v1
#min value   :              0
#max value   :              1

# Save the raster as 8-bit unsigned integer
aci_binary_out_path <- file.path(aci_aligned_dir, "aafc_aci_2025_merged_reproj_algn_binary_noshrub.tif")
writeRaster(aci_algn_binary, aci_binary_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# load ACI raster
#aci <- file.path(aci_aligned_dir, "aafc_aci_merged_reproj_algn.tif")
#aci_algn <- rast(aci)

# reclassify to native vs. tame (0 = not grass, 1 = tame, 2 = native)
aci_algn_nattam <- classify(aci_merged_reproj_algn, 
                            rcl = matrix(c(
                              10, 0, # cloud
                              20, 0, # water
                              30, 0, # exposed/barren land
                              34, 0, # urban/developed
                              35, 0, # greenhouses
                              50, 0, # shrubland
                              60, 0, # forest fire/burnt area
                              80, 0, # wetland
                              85, 0, # peatland
                              110, 2, # grassland
                              120, 0, # agriculture (undifferentiated)
                              121, 0, # cropland
                              122, 1, # pasture/forages
                              130, 0, # too wet to be seeded
                              131, 0, # fallow
                              132, 0, # cereals
                              133, 0, # barley
                              134, 0, # other grains
                              135, 0, # millet
                              136, 0, # oats
                              137, 0, # rye
                              138, 0, # spelt
                              139, 0, # triticale
                              140, 0, # wheat
                              141, 0, # switchgrass
                              142, 0, # sorghum
                              143, 0, # quinoa
                              145, 0, # winter wheat
                              146, 0, # spring wheat
                              147, 0, # corn for grain
                              148, 0, # tobacco
                              149, 0, # ginseng
                              150, 0, # oilseeds
                              151, 0, # borage
                              152, 0, # camelina
                              153, 0, # canola/rapeseed
                              154, 0, # flaxseed
                              155, 0, # mustard
                              156, 0, # safflower
                              157, 0, # sunflower
                              158, 0, # soybeans
                              159, 0, # other oilseeds
                              160, 0, # pulses
                              161, 0, # other pulses
                              162, 0, # peas
                              163, 0, # chickpeas
                              167, 0, # beans
                              168, 0, # fababeans
                              174, 0, # lentils
                              175, 0, # vegetables
                              176, 0, # tomatoes
                              177, 0, # potatoes
                              178, 0, # sugarbeets
                              179, 0, # other vegetables
                              180, 0, # fruits
                              181, 0, # berries
                              182, 0, # blueberry
                              183, 0, # cranberry
                              185, 0, # other berries
                              188, 0, # orchards
                              189, 0, # other fruits
                              190, 0, # vineyards
                              191, 0, # hops
                              192, 0, # sod
                              193, 0, # herbs
                              194, 0, # nursery
                              195, 0, # buckwheat
                              196, 0, # canaryseed
                              197, 0, # hemp
                              198, 0, # vetch
                              199, 0, # other crops
                              200, 0, # forest (undifferentiated)
                              210, 0, # coniferous
                              220, 0, # broadleaf
                              230, 0), # mixedwood
                              ncol = 2, byrow = TRUE))
plot(aci_algn_nattam)
aci_algn_nattam

# Save the raster as 8-bit unsigned integer
aci_nattam_out_path <- file.path(aci_aligned_dir, "aafc_aci_2025_merged_reproj_algn_natvtam_noshrub.tif")
writeRaster(aci_algn_nattam, aci_nattam_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# reclassify to isolate classes
aci_algn_grass <- classify(aci_merged_reproj_algn, 
                            rcl = matrix(c(
                              10, 0, # cloud
                              20, 0, # water
                              30, 0, # exposed/barren land
                              34, 0, # urban/developed
                              35, 0, # greenhouses
                              50, 0, # shrubland
                              60, 0, # forest fire/burnt area
                              80, 0, # wetland
                              85, 0, # peatland
                              110, 1, # grassland
                              120, 0, # agriculture (undifferentiated)
                              121, 0, # cropland
                              122, 0, # pasture/forages
                              130, 0, # too wet to be seeded
                              131, 0, # fallow
                              132, 0, # cereals
                              133, 0, # barley
                              134, 0, # other grains
                              135, 0, # millet
                              136, 0, # oats
                              137, 0, # rye
                              138, 0, # spelt
                              139, 0, # triticale
                              140, 0, # wheat
                              141, 0, # switchgrass
                              142, 0, # sorghum
                              143, 0, # quinoa
                              145, 0, # winter wheat
                              146, 0, # spring wheat
                              147, 0, # corn for grain
                              148, 0, # tobacco
                              149, 0, # ginseng
                              150, 0, # oilseeds
                              151, 0, # borage
                              152, 0, # camelina
                              153, 0, # canola/rapeseed
                              154, 0, # flaxseed
                              155, 0, # mustard
                              156, 0, # safflower
                              157, 0, # sunflower
                              158, 0, # soybeans
                              159, 0, # other oilseeds
                              160, 0, # pulses
                              161, 0, # other pulses
                              162, 0, # peas
                              163, 0, # chickpeas
                              167, 0, # beans
                              168, 0, # fababeans
                              174, 0, # lentils
                              175, 0, # vegetables
                              176, 0, # tomatoes
                              177, 0, # potatoes
                              178, 0, # sugarbeets
                              179, 0, # other vegetables
                              180, 0, # fruits
                              181, 0, # berries
                              182, 0, # blueberry
                              183, 0, # cranberry
                              185, 0, # other berries
                              188, 0, # orchards
                              189, 0, # other fruits
                              190, 0, # vineyards
                              191, 0, # hops
                              192, 0, # sod
                              193, 0, # herbs
                              194, 0, # nursery
                              195, 0, # buckwheat
                              196, 0, # canaryseed
                              197, 0, # hemp
                              198, 0, # vetch
                              199, 0, # other crops
                              200, 0, # forest (undifferentiated)
                              210, 0, # coniferous
                              220, 0, # broadleaf
                              230, 0), # mixedwood
                              ncol = 2, byrow = TRUE))
plot(aci_algn_grass)
aci_algn_grass

# Save the raster as 8-bit unsigned integer
aci_grass_out_path <- file.path(aci_aligned_dir, "aafc_aci_2025_merged_reproj_algn_grass.tif")
writeRaster(aci_algn_grass, aci_grass_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# reclassify to isolate classes
aci_algn_past <- classify(aci_merged_reproj_algn, 
                           rcl = matrix(c(
                             10, 0, # cloud
                             20, 0, # water
                             30, 0, # exposed/barren land
                             34, 0, # urban/developed
                             35, 0, # greenhouses
                             50, 0, # shrubland
                             60, 0, # forest fire/burnt area
                             80, 0, # wetland
                             85, 0, # peatland
                             110, 0, # grassland
                             120, 0, # agriculture (undifferentiated)
                             121, 0, # cropland
                             122, 1, # pasture/forages
                             130, 0, # too wet to be seeded
                             131, 0, # fallow
                             132, 0, # cereals
                             133, 0, # barley
                             134, 0, # other grains
                             135, 0, # millet
                             136, 0, # oats
                             137, 0, # rye
                             138, 0, # spelt
                             139, 0, # triticale
                             140, 0, # wheat
                             141, 0, # switchgrass
                             142, 0, # sorghum
                             143, 0, # quinoa
                             145, 0, # winter wheat
                             146, 0, # spring wheat
                             147, 0, # corn for grain
                             148, 0, # tobacco
                             149, 0, # ginseng
                             150, 0, # oilseeds
                             151, 0, # borage
                             152, 0, # camelina
                             153, 0, # canola/rapeseed
                             154, 0, # flaxseed
                             155, 0, # mustard
                             156, 0, # safflower
                             157, 0, # sunflower
                             158, 0, # soybeans
                             159, 0, # other oilseeds
                             160, 0, # pulses
                             161, 0, # other pulses
                             162, 0, # peas
                             163, 0, # chickpeas
                             167, 0, # beans
                             168, 0, # fababeans
                             174, 0, # lentils
                             175, 0, # vegetables
                             176, 0, # tomatoes
                             177, 0, # potatoes
                             178, 0, # sugarbeets
                             179, 0, # other vegetables
                             180, 0, # fruits
                             181, 0, # berries
                             182, 0, # blueberry
                             183, 0, # cranberry
                             185, 0, # other berries
                             188, 0, # orchards
                             189, 0, # other fruits
                             190, 0, # vineyards
                             191, 0, # hops
                             192, 0, # sod
                             193, 0, # herbs
                             194, 0, # nursery
                             195, 0, # buckwheat
                             196, 0, # canaryseed
                             197, 0, # hemp
                             198, 0, # vetch
                             199, 0, # other crops
                             200, 0, # forest (undifferentiated)
                             210, 0, # coniferous
                             220, 0, # broadleaf
                             230, 0), # mixedwood
                             ncol = 2, byrow = TRUE))
plot(aci_algn_past)
aci_algn_past

# Save the raster as 8-bit unsigned integer
aci_past_out_path <- file.path(aci_aligned_dir, "aafc_aci_2025_merged_reproj_algn_pastureforage.tif")
writeRaster(aci_algn_past, aci_past_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# reclassify to isolate classes
aci_algn_shrub <- classify(aci_merged_reproj_algn, 
                           rcl = matrix(c(
                             10, 0, # cloud
                             20, 0, # water
                             30, 0, # exposed/barren land
                             34, 0, # urban/developed
                             35, 0, # greenhouses
                             50, 1, # shrubland
                             60, 0, # forest fire/burnt area
                             80, 0, # wetland
                             85, 0, # peatland
                             110, 0, # grassland
                             120, 0, # agriculture (undifferentiated)
                             121, 0, # cropland
                             122, 0, # pasture/forages
                             130, 0, # too wet to be seeded
                             131, 0, # fallow
                             132, 0, # cereals
                             133, 0, # barley
                             134, 0, # other grains
                             135, 0, # millet
                             136, 0, # oats
                             137, 0, # rye
                             138, 0, # spelt
                             139, 0, # triticale
                             140, 0, # wheat
                             141, 0, # switchgrass
                             142, 0, # sorghum
                             143, 0, # quinoa
                             145, 0, # winter wheat
                             146, 0, # spring wheat
                             147, 0, # corn for grain
                             148, 0, # tobacco
                             149, 0, # ginseng
                             150, 0, # oilseeds
                             151, 0, # borage
                             152, 0, # camelina
                             153, 0, # canola/rapeseed
                             154, 0, # flaxseed
                             155, 0, # mustard
                             156, 0, # safflower
                             157, 0, # sunflower
                             158, 0, # soybeans
                             159, 0, # other oilseeds
                             160, 0, # pulses
                             161, 0, # other pulses
                             162, 0, # peas
                             163, 0, # chickpeas
                             167, 0, # beans
                             168, 0, # fababeans
                             174, 0, # lentils
                             175, 0, # vegetables
                             176, 0, # tomatoes
                             177, 0, # potatoes
                             178, 0, # sugarbeets
                             179, 0, # other vegetables
                             180, 0, # fruits
                             181, 0, # berries
                             182, 0, # blueberry
                             183, 0, # cranberry
                             185, 0, # other berries
                             188, 0, # orchards
                             189, 0, # other fruits
                             190, 0, # vineyards
                             191, 0, # hops
                             192, 0, # sod
                             193, 0, # herbs
                             194, 0, # nursery
                             195, 0, # buckwheat
                             196, 0, # canaryseed
                             197, 0, # hemp
                             198, 0, # vetch
                             199, 0, # other crops
                             200, 0, # forest (undifferentiated)
                             210, 0, # coniferous
                             220, 0, # broadleaf
                             230, 0), # mixedwood
                             ncol = 2, byrow = TRUE))
plot(aci_algn_shrub)
aci_algn_shrub

# Save the raster as 8-bit unsigned integer
aci_shrub_out_path <- file.path(aci_aligned_dir, "aafc_aci_2025_merged_reproj_algn_shrub.tif")
writeRaster(aci_algn_shrub, aci_shrub_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)


# Load study area raster
pa <- rast(pa_file)

# print raster details 
origin(pa)

# crop to pilot area
crop_ext <- ext(pa)
aci_pa <- crop(aci_merged_reproj_algn, crop_ext)
aci_pa
plot (aci_pa)

# mask to pilot area
aci_masked <- mask(aci_pa, pa)

# print raster details 
origin(aci_masked) # 
aci_masked
plot(aci_masked)

#class       : SpatRaster 
#dimensions  : 1367, 1692, 1  (nrow, ncol, nlyr)
#resolution  : 30, 30  (x, y)
#extent      : -828263.3, -777503.3, 1111437, 1152447  (xmin, xmax, ymin, ymax)
#coord. ref. : North_America_Albers_Equal_Area_Conic 
#source(s)   : memory
#color table : 1 
#varname     : class_study_area_aligned_cnvc_Dec2025 
#categories  : aci_2024_bc_v2 
#name        : aci_2024_bc_v2 
#min value   :             20 
#max value   :            220 

# Save the raster as 8-bit unsigned integer
aci_pa_out_path <- file.path(aci_pilot_dir, "aci_merged_reproj_algn_crop_pa.tif")
writeRaster(aci_masked, aci_pa_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# crop binary raster to pilot area
aci_bin_pa <- crop(aci_algn_binary, crop_ext)
aci_bin_pa
plot (aci_bin_pa)

# mask to pilot area
aci_bin_masked <- mask(aci_bin_pa, pa)
aci_bin_masked
plot(aci_bin_masked)

# Save the raster as 8-bit unsigned integer
aci_bin_pa_out_path <- file.path(aci_pilot_dir, "aci_merged_reproj_algn_crop_pa_binary_noshrub.tif")
writeRaster(aci_bin_masked, aci_bin_pa_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# crop native vs. tame raster to pilot area
aci_nt_pa <- crop(aci_algn_nattam, crop_ext)
aci_nt_pa
plot (aci_nt_pa)

# mask to pilot area
aci_nt_masked <- mask(aci_nt_pa, pa)
aci_nt_masked
plot(aci_nt_masked)

# Save the raster as 8-bit unsigned integer
aci_nt_pa_out_path <- file.path(aci_pilot_dir, "aci_merged_reproj_algn_crop_pa_natvtam_noshrub.tif")
writeRaster(aci_nt_masked, aci_nt_pa_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)
