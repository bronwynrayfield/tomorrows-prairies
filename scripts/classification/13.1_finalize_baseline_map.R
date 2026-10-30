
###################################################################################
##                    FINALIZE BASELINE MAP FOR MODELLING                        ##
###################################################################################

##----------------------- Alison Long - March 2026 -------------------------

# STEP 13: Fill in non-grassland pixels of the consensus map with other land cover
#          classes to generate baseline map for landscape modelling. Modification 5
#          was chosen as the final consensus map.

# load packages
library(terra)
library(sf)
library(tidyverse)

# Set the tempocmry directory to a drive with sufficient space
terraOptions(tempdir = "R:/CanadianGrasslands/WG_classification/scratch")

####--- DIRECTORIES ---####

# pilot area
#pa_dir <- "R:/CanadianGrasslands/processeddatasets/pilotarea"

# study area
sa_dir <- "R:/CanadianGrasslands/processeddatasets/studyarea"

# ecodistricts
eco_dir <- "R:/CanadianGrasslands/processeddatasets/boundaries/national-ecological-framework-for-canada/Ecodistricts"

# grasslands
raw_grass_data_dir  <- "R:/CanadianGrasslands/existingdatasets"
proc_grass_data_dir <- "R:/CanadianGrasslands/WG_Classification/processed"
consens_data_dir    <- "R:/CanadianGrasslands/WG_Classification/results/consensus_mapping"
out_dir             <- "R:/CanadianGrasslands/WG_Classification/results/consensus_mapping"

# LULC raster directory
luts_aligned_dir  <- file.path(proc_grass_data_dir, "aafc_luts_2020_30m_v5/aligned")
aci_aligned_dir   <- file.path(proc_grass_data_dir, "aafc_aci/aligned")

####--- FILES ---####

# Pilot area
#pa_file <- file.path(pa_dir, "pilot_area_aligned_cnvc_Jan2026.tif")

# Study area
sa_file <- file.path(sa_dir, "modelling_extent_aligned_cnvc_Jan2026.tif")

# ecodistricts
eco_file <- file.path(eco_dir, "ecodistricts_clean_class_sa.shp")

# grass raster files
luts_file <- file.path(luts_aligned_dir, "aafc_luts_merged_reproj_algn.tif")
aci_file  <- file.path(aci_aligned_dir, "aafc_aci_2025_merged_reproj_algn.tif")

# consensus map file
#cons_file <- file.path(consens_data_dir, "v4_temporal_modelling_extent_24Jul2026/grass_3class.tif")
#cons_file <- file.path(consens_data_dir, "new_fusion_global_norm_modext_26Aug2026/hard_class.tif")
#cons_file <- file.path(consens_data_dir, "new_fusion_global_norm_modext_classdifficultynorm_31Aug2026/hard_class.tif")
cons_file <- file.path(consens_data_dir, "new_fusion_global_norm_modext_classdifficultynorm_31Aug2026/hard_class_m5_boost_tame_grass.tif")

####--- PREP LC RASTER ---####

# Read rasters
sa    <- rast(sa_file)
luts  <- rast(luts_file)
aci   <- rast(aci_file)
cm    <- rast(cons_file)

luts
aci
cm

plot(luts)
plot(aci)
plot(cm)

####--- MASK TO MODELLING EXTENT ---####

# Mask to the extent of sa
# `mask()` keeps only cells where sa is not NA
# `crop()` ensures the output extent exactly matches sa
luts_masked <- mask(crop(luts, sa), sa)
aci_masked  <- mask(crop(aci, sa), sa)
cm_masked   <- mask(crop(cm, sa), sa)

# consensus map has the northern-most tip cut off
# will crop/mask to it instead of sa
#aci_masked  <- mask(crop(aci, cm), cm)
#plot(aci_masked)

# Get values and counts 
freq_table <- freq(aci_masked) 
print(freq_table)

####--- RECLASSIFY LULC RASTER ---####

# reclassify LUTS to modelling classes 
# 8 = tame grass
# 9 = native grass
# 4 = settlement
# 3 = annual crop
# 5 = other 
luts_masked_reclass <- classify(luts_masked, 
                             rcl = matrix(c(21, 4, # Settlement
                                            22, 4, # High Reflectance Settlement
                                            24, 4, # Settlement Forest
                                            25, 5, # Roads
                                            28, 4, # Vegetated Settlement
                                            29, 4, # Very High Reflectance Settlement
                                            31, 5, # Water
                                            41, 5, # Forest
                                            42, 5, # Forest Wetland
                                            43, 5, # Forest Regenerating after Harvest (<20 years)
                                            44, 5, # Forest Wetland Regenerating after Harvest (<20 years)
                                            47, 5, # Forest Regenerating after Harvest (20-29 years)
                                            48, 5, # Forest Wetland Regenerating after Harvest (20-29 years)
                                            49, 5, # Forest Regenerating after Fire (<20 years)
                                            51, 8, # Cropland
                                            52, 3, # Annual Cropland
                                            55, 8, # Land Converted to Cropland
                                            56, 3, # Land Converted to Annual Cropland
                                            63, 9, # Prairie Grassland
                                            64, 9, # Temperate Grassland
                                            65, 9, # Temperate Rainforest Grassland
                                            66, 9, # Boreal Grassland
                                            67, 9, # Montane Grassland
                                            68, 9, # Alpine Tundra Grassland
                                            69, 9, # Arctic Tundra Grassland
                                            73, 5, # Prairie Wetland
                                            74, 5, # Temperate Wetland
                                            75, 5, # Temperate Rainforest Wetland
                                            76, 5, # Boreal Wetland
                                            78, 5, # Alpine Tundra Wetland
                                            79, 5, # Arctic Tundra Wetland
                                            81, 4, # Newly-Detected Settlement (<10 years)
                                            82, 4, # Newly-Detected High Reflectance Settlement (<10 years)
                                            84, 4, # Newly-Detected Settlement Forest (<10 years)
                                            85, 5, # Newly-detected roads (<10 years)
                                            88, 4, # Newly-Detected Vegetated Settlement (<10 years)
                                            89, 4, # Newly-Detected Very High Reflectance Settlement (<10 years)
                                            91, 5, # Other Land
                                            92, 5), # Snow and Ice
                                          ncol = 2, byrow = TRUE))

luts_masked_reclass
plot(luts_masked_reclass)

# reclassify ACI to modelling classes
# 8 = tame grass (from ACI not consensus map)
# 9 = native grass (from ACI not consensus map)
# 4 = settlement
# 3 = annual crop
# 5 = other 
aci_masked_reclass <- classify(aci_masked, 
                            rcl = matrix(c(
                              10, 5, # cloud
                              20, 5, # water
                              30, 5, # exposed/barren land
                              34, 4, # urban/developed
                              35, 5, # greenhouses
                              50, 5, # shrubland
                              60, 5, # forest fire/burnt area
                              80, 5, # wetland
                              85, 5, # peatland
                              110, 9, # grassland
                              120, 3, # agriculture (undifferentiated)
                              121, 3, # cropland
                              122, 8, # pasture/forages
                              130, 5, # too wet to be seeded
                              131, 3, # fallow
                              132, 3, # cereals
                              133, 3, # barley
                              134, 3, # other grains
                              135, 3, # millet
                              136, 3, # oats
                              137, 3, # rye
                              138, 3, # spelt
                              139, 3, # triticale
                              140, 3, # wheat
                              141, 3, # switchgrass
                              142, 3, # sorghum
                              143, 3, # quinoa
                              145, 3, # winter wheat
                              146, 3, # spring wheat
                              147, 3, # corn for grain
                              148, 3, # tobacco
                              149, 3, # ginseng
                              150, 3, # oilseeds
                              151, 3, # borage
                              152, 3, # camelina
                              153, 3, # canola/rapeseed
                              154, 3, # flaxseed
                              155, 3, # mustard
                              156, 3, # safflower
                              157, 3, # sunflower
                              158, 3, # soybeans
                              159, 3, # other oilseeds
                              160, 3, # pulses
                              161, 3, # other pulses
                              162, 3, # peas
                              163, 3, # chickpeas
                              167, 3, # beans
                              168, 3, # fababeans
                              174, 3, # lentils
                              175, 3, # vegetables
                              176, 3, # tomatoes
                              177, 3, # potatoes
                              178, 3, # sugarbeets
                              179, 3, # other vegetables
                              180, 3, # fruits
                              181, 3, # berries
                              182, 3, # blueberry
                              183, 3, # cranberry
                              185, 3, # other berries
                              188, 5, # orchards
                              189, 3, # other fruits
                              190, 5, # vineyards
                              191, 3, # hops
                              192, 8, # sod - only 700 pixels
                              193, 3, # herbs
                              194, 5, # nursery
                              195, 3, # buckwheat
                              196, 3, # canaryseed
                              197, 3, # hemp
                              198, 3, # vetch
                              199, 3, # other crops
                              200, 5, # forest (undifferentiated)
                              210, 5, # coniferous
                              220, 5, # broadleaf
                              230, 5), # mixedwood
                              ncol = 2, byrow = TRUE))

aci_masked_reclass
plot(aci_masked_reclass)

####--- RECLASSIFY GRASS RASTER ---####

# Reclassify 0 values (non-grass) to NA
# 1 --> NA (nongrass)
# 2 --> 2  (tame)
# 3 --> 1  (native)
cm_masked_reclass <- classify(cm_masked, 
                                rcl = matrix(c(3, 1,  # Native grass
                                               2, 2,  # Tame grass
                                               1, NA), # Non-grass
                                             ncol = 2, byrow = TRUE))

cm_masked_reclass
plot(cm_masked_reclass)

####--- FILL NA PIXELS ---####

# Validate the rasters align 
stopifnot(
  all.equal(ext(cm_masked_reclass), ext(luts_reclass)),
  all.equal(res(cm_masked_reclass), res(luts_reclass)),
  crs(cm_masked_reclass) == crs(luts_reclass)
)

stopifnot(
  all.equal(ext(cm_masked_reclass), ext(aci_masked_reclass)),
  all.equal(res(cm_masked_reclass), res(aci_masked_reclass)),
  crs(cm_masked_reclass) == crs(aci_masked_reclass)
)

# Fill NA pixels in A with values from B
# cover(A, B) returns A where not NA, otherwise B
cm_filled <- cover(cm_masked_reclass, luts_masked_reclass)
plot(cm_filled)

cm_filled <- cover(cm_masked_reclass, aci_masked_reclass)
plot(cm_filled)

# Get values and counts 
freq_table <- freq(cm_filled) 
print(freq_table)

# LUTS version (consensus map 1.0 - v4)
#   layer value     count
# 1     1     1  97310792
# 2     1     2  64075412
# 3     1     3 272870339
# 4     1     4  11026463
# 5     1     5 177284642
# 6     1     8  55260431
# 7     1     9   7667540

# ACI version (consensus map 1.0 - v4)
#   layer value     count
# 1     1     1  97310792
# 2     1     2  64075412
# 3     1     3 321304392
# 4     1     4  16451855
# 5     1     5 160122868
# 6     1     8  15343363
# 7     1     9  10890541

# ACI version (consensus map 2.0 - new data fusion global norm)
#   layer value     count
# 1     1     1  98148288
# 2     1     2  53304622
# 3     1     3 327433111
# 4     1     4  17522116
# 5     1     5 162451540
# 6     1     8  14866664
# 7     1     9  10769043

# ACI version (consensus map 2.0 - new data fusion global norm with class difficulty norm)
#   layer value       count
# 1     1     1  97,916,599 - native
# 2     1     2  63,179,602 - tame
# 3     1     3 324,903,505 - annual crop
# 4     1     4  17,421,116 - settlement
# 5     1     5 161,288,034 - other
# 6     1     8  11,385,049 - ACI pasture/forages (tame)
# 7     1     9   9,418,767 - ACI grassland (native)

# LUTS version (consensus map 2.0 - new data fusion global norm with class difficulty norm)
#   layer value       count
# 1     1     1  97,916,599 - native
# 2     1     2  63,179,602 - tame
# 3     1     3 274,325,139 - annual crop
# 4     1     4  11,414,568 - settlement
# 5     1     5 178,582,139 - other
# 6     1     8  52,166,366 - LUTS cropland & land converted to cropland (tame)
# 7     1     9   7,923,916 - LUTS native

# ACI version (consensus map 2.0 - new data fusion global norm with class difficulty norm - mod5)
#   layer value       count
# 1     1     1  93,487,558
# 2     1     2  75,760,518
# 3     1     3 321,673,550
# 4     1     4  17,127,974
# 5     1     5 159,984,616
# 6     1     8   9,222,987
# 7     1     9   8,258,481

# what to do with pixels that WERE NOT classified as grass
# in consensus map but WERE in lulc dataset (classes 8 and 9)?

#	ACI says forage: assume that it is annual crop
#	ACI says native: could go to OTHER or NATIVE in filled consensus map

cm_filled_reclass <- classify(cm_filled, 
                              rcl = matrix(c(1, 1,  # Native grass
                                             2, 2,  # Tame grass
                                             3, 3,  # Annual crop
                                             4, 4,  # Settlement
                                             5, 5,  # Other
                                             8, 3,  # ACI pasture/forage --> annual crop
                                             9, 1), # ACI grassland --> native grass
                                           ncol = 2, byrow = TRUE))

cm_filled_reclass
plot(cm_filled_reclass)

# Get values and counts
freq_table <- freq(cm_filled_reclass) 
print(freq_table)

# LUTS version (consensus map 1.0 - v4)
#   layer value     count
# 1     1     1 104978332
# 2     1     2 119335843
# 3     1     3 272870339
# 4     1     4  11026463
# 5     1     5 177284642

# ACI version (consensus map 1.0 - v4)
#   layer value     count
# 1     1     1 108201333
# 2     1     2  79418775
# 3     1     3 321304392
# 4     1     4  16451855
# 5     1     5 160122868

####--- EXPORT ---####

# Write result efficiently for large rasters (chunked writing)
writeRaster(cm_filled_reclass,
            filename = file.path(out_dir, "TP_baseline_map_3Sep2026.tif"),
            overwrite = TRUE,
            wopt = list(datatype = "INT1U", gdal = c("COMPRESS=LZW", "TILED=YES")))


#----- ACI reclass for comparison to consensus map -----#

# reclassify ACI to modelling classes
# 2 = tame grass
# 1 = native grass
# 4 = settlement
# 3 = annual crop
# 5 = other 
aci_masked_reclass_compare <- classify(aci_masked, 
                               rcl = matrix(c(
                                 10, 5, # cloud
                                 20, 5, # water
                                 30, 5, # exposed/barren land
                                 34, 4, # urban/developed
                                 35, 5, # greenhouses
                                 50, 5, # shrubland
                                 60, 5, # forest fire/burnt area
                                 80, 5, # wetland
                                 85, 5, # peatland
                                 110, 1, # grassland
                                 120, 3, # agriculture (undifferentiated)
                                 121, 3, # cropland
                                 122, 2, # pasture/forages
                                 130, 5, # too wet to be seeded
                                 131, 3, # fallow
                                 132, 3, # cereals
                                 133, 3, # barley
                                 134, 3, # other grains
                                 135, 3, # millet
                                 136, 3, # oats
                                 137, 3, # rye
                                 138, 3, # spelt
                                 139, 3, # triticale
                                 140, 3, # wheat
                                 141, 3, # switchgrass
                                 142, 3, # sorghum
                                 143, 3, # quinoa
                                 145, 3, # winter wheat
                                 146, 3, # spring wheat
                                 147, 3, # corn for grain
                                 148, 3, # tobacco
                                 149, 3, # ginseng
                                 150, 3, # oilseeds
                                 151, 3, # borage
                                 152, 3, # camelina
                                 153, 3, # canola/rapeseed
                                 154, 3, # flaxseed
                                 155, 3, # mustard
                                 156, 3, # safflower
                                 157, 3, # sunflower
                                 158, 3, # soybeans
                                 159, 3, # other oilseeds
                                 160, 3, # pulses
                                 161, 3, # other pulses
                                 162, 3, # peas
                                 163, 3, # chickpeas
                                 167, 3, # beans
                                 168, 3, # fababeans
                                 174, 3, # lentils
                                 175, 3, # vegetables
                                 176, 3, # tomatoes
                                 177, 3, # potatoes
                                 178, 3, # sugarbeets
                                 179, 3, # other vegetables
                                 180, 3, # fruits
                                 181, 3, # berries
                                 182, 3, # blueberry
                                 183, 3, # cranberry
                                 185, 3, # other berries
                                 188, 5, # orchards
                                 189, 3, # other fruits
                                 190, 5, # vineyards
                                 191, 3, # hops
                                 192, 2, # sod
                                 193, 3, # herbs
                                 194, 5, # nursery
                                 195, 3, # buckwheat
                                 196, 3, # canaryseed
                                 197, 3, # hemp
                                 198, 3, # vetch
                                 199, 3, # other crops
                                 200, 5, # forest (undifferentiated)
                                 210, 5, # coniferous
                                 220, 5, # broadleaf
                                 230, 5), # mixedwood
                                 ncol = 2, byrow = TRUE))

aci_masked_reclass_compare
plot(aci_masked_reclass_compare)

# Write result efficiently for large rasters 
writeRaster(aci_masked_reclass_compare,
            filename = file.path(out_dir, "ACI_2025_reclass_for_baseline.tif"),
            overwrite = TRUE,
            wopt = list(datatype = "INT1U", gdal = c("COMPRESS=LZW", "TILED=YES")))
