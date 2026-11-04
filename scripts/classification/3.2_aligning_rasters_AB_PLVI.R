
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

# PLVI directory
plvi_root <- file.path(raw_grass_data_dir, "Landcover/AB Primary Land and Vegetation Inventory/PLVI_SHP/Data")
plvi_aligned_dir <- file.path(proc_grass_data_dir, "ab_plvi/aligned")
plvi_cropped_dir <- file.path(proc_grass_data_dir, "ab_plvi/cropped")
plvi_rast_dir <- file.path(proc_grass_data_dir, "ab_plvi/rasterized")
plvi_pilot_dir <- file.path(proc_grass_data_dir, "ab_plvi/pilot")

####--- Files ---####
# Pilot area
pa_file <- file.path(pa_dir, "pilot_area_aligned_cnvc_Jan2026.tif")

# Study area
sa_file <- file.path(sa_dir, "class_study_area_aligned_cnvc_Dec2025.tif")

# CNVC raster
#cnvc_file <- file.path(raw_grass_data_dir, "ns_cnvc/NatureServe_CNVC_COMPRESSED_20251023.tif")

# plvi raster
plvi_file <- file.path(plvi_root, "PrimaryLandAndVegetationInventory.shp")

##### use study area file as snap raster #####

# Load study area raster
sa <- rast(sa_file)

# print raster details 
origin(sa) # 6.696121 -2.550334
sa

####--- ALGIN RASTERS ---#### 

#----- AB Primary Land and Vegetation Inventory ----#

# Load rasters
v_plvi <- vect(plvi_file)
v_plvi

#class       : SpatVector 
#geometry    : polygons 
#dimensions  : 179087, 58  (geometries, attributes)
#extent      : 181310.5, 842339.2, 5756870, 6441873  (xmin, xmax, ymin, ymax)
#source      : PrimaryLandAndVegetationInventory.shp
#coord. ref. : NAD83 / Alberta 10-TM (Forest) (EPSG:3400) 
#names       : POLYGON_ID   NSR PRIME_CLS1 LAND_CLS1 SITE_TYP1 SITE_PCT1 RNG_SITE1 MAP_CODE1 ECO_PH1_ID
#type        :      <num> <chr>      <chr>     <chr>     <chr>     <num>     <chr>     <chr>      <num>
#values      :        212    CM        VEG       WET        TF         8        NA        9C       3071
#213    CM        VEG       NAW        FT         8        NA        5C       3030
#214    CM        VEG       NAW        FT         8        NA        5C       3030
#ECO_PHASE1 (and 48 more)
#<chr>              
#  j1              
#d1              
#d1  

# Filter plvi to include only records with 2015 and later
v_plvi_gte2015 <- v_plvi[v_plvi$IMAGE_YR >= 2015, ]
v_plvi_gte2015

# 57594, 58  (geometries, attributes)

##---- Rasterize vector data ----##

# ---- 1) Ensure both layers use the same CRS as 'sa' ----
if (!same.crs(v_plvi, sa)) {
  v_plvi <- project(v_plvi, sa)  # reproject polygons to match 'sa'
}

# ---- 2) Build a categorical code for SITE_TYP1 ----
# terra rasters can't store character values; we map each category to an integer
cats  <- sort(unique(v_plvi$SITE_TYP1))
codes <- seq_along(cats)
# attach integer code to the vector attributes
v_plvi$site_code <- match(v_plvi$SITE_TYP1, cats)

# ---- 3) Create an empty raster aligned to 'sa' ----
r_tmpl <- rast(sa)  # inherits extent, resolution, origin, and CRS from 'sa'

# ---- 4) Rasterize polygons using the integer code ----
# For polygons, 'touches=TRUE' helps fill edge pixels; you can set it to FALSE if undesired.
r_code <- rasterize(v_plvi, y = r_tmpl,
                    field     = "site_code",
                    touches   = TRUE,
                    background = NA)
                    

# (Optional) If 'sa' has NA areas you want to enforce, mask to 'sa'’s valid footprint:
# r_code <- mask(r_code, sa)

# ---- 5) Make the raster categorical (factor) and attach the lookup table (RAT) ----
r_cat <- as.factor(r_code)
# Set the level table so that cell values map back to the two-letter codes
levels(r_cat) <- data.frame(ID = codes, SITE_TYP1 = cats)

# ---- 6) Write to disk ----
# Save the raster as 16-bit unsigned integer
plvi_out_path <- file.path(plvi_aligned_dir, "plvi_rast_reproj_algn.tif")
writeRaster(r_cat, plvi_out_path, datatype = "INT2U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 65535,
            overwrite = TRUE)

# ---- 7) Quick check ----
levels(r_cat)       # view the RAT
plot(r_cat)         # quick look (legend will show SITE_TYP1 codes)
r_cat

#class       : SpatRaster 
#size        : 86539, 130223, 1  (nrow, ncol, nlyr)
#resolution  : 30, 30  (x, y)
#extent      : -2321003, 1585687, 270327.4, 2866497  (xmin, xmax, ymin, ymax)
#coord. ref. : North_America_Albers_Equal_Area_Conic 
#source      : spat_3edb4548973b2_257460_o8up7yuDdJ6XHBX.tif 
#varname     : class_study_area_aligned_cnvc_Dec2025 
#categories  : SITE_TYP1 
#name        : SITE_TYP1 
#min value   :       AIH 
#max value   :        WT 

# print frequency table
freq_table <- freq(r_cat) 
print(freq_table)

#   layer value    count
#1      1   AIH   612445 # Transportation Surface (AIH)
#2      1   AII   234748 # Industrial Tracts (AII)
#3      1   ASC   421150 # Settlement Tracts (ASC)
#4      1    CA 31186109 # Annual Crops (CA)
#5      1    CC    50999 # Recent Cut (CC)
#6      1   CIP   115656 # Non-treed Clearing (CIP)
#7      1   CIT    31722 # Treed Clearing (CIT)
#8      1    CP 10918545 # Tame Pasture (CP)
#9      1   CPR  3297379 # Rough Pasture (CPR)
#10     1    CR   596359 # Regenerated Cut Block (CR)
#11     1    CS   530903 # Closed Shrub (CS)
#12     1    FT 19278896 # Forested (FT)
#13     1    GF  1359202 # Grass Fen (GF)
#14     1    HG  2762417 # Herbaceous Grass (HG)
#15     1     M   527983 # Marsh (M)
#16     1    MS   431765 # Medial Shrub (MS)
#17     1   NMR    41100 # Rock (NMR)
#18     1   NMS    13113 # Sand (NMS)
#19     1    NW  3163761 # Water (NW)
#20     1   NWI      167 # Ice (NWI)
#21     1    OS  1631040 # Open Shrub (OS)
#22     1    PF    46041 # Patterned Fen (PF)
#23     1    SF  1012655 # Shrub Fen (SF)
#24     1    SW     5111 # Swamp (SW)
#25     1    TF  2457209 # Treed Fen (TF)
#26     1    WS    31592 # Shrub (WS)
#27     1    WT   387090 # Treed (WT)

##---- Rasterize filtered vector data ----##

# ---- 1) Ensure both layers use the same CRS as 'sa' ----
if (!same.crs(v_plvi_gte2015, sa)) {
  v_plvi_gte2015 <- project(v_plvi_gte2015, sa)  # reproject polygons to match 'sa'
}

# ---- 2) Build a categorical code for SITE_TYP1 ----
# terra rasters can't store character values; we map each category to an integer
cats_filt  <- sort(unique(v_plvi_gte2015$SITE_TYP1))
codes_filt <- seq_along(cats_filt)
# attach integer code to the vector attributes
v_plvi_gte2015$site_code <- match(v_plvi_gte2015$SITE_TYP1, cats_filt)

# ---- 3) Create an empty raster aligned to 'sa' ----
r_tmpl_filt <- rast(sa)  # inherits extent, resolution, origin, and CRS from 'sa'

# ---- 4) Rasterize polygons using the integer code ----
# For polygons, 'touches=TRUE' helps fill edge pixels; you can set it to FALSE if undesired.
r_code_filt <- rasterize(v_plvi_gte2015, y = r_tmpl_filt,
                    field     = "site_code",
                    touches   = TRUE,
                    background = NA)

plot(r_code_filt)

# ---- 5) Make the raster categorical (factor) and attach the lookup table (RAT) ----
r_cat_filt <- as.factor(r_code_filt)
# Set the level table so that cell values map back to the two-letter codes
levels(r_cat_filt) <- data.frame(ID = codes_filt, SITE_TYP1 = cats_filt)

# ---- 6) Write to disk ----
# Save the raster as 16-bit unsigned integer
plvi_out_path <- file.path(plvi_aligned_dir, "plvi_rast_reproj_algn_gte2015.tif")
writeRaster(r_cat_filt, plvi_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

## reclassify to binary (1 for grass)

# 1) Confirm the levels table and the column name that holds the labels
levels(r_cat)[[1]] |> head()
# Should see a column named "SITE_TYP1" with values like "AIH", "CP", "HG", etc.

# 2) Define which labels map to 1 (grass) and all others to 0
grass_codes <- c("CP", "CPR", "HG")  # Tame Pasture, Rough Pasture, Herbaceous Grass

# 3) Build a lookup table (label -> 0/1)
# Pull all labels that actually occur in your raster from the freq table (or from levels)
#labels_present <- unique(freq(r_cat)$value)

#lut_bin <- data.frame(
#  SITE_TYP1  = labels_present,
#  grass_bin  = ifelse(labels_present %in% grass_codes, 1, 0)
#)

# 4) Substitute by label column ("SITE_TYP1"), returning 0/1
#plvi_algn_binary <- terra::subs(
#  r_cat,
#  from = lut_bin,
#  by   = "SITE_TYP1",   # column in the raster's levels table that holds the labels
#  which= "grass_bin",   # the column from `lut` we want as output
#  subsWithNA = TRUE     # make unmatched labels NA; here all labels are covered
#)

# Binary raster: 1 for grass classes, 0 otherwise (except NA)
plvi_algn_binary_filt <- ifel(
  is.na(r_cat_filt),                     # if the input cell is NA
  NA,                               # keep NA
  ifel(r_cat_filt %in% grass_codes, 1, 0)  # otherwise 1 for grass, 0 for non-grass
)

# 5) Quick checks
plvi_algn_binary
freq(plvi_algn_binary)
plot(plvi_algn_binary) 

# Save the raster as 8-bit unsigned integer
plvi_binary_out_path <- file.path(plvi_aligned_dir, "plvi_reproj_algn_binary_gte2015.tif")
writeRaster(plvi_algn_binary_filt, plvi_binary_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# load raster
plvi <- file.path(plvi_aligned_dir, "plvi_rast_reproj_algn.tif")
plvi_algn <- rast(plvi)
plot(plvi_algn)
plvi_algn

#### reclassify to native vs. tame (0 = not grass, 1 = tame, 2 = native)
# Define mappings
tame_codes <- c("CP", "CPR")   # -> 1
native_codes <- c("HG")          # -> 2

plvi_algn_nt_filt <- ifel(
  is.na(r_cat_filt), 
  NA,
  ifel(r_cat_filt %in% native_codes, 2,    # HG -> 2
       ifel(r_cat_filt %in% tame_codes, 1,  # CP, CPR -> 1
            0)                           # all others -> 0
  )
)

plot(plvi_algn_nt_filt)

# Save the raster as 8-bit unsigned integer
plvi_nt_out_path <- file.path(plvi_aligned_dir, "plvi_reproj_algn_natvtam_gte2015.tif")
writeRaster(plvi_algn_nt_filt, plvi_nt_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# load raster
plvi <- file.path(plvi_aligned_dir, "plvi_rast_reproj_algn.tif")
plvi_algn <- rast(plvi)
plot(plvi_algn)
plvi_algn

#### reclassify to native vs. tame (0 = not grass, 1 = tame, 2 = native)

# 1) Build lookup table
lut_nt <- data.frame(
  SITE_TYP1 = labels_present,
  recode    = ifelse(labels_present == "HG", 2,          # HG -> 2
                     ifelse(labels_present %in% c("CP", "CPR"), 1,  # CP, CPR -> 1
                            0))                                  # all others -> 0
)

# 2) Substitute using SITE_TYP1
plvi_algn_nattam <- terra::subs(
  r_cat,
  from       = lut_nt,
  by         = "SITE_TYP1",
  which      = "recode",
  subsWithNA = TRUE   # unmatched labels become NA (but all are covered here)
)

## another way

# Define mappings
tame_codes <- c("CP", "CPR")   # -> 1
native_codes <- c("HG")          # -> 2

plvi_algn_nattam <- ifel(
  is.na(plvi_algn), 
  NA,
  ifel(plvi_algn %in% native_codes, 2,    # HG -> 2
       ifel(plvi_algn %in% tame_codes, 1,  # CP, CPR -> 1
            0)                           # all others -> 0
  )
)

plot(plvi_algn_nattam)
plvi_algn_nattam

# Save the raster as 8-bit unsigned integer
plvi_nattam_out_path <- file.path(plvi_aligned_dir, "plvi_reproj_algn_natvtam.tif")
writeRaster(plvi_algn_nattam, plvi_nattam_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

## not present in pilot area!!

################################################

# isolate shrub classes

# load raster
#plvi <- file.path(plvi_aligned_dir, "plvi_rast_reproj_algn.tif")
plvi <- file.path(plvi_aligned_dir, "plvi_rast_reproj_algn_gte2015.tif")
plvi <- rast(plvi)

## reclassify to shrub (1 for grass)

# Define which labels map to 1 (shrub) and all others to 0
shrub_codes <- c("OS", "MS", "CS")  # open shrub, medial shrub, closed shrub

# Reclassify
plvi_shrub <- ifel(
  is.na(plvi),      # if the input cell is NA
  NA,               # keep NA
  ifel(plvi %in% shrub_codes, 1, 0))  # otherwise 1 for shrub, 0 for non-grass

# Check
plot(plvi_shrub)
plvi_shrub

# Save the raster as 8-bit unsigned integer
plvi_shrub_out_path <- file.path(plvi_aligned_dir, "plvi_reproj_algn_shrub_gte2015.tif")
writeRaster(plvi_shrub, plvi_shrub_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

################################################

# isolate each class

# load raster
#plvi <- file.path(plvi_aligned_dir, "plvi_rast_reproj_algn.tif")
plvi <- file.path(plvi_aligned_dir, "plvi_rast_reproj_algn_gte2015.tif")
plvi <- rast(plvi)

## reclassify to individual class 

# Define which labels map to 1 (desired class) and all others to 0
hg_codes <- c("HG")  # herbaceous grass

# Reclassify
plvi_hg <- ifel(
  is.na(plvi),      # if the input cell is NA
  NA,               # keep NA
  ifel(plvi %in% hg_codes, 1, 0))  # otherwise 1 for desired class, 0 for everythning else

# Check
plot(plvi_hg)
plvi_hg

# Save the raster as 8-bit unsigned integer
plvi_hg_out_path <- file.path(plvi_aligned_dir, "plvi_reproj_algn_herbgrass_gte2015.tif")
writeRaster(plvi_hg, plvi_hg_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

## reclassify to individual class 

# Define which labels map to 1 (desired class) and all others to 0
cp_codes <- c("CP")  # tame pasture

# Reclassify
plvi_cp <- ifel(
  is.na(plvi),      # if the input cell is NA
  NA,               # keep NA
  ifel(plvi %in% cp_codes, 1, 0))  # otherwise 1 for desired class, 0 for everythning else

# Check
plot(plvi_cp)
plvi_cp

# Save the raster as 8-bit unsigned integer
plvi_cp_out_path <- file.path(plvi_aligned_dir, "plvi_reproj_algn_tamepasture_gte2015.tif")
writeRaster(plvi_cp, plvi_cp_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

## reclassify to individual class 

# Define which labels map to 1 (desired class) and all others to 0
cpr_codes <- c("CPR")  # rough pasture

# Reclassify
plvi_cpr <- ifel(
  is.na(plvi),      # if the input cell is NA
  NA,               # keep NA
  ifel(plvi %in% cpr_codes, 1, 0))  # otherwise 1 for desired class, 0 for everything else

# Check
plot(plvi_cpr)
plvi_cpr

# Save the raster as 8-bit unsigned integer
plvi_cpr_out_path <- file.path(plvi_aligned_dir, "plvi_reproj_algn_roughpasture_gte2015.tif")
writeRaster(plvi_cpr, plvi_cpr_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

## reclassify to individual class 

# Define which labels map to 1 (desired class) and all others to 0
os_codes <- c("OS")  # open shrub

# Reclassify
plvi_os <- ifel(
  is.na(plvi),      # if the input cell is NA
  NA,               # keep NA
  ifel(plvi %in% os_codes, 1, 0))  # otherwise 1 for desired class, 0 for everything else

# Check
plot(plvi_os)
plvi_os

# Save the raster as 8-bit unsigned integer
plvi_os_out_path <- file.path(plvi_aligned_dir, "plvi_reproj_algn_openshrub_gte2015.tif")
writeRaster(plvi_os, plvi_os_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

## reclassify to individual class 

# Define which labels map to 1 (desired class) and all others to 0
ms_codes <- c("MS")  # medial shrub

# Reclassify
plvi_ms <- ifel(
  is.na(plvi),      # if the input cell is NA
  NA,               # keep NA
  ifel(plvi %in% ms_codes, 1, 0))  # otherwise 1 for desired class, 0 for everythning else

# Check
plot(plvi_ms)
plvi_ms

# Save the raster as 8-bit unsigned integer
plvi_ms_out_path <- file.path(plvi_aligned_dir, "plvi_reproj_algn_medialshrub_gte2015.tif")
writeRaster(plvi_ms, plvi_ms_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

## reclassify to individual class 

# Define which labels map to 1 (desired class) and all others to 0
cs_codes <- c("CS")  # closed shrub

# Reclassify
plvi_cs <- ifel(
  is.na(plvi),      # if the input cell is NA
  NA,               # keep NA
  ifel(plvi %in% cs_codes, 1, 0))  # otherwise 1 for desired class, 0 for everythning else

# Check
plot(plvi_cs)
plvi_cs

# Save the raster as 8-bit unsigned integer
plvi_cs_out_path <- file.path(plvi_aligned_dir, "plvi_reproj_algn_closedshrub_gte2015.tif")
writeRaster(plvi_cs, plvi_cs_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)
