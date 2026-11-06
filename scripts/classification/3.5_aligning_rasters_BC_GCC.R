

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

# BC GCC directory
gcc_root <- file.path(raw_grass_data_dir, "Grassland/BC Grassland Conservation Council")
gcc_aligned_dir <- file.path(proc_grass_data_dir, "bc_gcc/aligned")
gcc_cropped_dir <- file.path(proc_grass_data_dir, "bc_gcc/cropped")
gcc_rast_dir <- file.path(proc_grass_data_dir, "bc_gcc/rasterized")
gcc_pilot_dir <- file.path(proc_grass_data_dir, "bc_gcc/pilot")

####--- Files ---####
# Pilot area
pa_file <- file.path(pa_dir, "pilot_area_aligned_cnvc_Jan2026.tif")

# Study area
sa_file <- file.path(sa_dir, "class_study_area_aligned_cnvc_Dec2025.tif")

# CNVC raster
#cnvc_file <- file.path(raw_grass_data_dir, "ns_cnvc/NatureServe_CNVC_COMPRESSED_20251023.tif")

# gcc raster
gcc_file <- file.path(gcc_root, "BC_grass_only_updated.shp")

##### use study area file as snap raster #####

# Load study area raster
sa <- rast(sa_file)

# print raster details 
origin(sa) # 6.696121 -2.550334
sa

####--- ALGIN RASTERS ---#### 

#----- Central Grasslands Roadmap ----#

# Load rasters
v_gcc <- vect(gcc_file)
v_gcc

#class       : SpatVector 
#geometry    : polygons 
#dimensions  : 14779, 13  (geometries, attributes)
#extent      : 1023000, 1847808, 462970.5, 1341179  (xmin, xmax, ymin, ymax)
#source      : BC_grass_only_updated.shp
#coord. ref. : PCS_Albers 
#names       : OBJECTID_1 NPFORESTDE FOREST_SPE      GRASSLANDS GRS_ASSO_V GRS_ASSO_S GRS_ASSO_C
#type        :      <num>      <chr>      <chr>           <chr>      <chr>      <chr>      <chr>
#values      :          1         NA         NA burned forest ~         NA         NA         NA

##---- Rasterize vector data ----##

# ---- 1) Ensure both layers use the same CRS as 'sa' ----
if (!same.crs(v_gcc, sa)) {
  v_gcc <- project(v_gcc, sa)  # reproject polygons to match 'sa'
}

# ---- 2) Build a categorical code for GRASSLANDS ----
# terra rasters can't store character values; we map each category to an integer
cats  <- sort(unique(v_gcc$GRASSLANDS))
codes <- seq_along(cats)
# attach integer code to the vector attributes
v_gcc$grass_type <- match(v_gcc$GRASSLANDS, cats)

# ---- 3) Create an empty raster aligned to 'sa' ----
r_tmpl <- rast(sa)  # inherits extent, resolution, origin, and CRS from 'sa'

# ---- 4) Rasterize polygons using the integer code ----
# For polygons, 'touches=TRUE' helps fill edge pixels; you can set it to FALSE if undesired.
r_code <- rasterize(v_gcc, y = r_tmpl,
                    field     = "grass_type",
                    touches   = TRUE,
                    background = NA)


# enforce NAs from sa by masking to 'sa'’s valid footprint:
r_code <- mask(r_code, sa)

# ---- 5) Make the raster categorical (factor) and attach the lookup table (RAT) ----
r_cat <- as.factor(r_code)
# Set the level table so that cell values map back to the two-letter codes
levels(r_cat) <- data.frame(ID = codes, GRASSLANDS = cats)

plot(r_cat)

# ---- 6) Write to disk ----
# Save the raster as 16-bit unsigned integer
gcc_out_path <- file.path(gcc_aligned_dir, "gcc_rast_reproj_algn.tif")
writeRaster(r_cat, gcc_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# ---- 7) Quick check ----
levels(r_cat)       # view the RAT
plot(r_cat)         # quick look (legend will show GRASSLANDS codes)
r_cat

#class       : SpatRaster 
#size        : 86539, 130223, 1  (nrow, ncol, nlyr)
#resolution  : 30, 30  (x, y)
#extent      : -2321003, 1585687, 270327.4, 2866497  (xmin, xmax, ymin, ymax)
#coord. ref. : North_America_Albers_Equal_Area_Conic (ESRI:102008) 
#source      : spat_a8c849a23260_43208_zxzxx9jyyx6UnnI.tif 
#varname     : class_study_area_aligned_cnvc_Dec2025 
#categories  : GRASSLANDS 
#name        :                           GRASSLANDS 
#min value   :   burned forest in PP or BG BGC zone 
#max value   : open grasslands with logging history 

# print frequency table
freq_table <- freq(r_cat) 
print(freq_table)

#  layer                                       value   count
#1     1          burned forest in PP or BG BGC zone   65924
#2     1              grassland associated clay bank  195792
#3     1 open dry forest adjacent to open grasslands  332452
#4     1                     open dry forest in NDT4  238015
#5     1                             open grasslands 7797733
#6     1        open grasslands with logging history   35575

## reclassify to binary (0 = not grass, 1 = grass)

# Confirm the levels table and the column name that holds the labels
levels(r_cat)[[1]] |> head()

# Define which labels map to 1 (grass) and all others to 0
grass_codes <- c("open grasslands")  

# Binary raster: 1 for grass classes, 0 otherwise (except NA)
gcc_algn_binary <- ifel(
  is.na(r_cat),                     # if the input cell is NA
  NA,                               # keep NA
  ifel(r_cat %in% grass_codes, 1, 0)  # otherwise 1 for grass, 0 for non-grass
)

# 5) Quick checks
gcc_algn_binary
freq(gcc_algn_binary)
plot(gcc_algn_binary) 

# Save the raster as 8-bit unsigned integer
gcc_binary_out_path <- file.path(gcc_aligned_dir, "gcc_reproj_algn_binary.tif")
writeRaster(gcc_algn_binary, gcc_binary_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

## not present in pilot area!!
