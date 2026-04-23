

###################################################################################
##                    FINALIZE BASELINE MAP FOR MODELLING                        ##
###################################################################################

##----------------------- Alison Long - March 2026 -------------------------

# step 9: Fill in non-grassland pixels of the consensus map with other land cover
#         classes to generate baseline map for landscape modelling

# load packages
library(terra)
library(sf)
library(tidyverse)

# Set the tempocmry directory to a drive with sufficient space
terraOptions(tempdir = "R:/CanadianGrasslands/WG_classification/scratch")

####--- DIRECTORIES ---####

# pilot area
pa_dir <- "R:/CanadianGrasslands/processeddatasets/pilotarea"

# study area
sa_dir <- "R:/CanadianGrasslands/processeddatasets/studyarea"

# grasslands
raw_grass_data_dir  <- "R:/CanadianGrasslands/existingdatasets"
proc_grass_data_dir <- "R:/CanadianGrasslands/WG_Classification/processed"
consens_data_dir    <- "R:/CanadianGrasslands/WG_Classification/results/consensus_mapping"
out_dir             <- "R:/CanadianGrasslands/WG_Classification/results/consensus_mapping"

# grass raster directories
luts_aligned_dir  <- file.path(proc_grass_data_dir, "aafc_luts/aligned")

####--- FILES ---####

# Pilot area
pa_file <- file.path(pa_dir, "pilot_area_aligned_cnvc_Jan2026.tif")

# Study area
sa_file <- file.path(sa_dir, "modelling_extent_aligned_cnvc_Jan2026.tif")

# ecodistricts
eco_file <- file.path(eco_dir, "ecodistricts_clean.shp")

# grass raster files
luts_file <- file.path(luts_aligned_dir, "aafc_luts_merged_reproj_algn_natvtam_noshrub.tif")

# consensus map file
cons_file <- file.path(consens_data_dir, "??")

####--- PREP LC cmSTER ---####

# Read cmsters
sa   <- rast(sa_file)
luts  <- rast(luts_file)
cm   <- rast(cons_file)

####--- FILL NA PIXELS ---####

# Validate the rasters align 
stopifnot(
  all.equal(ext(cm), ext(luts)),
  all.equal(res(cm), res(luts)),
  crs(cm) == crs(luts)
)

# Fill NA pixels in A with values from B
# cover(A, B) returns A where not NA, otherwise B
filled_cm <- cover(cm, luts)
plot(filled_cm)

####--- MASK TO MODELLING EXTENT ---####

# Mask the result to the extent of sa
# `mask()` keeps only cells where sa is not NA
# `crop()` ensures the output extent exactly matches sa
filled_masked_cm <- mask(crop(filled_cm, sa), sa)
plot(filled_masked_cm)

####--- EXPORT ---####

# Write result efficiently for large rasters (chunked writing)
writeRaster(filled_masked_cm,
            filename = "TP_baseline_map_Mar2026.tif",
            overwrite = TRUE,
            wopt = list(gdal = c("COMPRESS=LZW", "TILED=YES")))
