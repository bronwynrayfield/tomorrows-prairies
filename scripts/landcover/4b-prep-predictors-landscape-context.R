# =============================================================================
# Purpose:
#   Prepares intermediate binary and mask layers needed to derive spatial
#   predictors in subsequent scripts.
#   All outputs are at native 30m resolution and saved to processeddatasets/.
#   Resolution-dependent operations (focal statistics etc.) are in 5-*.R
#
# Inputs:
#   - AAFC LUTS 30m rasters (2010, 2015) : projected and cropped
#
# Outputs:
#   - AAFC binary rasters (native, tame, crop, settlement, forest, water)
#     for 2010 and 2015
#
# Notes:
#   - crop class (codes 52, 56) not present in 2010 AAFC LUTS;
#
# Dependencies:
#   - 0-setup.R (parameters, file paths, functions)
#   - 3-crop-data.R outputs
# =============================================================================

source(file.path("scripts/landcover", "0-setup.R"))

#--- Load files ----
# Study area
sa <- st_read(file.path(sa_filepath, sa_filename))
me <- st_read(file.path(sa_filepath, me_filename))
me_r <- rast(file.path(sa_filepath, me_raster_filename))

# Land cover (AAFC LUTS)
lulc2010 <- rast(tagged_path(file.path("./processeddatasets", aafc_filepath), "cropped", aafc2010_filename, "modelling-extent"))
lulc2015 <- rast(tagged_path(file.path("./processeddatasets", aafc_filepath), "cropped", aafc2015_filename, "modelling-extent"))
lulc2020 <- rast(tagged_path(file.path("./processeddatasets", aafc_filepath), "cropped", aafc2020_filename, "modelling-extent"))


#--- LULC binary layers ----
# --- 2015 binary layers ---
native2015_bin     <- reclass_aafc(lulc2015, scheme = "binary", binary_class = "native",     check = TRUE)
tame2015_bin       <- reclass_aafc(lulc2015, scheme = "binary", binary_class = "tame",       check = TRUE)
crop2015_bin       <- reclass_aafc(lulc2015, scheme = "binary", binary_class = "crop",       check = TRUE)
settlement2015_bin <- reclass_aafc(lulc2015, scheme = "binary", binary_class = "settlement", check = TRUE)
forest2015_bin     <- reclass_aafc(lulc2015, scheme = "binary", binary_class = "forest",     check = TRUE)
water2015_bin      <- reclass_aafc(lulc2015, scheme = "binary", binary_class = "water",      check = TRUE)

# --- 2010 binary layers ---
native2010_bin     <- reclass_aafc(lulc2010, scheme = "binary", binary_class = "native",     check = TRUE)
# Note that AAFC LUTS 2010 doesn't separate tame and annual cropland, both are lumped into "tame"
tame2010_bin       <- reclass_aafc(lulc2010, scheme = "binary", binary_class = "tame",       check = TRUE)
settlement2010_bin <- reclass_aafc(lulc2010, scheme = "binary", binary_class = "settlement", check = TRUE)
forest2010_bin     <- reclass_aafc(lulc2010, scheme = "binary", binary_class = "forest",     check = TRUE)
water2010_bin      <- reclass_aafc(lulc2010, scheme = "binary", binary_class = "water",      check = TRUE)


