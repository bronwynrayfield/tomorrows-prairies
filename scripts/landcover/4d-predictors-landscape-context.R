# =============================================================================
# Purpose:
#   Prepares intermediate binary and mask layers needed to derive spatial
#   predictors in subsequent scripts.
#   All outputs are at native 30m resolution and saved to processeddatasets/.
#   Resolution-dependent operations (focal statistics etc.) are in 5-*.R
#
# Inputs:
#   - AAFC LUTS 30m rasters (2010, 2015, 2020) : projected and cropped
#
# Outputs:
#   - AAFC binary rasters (native, tame, crop, settlement, forest, water)
#     for 2010, 2015, and 2020
#
# Notes:
#   - crop class (codes 52, 56) not present in 2010 AAFC LUTS, so "crop" is
#     omitted from that year's class list (see reclass_aafc() docs)
#
# Dependencies:
#   - 0-setup.R (parameters, file paths, functions)
#   - 3-crop-data.R outputs
# =============================================================================
source(file.path("scripts/landcover", "0-setup.R"))

# Output directory
out_dir <- file.path("./processeddatasets", aafc_filepath, "binary")

#--- Load files ----
sa   <- st_read(file.path(sa_filepath, sa_filename))
me   <- st_read(file.path(sa_filepath, me_filename))
me_r <- rast(file.path(sa_filepath, me_raster_filename))

lulc <- list(
  `2010` = rast(tagged_path(file.path("./processeddatasets", aafc_filepath), "cropped", aafc2010_filename, "modelling-extent")),
  `2015` = rast(tagged_path(file.path("./processeddatasets", aafc_filepath), "cropped", aafc2015_filename, "modelling-extent")),
  `2020` = rast(tagged_path(file.path("./processeddatasets", aafc_filepath), "cropped", aafc2020_filename, "modelling-extent"))
)

classes <- list(
  `2010` = c("native", "tame", "settlement", "forest", "water"),
  `2015` = c("native", "tame", "crop", "settlement", "forest", "water"),
  `2020` = c("native", "tame", "crop", "settlement", "forest", "water")
)

#--- LULC binary layers ----
bin_layers <- list()

# Reclass to binary layer and save result to disk
for (yr in names(lulc)) {
  bin_layers[[yr]] <- list()
  for (cls in classes[[yr]]) {
    r <- reclass_aafc(lulc[[yr]], scheme = "binary", binary_class = cls, check = TRUE)
    writeRaster(
      r, file.path(out_dir, sprintf("AAFC%s_%s.tif", yr, cls)),
      overwrite = TRUE, gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2")
    )
    bin_layers[[yr]][[cls]] <- r
  }
}


