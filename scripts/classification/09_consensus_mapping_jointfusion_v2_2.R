
###################################################################################
##                     GRASSLAND CONSENSUS MAPPING ANALYSIS                      ##
###################################################################################

##----------------------- Alison Long - May 2026 -------------------------

# step 8: use data fusion techniques to create a consensus grassland map 
#         to serve as baseline for modelling

# THIS VERSION IS OPTIMIZED FOR SAGEMAKER
# Next optimization to try: build stacks directly from S3 using GDAL VRT (zero downloads)

################################################################################
# STRATEGY A: JOINT 3‑CLASS FUSION WITH EMPIRICAL AGREEMENT PENALTIES
#
# Purpose:
#   - Produce joint probabilities for {not grass, tame grass, native grass}
#   - Use binary rasters (grass / not-grass) AND tri‑class rasters together
#   - Handle raster dependence explicitly using precomputed agreement metrics
#   - Non‑Bayesian: no priors, no likelihoods, no log‑odds
#
# Key Ideas:
#   1) Binary rasters inform grass presence only
#   2) Tri‑class rasters inform native vs tame directly
#   3) All votes are weighted by:
#        accuracy_weight × redundancy_penalty
#   4) Redundant rasters are downweighted but never removed
################################################################################

library(terra)

files <- list.files("tiles_out", full.names=TRUE, pattern="\\.tif$")
files <- Filter(Negate(is.null), results) # handle null outputs

rast_list <- lapply(files, rast)

mosaic_out <- do.call(mosaic, c(rast_list, list(
  fun="mean"
)))

writeRaster(
  mosaic_out,
  "joint_prob_full.tif",
  overwrite=TRUE,
  wopt=list(
    datatype="FLT4S",
    gdal=c("COMPRESS=ZSTD","TILED=YES","BIGTIFF=YES")
  )
)
