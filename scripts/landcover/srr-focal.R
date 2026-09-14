source(file.path("scripts/landcover", "0-setup.R"))
terraOptions(memfrac = 0.5, todisk = FALSE, tempdir = "D:/terra-temp")

dem     <- rast(file.path("./processeddatasets", dem_filepath,
                          "MRDEM Digital Terrain Model (COG)_modelling-extent.tif"))
out_dir <- file.path("./processeddatasets", dem_filepath)
w       <- 27
gopt    <- c("COMPRESS=DEFLATE", "PREDICTOR=3", "TILED=YES")
out <- function(nm) file.path(out_dir, paste0(nm, "_modelling-extent.tif"))

for (f in c("mean", "max", "min")) {
  cat("starting", f, format(Sys.time()), "\n")
  focal(dem, matrix(1, w, w), f, na.rm = TRUE,
        filename = out(paste0("terr_", f, "-27x27")),
        overwrite = TRUE, wopt = list(gdal = gopt))
  cat("done", f, format(Sys.time()), "\n")
}