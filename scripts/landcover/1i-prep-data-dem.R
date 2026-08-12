# session 1: fetch only
source(file.path("scripts/landcover", "0-setup.R"))
gopt    <- c("COMPRESS=DEFLATE", "PREDICTOR=3", "TILED=YES")

dem <- rast(file.path("./processeddatasets", dem_filepath,
                      "MRDEM Digital Terrain Model (COG)_modelling-extent.tif"))
e   <- ext(dem)
es  <- ext(e[1], e[2], e[3] - 2000, e[3])
src <- rast("/vsicurl/https://canelevation-dem.s3.ca-central-1.amazonaws.com/mrdem-30/mrdem-30-dtm.tif")

project(crop(src, project(es, crs(dem), crs(src))),
        rast(es, resolution = res(dem), crs = crs(dem)),
        method = "bilinear",
        filename = file.path("./processeddatasets", dem_filepath, "mrdem-strip.tif"),
        overwrite = TRUE,
        gdal = gopt)
