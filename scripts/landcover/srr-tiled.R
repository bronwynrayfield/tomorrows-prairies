source(file.path("scripts/landcover", "0-setup.R"))
terraOptions(memfrac = 0.5, tempdir = "D:/terra-temp")

dem     <- rast(file.path("./processeddatasets", dem_filepath,
                          "MRDEM Digital Terrain Model (COG)_modelling-extent.tif"))
out_dir <- file.path("./processeddatasets", dem_filepath)
w       <- 27
gopt    <- c("COMPRESS=DEFLATE", "PREDICTOR=3", "TILED=YES")
out <- function(nm) file.path(out_dir, paste0(nm, "_modelling-extent.tif"))

tmpl <- rast(nrows = 4, ncols = 4, extent = ext(dem), crs = crs(dem))

tiles <- makeTiles(dem, tmpl,
                   filename = "D:/terra-temp/dem_tile_.tif",
                   extend = TRUE, buffer = 13, overwrite = TRUE)

for (f in c("mean", "max", "min")) {
  outs <- character(0)
  for (i in seq_along(tiles)) {
    cat(f, "tile", i, format(Sys.time()), "\n")
    o  <- crop(focal(rast(tiles[i]), matrix(1, w, w), f, na.rm = TRUE), ext(tmpl, i))
    fn <- sprintf("D:/terra-temp/%s_tile_%02d.tif", f, i)
    writeRaster(o, fn, overwrite = TRUE)
    outs <- c(outs, fn)
  }
  mosaic(sprc(lapply(outs, rast)), fun = "first",
         filename = out(paste0("terr_", f, "-27x27")),
         overwrite = TRUE, gdal = gopt)
}

mn <- rast(out("terr_mean-27x27"))
hi <- rast(out("terr_max-27x27"))
lo <- rast(out("terr_min-27x27"))

lapp(c(mn, hi, lo), function(mn, hi, lo) {
  r <- (mn - lo) / (hi - lo)
  r[hi == lo] <- 0.5
  r
}, filename = out("terr_srr-27x27"), overwrite = TRUE, wopt = list(gdal = gopt))


srr <- rast(out("terr_srr-27x27"))
srr_sa <- mask(crop(srr, sa_r), sa_r)
minmax(srr_sa)                                    # expect within 0-1
global(is.na(srr_sa), "sum")[[1]]                 # vs 437,104,948


tri <- rast(out("terr_tri-27x27"))
srr <- rast(out("terr_srr-27x27"))
tpi <- rast(out("terr_tpi-27x27"))

s <- spatSample(mask(crop(c(tri, srr, tpi), sa_r), sa_r), 1e5,
                method = "regular", na.rm = TRUE)
names(s) <- c("tri", "srr", "tpi")

round(cor(s, method = "spearman"), 3)
round(sapply(s, function(x) c(mean = mean(x), sd = sd(x),
                              q01 = quantile(x, .01), q50 = median(x),
                              q99 = quantile(x, .99))), 3)

