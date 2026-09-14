# northness / eastness at quarter-section scale (27x27 = 810 m)
# tomorrows-prairies
#
# Plane fitted by least squares over the window; gradients obtained as weighted
# focal sums (one pass each) rather than via a quadratic fit, which is far cheaper.
# Sign convention verified against synthetic ramps: positive northness = the
# surface faces north. Aspect is atan2(-gx, -gy), both components negated so the
# angle describes the direction the slope faces rather than the direction of ascent.

source(file.path("scripts/landcover", "0-setup.R"))
terraOptions(memfrac = 0.5, tempdir = "D:/terra-temp")

sa_r <- rast(file.path(sa_filepath, sa_raster_filename))
dem     <- rast(file.path("./processeddatasets", dem_filepath,
                          "MRDEM Digital Terrain Model (COG)_modelling-extent.tif"))
out_dir <- file.path("./processeddatasets", dem_filepath)
w       <- 27
gopt    <- c("COMPRESS=DEFLATE", "PREDICTOR=3", "TILED=YES")
out <- function(nm) file.path(out_dir, paste0(nm, "_modelling-extent.tif"))

# ---- plane-fit weights ----
# dz/dx = sum(x_i * z_i) / sum(x_i^2), offsets in metres from the window centre.
# Denominator is a constant, so each gradient is a single weighted focal sum.
h   <- (w - 1) / 2
d   <- res(dem)[1]
off <- seq(-h, h) * d
wx  <- matrix(rep(off,      each  = w), w, w)   # varies across columns, west negative
wy  <- matrix(rep(rev(off), times = w), w, w)   # varies down rows, north positive
s   <- sum(off^2) * w

# ---- tile, compute, trim, mosaic ----
# focal() at this window crashes on the full extent; tiles of ~69M cells are fine.
# Buffer is (w-1)/2 = 13 cells. Tiles MUST be trimmed back to their unbuffered
# extent before mosaicking -- buffer cells have truncated windows, and letting
# them into the mosaic produces visible seams.
tmpl  <- rast(nrows = 4, ncols = 4, extent = ext(dem), crs = crs(dem))
tiles <- makeTiles(dem, tmpl, filename = "D:/terra-temp/dem_tile_.tif",
                   extend = TRUE, buffer = h, overwrite = TRUE)

for (g in c("x", "y")) {
  wt   <- if (g == "x") wx else wy
  outs <- character(0)
  for (i in seq_along(tiles)) {
    cat("dz", g, "tile", i, format(Sys.time()), "\n")
    o  <- crop(focal(rast(tiles[i]), wt / s, "sum", na.rm = TRUE), ext(tmpl, i))
    fn <- sprintf("D:/terra-temp/dz%s_tile_%02d.tif", g, i)
    writeRaster(o, fn, overwrite = TRUE)
    outs <- c(outs, fn)
  }
  mosaic(sprc(lapply(outs, rast)), fun = "first",
         filename = out(paste0("terr_dzd", g, "-27x27")),
         overwrite = TRUE, gdal = gopt)
  cat("done dz", g, format(Sys.time()), "\n")
}

dzdx <- rast(out("terr_dzdx-27x27"))
dzdy <- rast(out("terr_dzdy-27x27"))

# ---- Stage (1976) transformations ----
# slope * cos(aspect) and slope * sin(aspect). Self-damping: as slope goes to
# zero both go to zero regardless of the (ill-determined) aspect there, which is
# the correct answer for flat ground and avoids masking out most of the extent.
lapp(c(dzdx, dzdy), function(gx, gy) {
  sl <- atan(sqrt(gx^2 + gy^2))
  as <- atan2(-gx, -gy)
  sl * cos(as)
}, filename = out("terr_northness-27x27"), overwrite = TRUE, wopt = list(gdal = gopt))

lapp(c(dzdx, dzdy), function(gx, gy) {
  sl <- atan(sqrt(gx^2 + gy^2))
  as <- atan2(-gx, -gy)
  sl * sin(as)
}, filename = out("terr_eastness-27x27"), overwrite = TRUE, wopt = list(gdal = gopt))

# ---- optional: slope from the same fit, radians ----
# lapp(c(dzdx, dzdy), function(gx, gy) atan(sqrt(gx^2 + gy^2)),
#      filename = out("terr_slope-27x27"), overwrite = TRUE, wopt = list(gdal = gopt))



n <- rast(out("terr_northness-27x27"))
e <- rast(out("terr_eastness-27x27"))

n_sa <- mask(crop(n, sa_r), sa_r)
e_sa <- mask(crop(e, sa_r), sa_r)

# 1. range and NA count
minmax(n_sa); minmax(e_sa)
global(is.na(n_sa), "sum")[[1]]        # vs 437,104,948

# 2. distribution -- the key question is whether there's usable variance
s <- spatSample(c(n_sa, e_sa), 1e5, method = "regular", na.rm = TRUE)
names(s) <- c("north", "east")
round(sapply(s, function(x) c(mean = mean(x), sd = sd(x),
                              q01 = quantile(x, .01), q25 = quantile(x, .25),
                              q50 = median(x), q75 = quantile(x, .75),
                              q99 = quantile(x, .99))), 5)

# 3. independence
cor(s$north, s$east, method = "spearman")

# 4. seam check
xb <- ext(dem)[1] + (ext(dem)[2] - ext(dem)[1]) / 4
plot(terra::extract(n, cbind(seq(xb - 600, xb + 600, by = 30), 1403000))[,1], type = "l")



