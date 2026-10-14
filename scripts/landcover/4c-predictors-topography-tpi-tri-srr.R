# =============================================================================
# Purpose:
# TPI, TRI and SRR at quarter-section scale (27x27 = 810 m)
# AWS instance - r6i.2xlarge (8 vCPU, 64 GB)
#
# All three come from five focal passes per tile on the mean-centred DEM:
# n (valid-cell count), S1, S2, max, min. Every quantity uses the empirical n,
# so truncated windows (southern border, SE corner, DEM edges) are handled
# exactly.
#
# Centring: the TRI expansion S2 - 2*z*S1 + n*z^2 cancels below float32
# precision on flat ground (it produced interior NAs earlier). Centring plus
# double-precision temp files on the workers works
#
# Outputs are masked to the DEM so cells with NA elevation get no value.
# Within ~390 m (13 rows) of the southern border windows are truncated to their
# northern half (down to 14 x 27 cells)
# =============================================================================

source(file.path("scripts/landcover", "0a-setup.R"))
terraOptions(memfrac = 0.5)

sa_r    <- rast(file.path(sa_filepath, sa_raster_filename))
dem     <- rast(file.path("./processeddatasets", dem_filepath,
                          "MRDEM Digital Terrain Model (COG)_modelling-extent.tif"))
out_dir <- file.path("./processeddatasets", dem_filepath)
w       <- 27
gopt    <- c("COMPRESS=DEFLATE", "PREDICTOR=3", "TILED=YES")
out <- function(nm) file.path(out_dir, paste0(nm, "_modelling-extent.tif"))
tmp <- function(nm, i) sprintf("D:/terra-temp/%s_tile_%02d.tif", nm, i)

h    <- (w - 1) / 2
ones <- matrix(1, w, w)
mu   <- global(dem, "mean", na.rm = TRUE)[[1]]

# for parallelization
n_workers <- 8

# ---- tiles ----
# 8x8 = 64 tiles of ~17M cells; buffer h = 13 cells covers the window.
# Tiles are trimmed to their unbuffered extent before mosaicking (seams otherwise).
tmpl  <- rast(nrows = 8, ncols = 8, extent = ext(dem), crs = crs(dem))
tiles <- makeTiles(dem, tmpl, filename = "D:/terra-temp/dem_tile_.tif",
                   extend = TRUE, buffer = h, overwrite = TRUE)
bx    <- lapply(seq_along(tiles), function(i) as.vector(ext(tmpl, i)))   # SpatExtent can't cross processes

do_tile <- function(i) {
  cat("tile", i, format(Sys.time()), "\n")
  z  <- rast(tiles[i]) - mu
  n  <- focal(!is.na(z), ones, "sum", na.rm = TRUE)
  s1 <- focal(z,         ones, "sum", na.rm = TRUE)
  s2 <- focal(z^2,       ones, "sum", na.rm = TRUE)
  hi <- focal(z,         ones, "max", na.rm = TRUE)
  lo <- focal(z,         ones, "min", na.rm = TRUE)
  mn <- s1 / n
  
  # TPI (Weiss 2001): elevation minus window mean, m
  tpi <- z - mn
  # TRI (adapted from Riley et al. 1999): RMS deviation of the n - 1 neighbours
  # from the focal cell, m. clamp() rather than pmax(): pmax doesn't dispatch on SpatRaster.
  tri <- sqrt(clamp(s2 - 2 * z * s1 + n * z^2, lower = 0, values = TRUE) / (n - 1))
  # SRR (Pike & Wilson 1971): (mean - min) / (max - min), 0.5 on flat windows
  srr <- ifel(hi == lo, 0.5, (mn - lo) / (hi - lo))
  
  b  <- ext(bx[[i]])
  zb <- crop(z, b)
  writeRaster(mask(crop(tpi, b), zb), tmp("tpi", i), overwrite = TRUE, datatype = "FLT4S")
  writeRaster(mask(crop(tri, b), zb), tmp("tri", i), overwrite = TRUE, datatype = "FLT4S")
  writeRaster(mask(crop(srr, b), zb), tmp("srr", i), overwrite = TRUE, datatype = "FLT4S")
  
  rm(list = setdiff(ls(), "i")); gc()
  i
}

# ---- run ----
# 5 focal passes per tile 
cl <- makeCluster(n_workers, outfile = "D:/terra-temp/workers.log")
clusterEvalQ(cl, {
  library(terra)
  # 0.1 x 64 GB per worker; FLT8S keeps intermediates double if they spill to disk
  terraOptions(memfrac = 0.1, tempdir = "D:/terra-temp", datatype = "FLT8S")
  NULL
})
clusterExport(cl, c("tiles", "bx", "mu", "ones", "tmp"))
invisible(parLapplyLB(cl, seq_along(tiles), do_tile, chunk.size = 1))
stopCluster(cl)

for (nm in c("tpi", "tri", "srr")) {
  mosaic(sprc(lapply(tmp(nm, seq_along(tiles)), rast)), fun = "first",
         filename = out(paste0("terr_", nm, "-27x27")),
         overwrite = TRUE, gdal = gopt, datatype = "FLT4S")
  cat("done", nm, format(Sys.time()), "\n")
}


# ---- QA/QC ----
r    <- c(rast(out("terr_tpi-27x27")), rast(out("terr_tri-27x27")), rast(out("terr_srr-27x27")))
names(r) <- c("tpi", "tri", "srr")
r_sa <- mask(crop(r, sa_r), sa_r)

# 1. ranges: TPI centred near 0; TRI >= 0, max ~112 m (pre-artifact run); SRR within 0-1
minmax(r_sa)
#       tpi      tri         srr
# min -97.73172   0.0000 0.002731479
# max  94.39059 112.4405 0.992999971

# 2. NA count, each layer vs 437,104,948
sapply(1:3, function(k) global(is.na(r_sa[[k]]), "sum")[[1]])
# [1] 437104948 437104948 437104948 

# 3. where is the TRI max? should be an escarpment, not the SE corner
xyFromCell(r_sa, where.max(r_sa[["tri"]])[, "cell"])
crds(project(vect(xyFromCell(r_sa, where.max(r_sa[["tri"]])[, "cell"]), crs = crs(r_sa)), "EPSG:4326"))

# 4. distribution and redundancy
s <- spatSample(r_sa, 1e5, method = "regular", na.rm = TRUE)
round(sapply(s, function(x) c(mean = mean(x), sd = sd(x),
                              q01 = quantile(x, .01), q50 = median(x),
                              q99 = quantile(x, .99))), 3)
round(cor(s, method = "spearman"), 3)

# 5. southern border: edge cells should look like the terrain just north of them
p  <- crds(project(vect(cbind(-105, 49), crs = "EPSG:4326"), crs(dem)))
bz <- ext(p[1] - 10000, p[1] + 10000, p[2] - 3000, p[2] + 5000)
par(mfrow = c(1, 3)); for (k in 1:3) plot(crop(r_sa[[k]], bz), main = names(r)[k]); par(mfrow = c(1, 1))

# 6. tile seams: a vertical and a horizontal profile across tile boundaries
ed <- ext(dem)
xb <- ed[1] + (ed[2] - ed[1]) / 4
yb <- ed[3] + (ed[4] - ed[3]) / 4
plot(terra::extract(r[["tri"]], cbind(seq(xb - 600, xb + 600, by = 30), yb + 5000))[, 1], type = "l")
plot(terra::extract(r[["tri"]], cbind(xb + 5000, seq(yb - 600, yb + 600, by = 30)))[, 1], type = "l")