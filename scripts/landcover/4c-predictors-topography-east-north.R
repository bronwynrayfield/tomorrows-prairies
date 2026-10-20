# =============================================================================
# Purpose: 
# Calculate northness / eastness at quarter-section scale (27x27 = 810 m)
# AWS instance - r6i.2xlarge (8 vCPU, 64 GB)
#
# Plane fitted by least squares over the window, from centred focal sums, so
# truncated windows (southern border, DEM edges, interior NA) are fitted
# correctly.
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

# for parrallel processing
n_workers <- 7

# ---- window offsets (metres from centre) ----
h    <- (w - 1) / 2
d    <- res(dem)[1]
off  <- seq(-h, h) * d
ones <- matrix(1, w, w)
wx   <- matrix(rep(off,      each  = w), w, w)   # varies across columns, west negative
wy   <- matrix(rep(rev(off), times = w), w, w)   # varies down rows, north positive

# centring keeps the Sy*Sz/n cancellation well inside float32 temp-file precision;
# a plane fit is invariant to the offset
mu <- global(dem, "mean", na.rm = TRUE)[[1]]

# ---- tile, compute, trim, mosaic ----
# focal() at this window crashes on the full extent; tiles of ~69M cells are fine.
# Buffer is (w-1)/2 = 13 cells. Tiles are trimmed to their unbuffered extent
# before mosaicking -- buffer cells have truncated windows and produce seams.
tmpl  <- rast(nrows = 8, ncols = 8, extent = ext(dem), crs = crs(dem))
tiles <- makeTiles(dem, tmpl, filename = "D:/terra-temp/dem_tile_.tif",
                   extend = TRUE, buffer = h, overwrite = TRUE)

# 4 workers: 8 vCPU on r6i is 4 physical cores + hyperthreads, so expect ~3-4x, not 8x.
# Each tile holds ~69M cells; ~15 live double layers is roughly 8 GB per worker,
# so 4 workers x 8 GB fits in 64 GB without spilling to EBS (gp3 default 125 MB/s).
bx <- lapply(seq_along(tiles), function(i) as.vector(ext(tmpl, i)))   # SpatExtent can't cross processes

do_tile <- function(i) {
  fs <- function(r, k) focal(r, k, "sum", na.rm = TRUE)
  cat("tile", i, format(Sys.time()), "\n")
  z <- rast(tiles[i]) - mu
  v <- !is.na(z)
  
  n   <- fs(v, ones);  Sz  <- fs(z, ones)
  Sx  <- fs(v, wx);    Sxx <- fs(v, wx^2);  Sxz <- fs(z, wx)
  Sy  <- fs(v, wy);    Syy <- fs(v, wy^2);  Syz <- fs(z, wy)
  Sxy <- fs(v, wx * wy)
  
  cxx <- Sxx - Sx^2 / n;     cyy <- Syy - Sy^2 / n;     cxy <- Sxy - Sx * Sy / n
  cxz <- Sxz - Sx * Sz / n;  cyz <- Syz - Sy * Sz / n
  det <- cxx * cyy - cxy^2
  
  b <- ext(bx[[i]])
  writeRaster(crop(n, b),                             tmp("ncells", i), overwrite = TRUE)
  writeRaster(crop((cyy * cxz - cxy * cyz) / det, b), tmp("dzdx",   i), overwrite = TRUE)
  writeRaster(crop((cxx * cyz - cxy * cxz) / det, b), tmp("dzdy",   i), overwrite = TRUE)
  rm(list = setdiff(ls(), "i")); gc() 
  i
}

cl <- makeCluster(n_workers, outfile = "D:/terra-temp/workers.log")   # cat() output lands here
clusterExport(cl, c("tiles", "bx", "mu", "ones", "wx", "wy", "tmp"))
clusterEvalQ(cl, { library(terra); terraOptions(memfrac = 0.12, tempdir = "D:/terra-temp"); NULL })
invisible(parLapplyLB(cl, seq_along(tiles), do_tile, chunk.size = 1))
stopCluster(cl)

# ncells is the empirical window count TRI needs - reuse it there instead of n_r
for (nm in c("ncells", "dzdx", "dzdy")) {
  mosaic(sprc(lapply(tmp(nm, seq_along(tiles)), rast)), fun = "first",
         filename = out(paste0(nm, "-27x27")),
         overwrite = TRUE, gdal = gopt)
  cat("done", nm, format(Sys.time()), "\n")
}

# masked to the DEM so cells with NA elevation get no value
dzdx <- mask(rast(out("dzdx-27x27")), dem)
dzdy <- mask(rast(out("dzdy-27x27")), dem)


# ---- Stage (1976) transformations ----
# slope * cos(aspect) and slope * sin(aspect). Self-damping: as slope goes to
# zero both go to zero regardless of the (ill-determined) aspect there, which is
# the correct answer for flat ground and avoids masking out most of the extent.
lapp(c(dzdx, dzdy), function(gx, gy) {
  atan(sqrt(gx^2 + gy^2)) * cos(atan2(-gx, -gy))
}, filename = out("northness-27x27"), overwrite = TRUE, wopt = list(gdal = gopt))

lapp(c(dzdx, dzdy), function(gx, gy) {
  atan(sqrt(gx^2 + gy^2)) * sin(atan2(-gx, -gy))
}, filename = out("eastness-27x27"), overwrite = TRUE, wopt = list(gdal = gopt))

# not used: slope from the same fit, radians
# lapp(c(dzdx, dzdy), function(gx, gy) atan(sqrt(gx^2 + gy^2)),
#      filename = out("slope-27x27"), overwrite = TRUE, wopt = list(gdal = gopt))


# ---- QA/QC ----
n <- rast(out("northness-27x27"))
e <- rast(out("eastness-27x27"))
n_sa <- mask(crop(n, sa_r), sa_r)
e_sa <- mask(crop(e, sa_r), sa_r)

# 1. range and NA count
minmax(n_sa)
# lyr1
# min -0.2629849
# max  0.3287703
minmax(e_sa)
# lyr1
# min -0.2907141
# max  0.3460160
global(is.na(n_sa), "sum")[[1]]
# [1] 437104948 

# 2. window count: 729 interior, down to 14 * 27 = 378 on the border row
minmax(rast(out("ncells-27x27")))
# focal_sum
# min         0
# max       729

# 3. distribution - the key question is whether there's usable variance
s <- spatSample(c(n_sa, e_sa), 1e5, method = "regular", na.rm = TRUE)
names(s) <- c("north", "east")
round(sapply(s, function(x) c(mean = mean(x), sd = sd(x),
                              q01 = quantile(x, .01), q25 = quantile(x, .25),
                              q50 = median(x), q75 = quantile(x, .75),
                              q99 = quantile(x, .99))), 5)
#          north     east
# mean     0.00049  0.00105
# sd       0.01326  0.01425
# q01.1%  -0.03905 -0.04300
# q25.25% -0.00301 -0.00295
# q50      0.00007  0.00061
# q75.75%  0.00362  0.00487
# q99.99%  0.04409  0.04560

# 4. independence
cor(s$north, s$east, method = "spearman")
# [1] 0.08471683

# 5. southern border: should look like the strip north of it, not a positive band
ed <- ext(dem)
plot(crop(n, ext(ed[1], ed[2], ed[3], ed[3] + 1000)))

# 6. tile seam check
xb <- ed[1] + (ed[2] - ed[1]) / 4
plot(terra::extract(n, cbind(seq(xb - 600, xb + 600, by = 30), 1403000))[,1], type = "l")