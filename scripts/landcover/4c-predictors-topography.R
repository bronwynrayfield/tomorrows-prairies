# terrain covariates from MRDEM -- tomorrows-prairies

source(file.path("scripts/landcover", "0-setup.R"))

sa_r <- rast(file.path(sa_filepath, sa_raster_filename))
dem   <- rast(file.path("./processeddatasets", dem_filepath, "MRDEM Digital Terrain Model (COG)_modelling-extent.tif"))
wins    <- c(3, 27)                      # cells; at 30 m = 90 m and 810 m

# ---- Slope and Aspect ----
slope  <- terrain(dem, "slope",  unit="degrees")
aspect <- terrain(dem, "aspect", unit="degrees")
aspect_mask <- ifel(slope < 0.5, NA, aspect)

write_rast_tagged(slope, "./processeddatasets", dem_filepath, slope_filename, "modelling-extent")
write_rast_tagged(aspect_mask, "./processeddatasets", dem_filepath, aspect_filename, "modelling-extent")
rm(slope, aspect, aspect_mask)

# ---- Slope-aspect transformations ----
slp <- terrain(dem, v = "slope",  unit = "radians")
asp <- terrain(dem, v = "aspect", unit = "radians")

write_rast_tagged(slp * cos(asp), "./processeddatasets", dem_filepath, "terr_slope-cos-aspect", "modelling-extent")
write_rast_tagged(slp * sin(asp), "./processeddatasets", dem_filepath, "terr_slope-sin-aspect", "modelling-extent")
rm(slp, asp)

# ---- Hierarchical slope position (Murphy et al. 2010, per scale) ----
tpi <- function(z, w) z - focal(z, matrix(1, w, w), "mean", na.rm = TRUE)

for (w in wins) {
  tag <- sprintf("%02dx%02d", w, w)
  write_rast_tagged(tpi(dem, w), "./processeddatasets", dem_filepath, paste0("terr_tpi-", tag), "modelling-extent")
}

# QA/QC - TPI layers within study area
tpi3  <- rast(file.path("./processeddatasets/", dem_filepath, "terr_tpi-03x03_modelling-extent.tif"))
tpi27  <- rast(file.path("./processeddatasets/", dem_filepath, "terr_tpi-27x27_modelling-extent.tif"))

tpi3 <- mask(crop(tpi3, sa_r), sa_r)
tpi27 <- mask(crop(tpi27, sa_r), sa_r)
dem_sa <- mask(crop(dem, sa_r), sa_r)

# regular sample keeps this cheap; NAs retained so coverage can be assessed
s <- spatSample(c(dem_sa, tpi3, tpi27), 1e5, method = "regular",
                na.rm = FALSE, xy = TRUE)
names(s) <- c("x", "y", "dem", "tpi3", "tpi27")

# - 1. coverage: TPI NAs should track DEM NAs
colSums(is.na(s[, c("dem", "tpi3", "tpi27")])) / nrow(s)

k <- complete.cases(s)
d <- s[k, ]

# - 2. centring and spread
# TPI is a deviation from a local mean, so mean ~ 0 is definitional.
# A mean far off zero means the focal window is not doing what you think.
round(sapply(d[, c("tpi3", "tpi27")],
             function(x) c(mean = mean(x), median = median(x), sd = sd(x),
                           min = min(x), max = max(x),
                           q01 = quantile(x, .01), q99 = quantile(x, .99))), 3)

# - 3. scale ordering: sd(27) must exceed sd(3)
sd(d$tpi27) > sd(d$tpi3)

# - 4. cross-scale redundancy
cor(d$tpi3, d$tpi27, method = "spearman")

# - 5. quantisation: discrete banding means DEM precision is the limit
length(unique(round(d$tpi3, 3)))
hist(d$tpi3[abs(d$tpi3) < 2], breaks = 200, main = "TPI 3x3 (core)")

# - 6. TPI must be independent of elevation itself
cor(d$dem, d$tpi3,  method = "spearman")
cor(d$dem, d$tpi27, method = "spearman")

# - 7. visual: a tile with known relief, and one at the extent edge
e   <- ext(dem)
mid <- c(mean(e[1:2]), mean(e[3:4]))
sz  <- 500 * res(dem)[1]
tl  <- ext(mid[1] - sz, mid[1] + sz, mid[2] - sz, mid[2] + sz)

par(mfrow = c(1, 3))
plot(crop(dem,   tl), main = "DEM")
plot(crop(tpi3,  tl), main = "TPI 3x3",   range = c(-2, 2))
plot(crop(tpi27, tl), main = "TPI 27x27", range = c(-8, 8))
par(mfrow = c(1, 1))

# ---- Roughness, terrain ruggedness (Riley et al., 1999) ----
# sum_i (z_i - z_c)^2  expands to  S2 - 2*z_c*S1 + n*z_c^2, so three fast
# focal sums replace a per-window R callback (matters a lot at w = 27).
# rms = TRUE returns root-mean-square deviation, which is comparable across
# window sizes; Riley's original (rms = FALSE) scales with n.
tri <- function(z, w, rms = TRUE) {
  m  <- matrix(1, w, w)
  n  <- focal(!is.na(z), m, "sum", na.rm = TRUE)
  s1 <- focal(z,         m, "sum", na.rm = TRUE)
  s2 <- focal(z^2,       m, "sum", na.rm = TRUE)
  ss <- s2 - 2 * z * s1 + n * z^2
  if (rms) sqrt(ss / (n - 1)) else sqrt(ss)
}

# ---- Pike & Wilson (1971) surface relief ratio ----
# (mean - min) / (max - min); undefined on perfectly flat windows, which are
# common on the prairies with a quantised DEM -- set those to the flat case.
srr <- function(z, w) {
  m  <- matrix(1, w, w)
  hi <- focal(z, m, "max",  na.rm = TRUE)
  lo <- focal(z, m, "min",  na.rm = TRUE)
  mn <- focal(z, m, "mean", na.rm = TRUE)
  ifel(hi == lo, 0.5, (mn - lo) / (hi - lo))
}


#  write_rast_tagged(tri(dem, w), "./processeddatasets", dem_filepath, paste0("terr_tri-", tag), "modelling-extent")
#  wr(tri(dem, w), paste0("terr_tri-", tag))
#  wr(srr(dem, w), paste0("terr_srr-", tag))
#  wr(tpi(dem, w), paste0("terr_hsp-", tag))
}

# ---- curvature + heat load (spatialEco / MultiscaleDTM) ----
# spatialEco is Evans' own R port of the toolbox O&R cite, so closest parity.
# Not always on CRAN: remotes::install_github("jeffreyevans/spatialEco")
#
# wr(spatialEco::curvature(dem, type = "mcnab"),   "terr_curv-mcnab")
# wr(spatialEco::curvature(dem, type = "bolstad"), "terr_curv-bolstad")
#
# hli() derives one latitude for the whole raster by default. The Prairies
# Ecozone spans ~49-55 N, so tile the extent by latitude band or pass a
# per-cell latitude before trusting the result.
# wr(spatialEco::hli(dem), "terr_hli")