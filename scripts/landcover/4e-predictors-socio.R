# Purpose:
#   Distance, capacity and net-return predictors for the RF conversion risk
#   model: Canfax net returns (2005-2009, 2020-2024, standardised within
#   period), distance to elevators, large elevators, roads, rails, population
#   and regional centres (2011, 2021), capacity of nearest elevator.
#   Ends with a collinearity/range check and figures.
#
# Inputs:
#   - Modelling extent and study area raster templates (me_r, sa_r)
#   - Canfax polygons, elevator points, road/rail rasters, popctr polygons,
#     all clipped to extent by 4-prepare-spatial-inputs.R
#
# Outputs:
#   - Canfax rasters (raw and _z) to predictors dir, tagged {res_m}m
#   - dist_to_* and elevator capacity rasters to their source dirs,
#     tagged "modelling-extent"
#   - Predictor maps and correlation matrix to fig_dir
#
# Notes:
#   - Regional centres = PCTYPE 1, 6 (CMA/CA cores, >=10,000 in core);
#     distinct from popctr (Spearman ~0.54 within year)
#   - Same predictor across years is near-identical (0.94 popctr,
#     0.97 regctr) - don't use both years of one variable
#   - Large elevator threshold 50,000 t (n = 6), at the curvature change in
#     sorted CAP_T; Spearman 0.32 with distance to any elevator
#   - Canfax z-scored within period: 62% of raw 2024 cells above the 2009 max
#     where RF cannot extrapolate, 2.4% after. Assumes conversion responds to
#     relative, not absolute, profitability
#   - Distances in metres

source(file.path("scripts/landcover", "0a-setup.R"))

# ---- Load files ----
# Modelling extent
me_r <- rast(file.path(sa_filepath, me_raster_filename))
sa <- st_read(file.path(sa_filepath, sa_filename))
sa_r <- rast(file.path(sa_filepath, sa_raster_filename))
r <- nrow(sa_r) / ncol(sa_r)  # needed for plots
fig_dir <- "./output/landcover/figures"

# Ecodistricts
ecodistricts <- st_read(tagged_path("./processeddatasets", ecodistricts_filepath,
                                    ecodistricts_filename, "modelling-extent"))
# Canfax results
canfax2009 <- st_read(tagged_path("./processeddatasets", canfax_filepath, 
                                  canfax2009_filename, "modelling-extent"))
canfax2024 <- st_read(tagged_path("./processeddatasets", canfax_filepath, 
                                  canfax2024_filename, "modelling-extent"))

# Process elevators
elevator <- st_read(tagged_path("./processeddatasets", elevator_filepath, 
                              elevator_filename, "modelling-extent"))
# Roads
roads <- rast(tagged_path("./processeddatasets", road_filepath,
                              paste0(road_filename, ".tif"), "modelling-extent"))
# Rails
rails <- rast(tagged_path("./processeddatasets", rail_filepath,
                          paste0(rail_filename, ".tif"), "modelling-extent"))

# Population centres
popctr2011 <- st_read(tagged_path("./processeddatasets", popctr_filepath, 
                              paste0(popctr2011_filename, ".shp"), "modelling-extent"))
popctr2021 <- st_read(tagged_path("./processeddatasets", popctr_filepath, 
                                  paste0(popctr2021_filename, ".shp"), "modelling-extent"))

# ---- Rasterize Canfax results ----
# 2005 - 2009
canfax2009_v <- vect(canfax2009)
canfax2009_r <- rasterize(canfax2009_v, me_r, field = "nt_b_25")
canfax2009_r_sa <- mask(crop(canfax2009_r, sa_r), sa_r)
png(file.path(fig_dir, "canfax2005-2009.png"), width = 2400, height = round(2400 * r), res = 300)
  plot(canfax2009_r_sa, main = "Net returns of cropping 2005 - 2009 (CAD)", axes=FALSE, box=FALSE)
dev.off()

# 2020 - 2024
canfax2024_v <- vect(canfax2024)
canfax2024_r <- rasterize(canfax2024_v, me_r, field = "nt_b_25")
canfax2024_r_sa <- mask(crop(canfax2024_r, sa_r), sa_r)
png(file.path(fig_dir, "canfax2020-2024.png"), width = 2400, height = round(2400 * r), res = 300)
  plot(canfax2024_r_sa, main = "Net returns of cropping 2020 - 2024 (CAD)", axes=FALSE, box=FALSE)
dev.off()

### Standardize ###
# Canfax net returns standardised within period: 2005-2009 is entirely negative
# (95th pct -132) while 2020-2024 is mostly positive (median 814), so 62% of raw
# 2024 cells fell above the 2009 max, where the RF cannot extrapolate.
# Distributions have similar spread (IQR 1397 vs 1763), so z-scoring maps each
# period onto relative position; 2.4% of cells remain above the 2009 z-max
# (2.19 vs 2.76), 2024 being slightly more right-skewed. Assumes conversion
# responds to relative, not absolute, profitability.
z <- function(r) (r - global(r, "mean", na.rm = TRUE)[[1]]) /
  global(r, "sd", na.rm = TRUE)[[1]]

# standardize across modelling extent
canfax2009_z <- z(canfax2009_r)
canfax2024_z <- z(canfax2024_r)

# standardize across study area
canfax2009_z_sa <- z(canfax2009_r_sa)
canfax2024_z_sa <- z(canfax2024_r_sa)

png(file.path(fig_dir, "canfax2005-2009_z.png"), width = 2400, height = round(2400 * r), res = 300)
  plot(canfax2009_z_sa, main = "Net returns of cropping 2005 - 2009 (z-score)", axes=FALSE, box=FALSE)
dev.off()

png(file.path(fig_dir, "canfax2020-2024_z.png"), width = 2400, height = round(2400 * r), res = 300)
  plot(canfax2024_z_sa, main = "Net returns of cropping 2020 - 2024 (z-score)", axes=FALSE, box=FALSE)
dev.off()


# ---- Distance to process elevators ----
# Run this for plots
ev <- vect(elevator)
ev_sa <- ev[!is.na(terra::extract(sa_r, ev)[, 2]), ]

### Distance to all elevators ###
elev_d <- distance(me_r, vect(elevator))   # metres
elev_d <- mask(elev_d, me_r)

# Plot
elev_d_sa <- mask(crop(elev_d, sa_r), sa_r)
#plot(elev_d_sa / 1000, main = "Distance to process elevator (km)", axes=FALSE, box=FALSE)
lim <- quantile(spatSample(elev_d_sa / 1000, 2e4, method = "regular",
                           na.rm = TRUE)[, 1], c(.02, .98))
# clamp() piles the clipped tails onto lim exactly; plot() leaves cells at the
# range boundary outside the last colour interval and draws them as NA, so pad
# range by 0.1% of the span
rng <- lim + c(-1, 1) * diff(lim) * 1e-3
png(file.path(fig_dir, "dist_to_elev.png"), width = 2400, height = round(2400 * r), res = 300)
  plot(clamp(elev_d_sa / 1000, lim[1], lim[2], values = TRUE),
     range = rng, col = hcl.colors(100, "viridis"),
     main = "Distance to process elevator (km)", axes = FALSE, box = FALSE)
  points(ev_sa, pch = 21, bg = "white", col = "black", lwd = 0.8,
       cex = 0.4 + 2 * sqrt(ev_sa$CAP_T / max(ev_sa$CAP_T, na.rm = TRUE)))
dev.off()

### Distance to large elevators ###
# Check elevator capacity values
sort(elevator$CAP_T)
# [1]    200    370    650    800    820   1024   1465   1500   1950   2040   2100   2400   2500   2780   2850   2950   3000   3100   4150   4170   4500   4870   5000   5000   5254   5860   6000   6600   7694
# [30]   7900   8250   8400   8500   8500   8700   9390  10000  11030  11750  12000  12125  13000  13230  15390  16800  24000  31570  32000  35000  36252  38300  46800  48990  62700  63500  66470  89330  89640
# [59] 107380
plot(sort(elevator$CAP_T))
# Identify 90th percentile
quantile(elevator$CAP_T, 0.9, na.rm = TRUE)
# 90% 
# 51732
# Adjacent values are 48990 and 62700 so use 50000 as threshold

# Threshold from the 90th percentile of CAP_T; no external size class available,
# so this is a distributional cut, not a functional one. Fill in value + n here.
cap_thresh <- 50000 #quantile(elevator$CAP_T, 0.9, na.rm = TRUE)
elev_large <- elevator[!is.na(elevator$CAP_T) & elevator$CAP_T >= cap_thresh, ]
elev_large_d <- mask(distance(me_r, vect(elev_large)), me_r)

# Plot
elev_large_d_sa <- mask(crop(elev_large_d, sa_r), sa_r)
#plot(elev_large_d_sa / 1000, main = "Distance to large process elevator (km)", axes=FALSE, box=FALSE)

# elevator labels
lab <- sub(" (Company|Canada|Inc\\.?|Ltd\\.?|Limited).*", "", elev_large$LICENSEE)
p <- rep(3, nrow(elev_large))
p[c(2, 3, 5)] <- c(1, 3, 4)   # one below, one above

lim <- quantile(spatSample(elev_large_d_sa / 1000, 2e4, method = "regular",
                           na.rm = TRUE)[, 1], c(.02, .98))
# clamp() piles the clipped tails onto lim exactly; plot() leaves cells at the
# range boundary outside the last colour interval and draws them as NA, so pad
# range by 0.1% of the span
rng <- lim + c(-1, 1) * diff(lim) * 1e-3
png(file.path(fig_dir, "dist_to_elev_large.png"), width = 2400, height = round(2400 * r), res = 300)
  plot(clamp(elev_large_d_sa / 1000, lim[1], lim[2], values = TRUE),
       range = rng, col = hcl.colors(100, "viridis"),
       main = "Distance to large process elevator (km)", axes = FALSE, box = FALSE)
  points(vect(elev_large), pch = 21, bg = "white", col = "black", lwd = 0.8,
         cex = 0.4 + 2 * sqrt(elev_large$CAP_T / max(elev_large$CAP_T, na.rm = TRUE)))
  text(vect(elev_large), labels = lab, pos = p, cex = 0.45,
       halo = TRUE, hc = "white", hw = 0.15)
dev.off()


# does the threshold matter? - yes!
# th <- quantile(elevator$CAP_T, c(0.75, 0.9, 0.95), na.rm = TRUE)
# d <- rast(lapply(th, function(t) {
#   distance(me_r, vect(elevator[!is.na(elevator$CAP_T) & elevator$CAP_T >= t, ]))
# }))
# s <- spatSample(c(d, elev_d), 1e5, method = "regular", na.rm = TRUE, as.df = TRUE)
# round(cor(s, method = "spearman"), 3)


### Capacity of the nearest elevator: voronoi cells are exactly the nearest-point ###
# assignment used by distance(), so elev_cap and elev_d refer to the same elevator
elev_vor <- voronoi(vect(elevator), bnd = ext(vect(elevator)) + 1e5)
elev_cap <- mask(rasterize(elev_vor, me_r, field = "CAP_T"), me_r)

# Plot
elev_cap_sa <- mask(crop(elev_cap, sa_r), sa_r)
png(file.path(fig_dir, "cap_nearest_elev.png"), width = 2400, height = round(2400 * r), res = 300)
  plot(elev_cap_sa / 1000, main = "Capacity of nearest elevator (kt)", axes=FALSE, box=FALSE)
  points(ev_sa, pch = 21, bg = "white", col = "black", lwd = 0.8,
         cex = 0.4 + 2 * sqrt(ev_sa$CAP_T / max(ev_sa$CAP_T, na.rm = TRUE)))
dev.off()


# ---- Distance to roads ----
roads_for_dist <- ifel(roads == 1, 1, NA)
dist_to_road  <- distance(roads_for_dist)
dist_to_road  <- mask(dist_to_road, me_r)
# plot(dist_to_road / 1000, main = "Distance to roads (km)")
# Plot
dist_to_road_sa <- mask(crop(dist_to_road, sa_r), sa_r)
lim <- quantile(spatSample(dist_to_road_sa / 1000, 2e4, method = "regular",
                           na.rm = TRUE)[, 1], c(.02, .98))
# clamp() piles the clipped tails onto lim exactly; plot() leaves cells at the
# range boundary outside the last colour interval and draws them as NA, so pad
# range by 0.1% of the span
rng <- lim + c(-1, 1) * diff(lim) * 1e-3
png(file.path(fig_dir, "dist_to_road.png"), width = 2400, height = round(2400 * r), res = 300)
plot(clamp(dist_to_road_sa / 1000, lim[1], lim[2], values = TRUE),
     range = rng, col = hcl.colors(100, "viridis"),
     main = "Distance to roads (km)", axes = FALSE, box = FALSE)
dev.off()

v <- spatSample(dist_to_road_sa / 1000, 1e5, method = "regular", na.rm = TRUE)[, 1]
hist(v, breaks = 50, main = "", xlab = "Distance to road (km)")

# ---- Distance to railways ----
rails_for_dist <- ifel(rails == 1, 1, NA)
dist_to_rail  <- distance(rails_for_dist)
dist_to_rail  <- mask(dist_to_rail, me_r)
# plot(dist_to_rail / 1000, main = "Distance to railways (km)")
# Plot
dist_to_rail_sa <- mask(crop(dist_to_rail, sa_r), sa_r)
lim <- quantile(spatSample(dist_to_rail_sa / 1000, 2e4, method = "regular",
                           na.rm = TRUE)[, 1], c(.02, .98))
# clamp() piles the clipped tails onto lim exactly; plot() leaves cells at the
# range boundary outside the last colour interval and draws them as NA, so pad
# range by 0.1% of the span
rng <- lim + c(-1, 1) * diff(lim) * 1e-3
png(file.path(fig_dir, "dist_to_rail.png"), width = 2400, height = round(2400 * r), res = 300)
plot(clamp(dist_to_rail_sa / 1000, lim[1], lim[2], values = TRUE),
     range = rng, col = hcl.colors(100, "viridis"),
     main = "Distance to railways (km)", axes = FALSE, box = FALSE)
dev.off()


# ---- Population centre distance predictors ----
# POPCTR (>=1,000 people, >=400/km2), 2011: 175 polygons in extent.
# PCTYPE: 1 core (28), 2 fringe (30), 4 outside any CMA/CA (112), 6 secondary core (5).
# PCCLASS by population: 2 small (1,000-29,999), 3 medium (30,000-99,999), 4 large (100,000+).

# Subset population centres that serve as the core or secondary core of a 
# census metropolitan area (CMA) or census agglomeration (CA).
# A CMA needs a total population of at least 100,000 with 50,000 or more in the core
# A CA needs a core of at least 10,000.
pc_cma_ca_cores2011 <- popctr2011[popctr2011$PCTYPE %in% c(1, 6), ]
pc_cma_ca_cores2021 <- popctr2021[popctr2021$PCTYPE %in% c(1, 6), ]

# Distance to nearest population centre and distance to nearest regional centre
### 2011 ###
pc_d_all_2011 <- distance(me_r, vect(popctr2011), rasterize = TRUE, unit = "m")
pc_d_all_2011 <- mask(pc_d_all_2011, me_r)
#plot(pc_d_all_2011/1000, main = "Distance to nearest population centre (km)")
#plot(st_geometry(popctr2011), border = "white", lwd = 0.4, add = TRUE)
# Plot
pc_d_all_2011_sa <- mask(crop(pc_d_all_2011, sa_r), sa_r)
lim <- quantile(spatSample(pc_d_all_2011_sa / 1000, 2e4, method = "regular",
                           na.rm = TRUE)[, 1], c(.02, .98))
# clamp() piles the clipped tails onto lim exactly; plot() leaves cells at the
# range boundary outside the last colour interval and draws them as NA, so pad
# range by 0.1% of the span
rng <- lim + c(-1, 1) * diff(lim) * 1e-3
png(file.path(fig_dir, "dist_popctr_2011.png"), width = 2400, height = round(2400 * r), res = 300)
plot(clamp(pc_d_all_2011_sa / 1000, lim[1], lim[2], values = TRUE),
     range = rng, col = hcl.colors(100, "viridis"),
     main = "Distance to population centre 2011 (km)", axes = FALSE, box = FALSE)
plot(st_geometry(popctr2011), border = "white", lwd = 0.4, add = TRUE)
dev.off()

pc_d_core_2011 <- distance(me_r, vect(pc_cma_ca_cores2011), rasterize = TRUE, unit = "m")
pc_d_core_2011 <- mask(pc_d_core_2011, me_r)
#plot(pc_d_core_2011/1000, main = "Distance to nearest regional centre (km)")
#plot(st_geometry(pc_cma_ca_cores2011), border = "white", lwd = 0.4, add = TRUE)
# Plot
pc_d_core_2011_sa <- mask(crop(pc_d_core_2011, sa_r), sa_r)
lim <- quantile(spatSample(pc_d_core_2011_sa / 1000, 2e4, method = "regular",
                           na.rm = TRUE)[, 1], c(.02, .98))
# clamp() piles the clipped tails onto lim exactly; plot() leaves cells at the
# range boundary outside the last colour interval and draws them as NA, so pad
# range by 0.1% of the span
rng <- lim + c(-1, 1) * diff(lim) * 1e-3
png(file.path(fig_dir, "dist_regctr_2011.png"), width = 2400, height = round(2400 * r), res = 300)
  plot(clamp(pc_d_core_2011_sa / 1000, lim[1], lim[2], values = TRUE),
     range = rng, col = hcl.colors(100, "viridis"),
     main = "Distance to regional centre 2011 (km)", axes = FALSE, box = FALSE)
  plot(st_geometry(pc_cma_ca_cores2011), border = "white", lwd = 0.4, add = TRUE)
dev.off()

### 2021 ###
pc_d_all_2021 <- distance(me_r, vect(popctr2021), rasterize = TRUE, unit = "m")
pc_d_all_2021 <- mask(pc_d_all_2021, me_r)
# plot(pc_d_all_2021/1000, main = "Distance to nearest population centre (km)")
# plot(st_geometry(popctr2021), border = "white", lwd = 0.4, add = TRUE)
# Plot
pc_d_all_2021_sa <- mask(crop(pc_d_all_2021, sa_r), sa_r)
lim <- quantile(spatSample(pc_d_all_2021_sa / 1000, 2e4, method = "regular",
                           na.rm = TRUE)[, 1], c(.02, .98))
# clamp() piles the clipped tails onto lim exactly; plot() leaves cells at the
# range boundary outside the last colour interval and draws them as NA, so pad
# range by 0.1% of the span
rng <- lim + c(-1, 1) * diff(lim) * 1e-3
png(file.path(fig_dir, "dist_popctr_2021.png"), width = 2400, height = round(2400 * r), res = 300)
plot(clamp(pc_d_all_2021_sa / 1000, lim[1], lim[2], values = TRUE),
     range = rng, col = hcl.colors(100, "viridis"),
     main = "Distance to population centre 2021 (km)", axes = FALSE, box = FALSE)
plot(st_geometry(popctr2021), border = "white", lwd = 0.4, add = TRUE)
dev.off()

pc_d_core_2021 <- distance(me_r, vect(pc_cma_ca_cores2021), rasterize = TRUE, unit = "m")
pc_d_core_2021 <- mask(pc_d_core_2021, me_r)
# plot(pc_d_core_2021/1000, main = "Distance to nearest regional centre (km)")
# plot(st_geometry(pc_cma_ca_cores2021), border = "white", lwd = 0.4, add = TRUE)
# Plot
pc_d_core_2021_sa <- mask(crop(pc_d_core_2021, sa_r), sa_r)
lim <- quantile(spatSample(pc_d_core_2021_sa / 1000, 2e4, method = "regular",
                           na.rm = TRUE)[, 1], c(.02, .98))
# clamp() piles the clipped tails onto lim exactly; plot() leaves cells at the
# range boundary outside the last colour interval and draws them as NA, so pad
# range by 0.1% of the span
rng <- lim + c(-1, 1) * diff(lim) * 1e-3
png(file.path(fig_dir, "dist_regctr_2021.png"), width = 2400, height = round(2400 * r), res = 300)
plot(clamp(pc_d_core_2021_sa / 1000, lim[1], lim[2], values = TRUE),
     range = rng, col = hcl.colors(100, "viridis"),
     main = "Distance to regional centre 2021 (km)", axes = FALSE, box = FALSE)
plot(st_geometry(pc_cma_ca_cores2021), border = "white", lwd = 0.4, add = TRUE)
dev.off()

# ---- Save rasters ----
write_rast_tagged(canfax2009_r, "./processeddatasets", canfax_filepath, "canfax_nt_b_25_2005-2009.tif", "modelling-extent")
write_rast_tagged(canfax2024_r, "./processeddatasets", canfax_filepath, "canfax_nt_b_25_2020-2024.tif", "modelling-extent")
write_rast_tagged(canfax2009_z, "./processeddatasets", canfax_filepath, "canfax_nt_b_25_2005-2009_z.tif", "modelling-extent")
write_rast_tagged(canfax2024_z, "./processeddatasets", canfax_filepath, "canfax_nt_b_25_2020-2024_z.tif", "modelling-extent")
write_rast_tagged(canfax2009_z_sa, "./processeddatasets", canfax_filepath, "canfax_nt_b_25_2005-2009_z.tif", "study-area")
write_rast_tagged(canfax2024_z_sa, "./processeddatasets", canfax_filepath, "canfax_nt_b_25_2020-2024_z.tif", "study-area")
write_rast_tagged(elev_d, "./processeddatasets", elevator_filepath, "dist_to_elevator.tif", "modelling-extent")
write_rast_tagged(elev_large_d, "./processeddatasets", elevator_filepath, "dist_to_elevator_large.tif", "modelling-extent")
write_rast_tagged(elev_cap, "./processeddatasets", elevator_filepath, "cap_nearest_elevator.tif", "modelling-extent")
write_rast_tagged(dist_to_road,  "./processeddatasets", road_filepath, "dist_to_road.tif", "modelling-extent")
write_rast_tagged(dist_to_rail,  "./processeddatasets", rail_filepath, "dist_to_rail.tif", "modelling-extent")
write_rast_tagged(pc_d_all_2011,  "./processeddatasets", popctr_filepath, "dist_to_popctr2011.tif", "modelling-extent")
write_rast_tagged(pc_d_all_2021,  "./processeddatasets", popctr_filepath, "dist_to_popctr2021.tif", "modelling-extent")
write_rast_tagged(pc_d_core_2011,  "./processeddatasets", popctr_filepath, "dist_to_regctr2011.tif", "modelling-extent")
write_rast_tagged(pc_d_core_2021,  "./processeddatasets", popctr_filepath, "dist_to_regctr2021.tif", "modelling-extent")

# ---- Correlation ----
# Read files back in and check for correlation and ranges across study area extent
read_sa <- function(...) {
  r <- rast(tagged_path(...))
  mask(crop(r, sa_r), sa_r)
}

canfax2009_r <- read_sa("./processeddatasets", canfax_filepath, "canfax_nt_b_25_2005-2009.tif", "modelling-extent")
canfax2024_r <- read_sa("./processeddatasets", canfax_filepath, "canfax_nt_b_25_2020-2024.tif", "modelling-extent")
elev_d <- read_sa("./processeddatasets", elevator_filepath, "dist_to_elevator.tif", "modelling-extent")
elev_large_d <- read_sa("./processeddatasets", elevator_filepath, "dist_to_elevator_large.tif", "modelling-extent")
elev_cap <- read_sa("./processeddatasets", elevator_filepath, "cap_nearest_elevator.tif", "modelling-extent")
dist_to_road <- read_sa("./processeddatasets", road_filepath, "dist_to_road.tif", "modelling-extent")
dist_to_rail <- read_sa("./processeddatasets", rail_filepath, "dist_to_rail.tif", "modelling-extent")
pc_d_all_2011 <- read_sa("./processeddatasets", popctr_filepath, "dist_to_popctr2011.tif", "modelling-extent")
pc_d_all_2021 <- read_sa("./processeddatasets", popctr_filepath, "dist_to_popctr2021.tif", "modelling-extent")
pc_d_core_2011 <- read_sa("./processeddatasets", popctr_filepath, "dist_to_regctr2011.tif", "modelling-extent")
pc_d_core_2021 <- read_sa("./processeddatasets", popctr_filepath, "dist_to_regctr2021.tif", "modelling-extent")


socio <- c(canfax2009_r, elev_d,
           dist_to_road, dist_to_rail,
           pc_d_all_2011, pc_d_core_2011)
names(socio) <- c("canfax2009", "d_elev", 
                  "d_road", "d_rail", "popctr2011",
                  "regctr2011")

s <- spatSample(socio, 1e4, method = "regular", na.rm = TRUE, as.df = TRUE)
m <- round(cor(s, method = "spearman"), 3)


o <- hclust(as.dist(1 - m))$order
m <- m[o, o]
m[lower.tri(m)] <- NA

png(file.path(fig_dir, "socio_cor.png"), width = 1800, height = 1800, res = 300)
par(mar = c(8, 8, 1, 3))
image(seq_len(ncol(m)), seq_len(nrow(m)), t(m)[, nrow(m):1],
      col = rev(hcl.colors(200, "RdBu")), zlim = c(-1, 1),
      axes = FALSE, xlab = "", ylab = "")
axis(1, seq_len(ncol(m)), colnames(m), las = 2, tick = FALSE)
axis(2, seq_len(nrow(m)), rev(rownames(m)), las = 2, tick = FALSE)
text(rep(seq_len(ncol(m)), each = nrow(m)),
     rep(nrow(m):1, ncol(m)),
     round(t(m), 2), cex = 0.6)
dev.off()


# Check range of dynamic predictors
### population and regional centres ###
pc <- c(pc_d_all_2011, pc_d_all_2021, pc_d_core_2011, pc_d_core_2021)
names(pc) <- c("popctr2011", "popctr2021", "regctr2011", "regctr2021")
# range: 2021 values outside the 2011 range could be a problem for RF model
minmax(pc)
# popctr2011 popctr2021 regctr2011 regctr2021
# min        0.0        0.0        0.0        0.0
# max   118994.2   118994.4   187207.3   190408.7

### canfax ###
cf <- c(canfax2009_r, canfax2024_r)
names(cf) <- c("canfax2009", "canfax2024")
# range: 2021 values outside the 2011 range could be a problem for RF model
minmax(cf)
# canfax2009 canfax2024
# min -3895.5894  -1284.636
# max   341.7303   3743.167
global(cf[[2]] > minmax(cf)[2, 1], "mean", na.rm = TRUE)
