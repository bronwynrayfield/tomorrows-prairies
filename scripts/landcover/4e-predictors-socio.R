# =============================================================================
# Purpose:
#   Socio-economic predictors for the RF conversion risk model: 
#   Canfax net returns (2005-2009, 2020-2024, standardised within
#   period), distance to elevators, large elevators, roads, rails, population
#   and regional centres (2011, 2021), population density (2011, 2021), capacity of nearest elevator.
#   Ends with a collinearity/range check and figures.
#
# Inputs:
#   - Modelling extent and study area raster templates (me_r, sa_r)
#   - Canfax polygons, elevator points, road/rail rasters, popctr polygons,
#     all clipped to extent by 4-prepare-spatial-inputs.R
#
# Outputs:
#   - Predictor maps and correlation matrix to fig_dir
#
# Notes:
#   - Regional centres = PCTYPE 1, 6 (CMA/CA cores, >=10,000 in core);
#     distinct from popctr (Spearman ~0.54 within year)
#   - Same predictor across years is near-identical (0.94 popctr,
#     0.97 regctr) - don't use both years of one variable
#   - Large elevator threshold 50,000 t (n = 6), at 90th percentile; 
#     Spearman 0.32 with distance to any elevator
#   - Canfax z-scored within period: 62% of raw 2024 cells above the 2009 max
#     where RF cannot extrapolate, 2.4% after. Assumes conversion responds to
#     relative, not absolute, profitability
#   - Population density extrapolation not an issue
#   - Distances in metres
# =============================================================================
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
# Water
water2010 <- rast(file.path("./processeddatasets", aafc_filepath, "binary", "AAFC2010_water.tif"))
water2020 <- rast(file.path("./processeddatasets", aafc_filepath, "binary", "AAFC2020_water.tif"))

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

# Dissemination blocks
db2011 <- st_read(tagged_path("./processeddatasets", popdens_filepath,
                              paste0(db2011_filename, ".shp"), "modelling-extent"))
db2021 <- st_read(tagged_path("./processeddatasets", popdens_filepath,
                              paste0(db2021_filename, ".shp"), "modelling-extent"))
# Geographic Attribute Files: DB population and land area (km2, land only)
f11 <- file.path("./existingdatasets", capitalize_first(popdens_filepath), paste0(gaf2011_filename, ".txt"))
gaf2011 <- read.fwf(f11, widths = c(10, 8, 8, 8, 13, 1), fileEncoding = "latin1",
                    col.names = c("DBuid", "DBpop2011", "DBtdwell2011", "DBurdwell2011", "DBarea2011", "DBir2011"),
                    colClasses = c("character", "numeric", "numeric", "numeric", "numeric", "character"),
                    na.strings = c(".", "..", "..."), strip.white = TRUE)
f21 <- file.path("./existingdatasets", capitalize_first(popdens_filepath), paste0(gaf2021_filename, ".csv"))
keep <- c(DBUID_IDIDU = "character", DBPOP2021_IDPOP2021 = "numeric", DBAREA2021_IDSUP2021 = "numeric", DBIR2021_IDRI2021 = "logical")
nm <- names(read.csv(f21, nrows = 1, fileEncoding = "latin1"))
gaf2021 <- read.csv(f21, fileEncoding = "latin1", colClasses = ifelse(nm %in% names(keep), keep[nm], "NULL"))


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


# ---- Population density ----
# Persons per km2 of land within a square window: DB population spread evenly
# over each DB's land cells, then window sum of persons / window sum of land area.
# Land area from the Geographic Attribute File GAF (excludes water), not st_area(). Water-only DBs (area 0)
# and incompletely enumerated reserves (2011 refusal flag) set to NA.
i <- match(db2011$DBUID, as.character(gaf2011$DBuid))
sum(is.na(i))
db2011$dens <- ifelse(gaf2011$DBir2011[i] == "T" | gaf2011$DBarea2011[i] == 0, NA,
                      gaf2011$DBpop2011[i] / gaf2011$DBarea2011[i])
i <- match(db2021$DBUID, as.character(gaf2021$DBUID_IDIDU))
sum(is.na(i))
db2021$dens <- ifelse(gaf2021$DBIR2021_IDRI2021[i] | gaf2021$DBAREA2021_IDSUP2021[i] == 0, NA,
                      gaf2021$DBPOP2021_IDPOP2021[i] / gaf2021$DBAREA2021_IDSUP2021[i])

# Window side = sqrt(median land area of rural DBs, i.e. not touching a 2011 POPCTR),
# so each window spans more than one block.
rural <- lengths(st_intersects(db2011, popctr2011)) == 0
side_m <- sqrt(median(gaf2011$DBarea2011[match(db2011$DBUID[rural], as.character(gaf2011$DBuid))],
                      na.rm = TRUE) * 1e6)
fact <- 9
# Window in quarter-section units (27 cells = 810 m, as for terrain predictors): smallest odd number
# of quarters >= side of the median rural DB (3.98 km2 -> 1995 m -> 2.46 -> 3 quarters = 81 cells, 2430 m)
q <- 27
n_q <- 2 * ceiling((side_m / (q * res(me_r)[1]) - 1) / 2) + 1
w <- n_q * q / fact   # 9 cells at fact = 9
cell_km2 <- prod(res(me_r)) / 1e6

pop_dens <- function(db, water) {
  d    <- mask(rasterize(vect(db), me_r, field = "dens"), water, maskvalues = 1)   # GAF density is per land km2; drop water cells
  pop  <- aggregate(d * cell_km2, fact, sum, na.rm = TRUE)    # persons per 270 m cell
  land <- aggregate(!is.na(d), fact, sum) * cell_km2          # land km2 per 270 m cell
  dens <- focal(pop, w, sum, na.rm = TRUE) / focal(land, w, sum, na.rm = TRUE)
  mask(mask(resample(dens, me_r, "bilinear"), me_r), water, maskvalues = 1)
}
pd_2011 <- pop_dens(db2011, water2010)
pd_2021 <- pop_dens(db2021, water2020)

# Study area
pd_2011_sa <- mask(crop(pd_2011, sa_r), sa_r)
pd_2021_sa <- mask(crop(pd_2021, sa_r), sa_r)

# Check extrapolation
# How many cells in 2021 are outside the range of values in 2011?
mx <- global(pd_2011_sa, "max", na.rm = TRUE)[[1]]
global(pd_2021_sa > mx, "sum", na.rm = TRUE) 
#           sum                            
# focal_sum 2663
# only 2663 cells in 2021 exceed the max in 2011 (none overlapped with grass when I checked)

# Plot (raw)
lim <- quantile(spatSample(pd_2011_sa, 2e4, method = "regular", na.rm = TRUE)[, 1], c(.02, .98))
rng <- lim + c(-1, 1) * diff(lim) * 1e-3
png(file.path(fig_dir, "pop_density_2011.png"), width = 2400, height = round(2400 * r), res = 300)
plot(clamp(pd_2011_sa, lim[1], lim[2], values = TRUE),
     range = rng, col = hcl.colors(100, "viridis"),
     main = "Population density 2011 (persons/km2)", axes = FALSE, box = FALSE)
dev.off()

lim <- quantile(spatSample(pd_2021_sa, 2e4, method = "regular", na.rm = TRUE)[, 1], c(.02, .98))
rng <- lim + c(-1, 1) * diff(lim) * 1e-3
png(file.path(fig_dir, "pop_density_2021.png"), width = 2400, height = round(2400 * r), res = 300)
plot(clamp(pd_2021_sa, lim[1], lim[2], values = TRUE),
     range = rng, col = hcl.colors(100, "viridis"),
     main = "Population density 2021 (persons/km2)", axes = FALSE, box = FALSE)
dev.off()


# Plot (log10 scale for display only; RF uses raw density)
lim <- quantile(spatSample(log10(pd_2011_sa + 1), 2e4, method = "regular", na.rm = TRUE)[, 1], c(.02, .98))
rng <- lim + c(-1, 1) * diff(lim) * 1e-3
png(file.path(fig_dir, "pop_density_2011_log10.png"), width = 2400, height = round(2400 * r), res = 300)
plot(clamp(log10(pd_2011_sa + 1), lim[1], lim[2], values = TRUE),
     range = rng, col = hcl.colors(100, "viridis"),
     main = "Population density 2011 (log10 persons/km2 + 1)", axes = FALSE, box = FALSE)
dev.off()

lim <- quantile(spatSample(log10(pd_2021_sa + 1), 2e4, method = "regular", na.rm = TRUE)[, 1], c(.02, .98))
rng <- lim + c(-1, 1) * diff(lim) * 1e-3
png(file.path(fig_dir, "pop_density_2021_log10.png"), width = 2400, height = round(2400 * r), res = 300)
plot(clamp(log10(pd_2021_sa + 1), lim[1], lim[2], values = TRUE),
     range = rng, col = hcl.colors(100, "viridis"),
     main = "Population density 2021 (log10 persons/km2 + 1)", axes = FALSE, box = FALSE)
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
write_rast_tagged(pd_2011, "./processeddatasets", popdens_filepath, "pop_density2011.tif", "modelling-extent")
write_rast_tagged(pd_2021, "./processeddatasets", popdens_filepath, "pop_density2021.tif", "modelling-extent")

# ---- Correlation ----

# Read files back in and check for correlation and ranges across study area extent
# If reading modelling-extent, need to mask and crop
# read_sa <- function(...) {
#   r <- rast(tagged_path(...))
#   mask(crop(r, sa_r), sa_r)
# }
# 
# canfax2009_r <- read_sa("./processeddatasets", canfax_filepath, "canfax_nt_b_25_2005-2009.tif", "modelling-extent")
# elev_d <- read_sa("./processeddatasets", elevator_filepath, "dist_to_elevator.tif", "modelling-extent")
# dist_to_road <- read_sa("./processeddatasets", road_filepath, "dist_to_road.tif", "modelling-extent")
# dist_to_rail <- read_sa("./processeddatasets", rail_filepath, "dist_to_rail.tif", "modelling-extent")
# pc_d_all_2011 <- read_sa("./processeddatasets", popctr_filepath, "dist_to_popctr2011.tif", "modelling-extent")
# pc_d_core_2011 <- read_sa("./processeddatasets", popctr_filepath, "dist_to_regctr2011.tif", "modelling-extent")
# pd_2011 <- read_sa("./processeddatasets", popdens_filepath, "pop_density2011.tif", "modelling-extent")
# canfax2024_r <- read_sa("./processeddatasets", canfax_filepath, "canfax_nt_b_25_2020-2024.tif", "modelling-extent")
# pc_d_all_2021 <- read_sa("./processeddatasets", popctr_filepath, "dist_to_popctr2021.tif", "modelling-extent")
# pc_d_core_2021 <- read_sa("./processeddatasets", popctr_filepath, "dist_to_regctr2021.tif", "modelling-extent")
# pd_2021 <- read_sa("./processeddatasets", popdens_filepath, "pop_density2021.tif", "modelling-extent")
# elev_large_d <- read_sa("./processeddatasets", elevator_filepath, "dist_to_elevator_large.tif", "modelling-extent")
# elev_cap <- read_sa("./processeddatasets", elevator_filepath, "cap_nearest_elevator.tif", "modelling-extent")

# Read study area files directly
canfax2009_r <- rast(tagged_path("./processeddatasets", file.path(predictors_filepath, "30m"), "net_return.tif", "study-area"))
elev_d <- rast(tagged_path("./processeddatasets", file.path(predictors_filepath, "30m"), "dist_to_elevator.tif", "study-area"))
dist_to_road <- rast(tagged_path("./processeddatasets", file.path(predictors_filepath, "30m"), "dist_to_road.tif", "study-area"))
dist_to_rail <- rast(tagged_path("./processeddatasets", file.path(predictors_filepath, "30m"), "dist_to_rail.tif", "study-area"))
pc_d_all_2011 <- rast(tagged_path("./processeddatasets", file.path(predictors_filepath, "30m"), "dist_to_popctr2011.tif", "study-area"))
pc_d_core_2011 <- rast(tagged_path("./processeddatasets", file.path(predictors_filepath, "30m"), "dist_to_regctr2011.tif", "study-area"))
pd_2011 <- rast(tagged_path("./processeddatasets", file.path(predictors_filepath, "30m"), "pop_density2011.tif", "study-area"))
# canfax2024_r <- rast(tagged_path("./processeddatasets", file.path(predictors_filepath, "30m"), "canfax_nt_b_25_2020-2024.tif", "study-area"))
# pc_d_all_2021 <- rast(tagged_path("./processeddatasets", file.path(predictors_filepath, "30m"), "dist_to_popctr2021.tif", "study-area"))
# pc_d_core_2021 <- rast(tagged_path("./processeddatasets", file.path(predictors_filepath, "30m"), "dist_to_regctr2021.tif", "study-area"))
# pd_2021 <- rast(tagged_path("./processeddatasets", file.path(predictors_filepath, "30m"), "pop_density2021.tif", "study-area"))
# elev_large_d <- rast(tagged_path("./processeddatasets", elevator_filepath, "dist_to_elevator_large.tif", "study-area"))
# elev_cap <- rast(tagged_path("./processeddatasets", elevator_filepath, "cap_nearest_elevator.tif", "study-area"))

# Correlation between main s-e predictors
socio <- c(canfax2009_r, elev_d,
           dist_to_road, dist_to_rail,
           pc_d_all_2011, pc_d_core_2011, pd_2011)
names(socio) <- c("canfax2009", "d_elev", 
                  "d_road", "d_rail", "popctr2011",
                  "regctr2011", "popdens2011")

s <- spatSample(socio, 1e4, method = "regular", na.rm = TRUE, as.df = TRUE)
m <- round(cor(s, method = "spearman"), 3)
m[lower.tri(m)] <- NA
m
#             canfax2009 d_elev d_road d_rail popctr2011 regctr2011 popdens2011
# canfax2009           1 -0.056 -0.061  0.041     -0.126     -0.248       0.104
# d_elev              NA  1.000  0.140  0.334      0.526      0.534      -0.488
# d_road              NA     NA  1.000  0.137      0.156      0.130      -0.193
# d_rail              NA     NA     NA  1.000      0.418      0.239      -0.314
# popctr2011          NA     NA     NA     NA      1.000      0.516      -0.501
# regctr2011          NA     NA     NA     NA         NA      1.000      -0.391
# popdens2011         NA     NA     NA     NA         NA         NA       1.000

png(file.path(fig_dir, "socio_correlation.png"), width = 1800, height = 1800, res = 300)
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


# Correlation between elevator predictors
ev <- c(elev_d, elev_large_d, elev_cap)
names(ev) <- c("d_elev", "d_elev_large", "cap_elev")
s <- spatSample(ev, 1e4, method = "regular", na.rm = TRUE, as.df = TRUE)
round(cor(s, method = "spearman"), 3)
#               d_elev d_elev_large cap_elev
# d_elev        1.000        0.355    0.270
# d_elev_large  0.355        1.000   -0.284
# cap_elev      0.270       -0.284    1.000

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
