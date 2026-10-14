# =============================================================================
# Purpose: download railways
# =============================================================================
source(file.path("scripts/landcover", "0a-setup.R"))

# --- Parameters -------------------------------------------------------

rail_dir <- file.path("existingdatasets", capitalize_first(rail_filepath))

# Study area & Modelling extent
sa <- st_read(file.path(sa_filepath, sa_filename))
me <- st_read(file.path(sa_filepath, me_filename))
me_r <- rast(file.path(sa_filepath, me_raster_filename))

url  <- "https://ftp.maps.canada.ca/pub/nrcan_rncan/vector/geobase_nrwn_rfn/"
prov <- c("ab", "sk", "mb")
zips <- sprintf("%s%s/nrwn_rfn_%s_shp_en.zip", url, prov, prov)

for (z in zips) {
  f <- file.path(rail_dir, basename(z))
  download.file(z, f, mode = "wb")
  unzip(f, exdir = rail_dir)
}

shp <- list.files(rail_dir, pattern = "TRACK.*\\.shp$",
                  recursive = TRUE, full.names = TRUE, ignore.case = TRUE)

rail <- do.call(rbind, lapply(shp, vect))   # same product spec, so schemas should align
rail <- project(rail, crs(me))         # NRWN is geographic (NAD83-family) — verify with crs(vect(shp[1]))
rail <- crop(rail, ext(me_r))

rail_r <- rasterize(
  rail, me_r, 
  field = 1, 
  background = NA, 
  touches = TRUE
)
rail_r <- mask(rail_r, me_r)
write_rast_tagged(rail_r, "./processeddatasets", capitalize_first(rail_filepath), paste0(rail_filename, ".tif"), "modelling-extent")


