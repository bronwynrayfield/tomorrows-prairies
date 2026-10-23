# =============================================================================
# Purpose:
# Create protected areas layer
#
# Protected areas: CPCAD 2024 + SK Crown Conservation Easements
#
# Output: 2010 and 2024 rasters on me_r. 0 = unprotected, 1 = CPCAD,
# 2 = SK CCE not in CPCAD, NA outside me_r.
#
# Most CCEs are already in CPCAD.
# Both datasets are filtered by year. No filtering by IUCN category.
#
# Checked against NCC Tenure data layer
#
# Note: former PFRA pastures are already in CPCAD
# =============================================================================

source(file.path("scripts/landcover", "0a-setup.R"))

# Output directory
out_dir <- "./processeddatasets/management/Protected Areas Combined"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

#--- Load files ----
# Study area
sa <- st_read(file.path(sa_filepath, sa_filename))
sa_r <- rast(file.path(sa_filepath, sa_raster_filename))
me <- st_read(file.path(sa_filepath, me_filename))
me_r <- rast(file.path(sa_filepath, me_raster_filename))

# Protected areas (CPCAD)
pa <- st_read(tagged_path("./processeddatasets", pa_filepath, pa_filename, "modelling-extent", driver = "GPKG"),
              layer = paste0(tools::file_path_sans_ext(pa_filename), "_modelling-extent")
)

# NCC/NatureUnited ownership
ncc <- st_read(
  tagged_path("./processeddatasets", ncc_easement_filepath, ncc_easement_filename, "modelling-extent", driver = "GPKG"),
  layer = paste0(tools::file_path_sans_ext(ncc_easement_filename), "_modelling-extent")
)

# SK Crown Conservation Easements
sk_cce <- st_read(
  tagged_path("./processeddatasets", sk_cce_filepath, sk_cce_filename, "modelling-extent", driver = "GPKG"),
  layer = paste0(tools::file_path_sans_ext(sk_cce_filename), "_modelling-extent")
)

# Common geometry column name for bind_rows
st_geometry(pa) <- "geometry"
st_geometry(ncc) <- "geometry"
st_geometry(sk_cce) <- "geometry"


# --- protected areas from CPCAD -----------------------------------
pa_union <- st_union(pa)

# # --- per-parcel coverage fraction between ncc "PA - Public" lands and CPCAD -------------------------------------------
# Only needed to run this once to understand the overlap
# ncc <- ncc |> filter(OWNERSHIP_FINAL == "PA - Public") |>
#   mutate(pid = row_number(), 
#          area_ha = units::set_units(st_area(geometry), "ha"))
# 
# cov <- st_intersection(ncc, pa_union) |>
#   st_collection_extract("POLYGON") |>          # drop line artifacts at shared edges
#   mutate(inside_ha = units::set_units(st_area(geometry), "ha")) |>
#   st_drop_geometry() |>
#   summarise(inside_ha = sum(inside_ha), .by = pid)   # parcels split across >1 pa
# 
# ncc <- ncc |> left_join(cov, by = "pid") |>
#   mutate(inside_ha = units::set_units(replace_na(as.numeric(inside_ha), 0), "ha"),
#          frac_in = as.numeric(inside_ha / area_ha))
# 
# hist(ncc$frac_in, breaks = 50)
# 
# # --- outside-pa parcels, rasterized -----------------------------------------
# # 0.05 threshold: frac_in distribution is bimodal with 22,598 parcels >0.99 and
# # 380 <0.01; only 26 parcels fall between 0.01 and 0.95
# ncc_out <- ncc |> filter(frac_in < 0.05)
# 
# # Remove digitization slivers: parcels that vanish under a 3-cell (90 m) negative
# # buffer are narrower than ~180 m and are boundary remnants from ncc/pa geometry
# # mismatch, not real land units 
# buf <- -(terra::res(me_r)[1])*3
# ncc_out <- ncc_out |> filter(!st_is_empty(st_buffer(geometry, buf)))
# nrow(ncc_out)
# # [1] 37
# ncc_out |> st_drop_geometry()
# # All ncc features that did no overlap with CPCAD were in SK
# # 29 features had SOURCE = Crown Conservation Easements (2023-06-01)
# # 8 feature had SOURCE = Federal Pastures (2023-06-01)
# # 1 feature had SOURCE = Directory of Federal Real Property (2023-05-30)
# # Therefore, instead of using this compiled dataset from NCC, we directly downloaded the SK Crown Conservation Easements and Federal Pastures.
# 
# # rasterize
# ncc_out_r <- terra::rasterize(terra::vect(ncc_out), me_r, field = 1)
# cell_ha <- prod(terra::res(me_r)) / 1e4
# terra::global(ncc_out_r, "sum", na.rm = TRUE)[[1]] * cell_ha
# # [1] 14023.98
# 
# terra::writeRaster(ncc_out_r, file.path("processeddatasets", "management", "Protected Areas Combined", "ncc_pa_public_outside_cpcad.tif"),
#                    overwrite = TRUE, datatype = "INT1U", gdal = "COMPRESS=DEFLATE")


# --- per-parcel coverage fraction between ncc "PA - Public" lands and CPCAD -------------------------------------------
sk_cce <- sk_cce |> mutate(pid = row_number(), 
         area_ha = units::set_units(st_area(geometry), "ha"))

sk_cce_cov <- st_intersection(sk_cce, pa_union) |>
  st_collection_extract("POLYGON") |>          # drop line artifacts at shared edges
  mutate(inside_ha = units::set_units(st_area(geometry), "ha")) |>
  st_drop_geometry() |>
  summarise(inside_ha = sum(inside_ha), .by = pid)   # parcels split across >1 pa

sk_cce <- sk_cce |> left_join(sk_cce_cov, by = "pid") |>
  mutate(inside_ha = units::set_units(replace_na(as.numeric(inside_ha), 0), "ha"),
         frac_in = as.numeric(inside_ha / area_ha))

hist(sk_cce$frac_in, breaks = 50)
table(cut(sk_cce$frac_in, c(-Inf, 0.01, 0.5, 0.9, 0.95, 0.99, Inf)))
# (-Inf,0.01]  (0.01,0.5]   (0.5,0.9]  (0.9,0.95] (0.95,0.99] (0.99, Inf] 
# 138           1           3           2           3         924

# --- outside-pa parcels, rasterize -----------------------------------------
# 0.05 threshold: frac_in distribution is bimodal with 22,598 parcels >0.99 and
# 380 <0.01; only 26 parcels fall between 0.01 and 0.95
sk_cce_out <- sk_cce |> filter(frac_in < 0.05)

# Remove digitization slivers: parcels that vanish under a 3-cell (90 m) negative
# buffer are narrower than ~180 m and are boundary remnants from ncc/pa geometry
# mismatch, not real land units 
buf <- -(terra::res(me_r)[1])*3
sk_cce_out <- sk_cce_out |> filter(!st_is_empty(st_buffer(geometry, buf)))
nrow(sk_cce_out)
# [1] 106
sk_cce_out |> st_drop_geometry()

sk_cce_out_r <- terra::rasterize(terra::vect(sk_cce_out), me_r, field = 1)
cell_ha <- prod(terra::res(me_r)) / 1e4
terra::global(sk_cce_out_r2, "sum", na.rm = TRUE)[[1]] * cell_ha
# [1] 5063.4


# Create a protected areas layer for 2010 and 2025
pa_layer <- function(yr) {
  p <- pa |> filter(is.na(ESTYEAR) | ESTYEAR <= yr)
  c <- sk_cce_out |> filter(PROTDATE <= yr)
  rbind(st_sf(source = rep(1L, nrow(p)), geometry = st_geometry(p)),
        st_sf(source = rep(2L, nrow(c)), geometry = st_geometry(c))) |>
  st_transform(st_crs(ncc)) |> st_make_valid()
}

pa_2010 <- pa_layer(2010)
pa_2025 <- pa_layer(2024) # only keep pa's designated by the end of 2024

pa_2010_r <- terra::rasterize(terra::vect(pa_2010), me_r, field = "source", fun = "min")
pa_2025_r <- terra::rasterize(terra::vect(pa_2025), me_r, field = "source", fun = "min")
# 0 inside the modelling extent, NA outside — me_r defines the valid domain
pa_2010_r <- terra::mask(terra::subst(pa_2010_r, NA, 0L), me_r)
pa_2025_r <- terra::mask(terra::subst(pa_2025_r, NA, 0L), me_r)

terra::freq(pa_2025_r) |> transform(area_ha = count * cell_ha)
# layer value     count     area_ha
# 1     1     0 635056114 57155050.26
# 2     1     1  50744246  4566982.14
# 3     1     2     51023     4592.07

# Write to disk
terra::writeRaster(sk_cce_out_r, file.path("processeddatasets", "management", "Protected Areas Combined", "sk_cce_outside_cpcad.tif"),
                   overwrite = TRUE, datatype = "INT1U", gdal = "COMPRESS=DEFLATE")
terra::writeRaster(pa_2010_r, file.path("processeddatasets", "management", "Protected Areas Combined", "CPCAD_pa_2010.tif"),
                   overwrite = TRUE, datatype = "INT1U", gdal = "COMPRESS=DEFLATE")
terra::writeRaster(pa_2025_r, file.path("processeddatasets", "management", "Protected Areas Combined", "CPCAD-SK-CCE_pa_2025_source.tif"),
                   overwrite = TRUE, datatype = "INT1U", gdal = "COMPRESS=DEFLATE")
terra::writeRaster(pa_2025_r>0, file.path("processeddatasets", "management", "Protected Areas Combined", "CPCAD-SK-CCE_pa_2025.tif"),
                   overwrite = TRUE, datatype = "INT1U", gdal = "COMPRESS=DEFLATE")

