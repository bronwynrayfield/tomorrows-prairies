# =============================================================================
# Purpose:
#   Prepares protected area mask layers needed to derive spatial
#   predictors in subsequent scripts.
#   All outputs are at native 30m resolution and saved to processeddatasets/.
#
# Inputs:
#   - CPCAD 2024 protected areas : projected and cropped
#
# Outputs:
#   - Protected areas rasters:
#       pa_hard        : IUCN I-V established protected areas (binary mask)
#       pa_provpasture : Saskatchewan provincial pastures (binary)
#       pa_easements   : Conservation easements (saved for future use)
#
# Notes:
#   - pa_hard filter is interim (IUCN I-V, STATUS 1-2) pending expert review
#     of CPCAD classification conventions
#   - pa_easements coverage in CPCAD 2024 is incomplete (Alberta only);
#     saved for future use if DUC/NCC data becomes available
#
# Dependencies:
#   - 0-setup.R (parameters, file paths, functions)
#   - 3-crop-data.R outputs
# =============================================================================

source(file.path("scripts/landcover", "0-setup.R"))

# Output directory
out_dir <- "./processeddatasets/management/Canadian Protected and Conserved Areas Database"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

#--- Load files ----
# Study area
sa <- st_read(file.path(sa_filepath, sa_filename))
me <- st_read(file.path(sa_filepath, me_filename))
me_r <- rast(file.path(sa_filepath, me_raster_filename))

# Protected areas
pa <- st_read(tagged_path("./processeddatasets", pa_filepath, pa_filename, "study-area", driver = "GPKG"),
  layer = paste0(tools::file_path_sans_ext(pa_filename), "_study-area")
)

# Grassland
native <- rast(file.path("./processeddatasets", aafc_filepath, "binary", "AAFC2020_native.tif"))
tame <- rast(file.path("./processeddatasets", aafc_filepath, "binary", "AAFC2020_tame.tif"))

#--- Clean geometry ----
# One CPCAD feature (IUCN_CAT 9, Wildlife Habitat Protection) is a
# GEOMETRYCOLLECTION; extracting to polygons splits it into ~9,023 rows,
# so O_AREA_HA is no longer usable and area is recomputed from geometry
pa <- st_collection_extract(pa, "POLYGON")
pa <- st_cast(pa, "MULTIPOLYGON")
pa <- pa |>
  ungroup() |>
  mutate(area_ha = as.numeric(st_area(pa)) / 10000,
         ZONE_ID = row_number())

# validate before proceeding
stopifnot(length(unique(pa$ZONE_ID)) == nrow(pa))

pa_v <- vect(pa)

#--- Grassland zonal stats ----
# terra::extract(fun = sum) runs out of memory on this polygon set
# (one feature's 9,023 fragments blow up the per-polygon index), so
# rasterize zones and use zonal() instead
zone_r <- terra::rasterize(pa_v, native, field = "ZONE_ID")

grass_stack <- c(native, tame)
names(grass_stack) <- c("native", "tame")
cell_ha <- prod(res(native)) / 10000

zonal_sums <- terra::zonal(grass_stack, zone_r, fun = "sum", na.rm = TRUE) |>
  rename(ZONE_ID = 1)

# Join per-polygon grassland cell counts back onto pa, convert cell counts
# to hectares, and derive combined grassland total
pa <- pa |>
  left_join(zonal_sums, by = "ZONE_ID") |>
  mutate(
    native_ha    = coalesce(native, 0) * cell_ha,
    tame_ha      = coalesce(tame, 0) * cell_ha,
    grassland_ha = native_ha + tame_ha
  ) |>
  select(-ZONE_ID, -native, -tame)

# Check
# 1. Any polygons that didn't get a zone match?
sum(is.na(pa$grassland_ha))   # expect 0, or investigate why if not
# [1] 0

# 2. Sanity check: grassland can't exceed total polygon area
sum(pa$grassland_ha > pa$area_ha, na.rm = TRUE)   # expect 0
# [1] 171

# 3. Quick look at the distribution
summary(pa$grassland_ha)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.00      0.00      8.82    159.42     63.72 210301.20 
summary(pa$area_ha)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.00     52.04     64.83    600.88    191.78 507308.86 

# 4. Spot check problematic feature (the split GEOMETRYCOLLECTION)
# total_area_ha should be ~1,367,146 ha, matching the original pre-split O_AREA_HA
pa |>
  st_drop_geometry() |>
  filter(TYPE_E == "Wildlife Habitat Protection", IUCN_CAT == 9) |>
  summarise(n = n(), total_area_ha = sum(area_ha), total_grassland_ha = sum(grassland_ha))
# n total_area_ha total_grassland_ha
# 1 9023       1367146           741620.5

# Investigate rows that failed check #2. rasterize() assigns whole cells by
# cell-center containment, not proportional overlap, so a polygon can pick up
# slightly more grassland area than its own true geometric area near its edges
mismatch <- pa |>
  st_drop_geometry() |>
  filter(grassland_ha > area_ha)

summary(mismatch$area_ha)          # expect these to be small/sub-pixel
sum(mismatch$area_ha)               # total area involved
sum(mismatch$grassland_ha - mismatch$area_ha)   # total excess hectares from this issue

# Cap grassland at each polygon's own area to remove the sub-pixel edge
# effect confirmed above; only affects the 171 rows that failed check #2
pa <- pa |>
  mutate(
    native_ha    = pmin(native_ha, area_ha),
    tame_ha      = pmin(tame_ha, area_ha),
    grassland_ha = pmin(grassland_ha, area_ha)
  )


#--- Summary tables ----

# IUCN category summary
iucn_summary <- pa |>
  st_drop_geometry() |>
  group_by(IUCN_CAT) |>
  summarise(
    n_features    = n(),
    total_ha      = round(sum(area_ha, na.rm = TRUE)),
    native_ha     = round(sum(native_ha, na.rm = TRUE)),
    tame_ha       = round(sum(tame_ha, na.rm = TRUE)),
    grassland_ha  = round(sum(grassland_ha, na.rm = TRUE)),
    pct_grassland = round(100 * grassland_ha / total_ha, 1)
  ) |>
  arrange(IUCN_CAT) |>
  as.data.frame()

# Category IV TYPE_E breakdown
catIV_summary <- pa |>
  st_drop_geometry() |>
  filter(IUCN_CAT == 5) |>
  group_by(TYPE_E) |>
  summarise(
    n             = n(),
    area_ha       = round(sum(area_ha, na.rm = TRUE)),
    native_ha     = round(sum(native_ha, na.rm = TRUE)),
    tame_ha       = round(sum(tame_ha, na.rm = TRUE)),
    grassland_ha  = round(sum(grassland_ha, na.rm = TRUE)),
    pct_grassland = round(100 * grassland_ha / area_ha, 1)
  ) |>
  arrange(desc(grassland_ha)) |>
  as.data.frame()

# Category V TYPE_E breakdown
catV_summary <- pa |>
  st_drop_geometry() |>
  filter(IUCN_CAT == 6) |>
  group_by(TYPE_E) |>
  summarise(
    n             = n(),
    area_ha       = round(sum(area_ha, na.rm = TRUE)),
    native_ha     = round(sum(native_ha, na.rm = TRUE)),
    tame_ha       = round(sum(tame_ha, na.rm = TRUE)),
    grassland_ha  = round(sum(grassland_ha, na.rm = TRUE)),
    pct_grassland = round(100 * grassland_ha / area_ha, 1)
  ) |>
  arrange(desc(grassland_ha)) |>
  as.data.frame()

# Category VI TYPE_E breakdown
catVI_summary <- pa |>
  st_drop_geometry() |>
  filter(IUCN_CAT == 7) |>
  group_by(TYPE_E) |>
  summarise(
    n             = n(),
    area_ha       = round(sum(area_ha, na.rm = TRUE)),
    native_ha     = round(sum(native_ha, na.rm = TRUE)),
    tame_ha       = round(sum(tame_ha, na.rm = TRUE)),
    grassland_ha  = round(sum(grassland_ha, na.rm = TRUE)),
    pct_grassland = round(100 * grassland_ha / area_ha, 1)
  ) |>
  arrange(desc(grassland_ha)) |>
  as.data.frame()

# Category Not Applicable/Not Reported TYPE_E breakdown
catNANR_summary <- pa |>
  st_drop_geometry() |>
  filter(IUCN_CAT %in% c(8, 9)) |>
  group_by(TYPE_E) |>
  summarise(
    n             = n(),
    area_ha       = round(sum(area_ha, na.rm = TRUE)),
    native_ha     = round(sum(native_ha, na.rm = TRUE)),
    tame_ha       = round(sum(tame_ha, na.rm = TRUE)),
    grassland_ha  = round(sum(grassland_ha, na.rm = TRUE)),
    pct_grassland = round(100 * grassland_ha / area_ha, 1)
  ) |>
  arrange(desc(grassland_ha)) |>
  as.data.frame()

# Write out summaries
write.csv(iucn_summary,   file.path("./processeddatasets", pa_filepath, "cpcad_iucn_summary_study-area.csv"),        row.names = FALSE)
write.csv(catIV_summary,   file.path("./processeddatasets", pa_filepath, "cpcad_catIV_type_e_summary_study-area.csv"), row.names = FALSE)
write.csv(catV_summary,   file.path("./processeddatasets", pa_filepath, "cpcad_catV_type_e_summary_study-area.csv"), row.names = FALSE)
write.csv(catVI_summary,   file.path("./processeddatasets", pa_filepath, "cpcad_catVI_type_e_summary_study-area.csv"), row.names = FALSE)
write.csv(catNANR_summary,  file.path("./processeddatasets", pa_filepath, "cpcad_catNANR_type_e_summary_study-area.csv"),row.names = FALSE)

# Write cleaned PA layer
st_write(
  pa,
  tagged_path("./processeddatasets", pa_filepath, pa_filename, "study-area_cleaned-grassland", driver = "GPKG"),
  layer = paste0(tools::file_path_sans_ext(pa_filename), "_study-area_cleaned-grassland"),
  delete_layer = TRUE
)





















#--- Grassland zonal stats ----
grass_stack <- c(native, tame)
names(grass_stack) <- c("native", "tame")

cell_ha <- prod(res(native)) / 10000

#--- Clean geometry ----
# One CPCAD feature (IUCN_CAT 9, Wildlife Habitat Protection) is a
# GEOMETRYCOLLECTION; extracting to polygons splits it into ~9,023 rows,
# so O_AREA_HA is no longer usable and area is recomputed from geometry
pa <- st_collection_extract(pa, "POLYGON")
pa <- st_cast(pa, "MULTIPOLYGON")
pa <- pa |>
  mutate(area_ha = as.numeric(st_area(pa)) / 10000)

pa_v <- vect(pa)

#--- Grassland zonal stats ----
# terra::extract(fun = sum) runs out of memory on this polygon set
# (one feature's 9,023 fragments blow up the per-polygon index), so
# rasterize zones and use zonal() instead
zone_r <- terra::rasterize(pa_v, native, field = "ZONE_ID")

grass_stack <- c(native, tame)
names(grass_stack) <- c("native", "tame")

cell_ha <- prod(res(native)) / 10000

zonal_sums <- terra::zonal(grass_stack, zone_r, fun = "sum", na.rm = TRUE) |>
  rename(zone_id = 1)

pa <- pa |>
  left_join(zonal_sums, by = c("ZONE_ID" = "zone_id")) |>
  mutate(
    native_ha    = coalesce(native, 0) * cell_ha,
    tame_ha      = coalesce(tame, 0) * cell_ha,
    grassland_ha = native_ha + tame_ha
  ) |>
  select(-ZONE_ID, -native, -tame)

# Check
# 1. Any polygons that didn't get a zone match?
sum(is.na(pa$grassland_ha))   # expect 0, or investigate why if not

# 2. Sanity check: grassland can't exceed total polygon area
sum(pa$grassland_ha > pa$area_ha, na.rm = TRUE)   # expect 0

# 3. Quick look at the distribution
summary(pa$grassland_ha)
summary(pa$area_ha)

# 4. Spot check the problem feature from earlier (the split GEOMETRYCOLLECTION)
pa |>
  st_drop_geometry() |>
  filter(TYPE_E == "Wildlife Habitat Protection", IUCN_CAT == 9) |>
  summarise(n = n(), total_area_ha = sum(area_ha), total_grassland_ha = sum(grassland_ha))













pa_v <- vect(pa)


pa <- st_collection_extract(pa, "POLYGON")
pa <- st_cast(pa, "MULTIPOLYGON")

pa_v <- vect(pa)

pa |>
  st_drop_geometry() |>
  filter(row_number() == which(st_geometry_type(pa) == "GEOMETRYCOLLECTION")) |>
  select(any_of(c("IUCN_CAT", "TYPE_E", "O_AREA_HA")))


grass_counts <- terra::extract(grass_stack, pa_v, fun = sum, na.rm = TRUE)

pa <- pa |>
  mutate(
    native_ha    = grass_counts$native * cell_ha,
    tame_ha      = grass_counts$tame * cell_ha,
    grassland_ha = native_ha + tame_ha
  )

#--- Summary tables ----
# IUCN category summary
pa |>
  st_drop_geometry() |>
  group_by(IUCN_CAT) |>
  summarise(
    n_features    = n(),
    total_ha      = round(sum(O_AREA_HA, na.rm = TRUE)),
    native_ha     = round(sum(native_ha, na.rm = TRUE)),
    tame_ha       = round(sum(tame_ha, na.rm = TRUE)),
    grassland_ha  = round(sum(grassland_ha, na.rm = TRUE)),
    pct_grassland = round(100 * grassland_ha / total_ha, 1)
  ) |>
  arrange(IUCN_CAT) |>
  as.data.frame() |>
  print()

# Category IV (5) TYPE_E breakdown
pa |>
  st_drop_geometry() |>
  filter(IUCN_CAT == 5) |>
  group_by(TYPE_E) |>
  summarise(
    n             = n(),
    area_ha       = round(sum(O_AREA_HA, na.rm = TRUE)),
    native_ha     = round(sum(native_ha, na.rm = TRUE)),
    tame_ha       = round(sum(tame_ha, na.rm = TRUE)),
    grassland_ha  = round(sum(grassland_ha, na.rm = TRUE)),
    pct_grassland = round(100 * grassland_ha / area_ha, 1)
  ) |>
  arrange(desc(n)) |>
  as.data.frame() |>
  print()

# Category 8/9 TYPE_E breakdown
pa |>
  st_drop_geometry() |>
  filter(IUCN_CAT %in% c(8, 9)) |>
  group_by(TYPE_E) |>
  summarise(
    n             = n(),
    area_ha       = round(sum(O_AREA_HA, na.rm = TRUE)),
    native_ha     = round(sum(native_ha, na.rm = TRUE)),
    tame_ha       = round(sum(tame_ha, na.rm = TRUE)),
    grassland_ha  = round(sum(grassland_ha, na.rm = TRUE)),
    pct_grassland = round(100 * grassland_ha / area_ha, 1)
  ) |>
  arrange(desc(n)) |>
  as.data.frame() |>
  print()
























# easement_types <- c(
#   "Conservation Easement",
#   "Conservation easement",
#   "Conservation Easement, Ownership by Environmental Non-Governmental Organization",
#   "Conservation Agreement",
#   "Conservation covenant",
#   "Natural Area Under Private Stewardship",
#   "NGO/Land Trust/Agency - Fee Simple Property",
#   "Private Ecologically Significant Area"
# )
#
# ngo_ownership_types <- c(
#   "Ownership by Environmental Non-Governmental Organization"
# )
#
# additional_hard_types <- c(
#   "Fee Simple Property",
#   "Fee simple – Ownership by land trust or conservancy",
#   "National Wildlife Area",
#   "Wildlife Habitat Areas",
#   "Wildlife Habitat Protection",
#   "Conservation Area",
#   "Privately Owned Conservation Area",
#   "Privately Owned Conservation Area (INT)"
# )
#
# # Hard protection mask
# pa_hard <- pa |>
#   filter(
#     STATUS %in% c(1, 2) &
#       !TYPE_E %in% easement_types &
#       (IUCN_CAT %in% c(1, 2, 3, 4, 5) |
#          (IUCN_CAT %in% c(8, 9) & TYPE_E %in% c(
#            "A - Park", "C - Park",
#            "Ecological Reserve",
#            "Nature Reserve",
#            "Migratory Bird Sanctuary",
#            "Wildlife Management Area",
#            "Wildland Provincial Park",
#            "Provincial Park",
#            "Wilderness Area",
#            "Conservancy"
#          )) |
#          TYPE_E %in% ngo_ownership_types |
#          TYPE_E %in% additional_hard_types)
#   )
#
# # Easements
# pa_easements <- pa |>
#   filter(
#     STATUS %in% c(1, 2) &
#       TYPE_E %in% easement_types
#   )
#
# # Provincial pastures
# pa_provpasture <- pa |>
#   filter(
#     STATUS %in% c(1, 2) &
#       TYPE_E == "Provincial Pasture"
#   )

# Hard protection mask — IUCN I-V, established areas only
# -----------------------------------------------------------------------------
# TODO: Revisit filtering before final model run:
#   1. IUCN category filter — currently using I-V (codes 1-5); strict Barnes
#      et al. equivalent would be I-III only (codes 1-3) but this gives only
#      194 features, suggesting many Prairie protected areas are miscoded in
#      CPCAD. Pending expert review of CPCAD classification conventions.
#   2. Establishment year — consider filtering to ESTYEAR <= 2010 to match
#      predictor year. However 2,798 of 4,663 features have no ESTYEAR,
#      including NGO ownership (n=2,238) and Provincial Pastures (n=112)
#      which almost certainly predate 2010. Low practical impact but worth
#      noting as a limitation. 739 features confirmed post-2010.
#   3. TYPE_E inclusions — earlier analysis identified TYPE_E types coded as
#      IUCN 8/9 ("Not Applicable"/"Not Reported") that represent genuine hard
#      protection (Provincial Parks, Ecological Reserves, Migratory Bird
#      Sanctuaries, National Wildlife Areas etc.) and are currently excluded
#      by the IUCN I-V filter. Consider adding explicit TYPE_E inclusions for
#      these well-defined protection types regardless of IUCN coding.
#      See earlier exploration code for full list of TYPE_E values in codes 8/9.
# -----------------------------------------------------------------------------
# Explore pa data
# Crop pa to study area
pa_sa <- st_intersection(pa, sa)

# IUCN summary
pa_sa |>
  st_drop_geometry() |>
  group_by(IUCN_CAT) |>
  summarise(
    n_features = n(),
    total_ha   = round(sum(O_AREA_HA, na.rm = TRUE))
  ) |>
  arrange(IUCN_CAT) |>
  mutate(total_ha = format(total_ha, big.mark = ",")) |>
  as.data.frame() |>
  print()

# Category IV breakdown
pa_sa |>
  st_drop_geometry() |>
  filter(IUCN_CAT == 5) |>
  group_by(TYPE_E) |>
  summarise(n = n(), area_ha = round(sum(O_AREA_HA, na.rm = TRUE))) |>
  arrange(desc(n)) |>
  as.data.frame() |>
  print()

# Categories 8 and 9 breakdown
pa_sa |>
  st_drop_geometry() |>
  filter(IUCN_CAT %in% c(8, 9)) |>
  group_by(TYPE_E) |>
  summarise(n = n(), area_ha = round(sum(O_AREA_HA, na.rm = TRUE))) |>
  arrange(desc(n)) |>
  as.data.frame() |>
  print()

write_sf_tagged(pa_sa, "./processeddatasets", pa_filepath, pa_filename,
                "study-area", driver = "GPKG")





pa_hard <- pa |>
  filter(
    STATUS %in% c(1, 2) &
      IUCN_CAT %in% c(1, 2, 3, 4, 5)
  )

# Provincial pastures
pa_provpasture <- pa |>
  filter(
    STATUS %in% c(1, 2) &
      TYPE_E == "Provincial Pasture"
  )

# Easements — saved for future use if DUC/NCC data becomes available
pa_easements <- pa |>
  filter(
    STATUS %in% c(1, 2) &
      TYPE_E %in% c(
        "Conservation Easement",
        "Conservation easement",
        "Conservation Easement, Ownership by Environmental Non-Governmental Organization",
        "Conservation Agreement",
        "Conservation covenant",
        "Natural Area Under Private Stewardship",
        "NGO/Land Trust/Agency - Fee Simple Property",
        "Private Ecologically Significant Area"
      )
  )

# Check
nrow(pa_hard)
nrow(pa_provpasture)
nrow(pa_easements)

# How many are unaccounted for?
# ok, these are proposed areas or recreation sites, municipal parks etc.
nrow(pa) - nrow(pa_hard) - nrow(pa_easements) - nrow(pa_provpasture)

# Overlap between layers
# ok, there is some overlap but they'll just get written in binary rasters
any(st_intersects(pa_hard, pa_easements, sparse = FALSE))
any(st_intersects(pa_hard, pa_provpasture, sparse = FALSE))
any(st_intersects(pa_easements, pa_provpasture, sparse = FALSE))

# Check total area looks reasonable
sum(pa_hard$O_AREA_HA, na.rm = TRUE)
sum(pa_easements$O_AREA_HA, na.rm = TRUE)
sum(pa_provpasture$O_AREA_HA, na.rm = TRUE)


# Rasterize
pa_hard_rast <- rasterize(vect(pa_hard), me_r,
                          field = 1,
                          background = 0,
                          touches = TRUE)

# Mask to modelling extent
pa_hard_rast <- mask(pa_hard_rast, me_r)

# Rasterize
pa_provpasture_rast <- rasterize(vect(pa_provpasture), me_r,
                                 field = 1,
                                 background = 0,
                                 touches = TRUE)
# Mask to modelling extent
pa_provpasture_rast <- mask(pa_provpasture_rast, me_r)


#--- Save to disk ----
# LULC binary layers
# Binary AAFC filenames = cleaner than using the merged filename
aafc2015_binary_filename <- "AAFC2015.tif"
aafc2010_binary_filename <- "AAFC2010.tif"

# Binary filepath
aafc_binary_filepath <- file.path(aafc_filepath, "binary")

# 2015
write_rast_tagged(native2015_bin,     "./processeddatasets", aafc_binary_filepath, aafc2015_binary_filename, "native",     datatype = "INT1U")
write_rast_tagged(tame2015_bin,       "./processeddatasets", aafc_binary_filepath, aafc2015_binary_filename, "tame",       datatype = "INT1U")
write_rast_tagged(crop2015_bin,       "./processeddatasets", aafc_binary_filepath, aafc2015_binary_filename, "crop",       datatype = "INT1U")
write_rast_tagged(settlement2015_bin, "./processeddatasets", aafc_binary_filepath, aafc2015_binary_filename, "settlement", datatype = "INT1U")
write_rast_tagged(forest2015_bin,     "./processeddatasets", aafc_binary_filepath, aafc2015_binary_filename, "forest",     datatype = "INT1U")
write_rast_tagged(water2015_bin,      "./processeddatasets", aafc_binary_filepath, aafc2015_binary_filename, "water",      datatype = "INT1U")
# 2010
write_rast_tagged(native2010_bin,     "./processeddatasets", aafc_binary_filepath, aafc2010_binary_filename, "native",     datatype = "INT1U")
write_rast_tagged(tame2010_bin,       "./processeddatasets", aafc_binary_filepath, aafc2010_binary_filename, "tame",       datatype = "INT1U")
write_rast_tagged(settlement2010_bin, "./processeddatasets", aafc_binary_filepath, aafc2010_binary_filename, "settlement", datatype = "INT1U")
write_rast_tagged(forest2010_bin,     "./processeddatasets", aafc_binary_filepath, aafc2010_binary_filename, "forest",     datatype = "INT1U")
write_rast_tagged(water2010_bin,      "./processeddatasets", aafc_binary_filepath, aafc2010_binary_filename, "water",      datatype = "INT1U")


# Protected areas
# Save intermediate PA shapefiles
write_sf_tagged(pa_hard, "./processeddatasets", pa_filepath, pa_filename, "hard", driver = "GPKG")
write_sf_tagged(pa_provpasture, "./processeddatasets", pa_filepath, pa_filename, "provpasture", driver = "GPKG")
write_sf_tagged(pa_easements, "./processeddatasets", pa_filepath, pa_filename, "easements", driver = "GPKG")

# Save PA rasters
write_rast_tagged(pa_hard_rast, "./processeddatasets", pa_filepath, pa_filename, "hard", "binary", datatype = "INT1U")
write_rast_tagged(pa_provpasture_rast, "./processeddatasets", pa_filepath, pa_filename, "provpasture", "binary", datatype = "INT1U")


