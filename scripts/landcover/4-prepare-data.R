# =============================================================================
# 4-prepare-data.R
# =============================================================================
#
# Purpose:
#   Prepares intermediate binary and mask layers needed to derive spatial
#   predictors in subsequent scripts.
#   All outputs are at native 30m resolution and saved to processeddatasets/.
#   Resolution-dependent operations (focal statistics etc.) are in 5-*.R
#
# Inputs:
#   - AAFC LUTS 30m rasters (2010, 2015) : projected and cropped
#   - CPCAD 2024 protected areas : projected and cropped
#
# Outputs:
#   - AAFC binary rasters (native, tame, crop, settlement, forest, water)
#     for 2010 and 2015
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
#   - crop class (codes 52, 56) not present in 2010 AAFC LUTS;
#
# Dependencies:
#   - 0-setup.R (parameters, file paths, functions)
#   - 3-crop-data.R outputs
# =============================================================================

source(file.path("scripts/landcover", "0-setup.R"))
source(file.path("scripts/landcover/functions", "0-functions-general.R"))
source(file.path("scripts/landcover/functions", "reclass-aafc-luts.R"))

#--- Load files ----
# Study area
sa <- st_read(file.path(sa_filepath, sa_filename))
me <- st_read(file.path(sa_filepath, me_filename))
me_r <- rast(file.path(sa_filepath, me_raster_filename))

# Land cover (AAFC LUTS)
lulc2010 <- rast(tagged_path(file.path("./processeddatasets", aafc_filepath), "cropped", aafc2010_filename, "modelling-extent"))
lulc2015 <- rast(tagged_path(file.path("./processeddatasets", aafc_filepath), "cropped", aafc2015_filename, "modelling-extent"))
lulc2020 <- rast(tagged_path(file.path("./processeddatasets", aafc_filepath), "cropped", aafc2020_filename, "modelling-extent"))

# Protected areas
pa <- read_sf_tagged("./processeddatasets", pa_filepath, pa_filename, "modelling-extent", driver = "GPKG")

# Olimb & Robinson 2019 cropland suitability
crop_suit <- rast(tagged_path("./processeddatasets", crop_suit_2019_filepath, crop_suit_2019_filename, "modelling-extent"))



#--- LULC binary layers ----
# --- 2015 binary layers ---
native2015_bin     <- reclass_aafc(lulc2015, scheme = "binary", binary_class = "native",     check = TRUE)
tame2015_bin       <- reclass_aafc(lulc2015, scheme = "binary", binary_class = "tame",       check = TRUE)
crop2015_bin       <- reclass_aafc(lulc2015, scheme = "binary", binary_class = "crop",       check = TRUE)
settlement2015_bin <- reclass_aafc(lulc2015, scheme = "binary", binary_class = "settlement", check = TRUE)
forest2015_bin     <- reclass_aafc(lulc2015, scheme = "binary", binary_class = "forest",     check = TRUE)
water2015_bin      <- reclass_aafc(lulc2015, scheme = "binary", binary_class = "water",      check = TRUE)

# --- 2010 binary layers ---
native2010_bin     <- reclass_aafc(lulc2010, scheme = "binary", binary_class = "native",     check = TRUE)
# Note that AAFC LUTS 2010 doesn't separate tame and annual cropland, both are lumped into "tame"
tame2010_bin       <- reclass_aafc(lulc2010, scheme = "binary", binary_class = "tame",       check = TRUE)
settlement2010_bin <- reclass_aafc(lulc2010, scheme = "binary", binary_class = "settlement", check = TRUE)
forest2010_bin     <- reclass_aafc(lulc2010, scheme = "binary", binary_class = "forest",     check = TRUE)
water2010_bin      <- reclass_aafc(lulc2010, scheme = "binary", binary_class = "water",      check = TRUE)


#--- Protected area rasters ----

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
# Note: interim filter pending expert review of CPCAD classification conventions
# Other filtering options in the comments above
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


# Rasterize and save pa_hard
pa_hard_rast <- rasterize(vect(pa_hard), me_r,
                          field = 1,
                          background = 0,
                          touches = TRUE)

# Rasterize and save pa_provpasture
pa_provpasture_rast <- rasterize(vect(pa_provpasture), me_r,
                                 field = 1,
                                 background = 0,
                                 touches = TRUE)




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


