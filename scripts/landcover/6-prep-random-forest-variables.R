# Data pre-processing -> cropping

# Setup
source(file.path("scripts/landcover", "0-setup.R"))
source(file.path("scripts/landcover/functions", "0-functions-general.R"))
source(file.path("scripts/landcover/functions", "reclass-aafc-luts.R"))

# ---- Resolution setting ----
target_res <- 30
agg_fact   <- target_res / 30

outdir <- file.path("./output/landcover/conversion-predictors", paste0(target_res, "m"))

#--- Load files ----
# Study area
sa <- st_read(file.path(sa_filepath, sa_filename))
sa_r <- rast(file.path(sa_filepath, sa_raster_filename))
# me <- st_read(file.path(sa_filepath, me_filename))
me_r <- rast(file.path(sa_filepath, me_raster_filename))

ecozones <- st_read(tagged_path("./processeddatasets", ecozones_filepath, ecozones_filename, "modelling-extent"))
ecoregions <- st_read(tagged_path("./processeddatasets", ecoregions_filepath, ecoregions_filename, "modelling-extent"))
ecodistricts <- st_read(tagged_path("./processeddatasets", ecodistricts_filepath, ecodistricts_filename, "modelling-extent"))

lulc2020 <- rast(tagged_path(file.path("./processeddatasets", aafc_filepath), "cropped", aafc2020_filename, "modelling-extent"))
lulc2015 <- rast(tagged_path(file.path("./processeddatasets", aafc_filepath), "cropped", aafc2015_filename, "modelling-extent"))
lulc2010 <- rast(tagged_path(file.path("./processeddatasets", aafc_filepath), "cropped", aafc2010_filename, "modelling-extent"))


# Olimb & Robinson 2019 cropland suitability
crop_suit <- rast(tagged_path("./processeddatasets", crop_suit_2019_filepath, crop_suit_2019_filename, "modelling-extent"))


# ---- Resample to target resolution (skip if target_res == 30) ----
if (agg_fact > 1) {
  sa_r      <- aggregate(sa_r, fact = agg_fact, fun = "modal")
  lulc2010  <- resample(lulc2010, sa_r, method = "near")
  lulc2020  <- resample(lulc2020, sa_r, method = "near")
  crop_suit <- resample(crop_suit, sa_r, method = "bilinear")
}

# --- Response variable ----
# get numeric codes from lookup tables
lut_2010 <- levels(lulc2010)[[1]]
lut_2020 <- levels(lulc2020)[[1]]

native_labels    <- c("Prairie Grassland", "Temperate Grassland",
                      "Temperate Rainforest Grassland", "Boreal Grassland",
                      "Montane Grassland", "Alpine Tundra Grassland")
crop_labels_2020 <- c("Annual Cropland", "Land Converted to Annual Cropland")

native_codes <- lut_2010$Value[lut_2010$LU2010 %in% native_labels]
crop_codes   <- lut_2020$Value[lut_2020$LU2020 %in% crop_labels_2020]

# drop factor encoding
lulc2010_num <- as.numeric(lulc2010)
lulc2020_num <- as.numeric(lulc2020)

# create response variable
was_native <- app(lulc2010_num, fun = function(x) as.integer(x %in% native_codes))
is_crop    <- app(lulc2020_num, fun = function(x) as.integer(x %in% crop_codes))

response <- ifel(was_native == 1 & is_crop == 1, 1, 0)
response <- crop(response, sa_r) %>% mask(sa_r)

# check
# freq(response)
# layer value     count
# 1     1     0 685647465
# 2     1     1    203918

ecoregion_zones <- rasterize(ecoregions, sa_r, field = "ECOREGION1")

conversion_by_ecoregion <- zonal(response, ecoregion_zones, fun = "sum", na.rm = TRUE)
names(conversion_by_ecoregion) <- c("ecoregion", "n_conversions")
conversion_by_ecoregion <- conversion_by_ecoregion[order(conversion_by_ecoregion$n_conversions), ]

summary(conversion_by_ecoregion$n_conversions)
sum(conversion_by_ecoregion$n_conversions < 100)
sum(conversion_by_ecoregion$n_conversions == 0)

# also print the full table so we can see per-ecoregion counts
print(conversion_by_ecoregion)

# how many ecodistricts have fewer than 100 conversion events?
sum(conversion_by_ecoregion$n_conversions < 100, na.rm = TRUE)

# how many have zero?
sum(conversion_by_ecoregion$n_conversions == 0, na.rm = TRUE)


# Plot
# aggregate to lower resolution for plotting (factor = 10 means 10x coarser = 300m pixels)
response_plot <- aggregate(response, fact = 3, fun = "max")

# now convert - much more manageable size
response_df <- as.data.frame(response_plot, xy = TRUE) %>%
  filter(lyr.1 == 1)

# check how many rows
nrow(response_df)

ecoregions_sf <- st_as_sf(ecoregions)
ecozone_sf <- st_as_sf(ecozones)
ecoregions_sf <- st_intersection(ecoregions_sf, sa)
ecozone_sf    <- st_intersection(ecozone_sf, sa)

ecoregions_label <- ecoregions_sf %>%
  left_join(conversion_by_ecoregion, by = c("ECOREGION1" = "ecoregion")) %>%
  filter(n_conversions >= 1000)

map <- ggplot() +
  geom_tile(data = response_df, aes(x = x, y = y), fill = "darkred", alpha = 0.8) +
  geom_sf(data = ecoregions_sf, fill = NA, colour = "grey40", linewidth = 0.3) +
  geom_sf(data = ecozone_sf, fill = NA, colour = "black", linewidth = 0.7) +
  geom_sf_label(data = ecoregions_label, aes(label = ECOREGION1), 
                size = 2.5, fill = "white", alpha = 0.7,
                label.size = 0.1) +
  theme_void()

ggsave("grassland_conversion_map.png", plot = map,
       width = 10, height = 7, dpi = 300, bg = "white")

# End of response variable


#--- Prepare protected areas ----
pa <- read_sf_tagged("./processeddatasets", pa_filepath, pa_filename, "modelling-extent", driver = "GPKG")

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


write_sf_tagged(pa_hard, "./processeddatasets", pa_filepath, pa_filename,
                "modelling-extent", "hard", driver = "GPKG")

write_sf_tagged(pa_provpasture, "./processeddatasets", pa_filepath, pa_filename,
                "modelling-extent", "provpasture", driver = "GPKG")

write_sf_tagged(pa_easements, "./processeddatasets", pa_filepath, pa_filename,
                "modelling-extent", "easements", driver = "GPKG")


# Rasterize and save pa_hard
pa_hard_rast <- rasterize(vect(pa_hard), me_r,
                          field = 1,
                          background = 0,
                          touches = TRUE)

write_rast_tagged(pa_hard_rast, "./processeddatasets", pa_filepath, pa_filename,
                  "modelling-extent", "hard", "binary",
                  datatype = "INT1U")

# Rasterize and save pa_provpasture
pa_provpasture_rast <- rasterize(vect(pa_provpasture), me_r,
                                 field = 1,
                                 background = 0,
                                 touches = TRUE)

write_rast_tagged(pa_provpasture_rast, "./processeddatasets", pa_filepath, pa_filename,
                  "modelling-extent", "provpasture", "binary",
                  datatype = "INT1U")


# --- Proportion cropland ----
win_size_m     <- 3200
win_size_cells <- round(win_size_m / res(lulc2010)[1])
w              <- matrix(1, nrow = win_size_cells, ncol = win_size_cells)

#cropland_bin <- subst(lulc2010, 
#                      from = c("Cropland", "Land Converted to Cropland"),
#                      to   = c(1, 1),
#                      others = 0) %>%
#  crop(sa_r) %>%
#  mask(sa_r)

# cropland_codes <- lut_2010$Value[lut_2010$LU2010 %in% c("Cropland", "Land Converted to Cropland")]
# 
# cropland_bin <- app(lulc2010_num, fun = function(x) as.integer(x %in% cropland_codes)) %>%
#   crop(sa_r) %>%
#   mask(sa_r)

# create binary crop/non-crop raster
cropland2015_bin <- reclass_aafc(lulc2015, scheme = "binary", binary_class = "crop", check = TRUE)

# for a binary raster, the function mean is the same as proportion
cropland2015_prop <- focal(cropland2015_bin, w = w, fun = "mean",
                       na.policy = "omit", na.rm = TRUE,
                       filename = file.path(outdir, "cropland2015_prop.tif"),
                       overwrite = TRUE)

# QA/QC
# basic summary
summary(cropland_prop)

# should be between 0 and 1
minmax(cropland_prop)

# check for NAs
freq(is.na(cropland_prop))

# distribution of values
hist(cropland_prop, main = "Cropland proportion within 3.2km window",
     xlab = "Proportion", col = "grey80")

# visual check
plot(cropland_prop, main = "Cropland proportion")


# Cropland suitability
# rasterize ecodistricts to 300m template
ecodistrict_zones <- rasterize(ecodistricts, sa_r, field = "ECODISTRIC")

# calculate mean of crop_suit per ecodistrict
eco_means <- zonal(crop_suit, ecodistrict_zones, fun = "mean", na.rm = TRUE)
names(eco_means) <- c("ECODISTIC", "mean_suit")

# create raster of ecodistrict means
ecodistrict_means <- classify(ecodistrict_zones, eco_means)

# impute: replace NAs in crop_suit with ecodistrict mean, only within sa_r
crop_suit_filled <- ifel(is.na(crop_suit) & !is.na(sa_r), ecodistrict_means, crop_suit)

# quick check
freq(is.na(crop_suit_filled))

# visual check
plot(crop_suit_filled, main = "Cropland suitability")










