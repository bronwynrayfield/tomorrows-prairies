# =============================================================================
# presentation-candidate-locations.R
# =============================================================================
#
# Purpose:
#   Identifies and visualizes candidate grassland conversion locations for
#   presentation narrative. Finds large contiguous native-to-annual-cropland
#   conversion patches in Aspen Parkland with high surrounding cropland context.
#   Not part of the main analysis pipeline.
#
# Outputs:
#   - candidate_locations.csv — top 20 conversion patch candidates with
#     coordinates and cropland proportion
#   - fig_conversion_example.png — six-panel figure showing field and
#     landscape scale before/after maps for selected candidate
#
# Notes:
#   - Candidate 4 (clump 10787, 68 ha, 42% annual cropland) selected as
#     presentation example — central Saskatchewan, Aspen Parkland
#   - Candidate 5 excluded — possible tame grass misclassification (hay
#     bales visible in 2020 Google Earth imagery)
#   - Candidates 15, 17, 19 are in Alberta (lon ~-114) — outside typical
#     Aspen Parkland extent, likely edge effects
#   - This script is for presentation preparation only and is not part of
#     the main modelling pipeline
#
# Dependencies:
#   - 0-setup.R
#   - 6-response.R outputs (response_30m.tif)
#   - 5-spatial-predictors.R outputs (crop2015_prop_3200m_300m.tif)
#   - 3-crop-data.R outputs (AAFC LUTS 2010, 2020 cropped)
#
# Author: [your name]
# Date:   [date]
# =============================================================================

source(file.path("scripts/landcover", "0-setup.R"))
source(file.path("scripts/landcover/functions", "0-functions-general.R"))
source(file.path("scripts/landcover/functions", "reclass-aafc-luts.R"))

# ---- Load files ----
lulc2010_full <- rast(tagged_path(file.path("./processeddatasets", aafc_filepath),
                                  "cropped", aafc2010_filename, "modelling-extent"))
lulc2020_full <- rast(tagged_path(file.path("./processeddatasets", aafc_filepath),
                                  "cropped", aafc2020_filename, "modelling-extent"))
response_full <- rast(tagged_path("./processeddatasets", "response",
                                  "response.tif", "30m"))
crop_prop_annual <- rast(file.path("./processeddatasets", "predictors", "300m",
                                   "crop2015_prop_3200m_300m.tif"))
ecoregions <- st_read(tagged_path("./processeddatasets", ecoregions_filepath,
                                  ecoregions_filename, "modelling-extent"))

# ---- Crop to Aspen Parkland ----
ap_boundary  <- ecoregions |> filter(ECOREGION1 == "Aspen Parkland")
response_ap  <- mask(crop(response_full, ap_boundary), ap_boundary)

# ---- Find contiguous conversion patches ----
conversion_ap  <- ifel(response_ap == 1, 1, NA)
clumps_ap      <- patches(conversion_ap, directions = 8, zeroAsNA = TRUE)
clumps_ap_freq <- freq(clumps_ap) |> arrange(desc(count))

# ---- Get centroids and cropland proportion for top 20 clumps ----
crop_prop_ap <- mask(crop(crop_prop_annual, ap_boundary), ap_boundary)
top_ids      <- head(clumps_ap_freq, 20)$value

results <- lapply(top_ids, function(id) {
  cells <- which(values(clumps_ap) == id)
  xy    <- xyFromCell(clumps_ap, cells)
  data.frame(
    clump_id = id,
    x        = mean(xy[, 1]),
    y        = mean(xy[, 2]),
    n_pixels = length(cells),
    area_ha  = length(cells) * 900 / 10000
  )
}) |> bind_rows()

# Extract annual cropland proportion at centroids
cents_vect        <- vect(results, geom = c("x", "y"), crs = crs(clumps_ap))
results$crop_prop <- terra::extract(crop_prop_ap, cents_vect)[, 2]

# Add lat/long columns
results_latlong <- cents_vect |>
  project("EPSG:4326") |>
  as.data.frame(geom = "XY") |>
  rename(lon = x, lat = y)

results <- results |>
  mutate(lon        = results_latlong$lon,
         lat        = results_latlong$lat,
         google_maps = paste0(lat, ", ", lon)) |>
  arrange(desc(crop_prop))

# Save candidate locations
write.csv(results, "./output/presentation/candidate_locations.csv",
          row.names = FALSE)

print(results[, c("clump_id", "area_ha", "crop_prop", "lon", "lat", "google_maps")])
> print(results[, c("clump_id", "area_ha", "crop_prop", "lon", "lat", "google_maps")])
# clump_id area_ha  crop_prop        lon      lat                         google_maps
# 1     27635   31.77 0.62735540  -99.26557 49.97096 49.9709643791721, -99.2655714196983
# 2     26499   31.41 0.51074380  -99.35746 50.07349 50.0734918207551, -99.3574648499672
# 3     19731  112.68 0.48743802 -102.63745 50.54843 50.5484304471983, -102.637447965444
# 4     10787   68.13 0.42396694 -103.32356 51.47929 51.4792905866451, -103.323563144099
# 5     20253   21.78 0.39107439 -101.63476 50.56278 50.5627816921504, -101.634756100962
# 6     29801   52.92 0.38661158  -99.62677 49.69087 49.6908666030555, -99.6267693715946
# 7     29772   22.86 0.31206611  -99.00325 49.70692  49.706917117671, -99.0032487234483
# 8     29746   28.62 0.30768594  -99.59250 49.69763 49.6976267036341, -99.5924985784103
# 9     18756   76.59 0.30355373 -101.40569 50.65003 50.6500257540512, -101.405692826659
# 10     9716   50.22 0.27685949 -103.13313 51.64819 51.6481852365409, -103.133134393711
# 11    26335   20.61 0.21206611  -99.40186 50.09099 50.0909915368551, -99.4018588435319
# 12    28520   26.01 0.21165289 -100.74535 49.82071 49.8207066375532, -100.745346813434
# 13    17958   24.12 0.20181818 -101.42009 50.71704 50.7170424581194, -101.420087803852
# 14    29690   41.22 0.18041322  -99.58098 49.70491 49.7049144953746, -99.5809795236596
# 15     9062   43.20 0.15834710 -114.30942 50.97587 50.9758695620673, -114.309422866711
# 16    24567   47.97 0.13520661 -101.10187 50.21309  50.213092462602, -101.101872403594
# 17     8653   35.82 0.13256198 -114.49385 51.12904 51.1290386640466, -114.493850480166
# 18    13798   32.22 0.12793389 -102.51473 51.11844 51.1184410341716, -102.514726082442
# 19     8697   23.94 0.11545455 -114.27560 51.14434 51.1443389878409, -114.275599891529
# 20    20924   33.75 0.08314049 -101.64720 50.52330   50.5233041057238, -101.6471955173

# ---- Visualize selected candidate (Candidate 4) ----
# clump 10787, 68 ha, 42% annual cropland, central Saskatchewan
c4_pt  <- vect(data.frame(x = results$x[4], y = results$y[4]),
               geom = c("x", "y"), crs = crs(clumps_ap))

buf_field     <- buffer(c4_pt, width = 3000)
buf_landscape <- buffer(c4_pt, width = 15000)

# Crop to buffers
lulc2010_field <- mask(crop(lulc2010_full, buf_field),     buf_field)
lulc2020_field <- mask(crop(lulc2020_full, buf_field),     buf_field)
response_field <- mask(crop(response_ap,   buf_field),     buf_field)
lulc2010_land  <- mask(crop(lulc2010_full, buf_landscape), buf_landscape)
lulc2020_land  <- mask(crop(lulc2020_full, buf_landscape), buf_landscape)
response_land  <- mask(crop(response_ap,   buf_landscape), buf_landscape)

# Colour scheme
lulc_cols_4 <- c(
  "1" = "#5AAE61",  # native grassland
  "3" = "#FEA500",  # cropland (incl tame and annual)
  "4" = "#7B3294",  # settlement
  "5" = "#D3D3D3"   # other
)

# ---- Save figure ----
dir.create("./output/presentation", recursive = TRUE, showWarnings = FALSE)

png("./output/presentation/fig_conversion_example.png",
    width = 3000, height = 2000, res = 300)

par(mfrow = c(2, 3))

plot(reclass_aafc(lulc2010_land, scheme = "4class"),
     main = "Landscape 2010", col = lulc_cols_4,
     legend = FALSE, axes = FALSE, box = FALSE)
plot(reclass_aafc(lulc2020_land, scheme = "4class"),
     main = "Landscape 2020", col = lulc_cols_4,
     legend = FALSE, axes = FALSE, box = FALSE)
plot(response_land, main = "Conversions 2010-2020",
     col = c("grey90", "darkred"),
     legend = FALSE, axes = FALSE, box = FALSE)
plot(reclass_aafc(lulc2010_field, scheme = "4class"),
     main = "Field 2010", col = lulc_cols_4,
     legend = FALSE, axes = FALSE, box = FALSE)
plot(reclass_aafc(lulc2020_field, scheme = "4class"),
     main = "Field 2020", col = lulc_cols_4,
     legend = FALSE, axes = FALSE, box = FALSE)
plot(response_field, main = "Conversions 2010-2020",
     col = c("grey90", "darkred"),
     legend = FALSE, axes = FALSE, box = FALSE)

dev.off()

par(mfrow = c(1, 1))



# Response map
# Get label point for largest polygon per ecoregion
ecoregions_labels <- ecoregions |>
  filter(ECOZONE_ID == 10) |>
  group_by(ECOREGION1) |>
  slice_max(Shape_Area, n = 1) |>  # keep only largest polygon per ecoregion
  ungroup() |>
  st_centroid()

map_response <- ggplot() +
  geom_sf(data = ecoregions |> filter(ECOZONE_ID == 10),
          aes(fill = ECOREGION1), colour = "white", linewidth = 0.4, alpha = 0.2) +
  geom_tile(data = response_df, aes(x = x, y = y),
            fill = "darkred", alpha = 0.9) +
  geom_sf(data = ecoregions |> filter(ECOZONE_ID == 10),
          fill = NA, colour = "grey40", linewidth = 0.4) +
  geom_sf_label(data = ecoregions_labels,
                aes(label = ECOREGION1),
                size = 3, fill = "white", alpha = 0.8,
                label.size = 0.1, fontface = "bold") +
  scale_fill_brewer(palette = "Set3") +
  theme_void() +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white", colour = NA),
        plot.background  = element_rect(fill = "white", colour = NA))

ggsave("./output/presentation/fig_response_map.png",
       plot = map_response, width = 10, height = 7, dpi = 300, bg = "white")




# Add circle for candidate 4 location
c4_latlong <- vect(data.frame(x = results$x[4], y = results$y[4]),
                   geom = c("x", "y"), crs = crs(clumps_ap)) |>
  project("EPSG:4326") |>
  as.data.frame(geom = "XY")

# Create circle buffer around candidate 4
c4_sf <- st_as_sf(data.frame(x = results$x[4], y = results$y[4]),
                  coords = c("x", "y"), crs = crs(clumps_ap)) |>
  st_buffer(dist = 15000) |>  # 15km radius to match landscape buffer
  st_transform(st_crs(ecoregions))

map_response_circle <- ggplot() +
  geom_sf(data = ecoregions |> filter(ECOZONE_ID == 10),
          aes(fill = ECOREGION1), colour = "white", linewidth = 0.4, alpha = 0.2) +
  geom_tile(data = response_df, aes(x = x, y = y),
            fill = "darkred", alpha = 0.9) +
  geom_sf(data = ecoregions |> filter(ECOZONE_ID == 10),
          fill = NA, colour = "grey40", linewidth = 0.4) +
  geom_sf(data = c4_sf, fill = NA, colour = "black", linewidth = 1) +
  geom_sf_label(data = ecoregions_labels,
                aes(label = ECOREGION1),
                size = 3, fill = "white", alpha = 0.8,
                label.size = 0.1, fontface = "bold") +
  scale_fill_brewer(palette = "Set3") +
  theme_void() +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white", colour = NA),
        plot.background  = element_rect(fill = "white", colour = NA))

ggsave("./output/presentation/fig_response_map_circle.png",
       plot = map_response_circle, width = 10, height = 7, dpi = 300, bg = "white")



# Landscape scale inset — 2010 and 2020 side by side
png("./output/presentation/fig_landscape_inset.png",
    width = 2000, height = 1000, res = 300)

par(mfrow = c(1, 2), mar = c(0.5, 0.5, 1.5, 0.5))

plot(reclass_aafc(lulc2010_land, scheme = "4class"),
     main = "2010", col = lulc_cols_4,
     legend = FALSE, axes = FALSE, box = FALSE)

plot(reclass_aafc(lulc2020_land, scheme = "4class"),
     main = "2020", col = lulc_cols_4,
     legend = FALSE, axes = FALSE, box = FALSE)

dev.off()

par(mfrow = c(1, 1))

# ---- Figure 3: Landscape 2010 only with field circle ----
c4_field_sf <- st_as_sf(data.frame(x = results$x[4], y = results$y[4]),
                        coords = c("x", "y"), crs = crs(clumps_ap)) |>
  st_buffer(dist = 3000) |>
  st_transform(st_crs(ecoregions))

png("./output/presentation/fig_landscape_2010.png",
    width = 1200, height = 1200, res = 300)
par(mar = c(0.5, 0.5, 1.5, 0.5))
plot(reclass_aafc(lulc2010_land, scheme = "4class"),
     main = "2010", col = lulc_cols_4,
     legend = FALSE, axes = FALSE, box = FALSE)
plot(st_geometry(c4_field_sf), add = TRUE,
     border = "black", col = NA, lwd = 2)
dev.off()

# ---- Figure 4: Landscape 2010 and 2020 side by side with field circle ----
png("./output/presentation/fig_landscape_both.png",
    width = 2400, height = 1200, res = 300)
par(mfrow = c(1, 2), mar = c(0.5, 0.5, 1.5, 0.5))
plot(reclass_aafc(lulc2010_land, scheme = "4class"),
     main = "2010", col = lulc_cols_4,
     legend = FALSE, axes = FALSE, box = FALSE)
plot(st_geometry(c4_field_sf), add = TRUE,
     border = "black", col = NA, lwd = 2)
plot(reclass_aafc(lulc2020_land, scheme = "4class"),
     main = "2020", col = lulc_cols_4,
     legend = FALSE, axes = FALSE, box = FALSE)
plot(st_geometry(c4_field_sf), add = TRUE,
     border = "black", col = NA, lwd = 2)
dev.off()
par(mfrow = c(1, 1))

# ---- Figure 5: Field 2010 and 2020 side by side ----
png("./output/presentation/fig_field_both.png",
    width = 2400, height = 1200, res = 300)
par(mfrow = c(1, 2), mar = c(0.5, 0.5, 1.5, 0.5))
plot(reclass_aafc(lulc2010_field, scheme = "4class"),
     main = "2010", col = lulc_cols_4,
     legend = FALSE, axes = FALSE, box = FALSE)
plot(reclass_aafc(lulc2020_field, scheme = "4class"),
     main = "2020", col = lulc_cols_4,
     legend = FALSE, axes = FALSE, box = FALSE)
dev.off()
par(mfrow = c(1, 1))

# ---- Figure 8: Conversion panel ----
png("./output/presentation/fig_field_conversions.png",
    width = 1200, height = 1200, res = 300)
par(mar = c(0.5, 0.5, 1.5, 0.5))
plot(response_field, main = "Native → Crop 2010-2020",
     col = c("grey90", "darkred"),
     legend = FALSE, axes = FALSE, box = FALSE)
dev.off()




#study area map with two circles
map_response_circles <- ggplot() +
  geom_sf(data = ecoregions |> filter(ECOZONE_ID == 10),
          aes(fill = ECOREGION1), colour = "white", linewidth = 0.4, alpha = 0.2) +
  geom_tile(data = response_df, aes(x = x, y = y),
            fill = "darkred", alpha = 0.9) +
  geom_sf(data = ecoregions |> filter(ECOZONE_ID == 10),
          fill = NA, colour = "grey40", linewidth = 0.4) +
  # Landscape circle
  geom_sf(data = c4_landscape_sf, fill = NA,
          colour = "#FEA500", linewidth = 0.9) +
  # Field circle
  geom_sf(data = c4_field_sf, fill = NA,
          colour = "#FEA500", linewidth = 0.9) +
  geom_sf_label(data = ecoregions_labels,
                aes(label = ECOREGION1),
                size = 3, fill = "white", alpha = 0.8,
                label.size = 0.1, fontface = "bold") +
  scale_fill_brewer(palette = "Set3") +
  theme_void() +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white", colour = NA),
        plot.background  = element_rect(fill = "white", colour = NA))

ggsave("./output/presentation/fig_response_map_circles.png",
       plot = map_response_circles, width = 10, height = 7, dpi = 300, bg = "white")

