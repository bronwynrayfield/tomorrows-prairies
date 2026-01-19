library(terra)

# Example class codes
grassland_code <- 60
cropland_code <- 50
settlement_code <- 20

# Load rasters
r2010 <- rast("lulc_2010.tif")
r2020 <- rast("lulc_2020.tif")

# Identify transitions
grass_to_crop <- (r2010 == grassland_code) & (r2020 == cropland_code)
grass_to_settle <- (r2010 == grassland_code) & (r2020 == settlement_code)

# Combine into one raster
transition_map <- classify(grass_to_crop, cbind(1, 1))  # Grass → Crop = 1
transition_map[grass_to_settle] <- 2                    # Grass → Settle = 2

# Plot
plot(transition_map, col = c("white", "black", "red"),
     legend = FALSE, main = "Grassland Transitions")
legend("topright", legend = c("Grass → Cropland", "Grass → Settlement"),
       fill = c("black", "red"))


# Aggregate by factor of 5 (i.e., coarser 150 m resolution)
transition_map_small <- aggregate(transition_map, fact = 5, fun = "modal")

transition_df <- as.data.frame(transition_map_small, xy = TRUE, na.rm = TRUE)
colnames(transition_df)[3] <- "transition"
transition_df <- transition_df[transition_df$transition > 0, ]




# Convert to data frame for ggplot
transition_df <- as.data.frame(transition_map, xy = TRUE)
colnames(transition_df)[3] <- "transition"

# Drop NA and unneeded zeros
transition_df <- transition_df[!is.na(transition_df$transition) & transition_df$transition > 0, ]

# Plot
ggplot(transition_df) +
  geom_raster(aes(x = x, y = y, fill = factor(transition))) +
  scale_fill_manual(
    name = "Transition",
    values = c("1" = "forestgreen", "2" = "firebrick"),
    labels = c("Grass → Cropland", "Grass → Settlement")
  ) +
  coord_equal() +
  theme_minimal() +
  labs(title = "Grassland Transitions (2010–2020)", x = "Longitude", y = "Latitude")


ggplot() +
  geom_raster(data = transition_df, aes(x = x, y = y, fill = factor(transition))) +
  geom_sf(data = veg_zones, fill = NA, color = "black", size = 0.3) +  # outlines only
  scale_fill_manual(
    name = "Transition",
    values = c("1" = "forestgreen", "2" = "firebrick"),
    labels = c("Grass → Cropland", "Grass → Settlement")
  ) +
  coord_sf() +
  theme_minimal() +
  labs(
    title = "Grassland Transitions with Vegetation Zones",
    x = "Longitude", y = "Latitude"
  )




####################
# 1. Filter to just cells with transitions
transition_df <- transition_df[transition_df$transition > 0, ]

# 2. Get bounding box of just the transition pixels
bbox_transitions <- st_as_sf(transition_df, coords = c("x", "y"), crs = crs(veg_zones)) |>
  st_bbox() |> st_as_sfc()

# 3. Crop vegetation zones to that extent (with buffer)
buffered_extent <- st_buffer(bbox_transitions, dist = 10000)  # 10 km buffer
veg_zones_crop <- st_crop(veg_zones, buffered_extent)

# 4. Plot
ggplot() +
  geom_raster(data = transition_df, aes(x = x, y = y, fill = factor(transition))) +
  geom_sf(data = veg_zones_crop, fill = NA, color = "black", size = 0.3) +
  scale_fill_manual(
    name = "Transition",
    values = c("1" = "forestgreen", "2" = "firebrick"),
    labels = c("Grass → Cropland", "Grass → Settlement")
  ) +
  coord_sf(xlim = st_bbox(buffered_extent)[c("xmin", "xmax")],
           ylim = st_bbox(buffered_extent)[c("ymin", "ymax")]) +
  theme_minimal() +
  labs(
    title = "Zoomed-in Grassland Transitions",
    x = "Longitude", y = "Latitude"
  )



##############
library(sf)
library(dplyr)

# City coordinates (decimal degrees)
cities <- tribble(
  ~name,          ~lon,        ~lat,
  "Edmonton",     -113.4909,   53.5461,
  "Calgary",      -114.0719,   51.0447,
  "Regina",       -104.6189,   50.4452,
  "Saskatoon",    -106.6702,   52.1332,
  "Prince Albert",-105.7513,   53.2032,
  "Winnipeg",     -97.1384,    49.8951,
  "Brandon",      -99.9539,    49.8485
)

cities_sf <- cities %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  st_transform(st_crs(veg_zones))

ggplot() +
  geom_point(
    data = transition_df,
    aes(x = X, y = Y, color = factor(transition)),
    shape = 15, size = 0.25
  ) +
  geom_sf(data = veg_zones, fill = NA, color = "black", linewidth = 0.3) +
  geom_sf(data = cities_sf, size = 4, shape = 21, fill = "yellow", color = "black") +
  geom_sf_text(data = cities_sf, aes(label = name),
               size = 5, fontface = "bold", nudge_y = 30000) +
  scale_color_manual(
    name = "Transition",
    values = c("1" = "black", "2" = "firebrick"),
    labels = c("Grass → Cropland", "Grass → Settlement")
  ) +
  coord_sf(
    xlim = st_bbox(bbox_buffered)[c("xmin", "xmax")],
    ylim = st_bbox(bbox_buffered)[c("ymin", "ymax")]
  ) +
  theme_minimal() +
  labs(
    title = "Grassland Transitions (2010–2020) in Prairie Ecozone",
    x = "Longitude", y = "Latitude"
  )


library(rnaturalearth)
library(rnaturalearthdata)
library(sf)

# Get rivers and lakes as sf
rivers <- ne_download(scale = 50, type = "rivers_lake_centerlines", category = "physical", returnclass = "sf")
lakes <- ne_download(scale = 50, type = "lakes", category = "physical", returnclass = "sf")

# Transform to match your map's CRS
rivers <- st_transform(rivers, crs = st_crs(veg_zones))
lakes <- st_transform(lakes, crs = st_crs(veg_zones))

# Use veg_zones_crop extent for clipping
study_extent <- st_bbox(veg_zones) |> st_as_sfc() |> st_as_sf()
study_extent <- st_transform(study_extent, crs = st_crs(rivers))


rivers_clipped <- st_intersection(rivers, study_extent)
lakes_clipped <- st_intersection(lakes, study_extent)


ggplot() +
  geom_point(data = transition_df, aes(x = X, y = Y, color = factor(transition)),
             shape = 15, size = 0.25) +
  
  # Add cropped rivers and lakes
  geom_sf(data = lakes_clipped, fill = "lightblue", color = NA, alpha = 0.4) +
  geom_sf(data = rivers_clipped, color = "blue", size = 0.3) +
  
  # Vegetation zones
  geom_sf(data = veg_zones, fill = NA, color = "darkgrey", linewidth = 0.3) +
  
  # Cities
  geom_sf(data = cities_sf, shape = 21, fill = "yellow", color = "black", size = 2) +
  geom_sf_text(data = cities_sf, aes(label = name), size = 3, fontface = "bold", nudge_y = 30000) +
  
  # Color scale
  scale_color_manual(
    name = "Transition",
    values = c("1" = "black", "2" = "firebrick"),
    labels = c("Grass → Cropland", "Grass → Settlement")
  ) +
  
  theme_minimal() +
  labs(
    title = "Grassland Transitions (2010–2020) with Lakes, Rivers, and Vegetation Zones",
    x = "Longitude", y = "Latitude"
  )




