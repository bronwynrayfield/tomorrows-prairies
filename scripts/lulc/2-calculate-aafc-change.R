##################################################
# Quantify LULC change from from_year to to_year
# Summarize LULC change and save tables

####--- Setup ---####
source(file.path("scripts", "0-setup.R"))

####--- Directories ---#### 
aafc_cropped_dir <- file.path(proc_lulc_data_dir, "aafc/cropped")

####--- Load rasters ---
# Load cropped AAFC rasters
r2000 <- rast(file.path(aafc_cropped_dir, "LU2000_prairie.tif"))
r2010 <- rast(file.path(aafc_cropped_dir, "LU2010_prairie.tif"))
r2020 <- rast(file.path(aafc_cropped_dir, "LU2020_prairie.tif"))

# Load vegetation zones
veg_zones <- st_read(veg_zones_file)

####--- Prep spatial data ---####
# Make sure AAFC rasters are aligned (should be if created by same process)
stopifnot(compareGeom(r2000, r2010), compareGeom(r2010, r2020))

# Give vegetation zones Same CRS as AAFC rasters
veg_zones <- st_transform(veg_zones, crs(r2000))

# Rasterize vegetation zones to match IPCC raster resolution and extent
veg_rast <- terra::rasterize(vect(veg_zones), ipcc_rasters[[1]], field = "level_2")

####--- Reclass rasters ---####
# Make reclass table to go from AAFC classes to IUCN classes
rat <- cats(r2020)[[1]]
ipcc_lookup <- rat[, c("Value", "IPCC_Class")]
# str(ipcc_lookup)
ipcc_classes <- unique(ipcc_lookup$IPCC_Class)
ipcc_ids <- seq_along(ipcc_classes)

ipcc_reclass <- data.frame(
  Value = ipcc_lookup$Value,
  IPCC_ID = match(ipcc_lookup$IPCC_Class, ipcc_classes)
)

# Keep class labels handy for later lookup
ipcc_labels <- data.frame(
  IPCC_ID = ipcc_ids,
  IPCC_Class = ipcc_classes
)

# Reclass
r2000_ipcc <- classify(as.numeric(r2000), rcl = ipcc_reclass, others = NA)
r2010_ipcc <- classify(as.numeric(r2010), rcl = ipcc_reclass, others = NA)
r2020_ipcc <- classify(as.numeric(r2020), rcl = ipcc_reclass, others = NA)

# Create named list of rasters
ipcc_rasters <- list(
  "2000" = r2000_ipcc,
  "2010" = r2010_ipcc,
  "2020" = r2020_ipcc
)


####--- Quantify LULC change from from_year to to_year ---####
# Define years
from_year <- 2010
to_year <- 2020
n_years <- to_year - from_year

## Crosstabulate grassland transitions by vegetation zone
# Stack the from, to, and vegetation zone rasters
x <- c(ipcc_rasters[[as.character(from_year)]],
       ipcc_rasters[[as.character(to_year)]],
       veg_rast)

# Crosstab across all 3 layers
transitions_by_zone <- terra::crosstab(x, useNA = FALSE, long = TRUE)

# Rename dynamically
colnames(transitions_by_zone) <- c(
  paste0("LU", from_year),
  paste0("LU", to_year),
  "VegZone",
  "n")

# Add class names
trans_by_zone_labelled <- transitions_by_zone %>%
  mutate(across(starts_with("LU"), as.integer)) %>%
  left_join(ipcc_labels, by = setNames("IPCC_ID", paste0("LU", from_year))) %>%
  rename(From_Class = IPCC_Class) %>%
  left_join(ipcc_labels, by = setNames("IPCC_ID", paste0("LU", to_year))) %>%
  rename(To_Class = IPCC_Class)

####--- Summarize results ---####
# Set output unit: "ha" for hectares, "km2" for square kilometres
output_unit <- "km2"  # or "km2"
unit_conversion <- if (output_unit == "ha") 100 else 1e6

# Summarize all transitions across the full study area
trans_summary <- trans_by_zone_labelled %>%
  group_by(From_Class, To_Class) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  mutate(
    area = n * res^2 / unit_conversion,
    area_annual = area / n_years
  ) %>%
  mutate(
    total_area = sum(area)  # Total area of all transitions
  ) %>%
  group_by(From_Class) %>%
  mutate(
    from_total = sum(area),                       # Total area from each From_Class
    prop_of_from = area / from_total     # proportion of From_Class transitions
  ) %>%
  ungroup() %>%
  mutate(
    prop_of_studyarea = area / total_area
  ) %>%
  filter(From_Class != To_Class) %>%
  arrange(desc(area))

# Round & prepare clean version
trans_summary_clean <- trans_summary %>%
  select(From_Class, To_Class, area, prop_of_studyarea, prop_of_from, area_annual) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))


# Summarize Grassland change by vegetation zone
veg_zone_grass <- trans_by_zone_labelled %>%
  mutate(area = n * (res^2) / unit_conversion) %>%
  group_by(VegZone) %>%
  summarise(
    VegZone_Area = sum(area),
    Grass = sum(area[From_Class == "Grassland" & To_Class == "Grassland"]),
    Loss_to_Crop = sum(area[From_Class == "Grassland" & To_Class == "Cropland"]),
    Loss_to_Dev = sum(area[From_Class == "Grassland" & To_Class == "Settlement"]),
    Loss_to_Crop_or_Dev = sum(area[From_Class == "Grassland" & To_Class %in% c("Cropland", "Settlement")])
  ) %>%
  mutate(
    Prop_Loss_to_Crop = Loss_to_Crop / Grass,
    Prop_Loss_to_Dev = Loss_to_Dev / Grass,
    Prop_Loss_to_Crop_or_Dev = Loss_to_Crop_or_Dev / Grass
  )

# Round & prepare clean version
veg_zone_grass_clean <- veg_zone_grass %>%
   mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
   arrange(desc(Loss_to_Crop_or_Dev))


####--- Write outputs ---####
# Write to disk with dynamic filenames
# Raw transitions
write.csv(
  trans_by_zone_labelled,
  file.path(out_lulc_historic_dir, paste0("aafc-raw-transitions-by-zone-", from_year, "-", to_year, ".csv")),
  row.names = FALSE
)
# Study area - full summary
write.csv(
  trans_summary,
  file.path(out_lulc_historic_dir, paste0("aafc-lulc-change-", from_year, "-", to_year, "-full-", output_unit, ".csv")),
  row.names = FALSE
)
# Study area - clean summary
write.csv(
  trans_summary_clean,
  file.path(out_lulc_historic_dir, paste0("aafc-lulc-change-", from_year, "-", to_year, "-", output_unit, ".csv")),
  row.names = FALSE
)
# By vegetation zone - clean summary
write.csv(
  veg_zone_grass_clean,
  file.path(out_lulc_historic_dir, paste0("aafc-grass-change-by-zone-", from_year, "-", to_year, "-", output_unit, ".csv")),
  row.names = FALSE
)
