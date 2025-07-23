##################################################
# Quantify LULC change from from_year to to_year
# Summarize LULC change and save tables

####--- Setup ---####
source(file.path("scripts", "0-setup.R"))

####--- Directories ---#### 
aafc_cropped_dir <- file.path(proc_lulc_data_dir, "aafc/cropped")

####--- Definitions ----####
# Define years
from_year <- 2010
to_year <- 2020
n_years <- to_year - from_year

####--- Load data and prepare ---
# Load cropped AAFC rasters
#r2000_aafc <- rast(file.path(aafc_cropped_dir, "AAFC2000_sa_nobuffer.tif"))
r2010_aafc <- rast(file.path(aafc_cropped_dir, "AAFC2010_sa_nobuffer.tif"))
#r2015_aafc <- rast(file.path(aafc_cropped_dir, "AAFC2015_sa_nobuffer.tif"))
r2020_aafc <- rast(file.path(aafc_cropped_dir, "AAFC2020_sa_nobuffer.tif"))

# Load vegetation zones
veg_zones <- st_read(veg_zones_file)

# Load AAFC crosswalk
aafc_crosswalk_2010 <- read.csv(aafc_crosswalk_file_2010)
aafc_crosswalk_2015 <- read.csv(aafc_crosswalk_file_2015)

if(from_year == 2010){ 
  aafc_crosswalk <- aafc_crosswalk_2010
} else if(from_year == 2015){ 
    aafc_crosswalk <- aafc_crosswalk_2015
    }

# Create reclass table with IDs
aafc_reclass <- aafc_crosswalk %>%
  select(Value, Our_ID)

# Keep class labels (remove duplicates)
labels <- aafc_crosswalk %>%
  transmute(ID = Our_ID, Class = Our_Class) %>%
  distinct()


####--- Prep spatial data ---####
# Make sure AAFC rasters are aligned (should be if created by same process)
#stopifnot(compareGeom(r2000, r2010), compareGeom(r2010, r2020))


####--- Reclass rasters ---####
#r2000 <- classify(as.numeric(r2000_aafc), rcl = aafc_reclass, others = NA)
r2010 <- classify(as.numeric(r2010_aafc), rcl = aafc_reclass, others = NA)
#r2015 <- classify(as.numeric(r2015_aafc), rcl = ipcc_reclass, others = NA)
r2020 <- classify(as.numeric(r2020_aafc), rcl = aafc_reclass, others = NA)

# Create named list of rasters
lulc_rasters <- list(
  "2010" = r2010,
#  "2015" = r2015,
  "2020" = r2020
)

# Give vegetation zones Same CRS as AAFC rasters
veg_zones <- st_transform(veg_zones, crs(aafc_rasters[[1]]))

# Rasterize vegetation zones to match LULC raster resolution and extent
veg_rast <- terra::rasterize(vect(veg_zones), aafc_rasters[[1]], field = "id_2")

####--- Quantify LULC change from from_year to to_year ---####

## Crosstabulate grassland transitions by vegetation zone
# Stack the from, to, and vegetation zone rasters
x <- c(lulc_rasters[[as.character(from_year)]],
       lulc_rasters[[as.character(to_year)]],
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
  left_join(labels, by = setNames("ID", paste0("LU", from_year))) %>%
  rename(From_Class = Class) %>%
  left_join(labels, by = setNames("ID", paste0("LU", to_year))) %>%
  rename(To_Class = Class)

trans_by_zone_labelled <- read.csv(file.path(out_lulc_historic_dir, paste0("aafc-raw-transitions-by-zone-", from_year, "-", to_year, ".csv")))

####--- Summarize results ---####
# Set output unit: "ha" for hectares, "km2" for square kilometres
output_unit <- "km2"  # or "km2"
unit_conversion <- if (output_unit == "ha") 1e4 else 1e6

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
    Grass = sum(area[From_Class == "Grassland"]),
    Loss_to_Crop = sum(area[From_Class == "Grassland" & To_Class == "Cropland"])/n_years,
    Loss_to_Dev = sum(area[From_Class == "Grassland" & To_Class == "Settlement"]/n_years),
    Loss_to_Crop_or_Dev = sum(area[From_Class == "Grassland" & To_Class %in% c("Cropland", "Settlement")]/n_years)
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
