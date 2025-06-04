library(terra)
library(dplyr)

# Parameters
aafc_clipped_dir <- "D:/Projects/E104 - NU Grasslands/aafc_clipped"
res <- 30

# Load rasters
r2000 <- rast(file.path(aafc_clipped_dir, "LU2000_prairie_clipped.tif"))
r2010 <- rast(file.path(aafc_clipped_dir, "LU2010_prairie_clipped.tif"))
r2020 <- rast(file.path(aafc_clipped_dir, "LU2020_prairie_clipped.tif"))

r2000_raw <- as.numeric(r2000)
r2010_raw <- as.numeric(r2010)
r2020_raw <- as.numeric(r2020)

# Make sure all rasters are aligned (should be if created by same process)
stopifnot(compareGeom(r2000, r2010), compareGeom(r2010, r2020))

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

# Reclass rasters
r2000_ipcc <- classify(r2000_raw, rcl = ipcc_reclass, others = NA)
r2010_ipcc <- classify(r2010_raw, rcl = ipcc_reclass, others = NA)
r2020_ipcc <- classify(r2020_raw, rcl = ipcc_reclass, others = NA)

# Not sure if this is necessary
r2000_ipcc <- as.factor(r2000_ipcc)
r2010_ipcc <- as.factor(r2010_ipcc)
r2020_ipcc <- as.factor(r2020_ipcc)

# Create names list of rasters
ipcc_rasters <- list(
  "2000" = r2000_ipcc,
  "2010" = r2010_ipcc,
  "2020" = r2020_ipcc
)

###################################################
# Quantify LULC change from from_year to to_year ##
###################################################

# Define years
from_year <- 2000
to_year <- 2020
n_years <- to_year - from_year

# Stack input rasters
x <- c(ipcc_rasters[[as.character(from_year)]],
       ipcc_rasters[[as.character(to_year)]])
#x <- c(r2000_ipcc, r2020_ipcc)  # update later to use dynamic input if needed
trans_id <- terra::crosstab(x, useNA = FALSE, long = TRUE)

# Rename dynamically
colnames(trans_id) <- c(
  paste0("LU", from_year),
  paste0("LU", to_year),
  "n"
)

# Add class names
trans_labeled <- trans_id %>%
  mutate(across(starts_with("LU"), as.integer)) %>%
  left_join(ipcc_labels, by = setNames("IPCC_ID", paste0("LU", from_year))) %>%
  rename(From_Class = IPCC_Class) %>%
  left_join(ipcc_labels, by = setNames("IPCC_ID", paste0("LU", to_year))) %>%
  rename(To_Class = IPCC_Class)

# Summarize transitions
trans_summary <- trans_labeled %>%
  select(From_Class, To_Class, n) %>%
  mutate(
    km2 = n * res^2 / 1e6,
    km2_annual = km2 / n_years,
    total_area_km2 = sum(km2),
    percent_of_studyarea = km2 / total_area_km2 * 100
  ) %>%
  group_by(From_Class) %>%
  mutate(
    from_total_km2 = sum(km2),
    percent_of_from = km2 / from_total_km2 * 100
  ) %>%
  ungroup() %>%
  filter(From_Class != To_Class) %>%
  arrange(desc(km2))

# Summarize transitions 2
summary_barnes <- trans_labeled %>%
  select(From_Class, To_Class, n) %>%
  mutate(
    km2 = n * res^2 / 1e6,
    km2_annual = km2 / n_years,
    total_area_km2 = sum(km2),
    percent_of_studyarea = km2 / total_area_km2 * 100
  ) %>%
  group_by(From_Class) %>%
  mutate(
    from_total_km2 = sum(km2),
    percent_of_from = km2 / from_total_km2 * 100
  ) %>%
  ungroup() %>%
  filter(From_Class != To_Class) %>%
  arrange(desc(km2))


# Round & prepare clean version
trans_summary_clean <- trans_summary %>%
  select(From_Class, To_Class, km2, percent_of_studyarea, percent_of_from, km2_annual) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))

# Write to disk with dynamic filenames
write.csv(
  trans_summary,
  file.path(aafc_clipped_dir, paste0("lulc-change-", from_year, "-", to_year, "-full.csv")),
  row.names = FALSE
)

write.csv(
  trans_summary_clean,
  file.path(aafc_clipped_dir, paste0("lulc-change-", from_year, "-", to_year, ".csv")),
  row.names = FALSE
)