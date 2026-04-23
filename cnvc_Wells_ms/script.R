
####################################################################################
##         ANALYSIS FOR CNVC MANUSCRIPT       ##
####################################################################################

##----------------------- Alison Long - April 2026 -------------------------

# Product 1. Separate distribution maps of each CNVC group, where CNVC plots for a group 
#            are shown on a map of the group.  We'd list to see if coding the plots for 
#            the group by alliance would highlight variation within the group.
# Product 2. Accuracy assessment on the CNVC groups. A 2000 cutoff would be appropriate 
#            (assuming that after that date the locations would be more accurate to due to 
#            removal of selective availability). However, there are many plots sampled in the 
#            1990s, so run the accuracy assessment using plots with 2 cut-off dates, 2000 and 
#            1990. Depending on how that works out, we could consider running the AA at the 
#            macrogroup level (most relevant to M051).


# ----------------------------
# 1) Setup
# ----------------------------
# load packages
library(tidyverse)
library(sf)
library(terra)
library(ggplot2)
library(scales)

# ----------------------------
# 2) Directories and file paths 
# ----------------------------

# Directories
map_dir    <- "R:/CanadianGrasslands/existingdatasets/Grassland/NatureServe_CNVC/updated version 23Oct2025/NS_CNVC_TIFF_COMPRESSED"
pts_dir    <- "R:/CanadianGrasslands/WG_Classification/CNVC_veg_data_analysis_Wells/veg_plots_Mar2026_v3/Prairie_Province_Veg_Plots_2026-03-31"
ecoreg_dir <- "R:/CanadianGrasslands/processeddatasets/boundaries/national-ecological-framework-for-canada/Ecoregions"
pe_dir     <- "R:/CanadianGrasslands/processeddatasets/boundaries/national-ecological-framework-for-canada/Ecozones"
out_dir    <- "R:/CanadianGrasslands/WG_Classification/CNVC_veg_data_analysis_Wells/results"

# File paths
map_file    <- file.path(map_dir, "NatureServe_CNVC_COMPRESSED_20251023.tif")
ecoreg_file <- file.path(ecoreg_dir, "ecoregions_clean.shp")
pe_file     <- file.path(pe_dir, "ecozone_clean_prairies.shp")

# ----------------------------
# 3) Load and inspect data
# ----------------------------
map <- rast(map_file)
map
plot(map)

ecoreg <- vect(ecoreg_file)
ecoreg
plot(ecoreg)

pe <- vect(pe_file)
pe
plot(pe)

points <- st_read(file.path(plots_dir, "Prairie_Province_Veg_Plots.gdb"),
                 layer = "Prairie_Province_Veg_Plot_Locations_wgs84")

# ----------------------------
# 4) Configuration
# ----------------------------
# Reproject plots to match raster CRS
pts_proj <- st_transform(points, crs(map))

# Reproject ecoregions to match raster CRS
ecoreg_proj <- project(ecoreg, crs(map))

# Reproject prairie ecozone to match raster CRS
pe_proj <- project(pe, crs(map))

# Create bounding box for area of interest

# Buffer pe by 100 km (100,000 meters)
pe_buf <- buffer(pe_proj, width = 100000)

# Create bounding box (extent) around buffered pe
pe_bbox <- as.polygons(ext(pe_buf), crs = crs(pe_buf))

# Plot
plot(pe_bbox, border = "red", lwd = 2)
plot(pe_buf, add = TRUE, col = NA, border = "blue")
plot(pe_proj, add = TRUE, col = NA, border = "black")

#crs(map)
#crs(pe)
#crs(pe_bbox)

#ext(map)
#ext(pe_bbox)

plot(map)
plot(pe_bbox, add = TRUE, border = "red", lwd = 2)

# Clip and mask to area of interest
map_crop <- crop(map, pe_bbox)
map_mask <- mask(map_crop, pe_bbox)
plot(map_mask)

# Fill in group code for NGP solonetzic prairie

# 1) Extract the raster attribute (category) table
cats <- levels(map_mask)[[1]]

# 2) Identify rows where Value == 1215 and ELCODE is blank or NA
cats$ELCODE[cats$Value == 1215 & (is.na(cats$ELCODE) | cats$ELCODE == "")] <- "G1215"

# 3) Write the modified table back to the raster
levels(map_mask) <- cats

# 4) Verify the fix
levels(map_mask)[[1]][levels(map_mask)[[1]]$Value == 1215, ]

# how to handle rows without group code?

# ----------------------------
# 5) Extract raster group codes
# ----------------------------
# Convert plots to SpatVector
pts_vect <- vect(pts_proj)

# Extract ELCODE values at point locations
vals <- extract(map, pts_vect)

# Add ELCODE to sf object
pts_proj$ELCODE <- vals$ELCODE

# ----------------------------
# 5) Join raster group codes
# ----------------------------
 
cats <- levels(map_mask)[[1]]

names(cats)

pts_proj <- pts_proj |>
  left_join(
    cats[, c("ELCODE", "GroupKey")],
    by = "ELCODE"
  ) |>
  rename(map_group = GroupKey)

# How many points fell outside the raster?
sum(is.na(pts_proj$GroupKey))

# Inspect results
head(pts_proj[, c("ELCODE", "GroupKey")])

# ----------------------------
# 6) Filter out old plots
# ----------------------------

# how should I handle plots with NULL for year?

# Filter by sample_year >= 1990
pts_1990plus <- pts_proj |>
  dplyr::filter(sample_year >= 1990)

# Filter by sample_year >= 2000
pts_2000plus <- pts_proj |>
  dplyr::filter(sample_year >= 2000)

# ----------------------------
# 7) Binary majority-rule accuracy assessment
# ----------------------------

# Assumptions:
#  map_mask is a SpatRaster
#  Raster values are group codes (e.g., GroupKey or ELCODE already translated)
#  pts_proj is an sf object with a column:
#      Rpts_proj$map_group   # reference group for the pointShow more lines
#  CRS of points and raster already match
#  Raster is square (e.g., 30 m × 30 m)

# ------------------------------------------------------------
# 7.1) Convert points to SpatVector
# ------------------------------------------------------------
pts_vect <- vect(pts_proj)

# ------------------------------------------------------------
# 7.2) Get raster cell index for each point
# ------------------------------------------------------------
center_cells <- cellFromVect(map_mask, pts_vect)

# ------------------------------------------------------------
# 7.3) Get 3x3 neighborhood (center + 8 neighbors)
#    directions = 8 gives a Moore neighborhood
# ------------------------------------------------------------
nbr_cells <- adjacent(
  map_mask,
  cells = center_cells,
  directions = 8,
  include = TRUE   # include center pixel
)

# nbr_cells is a data.frame:
#   from   to
# where:
#   from = index of center cell
#   to   = index of neighboring (or same) cell

# ------------------------------------------------------------
# 7.4) Extract raster values for all neighborhood cells
# ------------------------------------------------------------
nbr_values <- extract(
  map_mask,
  nbr_cells$to
)[,2]   # second column is raster value

# Combine center cell index, raster value
nbr_df <- data.frame(
  center_cell = nbr_cells$from,
  raster_group = nbr_values
)

# ------------------------------------------------------------
# 7.5) Add reference group from points
# ------------------------------------------------------------
pts_lookup <- data.frame(
  center_cell = center_cells,
  point_group = as.character(pts_proj$map_group)
)

nbr_df <- nbr_df |>
  left_join(pts_lookup, by = "center_cell")

# ------------------------------------------------------------
# 7.6) Compute majority match per point
# ------------------------------------------------------------
majority_result <- nbr_df |>
  group_by(center_cell) |>
  summarise(
    n_total = sum(!is.na(raster_group)),
    n_match = sum(raster_group == point_group, na.rm = TRUE),
    majority_match = n_match >= (n_total / 2)
  )

# ------------------------------------------------------------
# 7.7) Attach result back to points
# ------------------------------------------------------------
pts_proj <- pts_proj |>
  mutate(center_cell = center_cells) |>
  left_join(
    majority_result |> select(center_cell, majority_match),
    by = "center_cell"
  ) |>
  select(-center_cell)

# output is binary accuracy field in points table:
# TRUE  = Majority of the 3×3 pixels match the point’s group
# FALSE = Majority of neighbors are different
pts_proj$majority_match

######################################################
######################################################

############### DISTRIBUTION MAPS ####################

# ============================================================
# Multi-panel map:
# - One panel per target group
# - Featured group highlighted in color
# - All other groups (including other target groups) shown in gray
# - Points overlaid and panels labeled using group_title
# ============================================================

# ------------------------------------------------------------
# 1) Define target groups (GroupKey values you want panels for)
# ------------------------------------------------------------

target_groups <- c(
  "G075",  # Northern tallgrass prairie
  "G141",  # NGP mesic mixedgrass prairie
  "G146",  # Northern tallgrass oak-aspen woodland
  "G273",  # Central rocky mountain lower montane, foothill & valley grassland
  "G328",  # Northwestern great plains aspen woodland
  "G331",  # NGP dry mixedgrass prairie
  "G332",  # NGP rough fescue prairie
  "G679",  # NGP ruderal exotic grassland & shrubland
  "G889",  # NGP sand prairie
  "G1215"  # NGP solonetzic prairie
)

# ------------------------------------------------------------
# 2) Prepare raster for plotting
#    NOTE:
#    - The raster only contains ELCODE values
#    - We do NOT filter the raster
#    - We create one logical copy per panel using expand_grid()
# ------------------------------------------------------------

r_df <- as.data.frame(map, xy = TRUE, na.rm = TRUE)

r_plot <- r_df |>
  # Create a panel for each target group
  expand_grid(panel_group = target_groups) |>
  # Highlight only the featured group in each panel; everything else = "Other"
  mutate(
    fill_class = ifelse(ELCODE == panel_group, panel_group, "Other")
  )

# ------------------------------------------------------------
# 3) Prepare point data
#    - pts_proj already contains map_group (= GroupKey)
#    - We align panel_group with map_group
# ------------------------------------------------------------

p_plot <- as.data.frame(pts_proj, geom = "xy") |>
  mutate(
    map_group   = as.character(map_group),
    panel_group = map_group
  ) |>
  # Only keep points belonging to groups that have panels
  filter(panel_group %in% target_groups)

# ------------------------------------------------------------
# 4) Build facet labels using group_title from points
# ------------------------------------------------------------

group_labels <- p_plot |>
  distinct(panel_group, group_title) |>
  deframe()   # named vector: names = GroupKey, values = group_title

# ------------------------------------------------------------
# 5) Define colors
#    - Each target group gets a distinct color
#    - All non-featured groups are gray
# ------------------------------------------------------------

group_colors <- setNames(
  hue_pal()(length(target_groups)),
  target_groups
)

fill_colors <- c(group_colors, Other = "grey80")

# ------------------------------------------------------------
# 6) Create the multi-panel figure
# ------------------------------------------------------------

ggplot() +
  # Raster layer: highlighted group in color, all others gray
  geom_raster(
    data = r_plot,
    aes(x = x, y = y, fill = fill_class),
    show.legend = FALSE
  ) +
  # Overlay points
  geom_point(
    data = p_plot,
    aes(x = x, y = y),
    size = 0.4,
    colour = "black"
  ) +
  # One panel per target group, labeled using group_title
  facet_wrap(
    ~ panel_group,
    labeller = labeller(panel_group = group_labels)
  ) +
  # Apply colors
  scale_fill_manual(values = fill_colors) +
  # Preserve map aspect ratio
  coord_equal() +
  # Clean cartographic theme
  theme_minimal() +
  theme(
    axis.text   = element_blank(),
    axis.title  = element_blank(),
    panel.grid  = element_blank(),
    strip.text  = element_text(size = 10)
  ) +
  labs(
    title    = "Spatial distribution of target CNVC groups",
    subtitle = "Featured group highlighted; all other groups shown in gray"
  )