
###################################################################################
##                       PRELIMINARY ACCURACY ASSESSMENT                         ##
###################################################################################

##----------------------- Alison Long - February 2026 -------------------------

# step 6: assess accuracy of datasets mapping grasslands by comparing them
#         to ground-truthed point data

# load packages
library(terra)
library(sf)
library(tidyverse)
library(tools)   # for file_path_sans_ext
library(ggplot2)
library(viridisLite)
library(ggspatial)

# Set the temporary directory to a drive with sufficient space
terraOptions(tempdir = "R:/CanadianGrasslands/WG_classification/scratch")

####--- DIRECTORIES ---####

# pilot area
pa_dir <- "R:/CanadianGrasslands/processeddatasets/pilotarea"

# study area
sa_dir <- "R:/CanadianGrasslands/processeddatasets/studyarea"

# grasslands
raw_grass_data_dir <- "R:/CanadianGrasslands/existingdatasets"
proc_grass_data_dir <- "R:/CanadianGrasslands/WG_Classification/processed"
out_dir <- "R:/CanadianGrasslands/WG_Classification/results/accuracy_assessment"

# ground-truthed data
gt_dir <- "R:/CanadianGrasslands/processeddatasets/groundtruth_data"

# ecodistricts 
eco_dir <- "R:/CanadianGrasslands/processeddatasets/boundaries"

# grass raster directories
aci_aligned_dir <- file.path(proc_grass_data_dir, "aafc_aci/aligned")
luts_aligned_dir <- file.path(proc_grass_data_dir, "aafc_luts_2020_30m_v5/aligned")
cec_aligned_dir <- file.path(proc_grass_data_dir, "cec_nalcms/aligned")
ngi_aligned_dir <- file.path(proc_grass_data_dir, "cfga_ngi/aligned")
cgr_aligned_dir <- file.path(proc_grass_data_dir, "cgr_gam_adv/aligned")
eccc_aligned_dir <- file.path(proc_grass_data_dir, "eccc_grass/aligned")
esa_aligned_dir <- file.path(proc_grass_data_dir, "esa_worldcover/aligned")
glc_aligned_dir <- file.path(proc_grass_data_dir, "glc_fcs30d/aligned")
cnvc_aligned_dir <- file.path(proc_grass_data_dir, "natureserve_cnvc/aligned")
lcc_aligned_dir <- file.path(proc_grass_data_dir, "nrcan_lcc/aligned")
pudl_aligned_dir <- file.path(proc_grass_data_dir, "phjv_pudl/aligned")
pli_aligned_dir <- file.path(proc_grass_data_dir, "sk_pli/aligned")
lcr_aligned_dir <- file.path(proc_grass_data_dir, "statcan_lcr/aligned")
gpw_aligned_dir <- file.path(proc_grass_data_dir, "wri_gpw_ann/aligned")
plow_aligned_dir <- file.path(proc_grass_data_dir, "wwf_plowprint/aligned")

####--- FILES ---####
# Pilot area
pa_file <- file.path(pa_dir, "pilot_area_aligned_cnvc_Jan2026.tif")

# Study area
sa_file <- file.path(sa_dir, "class_study_area_aligned_cnvc_Dec2025.tif")

# ground-truthed data
gt_file <- file.path(gt_dir, "gt_points_mod_w_ecodist.shp")

# ecodistricts
eco_file <- file.path(eco_dir, "ecodistricts_clean.shp")

# grass raster files
aci_file <- file.path(aci_aligned_dir, "aafc_aci_merged_reproj_algn_binary_noshrub.tif")
luts_file <- file.path(luts_aligned_dir, "aafc_luts_merged_reproj_algn_binary.tif")
cec_file <- file.path(cec_aligned_dir, "nalcms_reproj_algn_binary_noshrub.tif")
ngi_file <- file.path(ngi_aligned_dir, "cfga_ngi_merged_reproj_algn_binary.tif")
cgr_file <- file.path(cgr_aligned_dir, "cgr_reproj_algn_binary.tif")
eccc_file <- file.path(eccc_aligned_dir, "eccc_reclass_reproj_algn_binary.tif")
esa_file <- file.path(esa_aligned_dir, "esa_reproj_algn_binary_noshrub.tif")
glc_file <- file.path(glc_aligned_dir, "glc_reproj_algn_binary_noshrub.tif")
cnvc_file <- file.path(cnvc_aligned_dir, "cnvc_algn_binary.tif")
lcc_file <- file.path(lcc_aligned_dir, "lcc_reproj_algn_binary_noshrub.tif")
pudl_file <- file.path(pudl_aligned_dir, "pudl_reclass_reproj_algn_binary_noshrub.tif")
pli_file <- file.path(pli_aligned_dir, "sk_pli_merged_reproj_algn_binary_noshrub.tif")
lcr_file <- file.path(lcr_aligned_dir, "lcr_reproj_algn_binary.tif")
gpw_file <- file.path(gpw_aligned_dir, "gpw_reclass_reproj_algn_binary_noshrub.tif")
plow_file <- file.path(plow_aligned_dir, "plow_reproj_algn_binary.tif")

####--- READ INTO TERRA ---####

gt <- vect(gt_file)
gt

aci <- rast(aci_file)
luts <- rast(luts_file)
cec <- rast(cec_file)
ngi <- rast(ngi_file)
cgr <- rast(cgr_file)
eccc <- rast(eccc_file)
esa <- rast(esa_file)
glc <- rast(glc_file)
cnvc <- rast(cnvc_file)
lcc <- rast(lcc_file)
pudl <- rast(pudl_file)
pli <- rast(pli_file)
lcr <- rast(lcr_file)
gpw <- rast(gpw_file)
plow <- rast(plow_file)

files <- c(
  aci_file, luts_file, cec_file, ngi_file, cgr_file, eccc_file,
  esa_file, glc_file, cnvc_file, lcc_file, pudl_file, pli_file,
  lcr_file, gpw_file, plow_file
)

keys <- c(
  "aci","luts","cec","ngi","cgr","eccc", "esa","glc",
  "cnvc","lcc","pudl","pli","lcr","gpw","plow"
)

names(files) <- keys
rlist <- lapply(files, rast)  # each is single-layer SpatRaster

####--- PREP DATA ---####

# Inspect key grid properties across rasters
crs_all <- sapply(rlist, \(r) as.character(crs(r)))
res_all <- t(sapply(rlist, res))
org_all <- t(sapply(rlist, origin))
ext_all <- t(sapply(rlist, \(r) as.vector(ext(r))))
dim_all <- t(sapply(rlist, \(r) c(ncol(r), nrow(r))))

table(crs_all)          # should be one unique CRS
print(res_all)           # should be identical across rows
print(org_all)           # should be identical
print(ext_all)           # tiny differences here cause the error
print(dim_all)           # ncol/nrow must match to c()

# ---- Combine into a single multi-layer SpatRaster ----
# Because all rasters are aligned, this is safe and fast
r_all <- do.call(c, unname(rlist))  # c() concatenates layers

# Set clean layer names from the keys in order
names(r_all) <- keys

# ---- Ensure points are terra SpatVector ----
gt_vec <- if (inherits(gt, "sf")) terra::vect(gt) else gt
stopifnot(inherits(gt_vec, "SpatVector"))
stopifnot(geomtype(gt_vec) == "points")

# ---- Reproject points once if needed ----
gt_in_r_crs <- if (!identical(crs(gt_vec), crs(r_all))) {
  project(gt_vec, crs(r_all))
} else {
  gt_vec
}

####--- EXTRACT AND BIND DATA ---####

# ---- Single extract for all layers ----
vals_df <- terra::extract(r_all, gt_in_r_crs, ID = TRUE)

nrow(gt_in_r_crs)      # should be 196931
nrow(vals_df)         # should be 196931

# Basic counts
n_pt   <- nrow(gt_in_r_crs)
n_ext  <- nrow(vals_df)
n_pt; n_ext

# ID problems (should be 0 in each ideally)
dups   <- vals_df$ID[duplicated(vals_df$ID)]
miss   <- setdiff(seq_len(n_pt), vals_df$ID)        # IDs present in points but absent in extract
extra  <- setdiff(vals_df$ID, seq_len(n_pt))        # IDs not in 1..nrow(points)

length(dups); length(miss); length(extra)


# ---- Bind values back to your points ----
# vals_df: result of terra::extract(r_all, gt_in_r_crs, ID = TRUE)
# gt_in_r_crs: SpatVector points used in the extract()
# r_all: your multi-layer SpatRaster; names(r_all) are the desired column names

# 1) Make a full, ordered value table (one row per point, in 1..n order)
n_pt     <- nrow(gt_in_r_crs)
val_cols <- setdiff(names(vals_df), "ID")

row_map   <- match(seq_len(n_pt), vals_df$ID)     # which extract row corresponds to each point row
vals_only <- as.data.frame(
  matrix(NA, nrow = n_pt, ncol = length(val_cols)),
  stringsAsFactors = FALSE
)
names(vals_only) <- val_cols
good <- which(!is.na(row_map))
vals_only[good, ] <- vals_df[row_map[good], val_cols, drop = FALSE]

stopifnot(nrow(vals_only) == n_pt)

# 2) Assign directly into the attribute table of the SAME object (gt_in_r_crs)
gt_out <- gt_in_r_crs                           # keep geometry and row order
terra::values(gt_out)[, names(r_all)] <- vals_only

# How many geometries vs. how many rows in the attribute table?
n_geom <- nrow(gt_out)
n_attr <- nrow(terra::values(gt_out))
n_geom; n_attr

# Find attribute columns whose length != number of geometries
attr_len <- sapply(terra::values(gt_out), NROW)
bad_cols <- names(attr_len)[attr_len != n_geom]
good_cols <- names(attr_len)[attr_len == n_geom]

length(bad_cols); bad_cols[1:min(10, length(bad_cols))]  

# 1) Attributes (already in gt_out after your binding step)
vals <- terra::values(gt_out)               # data.frame with 196,931 rows

# 2) Coordinates for a SpatVector (points) as a data.frame
g   <- terra::geom(gt_out, df = TRUE)       # has columns 'geom','part','x','y','hole'
xy  <- g[, c("x", "y")]

# 3) Sanity check and combine
stopifnot(nrow(vals) == nrow(xy))           # should both be 196,931 for points
attr_df <- cbind(vals, xy, stringsAsFactors = FALSE)

# rename coord columns
names(attr_df)[(ncol(attr_df)-1):ncol(attr_df)] <- c("x", "y")

stopifnot(nrow(attr_df) == nrow(gt_out))
# Inspect NA shares per layer (helpful for coverage/masking)
sort(sapply(attr_df[names(r_all)], function(x) mean(is.na(x))), decreasing = TRUE)[1:10]

# add a field that categories the ground-truthed data into grass vs not (only "not" is shrub)

# Normalize and derive the new field
norm <- function(x) tolower(trimws(as.character(x)))
vals <- norm(attr_df$LC_Cl_N)

attr_df$LC_grass <- NA_integer_
attr_df$LC_grass[vals %in% c("mixed","native","tame")] <- 1L
attr_df$LC_grass[vals == "shrubs"] <- 0L

# Quick QA
table(attr_df$LC_grass, useNA = "ifany")


####--- EXPORT ---####

# ---- Write result ----

sf_use_s2(FALSE)  # for safety with non-geodesic ops (not strictly needed for write)

gt_sf <- st_as_sf(gt_out)
st_write(
  gt_sf,
  dsn    = "groundtruth_points_with_raster_values.gpkg",
  layer  = "groundtruth_points",        # explicit layer name avoids clashes
  delete_dsn = TRUE                     # overwrite dataset
)

# modified version to avoid session abort

gt_sf <- st_as_sf(gt_out)               # SpatVector -> sf
gt_sf <- st_zm(gt_sf, drop = TRUE)      # drop Z/M dimensions if any
gt_sf <- gt_sf[!st_is_empty(gt_sf), ]   # drop empty geometries (if any)

# If anything ended up not being POINT, cast it (should be no-op for points)
if (!all(st_geometry_type(gt_sf) == "POINT")) {
  gt_sf <- st_cast(gt_sf, "POINT", warn = FALSE)
}

st_write(
  gt_sf,
  dsn    = "groundtruth_points_with_raster_values2.gpkg",
  layer  = "groundtruth_points",
  delete_dsn       = TRUE,                     # overwrite file
  layer_options    = "SPATIAL_INDEX=NO",       # <-- avoid R-Tree build path
  dataset_options  = "VERSION=1.2",            # optional, avoids 1.4 trigger quirks
  quiet  = FALSE
)

# Create a spatial index on layer 'groundtruth_points' and geometry column 'geom'
# (If your geometry column has a different name, change 'geom' accordingly.)
cmd <- sprintf(
  "SELECT CreateSpatialIndex('%s','%s')",
  "groundtruth_points", "geom"
)
sf::gdal_utils(
  util = "vectortranslate",       # or "ogrinfo" with -sql (both work)
  source = "groundtruth_points_with_raster_values2.gpkg",
  options = c("-sql", cmd, "-dialect", "SQLITE")
)

# write to csv
csv_out_path <- file.path(out_dir, "groundtruth_points_with_grass_raster_values_ecodist.csv")
write.csv(attr_df, csv_out_path, row.names = FALSE)

## transform to EPSG:4326 so ArcGIS Pro can read gpkg

library(sf)

dsn   <- "groundtruth_points_with_raster_values3.gpkg"
layer <- "groundtruth_points"

# 1) Start from current sf (points) and transform to a standard EPSG
gt_ok <- gt_sf |>
  st_zm(drop = TRUE) |>        # drop Z/M if present
  (\(x) x[!st_is_empty(x), ])() |>
  st_transform(4326)           # or st_transform(3978)

# 2) Write GPKG, turning OFF spatial index to avoid the R-Tree fast-path
st_write(
  gt_ok, dsn = dsn, layer = layer,
  delete_dsn      = TRUE,
  layer_options   = "SPATIAL_INDEX=NO",
  dataset_options = "VERSION=1.2",  # conservative GPKG version for compatibility
  quiet = FALSE
)

# add the spatial index after the write

# Find the geometry column name (often "geom")
gcol <- st_geometry_name(gt_ok)

# Build the spatial index using the stable helper
cmd <- sprintf("SELECT CreateSpatialIndex('%s','%s')", layer, gcol)
sf::gdal_utils(util="ogrinfo", source=dsn, options=c("-sql", cmd, "-dialect", "SQLITE"))

#########################################################################################

# count number of points per ecoregion

# --- Read ecodistricts as SpatVector ---
eco <- vect(eco_file)  

# Inspect fields to choose a key for grouping
names(eco)

# >>> Set this to the correct field in your ecodistrict layer <<<
ECO_KEY <- "ECODISTRIC"   

# 1) Use a template raster that matches analysis grid
tmpl <- cnvc  # or luts, aci, etc.

# --- 2) Reproject ecodistricts to the template CRS ---
if (!identical(crs(eco), crs(tmpl))) {
  message("Projecting ecodistricts to the template raster CRS...")
  eco_proj <- terra::project(eco, tmpl)  # vector -> raster’s CRS
} else {
  eco_proj <- eco
}

# --- 3) Quick extent overlap check ---
ov <- try(terra::intersect(ext(eco_proj), ext(tmpl)), silent = TRUE)
if (inherits(ov, "try-error") || is.null(ov)) {
  stop("The extents of 'eco_proj' and 'tmpl' do not overlap. ",
       "Verify you selected the right template and that CRS/areas match.")
}

# --- 4) Build an integer field for rasterization + lookup (LUT) ---
eco_df <- as.data.frame(eco_proj)
eco_key_raw <- eco_df[[ECO_KEY]]

if (is.numeric(eco_key_raw)) {
  eco_proj$eco_id_int <- as.integer(eco_key_raw)
  lut <- tibble(eco_id_int = eco_proj$eco_id_int,
                !!ECO_KEY := as.character(eco_key_raw)) |>
    distinct(eco_id_int, .keep_all = TRUE)
} else {
  eco_proj$eco_id_int <- as.integer(factor(trimws(as.character(eco_key_raw))))
  lut <- tibble(eco_id_int = eco_proj$eco_id_int,
                !!ECO_KEY := trimws(as.character(eco_key_raw))) |>
    distinct(eco_id_int, .keep_all = TRUE)
}
id_field <- "eco_id_int"

# --- 5) Rasterize ecodistricts onto the template grid ---
# 'touches=TRUE' assigns a cell if polygon touches it (good for full coverage)
eco_r <- terra::rasterize(eco_proj, tmpl, field = id_field, touches = TRUE, background = NA)
names(eco_r) <- id_field

# --- 6) Ensure points are also in the template CRS ---
gt_in_tmpl_crs <- if (!identical(crs(gt), crs(tmpl))) project(gt, crs(tmpl)) else gt

# --- 7) Extract ecodistrict id per point (raster -> points) ---
pt_vals <- terra::extract(eco_r, gt_in_tmpl_crs, ID = FALSE)  # returns data.frame with 'eco_id_int'
pt_keys <- tibble(eco_id_int = pt_vals[[id_field]]) |>
  left_join(lut, by = "eco_id_int") |>
  transmute(!!ECO_KEY := .data[[ECO_KEY]])

# --- 8) Count points per ecodistrict, include zeros ---
eco_counts <- pt_keys |>
  filter(!is.na(.data[[ECO_KEY]])) |>
  count(.data[[ECO_KEY]], name = "n_points") |>
  arrange(desc(n_points))

all_keys <- lut |>
  select(all_of(ECO_KEY)) |>
  distinct()

eco_counts_full <- all_keys |>
  left_join(eco_counts, by = ECO_KEY) |>
  mutate(n_points = replace_na(n_points, 0L))

# --- 9) Export CSV ---
eco_counts_csv <- file.path(out_dir, "ecodistrict_point_counts.csv")
write.csv(eco_counts_full, eco_counts_csv, row.names = FALSE)

# --- 10) (Optional) Attach counts to polygons and write GPKG for mapping ---
eco_with_counts <- eco_proj
match_idx <- match(trimws(as.character(as.data.frame(eco_proj)[[ECO_KEY]])),
                   eco_counts_full[[ECO_KEY]])
eco_with_counts$n_points <- eco_counts_full$n_points[match_idx]

gpkg_out <- file.path(out_dir, "ecodistrict_point_counts.gpkg")
eco_sf_out <- st_as_sf(eco_with_counts) |>
  st_zm(drop = TRUE)

st_write(eco_sf_out, dsn = gpkg_out, layer = "ecodistricts_with_counts",
         delete_dsn = TRUE, layer_options = "SPATIAL_INDEX=NO", quiet = FALSE)

# Optional: build spatial index after write
gcol <- st_geometry_name(eco_sf_out)
cmd  <- sprintf("SELECT CreateSpatialIndex('%s','%s')", "ecodistricts_with_counts", gcol)
sf::gdal_utils(util = "ogrinfo", source = gpkg_out, options = c("-sql", cmd, "-dialect", "SQLITE"))

# --- 11) Sanity checks ---
n_total  <- nrow(gt_in_tmpl_crs)
n_inside <- sum(!is.na(pt_keys[[ECO_KEY]]))
message(sprintf("Points inside ecodistricts: %s / %s (%.1f%%)",
                format(n_inside, big.mark=","), format(n_total, big.mark=","),
                100 * n_inside / n_total))

# Top few counts for a quick look
print(head(eco_counts, 10))

######## plot

# ---- Read data ----
gt_sf  <- st_read(gt_file, quiet = TRUE)
eco_sf <- st_read(eco_file, quiet = TRUE)

# ---- CRS handling ----
if (st_crs(gt_sf) != st_crs(eco_sf)) {
  gt_sf <- st_transform(gt_sf, st_crs(eco_sf))
}

# ---- Optional: ensure validity ----
if (any(!st_is_valid(eco_sf))) {
  eco_sf <- st_make_valid(eco_sf)
}

# ---- Normalize key ----
stopifnot(ECO_KEY %in% names(eco_sf))
eco_sf[[ECO_KEY]] <- trimws(as.character(eco_sf[[ECO_KEY]]))

# ---- Join points -> ecodistrict (1:1, deterministic) ----
sf_use_s2(FALSE)  # planar ops for projected data
gt_joined <- st_join(gt_sf, eco_sf[, ECO_KEY, drop = FALSE], left = TRUE, largest = TRUE)

# ---- Count per ecodistrict ----
eco_counts <- gt_joined |>
  st_drop_geometry() |>
  filter(!is.na(.data[[ECO_KEY]])) |>
  count(.data[[ECO_KEY]], name = "n_points") |>
  arrange(desc(n_points))

# ---- Attach counts back to polygons ----
eco_map <- eco_sf |>
  left_join(eco_counts, by = ECO_KEY)

# ---- EXCLUDE ecodistricts with zero (or missing) points ----
eco_map_nonzero <- eco_map |>
  filter(!is.na(n_points) & n_points > 0)

# ---- Sanity message ----
n_total  <- nrow(gt_sf)
n_inside <- sum(!is.na(gt_joined[[ECO_KEY]]))
n_mapped <- nrow(eco_map_nonzero)
message(sprintf("Points assigned: %s / %s (%.1f%%) • Ecodistricts mapped (n_points > 0): %s",
                format(n_inside, big.mark=","), format(n_total, big.mark=","),
                100 * n_inside / n_total, n_mapped))

# ---- Plot ----
title_txt <- "Ground-Truth Points per Ecodistrict"
subtitle_txt <- sprintf("Only ecodistricts with > 0 points (Assigned: %s of %s)",
                        format(n_inside, big.mark=","), format(n_total, big.mark=","))

p <- ggplot() +
  geom_sf(data = eco_map_nonzero, aes(fill = n_points), color = "grey30", size = 0.1) +
  scale_fill_viridis_c(
    option = "C", direction = 1,
    trans = "sqrt",   # use "identity" if you want raw (untransformed) scale
    name = "Points"
  ) +
  coord_sf(datum = NA) +
  labs(title = title_txt, subtitle = subtitle_txt, caption = "Ecodistricts with > 0 points") +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "right",
    panel.grid.major = element_line(color = "grey85", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    axis.title = element_blank()
  ) #+
  #annotation_scale(location = "bl", width_hint = 0.2, line_width = 0.3, text_cex = 0.7) +
  #annotation_north_arrow(location = "bl", which_north = "true",
  #                       style = north_arrow_fancy_orienteering)

print(p)

# ---- Save PNG ----
png_out <- file.path(out_dir, "map_ecodistricts_point_counts_nonzero.png")
ggsave(filename = png_out, plot = p, width = 10, height = 7, dpi = 300)
message("Map saved: ", png_out)

# ---- (Optional) Save only the non-zero polygons to GPKG ----
gpkg_out <- file.path(out_dir, "ecodistricts_with_counts_nonzero.gpkg")
st_write(eco_map_nonzero, dsn = gpkg_out, layer = "ecodistricts_with_counts_nonzero",
         delete_dsn = TRUE, layer_options = "SPATIAL_INDEX=NO", quiet = FALSE)

# Optional: build spatial index after write
gcol <- st_geometry_name(eco_map_nonzero)
cmd  <- sprintf("SELECT CreateSpatialIndex('%s','%s')", "ecodistricts_with_counts_nonzero", gcol)
sf::gdal_utils(util = "ogrinfo", source = gpkg_out, options = c("-sql", cmd, "-dialect", "SQLITE"))

## symbolize ecodistricts with 5 or less points differently

# Split into "low" (<=5) and "high" (>5) groups
eco_low  <- eco_map_nonzero |> dplyr::filter(n_points <= 5)
eco_high <- eco_map_nonzero |> dplyr::filter(n_points > 5)

title_txt <- "Ground-Truth Points per Ecodistrict"
subtitle_txt <- "Green = 5 or less points; others colored by count"

p <- ggplot() +
  # Draw >5 first with a continuous scale
  geom_sf(data = eco_high, aes(fill = n_points), color = "grey30", size = 0.1) +
  scale_fill_viridis_c(
    option = "C", direction = 1,
    trans = "sqrt",        # change to "identity" if you want raw counts
    name = "Points (>5)"
  ) +
  # Overlay ≤5 in green (kept out of the continuous legend)
  geom_sf(data = eco_low, fill = "#2ca25f", color = "grey30", size = 0.1, inherit.aes = FALSE) +
  coord_sf(datum = NA) +
  labs(title = title_txt, subtitle = subtitle_txt) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "right",
    panel.grid.major = element_line(color = "grey85", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    axis.title = element_blank()
  ) #+
  #annotation_scale(location = "bl", width_hint = 0.2, line_width = 0.3, text_cex = 0.7) +
  #annotation_north_arrow(location = "bl", which_north = "true",
  #                       style = north_arrow_fancy_orienteering)

print(p)

# Save
png_out <- file.path(out_dir, "map_ecodistricts_point_counts_nonzero_green_le5.png")
ggsave(filename = png_out, plot = p, width = 10, height = 7, dpi = 300)
message("Map saved: ", png_out)

##########################################################################################

# version that retains only groundtruthed data for past 10 years
attr_df_2015 <- subset(attr_df, Yer_Cll >= 2015)

# Confirm the filter worked
summary(attr_df_2015$Yer_Cll)

# Count how many rows remained
nrow(attr_df_2015)

# also filter spatial object
keep_rows <- which(attr_df$Yer_Cll >= 2015)
gt_sf_2015 <- gt_sf[keep_rows, ]

# =========================
# Use ONLY ground-truth points from 2015+
# =========================

# ---- Assumptions / Objects from above ----
# - gt: SpatVector of ALL ground-truth points (same order as attr_df)
# - attr_df: attributes for gt (same row order)
# - attr_df_2015: subset(attr_df, Yer_Cll >= 2015)
# - eco_file, cnvc (template raster), out_dir already defined
# - ECO_KEY is the ecodistrict key ("ECODISTRIC")

#ECO_KEY <- "ECODISTRIC"

# ---- Build row index for points since 2015 ----
stopifnot(nrow(attr_df) == nrow(gt))
idx_2015 <- which(attr_df$Yer_Cll >= 2015)
if (length(idx_2015) == 0L) {
  stop("No points with Yer_Cll >= 2015 found in attr_df.")
}
message("Points since 2015: ", length(idx_2015), " of ", nrow(gt))

# ---- Read ecodistricts as SpatVector ----
eco <- vect(eco_file)

# ---- 1) Template raster (aligned grid) ----
tmpl <- cnvc  # or luts/aci etc.

# ---- 2) Reproject ecodistricts to template CRS (if needed) ----
if (!identical(crs(eco), crs(tmpl))) {
  message("Projecting ecodistricts to the template raster CRS...")
  eco_proj <- terra::project(eco, tmpl)
} else {
  eco_proj <- eco
}

# ---- 3) Quick extent overlap check ----
ov <- try(terra::intersect(ext(eco_proj), ext(tmpl)), silent = TRUE)
if (inherits(ov, "try-error") || is.null(ov)) {
  stop("The extents of 'eco_proj' and 'tmpl' do not overlap. ",
       "Verify the template raster and ecodistrict AOI.")
}

# ---- 4) Build integer field for rasterization + LUT ----
eco_df <- as.data.frame(eco_proj)
eco_key_raw <- eco_df[[ECO_KEY]]

if (is.numeric(eco_key_raw)) {
  eco_proj$eco_id_int <- as.integer(eco_key_raw)
  lut <- tibble(eco_id_int = eco_proj$eco_id_int,
                !!ECO_KEY := as.character(eco_key_raw)) |>
    distinct(eco_id_int, .keep_all = TRUE)
} else {
  eco_proj$eco_id_int <- as.integer(factor(trimws(as.character(eco_key_raw))))
  lut <- tibble(eco_id_int = eco_proj$eco_id_int,
                !!ECO_KEY := trimws(as.character(eco_key_raw))) |>
    distinct(eco_id_int, .keep_all = TRUE)
}
id_field <- "eco_id_int"

# ---- 5) Rasterize ecodistricts to template grid ----
eco_r <- terra::rasterize(eco_proj, tmpl, field = id_field, touches = TRUE, background = NA)
names(eco_r) <- id_field

# ---- 6) Ensure points are in template CRS, THEN subset to 2015+ ----
gt_in_tmpl_crs <- if (!identical(crs(gt), crs(tmpl))) project(gt, crs(tmpl)) else gt

# Safety check row counts still align after projection
stopifnot(nrow(gt_in_tmpl_crs) == nrow(attr_df))

# Use only the points from 2015+
gt_in_tmpl_crs_2015 <- gt_in_tmpl_crs[idx_2015, ]

# ---- 7) Extract ecodistrict id per (filtered) point ----
pt_vals <- terra::extract(eco_r, gt_in_tmpl_crs_2015, ID = FALSE)  # data.frame with 'eco_id_int'

pt_keys <- tibble(eco_id_int = pt_vals[[id_field]]) |>
  left_join(lut, by = "eco_id_int") |>
  transmute(!!ECO_KEY := .data[[ECO_KEY]])

# ---- 8) Count points per ecodistrict (for filtered points), include zeros for reference table ----
eco_counts <- pt_keys |>
  filter(!is.na(.data[[ECO_KEY]])) |>
  count(.data[[ECO_KEY]], name = "n_points") |>
  arrange(desc(n_points))

all_keys <- lut |>
  select(all_of(ECO_KEY)) |>
  distinct()

eco_counts_full <- all_keys |>
  left_join(eco_counts, by = ECO_KEY) |>
  mutate(n_points = replace_na(n_points, 0L))

# ---- 9) Export CSV (tagged since2015) ----
eco_counts_csv <- file.path(out_dir, "ecodistrict_point_counts_since2015.csv")
write.csv(eco_counts_full, eco_counts_csv, row.names = FALSE)

# ---- 10) Attach counts to polygons & write GPKG (tagged since2015) ----
eco_with_counts <- eco_proj
match_idx <- match(trimws(as.character(as.data.frame(eco_proj)[[ECO_KEY]])),
                   eco_counts_full[[ECO_KEY]])
eco_with_counts$n_points <- eco_counts_full$n_points[match_idx]

gpkg_out <- file.path(out_dir, "ecodistrict_point_counts_since2015.gpkg")
eco_sf_out <- st_as_sf(eco_with_counts) |>
  st_zm(drop = TRUE)

st_write(eco_sf_out, dsn = gpkg_out, layer = "ecodistricts_with_counts_since2015",
         delete_dsn = TRUE, layer_options = "SPATIAL_INDEX=NO", quiet = FALSE)

# Optional spatial index
gcol <- st_geometry_name(eco_sf_out)
cmd  <- sprintf("SELECT CreateSpatialIndex('%s','%s')", "ecodistricts_with_counts_since2015", gcol)
sf::gdal_utils(util = "ogrinfo", source = gpkg_out, options = c("-sql", cmd, "-dialect", "SQLITE"))

# ---- 11) Sanity checks using filtered points ----
n_total_2015  <- nrow(gt_in_tmpl_crs_2015)
n_inside_2015 <- sum(!is.na(pt_keys[[ECO_KEY]]))
message(sprintf("[2015+] Points inside ecodistricts: %s / %s (%.1f%%)",
                format(n_inside_2015, big.mark=","), format(n_total_2015, big.mark=","),
                100 * n_inside_2015 / n_total_2015))
print(head(eco_counts, 10))


# =========================
# PLOTTING (filtered points only)
# =========================

# ---- Read with sf (polygons), build counts from filtered points again for mapping ----
eco_sf <- st_read(eco_file, quiet = TRUE)
if (any(!st_is_valid(eco_sf))) eco_sf <- st_make_valid(eco_sf)
stopifnot(ECO_KEY %in% names(eco_sf))
eco_sf[[ECO_KEY]] <- trimws(as.character(eco_sf[[ECO_KEY]]))

# Build an sf for the filtered points in the polygon CRS
gt_sf_all <- st_as_sf(gt)                # all points
if (st_crs(gt_sf_all) != st_crs(eco_sf)) {
  gt_sf_all <- st_transform(gt_sf_all, st_crs(eco_sf))
}

# Subset to the same 2015+ rows by index to ensure perfect alignment
gt_sf_2015 <- gt_sf_all[idx_2015, ]

# Deterministic 1:1 assignment for mapping (vector join on filtered points)
sf_use_s2(FALSE)
gt_joined_2015 <- st_join(gt_sf_2015, eco_sf[, ECO_KEY, drop = FALSE], left = TRUE, largest = TRUE)

eco_counts_2015 <- gt_joined_2015 |>
  st_drop_geometry() |>
  filter(!is.na(.data[[ECO_KEY]])) |>
  count(.data[[ECO_KEY]], name = "n_points") |>
  arrange(desc(n_points))

# Attach counts to polygons (note: for the MAP we will exclude zeros)
eco_map <- eco_sf |>
  left_join(eco_counts_2015, by = ECO_KEY)

# EXCLUDE ecodistricts with zero (or missing) points
eco_map_nonzero <- eco_map |>
  filter(!is.na(n_points) & n_points > 0)

# Split for ≤5 green vs >5 with continuous scale
eco_low  <- eco_map_nonzero |> filter(n_points <= 5)
eco_high <- eco_map_nonzero |> filter(n_points > 5)

# Plot
title_txt <- "Ground-Truth Points per Ecodistrict (2015+ Only)"
subtitle_txt <- "Green = 5 or fewer points; others colored by count"

p <- ggplot() +
  # draw >5 first with a continuous scale
  geom_sf(data = eco_high, aes(fill = n_points), color = "grey30", size = 0.1) +
  scale_fill_viridis_c(
    option = "C", direction = 1,
    trans = "sqrt",
    name  = "Points (>5)"
  ) +
  # overlay <=5 in green
  geom_sf(data = eco_low, fill = "#2ca25f", color = "grey30", size = 0.1, inherit.aes = FALSE) +
  coord_sf(datum = NA) +
  labs(title = title_txt, subtitle = subtitle_txt, caption = "Ecodistricts with > 0 points (2015+)") +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "right",
    panel.grid.major = element_line(color = "grey85", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    axis.title = element_blank()
  )

print(p)

# Save PNG (tagged since2015)
png_out <- file.path(out_dir, "map_ecodistricts_point_counts_nonzero_since2015_green_le5.png")
ggsave(filename = png_out, plot = p, width = 10, height = 7, dpi = 300)
message("Map saved: ", png_out)

# Optional: Save only non-zero polygons (2015+) to GPKG
gpkg_out2 <- file.path(out_dir, "ecodistricts_with_counts_nonzero_since2015.gpkg")
st_write(eco_map_nonzero, dsn = gpkg_out2, layer = "ecodistricts_with_counts_nonzero_since2015",
         delete_dsn = TRUE, layer_options = "SPATIAL_INDEX=NO", quiet = FALSE)

# Optional spatial index
gcol2 <- st_geometry_name(eco_map_nonzero)
cmd2  <- sprintf("SELECT CreateSpatialIndex('%s','%s')", "ecodistricts_with_counts_nonzero_since2015", gcol2)
sf::gdal_utils(util = "ogrinfo", source = gpkg_out2, options = c("-sql", cmd2, "-dialect", "SQLITE"))
