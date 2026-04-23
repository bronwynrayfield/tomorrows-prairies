
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
out_dir <- "R:/CanadianGrasslands/WG_Classification/processed/groundtruthed_data_points"

# ground-truthed data
gt_dir <- "R:/CanadianGrasslands/processeddatasets/groundtruth_data"

# ecodistricts 
eco_dir <- "R:/CanadianGrasslands/processeddatasets/boundaries/national-ecological-framework-for-canada/Ecodistricts"

# grass raster directories
aci_aligned_dir   <- file.path(proc_grass_data_dir, "aafc_aci/aligned")
luts30_masked_dir <- file.path(proc_grass_data_dir, "aafc_luts_2020_30m_v5/masked")
luts10_masked_dir <- file.path(proc_grass_data_dir, "aafc_luts_2020_10m_v1/masked")
cec_masked_dir    <- file.path(proc_grass_data_dir, "cec_nalcms/masked")
ngi_aligned_dir   <- file.path(proc_grass_data_dir, "cfga_ngi/aligned")
cgr_masked_dir    <- file.path(proc_grass_data_dir, "cgr_gam_adv/masked")
eccc_masked_dir   <- file.path(proc_grass_data_dir, "eccc_grass/masked")
esa_masked_dir    <- file.path(proc_grass_data_dir, "esa_worldcover/masked")
glc_masked_dir    <- file.path(proc_grass_data_dir, "glc_fcs30d/masked")
cnvc_masked_dir   <- file.path(proc_grass_data_dir, "natureserve_cnvc/masked")
lcc_masked_dir    <- file.path(proc_grass_data_dir, "nrcan_lcc/masked")
pudl_aligned_dir  <- file.path(proc_grass_data_dir, "phjv_pudl/aligned")
pli_aligned_dir   <- file.path(proc_grass_data_dir, "sk_pli/aligned")
lcr_masked_dir    <- file.path(proc_grass_data_dir, "statcan_lcr/masked")
gpw_masked_dir    <- file.path(proc_grass_data_dir, "wri_gpw_ann/masked")
plow_masked_dir   <- file.path(proc_grass_data_dir, "wwf_plowprint/masked")
plvi_aligned_dir  <- file.path(proc_grass_data_dir, "ab_plvi/aligned")
hfi_aligned_dir   <- file.path(proc_grass_data_dir, "abmi_hfi/aligned")
gcc_aligned_dir   <- file.path(proc_grass_data_dir, "bc_gcc/aligned")
dw_masked_dir     <- file.path(proc_grass_data_dir, "dw_2025_grass/masked")
esri_masked_dir   <- file.path(proc_grass_data_dir, "esri_rangelands/masked")
mgi_aligned_dir   <- file.path(proc_grass_data_dir, "mb_mgi/aligned")

####--- FILES ---####
# Pilot area
pa_file <- file.path(pa_dir, "pilot_area_aligned_cnvc_Jan2026.tif")

# Study area
sa_file <- file.path(sa_dir, "class_study_area_aligned_cnvc_Dec2025.tif")

# Modelling study area shapefile for plotting
sap_file <- file.path(sa_dir, "modelling_study_area_aligned_cnvc_Jan2026.shp")

# ground-truthed data
gt_file <- file.path(gt_dir, "gt_points_modified_Mar2026.shp")

# ecodistricts
eco_file <- file.path(eco_dir, "ecodistricts_clean.shp")

# grass raster files
aci_file    <- file.path(aci_aligned_dir, "aafc_aci_merged_reproj_algn_binary_noshrub.tif")
luts30_file <- file.path(luts30_masked_dir, "luts30_reproj_algn_binary_masked.tif")
luts10_file <- file.path(luts10_masked_dir, "luts10_reproj_algn_binary_masked.tif")
cec_file    <- file.path(cec_masked_dir, "nalcms_reproj_algn_binary_noshrub_masked.tif")
ngi_file    <- file.path(ngi_aligned_dir, "cfga_ngi_merged_reproj_algn_binary.tif")
cgr_file    <- file.path(cgr_masked_dir, "cgr_reproj_algn_binary_masked.tif")
eccc_file   <- file.path(eccc_masked_dir, "eccc_reclass_reproj_algn_binary_masked.tif")
esa_file    <- file.path(esa_masked_dir, "esa_reproj_algn_binary_noshrub_masked.tif")
glc_file    <- file.path(glc_masked_dir, "glc_reproj_algn_binary_noshrub_masked.tif")
cnvc_file   <- file.path(cnvc_masked_dir, "cnvc_algn_binary_masked.tif")
lcc_file    <- file.path(lcc_masked_dir, "lcc_reproj_algn_binary_noshrub_masked.tif")
pudl_file   <- file.path(pudl_aligned_dir, "pudl_reclass_reproj_algn_binary_noshrub.tif")
pli_ap_file <- file.path(pli_aligned_dir, "sk_pli_ap_reproj_algn_binary_noshrub.tif")
pli_cu_file <- file.path(pli_aligned_dir, "sk_pli_cu_reproj_algn_binary_noshrub.tif")
pli_mg_file <- file.path(pli_aligned_dir, "sk_pli_mg_reproj_algn_binary_noshrub.tif")
pli_mm_file <- file.path(pli_aligned_dir, "sk_pli_mm_reproj_algn_binary_noshrub.tif")
lcr_file    <- file.path(lcr_masked_dir, "lcr_reproj_algn_binary_masked.tif")
gpw_file    <- file.path(gpw_masked_dir, "gpw_reproj_algn_binary_noshrub_masked.tif")
plow_file   <- file.path(plow_masked_dir, "plow_reproj_algn_intact_masked.tif")
plvi_file   <- file.path(plvi_aligned_dir, "plvi_reproj_algn_binary_gte2015.tif")
hfi_file    <- file.path(hfi_aligned_dir, "hfi_reproj_algn_binary_gte2015.tif")
gcc_file    <- file.path(gcc_aligned_dir, "gcc_reproj_algn_binary.tif")
dw_file     <- file.path(dw_masked_dir, "dw_reproj_algn_masked_binary.tif")
esri_file   <- file.path(esri_masked_dir, "esri_reproj_algn_masked_binary.tif")
mgi_file    <- file.path(mgi_aligned_dir, "mb_mgi_reproj_algn_binary_noshrub.tif")

####--- READ INTO TERRA ---####

gt <- vect(gt_file)
gt

cnvc <- rast(cnvc_file) # for plotting
sap  <- vect(sap_file)  # for plotting

files <- c(
  aci_file, luts10_file, luts30_file, cec_file, ngi_file, cgr_file, eccc_file,
  esa_file, glc_file, cnvc_file, lcc_file, pudl_file, pli_ap_file,
  pli_cu_file, pli_mg_file, pli_mm_file, lcr_file, gpw_file, plow_file,
  plvi_file, hfi_file, gcc_file, dw_file, esri_file, mgi_file
)

keys <- c(
  "aci","luts10","luts30","cec","ngi","cgr","eccc", "esa","glc",
  "cnvc","lcc","pudl","pli_ap","pli_cu","pli_mg","pli_mm","lcr",
  "gpw","plow","plvi","hfi","gcc","dw","esri","mgi"
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

nrow(gt_in_r_crs)      # 200987
nrow(vals_df)          # 200987

# Basic counts
n_pt   <- nrow(gt_in_r_crs)
n_ext  <- nrow(vals_df)
n_pt; n_ext

# ID problems (should be 0 in each ideally)
dups   <- vals_df$ID[duplicated(vals_df$ID)]
miss   <- setdiff(seq_len(n_pt), vals_df$ID)        # IDs present in points but absent in extract
extra  <- setdiff(vals_df$ID, seq_len(n_pt))        # IDs not in 1..nrow(points)

length(dups); length(miss); length(extra)


# ---- Bind values back to points ----
# vals_df: result of terra::extract(r_all, gt_in_r_crs, ID = TRUE)
# gt_in_r_crs: SpatVector points used in the extract()
# r_all: multi-layer SpatRaster; names(r_all) are the desired column names

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

# 1) Attributes (already in gt_out after binding step)
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

# add a field that categorizes the ground-truthed data into grass vs not ("not" is shrub/woody and annual cropland)

# Normalize and derive the new field
norm <- function(x) tolower(trimws(as.character(x)))
vals <- norm(attr_df$LC_Cl_N)

attr_df$LC_grass <- NA_integer_
attr_df$LC_grass[vals %in% c("mixed","native","tame")] <- 1L
attr_df$LC_grass[vals == "woody"] <- 0L
attr_df$LC_grass[vals == "annual cropland"] <- 0L

# Quick QA
table(attr_df$LC_grass, useNA = "ifany")

head(attr_df)


####--- EXPORT ---####

# write to csv
csv_out_path <- file.path(out_dir, "groundtruth_points_with_grass_raster_values_binary_ecodist_pli_sep_Apr2026.csv")
write.csv(attr_df, csv_out_path, row.names = FALSE)

# filter attr_df to only include rows where Yer_Cll is 2015+
attr_df_2015 <- attr_df %>%
  filter(Yer_Cll >= 2015)

head(attr_df_2015)

csv_out_path <- file.path(out_dir, "groundtruth_points_with_grass_raster_values_binary_pli_sep_ge2015_Apr2026.csv")
write.csv(attr_df_2015, csv_out_path, row.names = FALSE)

####--- CLEAN FOR NEXT STEP ---####

# filter out any non-grass points (woody and annual cropland)
attr_df_2015 <- attr_df_2015 %>%
  filter(LC_grass == 1)

# add source column and fill when points were used as training data 
attr_df_2015 <- attr_df_2015 %>%
  mutate(
    Source = case_when(
      Sorc_ID == "PLI_SKMoE" ~ "pli",
      Sorc_ID == "ABMI_HFI" ~ "hfi",
      Sorc_ID == "CFGA_NGI" & Yer_Cll %in% c(2022, 2023) ~ "ngi",
      Sorc_ID == "ACI_field_2024" ~ "aci",
      Sorc_ID == "BC_ACI_2025" ~ "aci",
      TRUE ~ ""    # keep unmatched cases empty
    )
  )

# write cleaned dataset to csv
csv_out_path <- file.path(out_dir, "groundtruth_points_with_grass_raster_values_binary_pli_sep_ge2015_Apr2026_clean.csv")
write.csv(attr_df_2015, csv_out_path, row.names = FALSE)

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
eco_counts_csv <- file.path(out_dir, "ecodistrict_point_counts_Mar2026.csv")
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
sap_sf <- sf::st_as_sf(sap)

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
  geom_sf(
    data = sap_sf,
    fill  = NA,           # no fill
    color = "lightblue",     # outline color
    size  = 0.7,          # adjust thickness as needed
    linewidth = 0.7       # for newer ggplot versions
  ) +
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
png_out <- file.path(out_dir, "map_ecodistricts_point_counts_nonzero_Mar2026.png")
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
  geom_sf(
    data = sap_sf,
    fill  = NA,           # no fill
    color = "lightblue",     # outline color
    size  = 0.7,          # adjust thickness as needed
    linewidth = 0.7       # for newer ggplot versions
  ) +
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
png_out <- file.path(out_dir, "map_ecodistricts_point_counts_nonzero_green_le5_Mar2026.png")
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
eco_counts_csv <- file.path(out_dir, "ecodistrict_point_counts_since2015_Mar2026.csv")
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
  # Prairies Ecozone boundary
  geom_sf(
    data = sap_sf,
    fill  = NA,           # no fill
    color = "lightblue",     # outline color
    size  = 0.7,          # adjust thickness as needed
    linewidth = 0.7       # for newer ggplot versions
  ) +
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
png_out <- file.path(out_dir, "map_ecodistricts_point_counts_nonzero_since2015_green_le5_Mar2026.png")
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