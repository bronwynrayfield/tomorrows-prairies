
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
eco_dir <- "R:/CanadianGrasslands/processeddatasets/boundaries/national-ecological-framework-for-canada/Ecodistricts"

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

# grass raster files - NEED TO UPDATE WITH RASTERS THAT HAVE NATIVE (3), MIXED (2), AND TAME (1)!!
aci_file <- file.path(aci_aligned_dir, "aafc_aci_merged_reproj_algn_natvtam_noshrub.tif")
luts_file <- file.path(luts_aligned_dir, "aafc_luts_merged_reproj_algn_natvtam.tif")
cec_file <- file.path(cec_aligned_dir, "nalcms_reproj_algn_natvtam_noshrub.tif")
ngi_file <- file.path(ngi_aligned_dir, "cfga_ngi_merged_reproj_algn_natvtam.tif")
cgr_file <- file.path(cgr_aligned_dir, "cgr_reproj_algn_natvtam.tif")
eccc_file <- file.path(eccc_aligned_dir, "eccc_reclass_reproj_algn_natvtam.tif")
esa_file <- file.path(esa_aligned_dir, "esa_reproj_algn_natvtam_noshrub.tif")
glc_file <- file.path(glc_aligned_dir, "glc_reproj_algn_natvtam_noshrub.tif")
cnvc_file <- file.path(cnvc_aligned_dir, "cnvc_algn_natvtam.tif")
lcc_file <- file.path(lcc_aligned_dir, "lcc_reproj_algn_natvtam_noshrub.tif")
pudl_file <- file.path(pudl_aligned_dir, "pudl_reclass_reproj_algn_natvtam_noshrub.tif")
pli_ap_file <- file.path(pli_aligned_dir, "sk_pli_ap_reproj_algn_natvtam_noshrub.tif")
pli_cu_file <- file.path(pli_aligned_dir, "sk_pli_cu_reproj_algn_natvtam_noshrub.tif")
pli_mg_file <- file.path(pli_aligned_dir, "sk_pli_mg_reproj_algn_natvtam_noshrub.tif")
pli_mm_file <- file.path(pli_aligned_dir, "sk_pli_mm_reproj_algn_natvtam_noshrub.tif")
lcr_file <- file.path(lcr_aligned_dir, "lcr_reproj_algn_natvtam.tif")
gpw_file <- file.path(gpw_aligned_dir, "gpw_reclass_reproj_algn_natvtam_noshrub.tif")
plow_file <- file.path(plow_aligned_dir, "plow_reproj_algn_natvtam.tif")

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
pli_ap <- rast(pli_ap_file)
pli_cu <- rast(pli_cu_file)
pli_mg <- rast(pli_mg_file)
pli_mm <- rast(pli_mm_file)
lcr <- rast(lcr_file)
gpw <- rast(gpw_file)
plow <- rast(plow_file)

files <- c(
  aci_file, luts_file, cec_file, ngi_file, cgr_file, eccc_file,
  esa_file, glc_file, cnvc_file, lcc_file, pudl_file, pli_ap_file,
  pli_cu_file, pli_mg_file, pli_mm_file, lcr_file, gpw_file, plow_file
)

keys <- c(
  "aci","luts","cec","ngi","cgr","eccc", "esa","glc",
  "cnvc","lcc","pudl","pli_ap","pli_cu","pli_mg","pli_mm","lcr","gpw","plow"
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

# add a field that categorizes the ground-truthed data into native vs tame vs not (only "not" is shrub)

# Normalize and derive the new field
norm <- function(x) tolower(trimws(as.character(x)))
vals <- norm(attr_df$LC_Cl_N)

# ---- 3-class ground truth: 0 = not grass, 1 = tame, 2 = native ----
attr_df$LC_class4 <- NA_integer_

# not grass (shrubs only)
attr_df$LC_class4[vals == "shrubs"] <- 0L

# tame grass
attr_df$LC_class4[vals == "tame"] <- 1L

# mixed grass
attr_df$LC_class4[vals == "mixed"] <- 2L   

# native grass
attr_df$LC_class4[vals == "native"] <- 3L

# Quick QA
table(attr_df$LC_class4, useNA = "ifany")

head(attr_df)


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
csv_out_path <- file.path(out_dir, "groundtruth_points_with_grass_raster_values_natvtamvmix_ecodist_pli_sep.csv")
write.csv(attr_df, csv_out_path, row.names = FALSE)

# filter attr_df to only include rows where Yer_Cll is 2015+
attr_df_2015 <- attr_df %>%
  filter(Yer_Cll >= 2015)

head(attr_df_2015)

csv_out_path <- file.path(out_dir, "groundtruth_points_with_grass_raster_values_natvtamvmix_pli_sep_ge2015.csv")
write.csv(attr_df_2015, csv_out_path, row.names = FALSE)

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


