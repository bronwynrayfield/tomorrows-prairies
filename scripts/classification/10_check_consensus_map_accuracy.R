
###################################################################################
##                       CHECK ACCURACY OF CONSENSUS MAP                         ##
###################################################################################

##----------------------- Alison Long - July 2026 -------------------------

# Step 10: assess accuracy of consensus map by comparing to ground-truthed point data

# load packages
library(terra)
library(sf)
library(tidyverse)
library(tools)   # for file_path_sans_ext
library(viridisLite)
library(ggspatial)

# Set the temporary directory to a drive with sufficient space
terraOptions(tempdir = "R:/CanadianGrasslands/WG_Classification/scratch")

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
cnvc_aligned_dir  <- file.path(proc_grass_data_dir, "natureserve_cnvc/aligned")
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
pglc_aligned_dir  <- file.path(proc_grass_data_dir, "aafc_pglc/aligned")

####--- FILES ---####
# Pilot area
pa_file <- file.path(pa_dir, "pilot_area_aligned_cnvc_Jan2026.tif")

# Study area
sa_file <- file.path(sa_dir, "class_study_area_aligned_cnvc_Dec2025.tif")

# Modelling study area shapefile for plotting
sap_file <- file.path(sa_dir, "modelling_study_area_aligned_cnvc_Jan2026.shp")

# ground-truthed data
gt_file <- file.path(gt_dir, "gt_points_mod_Mar26_cleanJul26.shp")

# ecodistricts
eco_file <- file.path(eco_dir, "ecodistricts_clean_class_sa.shp")

# grass raster files
aci2024_file <- file.path(aci_aligned_dir, "aafc_aci_merged_reproj_algn_binary_noshrub.tif")
aci2025_file <- file.path(aci_aligned_dir, "aafc_aci_2025_merged_reproj_algn_binary_noshrub.tif")
luts30_file  <- file.path(luts30_masked_dir, "luts30_reproj_algn_binary_masked.tif")
luts10_file  <- file.path(luts10_masked_dir, "luts10_reproj_algn_binary_masked.tif")
cec_file     <- file.path(cec_masked_dir, "nalcms_reproj_algn_binary_noshrub_masked.tif")
ngi_file     <- file.path(ngi_aligned_dir, "cfga_ngi_merged_reproj_algn_binary.tif")
cgr_file     <- file.path(cgr_masked_dir, "cgr_reproj_algn_binary_masked.tif")
eccc_file    <- file.path(eccc_masked_dir, "eccc_reclass_reproj_algn_binary_masked.tif")
esa_file     <- file.path(esa_masked_dir, "esa_reproj_algn_binary_noshrub_masked.tif")
glc_file     <- file.path(glc_masked_dir, "glc_reproj_algn_binary_noshrub_masked.tif")
cnvc_file    <- file.path(cnvc_aligned_dir, "cnvc_algn_binary.tif")
lcc_file     <- file.path(lcc_masked_dir, "lcc_reproj_algn_binary_noshrub_masked.tif")
pudl_file    <- file.path(pudl_aligned_dir, "pudl_reclass_reproj_algn_binary_noshrub.tif")
pli_ap_file  <- file.path(pli_aligned_dir, "sk_pli_ap_reproj_algn_binary_noshrub.tif")
pli_cu_file  <- file.path(pli_aligned_dir, "sk_pli_cu_reproj_algn_binary_noshrub.tif")
pli_mg_file  <- file.path(pli_aligned_dir, "sk_pli_mg_reproj_algn_binary_noshrub.tif")
pli_mm_file  <- file.path(pli_aligned_dir, "sk_pli_mm_reproj_algn_binary_noshrub.tif")
lcr_file     <- file.path(lcr_masked_dir, "lcr_reproj_algn_binary_masked.tif")
gpw_file     <- file.path(gpw_masked_dir, "gpw_reproj_algn_binary_noshrub_masked.tif")
plow_file    <- file.path(plow_masked_dir, "plow_reproj_algn_intact_masked.tif")
plvi_file    <- file.path(plvi_aligned_dir, "plvi_reproj_algn_binary_gte2015.tif")
hfi_file     <- file.path(hfi_aligned_dir, "hfi_reproj_algn_binary_gte2015.tif")
gcc_file     <- file.path(gcc_aligned_dir, "gcc_reproj_algn_binary.tif")
dw_file      <- file.path(dw_masked_dir, "dw_reproj_algn_masked_binary.tif")
esri_file    <- file.path(esri_masked_dir, "esri_reproj_algn_masked_binary.tif")
mgi_file     <- file.path(mgi_aligned_dir, "mb_mgi_reproj_algn_binary_noshrub.tif")
pglc_file    <- file.path(pglc_aligned_dir, "pglc_reproj_algn_binary.tif")

####--- READ INTO TERRA ---####

#gt <- vect(gt_file)
gt_sf <- sf::st_read(gt_file)

# keep only point geometries and remove empty geometries
#gt_sf <- gt_sf[st_geometry_type(gt_sf) == "POINT", ]
gt_sf <- gt_sf[st_geometry_type(gt_sf) == "POINT" &
                 !sf::st_is_empty(gt_sf),]

# clean
coords <- st_coordinates(gt_sf)

summary(coords)

bad <- which(
  coords[,1] < -180 | coords[,1] > 180 |
    coords[,2] < -90  | coords[,2] > 90
)

length(bad)

gt_sf[bad, ]

keep <- coords[,1] >= -180 &
  coords[,1] <= 180 &
  coords[,2] >= -90 &
  coords[,2] <= 90

gt_sf <- gt_sf[keep, ]

# convert to terra
gt <- terra::vect(gt_sf)
gt # 200,985 points

cnvc <- rast(cnvc_file) # for plotting
sap  <- vect(sap_file)  # for plotting

files <- c(
  aci2024_file, aci2025_file, luts10_file, luts30_file, cec_file, ngi_file, cgr_file, 
  eccc_file, esa_file, glc_file, cnvc_file, lcc_file, pudl_file, pli_ap_file,
  pli_cu_file, pli_mg_file, pli_mm_file, lcr_file, gpw_file, plow_file,
  plvi_file, hfi_file, gcc_file, dw_file, esri_file, mgi_file, pglc_file
)

keys <- c(
  "aci2024","aci2025","luts10","luts30","cec","ngi","cgr","eccc", "esa","glc",
  "cnvc","lcc","pudl","pli_ap","pli_cu","pli_mg","pli_mm","lcr",
  "gpw","plow","plvi","hfi","gcc","dw","esri","mgi","pglc"
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

stopifnot(length(unique(crs_all)) == 1)

stopifnot(
  nrow(unique(as.data.frame(res_all))) == 1,
  nrow(unique(as.data.frame(org_all))) == 1,
  nrow(unique(as.data.frame(dim_all))) == 1
)

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

nrow(gt_in_r_crs)      # 200985
nrow(vals_df)          # 200985

# Basic counts
n_pt   <- nrow(gt_in_r_crs)
n_ext  <- nrow(vals_df)
n_pt; n_ext

# ID problems (should be 0 in each ideally)
dups   <- vals_df$ID[duplicated(vals_df$ID)]
miss   <- setdiff(seq_len(n_pt), vals_df$ID)        # IDs present in points but absent in extract
extra  <- setdiff(vals_df$ID, seq_len(n_pt))        # IDs not in 1..nrow(points)

length(dups); length(miss); length(extra)

stopifnot(length(dups) == 0)
stopifnot(length(miss) == 0)
stopifnot(length(extra) == 0)

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
#terra::values(gt_out)[, names(r_all)] <- vals_only

v <- terra::values(gt_out)

# remove any existing columns that will be replaced
v <- v[, setdiff(names(v), names(r_all)), drop = FALSE]

#v <- cbind(v, vals_only[, names(r_all), drop = FALSE])
#terra::values(gt_out) <- v

# How many geometries vs. how many rows in the attribute table?
n_geom <- nrow(gt_out)
n_attr <- nrow(terra::values(gt_out))
n_geom; n_attr

stopifnot(nrow(gt_in_r_crs) == nrow(vals_df))
stopifnot(nrow(gt_out) == nrow(gt_in_r_crs))

# Find attribute columns whose length != number of geometries
attr_len <- sapply(terra::values(gt_out), NROW)
bad_cols <- names(attr_len)[attr_len != n_geom]
good_cols <- names(attr_len)[attr_len == n_geom]

length(bad_cols); bad_cols[1:min(10, length(bad_cols))]  

# 1) Attributes (already in gt_out after binding step)
vals <- terra::values(gt_out)               # data.frame with 200,985 rows

# 2) Coordinates for a SpatVector (points) as a data.frame
g   <- terra::geom(gt_out, df = TRUE)       # has columns 'geom','part','x','y','hole'
stopifnot(nrow(g) == nrow(gt_out))
#xy  <- g[, c("x", "y")]

xy <- terra::crds(gt_out, df = TRUE)

#xy <- terra::crds(gt_out)
#xy <- as.data.frame(xy)

dim(vals)
dim(xy)

nrow(gt_out)

# 3) Sanity check and combine
stopifnot(nrow(vals) == nrow(xy))           # should both be 200,985 for points
#attr_df <- cbind(vals, xy, stringsAsFactors = FALSE)

attr_df <- cbind(
  vals,       # 14 original attributes
  vals_only,  # 27 raster layers
  xy          # coordinates
)

dim(attr_df)

# rename coord columns
names(attr_df)[(ncol(attr_df)-1):ncol(attr_df)] <- c("x", "y")

names(attr_df)

stopifnot(nrow(attr_df) == nrow(gt_out))
# Inspect NA shares per layer (helpful for coverage/masking)
sort(sapply(attr_df[names(r_all)], function(x) mean(is.na(x))), decreasing = TRUE)[1:10]

# add a field that categorizes the ground-truthed data into grass vs not ("not" is shrub/woody and annual cropland)

# Normalize and derive the new field
norm <- function(x) tolower(trimws(as.character(x)))
stopifnot("LC_Cl_N" %in% names(attr_df))
stopifnot("Sorc_ID" %in% names(attr_df))
vals <- norm(attr_df$LC_Cl_N)

attr_df$LC_grass <- NA_integer_
attr_df$LC_grass[vals %in% c("mixed","native","tame")] <- 1L
attr_df$LC_grass[vals == "woody"] <- 0L
attr_df$LC_grass[vals == "annual cropland"] <- 0L
attr_df$LC_grass[vals == "wetland"] <- 0L
attr_df$LC_grass[vals == "open water"] <- 0L

# Quick QA
table(attr_df$LC_grass, useNA = "ifany")
#    0      1 
# 4111 196874

table(vals, useNA = "ifany")
#vals
#annual cropland           mixed          native      open water            tame         wetland 
#           1355            4112           18181              87          174581              22 
#woody 
#2647 

head(attr_df)

####--- EXPORT ---####

# write to csv
csv_out_path <- file.path(out_dir, "groundtruth_points_with_grass_raster_values_binary_ecodist_pli_sep_18Jul2026.csv")
write.csv(attr_df, csv_out_path, row.names = FALSE)

# filter attr_df to only include rows where Yer_Cll is 2015+
stopifnot("Yer_Cll" %in% names(attr_df))
attr_df_2015 <- attr_df %>%
  filter(Yer_Cll >= 2015)

head(attr_df_2015)

csv_out_path <- file.path(out_dir, "groundtruth_points_with_grass_raster_values_binary_pli_sep_ge2015_18Jul2026.csv")
write.csv(attr_df_2015, csv_out_path, row.names = FALSE)

####--- CLEAN FOR NEXT STEP ---####

# filter out any non-grass points (woody and annual cropland)
#attr_df_2015 <- attr_df_2015 %>%
#  filter(LC_grass == 1)

# add source column and fill when points were used as training data 
attr_df_2015 <- attr_df_2015 %>%
  mutate(
    Training = case_when(
      Sorc_ID == "PLI_SKMoE" ~ "pli",
      Sorc_ID == "ABMI_HFI" ~ "hfi",
      Sorc_ID == "ABMI_THF" ~ "hfi",
      Sorc_ID == "CFGA_NGI" & Yer_Cll %in% c(2022, 2023) ~ "ngi",
      Sorc_ID == "ACI_field_2024" ~ "aci2024", 
      Sorc_ID == "ACI_field_2025" ~ "aci2025", 
      Sorc_ID == "BC_ACI_2025" & LC_Cl_N == "TAME" ~ "aci2025", 
      TRUE ~ ""    # keep unmatched cases empty
    )
  )

# Move new columns to left
attr_df_2015 <- attr_df_2015 %>%
  relocate(LC_grass, Training)

# write cleaned dataset to csv
csv_out_path <- file.path(out_dir, "groundtruth_points_with_grass_raster_values_binary_pli_sep_ge2015_18Jul2026_clean.csv")
write.csv(attr_df_2015, csv_out_path, row.names = FALSE)

# load packages
library(terra)
library(sf)
library(tidyverse)
library(tools)   # for file_path_sans_ext

# Set the temporary directory to a drive with sufficient space
terraOptions(tempdir = "R:/CanadianGrasslands/WG_classification/scratch")

####--- DIRECTORIES ---####

# pilot area
#pa_dir <- "R:/CanadianGrasslands/processeddatasets/pilotarea"

# study area
#sa_dir <- "R:/CanadianGrasslands/processeddatasets/studyarea"

# grasslands
#raw_grass_data_dir <- "R:/CanadianGrasslands/existingdatasets"
proc_data_dir <- "R:/CanadianGrasslands/WG_Classification/processed"
results_dir <- "R:/CanadianGrasslands/WG_Classification/results"
out_dir <- "R:/CanadianGrasslands/WG_Classification/results/consensus_mapping"

# ground-truthed data
gt_dir <- "R:/CanadianGrasslands/processeddatasets/groundtruth_data"

# ecodistricts 
eco_dir <- "R:/CanadianGrasslands/processeddatasets/boundaries/national-ecological-framework-for-canada/Ecodistricts"

# consensus map directory
cm_dir  <- file.path(results_dir, "consensus_mapping/v4_temporal_modelling_extent_24Jul2026")

####--- FILES ---####
# Pilot area
#pa_file <- file.path(pa_dir, "pilot_area_aligned_cnvc_Jan2026.tif")

# Study area
#sa_file <- file.path(sa_dir, "class_study_area_aligned_cnvc_Dec2025.tif")

# ground-truthed data
gt_file <- file.path(gt_dir, "gt_points_mod_Mar26_cleanJul26.shp")

# ecodistricts
eco_file <- file.path(eco_dir, "ecodistricts_clean_class_sa.shp")

# consensus map
cm_file <- file.path(cm_dir, "grass_3class.tif")

####--- READ IN DATA ---####

#gt <- vect(gt_file)
#gt
gt_sf <- sf::st_read(gt_file)

cm_rast <- rast(cm_file)

####--- PREP DATA ---####

# keep only point geometries and remove empty geometries
#gt_sf <- gt_sf[st_geometry_type(gt_sf) == "POINT", ]
gt_sf <- gt_sf[st_geometry_type(gt_sf) == "POINT" &
                 !sf::st_is_empty(gt_sf),]

# clean
coords <- st_coordinates(gt_sf)

summary(coords)

bad <- which(
  coords[,1] < -180 | coords[,1] > 180 |
    coords[,2] < -90  | coords[,2] > 90
)

length(bad)

gt_sf[bad, ]

keep <- coords[,1] >= -180 &
  coords[,1] <= 180 &
  coords[,2] >= -90 &
  coords[,2] <= 90

gt_sf <- gt_sf[keep, ]

# convert to terra
gt <- terra::vect(gt_sf)
gt # 200,985 points

# ---- Ensure points are terra SpatVector ----
gt_vec <- if (inherits(gt, "sf")) terra::vect(gt) else gt
stopifnot(inherits(gt_vec, "SpatVector"))
stopifnot(geomtype(gt_vec) == "points")

# ---- Reproject points once if needed ----
gt_in_r_crs <- if (!identical(crs(gt_vec), crs(cm_rast))) {
  project(gt_vec, crs(cm_rast))
} else {
  gt_vec
}

####--- EXTRACT AND BIND DATA ---####

# ---- Single extract for all layers ----
vals_df <- terra::extract(cm_rast, gt_in_r_crs, ID = TRUE)

nrow(gt_in_r_crs)      # 200985
nrow(vals_df)          # 200985

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
# vals_df: result of terra::extract(cm_rast, gt_in_r_crs, ID = TRUE)
# gt_in_r_crs: SpatVector points used in the extract()
# cm_rast: multi-layer SpatRaster; names(cm_rast) are the desired column names

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
#terra::values(gt_out)[, names(cm_rast)] <- vals_only

v <- terra::values(gt_out)

# remove any existing columns that will be replaced
v <- v[, setdiff(names(v), names(cm_rast)), drop = FALSE]

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
stopifnot(nrow(g) == nrow(gt_out))
#xy  <- g[, c("x", "y")]

xy <- terra::crds(gt_out, df = TRUE)

dim(vals)
dim(xy)

nrow(gt_out)

# 3) Sanity check and combine
stopifnot(nrow(vals) == nrow(xy))           # should both be 196,931 for points
attr_df <- cbind(vals, xy, stringsAsFactors = FALSE)

# rename coord columns
names(attr_df)[(ncol(attr_df)-1):ncol(attr_df)] <- c("x", "y")

stopifnot(nrow(attr_df) == nrow(gt_out))
# Inspect NA shares per layer (helpful for coverage/masking)
#sort(sapply(attr_df[names(cm_rast)], function(x) mean(is.na(x))), decreasing = TRUE)[1:10]

# add a field that categorizes the ground-truthed data into native vs tame vs not (only "not" is shrub)

# Normalize and derive the new field
norm <- function(x) tolower(trimws(as.character(x)))
vals <- norm(attr_df$LC_Cl_N)

# ---- 3-class ground truth: 0 = not grass, 1 = tame, 2 = native ----
attr_df$LC_class3 <- NA_integer_

# not grass (shrubs/woody)
attr_df$LC_class3[vals == "woody"] <- 0L

# not grass (cropland)
attr_df$LC_class3[vals == "annual cropland"] <- 0L

# not grass (open water)
attr_df$LC_class3[vals == "open water"] <- 0L

# not grass (wetland)
attr_df$LC_class3[vals == "wetland"] <- 0L

# tame grass
attr_df$LC_class3[vals == "tame"] <- 1L

# treat "mixed" as tame
attr_df$LC_class3[vals == "mixed"] <- 1L   

# native grass
attr_df$LC_class3[vals == "native"] <- 2L

# Quick QA
table(attr_df$LC_class3, useNA = "ifany")

head(attr_df)

attr_df %>%
  filter(is.na(LC_class3)) %>%
  print()


####--- EXPORT ---####

# write to csv
csv_out_path <- file.path(out_dir, "groundtruth_points_with_grass_raster_values_natvtam_ecodist_pli_sep_consmap_24Jul2026.csv")
write.csv(attr_df, csv_out_path, row.names = FALSE)

# filter attr_df to only include rows where Yer_Cll is 2015+
attr_df_2015 <- attr_df %>%
  filter(Yer_Cll >= 2015)

head(attr_df_2015)

csv_out_path <- file.path(out_dir, "groundtruth_points_with_grass_raster_values_natvtam_pli_sep_ge2015_consmap_24Jul2026.csv")
write.csv(attr_df_2015, csv_out_path, row.names = FALSE)

####--- CLEAN FOR NEXT STEP ---####

# filter out any non-grass points (woody and annual cropland)
#attr_df_2015_clean <- attr_df_2015 %>%
#  filter(LC_class3 %in% c(1, 2))

# add source column and fill when points were used as training data 
# attr_df_2015 <- attr_df_2015 %>%
#   mutate(
#     Training = case_when(
#       Sorc_ID == "PLI_SKMoE" ~ "pli",
#       Sorc_ID == "ABMI_HFI" ~ "hfi",
#       Sorc_ID == "ABMI_THF" ~ "hfi",
#       Sorc_ID == "CFGA_NGI" & Yer_Cll %in% c(2022, 2023) ~ "ngi",
#       Sorc_ID == "ACI_field_2024" ~ "aci2024", 
#       Sorc_ID == "ACI_field_2025" ~ "aci2025", 
#       Sorc_ID == "BC_ACI_2025" & LC_Cl_N == "TAME" ~ "aci2025", 
#       TRUE ~ ""    # keep unmatched cases empty
#     )
#   )

# Move new columns to left
attr_df_2015 <- attr_df_2015 %>%
  relocate(LC_class3)

# write cleaned dataset to csv
csv_out_path <- file.path(out_dir, "groundtruth_points_with_grass_raster_values_natvtam_pli_sep_ge2015_Mar2026_clean_v2.csv")
write.csv(attr_df_2015, csv_out_path, row.names = FALSE)

#############################################################################################
#############################################################################################

# load packages
library(readxl)
library(writexl)
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(purrr)
library(tibble)
library(sf)
library(scales)
library(mltools)

# Small utility: "x %||% y" returns y if x is NULL
`%||%` <- function(x, y) if (is.null(x)) y else x

# ----------------------------
# 2) File paths 
# ----------------------------
#points_excel_path <- "R:/CanadianGrasslands/WG_Classification/processed/groundtruthed_data_points/groundtruth_points_with_grass_raster_values_natvtam_pli_sep_ge2015_18Jul2026_clean_withnongrass.xlsx"
metadata_csv_path <- "R:/CanadianGrasslands/WG_Classification/github/tomorrows-prairies/processeddatasets/Grass/metadata_full_extent_nattam_pli_sep.csv"

#out_dir <- "R:/CanadianGrasslands/WG_Classification/results/accuracy_assessment/native_v_tame/compiled_derived"
#dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ----------------------------
# 3) Configuration
# ----------------------------
# IMPORTANT: Handle PLI as four separate rasters in the points file and metadata.
pli_variants <- c("pli_ap", "pli_cu", "pli_mg", "pli_mm")

# Expected raster columns (3-class values 0/1/2)
raster_cols <- c("which.max")

# Grouping helper: collapse pli_* -> "pli"; all others remain themselves
raster_to_group <- function(x) ifelse(x %in% pli_variants, "pli", x)

# ----------------------------
# 4) Load points check for required columns
# ----------------------------
#points <- readxl::read_excel(points_excel_path)
# points <- readxl::read_excel(
#   points_excel_path,
#   col_types = c("numeric", "text", "text", "numeric", rep("guess", 10), rep("numeric", 20))
# )

points <- attr_df_2015

names(points)
str(points)
unique(points$LC_class3)
table(points$LC_class3, useNA = "ifany")

# Add point ID column
points <- points %>%
  mutate(Pt_ID = row_number()) %>%
  relocate(Pt_ID)

# Drop rows where ecodistrict ID is 0
points <- points %>%
  filter(Ecds_ID != 0)

unique(points$Ecrgn_n)

# Create column for ecoprovince
points <- points %>%
  mutate(
    Ecprv = case_when(
      Ecrgn_n %in% c("Abitibi Plains", 
                     "Lake Nipigon",
                     "Rupert River Plateau",
                     "Lac Seul Upland") ~ "6.2",
      Ecrgn_n %in% c("Thunder Bay-Quetico", 
                     "Lake of the Woods",
                     "Rainy River") ~ "6.5",
      Ecrgn_n %in% c("Algonquin-Lake Nipissing", 
                     "Southern Laurentians",
                     "Lac Temiscamingue Lowland") ~ "6.6",
      Ecrgn_n %in% c("St-Laurent Lowlands", 
                     "St. Lawrence Lowlands", 
                     "Manitoulin-Lake Simcoe",
                     "Frontenac Axis") ~ "8.1",
      Ecrgn_n == "Lake Erie Lowland" ~ "8.2",
      Ecrgn_n %in% c("Western Alberta Upland", 
                     "Clear Hills Upland") ~ "9.1",
      Ecrgn_n %in% c("Wabasca Lowland", 
                     "Boreal Transition",
                     "Western Boreal",
                     "Peace Lowland",
                     "Mid-Boreal Uplands") ~ "9.2",
      Ecrgn_n %in% c("Interlake Plain", 
                     "Mid-Boreal Lowland") ~ "9.3",
      Ecrgn_n == "Lake Manitoba Plain" ~ "10.1",
      Ecrgn_n %in% c("Aspen Parkland", 
                     "Southwest Manitoba Uplands",
                     "Southwest Manotoba Uplands") ~ "10.2",
      Ecrgn_n %in% c("Mixed Grassland", 
                     "Moist Mixed Grassland",
                     "Fescue Grassland",
                     "Cypress Upland") ~ "10.3",
      Ecrgn_n %in% c("Lower Mainland", 
                     "Eastern Vancouver Island") ~ "13.1",
      Ecrgn_n %in% c("Pacific Ranges", 
                     "Cascade Ranges") ~ "13.2",
      Ecrgn_n %in% c("Central Canadian Rocky Mountains", 
                     "Skeena Mountains",
                     "Fraser Basin") ~ "14.1",
      Ecrgn_n %in% c("Fraser Plateau",
                     "Bulkley Ranges") ~ "14.2",
      Ecrgn_n %in% c("Thompson-Okanagan Plateau", 
                     "Okanogan Highland",
                     "Okanogan Range",
                     "Okanagan Highland",
                     "Okanagan Range",
                     "Interior Transition Ranges") ~ "14.3",
      Ecrgn_n %in% c("Columbia Mountains and Highlands", 
                     "Southern Rocky Mountain Trench",
                     "Northern Continental Divide",
                     "Eastern Continental Ranges",
                     "Western Continental Ranges",
                     "Selkirk-Bitterroot Foothills") ~ "14.4"))

unique(points$Ecprv)

# Define required columns
required_cols <- c("Pt_ID", "LC_class3", "Sorc_ID", "Ecds_ID", "Ecrg_ID", "Ecrgn_n", "Ecprv", "Yer_Cll", "LC_Cl_N")
missing_req <- setdiff(required_cols, names(points))
if (length(missing_req) > 0) {
  stop("Missing required columns in points file: ", paste(missing_req, collapse = ", "))
}

#Check LC types present
table(points$LC_Cl_N, useNA = "ifany")

# Rename all non-grass points as OTHER
points <- points %>%
  mutate(
    LC_Cl_N = if_else(
      LC_Cl_N %in% c("ANNUAL CROPLAND", "WOODY", "WETLAND", "OPEN WATER"),
      "OTHER",
      LC_Cl_N))

# Rename all MIXED points as TAME
points <- points %>%
  mutate(
    LC_Cl_N = if_else(
      LC_Cl_N == "MIXED",
      "TAME",
      LC_Cl_N))

#Check LC types present
table(points$LC_Cl_N, useNA = "ifany")

# How many OTHER points for each year?
pts.oth <- points %>%
  filter(LC_Cl_N == "OTHER")

pts.oth %>%
  group_by(Yer_Cll) %>%
  summarise(total_count = n(), .groups = "drop")

# How many TAME points for each year?
pts.tam <- points %>%
  filter(LC_Cl_N == "TAME")

pts.tam %>%
  group_by(Yer_Cll) %>%
  summarise(total_count = n(), .groups = "drop")

# How many NATIVE points for each year?
pts.nat <- points %>%
  filter(LC_Cl_N == "NATIVE")

pts.nat %>%
  group_by(Yer_Cll) %>%
  summarise(total_count = n(), .groups = "drop")

#Check LC types present
table(points$LC_class3, useNA = "ifany")

# Identify which raster columns actually exist in the points file
present_rasters <- intersect(raster_cols, names(points))
if (length(present_rasters) == 0) {
  stop("None of the expected raster columns were found in the points file.")
}

# ----------------------------
# 5) Load metadata and allow MULTIPLE years per raster
# ----------------------------
meta_raw <- readr::read_csv(metadata_csv_path, show_col_types = FALSE)

parse_year_vector <- function(y) {
  # Accept numbers or delimited strings "2020|2021" or "2020, 2021"
  if (is.null(y) || (length(y) == 1 && is.na(y))) return(integer(0))
  if (is.numeric(y)) return(as.integer(stats::na.omit(y)))
  y_chr <- as.character(y)
  y_chr <- stringr::str_trim(y_chr)
  y_chr[y_chr %in% c("", "NA", "NaN", "null", "NULL", "missing")] <- NA_character_
  y_chr <- y_chr[!is.na(y_chr)]
  if (length(y_chr) == 0) return(integer(0))
  parts <- unlist(strsplit(y_chr, "[|,]"))
  parts <- stringr::str_trim(parts)
  suppressWarnings(as.integer(stats::na.omit(as.numeric(parts))))
}

if (all(c("raster", "year") %in% names(meta_raw))) {
  # Long format
  meta_norm <- meta_raw %>%
    transmute(
      raster = str_trim(tolower(as.character(raster))),
      year_raw = year
    ) %>%
    group_by(raster) %>%
    summarise(years = list(parse_year_vector(year_raw)), .groups = "drop")
} else if (nrow(meta_raw) == 1 && any(intersect(raster_cols, names(meta_raw)) %in% names(meta_raw))) {
  # Wide single-row
  meta_norm <- tibble(
    raster = raster_cols,
    years  = map(raster_cols, ~ parse_year_vector(meta_raw[[.x]]))
  ) %>%
    mutate(raster = tolower(raster))
} else {
  stop(
    "Metadata format not recognized. Provide either:\n",
    "  (A) Long CSV with columns: raster, year (multiple rows allowed per raster), OR\n",
    "  (B) Wide 1-row CSV with raster columns; cells may be single years or delimited lists (e.g., '2020|2021')."
  )
}

# Build named list: raster -> unique sorted integer vector of years
raster_years_map <- meta_norm %>%
  mutate(years = map(years, ~ unique(sort(as.integer(.x))))) %>%
  deframe()

# Warn if any present raster has no metadata years
rasters_missing_years <- setdiff(present_rasters, names(raster_years_map))
if (length(rasters_missing_years) > 0) {
  warning("These rasters are present in points but have no years in metadata (they will be excluded): ",
          paste(rasters_missing_years, collapse = ", "))
}

# ----------------------------
# 6) Reshape points to long and attach allowed years (list-column)
# ----------------------------
# Convert raster columns to integer 0/1/2/NA. Anything not 0/1/2 becomes NA.
points_norm <- points %>%
  mutate(
    across(
      all_of(present_rasters),
      ~ {
        x <- as.character(.)
        x <- str_trim(x)
        x[x %in% c("", "NA", "NaN", "null", "NULL")] <- NA_character_
        ifelse(x %in% c("0","1","2"), as.integer(x), NA_integer_)
      }
    )
  )

# Long format: one row per point × raster
# Also attach:
#   - raster_group (pli_* -> pli)
#   - allowed_years (list of integers from metadata)
long <- points_norm %>%
  select(all_of(required_cols), all_of(present_rasters)) %>%
  pivot_longer(
    cols = all_of(present_rasters),
    names_to = "raster",
    values_to = "pred"
  ) %>%
  mutate(
    raster = tolower(raster),
    raster_group = raster_to_group(raster),
    allowed_years = map(raster, ~ raster_years_map[[.x]] %||% integer(0))
  )

# Collapse lists to strings and export
long_out <- long %>%
  mutate(allowed_years = map_chr(allowed_years, ~ paste(.x, collapse = ";")))

write.csv(long_out, file.path(out_dir, "long_prefiltering_27Jul2026_natvtam_withnongrass_compiled.csv"))

# ----------------------------
# 7) Filter rows via class-based time rules using MULTIPLE allowed years
# ----------------------------
# Set up class-dependent year matching rules:
#   Native grass (lc3 == 2L): accept any year ≤ reference year (choose most recent such year)
#   Tame grass (lc3 == 1L): accept only if year matches exactly
#   Non-grass (lc3 == 0L): can configure rule to be like native or tame rule 
#compute_matched_year_class3 <- function(yrs, yer, lc3, notgrass_rule = c("exact","lte")) {
#  notgrass_rule <- match.arg(notgrass_rule)
#  yrs <- as.integer(yrs)
#  if (length(yrs) == 0L || is.na(yer) || is.na(lc3)) return(NA_integer_)
#  if (lc3 == 2L) {                # native
#    ok <- yrs[yrs <= yer]
#    if (length(ok) == 0L) return(NA_integer_)
#    return(max(ok))
#  } else if (lc3 == 1L) {         # tame
#    if (yer %in% yrs) return(yer) else return(NA_integer_)
#  } else if (lc3 == 0L) {         # not grass
#    if (notgrass_rule == "exact") {
#      if (yer %in% yrs) return(yer) else return(NA_integer_)
#    } else {"lte"
#      ok <- yrs[yrs <= yer]
#      if (length(ok) == 0L) return(NA_integer_)
#      return(max(ok))
#    }
#  } else {
#    return(NA_integer_)
#  }
#}

compute_matched_year_class3 <- function(yrs, yer, lc3) {
  yrs <- as.integer(unlist(yrs))
  if (length(yrs) == 0L || is.na(yer) || is.na(lc3)) return(NA_integer_)
  if (lc3 == 2L) { # native
    ok <- yrs[yrs <= yer]
    if (length(ok) == 0L) return(NA_integer_)
    return(max(ok))
  } else if (lc3 == 1L) { # tame
    if (yer %in% yrs) return(yer) else return(NA_integer_)
  } else if (lc3 == 0L) { # nongrass
    if (yer %in% yrs) return(yer) else return(NA_integer_)
  } else {
    return(NA_integer_)
  }
}

# Apply filtering:
#   - require LC_class3 ∈ {0,1,2}, Yer_Cll present, pred present
#   - DO NOT filter out LC_class3 == 0
#   - DO NOT exclude LC_Cl_N == "SHRUB"
#   - require allowed_years matching via matched_year
#   - Exclude training sources (updated filtering)
eval_data <- long %>%
  filter(
    !is.na(LC_class3),
    LC_class3 %in% c(0L, 1L, 2L),  # keep all 3 classes
    !is.na(Yer_Cll)) %>%
  filter(!is.na(pred)) %>%
  rowwise() %>%
  mutate(matched_year = compute_matched_year_class3(allowed_years, Yer_Cll, LC_class3)) %>%
  ungroup() %>%
  filter(!is.na(matched_year)) %>%
  filter(is.na(Training) | (Training != raster & Training != raster_group)) # updated

# Collapse lists to strings and export
eval_data_out <- eval_data %>%
  mutate(allowed_years = map_chr(allowed_years, ~ paste(.x, collapse = ";")))

write.csv(eval_data_out, file.path(out_dir, "eval_data_27Jul2026_natvtam_withnongrass_compiled.csv"))

# Isolate datasets that need to be handled differently
#   Datasets with training data need to be filtered separately
#   so that data is not lost

# ACI
#eval_data_aci <- eval_data %>%
#  filter(raster_group == "aci")

# NGI
#eval_data_ngi <- eval_data %>%
#  filter(raster_group == "ngi")

# PLI
#eval_data_pli <- eval_data %>%
#  filter(raster_group == "pli")

# Filter out training data
#   Remove points used to train the model

# ACI
#eval_data_aci_filt <- eval_data_aci %>%
#  filter(is.na(Source) | Source != "aci")

# HFI
#eval_data_hfi_filt <- eval_data_hfi %>%
#  filter(is.na(Source) | Source != "hfi")

# NGI
#eval_data_ngi_filt <- eval_data_ngi %>%
#  filter(is.na(Source) | Source != "ngi")

# PLI
#eval_data_pli_filt <- eval_data_pli %>%
#  filter(is.na(Source) | Source != "pli")

# Filter ACI, HFI, NGI, and PLI out of eval_data
#eval_data_reg <- eval_data %>%
#  filter(!(raster_group %in% c("aci", "pli", "ngi", "hfi")))

# Combine back together
#eval_data_combined <- bind_rows(
#  eval_data_reg,
#  eval_data_aci_filt,
#  eval_data_hfi_filt,
#  eval_data_ngi_filt,
#  eval_data_pli_filt)

# Check number of rows - 60,049 total
#sapply(
#  list(
#    table1 = eval_data_reg,
#    table2 = eval_data_aci_filt,
#    table3 = eval_data_hfi_filt,
#    table4 = eval_data_ngi_filt,
#    table5 = eval_data_pli_filt),
#  nrow)

# Collapse lists to strings and export
#eval_data_combined_out <- eval_data_combined %>%
#  mutate(allowed_years = map_chr(allowed_years, ~ paste(.x, collapse = ";")))

#write.csv(eval_data_combined_out, file.path(out_dir, "eval_data_29Apr2026_withnongrass_good.csv"))

# ----------------------------
# 8) Multiclass metrics helpers
# ----------------------------

safe_div <- function(n, d) ifelse(d == 0 | is.na(d), NA_real_, n / d)

# ----------------------------
# MCC helpers
# ----------------------------
binary_mcc <- function(tp, tn, fp, fn) {
  
  # force numeric (double) to avoid integer overflow
  tp <- as.numeric(tp)
  tn <- as.numeric(tn)
  fp <- as.numeric(fp)
  fn <- as.numeric(fn)
  
  num <- tp * tn - fp * fn
  
  den <- sqrt(
    (tp + fp) *
      (tp + fn) *
      (tn + fp) *
      (tn + fn)
  )
  
  # robust guard
  if (is.na(den) || den == 0) return(NA_real_)
  
  num / den
}

multiclass_mcc <- function(cm) {
  
  cm <- as.matrix(cm)
  if (!is.numeric(cm)) stop("Confusion matrix must be numeric")
  
  n <- sum(cm)
  if (n == 0) return(NA_real_)
  
  r <- rowSums(cm)
  c <- colSums(cm)
  
  sum_diag <- sum(diag(cm))
  sum_row_col <- sum(r * c)
  
  numerator <- sum_diag * n - sum_row_col
  
  denom_row <- n^2 - sum(r^2)
  denom_col <- n^2 - sum(c^2)
  
  if (denom_row <= 0 || denom_col <= 0) return(NA_real_)
  
  numerator / sqrt(denom_row * denom_col)
}

# ----------------------------
# Build 3x3 confusion matrix
# ----------------------------
compute_confusion_3 <- function(df,
                                truth_col = "LC_class3",
                                pred_col  = "pred") {
  
  lvls <- c(0L, 1L, 2L)
  
  truth <- factor(df[[truth_col]], levels = lvls)
  pred  <- factor(df[[pred_col]],  levels = lvls)
  
  tab <- table(truth, pred)
  tab <- as.matrix(tab)[as.character(lvls), as.character(lvls), drop = FALSE]
  
  counts_long <- as_tibble(as.data.frame(as.table(tab))) %>%
    transmute(
      Truth = as.integer(as.character(truth)),
      Pred  = as.integer(as.character(pred)),
      N     = as.integer(Freq)
    )
  
  list(tab = tab, n = sum(tab), counts_long = counts_long)
}

# ----------------------------
# Compute metrics from 3x3 CM
# ----------------------------
compute_metrics_3 <- function(tab) {
  
  n     <- sum(tab)
  diagv <- diag(tab)
  rows  <- rowSums(tab)   # support (truth)
  cols  <- colSums(tab)   # predicted totals
  
  # ---- Accuracy ----
  overall <- safe_div(sum(diagv), n)
  
  # ---- Precision / Recall / F1 ----
  prec <- safe_div(diagv, cols)
  rec  <- safe_div(diagv, rows)
  f1   <- ifelse(
    is.na(prec) | is.na(rec) | (prec + rec) == 0,
    NA_real_,
    2 * prec * rec / (prec + rec)
  )
  
  macro_prec <- mean(prec, na.rm = TRUE)
  macro_rec  <- mean(rec,  na.rm = TRUE)
  macro_f1   <- mean(f1,   na.rm = TRUE)
  
  balanced_acc <- macro_rec
  
  weights     <- safe_div(rows, sum(rows))
  weighted_f1 <- sum(weights * f1, na.rm = TRUE)
  
  # ---- Cohen's kappa ----
  pe     <- safe_div(sum(rows * cols), n * n)
  kappa <- ifelse(is.na(pe) | (1 - pe) == 0, NA_real_, (overall - pe) / (1 - pe))
  
  # ---- Multiclass MCC ----
  mcc_multiclass <- multiclass_mcc(tab)
  
  # ---- One-vs-rest MCC ----
  classes <- seq_len(nrow(tab))
  
  ovr_mcc <- sapply(classes, function(i) {
    tp <- tab[i, i]
    fn <- sum(tab[i, -i])
    fp <- sum(tab[-i, i])
    tn <- sum(tab[-i, -i])
    binary_mcc(tp, tn, fp, fn)
  })
  
  mcc_ovr_macro <- mean(ovr_mcc, na.rm = TRUE)
  
  tibble(
    N_Evaluated        = as.integer(n),
    Overall_Accuracy  = overall,
    Balanced_Accuracy = balanced_acc,
    Kappa             = kappa,
    
    MCC_Multiclass    = mcc_multiclass,
    MCC_OvR_Macro     = mcc_ovr_macro,
    MCC_OvR_0         = ovr_mcc[1],
    MCC_OvR_1         = ovr_mcc[2],
    MCC_OvR_2         = ovr_mcc[3],
    
    Macro_Precision   = macro_prec,
    Macro_Recall      = macro_rec,
    Macro_F1          = macro_f1,
    Weighted_F1       = weighted_f1,
    
    Precision_0 = prec[1], Recall_0 = rec[1], F1_0 = f1[1], Support_0 = rows[1],
    Precision_1 = prec[2], Recall_1 = rec[2], F1_1 = f1[2], Support_1 = rows[2],
    Precision_2 = prec[3], Recall_2 = rec[3], F1_2 = f1[3], Support_2 = rows[3]
  )
}

# ----------------------------
# Summarize metrics per group
# ----------------------------
summarize_block_multiclass <- function(df, key_cols) {
  
  df %>%
    group_by(across(all_of(key_cols))) %>%
    group_modify(~ {
      c3  <- compute_confusion_3(.x)
      met <- compute_metrics_3(c3$tab)
      met
    }) %>%
    ungroup()
}

# ----------------------------
# Return long-form confusion
# ----------------------------
confusion_block_multiclass <- function(df, key_cols) {
  
  df %>%
    group_by(across(all_of(key_cols))) %>%
    group_modify(~ compute_confusion_3(.x)$counts_long) %>%
    ungroup()
}

# ----------------------------
# Create output tables
# ----------------------------
metrics_by_dataset <- summarize_block_multiclass(
  eval_data,
  key_cols = "raster_group"
)

confusions_by_dataset <- confusion_block_multiclass(
  eval_data,
  key_cols = "raster_group"
)

metrics_by_dataset

# ----------------------------
# Export output tables
# ----------------------------
out_xlsx <- file.path(out_dir, "nat_v_tame_v_nongrass_accuracy_assessment_results_27Jul2026.xlsx")

write_xlsx(
  list(
    "Metrics_by_Dataset"   = metrics_by_dataset,
    "Confusion_Long_Form"  = confusions_by_dataset
  ),
  path = out_xlsx
)

#####################################################################################
#####################################################################################

# ----------------------------
# 13) Accuracy by Ecodistrict (TRICLASS)
# ----------------------------

if (!"Ecds_ID" %in% names(eval_data)) {
  stop("Ecds_ID column not found in eval_data.")
}

# ---- ATOMIC: raster × ecodistrict ----
accuracy_by_ecodistrict_atomic <- summarize_block_multiclass(
  eval_data,
  key_cols = c("raster", "Ecds_ID")
) %>%
  arrange(raster, Ecds_ID)

# ---- GROUPED: collapse pli_* -> pli ----
accuracy_by_ecodistrict_grouped <- eval_data %>%
  mutate(raster = raster_group) %>%
  summarize_block_multiclass(
    key_cols = c("raster", "Ecds_ID")
  ) %>%
  arrange(raster, Ecds_ID)

# ---- OPTIONAL: atomic with year ----
accuracy_by_ecodistrict_year_atomic <- summarize_block_multiclass(
  eval_data,
  key_cols = c("raster", "matched_year", "Ecds_ID")
) %>%
  rename(raster_year = matched_year) %>%
  arrange(raster, raster_year, Ecds_ID)

# ---- Save ----
readr::write_csv(accuracy_by_ecodistrict_grouped,
                 file.path(out_dir, "accuracy_by_ecodistrict_grouped_compiled.csv"))

readr::write_csv(accuracy_by_ecodistrict_atomic,
                 file.path(out_dir, "accuracy_by_ecodistrict_atomic_compiled.csv"))

readr::write_csv(accuracy_by_ecodistrict_year_atomic,
                 file.path(out_dir, "accuracy_by_ecodistrict_year_atomic_compiled.csv"))

# ----------------------------
# 14) Multiclass fallback hierarchy
# ecodistrict -> ecoregion -> ecoprovince
# ----------------------------

min_n <- 20

# ---- Ecoregion metrics ----
accuracy_by_ecoregion <- summarize_block_multiclass(
  eval_data,
  key_cols = c("raster", "Ecrgn_n")
)

# ---- Ecoprovince metrics ----
accuracy_by_ecoprovince <- summarize_block_multiclass(
  eval_data,
  key_cols = c("raster", "Ecprv")
) %>%
  rename_with(
    ~ paste0(.x, "_ecoP"),
    -c(raster, Ecprv)
  )

# ---- Mapping table ----
ecodist_to_regions <- eval_data %>%
  distinct(Ecds_ID, Ecrgn_n, Ecprv)

# --------------------------------------------------
# ATOMIC
# --------------------------------------------------

accuracy_by_ecodistrict_atomic2 <- accuracy_by_ecodistrict_atomic %>%
  left_join(ecodist_to_regions, by = "Ecds_ID")

fallback_joined_atomic <- accuracy_by_ecodistrict_atomic2 %>%
  left_join(
    accuracy_by_ecoregion,
    by = c("raster", "Ecrgn_n"),
    suffix = c("_ecoD", "_ecoR")
  ) %>%
  left_join(
    accuracy_by_ecoprovince,
    by = c("raster", "Ecprv")
  )

# Make sure suffixes are present 
names(fallback_joined_atomic)

metric_cols <- setdiff(
  names(accuracy_by_ecodistrict_atomic),
  c("raster", "Ecds_ID")
)

all(paste0(metric_cols, "_ecoD") %in% names(fallback_joined_atomic))
all(paste0(metric_cols, "_ecoR") %in% names(fallback_joined_atomic))
all(paste0(metric_cols, "_ecoP") %in% names(fallback_joined_atomic))

accuracy_with_fallback_atomic <- fallback_joined_atomic %>%
  mutate(
    source_level = case_when(
      N_Evaluated_ecoD >= min_n ~ "ecodistrict",
      N_Evaluated_ecoR >= min_n ~ "ecoregion",
      TRUE ~ "ecoprovince"
    )
  )

for (m in metric_cols) {
  accuracy_with_fallback_atomic[[m]] <- case_when(
    accuracy_with_fallback_atomic$source_level == "ecodistrict" ~
      accuracy_with_fallback_atomic[[paste0(m, "_ecoD")]],
    accuracy_with_fallback_atomic$source_level == "ecoregion" ~
      accuracy_with_fallback_atomic[[paste0(m, "_ecoR")]],
    TRUE ~
      accuracy_with_fallback_atomic[[paste0(m, "_ecoP")]]
  )
}

#accuracy_with_fallback_atomic <- fallback_joined_atomic %>%
#  mutate(
#    source_level = case_when(
#      N_Evaluated_ecoD >= min_n ~ "ecodistrict",
#      N_Evaluated_ecoR >= min_n ~ "ecoregion",
#      TRUE ~ "ecoprovince"
#    )
#  ) %>%
#  mutate(
#    across(
#      all_of(metric_cols),
#      ~ case_when(
#        source_level == "ecodistrict" ~ .data[[paste0(cur_column(), "_ecoD")]],
#        source_level == "ecoregion" ~ .data[[paste0(cur_column(), "_ecoR")]],
#        TRUE ~ .data[[paste0(cur_column(), "_ecoP")]]
#      ),
#      .names = "{.col}"
#    )
#  )

# --------------------------------------------------
# GROUPED
# --------------------------------------------------

accuracy_by_ecodistrict_grouped2 <- accuracy_by_ecodistrict_grouped %>%
  left_join(ecodist_to_regions, by = "Ecds_ID")

fallback_joined_grouped <- accuracy_by_ecodistrict_grouped2 %>%
  left_join(
    accuracy_by_ecoregion,
    by = c("raster", "Ecrgn_n"),
    suffix = c("_ecoD", "_ecoR")
  ) %>%
  left_join(
    accuracy_by_ecoprovince,
    by = c("raster", "Ecprv")
  )

# Make sure suffixes are present 
names(fallback_joined_grouped)

#accuracy_with_fallback_grouped <- fallback_joined_grouped %>%
#  mutate(
#    source_level = case_when(
#      N_Evaluated_ecoD >= min_n ~ "ecodistrict",
#      N_Evaluated_ecoR >= min_n ~ "ecoregion",
#      TRUE ~ "ecoprovince"
#    )
#  ) %>%
#  mutate(
#    across(
#      all_of(metric_cols),
#      ~ case_when(
#        source_level == "ecodistrict" ~ .data[[paste0(cur_column(), "_ecoD")]],
#        source_level == "ecoregion" ~ .data[[paste0(cur_column(), "_ecoR")]],
#        TRUE ~ .data[[paste0(cur_column(), "_ecoP")]]
#      ),
#      .names = "{.col}"
#    )
#  )

accuracy_with_fallback_grouped <- fallback_joined_grouped %>%
  mutate(
    source_level = case_when(
      N_Evaluated_ecoD >= min_n ~ "ecodistrict",
      N_Evaluated_ecoR >= min_n ~ "ecoregion",
      TRUE ~ "ecoprovince"
    )
  )

for (m in metric_cols) {
  accuracy_with_fallback_grouped[[m]] <- case_when(
    accuracy_with_fallback_grouped$source_level == "ecodistrict" ~
      accuracy_with_fallback_grouped[[paste0(m, "_ecoD")]],
    accuracy_with_fallback_grouped$source_level == "ecoregion" ~
      accuracy_with_fallback_grouped[[paste0(m, "_ecoR")]],
    TRUE ~
      accuracy_with_fallback_grouped[[paste0(m, "_ecoP")]]
  )
}

# ----------------------------
# Export multiclass fallback results
# ----------------------------

readr::write_csv(
  accuracy_with_fallback_atomic,
  file.path(
    out_dir,
    "accuracy_by_ecodistrict_atomic_nattam_ecoregprovfallback_compiled_27Jul2026.csv"
  )
)

readr::write_csv(
  accuracy_with_fallback_grouped,
  file.path(
    out_dir,
    "accuracy_by_ecodistrict_grouped_nattam_ecoregprovfallback_compiled_18Jul2026.csv"
  )
)


