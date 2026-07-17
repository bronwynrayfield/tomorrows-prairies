
###################################################################################
##               PRELIMINARY ACCURACY ASSESSMENT  (NATIVE V TAME)                ##
###################################################################################

##----------------------- Alison Long - February 2026 -------------------------

# step 7: set up data to assess accuracy of datasets mapping grasslands by comparing 
#         them to ground-truthed point data

# load packages
library(terra)
library(sf)
library(tidyverse)
library(tools)   # for file_path_sans_ext

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
aci_aligned_dir  <- file.path(proc_grass_data_dir, "aafc_aci/aligned")
ngi_aligned_dir  <- file.path(proc_grass_data_dir, "cfga_ngi/aligned")
cgr_masked_dir   <- file.path(proc_grass_data_dir, "cgr_gam_adv/masked")
eccc_masked_dir  <- file.path(proc_grass_data_dir, "eccc_grass/masked")
cnvc_aligned_dir <- file.path(proc_grass_data_dir, "natureserve_cnvc/aligned")
pudl_aligned_dir <- file.path(proc_grass_data_dir, "phjv_pudl/aligned")
pli_aligned_dir  <- file.path(proc_grass_data_dir, "sk_pli/aligned")
gpw_masked_dir   <- file.path(proc_grass_data_dir, "wri_gpw_ann/masked")
plvi_aligned_dir <- file.path(proc_grass_data_dir, "ab_plvi/aligned")
hfi_aligned_dir  <- file.path(proc_grass_data_dir, "abmi_hfi/aligned")
mgi_aligned_dir  <- file.path(proc_grass_data_dir, "mb_mgi/aligned")
pglc_aligned_dir <- file.path(proc_grass_data_dir, "aafc_pglc/aligned")

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
eco_file <- file.path(eco_dir, "ecodistricts_clean_class_sa.shp")

# grass raster files
aci2024_file <- file.path(aci_aligned_dir, "aafc_aci_merged_reproj_algn_natvtam_noshrub.tif")
aci2025_file <- file.path(aci_aligned_dir, "aafc_aci_2025_merged_reproj_algn_natvtam_noshrub.tif")
ngi_file     <- file.path(ngi_aligned_dir, "cfga_ngi_merged_reproj_algn_natvtam_good.tif")
cgr_file     <- file.path(cgr_masked_dir, "cgr_reproj_algn_natvtam_masked.tif")
eccc_file    <- file.path(eccc_masked_dir, "eccc_reclass_reproj_algn_natvtam_masked.tif")
cnvc_file    <- file.path(cnvc_aligned_dir, "cnvc_algn_natvtam.tif")
pudl_file    <- file.path(pudl_aligned_dir, "pudl_reclass_reproj_algn_natvtam_noshrub.tif")
pli_ap_file  <- file.path(pli_aligned_dir, "sk_pli_ap_reproj_algn_natvtam_noshrub.tif")
pli_cu_file  <- file.path(pli_aligned_dir, "sk_pli_cu_reproj_algn_natvtam_noshrub.tif")
pli_mg_file  <- file.path(pli_aligned_dir, "sk_pli_mg_reproj_algn_natvtam_noshrub.tif")
pli_mm_file  <- file.path(pli_aligned_dir, "sk_pli_mm_reproj_algn_natvtam_noshrub.tif")
gpw_file     <- file.path(gpw_masked_dir, "gpw_reproj_algn_natvtam_noshrub_masked.tif")
plvi_file    <- file.path(plvi_aligned_dir, "plvi_reproj_algn_natvtam_gte2015.tif")
hfi_file     <- file.path(hfi_aligned_dir, "hfi_reproj_algn_binary_gte2015.tif") # tame only - same as binary
mgi_file     <- file.path(mgi_aligned_dir, "mb_mgi_reproj_algn_natvtam_noshrub.tif")
pglc_file    <- file.path(pglc_aligned_dir, "pglc_reproj_algn_natvtam.tif")

####--- READ INTO TERRA ---####

#gt <- vect(gt_file)
gt_sf <- sf::st_read(gt_file)

# keep only point geometries and remove empty geometries
#gt_sf <- gt_sf[st_geometry_type(gt_sf) == "POINT", ]
gt_sf <- gt_sf[st_geometry_type(gt_sf) == "POINT" &
                 !sf::st_is_empty(gt_sf),]

# convert to terra
gt <- terra::vect(gt_sf)
gt # 200,985 points

cnvc <- rast(cnvc_file) # for plotting
sap  <- vect(sap_file)  # for plotting

files <- c(
  aci2024_file, aci2025_file, ngi_file, cgr_file, eccc_file,
  cnvc_file, pudl_file, pli_ap_file,
  pli_cu_file, pli_mg_file, pli_mm_file, gpw_file,
  plvi_file, hfi_file, mgi_file, pglc_file
)

keys <- c(
  "aci2024","aci2025","ngi","cgr","eccc","cnvc","pudl","pli_ap",
  "pli_cu","pli_mg","pli_mm","gpw","plvi","hfi","mgi","pglc"
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

# 1) Attributes (already in gt_out after your binding step)
vals <- terra::values(gt_out)               # data.frame with 200,985 rows

# 2) Coordinates for a SpatVector (points) as a data.frame
g   <- terra::geom(gt_out, df = TRUE)       # has columns 'geom','part','x','y','hole'
stopifnot(nrow(g) == nrow(gt_out))

#xy  <- g[, c("x", "y")]

xy <- terra::crds(gt_out, df = TRUE)

dim(vals)
dim(xy)

nrow(gt_out)

# 3) Sanity check and combine
stopifnot(nrow(vals) == nrow(xy))           
#attr_df <- cbind(vals, xy, stringsAsFactors = FALSE)

attr_df <- cbind(
  vals,       # 14 original attributes
  vals_only,  # 27 raster layers
  xy          # coordinates
)

dim(attr_df)

# rename coord columns
names(attr_df)[(ncol(attr_df)-1):ncol(attr_df)] <- c("x", "y")

stopifnot(nrow(attr_df) == nrow(gt_out))
# Inspect NA shares per layer (helpful for coverage/masking)
sort(sapply(attr_df[names(r_all)], function(x) mean(is.na(x))), decreasing = TRUE)[1:10]

# add a field that categorizes the ground-truthed data into native vs tame vs not (only "not" is shrub)

# Normalize and derive the new field
norm <- function(x) tolower(trimws(as.character(x)))
stopifnot("LC_Cl_N" %in% names(attr_df))
stopifnot("Sorc_ID" %in% names(attr_df))
vals <- norm(attr_df$LC_Cl_N)

# ---- 3-class ground truth: 0 = not grass, 1 = tame, 2 = native ----
attr_df$LC_class3 <- NA_integer_

# not grass 
attr_df$LC_class3[vals == "woody"] <- 0L
attr_df$LC_class3[vals == "wetland"] <- 0L
attr_df$LC_class3[vals == "open water"] <- 0L
attr_df$LC_class3[vals == "annual cropland"] <- 0L

# tame grass
attr_df$LC_class3[vals == "tame"] <- 1L

# treat "mixed" as tame
attr_df$LC_class3[vals == "mixed"] <- 1L   

# native grass
attr_df$LC_class3[vals == "native"] <- 2L

# Quick QA
table(attr_df$LC_class3, useNA = "ifany")
#    0      1      2 
# 4111 178693  18181

table(vals, useNA = "ifany")
#vals
#annual cropland           mixed          native      open water            tame         wetland 
#1355            4112           18181              87          174581              22 
#woody 
#2647

head(attr_df)


####--- EXPORT ---####

# write to csv
csv_out_path <- file.path(out_dir, "groundtruth_points_with_grass_raster_values_natvtam_ecodist_pli_sep_17Jul2026.csv")
write.csv(attr_df, csv_out_path, row.names = FALSE)

# filter attr_df to only include rows where Yer_Cll is 2015+
stopifnot("Yer_Cll" %in% names(attr_df))
attr_df_2015 <- attr_df %>%
  filter(Yer_Cll >= 2015)

head(attr_df_2015)

csv_out_path <- file.path(out_dir, "groundtruth_points_with_grass_raster_values_natvtam_pli_sep_ge2015_17Jul2026.csv")
write.csv(attr_df_2015, csv_out_path, row.names = FALSE)

####--- CLEAN FOR NEXT STEP ---####

# filter out any non-grass points (woody and annual cropland)
#attr_df_2015_clean <- attr_df_2015 %>%
#  filter(LC_class3 %in% c(1, 2))

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
      Sorc_ID == "BC_ACI_2025" ~ "aci2024", 
      TRUE ~ ""    # keep unmatched cases empty
    )
  )

# Move new columns to left
attr_df_2015 <- attr_df_2015 %>%
  relocate(LC_class3, Training)

# write cleaned dataset to csv
csv_out_path <- file.path(out_dir, "groundtruth_points_with_grass_raster_values_natvtam_pli_sep_ge2015_17Jul2026_clean_withnongrass.csv")
write.csv(attr_df_2015, csv_out_path, row.names = FALSE)

