
###################################################################################
##                VALIDATION: CHECK ACCURACY OF CONSENSUS MAP                    ##
###################################################################################

##----------------------- Alison Long - July 2026 -------------------------

# STEP 11: Assess accuracy of consensus map by comparing to ground-truthed point 
#          data derived from crop insurance data (not same points used to generate
#          accuracy metrics for input rasters). Generates accuracy metrics for 
#          consensus map and modified versions of the map to choose final version.

# load packages
library(terra)
library(sf)
library(tidyverse)
library(tools)   # for file_path_sans_ext
library(readxl)
# library(rnaturalearth)
# library(viridisLite)
# library(ggspatial)

# Set the temporary directory to a drive with sufficient space
terraOptions(tempdir = "R:/CanadianGrasslands/WG_Classification/scratch")

####--- DIRECTORIES ---####

# pilot area
pa_dir <- "R:/CanadianGrasslands/processeddatasets/pilotarea"

# study area
sa_dir <- "R:/CanadianGrasslands/processeddatasets/studyarea"

# grasslands
raw_grass_data_dir <- "R:/CanadianGrasslands/existingdatasets"
cons_map_dir <- "R:/CanadianGrasslands/WG_Classification/results/consensus_mapping"
out_dir <- "R:/CanadianGrasslands/WG_Classification/processed/groundtruthed_data_points"

# ground-truthed data
gt_dir <- "R:/CanadianGrasslands/processeddatasets/groundtruth_data"

# ecodistricts 
eco_dir <- "R:/CanadianGrasslands/processeddatasets/boundaries/national-ecological-framework-for-canada/Ecodistricts"

# grass raster directories
cm01_dir <- file.path(cons_map_dir, "v4_temporal_modelling_extent_24Jul2026")
cm02_dir <- file.path(cons_map_dir, "new_fusion_global_norm_modext_26Aug2026")
cm03_dir <- file.path(cons_map_dir, "new_fusion_global_norm_modext_classdifficultynorm_31Aug2026")

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
eco_vect_file <- file.path(eco_dir, "ecodistricts_clean_class_sa.shp")
eco_rast_file <- file.path(eco_dir, "ecodistricts_clean_rast.tif")

# grass raster files
cm01_file <- file.path(cm01_dir, "grass_3class.tif")
cm02_file <- file.path(cm02_dir, "hard_class.tif")

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

# ignore error:
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

files <- c(cm01_file, cm02_file)

keys <- c("cm01", "cm02")

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

# crop
e <- terra::ext(rlist$cm02)

rlist <- lapply(rlist, function(r) terra::crop(r, e))

s <- rast(rlist)


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

# add a field that categorizes the ground-truthed data 

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

# fix values for cm02 (0 to 2 instead of 1 to 3)
attr_df <- attr_df %>%
  mutate(cm02 = cm02 - 1)

####--- EXPORT ---####

# write to csv
csv_out_path <- file.path(out_dir, "groundtruth_points_with_consensus_map_raster_values_26Aug2026.csv")
write.csv(attr_df, csv_out_path, row.names = FALSE)

# filter attr_df to only include rows where Yer_Cll is 2015+
stopifnot("Yer_Cll" %in% names(attr_df))
attr_df_2015 <- attr_df %>%
  filter(Yer_Cll >= 2015)

head(attr_df_2015)

csv_out_path <- file.path(out_dir, "groundtruth_points_with_consensus_map_raster_values_ge2015_26Aug2026.csv")
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
#metadata_csv_path <- "R:/CanadianGrasslands/WG_Classification/github/tomorrows-prairies/processeddatasets/Grass/metadata_full_extent_nattam_pli_sep.csv"

out_dir <- "R:/CanadianGrasslands/WG_Classification/results/consensus_mapping"
#dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# crop insurance-derived points directory 
points_dir <- "R:/CanadianGrasslands/existingdatasets/Grassland/groundtruth_points/from_EL/crop_ins_derived"

# ----------------------------
# 3) Configuration
# ----------------------------

# 

# ----------------------------
# 4) Load points check for required columns
# ----------------------------
# points  from Emily 3Sep2026 to validate consensus maps
points <- read_csv(file.path(points_dir, "point_raster_extract_native_tame_other_2025_v3.csv"))
points

names(points)
str(points)
unique(points$Source_ID)
table(points$LC_Class_N, useNA = "ifany")

# NATIVE  OTHER   TAME 
# 47197 165653  42042

# Rename columns
points <- points %>%
  rename(
    Pt_ID_O    = Point_ID,
    Sorc_ID  = Source_ID,
    Yer_Cll  = Year_Colle,
    LC_Cl_N  = LC_Class_N,
    Ecds_ID = ecodist,
    Provinc  = PROV,
    Ecrgn_n  = ecoreg)

# correct values in consensus maps (should be 0-2 not 1-3)
points <- points %>%
  mutate(across(starts_with("cm2"),
      ~ ifelse(is.na(.x), NA, .x - 1)))

# Create column with numerical class for ground-truthed cover type (binary)
points <- points %>%
  mutate(
    LC_grass = case_when(
      LC_Cl_N %in% c("TAME", "NATIVE") ~ 1,
      LC_Cl_N == "OTHER" ~ 0,
      TRUE ~ NA_real_))

unique(points$LC_grass)
table(points$LC_grass, useNA = "ifany")

# Create column with numerical class for ground-truthed cover type (triclass)
points <- points %>%
  mutate(
    LC_tri = case_when(
      LC_Cl_N == "NATIVE" ~ 2,
      LC_Cl_N == "TAME"   ~ 1,
      LC_Cl_N == "OTHER"  ~ 0,
      TRUE ~ NA_real_))

unique(points$LC_tri)
table(points$LC_tri, useNA = "ifany")

# Check for NAs in Pt_ID
points %>%
  summarise(n_na = sum(is.na(Pt_ID_O)))

# Create new Pt_ID column 
points <- points %>%
  mutate(Pt_ID = row_number()) %>%
  relocate(Pt_ID)

# Drop rows where ecodistrict ID is 0
points <- points %>%
  filter(Ecds_ID != 0)

# Fix ecoregions
points <- points %>%
  mutate(
    Ecrgn_n = case_when(
      Ecds_ID == 618 ~ "Western Alberta Upland",
      Ecds_ID == 815 ~ "Mixed Grassland",
      Ecds_ID == 828 ~ "Mixed Grassland",
      Ecds_ID == 833 ~ "Mixed Grassland",
      Ecds_ID == 838 ~ "Cypress Upland",
      TRUE ~ Ecrgn_n
    )
  )

points %>%
  distinct(Ecds_ID, Ecrgn_n) %>%
  count(Ecds_ID) %>%
  filter(n > 1)


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
                     "Interior Transition Ranges") ~ "14.3",
      Ecrgn_n %in% c("Columbia Mountains and Highlands", 
                     "Southern Rocky Mountain Trench",
                     "Northern Continental Divide",
                     "Western Continental Ranges",
                     "Selkirk-Bitterroot Foothills") ~ "14.4"))

unique(points$Ecprv)

unique(points$Ecrgn_n)

# Define required columns
required_cols <- c("Pt_ID", "LC_tri", "LC_grass", "Sorc_ID", "Ecds_ID", 
                   "Ecrgn_n", "Ecprv", "Yer_Cll", "LC_Cl_N")
missing_req <- setdiff(required_cols, names(points))
if (length(missing_req) > 0) {
  stop("Missing required columns in points file: ", paste(missing_req, collapse = ", "))
}

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
table(points$LC_tri, useNA = "ifany")

#----------------------------

# filter to points collected in 2025
pts.2025 <- points %>%
  filter(Yer_Cll == 2025)

#Check LC types present
table(pts.2025$LC_tri, useNA = "ifany")
#      1      2      3 
# 145945  16001  23303

# # plot 
# pts_sf <- st_as_sf(
#   pts.2025,
#   coords = c("POINT_X", "POINT_Y"),
#   crs = 4326
# )
# 
# sap <- vect(sap_file)
# sap_sf <- st_as_sf(sap)
# 
# ggplot() +
#   geom_sf(data = sap_sf,
#           fill = NA,
#           colour = "grey40",
#           linewidth = 0.3) +
#   geom_sf(data = pts_sf,
#           aes(colour = factor(LC_class3)),
#           size = 0.5) +
#   theme_minimal()
# 
# # filter to points collected in 2024 or 2025
# pts.2024.2025 <- points %>%
#   filter(Yer_Cll %in% c(2024, 2025))
# 
# #Check LC types present
# table(pts.2024.2025$LC_class3, useNA = "ifany")
# #   0     1     2 
# # 315 10016   762 but consensus maps are 93% NA
# 
# # plot 
# pts_sf2 <- st_as_sf(
#   pts.2024.2025,
#   coords = c("POINT_X", "POINT_Y"),
#   crs = 4326
# )
# 
# sap <- vect(sap_file)
# sap_sf <- st_as_sf(sap)
# 
# ggplot() +
#   geom_sf(data = sap_sf,
#           fill = NA,
#           colour = "grey40",
#           linewidth = 0.3) +
#   geom_sf(data = pts_sf2,
#           aes(colour = factor(LC_class3)),
#           size = 0.5) +
#   theme_minimal()
# 
# # filter points to those where consensus maps are not NA
# pts.2025.valid <- pts.2025 %>%
#   filter(!is.na(cm01) & !is.na(cm02))
# 
# pts.2024.2025.valid <- pts.2024.2025 %>%
#   filter(!is.na(cm01) & !is.na(cm02))
# 
# #Check LC types present
# table(pts.2025.valid$LC_class3, useNA = "ifany")
# # 0   1   2 
# # 55 320 231 
# 
# table(pts.2024.2025.valid$LC_class3, useNA = "ifany")
# # 0   1   2 
# # 55 320 359

# ----------------------------
# 6) Reshape points to long and attach allowed years (list-column)
# ----------------------------
names(pts.2025)

eval_data <- pts.2025 %>%
  pivot_longer(
    cols = starts_with("cm"),
    names_to = "raster_group",
    values_to = "pred"
  ) %>%
  mutate(
    LC_class3 = as.integer(LC_tri),
    pred      = as.integer(pred)
  )

# rename raster_group column to raster
eval_data <- eval_data %>%
  rename(raster = raster_group)

write.csv(eval_data, file.path(out_dir, "consensus_maps_accuracy_check_2025only_3Sep2026.csv"))

# make sure classes are integers
eval_data <- eval_data %>%
  mutate(
    LC_class3 = as.integer(LC_class3),
    pred      = as.integer(pred)
  ) %>%
  filter(
    LC_class3 %in% c(0, 1, 2),
    pred %in% c(0, 1, 2)
  )

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
  
  #cm <- as.matrix(cm)
  cm <- matrix(as.numeric(cm),
               nrow = nrow(cm),
               ncol = ncol(cm),
               dimnames = dimnames(cm))
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
  
  # n     <- sum(tab)
  # diagv <- diag(tab)
  # rows  <- rowSums(tab)   # support (truth)
  # cols  <- colSums(tab)   # predicted totals
  
  n <- as.numeric(sum(tab))
  diagv <- as.numeric(diag(tab))
  rows <- as.numeric(rowSums(tab))
  cols <- as.numeric(colSums(tab))
  
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
  key_cols = "raster"
)

confusions_by_dataset <- confusion_block_multiclass(
  eval_data,
  key_cols = "raster"
)

metrics_by_dataset

# ----------------------------
# Export output tables
# ----------------------------
out_xlsx <- file.path(out_dir, "consensus_maps_accuracy_check_metrics_2025only_3Sep2026.xlsx")

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


# ---- GROUPED: collapse pli_* -> pli ----
accuracy_by_ecodistrict <- eval_data %>%
  summarize_block_multiclass(
    key_cols = c("raster", "Ecds_ID")
  ) %>%
  arrange(raster, Ecds_ID)


# ---- Save ----
readr::write_csv(accuracy_by_ecodistrict,
                 file.path(out_dir, "accuracy_by_ecodistrict_consensus_maps_3Sep2026.csv"))


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
# GROUPED
# --------------------------------------------------

eval_data %>%
  distinct(Ecds_ID, Ecrgn_n, Ecprv) %>%
  count(Ecds_ID) %>%
  filter(n > 1)

eval_data %>%
  distinct(Ecds_ID, Ecrgn_n, Ecprv) %>%
  filter(Ecds_ID %in% c(618, 815, 828, 833, 838)) %>%
  arrange(Ecds_ID)

#   Ecds_ID Ecrgn_n                Ecprv
# 1     618 Western Alberta Upland 9.1  
# 2     618 Peace Lowland          9.2  
# 3     815 Mixed Grassland        10.3 
# 4     815 Cypress Upland         10.3 
# 5     828 Mixed Grassland        10.3 
# 6     828 Moist Mixed Grassland  10.3 
# 7     833 Mixed Grassland        10.3 
# 8     833 Cypress Upland         10.3 
# 9     838 Cypress Upland         10.3 
# 10     838 Mixed Grassland        10.3 


accuracy_by_ecodistrict2 <- accuracy_by_ecodistrict %>%
  left_join(ecodist_to_regions, by = "Ecds_ID")

fallback_joined <- accuracy_by_ecodistrict2 %>%
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
names(fallback_joined)

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

accuracy_with_fallback <- fallback_joined %>%
  mutate(
    source_level = case_when(
      N_Evaluated_ecoD >= min_n ~ "ecodistrict",
      N_Evaluated_ecoR >= min_n ~ "ecoregion",
      TRUE ~ "ecoprovince"
    )
  )

metric_cols <- setdiff(
  names(accuracy_by_ecodistrict),
  c("raster", "Ecds_ID")
)

for (m in metric_cols) {
  accuracy_with_fallback[[m]] <- case_when(
    accuracy_with_fallback$source_level == "ecodistrict" ~
      accuracy_with_fallback[[paste0(m, "_ecoD")]],
    accuracy_with_fallback$source_level == "ecoregion" ~
      accuracy_with_fallback[[paste0(m, "_ecoR")]],
    TRUE ~
      accuracy_with_fallback[[paste0(m, "_ecoP")]]
  )
}

# ----------------------------
# Export multiclass fallback results
# ----------------------------

readr::write_csv(
  accuracy_with_fallback,
  file.path(
    out_dir,
    "accuracy_by_ecodistrict_nattam_ecoregprovfallback_consensus_maps_3Sep2026.csv"
  )
)

#########################################################################
#########################################################################

# ----------------------------
# 15) Diagnostics
# ----------------------------
table(eval_data$raster)

sum(is.na(eval_data$raster))

fallback_joined %>%
  summarise(
    ecoR_missing =
      sum(is.na(N_evaluated_ecoR)),
    ecoP_missing =
      sum(is.na(N_evaluated_ecoP))
  )

fallback_joined %>%
  filter(
    is.na(N_evaluated_ecoR) |
      is.na(N_evaluated_ecoP)
  ) %>%
  count(raster)

eval_data %>%
  count(raster) %>%
  ggplot(aes(x = reorder(raster, n), y = n)) +
  geom_col() +
  coord_flip() +
  scale_y_continuous(labels = scales::label_comma()) +
  labs(
    x = "Dataset",
    y = "Count of validation points",
    title = "Number of validation points by dataset"
  ) +
  theme_bw()

#------ NA rates for major metrics at each level ------#

fallback_joined %>%
  summarise(
    n_rows = n(),
    
    OA_ecoD_NA = sum(is.na(Overall_Accuracy_ecoD)),
    OA_ecoR_NA = sum(is.na(Overall_Accuracy_ecoR)),
    OA_ecoP_NA = sum(is.na(Overall_Accuracy_ecoP)),
    
    BA_ecoD_NA = sum(is.na(Balanced_Accuracy_ecoD)),
    BA_ecoR_NA = sum(is.na(Balanced_Accuracy_ecoR)),
    BA_ecoP_NA = sum(is.na(Balanced_Accuracy_ecoP)),
    
    Kappa_ecoD_NA = sum(is.na(Kappa_ecoD)),
    Kappa_ecoR_NA = sum(is.na(Kappa_ecoR)),
    Kappa_ecoP_NA = sum(is.na(Kappa_ecoP)),
    
    MCC_ecoD_NA = sum(is.na(MCC_Multiclass_ecoD)),
    MCC_ecoR_NA = sum(is.na(MCC_Multiclass_ecoR)),
    MCC_ecoP_NA = sum(is.na(MCC_Multiclass_ecoP))
  ) %>%
  print(width = Inf)

# n_rows OA_ecoD_NA OA_ecoR_NA OA_ecoP_NA BA_ecoD_NA BA_ecoR_NA BA_ecoP_NA Kappa_ecoD_NA Kappa_ecoR_NA
#   1816        480        264         80        480        264         80           527           288

# Kappa_ecoP_NA MCC_ecoD_NA MCC_ecoR_NA MCC_ecoP_NA
#           104         562         288         104

#------ Final fallback NA rates ------#

accuracy_with_fallback %>%
  summarise(
    n_rows = n(),
    
    OA_NA = sum(is.na(Overall_Accuracy)),
    OA_pct_NA = 100 * mean(is.na(Overall_Accuracy)),
    
    Macro_F1_NA = sum(is.na(Macro_F1)),
    Macro_F1_pct_NA = 100 * mean(is.na(Macro_F1)),
    
    Weighted_F1_NA = sum(is.na(Weighted_F1)),
    Weighted_F1_pct_NA = 100 * mean(is.na(Weighted_F1)),
    
    BA_NA = sum(is.na(Balanced_Accuracy)),
    BA_pct_NA = 100 * mean(is.na(Balanced_Accuracy)),
    
    Kappa_NA = sum(is.na(Kappa)),
    Kappa_pct_NA = 100 * mean(is.na(Kappa)),
    
    MCC_NA = sum(is.na(MCC_Multiclass)),
    MCC_pct_NA = 100 * mean(is.na(MCC_Multiclass))
  ) %>%
  print(width = Inf)

# n_rows OA_NA OA_pct_NA Macro_F1_NA Macro_F1_pct_NA Weighted_F1_NA Weighted_F1_pct_NA BA_NA BA_pct_NA
#   1816    80      4.41          80            4.41              0                  0    80      4.41

# Kappa_NA Kappa_pct_NA MCC_NA MCC_pct_NA
#      104         5.73    129       7.10

accuracy_with_fallback %>%
  filter(is.na(Overall_Accuracy)) %>%
  count(raster)

#   raster        n
# 1 cm1          10
# 2 cm2          10
# 3 cm2_cd       10
# 4 cm2_cd_m1    10
# 5 cm2_cd_m2    10
# 6 cm2_cd_m3    10
# 7 cm2_cd_m4    10
# 8 cm2_cd_m5    10

accuracy_with_fallback %>%
  filter(is.na(Overall_Accuracy)) %>%
  summarise(
    min_n = min(N_Evaluated, na.rm = TRUE),
    max_n = max(N_Evaluated, na.rm = TRUE)
  )

# min_n max_n
#     0     0

#------ Where are MCC/Kappa missing? ------#

accuracy_with_fallback %>%
  count(Truth_Classes_Present)

accuracy_with_fallback %>%
  count(Pred_Classes_Present)

accuracy_with_fallback %>%
  summarise(
    one_truth_class =
      sum(Truth_Classes_Present == 1),
    
    two_truth_classes =
      sum(Truth_Classes_Present == 2),
    
    three_truth_classes =
      sum(Truth_Classes_Present == 3),
    
    one_pred_class =
      sum(Pred_Classes_Present == 1),
    
    two_pred_classes =
      sum(Pred_Classes_Present == 2),
    
    three_pred_classes =
      sum(Pred_Classes_Present == 3),
    
    Kappa_NA =
      sum(is.na(Kappa)),
    
    MCC_NA =
      sum(is.na(MCC_Multiclass))
  ) %>%
  print(width = Inf)

# Which class is driving the 346/451 single-class units?
accuracy_with_fallback %>%
  summarise(
    only_class0 =
      sum(
        Support_0 > 0 &
          Support_1 == 0 &
          Support_2 == 0
      ),
    
    only_class1 =
      sum(
        Support_0 == 0 &
          Support_1 > 0 &
          Support_2 == 0
      ),
    
    only_class2 =
      sum(
        Support_0 == 0 &
          Support_1 == 0 &
          Support_2 > 0
      )
  )

# only_class0 only_class1 only_class2
#          32           0           0

# where are these happening?
accuracy_with_fallback %>%
  filter(
    Support_0 > 0,
    Support_1 == 0,
    Support_2 == 0
  )

# Which class pairings dominate 
accuracy_with_fallback %>%
  mutate(
    combo = case_when(
      Support_0 > 0 & Support_1 > 0 & Support_2 == 0 ~ "0+1",
      Support_0 > 0 & Support_1 == 0 & Support_2 > 0 ~ "0+2",
      Support_0 == 0 & Support_1 > 0 & Support_2 > 0 ~ "1+2",
      TRUE ~ NA_character_
    )
  ) %>%
  count(combo)

#   combo     n
# 1 0+1      64
# 2 NA     1752

#------ Which classes are absent? ------#

accuracy_with_fallback %>%
  summarise(
    no_class0 = sum(Support_0 == 0, na.rm = TRUE),
    no_class1 = sum(Support_1 == 0, na.rm = TRUE),
    no_class2 = sum(Support_2 == 0, na.rm = TRUE)
  )

# no_class0 no_class1 no_class2
#        80       112       176

# ----- Missing one-vs-rest MCCs ------#

accuracy_with_fallback %>%
  summarise(
    MCC0_NA = sum(is.na(MCC_OvR_0)),
    MCC1_NA = sum(is.na(MCC_OvR_1)),
    MCC2_NA = sum(is.na(MCC_OvR_2))
  )

# MCC0_NA MCC1_NA MCC2_NA
#     129     133     346

accuracy_with_fallback %>%
  summarise(
    no_class0_truth = sum(Support_0 == 0),
    no_class1_truth = sum(Support_1 == 0),
    no_class2_truth = sum(Support_2 == 0),
    
    no_class0_pred = sum(Precision_0 %>% is.na()),
    no_class1_pred = sum(Precision_1 %>% is.na()),
    no_class2_pred = sum(Precision_2 %>% is.na())
  )

# no_class0_truth no_class1_truth no_class2_truth no_class0_pred no_class1_pred no_class2_pred
#              80             112             176             80            132            324

#------ Grouped vs atomic results ------#

accuracy_with_fallback %>%
  summarise(
    OA_NA = sum(is.na(Overall_Accuracy)),
    MCC_NA = sum(is.na(MCC_Multiclass))
  )

# OA_NA MCC_NA
#    80    129

#------ Which classes are often missing ------#

accuracy_with_fallback %>%
  count(
    Truth_Classes_Present,
    Pred_Classes_Present
  ) %>%
  arrange(
    Truth_Classes_Present,
    Pred_Classes_Present
  )


ggplot(accuracy_with_fallback,
       aes(
         factor(Truth_Classes_Present),
         factor(Pred_Classes_Present))) +
  geom_count()


#------ Check correlations among metrics ------#

accuracy_with_fallback %>%
  select(
    Overall_Accuracy,
    Balanced_Accuracy,
    Weighted_F1,
    Macro_F1
  ) %>%
  cor(use = "pairwise.complete.obs")

#                   Overall_Accuracy Balanced_Accuracy Weighted_F1  Macro_F1
# Overall_Accuracy         1.0000000         0.1999424   0.7719941 0.4931031
# Balanced_Accuracy        0.1999424         1.0000000   0.1875836 0.5534888
# Weighted_F1              0.7719941         0.1875836   1.0000000 0.3002021
# Macro_F1                 0.4931031         0.5534888   0.3002021 1.0000000


#------ Check distribution of sample sizes --------#

accuracy_with_fallback %>%
  ggplot(aes(x = N_Evaluated)) +
  geom_histogram(
    binwidth = 5,
    color = "black",
    fill = "steelblue"
  ) +
  geom_vline(
    xintercept = 20,
    linetype = "dashed",
    color = "red",
    linewidth = 1
  ) +
  labs(
    title = "Distribution of Validation Sample Sizes",
    x = "n evaluated",
    y = "Count of ecodistricts"
  ) +
  theme_bw()

# log-scaled version to address right skew
accuracy_with_fallback %>%
  ggplot(aes(x = N_Evaluated)) +
  geom_histogram(
    bins = 50,
    color = "black",
    fill = "steelblue"
  ) +
  scale_x_log10() +
  labs(
    title = "Distribution of Validation Sample Sizes",
    x = "n evaluated (log10 scale)",
    y = "Count of ecodistricts"
  ) +
  theme_bw()

#------ Number of rows with n less than 20 ------#

accuracy_with_fallback %>%
  summarise(
    n_below_20 = sum(N_Evaluated < 20),
    pct_below_20 = 100 * mean(N_Evaluated < 20)
  )

# n_below_20 pct_below_20
#        104         5.73

# Examine low-sample rows
accuracy_with_fallback %>%
  filter(N_Evaluated < 20) %>%
  select(
    raster,
    Ecds_ID,
    Ecrgn_n,
    Ecprv,
    source_level,
    N_Evaluated,
    Truth_Classes_Present,
    Pred_Classes_Present,
    Macro_F1,
    Balanced_Accuracy,
    Weighted_F1,
    Overall_Accuracy,
    Support_0,
    Support_1,
    Support_2
  ) %>%
  arrange(N_evaluated, raster, Ecds_ID) %>%
  print(width = Inf, n = Inf)

#------ Examine rows with only 1 or 2 points ------#

accuracy_with_fallback %>%
  filter(N_Evaluated %in% c(1, 2)) %>%
  select(
    raster,
    Ecds_ID,
    Ecrgn_n,
    Ecprv,
    source_level,
    N_Evaluated,
    Truth_Classes_Present,
    Pred_Classes_Present,
    Macro_F1,
    Balanced_Accuracy,
    Weighted_F1,
    Overall_Accuracy
  ) %>%
  arrange(N_Evaluated, raster, Ecds_ID) %>%
  print(width = Inf)

# ----------------------------
# 16) Export file for consensus mapping
# ----------------------------

# Primary metric: Macro F1 
# Fallback metric: Balanced Accuracy 

# Grouped
accuracy_with_fallback_cons <- accuracy_with_fallback %>%
  mutate(
    acc = coalesce(Macro_F1, Balanced_Accuracy),
    acc_source = if_else(
      !is.na(Macro_F1),
      "Macro_F1",
      "Balanced_Accuracy"
    )
  )

# Add quality flag for rows with N less than 20
accuracy_with_fallback_cons <- accuracy_with_fallback_cons %>%
  mutate(
    low_confidence_acc =
      #acc_source == "Overall_Accuracy" &
      N_Evaluated < 20
  )

# Atomic
accuracy_with_fallback_atomic_cons <- accuracy_with_fallback_atomic %>%
  mutate(
    acc = coalesce(Macro_F1, Balanced_Accuracy),
    acc_source = if_else(
      !is.na(Macro_F1),
      "Macro_F1",
      "Balanced_Accuracy"
    )
  )

# Add quality flag for rows with N less than 20
accuracy_with_fallback_atomic_cons <- accuracy_with_fallback_atomic_cons %>%
  mutate(
    low_confidence_acc =
      #acc_source == "Overall_Accuracy" &
      N_Evaluated < 20
  )

#---- Export ----

# GROUPED
readr::write_csv(
  accuracy_with_fallback_cons,
  file.path(out_dir, "nattam_accuracy_combinedCIandCM_forconsensus_29Jul2026.csv")
)

# ATOMIC
readr::write_csv(
  accuracy_with_fallback_atomic_cons,
  file.path(out_dir, "nattam_accuracy_atomic_combinedCIandCM_forconsensus_29Jul2026.csv")
)

#------ How often is fallback used? ------#

accuracy_with_fallback_cons %>%
  count(acc_source)

accuracy_with_fallback_cons %>%
  filter(acc_source == "Balanced_Accuracy") %>%
  summarise(
    n = n(),
    pct = 100 * n() / nrow(accuracy_with_fallback_cons)
  )

#------ Why is Macro_F1 undefined? ------#

accuracy_with_fallback %>%
  filter(is.na(Macro_F1)) %>%
  summarise(
    n_rows = n(),
    one_truth_class = sum(Truth_Classes_Present == 1),
    two_truth_classes = sum(Truth_Classes_Present == 2),
    three_truth_classes = sum(Truth_Classes_Present == 3),
    
    one_pred_class = sum(Pred_Classes_Present == 1),
    two_pred_classes = sum(Pred_Classes_Present == 2),
    three_pred_classes = sum(Pred_Classes_Present == 3)
  ) %>%
  print(width = Inf)

#------ How many fallback rows have extreme BA values? ------#

accuracy_with_fallback_cons %>%
  filter(
    acc_source == "Balanced_Accuracy",
    acc %in% c(0, 1)
  ) %>%
  count(acc)

#------ Are extreme fallback values concentrated in small n? ------#

accuracy_with_fallback_cons %>%
  filter(
    acc_source == "Balanced_Accuracy",
    acc %in% c(0, 1)
  ) %>%
  summarise(
    n = n(),
    min_n = min(N_Evaluated),
    median_n = median(N_Evaluated),
    max_n = max(N_Evaluated)
  )

#------ Distribution of n among fallback rows ------#

accuracy_with_fallback_cons %>%
  filter(
    acc_source == "Balanced_Accuracy",
    acc %in% c(0, 1)
  ) %>%
  count(N_Evaluated, sort = TRUE) %>%
  print(n = Inf)

#------ Inspect all fallback rows ------#

accuracy_with_fallback_cons %>%
  filter(acc_source == "Balanced_Accuracy") %>%
  select(
    raster,
    Ecds_ID,
    Ecrgn_n,
    Ecprv,
    source_level,
    N_Evaluated,
    Truth_Classes_Present,
    Pred_Classes_Present,
    acc,
    acc_source,
    Macro_F1,
    Balanced_Accuracy,
    Overall_Accuracy,
    Weighted_F1,
    Support_0,
    Support_1,
    Support_2,
    F1_0,
    F1_1,
    F1_2
  ) %>%
  arrange(
    desc(N_Evaluated),
    raster,
    Ecds_ID
  ) %>%
  print(n = Inf, width = Inf)

accuracy_with_fallback_cons %>%
  filter(acc_source == "Balanced_Accuracy") %>%
  summarise(
    no_class0 = sum(Support_0 == 0),
    no_class1 = sum(Support_1 == 0),
    no_class2 = sum(Support_2 == 0)
  )

#------ CLASS-SPECIFIC METRICS ------#

# Summary by class
class_metrics_long <- accuracy_with_fallback %>%
  select(
    raster,
    Ecds_ID,
    N_Evaluated,
    starts_with("Precision_"),
    starts_with("Recall_"),
    starts_with("F1_"),
    starts_with("MCC_OvR_")
  ) %>%
  pivot_longer(
    cols = -c(raster, Ecds_ID, N_Evaluated),
    names_to = c("metric", "class"),
    names_pattern = "(.*)_(\\d)$"
  )

class_metrics_long %>%
  group_by(class, metric) %>%
  summarise(
    n = sum(!is.na(value)),
    NAs = sum(is.na(value)),
    mean = mean(value, na.rm = TRUE),
    median = median(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    q25 = quantile(value, 0.25, na.rm = TRUE),
    q75 = quantile(value, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

# Metric distributions among classes
class_metrics_long %>%
  filter(metric == "F1") %>%
  ggplot(aes(class, value)) +
  geom_boxplot() +
  labs(
    x = "Class",
    y = "F1 score",
    title = "Distribution of F1 by class"
  ) +
  theme_bw()

class_metrics_long %>%
  filter(metric == "Recall") %>%
  ggplot(aes(class, value)) +
  geom_boxplot() +
  labs(
    x = "Class",
    y = "Recall",
    title = "Distribution of recall by class"
  ) +
  theme_bw()

class_metrics_long %>%
  filter(metric == "Precision") %>%
  ggplot(aes(class, value)) +
  geom_boxplot() +
  labs(
    x = "Class",
    y = "Precision",
    title = "Distribution of precision by class"
  ) +
  theme_bw()

# How much support does each class have?
accuracy_with_fallback %>%
  summarise(
    median_support0 = median(Support_0),
    median_support1 = median(Support_1),
    median_support2 = median(Support_2),
    
    mean_support0 = mean(Support_0),
    mean_support1 = mean(Support_1),
    mean_support2 = mean(Support_2)
  )

accuracy_with_fallback %>%
  select(
    Support_0,
    Support_1,
    Support_2
  ) %>%
  pivot_longer(
    everything(),
    names_to = "class",
    values_to = "support"
  ) %>%
  ggplot(aes(support)) +
  geom_histogram(bins = 40) +
  facet_wrap(~ class, scales = "free_y") +
  theme_bw()

# How often does each class drive Macro F1 failures?
accuracy_with_fallback %>%
  filter(is.na(Macro_F1)) %>%
  summarise(
    no_class0 = sum(Support_0 == 0),
    no_class1 = sum(Support_1 == 0),
    no_class2 = sum(Support_2 == 0)
  )

accuracy_with_fallback %>%
  filter(is.na(Macro_F1)) %>%
  mutate(
    missing_pattern = paste0(
      ifelse(Support_0 == 0, "0", ""),
      ifelse(Support_1 == 0, "1", ""),
      ifelse(Support_2 == 0, "2", "")
    )
  ) %>%
  count(missing_pattern, sort = TRUE)

# Which rasters struggle with native grassland?

# Average native-class F1 by source
accuracy_with_fallback %>%
  group_by(raster) %>%
  summarise(
    mean_F1_native = mean(F1_2, na.rm = TRUE),
    median_F1_native = median(F1_2, na.rm = TRUE),
    n_valid = sum(!is.na(F1_2))
  ) %>%
  arrange(mean_F1_native)

# Average tame-class F1 by source
accuracy_with_fallback %>%
  group_by(raster) %>%
  summarise(
    mean_F1_tame = mean(F1_1, na.rm = TRUE),
    median_F1_tame = median(F1_1, na.rm = TRUE),
    n_valid = sum(!is.na(F1_1))
  ) %>%
  arrange(mean_F1_tame)

# Identify systematic omission vs commission errors
# Recall < Precision = class is being omitted
# Precision < Recall = class is being over-predicted
accuracy_with_fallback %>%
  summarise(
    native_precision = mean(Precision_2, na.rm = TRUE),
    native_recall = mean(Recall_2, na.rm = TRUE),
    
    tame_precision = mean(Precision_1, na.rm = TRUE),
    tame_recall = mean(Recall_1, na.rm = TRUE),
    
    not_precision = mean(Precision_0, na.rm = TRUE),
    not_recall = mean(Recall_0, na.rm = TRUE)
  )

# Relationship between sample size and class-specific performance
accuracy_with_fallback %>%
  ggplot(
    aes(
      N_Evaluated,
      F1_2
    )
  ) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "loess") +
  scale_x_log10() +
  labs(
    x = "N evaluated",
    y = "F1 Native",
    title = "Native F1 versus sample size"
  ) +
  theme_bw()

# Quantify whether native grasslands are intrinsically harder
# If class 2 remains substantially lower even after excluding NA values and accounting 
# for sample size, that would support a conclusion that native grassland is genuinely 
# the most difficult class for the underlying datasets.
accuracy_with_fallback %>%
  summarise(
    mean_F1_0 = mean(F1_0, na.rm = TRUE),
    mean_F1_1 = mean(F1_1, na.rm = TRUE),
    mean_F1_2 = mean(F1_2, na.rm = TRUE),
    
    median_F1_0 = median(F1_0, na.rm = TRUE),
    median_F1_1 = median(F1_1, na.rm = TRUE),
    median_F1_2 = median(F1_2, na.rm = TRUE)
  )

# Sample size plot
accuracy_with_fallback %>%
  summarise(
    cor_native =
      cor(log10(N_Evaluated), F1_2,
          use = "complete.obs"),
    
    cor_tame =
      cor(log10(N_Evaluated), F1_1,
          use = "complete.obs"),
    
    cor_not =
      cor(log10(N_Evaluated), F1_0,
          use = "complete.obs")
  )

#---- Class-specific metrics (class-specific F1 or fallbacks if NA) ----
accuracy_with_fallback_cons2 <- accuracy_with_fallback %>%
  mutate(
    native_metric = coalesce(F1_2, Macro_F1, Balanced_Accuracy),
    tame_metric   = coalesce(F1_1, Macro_F1, Balanced_Accuracy),
    not_metric    = coalesce(F1_0, Macro_F1, Balanced_Accuracy))

accuracy_with_fallback_atomic_cons2 <- accuracy_with_fallback_atomic %>%
  mutate(
    native_metric = coalesce(F1_2, Macro_F1, Balanced_Accuracy),
    tame_metric   = coalesce(F1_1, Macro_F1, Balanced_Accuracy),
    not_metric    = coalesce(F1_0, Macro_F1, Balanced_Accuracy))

#---- Export ----

# GROUPED
readr::write_csv(
  accuracy_with_fallback_cons2,
  file.path(out_dir, "nattam_accuracy_combinedCIandCM_forconsensus_29Jul2026_final.csv")
)

# ATOMIC
readr::write_csv(
  accuracy_with_fallback_atomic_cons2,
  file.path(out_dir, "nattam_accuracy_atomic_combinedCIandCM_forconsensus_29Jul2026_final.csv")
)

#---- Diagnostics for class-specific metrics with fallbacks ----#

# Summarize distributions 
# Does one class systematically have higher/lower weights?
accuracy_with_fallback_cons2 %>%
  summarise(
    mean_native = mean(native_metric, na.rm = TRUE),
    mean_tame   = mean(tame_metric, na.rm = TRUE),
    mean_not    = mean(not_metric, na.rm = TRUE),
    
    median_native = median(native_metric, na.rm = TRUE),
    median_tame   = median(tame_metric, na.rm = TRUE),
    median_not    = median(not_metric, na.rm = TRUE)
  )

accuracy_with_fallback_cons2 %>%
  select(native_metric, tame_metric, not_metric) %>%
  pivot_longer(everything(),
               names_to = "class",
               values_to = "acc") %>%
  ggplot(aes(class, acc)) +
  geom_boxplot() +
  theme_bw()

accuracy_with_fallback_cons2 %>%
  group_by(raster) %>%
  summarise(
    native = mean(native_metric, na.rm = TRUE),
    tame   = mean(tame_metric, na.rm = TRUE),
    not    = mean(not_metric, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(
    c(native, tame, not),
    names_to = "class",
    values_to = "metric"
  ) %>%
  mutate(
    class = factor(
      class,
      levels = c("not", "tame", "native")
    )
  ) %>%
  ggplot(aes(
    x = reorder(raster, metric),
    y = metric,
    fill = class
  )) +
  geom_col(position = "dodge") +
  coord_flip() +
  scale_fill_manual(
    values = c(
      native = "forestgreen",
      tame   = "goldenrod",
      not    = "grey50"
    ),
    breaks = c("native", "tame", "not"),
    #guide = guide_legend(reverse = TRUE)
  ) +
  labs(
    x = "Dataset",
    y = "Mean accuracy metric",
    fill = "Class",
    title = "Class-specific accuracy metrics by dataset"
  ) +
  theme_bw()

# with number of data points

# Counts per raster
count_df <- eval_data %>%
  count(raster, name = "n") %>%
  rename(raster = raster)

# Data for plotting
plot_df <- accuracy_with_fallback_cons2 %>%
  group_by(raster) %>%
  summarise(
    native = mean(native_metric, na.rm = TRUE),
    tame   = mean(tame_metric, na.rm = TRUE),
    not    = mean(not_metric, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(count_df, by = "raster") %>%
  pivot_longer(
    c(native, tame, not),
    names_to = "class",
    values_to = "metric"
  ) %>%
  mutate(
    class = factor(class,
                   levels = c("not", "tame", "native")
    )
  )

# One label position per raster
label_df <- plot_df %>%
  group_by(raster, n) %>%
  summarise(y = max(metric, na.rm = TRUE) + 0.06, .groups = "drop")

p <- ggplot(plot_df, aes(
  x = reorder(raster, metric),
  y = metric,
  fill = class
)) +
  geom_col(position = position_dodge()) +
  geom_text(
    data = label_df,
    aes(x = raster, y = y, 
        label = paste0("n=",scales::comma(n))),
    inherit.aes = FALSE,
    size = 4
  ) +
  coord_flip() +
  scale_fill_manual(
    values = c(
      native = "forestgreen",
      tame   = "goldenrod",
      not    = "grey50"
    ),
    breaks = c("native", "tame", "not")
  ) +
  labs(
    x = "Dataset",
    y = "Mean accuracy metric",
    fill = "Class",
    title = "Class-specific accuracy metrics by dataset"
  ) +
  theme_bw()

p

ggsave(file.path(out_dir, "accuracy_by_dataset_triclass.png"),
       plot = p,
       width = 10,
       height = 8,
       dpi = 600
)

# How often are fallbacks used?
accuracy_with_fallback %>%
  mutate(
    native_source = case_when(
      !is.na(F1_2) ~ "F1_2",
      !is.na(Macro_F1) ~ "Macro_F1",
      TRUE ~ "Balanced_Accuracy"
    ),
    tame_source = case_when(
      !is.na(F1_1) ~ "F1_1",
      !is.na(Macro_F1) ~ "Macro_F1",
      TRUE ~ "Balanced_Accuracy"
    ),
    not_source = case_when(
      !is.na(F1_0) ~ "F1_0",
      !is.na(Macro_F1) ~ "Macro_F1",
      TRUE ~ "Balanced_Accuracy"
    )
  ) %>%
  summarise(
    native_F1 = sum(native_source == "F1_2"),
    native_macro = sum(native_source == "Macro_F1"),
    native_BA = sum(native_source == "Balanced_Accuracy"),
    
    tame_F1 = sum(tame_source == "F1_1"),
    tame_macro = sum(tame_source == "Macro_F1"),
    tame_BA = sum(tame_source == "Balanced_Accuracy"),
    
    not_F1 = sum(not_source == "F1_0"),
    not_macro = sum(not_source == "Macro_F1"),
    not_BA = sum(not_source == "Balanced_Accuracy")
  )

# Where is fallback occurring?
accuracy_with_fallback %>%
  filter(is.na(F1_2)) %>%
  count(raster, sort = TRUE)

accuracy_with_fallback %>%
  filter(is.na(F1_2)) %>%
  count(Ecrgn_n, sort = TRUE) %>%
  print(n = Inf)

# Compare class-specific metrics against macro F1
accuracy_with_fallback_cons2 %>%
  summarise(
    cor_native =
      cor(native_metric,
          Macro_F1,
          use = "complete.obs"),
    
    cor_tame =
      cor(tame_metric,
          Macro_F1,
          use = "complete.obs"),
    
    cor_not =
      cor(not_metric,
          Macro_F1,
          use = "complete.obs")
  )

# How often do native and tame disagree?
accuracy_with_fallback_cons2 %>%
  summarise(
    mean_diff =
      mean(native_metric - tame_metric,
           na.rm = TRUE),
    
    median_diff =
      median(native_metric - tame_metric,
             na.rm = TRUE)
  )

accuracy_with_fallback_cons2 %>%
  ggplot(
    aes(native_metric, tame_metric)
  ) +
  geom_point(alpha = 0.3) +
  geom_abline(
    slope = 1,
    intercept = 0,
    colour = "red"
  ) +
  theme_bw()

# Relationship with sample size 
accuracy_with_fallback_cons2 %>%
  summarise(
    native =
      cor(log10(N_Evaluated),
          native_metric,
          use = "complete.obs"),
    
    tame =
      cor(log10(N_Evaluated),
          tame_metric,
          use = "complete.obs"),
    
    not =
      cor(log10(N_Evaluated),
          not_metric,
          use = "complete.obs")
  )

# Do datasets tend to be conservative about predicting native grass?
# If native metrics are frequently higher despite lower recall, the
# consensus weighting scheme could favor datasets that are very
# selective about calling pixels native
accuracy_with_fallback_cons2 %>%
  summarise(
    native_gt_tame =
      mean(native_metric > tame_metric,
           na.rm = TRUE),
    
    tame_gt_native =
      mean(tame_metric > native_metric,
           na.rm = TRUE)
  )

#----- Ecodistrict diagnostics ----

# Which ecodistricts consistently have low accuracy?
accuracy_with_fallback_cons2 %>%
  group_by(Ecds_ID) %>%
  summarise(
    native = mean(native_metric, na.rm = TRUE),
    tame   = mean(tame_metric, na.rm = TRUE),
    not    = mean(not_metric, na.rm = TRUE),
    n_rows = n(),
    .groups = "drop"
  ) %>%
  arrange(native) %>%
  print(n = 40)

accuracy_with_fallback_cons2 %>%
  group_by(Ecds_ID) %>%
  summarise(
    native = mean(native_metric, na.rm = TRUE),
    tame   = mean(tame_metric, na.rm = TRUE),
    not    = mean(not_metric, na.rm = TRUE),
    n_rows = n(),
    .groups = "drop"
  ) %>%
  arrange(Ecds_ID) %>%
  print(n = Inf)

# Which ecodistricts show the greatest disagreement among datasets?
accuracy_with_fallback_cons2 %>%
  group_by(Ecds_ID) %>%
  summarise(
    native_sd = sd(native_metric, na.rm = TRUE),
    tame_sd   = sd(tame_metric, na.rm = TRUE),
    not_sd    = sd(not_metric, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(native_sd)) %>%
  print(n = 20)

accuracy_with_fallback_cons2 %>%
  group_by(Ecds_ID) %>%
  summarise(
    native_sd = sd(native_metric, na.rm = TRUE),
    tame_sd   = sd(tame_metric, na.rm = TRUE),
    not_sd    = sd(not_metric, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Ecds_ID) %>%
  print(n = Inf)

# Are accuracy estimates dominated by fallback values?
accuracy_with_fallback %>%
  mutate(
    native_is_fallback = is.na(F1_2),
    tame_is_fallback   = is.na(F1_1),
    not_is_fallback    = is.na(F1_0)
  ) %>%
  group_by(Ecds_ID) %>%
  summarise(
    native_fallback_pct =
      mean(native_is_fallback) * 100,
    tame_fallback_pct =
      mean(tame_is_fallback) * 100,
    not_fallback_pct =
      mean(not_is_fallback) * 100) %>%
  print(n = Inf)

# Ecodistricts with poor class representation
accuracy_with_fallback %>%
  group_by(Ecds_ID) %>%
  summarise(
    mean_classes =
      mean(Truth_Classes_Present)) %>%
  arrange(mean_classes) %>%
  print(n = Inf)

# Support by ecodistrict 
accuracy_with_fallback %>%
  group_by(Ecds_ID) %>%
  summarise(
    support0 = sum(Support_0),
    support1 = sum(Support_1),
    support2 = sum(Support_2)
  ) %>%
  print(n = Inf)

# Native vs tame by ecodistrict 
accuracy_with_fallback_cons2 %>%
  group_by(Ecds_ID) %>%
  summarise(
    mean_diff =
      mean(native_metric - tame_metric,
           na.rm = TRUE)) %>%
  print(n = Inf)

# Relationship to sample size by ecodistrict 
accuracy_with_fallback_cons2 %>%
  group_by(Ecds_ID) %>%
  summarise(mean_acc = mean(native_metric, na.rm = TRUE),
            mean_n = mean(N_Evaluated)) %>%
  print(n = Inf)

ecodistrict_summary <- accuracy_with_fallback_cons2 %>%
  group_by(Ecds_ID) %>%
  summarise(
    mean_acc = mean(native_metric, na.rm = TRUE),
    mean_n = mean(N_Evaluated),
    .groups = "drop"
  )

ggplot(
  ecodistrict_summary,
  aes(mean_n, mean_acc)
) +
  geom_point() +
  scale_x_log10() +
  theme_bw() +
  labs(
    x = "Mean N evaluated (log10 scale)",
    y = "Mean native metric",
    title = "Native metric versus sample size by ecodistrict"
  )

# Which ecodistricts drive native omission?
accuracy_with_fallback_cons2 %>%
  group_by(Ecds_ID) %>%
  summarise(native_precision = mean(Precision_2, na.rm=TRUE),
            native_recall = mean(Recall_2, na.rm=TRUE)) %>%
  mutate(omission_gap = native_precision - native_recall) %>%
  arrange(omission_gap) %>%
  print(n = Inf)

# Spatial autocorrelation
# Join accuracy table to ecodistrict boundaries
# Map:
# native_metric
# tame_metric
# not_metric
# native–tame difference
# fallback frequency

# ----------------------------
# 17) Create rasters of accuracy by ecodistrict
# ----------------------------

#----- create accuracy metrics rasters -----#

# load accuracy data if needed
cm03_file <- file.path(cm03_dir, "accuracy/accuracy_by_ecodistrict_nattam_ecoregprovfallback_consensus_maps_3Sep2026_for_raster.xlsx") 
acc <- read_xlsx(cm03_file, sheet = "clean")

metric_cols <- c("Balanced_Accuracy", "MCC_Multiclass", "F1_0", "F1_1", "F1_2")

acc <- acc %>%
  mutate(
    across(all_of(metric_cols), as.numeric),
    Ecds_ID = as.integer(Ecds_ID)
  )

eco_rast <- rast(eco_file)
eco_rast
plot(eco_rast)

# set up output directory
out_dir <- "R:/CanadianGrasslands/WG_Classification/results/consensus_mapping/new_fusion_global_norm_modext_classdifficultynorm_31Aug2026/accuracy"

#define function
raster_names <- unique(acc$raster)

for (r_name in raster_names) {

  cat("Processing:", r_name, "\n")

  acc_sub <- acc %>%
    filter(raster == r_name)

  # create one raster per metric
  for (metric in metric_cols) {

    lookup <- acc_sub[, c("Ecds_ID", metric)]

    out_rast <- classify(
      eco_rast,
      as.matrix(lookup),
      others = NA
    )

    writeRaster(
      out_rast,
      file.path(
        out_dir,
        paste0(r_name, "_", metric, ".tif")
      ),
      overwrite = TRUE
    )
  }
}

plot(out_rast)

#----- create accuracy metrics shapefiles -----#

# read ecodistrict polygons
eco_vect <- vect(eco_vect_file)

# ensure matching data types
eco_vect$ECODISTRIC <- as.integer(as.character(eco_vect$ECODISTRIC))

# run function
raster_names <- unique(acc$raster)

for (r_name in raster_names) {
  
  cat("Processing:", r_name, "\n")
  
  acc_sub <- acc %>%
    filter(raster == r_name)
  
  for (metric in metric_cols) {
    
    # make a fresh copy of polygons
    out_vect <- eco_vect
    
    # join metric values
    out_vect[[metric]] <- acc_sub[[metric]][
      match(out_vect$ECODISTRIC, acc_sub$Ecds_ID)
    ]
    
    # export
    writeVector(
      out_vect,
      file.path(
        out_dir,
        paste0(r_name, "_", metric, ".shp")
      ),
      overwrite = TRUE
    )
  }
}

plot(out_vect)

#----- create source_level categorical raster -----#
template <- eco_rast

# data check
cats(eco_rast)
#range(values(eco_rast), na.rm = TRUE)
range(acc$Ecds_ID, na.rm = TRUE)

# create source_level lookup table
source_df <- acc %>%
  select(Ecds_ID, source_level) %>%
  distinct()

# check that each ecodistrict has only one source level
source_df %>%
  count(Ecds_ID) %>%
  filter(n > 1)

# convert categories to integer codes
source_levels <- sort(unique(source_df$source_level))

source_lut <- data.frame(
  value = seq_along(source_levels),
  source_level = source_levels
)

source_df <- source_df %>%
  left_join(source_lut, by = "source_level")

# create the new raster
rcl <- as.matrix(
  source_df[, c("Ecds_ID", "value")]
)

source_rast <- classify(
  eco_rast,
  rcl = rcl,
  others = NA
)

# attach category labels
#levels(source_rast) <- source_lut

#visualize
plot(source_rast)

#export
writeRaster(source_rast, file.path(out_dir, "source_level.tif"), overwrite = TRUE)

#----- create source_level categorical shapefile -----#

# read ecodistrict polygons
eco_vect <- vect(eco_vect_file)

# create source_level lookup table
source_df <- acc %>%
  select(Ecds_ID, source_level) %>%
  distinct()

# check that each ecodistrict has only one source level
source_df %>%
  count(Ecds_ID) %>%
  filter(n > 1)

# make sure join fields have same type
source_df$Ecds_ID <- as.integer(source_df$Ecds_ID)
eco_vect$ECODISTRIC <- as.integer(as.character(eco_vect$ECODISTRIC))

# join source_level to polygons
eco_vect$source_level <- source_df$source_level[
  match(eco_vect$ECODISTRIC, source_df$Ecds_ID)
]

# visualize
plot(eco_vect, "source_level")

# export shapefile
writeVector(eco_vect, file.path(out_dir, "source_level.shp"), overwrite = TRUE)

#----- create N_Evaluated raster -----#

# one record per ecodistrict
n_eval_df <- acc %>%
  select(Ecds_ID, N_Evaluated) %>%
  distinct()

# check for duplicate values per ecodistrict
n_eval_df %>%
  count(Ecds_ID) %>%
  filter(n > 1)

# define lookup
rcl <- as.matrix(
  n_eval_df[, c("Ecds_ID", "N_Evaluated")]
)

n_eval_rast <- classify(
  eco_rast,
  rcl,
  others = NA
)

# make numeric
n_eval_rast_num <- as.numeric(n_eval_rast)
plot(n_eval_rast_num)

freq(eco_rast)[1:10, ]

# inspect
plot(as.numeric(n_eval_rast))
plot(n_eval_rast)
n_eval_rast_num
n_eval_rast

plot(n_eval_df)

# export
writeRaster(n_eval_rast, file.path(out_dir, "N_evaluated.tif"), overwrite = TRUE)

#----- create N_Evaluated shapefile -----#

# read ecodistrict polygons
eco_vect <- vect(eco_vect_file)

# create lookup table
n_eval_df <- acc %>%
  select(Ecds_ID, N_Evaluated) %>%
  distinct()

# make sure IDs have same type
eco_vect$ECODISTRIC <- as.integer(as.character(eco_vect$ECODISTRIC))

# join N_Evaluated to polygons
eco_vect$N_Evaluated <- n_eval_df$N_Evaluated[
  match(eco_vect$ECODISTRIC, n_eval_df$Ecds_ID)
]

# inspect
summary(eco_vect$N_Evaluated)

plot(eco_vect, "N_Evaluated")

# export
writeVector(eco_vect, file.path(out_dir, "N_Evaluated.shp"), overwrite = TRUE)

# acc <- acc %>%
#   left_join(source_lookup, by = "source_level")
# 
# eco_vect <- vect("R:/CanadianGrasslands/processeddatasets/boundaries/national-ecological-framework-for-canada/Ecodistricts/ecodistricts_clean_class_sa.shp")
# eco_vect$Ecds_ID <- as.character(eco_vect$ECODISTRIC)
# 
# weights_ecodist_mcc$Ecds_ID <- as.character(weights_ecodist_mcc$Ecds_ID)
# 
# template <- rast("R:/CanadianGrasslands/WG_Classification/processed/aafc_aci/aligned/aafc_aci_merged_reproj_algn_binary_noshrub.tif")
# 
# raster_names <- unique(weights_ecodist_mcc$raster)
# 
# mcc_rasters <- vector("list", length(raster_names))
# names(mcc_rasters) <- raster_names
# 
# for (r_name in raster_names) {
#   
#   w_df <- weights_ecodist_mcc %>% filter(raster == r_name)
#   
#   eco_joined <- merge(eco_vect, w_df, by = "Ecds_ID", all.x = TRUE)
#   
#   fallback_val <- global_mcc$MCC_global[global_mcc$raster == r_name]
#   fallback_val <- (fallback_val + 1) / 2
#   
#   eco_joined$MCC_final[is.na(eco_joined$MCC_final)] <- fallback_val
#   
#   r_mcc <- rasterize(
#     eco_joined,
#     template,
#     field = "MCC_final",
#     touches = TRUE
#   )
#   
#   names(r_mcc) <- r_name
#   
#   mcc_rasters[[r_name]] <- r_mcc
# }
# 
# mcc_stack <- rast(mcc_rasters)
# 
# writeRaster(
#   mcc_stack,
#   file.path(out_dir, "mcc_weight_stack_triclass.tif"),
#   overwrite = TRUE,
#   wopt = list(datatype = "FLT4S",
#               gdal = c("COMPRESS=ZSTD", "TILED=YES", "BIGTIFF=YES"))
# )

# summary_by_ecodistrict <- summarize_block_multiclass(
#   eval_data,
#   key_cols = c("raster", "Ecds_ID")
# )
# 
# # ---- global MCC fallback (FIXED) ----
# global_mcc <- metrics_by_dataset %>%
#   select(raster, MCC_global = MCC_Multiclass) %>%
#   rename(raster = raster)
# 
# # ---- weights ----
# weights_ecodist_mcc <- summary_by_ecodistrict %>%
#   select(raster, Ecds_ID, MCC_Multiclass, N_Evaluated) %>%
#   
#   rename(
#     MCC = MCC_Multiclass,
#     N_Evaluated = N_Evaluated
#   ) %>%
#   
#   mutate(
#     MCC = ifelse(N_Evaluated < 30, NA, MCC)
#   ) %>%
#   
#   left_join(global_mcc, by = "raster") %>%
#   
#   mutate(
#     MCC_final = case_when(
#       is.na(MCC) ~ MCC_global,
#       TRUE ~ 0.7 * MCC + 0.3 * MCC_global
#     ),
#     
#     MCC_final = (MCC_final + 1) / 2,
#     MCC_final = pmin(pmax(MCC_final, 0), 1)
#   ) %>%
#   
#   select(raster, Ecds_ID, MCC_final)

# ----------------------------
# 15) Plots: Map by ecodistrict
# ----------------------------

# ---- Match CRS to template raster ----
ecodistricts <- terra::project(eco_vect, template)

# ---- Convert to sf for ggplot ----
ecodistricts_sf <- sf::st_as_sf(ecodistricts)

# ---- Ensure join key is character ----
ecodistricts_sf <- ecodistricts_sf %>%
  mutate(Ecds_ID = as.character(Ecds_ID))

accuracy_by_ecodistrict <- accuracy_by_ecodistrict %>%
  mutate(Ecds_ID = as.character(Ecds_ID))

# ----------------------------
# Join spatial + accuracy data
# ----------------------------
map_data <- ecodistricts_sf %>%
  left_join(accuracy_by_ecodistrict, by = "Ecds_ID")

# ----------------------------
# Metrics to map
# ----------------------------
metrics_to_map <- c(
  "Overall_Accuracy",
  "Balanced_Accuracy",
  "Kappa",
  "MCC_Multiclass",
  "N_Evaluated"
)

# ----------------------------
# Function to make maps
# ----------------------------
make_map <- function(data, raster_name, metric_name) {
  
  df <- data %>%
    filter(raster == raster_name)
  
  # Define fill scale 
  if (metric_name == "N_Evaluated") {
    fill_scale <- scale_fill_viridis_c(
      na.value = "grey90",
      option = "C",
      name = "N evaluated"
    )
  } else {
    fill_scale <- scale_fill_viridis_c(
      na.value = "grey90",
      option = "C",
      labels = scales::percent_format(accuracy = 1),
      name = metric_name
    )
  }
  
  # Build plot
  p <- ggplot(df) +
    geom_sf(aes(fill = .data[[metric_name]]), color = NA) +
    fill_scale +
    ggtitle(paste0("Raster: ", raster_name, " | ", metric_name)) +
    theme_minimal()
  
  return(p)
}

# ----------------------------
# Generate and save maps
# ----------------------------
unique_rasters <- sort(unique(map_data$raster))

for (r in unique_rasters) {
  
  for (m in metrics_to_map) {
    
    p <- make_map(map_data, r, m)
    
    outfile <- file.path(
      out_dir,
      paste0("map_", r, "_", m, ".png")
    )
    
    ggsave(
      filename = outfile,
      plot = p,
      width = 8,
      height = 6,
      dpi = 300
    )
  }
}

# ------ Mean N_Evaluated per ecodistrict (across rasters)------

# ---- Compute mean N across rasters ----
mean_N_by_ecodistrict <- accuracy_by_ecodistrict_atomic %>%
  group_by(Ecds_ID) %>%
  summarise(
    mean_N_Evaluated = mean(N_Evaluated, na.rm = TRUE),
    n_rasters = sum(!is.na(N_Evaluated)),
    .groups = "drop"
  )

# ---- Join to spatial data ----
mean_N_by_ecodistrict <- mean_N_by_ecodistrict %>%
  mutate(Ecds_ID = as.character(Ecds_ID))

mean_N_map <- ecodistricts_sf %>%
  left_join(mean_N_by_ecodistrict, by = "Ecds_ID")

# --- Plot map ------------
p_mean_N <- ggplot(mean_N_map) +
  geom_sf(aes(fill = mean_N_Evaluated), color = NA) +
  scale_fill_viridis_c(
    option = "C",
    na.value = "grey90",
    name = "Mean N evaluated"
  ) +
  ggtitle("Mean Sample Size (N_Evaluated) Across Rasters by Ecodistrict") +
  theme_minimal()

# ---- Save ----
ggsave(
  filename = file.path(out_dir, "map_mean_N_Evaluated_ecodistrict.png"),
  plot = p_mean_N,
  width = 8,
  height = 6,
  dpi = 300
)

# ----------------------------
# 16) Plots: Overall accuracy (GROUPED and ATOMIC)
# ----------------------------

# grouped summary 
summary_by_raster <- metrics_by_dataset

# atomic summary
summary_by_raster_atomic <- summarize_block_multiclass(
  eval_data,
  key_cols = "raster"
)

# GROUPED plot
p <- ggplot(summary_by_raster,
                    aes(x = reorder(raster, Overall_Accuracy), y = Overall_Accuracy)) +
  geom_col(fill = "#2d7f5e") +
  geom_text(aes(label = paste0("n = ", N_Evaluated)),
            hjust = 1.1,
            color = "white",
            size = 3.5) +
  coord_flip() +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.02))) +
  labs(title = "Overall Accuracy by raster",
       x = "Raster", y = "Overall Accuracy") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(out_dir, "raster_accuracy_plot_nattam_cropsins.png"),
       p, width = 7, height = 5, dpi = 200)

# ATOMIC plot
p_atomic <- ggplot(summary_by_raster_atomic,
                   aes(x = reorder(raster, Overall_Accuracy), y = Overall_Accuracy)) +
  geom_col(fill = "#2d7f5e") +
  geom_text(aes(label = paste0("n = ", N_Evaluated)),
            hjust = 1.1,
            color = "white",
            size = 3.5) +
  coord_flip() +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.02))) +
  labs(title = "Overall Accuracy by Raster",
       x = "Raster", y = "Overall Accuracy") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(out_dir, "raster_accuracy_plot_atomic_nattam_cropins.png"),
       p_atomic, width = 7, height = 5, dpi = 200)

#################################################################
#################################################################

# Compare macro metrics
summary_long <- summary_by_raster %>%
  select(raster,
         Overall_Accuracy,
         Balanced_Accuracy,
         Macro_F1,
         MCC_Multiclass) %>%
  pivot_longer(-raster, names_to = "metric", values_to = "value")

ggplot(summary_long,
       aes(x = reorder(raster, value), y = value, fill = metric)) +
  geom_col(position = "dodge") +
  coord_flip() +
  theme_minimal()
