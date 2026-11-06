

#################################################################################
##                 RESULTS: CALCULATE AREA OF BASELINE MAP                     ##
#################################################################################

##----------------------- Alison Long - March 2026 -------------------------

# STEP 14: Calculate grass area of baseline map in Prairies Ecozone 

library(tidyverse)
library(terra)
library(sf)

results_dir   <- "R:/CanadianGrasslands/WG_Classification/results/consensus_mapping"
proc_data_dir <- "R:/CanadianGrasslands/processeddatasets"

m5_file       <- file.path(results_dir, "new_fusion_global_norm_modext_classdifficultynorm_31Aug2026/hard_class_m5_boost_tame_grass.tif")
baseline_file <- file.path(results_dir, "new_fusion_global_norm_modext_classdifficultynorm_31Aug2026/TP_baseline_map_3Sep2026.tif")
ent_file      <- file.path(results_dir, "new_fusion_global_norm_modext_classdifficultynorm_31Aug2026/entropy_norm.tif")
acc_file      <- file.path(results_dir, "new_fusion_global_norm_modext_classdifficultynorm_31Aug2026/accuracy/cm2_cd_m5_Balanced_Accuracy.tif")
pe_file       <- file.path(proc_data_dir, "studyarea/modelling_study_area_aligned_cnvc_Jan2026.tif")
ecodist_file  <- file.path(proc_data_dir, "boundaries/national-ecological-framework-for-canada/Ecodistricts/ecodistricts_clean_rast_Oct2026.tif")
prov_file     <- file.path(proc_data_dir, "boundaries/provinces_prairiesecozone_aligned.shp")

m5       <- rast(m5_file)
baseline <- rast(baseline_file)
ent      <- rast(ent_file)
acc      <- rast(acc_file)
pe       <- rast(pe_file)
ecodist  <- rast(ecodist_file)
prov     <- vect(prov_file)

#---- crop and mask m5 map to Prairies Ecozone ----#
m5_crop <- crop(m5, pe)
m5_mask <- mask(m5_crop, pe)

# verify alignment
compareGeom(m5_mask, pe, stopOnError = TRUE)

plot(m5_mask)
freq(m5_mask)

#   layer value     count
# 1     1     1 373,138,249 = non-grass   = 33,582,442 ha
# 2     1     2  57,628,305 = tame grass   = 5,186,547 ha
# 3     1     3  88,259,473 = native grass = 7,943,353 ha

#---- crop and mask baseline map to Prairies Ecozone ----#
baseline_crop <- crop(baseline, pe)
baseline_mask <- mask(baseline_crop, pe)

# verify alignment
compareGeom(baseline_mask, pe, stopOnError = TRUE)

plot(baseline_mask)
freq(baseline_mask)

#   layer value     count
# 1     1     1  95,773,537 = native grass = 8,619,618 ha
# 2     1     2  57,628,305 = tame grass   = 5,186,547 ha
# 3     1     3 291,638,177 = annual crop = 26,247,436 ha
# 4     1     4  14,422,618 = settlement   = 1,298,036 ha
# 5     1     5  59,280,901 = other        = 5,335,281 ha

#---- baseline map area in Prairies Ecozone by province ----#

prov

prov_rast <- rasterize(prov, baseline_mask, field = "PRNAME")
prov_rast
plot(prov_rast)

tab <- crosstab(c(baseline_mask, prov_rast), long = TRUE)
tab

#    which.max       PRNAME         n
# 1          1      Alberta  50,477,543 = 
# 2          1     Manitoba   3,617,036 = 
# 3          1 Saskatchewan  41,518,673 =
# 4          2      Alberta  20,337,767 =
# 5          2     Manitoba   6,278,251 =
# 6          2 Saskatchewan  30,988,838 =
# 7          3      Alberta  81,599,651 =
# 8          3     Manitoba  45,745,576 =
# 9          3 Saskatchewan 164,291,145 =
# 10         4      Alberta   5,507,374 =
# 11         4     Manitoba   2,450,280 =
# 12         4 Saskatchewan   6,464,188 = 
# 13         5      Alberta  15,831,199 =
# 14         5     Manitoba  19,190,054 =
# 15         5 Saskatchewan  24,253,024 =

tab <- tab %>%
  mutate(
    hectares = n * 0.09,
    acres = hectares * 2.47105381)

tab <- tab %>%
  mutate(
    class = recode(
      as.character(which.max),
      "1" = "Native grass",
      "2" = "Tame grass",
      "3" = "Annual crop",
      "4" = "Settlement",
      "5" = "Other"))

tab

#    which.max       PRNAME         n   hectares      acres        class
# 1          1      Alberta  50477543  4542978.9 11225945.2 Native grass
# 2          1     Manitoba   3617036   325533.2   804410.2 Native grass
# 3          1 Saskatchewan  41518673  3736680.6  9233538.8 Native grass
# 4          2      Alberta  20337767  1830399.0  4523014.5   Tame grass
# 5          2     Manitoba   6278251   565042.6  1396250.6   Tame grass
# 6          2 Saskatchewan  30988838  2788995.4  6891757.8   Tame grass
# 7          3      Alberta  81599651  7343968.6 18147341.6  Annual crop
# 8          3     Manitoba  45745576  4117101.8 10173580.2  Annual crop
# 9          3 Saskatchewan 164291145 14786203.0 36537503.4  Annual crop
# 10         4      Alberta   5507374   495663.7  1224811.6   Settlement
# 11         4     Manitoba   2450280   220525.2   544929.6   Settlement
# 12         4 Saskatchewan   6464188   581776.9  1437602.1   Settlement
# 13         5      Alberta  15831199  1424807.9  3520777.0        Other
# 14         5     Manitoba  19190054  1727104.9  4267769.0        Other
# 15         5 Saskatchewan  24253024  2182772.2  5393747.5        Other

# export
tab_df <- as.data.frame(tab)

write.csv(tab_df,
  file.path(results_dir, "new_fusion_global_norm_modext_classdifficultynorm_31Aug2026/baseline_area_by_province_pe.csv"),
  row.names = FALSE)

#------------ another way

# stack class and province rasters
x <- c(baseline_mask, prov_rast)
names(x) <- c("class", "province")

# convert to data.frame and count combinations
df <- as.data.frame(x, na.rm = TRUE)

pixel_counts <- df %>%
  count(province, class, name = "pixels") %>%
  arrange(province, class)

pixel_counts

#        province class    pixels
# 1       Alberta     1  50477543
# 2       Alberta     2  20337767
# 3       Alberta     3  81599651
# 4       Alberta     4   5507374
# 5       Alberta     5  15831199
# 6      Manitoba     1   3617036
# 7      Manitoba     2   6278251
# 8      Manitoba     3  45745576
# 9      Manitoba     4   2450280
# 10     Manitoba     5  19190054
# 11 Saskatchewan     1  41518673
# 12 Saskatchewan     2  30988838
# 13 Saskatchewan     3 164291145
# 14 Saskatchewan     4   6464188
# 15 Saskatchewan     5  24253024

#---- baseline map area in Prairies Ecozone by ecodistrict ----#

ecodist

# crop and mask to Prairies Ecozone
ecodist_crop <- crop(ecodist, pe)
ecodist_mask <- mask(ecodist_crop, pe)

# verify alignment
compareGeom(ecodist_mask, pe, stopOnError = TRUE)

plot(ecodist_mask)

# filter out non-grass and combine grass into single class
baseline_grass <- ifel(baseline_mask %in% c(1, 2), 1, NA)
plot(baseline_grass)
freq(baseline_grass)

# calculate pixels per ecodistrict per class
tab_grass <- crosstab(c(baseline_grass, ecodist_mask), long = TRUE)
tab_grass

tab_grass <- tab_grass %>%
  mutate(
    hectares = n * 0.09,
    acres = hectares * 2.47105381)

tab_grass

# export
tab_grass_df <- as.data.frame(tab_grass)

write.csv(tab_grass_df,
          file.path(results_dir, "new_fusion_global_norm_modext_classdifficultynorm_31Aug2026/baseline_grass_area_by_ecodist_pe.csv"),
          row.names = FALSE)

# plot as raster

# create lookup table: ecodistrict -> grass pixel count
lookup <- tab_grass[, c("ECODISTRIC", "n")]

# convert to reclassification matrix
rcl <- as.matrix(lookup)

# create raster of grass pixel counts by ecodistrict
grass_n_rast <- classify(
  ecodist_mask,
  rcl = rcl,
  others = NA)

plot(
  grass_n_rast,
  col = rev(hcl.colors(100, "YlGn")),
  main = "Grass pixels by ecodistrict")

# by hectares
grass_ha_rast <- subst(
  ecodist_mask,
  from = tab_grass$ECODISTRIC,
  to   = tab_grass$hectares)

plot(
  grass_ha_rast,
  col = rev(hcl.colors(100, "YlGn")),
  main = "Grass area (ha) by ecodistrict")

### normalize by ecodistrict area

# total area of each ecodistrict (all pixels)
tab_total <- freq(ecodist_mask)
tab_total

# join with grass counts
tab_grass_norm <- tab_grass %>%
  select(ECODISTRIC, grass_pixels = n) %>%
  left_join(
    tab_total %>%
      select(ECODISTRIC = value,
             total_pixels = count),
    by = "ECODISTRIC"
  ) %>%
  mutate(
    prop_grass = grass_pixels / total_pixels,
    pct_grass = 100 * prop_grass
  )

# create the raster
grass_pct_rast <- classify(
  ecodist_mask,
  rcl = as.matrix(tab_grass_norm[, c("ECODISTRIC", "pct_grass")]))

grass_pct_rast

plot(
  grass_pct_rast,
  col = rev(hcl.colors(100, "YlGn")),
  main = "Percent grassland by ecodistrict")

plot(
  grass_pct_rast,
  col = rev(hcl.colors(100, "Greens")),
  main = "Percent grassland by ecodistrict")

# library(tidyterra)
# library(ggplot2)
# 
# ggplot() +
#   geom_spatraster(data = grass_pct_rast) +
#   scale_fill_viridis_c(
#     direction = -1,
#     name = "% Grassland"
#   ) +
#   theme_minimal()

# export raster
writeRaster(
  grass_pct_rast,
  filename = file.path(results_dir, "grassland_pct_by_ecodistrict.tif"),
  datatype = "FLT4S",
  gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
  NAflag = 255,
  overwrite = TRUE
)

# export to PNG
png(
  filename = file.path(results_dir, "grassland_pct_by_ecodistrict.png"),
  width = 4000,
  height = 3000,
  res = 400
)

plot(
  grass_pct_rast,
  col = rev(hcl.colors(100, "Greens")),
  main = "Percent grassland by ecodistrict"
)

dev.off()

# export table of values
write.csv(
  tab_grass_norm,
  file.path(results_dir, "percent_grassland_by_ecodistrict.csv"),
  row.names = FALSE)

#---- entropy in Prairies Ecozone by ecodistrict ----#

# inspect
ent
plot(ent)

# crop and mask to Prairies Ecozone
ent_crop <- crop(ent, pe)
ent_mask <- mask(ent_crop, pe)

# verify alignment
compareGeom(ent_mask, pe, stopOnError = TRUE)

plot(ent_mask)

plot(
  ent_mask,
  col = rev(hcl.colors(100, "Grays")),
  main = "Entropy (normalized)")

# calculate mean entropy by ecodistrict
tab_ent <- zonal(ent_mask, ecodist_mask, fun = "mean", na.rm = TRUE)
head(tab_ent)

# create raster of ecodistict mean entropy values
ent_ecodist_rast <- classify(
  ecodist_mask,
  rcl = as.matrix(tab_ent[, c("ECODISTRIC", "entropy_norm")]))

ent_ecodist_rast

plot(
  ent_ecodist_rast,
  col = rev(hcl.colors(100, "Grays")),
  main = "Mean entropy by ecodistrict")

# export raster
writeRaster(
  ent_ecodist_rast,
  filename = file.path(results_dir, "mean_entropy_norm_by_ecodistrict.tif"),
  datatype = "FLT4S",
  gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
  NAflag = 255,
  overwrite = TRUE)

writeRaster(
  ent_mask,
  filename = file.path(results_dir, "entropy_norm_pe.tif"),
  datatype = "FLT4S",
  gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
  NAflag = 255,
  overwrite = TRUE)

# export to PNG
png(
  filename = file.path(results_dir, "mean_entropy_norm_by_ecodistrict.png"),
  width = 4000,
  height = 3000,
  res = 400)

plot(
  grass_pct_rast,
  col = rev(hcl.colors(100, "Grays")),
  main = "Mean entropy by ecodistrict")

dev.off()

# export to PNG
png(
  filename = file.path(results_dir, "entropy_norm_pe.png"),
  width = 4000,
  height = 3000,
  res = 400)

plot(
  ent_mask,
  col = rev(hcl.colors(100, "Grays")),
  main = "Entropy (normalized)")

dev.off()

# export table of values
write.csv(
  tab_ent,
  file.path(results_dir, "mean_entropy_norm_by_ecodistrict.csv"),
  row.names = FALSE)

#---- balanced accuracy in Prairies Ecozone by ecodistrict ----#

# inspect
acc
plot(acc)

# crop and mask to Prairies Ecozone
acc_crop <- crop(acc, pe)
acc_mask <- mask(acc_crop, pe)

# verify alignment
compareGeom(acc_mask, pe, stopOnError = TRUE)

plot(acc_mask)

plot(
  acc_mask,
  col = rev(hcl.colors(100, "Blues")),
  main = "Balanced accuracy")

writeRaster(
  acc_mask,
  filename = file.path(results_dir, "bal_acc_m5_ecodist_pe.tif"),
  datatype = "FLT4S",
  gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
  NAflag = 255,
  overwrite = TRUE)

# export to PNG
png(
  filename = file.path(results_dir, "bal_acc_m5_ecodist_pe.png"),
  width = 4000,
  height = 3000,
  res = 400)

plot(
  acc_mask,
  col = rev(hcl.colors(100, "Blues")),
  main = "Balanced accuracy")

dev.off()

# # extract and export table of values
# levels(acc_mask)
# str(levels(acc_mask))
# is.factor(acc_mask)
# minmax(acc_mask)
# length(unique(values(acc_mask)))
# 
# acc_tab <- zonal(
#   acc_mask,
#   ecodist_mask,
#   fun = "mean",
#   na.rm = TRUE
# )
# 
# names(acc_tab) <- c("ECODISTRIC", "balanced_accuracy")
# acc_tab
# 
# write.csv(
#   acc_tab,
#   file.path(results_dir, "balanced_accuracy_by_ecodistrict.csv"),
#   row.names = FALSE
# )
