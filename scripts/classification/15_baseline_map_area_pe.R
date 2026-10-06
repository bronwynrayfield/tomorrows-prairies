
# calculate grass area of baseline map in Prairies Ecozone for Bronwyn

library(tidyverse)
library(terra)
library(sf)

results_dir   <- "R:/CanadianGrasslands/WG_Classification/results/consensus_mapping"
proc_data_dir <- "R:/CanadianGrasslands/processeddatasets"

m5_file       <- file.path(results_dir, "new_fusion_global_norm_modext_classdifficultynorm_31Aug2026/hard_class_m5_boost_tame_grass.tif")
baseline_file <- file.path(results_dir, "new_fusion_global_norm_modext_classdifficultynorm_31Aug2026/TP_baseline_map_3Sep2026.tif")
pe_file       <- file.path(proc_data_dir, "studyarea/modelling_study_area_aligned_cnvc_Jan2026.tif")

m5       <- rast(m5_file)
baseline <- rast(baseline_file)
pe       <- rast(pe_file)

#---- crop and mask m5 map to Prairies Ecozone ----#
m5_crop <- crop(m5, pe)
m5_mask <- mask(m5_crop, pe)

# verify alignment
compareGeom(m5_mask, pe, stopOnError = TRUE)

plot(m5_mask)
freq(m5_mask)

#   layer value     count
# 1     1     1 373,138,249 = non-grass
# 2     1     2  57,628,305 = tame grass
# 3     1     3  88,259,473 = native grass

#---- crop and mask baseline map to Prairies Ecozone ----#
baseline_crop <- crop(baseline, pe)
baseline_mask <- mask(baseline_crop, pe)

# verify alignment
compareGeom(baseline_mask, pe, stopOnError = TRUE)

plot(baseline_mask)
freq(baseline_mask)

#   layer value     count
# 1     1     1  95,773,537 = native grass
# 2     1     2  57,628,305 = tame grass
# 3     1     3 291,638,177 = annual crop
# 4     1     4  14,422,618 = settlement
# 5     1     5  59,280,901 = other


