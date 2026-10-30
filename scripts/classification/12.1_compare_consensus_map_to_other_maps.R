
###################################################################################
##             VALIDATION: COMPARE CONSENSUS MAP TO OTHER PRODUCTS               ##
###################################################################################

##----------------------- Alison Long - July 2026 -------------------------

# STEP 12: Assess accuracy of consensus map and modified versions by comparing 
#          to existing grassland maps

library(tidyverse)
library(terra)
library(sf)

results_dir   <- "R:/CanadianGrasslands/WG_Classification/results/consensus_mapping"
proc_data_dir <- "R:/CanadianGrasslands/WG_Classification/processed"

hard_class_file    <- file.path(results_dir, "new_fusion_global_norm_modext_classdifficultynorm_31Aug2026/hard_class.tif")
hard_class_m1_file <- file.path(results_dir, "new_fusion_global_norm_modext_classdifficultynorm_31Aug2026/hard_class_m1_mild_nongrass_penalty.tif")
hard_class_m2_file <- file.path(results_dir, "new_fusion_global_norm_modext_classdifficultynorm_31Aug2026/hard_class_m2_moderate_nongrass_penalty.tif")
hard_class_m3_file <- file.path(results_dir, "new_fusion_global_norm_modext_classdifficultynorm_31Aug2026/hard_class_m3_strong_nongrass_penalty.tif")
hard_class_m4_file <- file.path(results_dir, "new_fusion_global_norm_modext_classdifficultynorm_31Aug2026/hard_class_m4_boost_all_grass.tif")
hard_class_m5_file <- file.path(results_dir, "new_fusion_global_norm_modext_classdifficultynorm_31Aug2026/hard_class_m5_boost_tame_grass.tif")

aci_file  <- file.path(proc_data_dir, "aafc_aci/aligned/aafc_aci_2025_merged_reproj_algn_natvtam_noshrub.tif")
luts_file <- file.path(proc_data_dir, "aafc_luts_2020_30m_v5/masked/luts30_reproj_algn_binary_masked.tif")
pglc_file <- file.path(proc_data_dir, "aafc_pglc/aligned/pglc_reproj_algn_natvtam.tif")
ngi_file  <- file.path(proc_data_dir, "cfga_ngi/aligned/cfga_ngi_merged_reproj_algn_natvtam_good.tif")
gpw_file  <- file.path(proc_data_dir, "wri_gpw_ann/masked/gpw_reproj_algn_natvtam_noshrub_masked.tif")
eccc_file <- file.path(proc_data_dir, "eccc_grass/masked/eccc_reclass_reproj_algn_natvtam_masked.tif")

hard_class    <- rast(hard_class_file)
hard_class_m1 <- rast(hard_class_m1_file)
hard_class_m2 <- rast(hard_class_m1_file)
hard_class_m3 <- rast(hard_class_m1_file)
hard_class_m4 <- rast(hard_class_m1_file)
hard_class_m5 <- rast(hard_class_m1_file)

aci  <- rast(aci_file)
luts <- rast(luts_file)
pglc <- rast(pglc_file)
ngi  <- rast(ngi_file)
gpw  <- rast(gpw_file)
eccc <- rast(eccc_file)

hard_class
aci

freq(hard_class)
#   layer value     count
# 1     1     1 524,751,637 = 76.51% non
# 2     1     2  63,179,602 =  9.21% tame
# 3     1     3  97,916,599 = 14.28% native

#------- compare ACI and consensus map -------#

# reclassify ACI to match
aci_reclass <- classify(aci, rcl = matrix(c(
                              0, 1, # Non-grass
                              1, 2, # Tame grassland
                              2, 3), # Native grass
                              ncol = 2, byrow = TRUE))

# crop ACI to extent of consensus map
aci_reclass_crop <- crop(aci_reclass, hard_class)
aci_reclass_mask <- mask(aci_reclass_crop, hard_class)

# verify alignment
compareGeom(aci_reclass_mask, hard_class, stopOnError = TRUE)

plot(aci_reclass_mask)
freq(aci_reclass_mask)
#   layer value     count
# 1     1     1 531,272,978 = 77.52% non-grass
# 2     1     2  49,737,125 =  7.26% tame
# 3     1     3 104,317,609 = 15.22% native

# define comparison
comparison_aci <- ((aci_reclass_mask * 10) + hard_class)
names(comparison_aci) <- "aci_vs_consensus"

# 11 = LULC non-grass, hard_class non-grass
# 12 = LULC non-grass, hard_class tame
# 13 = LULC non-grass, hard_class native
# 
# 21 = LULC tame,      hard_class non-grass
# 22 = LULC tame,      hard_class tame
# 23 = LULC tame,      hard_class native
# 
# 31 = LULC native,    hard_class non-grass
# 32 = LULC native,    hard_class tame
# 33 = LULC native,    hard_class native

levels(comparison_aci) <- data.frame(
  value = c(11,12,13,21,22,23,31,32,33),
  class = c(
    "both_non_grass",
    "ACI_non_CM_tame",
    "ACI_non_CM_native",
    "ACI_tame_CM_non",
    "both_tame",
    "ACI_tame_CM_native",
    "ACI_native_CM_non",
    "ACI_native_CM_tame",
    "both_native"
  )
)

# export
writeRaster(comparison_aci, file.path(results_dir, "aci_vs_hardclass.tif"),
  overwrite = TRUE,
  wopt = list(
    datatype = "INT2U",
    gdal = c("COMPRESS=ZSTD", "TILED=YES", "BIGTIFF=YES")))

comparison_aci_export <- comparison_aci

levels(comparison_aci_export) <- NULL

writeRaster(
  comparison_aci_export,
  file.path(results_dir, "aci_vs_hardclass_numeric.tif"),
  overwrite = TRUE,
  wopt = list(
    datatype = "INT2U",
    gdal = c(
      "COMPRESS=ZSTD",
      "TILED=YES",
      "BIGTIFF=YES"
    )
  )
)

# examine
plot(comparison_aci)
freq(comparison_aci)
#   layer              value       count
# 1     1     both_non_grass 503,611,127
# 2     1    ACI_non_CM_tame  18,432,970
# 3     1  ACI_non_CM_native   9,228,881
# 4     1    ACI_tame_CM_non  11,384,341
# 5     1          both_tame  33,787,437
# 6     1 ACI_tame_CM_native   4,565,347
# 7     1  ACI_native_CM_non   9,418,765
# 8     1 ACI_native_CM_tame  10,940,853
# 9     1        both_native  83,957,991

#------- compare LUTS and consensus map -------#

# BUT can I do it? LUTS classes are wrong

# reclassify LUTS to match
luts_reclass <- classify(luts, rcl = matrix(c(
  0, 1, # Non-grass
  1, 2, # Tame grassland
  2, 3), # Native grass
  ncol = 2, byrow = TRUE))

# crop LUTS to extent of consensus map
luts_reclass_crop <- crop(luts_reclass, hard_class)
luts_reclass_mask <- mask(luts_reclass_crop, hard_class)

# verify alignment
compareGeom(luts_reclass_mask, hard_class, stopOnError = TRUE)

plot(luts_reclass_mask)
freq(luts_reclass_mask)

# define comparison
comparison_luts <- ((luts_reclass_mask * 10) + hard_class)
names(comparison_luts) <- "luts_vs_consensus"

# 11 = LULC non-grass, hard_class non-grass
# 12 = LULC non-grass, hard_class tame
# 13 = LULC non-grass, hard_class native
# 
# 21 = LULC tame,      hard_class non-grass
# 22 = LULC tame,      hard_class tame
# 23 = LULC tame,      hard_class native
# 
# 31 = LULC native,    hard_class non-grass
# 32 = LULC native,    hard_class tame
# 33 = LULC native,    hard_class native

levels(comparison_luts) <- data.frame(
  value = c(11,12,13,21,22,23,31,32,33),
  class = c(
    "both_non_grass",
    "luts_non_CM_tame",
    "luts_non_CM_native",
    "luts_tame_CM_non",
    "both_tame",
    "luts_tame_CM_native",
    "luts_native_CM_non",
    "luts_native_CM_tame",
    "both_native"
  )
)

# export
writeRaster(comparison_luts, file.path(results_dir, "luts_vs_hardclass.tif"),
            overwrite = TRUE,
            wopt = list(
              datatype = "INT2U",
              gdal = c("COMPRESS=ZSTD", "TILED=YES", "BIGTIFF=YES")))

comparison_luts_export <- comparison_luts

levels(comparison_luts_export) <- NULL

writeRaster(
  comparison_luts_export,
  file.path(results_dir, "luts_vs_hardclass_numeric.tif"),
  overwrite = TRUE,
  wopt = list(
    datatype = "INT2U",
    gdal = c(
      "COMPRESS=ZSTD",
      "TILED=YES",
      "BIGTIFF=YES"
    )
  )
)

# examine
plot(comparison_luts)
freq(comparison_luts)

#------- compare PGLC and consensus map -------#

pglc

# reclassify PGLC to match
pglc_reclass <- classify(pglc, rcl = matrix(c(
  0, 1, # Non-grass
  1, 2, # Tame grassland
  2, 3), # Native grass
  ncol = 2, byrow = TRUE))

# crop PGLC to extent of consensus map
pglc_reclass_crop <- crop(pglc_reclass, hard_class)
pglc_reclass_mask <- mask(pglc_reclass_crop, hard_class)

# verify alignment
compareGeom(pglc_reclass_mask, hard_class, stopOnError = TRUE)

plot(pglc_reclass_mask)
freq(pglc_reclass_mask)
#   layer value     count
# 1     1     1 179,307,864 = 63.71%
# 2     1     2  50,508,298 = 17.95%
# 3     1     3  51,634,502 = 18.35%

# define comparison
comparison_pglc <- ((pglc_reclass_mask * 10) + hard_class)
names(comparison_pglc) <- "pglc_vs_consensus"

# 11 = LULC non-grass, hard_class non-grass
# 12 = LULC non-grass, hard_class tame
# 13 = LULC non-grass, hard_class native
# 
# 21 = LULC tame,      hard_class non-grass
# 22 = LULC tame,      hard_class tame
# 23 = LULC tame,      hard_class native
# 
# 31 = LULC native,    hard_class non-grass
# 32 = LULC native,    hard_class tame
# 33 = LULC native,    hard_class native

levels(comparison_pglc) <- data.frame(
  value = c(11,12,13,21,22,23,31,32,33),
  class = c(
    "both_non_grass",
    "pglc_non_CM_tame",
    "pglc_non_CM_native",
    "pglc_tame_CM_non",
    "both_tame",
    "pglc_tame_CM_native",
    "pglc_native_CM_non",
    "pglc_native_CM_tame",
    "both_native"
  )
)

# export
writeRaster(comparison_pglc, file.path(results_dir, "pglc_vs_hardclass.tif"),
            overwrite = TRUE,
            wopt = list(
              datatype = "INT2U",
              gdal = c("COMPRESS=ZSTD", "TILED=YES", "BIGTIFF=YES")))

comparison_pglc_export <- comparison_pglc

levels(comparison_pglc_export) <- NULL

writeRaster(
  comparison_pglc_export,
  file.path(results_dir, "pglc_vs_hardclass_numeric.tif"),
  overwrite = TRUE,
  wopt = list(
    datatype = "INT2U",
    gdal = c(
      "COMPRESS=ZSTD",
      "TILED=YES",
      "BIGTIFF=YES")))

# examine
plot(comparison_pglc)
freq(comparison_pglc)
#   layer               value     count
# 1     1      both_non_grass 166,140,458
# 2     1    pglc_non_CM_tame   4,657,798
# 3     1  pglc_non_CM_native   8,509,608
# 4     1    pglc_tame_CM_non   9,462,169
# 5     1           both_tame  15,045,319
# 6     1 pglc_tame_CM_native  26,000,810
# 7     1  pglc_native_CM_non   2,818,416
# 8     1 pglc_native_CM_tame   2,153,106
# 9     1         both_native  46,662,980

#------- compare NGI and consensus map -------#

# reclassify NGI to match
ngi_reclass <- classify(ngi, rcl = matrix(c(
  0, 1, # Non-grass
  1, 2, # Tame grassland
  2, 3), # Native grass
  ncol = 2, byrow = TRUE))

# crop NGI to extent of consensus map
ngi_reclass_crop <- crop(ngi_reclass, hard_class)
ngi_reclass_mask <- mask(ngi_reclass_crop, hard_class)

# verify alignment
compareGeom(ngi_reclass_mask, hard_class, stopOnError = TRUE)

plot(ngi_reclass_mask)
freq(ngi_reclass_mask)
#   layer value       count
# 1     1     1 310,074,752 = 58.34% non-grass
# 2     1     2 124,425,184 = 23.41% tame grass
# 3     1     3  97,023,669 = 18.25% native grass

# define comparison
comparison_ngi <- ((ngi_reclass_mask * 10) + hard_class)
names(comparison_ngi) <- "ngi_vs_consensus"

# 11 = LULC non-grass, hard_class non-grass
# 12 = LULC non-grass, hard_class tame
# 13 = LULC non-grass, hard_class native
# 
# 21 = LULC tame,      hard_class non-grass
# 22 = LULC tame,      hard_class tame
# 23 = LULC tame,      hard_class native
# 
# 31 = LULC native,    hard_class non-grass
# 32 = LULC native,    hard_class tame
# 33 = LULC native,    hard_class native

levels(comparison_ngi) <- data.frame(
  value = c(11,12,13,21,22,23,31,32,33),
  class = c(
    "both_non_grass",
    "ngi_non_CM_tame",
    "ngi_non_CM_native",
    "ngi_tame_CM_non",
    "both_tame",
    "ngi_tame_CM_native",
    "ngi_native_CM_non",
    "ngi_native_CM_tame",
    "both_native"
  )
)

# export
writeRaster(comparison_ngi, file.path(results_dir, "ngi_vs_hardclass.tif"),
            overwrite = TRUE,
            wopt = list(
              datatype = "INT2U",
              gdal = c("COMPRESS=ZSTD", "TILED=YES", "BIGTIFF=YES")))

comparison_ngi_export <- comparison_ngi

levels(comparison_ngi_export) <- NULL

writeRaster(
  comparison_ngi_export,
  file.path(results_dir, "ngi_vs_hardclass_numeric.tif"),
  overwrite = TRUE,
  wopt = list(
    datatype = "INT2U",
    gdal = c("COMPRESS=ZSTD", "TILED=YES", "BIGTIFF=YES")))

# examine
plot(comparison_ngi)
freq(comparison_ngi)
#   layer              value       count
# 1     1     both_non_grass 300,899,025
# 2     1    ngi_non_CM_tame   5,765,740
# 3     1  ngi_non_CM_native   3,409,987
# 4     1    ngi_tame_CM_non  66,339,778
# 5     1          both_tame  40,836,978
# 6     1 ngi_tame_CM_native  17,248,428
# 7     1  ngi_native_CM_non  20,971,170
# 8     1 ngi_native_CM_tame   4,273,727
# 9     1        both_native  71,778,772

#------- compare GPW and consensus map -------#

# reclassify GPW to match
gpw_reclass <- classify(gpw, rcl = matrix(c(
  0, 1, # Non-grass
  1, 2, # Tame grassland
  2, 3), # Native grass
  ncol = 2, byrow = TRUE))

# crop GPW to extent of consensus map
gpw_reclass_crop <- crop(gpw_reclass, hard_class)
gpw_reclass_mask <- mask(gpw_reclass_crop, hard_class)

# verify alignment
compareGeom(gpw_reclass_mask, hard_class, stopOnError = TRUE)

plot(gpw_reclass_mask)
freq(gpw_reclass_mask)
#   layer value     count
# 1     1     1 401,695,078 = 64.49%
# 2     1     2  79,645,851 = 12.79%
# 3     1     3 141,494,500 = 22.72%

# define comparison
comparison_gpw <- ((gpw_reclass_mask * 10) + hard_class)
names(comparison_gpw) <- "gpw_vs_consensus"

# 11 = LULC non-grass, hard_class non-grass
# 12 = LULC non-grass, hard_class tame
# 13 = LULC non-grass, hard_class native
# 
# 21 = LULC tame,      hard_class non-grass
# 22 = LULC tame,      hard_class tame
# 23 = LULC tame,      hard_class native
# 
# 31 = LULC native,    hard_class non-grass
# 32 = LULC native,    hard_class tame
# 33 = LULC native,    hard_class native

levels(comparison_gpw) <- data.frame(
  value = c(11,12,13,21,22,23,31,32,33),
  class = c(
    "both_non_grass",
    "gpw_non_CM_tame",
    "gpw_non_CM_native",
    "gpw_tame_CM_non",
    "both_tame",
    "gpw_tame_CM_native",
    "gpw_native_CM_non",
    "gpw_native_CM_tame",
    "both_native"
  )
)

# export
writeRaster(comparison_gpw, file.path(results_dir, "gpw_vs_hardclass.tif"),
            overwrite = TRUE,
            wopt = list(
              datatype = "INT2U",
              gdal = c("COMPRESS=ZSTD", "TILED=YES", "BIGTIFF=YES")))

comparison_gpw_export <- comparison_gpw

levels(comparison_gpw_export) <- NULL

writeRaster(
  comparison_gpw_export,
  file.path(results_dir, "gpw_vs_hardclass_numeric.tif"),
  overwrite = TRUE,
  wopt = list(
    datatype = "INT2U",
    gdal = c("COMPRESS=ZSTD", "TILED=YES", "BIGTIFF=YES")))

# examine
plot(comparison_gpw)
freq(comparison_gpw)
#   layer              value       count
# 1     1     both_non_grass 399,116,938
# 2     1    gpw_non_CM_tame   1,258,095
# 3     1  gpw_non_CM_native   1,320,045
# 4     1    gpw_tame_CM_non  32,050,821
# 5     1          both_tame  42,256,072
# 6     1 gpw_tame_CM_native   5,338,958
# 7     1  gpw_native_CM_non  30,652,133
# 8     1 gpw_native_CM_tame  19,608,576
# 9     1        both_native  91,233,791

#------- compare ECCC and consensus map -------#

# reclassify ECCC to match
eccc_reclass <- classify(eccc, rcl = matrix(c(
  0, 1, # Non-grass
  1, 2, # Tame grassland
  2, 3), # Native grass
  ncol = 2, byrow = TRUE))

# crop ECCC to extent of consensus map
eccc_reclass_crop <- crop(eccc_reclass, hard_class)
eccc_reclass_mask <- mask(eccc_reclass_crop, hard_class)

# verify alignment
compareGeom(eccc_reclass_mask, hard_class, stopOnError = TRUE)

plot(eccc_reclass_mask)
freq(eccc_reclass_mask)
# layer value     count
# 1     1     1 512,288,174 = 75.17% non
# 2     1     2  97,398,954 = 14.29% tame
# 3     1     3  71,853,680 = 10.54% native

# define comparison
comparison_eccc <- ((eccc_reclass_mask * 10) + hard_class)
names(comparison_eccc) <- "eccc_vs_consensus"

# 11 = LULC non-grass, hard_class non-grass
# 12 = LULC non-grass, hard_class tame
# 13 = LULC non-grass, hard_class native
# 
# 21 = LULC tame,      hard_class non-grass
# 22 = LULC tame,      hard_class tame
# 23 = LULC tame,      hard_class native
# 
# 31 = LULC native,    hard_class non-grass
# 32 = LULC native,    hard_class tame
# 33 = LULC native,    hard_class native

levels(comparison_eccc) <- data.frame(
  value = c(11,12,13,21,22,23,31,32,33),
  class = c(
    "both_non_grass",
    "eccc_non_CM_tame",
    "eccc_non_CM_native",
    "eccc_tame_CM_non",
    "both_tame",
    "eccc_tame_CM_native",
    "eccc_native_CM_non",
    "eccc_native_CM_tame",
    "both_native"
  )
)

# export
writeRaster(comparison_eccc, file.path(results_dir, "eccc_vs_hardclass.tif"),
            overwrite = TRUE,
            wopt = list(
              datatype = "INT2U",
              gdal = c("COMPRESS=ZSTD", "TILED=YES", "BIGTIFF=YES")))

comparison_eccc_export <- comparison_eccc

levels(comparison_eccc_export) <- NULL

writeRaster(
  comparison_eccc_export,
  file.path(results_dir, "eccc_vs_hardclass_numeric.tif"),
  overwrite = TRUE,
  wopt = list(
    datatype = "INT2U",
    gdal = c("COMPRESS=ZSTD", "TILED=YES", "BIGTIFF=YES")))

# examine
plot(comparison_eccc)
freq(comparison_eccc)
# layer               value     count
# 1     1      both_non_grass 494,766,365
# 2     1    eccc_non_CM_tame   6,131,268
# 3     1  eccc_non_CM_native  11,390,541
# 4     1    eccc_tame_CM_non  24,202,807
# 5     1           both_tame  54,982,834
# 6     1 eccc_tame_CM_native  18,213,313
# 7     1  eccc_native_CM_non   1,859,990
# 8     1 eccc_native_CM_tame   2,053,615
# 9     1         both_native  67,940,075

#------ grassland retention in consensus map ------#

# what percentage of each grass dataset is included in the consensus layer? 

# focused on grassland inclusion rather than exact class agreement:
# grass_retention_pct: percentage of all grassland pixels (tame + native) in the source 
#                      dataset that are still grassland in the consensus.
# tame_retention_pct: percentage of source tame pixels that remain grassland (either 
#                     tame or native) in the consensus.
# native_retention_pct: percentage of source native pixels that remain grassland (either 
#                      tame or native) in the consensus.

#define function
calc_retention <- function(comp_raster, dataset_name){
  
  f <- freq(comp_raster)
  
  counts <- setNames(
    f$count,
    as.character(f$value)
  )
  
  # helper
  get_count <- function(x){
    ifelse(x %in% names(counts), counts[x], 0)
  }
  
  tame_non    <- get_count(paste0(dataset_name, "_tame_CM_non"))
  tame_native <- get_count(paste0(dataset_name, "_tame_CM_native"))
  
  native_non  <- get_count(paste0(dataset_name, "_native_CM_non"))
  native_tame <- get_count(paste0(dataset_name, "_native_CM_tame"))
  
  both_tame   <- get_count("both_tame")
  both_native <- get_count("both_native")
  
  # totals in source dataset
  tame_total <-
    tame_non +
    both_tame +
    tame_native
  
  native_total <-
    native_non +
    native_tame +
    both_native
  
  grass_total <-
    tame_total +
    native_total
  
  # retained as grass in consensus
  tame_retained <-
    both_tame +
    tame_native
  
  native_retained <-
    native_tame +
    both_native
  
  grass_retained <-
    tame_retained +
    native_retained
  
  data.frame(
    dataset = dataset_name,
    
    grass_total = grass_total,
    grass_retained = grass_retained,
    grass_retention_pct =
      100 * grass_retained / grass_total,
    
    tame_total = tame_total,
    tame_retained = tame_retained,
    tame_retention_pct =
      100 * tame_retained / tame_total,
    
    native_total = native_total,
    native_retained = native_retained,
    native_retention_pct =
      100 * native_retained / native_total
  )
}

# run
retention_summary <- rbind(
  calc_retention(comparison_aci,  "ACI"),
  calc_retention(comparison_pglc, "pglc"),
  calc_retention(comparison_ngi,  "ngi"),
  calc_retention(comparison_gpw,  "gpw"),
  calc_retention(comparison_eccc, "eccc")
)

retention_summary


# cleaner output
retention_summary <- retention_summary |>
  transform(
    grass_retention_pct  = round(grass_retention_pct, 1),
    tame_retention_pct   = round(tame_retention_pct, 1),
    native_retention_pct = round(native_retention_pct, 1)
  )

retention_summary

#   dataset grass_total grass_retained grass_retention_pct tame_total tame_retained tame_retention_pct
# 1     ACI   154054734      133251628                86.5   49737125      38352784               77.1
# 2    pglc   102142800       89862215                88.0   50508298      41046129               81.3
# 3     ngi   221448853      134137905                60.6  124425184      58085406               46.7
# 4     gpw   221140351      158437397                71.6   79645851      47595030               59.8
# 5    eccc   169252634      143189837                84.6   97398954      73196147               75.2

#   native_total native_retained native_retention_pct
# 1    104317609        94898844                 91.0
# 2     51634502        48816086                 94.5
# 3     97023669        76052499                 78.4
# 4    141494500       110842367                 78.3
# 5     71853680        69993690                 97.4

write.csv(
  retention_summary,
  file.path(results_dir, "grassland_retention_summary.csv"),
  row.names = FALSE
)

######################################################################################
######################################################################################

# same analysis but with chosen baseline map (modification 5)

#------- compare ACI and consensus map mod 5 -------#

# reclassify ACI to match
aci_reclass <- classify(aci, rcl = matrix(c(
  0, 1, # Non-grass
  1, 2, # Tame grassland
  2, 3), # Native grass
  ncol = 2, byrow = TRUE))

# crop ACI to extent of consensus map
aci_reclass_crop <- crop(aci_reclass, hard_class)
aci_reclass_mask <- mask(aci_reclass_crop, hard_class)

# verify alignment
compareGeom(aci_reclass_mask, hard_class, stopOnError = TRUE)

plot(aci_reclass_mask)
freq(aci_reclass_mask)
#   layer value     count
# 1     1     1 531,272,978 = 77.52% non-grass
# 2     1     2  49,737,125 =  7.26% tame
# 3     1     3 104,317,609 = 15.22% native

# define comparison
comparison_aci_m5 <- ((aci_reclass_mask * 10) + hard_class_m5)
names(comparison_aci_m5) <- "aci_vs_consensus_m5"

# 11 = LULC non-grass, hard_class non-grass
# 12 = LULC non-grass, hard_class tame
# 13 = LULC non-grass, hard_class native
# 
# 21 = LULC tame,      hard_class non-grass
# 22 = LULC tame,      hard_class tame
# 23 = LULC tame,      hard_class native
# 
# 31 = LULC native,    hard_class non-grass
# 32 = LULC native,    hard_class tame
# 33 = LULC native,    hard_class native

levels(comparison_aci_m5) <- data.frame(
  value = c(11,12,13,21,22,23,31,32,33),
  class = c(
    "both_non_grass",
    "ACI_non_CM_tame",
    "ACI_non_CM_native",
    "ACI_tame_CM_non",
    "both_tame",
    "ACI_tame_CM_native",
    "ACI_native_CM_non",
    "ACI_native_CM_tame",
    "both_native"
  )
)

# export
writeRaster(comparison_aci_m5, file.path(results_dir, "aci_vs_hardclass_m5.tif"),
            overwrite = TRUE,
            wopt = list(
              datatype = "INT2U",
              gdal = c("COMPRESS=ZSTD", "TILED=YES", "BIGTIFF=YES")))

comparison_aci_m5_export <- comparison_aci_m5

levels(comparison_aci_m5_export) <- NULL

writeRaster(
  comparison_aci_m5_export,
  file.path(results_dir, "aci_vs_hardclassm5_numeric.tif"),
  overwrite = TRUE,
  wopt = list(
    datatype = "INT2U",
    gdal = c(
      "COMPRESS=ZSTD",
      "TILED=YES",
      "BIGTIFF=YES"
    )
  )
)

# examine
plot(comparison_aci_m5)
freq(comparison_aci_m5)


#------- compare LUTS and consensus map mod 5 -------#

# BUT can I do it? LUTS classes are wrong

# reclassify LUTS to match
luts_reclass <- classify(luts, rcl = matrix(c(
  0, 1, # Non-grass
  1, 2, # Tame grassland
  2, 3), # Native grass
  ncol = 2, byrow = TRUE))

# crop LUTS to extent of consensus map
luts_reclass_crop <- crop(luts_reclass, hard_class)
luts_reclass_mask <- mask(luts_reclass_crop, hard_class)

# verify alignment
compareGeom(luts_reclass_mask, hard_class, stopOnError = TRUE)

plot(luts_reclass_mask)
freq(luts_reclass_mask)

# define comparison
comparison_luts_m5 <- ((luts_reclass_mask * 10) + hard_class_m5)
names(comparison_luts_m5) <- "luts_vs_consensus_m5"

# 11 = LULC non-grass, hard_class non-grass
# 12 = LULC non-grass, hard_class tame
# 13 = LULC non-grass, hard_class native
# 
# 21 = LULC tame,      hard_class non-grass
# 22 = LULC tame,      hard_class tame
# 23 = LULC tame,      hard_class native
# 
# 31 = LULC native,    hard_class non-grass
# 32 = LULC native,    hard_class tame
# 33 = LULC native,    hard_class native

levels(comparison_luts_m5) <- data.frame(
  value = c(11,12,13,21,22,23,31,32,33),
  class = c(
    "both_non_grass",
    "luts_non_CM_tame",
    "luts_non_CM_native",
    "luts_tame_CM_non",
    "both_tame",
    "luts_tame_CM_native",
    "luts_native_CM_non",
    "luts_native_CM_tame",
    "both_native"
  )
)

# export
writeRaster(comparison_luts_m5, file.path(results_dir, "luts_vs_hardclass_m5.tif"),
            overwrite = TRUE,
            wopt = list(
              datatype = "INT2U",
              gdal = c("COMPRESS=ZSTD", "TILED=YES", "BIGTIFF=YES")))

comparison_luts_m5_export <- comparison_luts_m5

levels(comparison_luts_m5_export) <- NULL

writeRaster(
  comparison_luts_m5_export,
  file.path(results_dir, "luts_vs_hardclass_m5_numeric.tif"),
  overwrite = TRUE,
  wopt = list(
    datatype = "INT2U",
    gdal = c(
      "COMPRESS=ZSTD",
      "TILED=YES",
      "BIGTIFF=YES"
    )
  )
)

# examine
plot(comparison_luts_m5)
freq(comparison_luts_m5)

#------- compare PGLC and consensus map -------#

pglc

# reclassify PGLC to match
pglc_reclass <- classify(pglc, rcl = matrix(c(
  0, 1, # Non-grass
  1, 2, # Tame grassland
  2, 3), # Native grass
  ncol = 2, byrow = TRUE))

# crop PGLC to extent of consensus map
pglc_reclass_crop <- crop(pglc_reclass, hard_class)
pglc_reclass_mask <- mask(pglc_reclass_crop, hard_class)

# verify alignment
compareGeom(pglc_reclass_mask, hard_class, stopOnError = TRUE)

plot(pglc_reclass_mask)
freq(pglc_reclass_mask)
#   layer value     count
# 1     1     1 179,307,864 = 63.71%
# 2     1     2  50,508,298 = 17.95%
# 3     1     3  51,634,502 = 18.35%

# define comparison
comparison_pglc_m5 <- ((pglc_reclass_mask * 10) + hard_class_m5)
names(comparison_pglc_m5) <- "pglc_vs_consensus_m5"

# 11 = LULC non-grass, hard_class non-grass
# 12 = LULC non-grass, hard_class tame
# 13 = LULC non-grass, hard_class native
# 
# 21 = LULC tame,      hard_class non-grass
# 22 = LULC tame,      hard_class tame
# 23 = LULC tame,      hard_class native
# 
# 31 = LULC native,    hard_class non-grass
# 32 = LULC native,    hard_class tame
# 33 = LULC native,    hard_class native

levels(comparison_pglc_m5) <- data.frame(
  value = c(11,12,13,21,22,23,31,32,33),
  class = c(
    "both_non_grass",
    "pglc_non_CM_tame",
    "pglc_non_CM_native",
    "pglc_tame_CM_non",
    "both_tame",
    "pglc_tame_CM_native",
    "pglc_native_CM_non",
    "pglc_native_CM_tame",
    "both_native"
  )
)

# export
writeRaster(comparison_pglc_m5, file.path(results_dir, "pglc_vs_hardclass_m5.tif"),
            overwrite = TRUE,
            wopt = list(
              datatype = "INT2U",
              gdal = c("COMPRESS=ZSTD", "TILED=YES", "BIGTIFF=YES")))

comparison_pglc_m5_export <- comparison_pglc_m5

levels(comparison_pglc_m5_export) <- NULL

writeRaster(
  comparison_pglc_m5_export,
  file.path(results_dir, "pglc_vs_hardclass_m5_numeric.tif"),
  overwrite = TRUE,
  wopt = list(
    datatype = "INT2U",
    gdal = c(
      "COMPRESS=ZSTD",
      "TILED=YES",
      "BIGTIFF=YES")))

# examine
plot(comparison_pglc_m5)
freq(comparison_pglc_m5)

#------- compare NGI and consensus map mod 5 -------#

# reclassify NGI to match
ngi_reclass <- classify(ngi, rcl = matrix(c(
  0, 1, # Non-grass
  1, 2, # Tame grassland
  2, 3), # Native grass
  ncol = 2, byrow = TRUE))

# crop NGI to extent of consensus map
ngi_reclass_crop <- crop(ngi_reclass, hard_class)
ngi_reclass_mask <- mask(ngi_reclass_crop, hard_class)

# verify alignment
compareGeom(ngi_reclass_mask, hard_class, stopOnError = TRUE)

plot(ngi_reclass_mask)
freq(ngi_reclass_mask)
#   layer value       count
# 1     1     1 310,074,752 = 58.34% non-grass
# 2     1     2 124,425,184 = 23.41% tame grass
# 3     1     3  97,023,669 = 18.25% native grass

# define comparison
comparison_ngi_m5 <- ((ngi_reclass_mask * 10) + hard_class_m5)
names(comparison_ngi_m5) <- "ngi_vs_consensus_m5"

# 11 = LULC non-grass, hard_class non-grass
# 12 = LULC non-grass, hard_class tame
# 13 = LULC non-grass, hard_class native
# 
# 21 = LULC tame,      hard_class non-grass
# 22 = LULC tame,      hard_class tame
# 23 = LULC tame,      hard_class native
# 
# 31 = LULC native,    hard_class non-grass
# 32 = LULC native,    hard_class tame
# 33 = LULC native,    hard_class native

levels(comparison_ngi_m5) <- data.frame(
  value = c(11,12,13,21,22,23,31,32,33),
  class = c(
    "both_non_grass",
    "ngi_non_CM_tame",
    "ngi_non_CM_native",
    "ngi_tame_CM_non",
    "both_tame",
    "ngi_tame_CM_native",
    "ngi_native_CM_non",
    "ngi_native_CM_tame",
    "both_native"
  )
)

# export
writeRaster(comparison_ngi_m5, file.path(results_dir, "ngi_vs_hardclass_m5.tif"),
            overwrite = TRUE,
            wopt = list(
              datatype = "INT2U",
              gdal = c("COMPRESS=ZSTD", "TILED=YES", "BIGTIFF=YES")))

comparison_ngi_m5_export <- comparison_ngi_m5

levels(comparison_ngi_m5_export) <- NULL

writeRaster(
  comparison_ngi_m5_export,
  file.path(results_dir, "ngi_vs_hardclass_m5_numeric.tif"),
  overwrite = TRUE,
  wopt = list(
    datatype = "INT2U",
    gdal = c("COMPRESS=ZSTD", "TILED=YES", "BIGTIFF=YES")))

# examine
plot(comparison_ngi_m5)
freq(comparison_ngi_m5)

#------- compare GPW and consensus map mod 5 -------#

# reclassify GPW to match
gpw_reclass <- classify(gpw, rcl = matrix(c(
  0, 1, # Non-grass
  1, 2, # Tame grassland
  2, 3), # Native grass
  ncol = 2, byrow = TRUE))

# crop GPW to extent of consensus map
gpw_reclass_crop <- crop(gpw_reclass, hard_class)
gpw_reclass_mask <- mask(gpw_reclass_crop, hard_class)

# verify alignment
compareGeom(gpw_reclass_mask, hard_class, stopOnError = TRUE)

plot(gpw_reclass_mask)
freq(gpw_reclass_mask)
#   layer value     count
# 1     1     1 401,695,078 = 64.49%
# 2     1     2  79,645,851 = 12.79%
# 3     1     3 141,494,500 = 22.72%

# define comparison
comparison_gpw_m5 <- ((gpw_reclass_mask * 10) + hard_class_m5)
names(comparison_gpw_m5) <- "gpw_vs_consensus_m5"

# 11 = LULC non-grass, hard_class non-grass
# 12 = LULC non-grass, hard_class tame
# 13 = LULC non-grass, hard_class native
# 
# 21 = LULC tame,      hard_class non-grass
# 22 = LULC tame,      hard_class tame
# 23 = LULC tame,      hard_class native
# 
# 31 = LULC native,    hard_class non-grass
# 32 = LULC native,    hard_class tame
# 33 = LULC native,    hard_class native

levels(comparison_gpw_m5) <- data.frame(
  value = c(11,12,13,21,22,23,31,32,33),
  class = c(
    "both_non_grass",
    "gpw_non_CM_tame",
    "gpw_non_CM_native",
    "gpw_tame_CM_non",
    "both_tame",
    "gpw_tame_CM_native",
    "gpw_native_CM_non",
    "gpw_native_CM_tame",
    "both_native"
  )
)

# export
writeRaster(comparison_gpw_m5, file.path(results_dir, "gpw_vs_hardclass_m5.tif"),
            overwrite = TRUE,
            wopt = list(
              datatype = "INT2U",
              gdal = c("COMPRESS=ZSTD", "TILED=YES", "BIGTIFF=YES")))

comparison_gpw_m5_export <- comparison_gpw_m5

levels(comparison_gpw_m5_export) <- NULL

writeRaster(
  comparison_gpw_m5_export,
  file.path(results_dir, "gpw_vs_hardclass_m5_numeric.tif"),
  overwrite = TRUE,
  wopt = list(
    datatype = "INT2U",
    gdal = c("COMPRESS=ZSTD", "TILED=YES", "BIGTIFF=YES")))

# examine
plot(comparison_gpw_m5)
freq(comparison_gpw_m5)

#------- compare ECCC and consensus map -------#

# reclassify ECCC to match
eccc_reclass <- classify(eccc, rcl = matrix(c(
  0, 1, # Non-grass
  1, 2, # Tame grassland
  2, 3), # Native grass
  ncol = 2, byrow = TRUE))

# crop ECCC to extent of consensus map
eccc_reclass_crop <- crop(eccc_reclass, hard_class)
eccc_reclass_mask <- mask(eccc_reclass_crop, hard_class)

# verify alignment
compareGeom(eccc_reclass_mask, hard_class, stopOnError = TRUE)

plot(eccc_reclass_mask)
freq(eccc_reclass_mask)
# layer value     count
# 1     1     1 512,288,174 = 75.17% non
# 2     1     2  97,398,954 = 14.29% tame
# 3     1     3  71,853,680 = 10.54% native

# define comparison
comparison_eccc_m5 <- ((eccc_reclass_mask * 10) + hard_class_m5)
names(comparison_eccc_m5) <- "eccc_vs_consensus_m5"

# 11 = LULC non-grass, hard_class non-grass
# 12 = LULC non-grass, hard_class tame
# 13 = LULC non-grass, hard_class native
# 
# 21 = LULC tame,      hard_class non-grass
# 22 = LULC tame,      hard_class tame
# 23 = LULC tame,      hard_class native
# 
# 31 = LULC native,    hard_class non-grass
# 32 = LULC native,    hard_class tame
# 33 = LULC native,    hard_class native

levels(comparison_eccc_m5) <- data.frame(
  value = c(11,12,13,21,22,23,31,32,33),
  class = c(
    "both_non_grass",
    "eccc_non_CM_tame",
    "eccc_non_CM_native",
    "eccc_tame_CM_non",
    "both_tame",
    "eccc_tame_CM_native",
    "eccc_native_CM_non",
    "eccc_native_CM_tame",
    "both_native"
  )
)

# export
writeRaster(comparison_eccc_m5, file.path(results_dir, "eccc_vs_hardclass_m5.tif"),
            overwrite = TRUE,
            wopt = list(
              datatype = "INT2U",
              gdal = c("COMPRESS=ZSTD", "TILED=YES", "BIGTIFF=YES")))

comparison_eccc_m5_export <- comparison_eccc_m5

levels(comparison_eccc_m5_export) <- NULL

writeRaster(
  comparison_eccc_m5_export,
  file.path(results_dir, "eccc_vs_hardclass_m5_numeric.tif"),
  overwrite = TRUE,
  wopt = list(
    datatype = "INT2U",
    gdal = c("COMPRESS=ZSTD", "TILED=YES", "BIGTIFF=YES")))

# examine
plot(comparison_eccc_m5)
freq(comparison_eccc_m5)

#------ grassland retention in consensus map ------#

# what percentage of each grass dataset is included in the consensus layer? 

# focused on grassland inclusion rather than exact class agreement:
# grass_retention_pct: percentage of all grassland pixels (tame + native) in the source 
#                      dataset that are still grassland in the consensus.
# tame_retention_pct: percentage of source tame pixels that remain grassland (either 
#                     tame or native) in the consensus.
# native_retention_pct: percentage of source native pixels that remain grassland (either 
#                      tame or native) in the consensus.

#define function
calc_retention <- function(comp_raster, dataset_name){
  
  f <- freq(comp_raster)
  
  counts <- setNames(
    f$count,
    as.character(f$value)
  )
  
  # helper
  get_count <- function(x){
    ifelse(x %in% names(counts), counts[x], 0)
  }
  
  tame_non    <- get_count(paste0(dataset_name, "_tame_CM_non"))
  tame_native <- get_count(paste0(dataset_name, "_tame_CM_native"))
  
  native_non  <- get_count(paste0(dataset_name, "_native_CM_non"))
  native_tame <- get_count(paste0(dataset_name, "_native_CM_tame"))
  
  both_tame   <- get_count("both_tame")
  both_native <- get_count("both_native")
  
  # totals in source dataset
  tame_total <-
    tame_non +
    both_tame +
    tame_native
  
  native_total <-
    native_non +
    native_tame +
    both_native
  
  grass_total <-
    tame_total +
    native_total
  
  # retained as grass in consensus
  tame_retained <-
    both_tame +
    tame_native
  
  native_retained <-
    native_tame +
    both_native
  
  grass_retained <-
    tame_retained +
    native_retained
  
  data.frame(
    dataset = dataset_name,
    
    grass_total = grass_total,
    grass_retained = grass_retained,
    grass_retention_pct =
      100 * grass_retained / grass_total,
    
    tame_total = tame_total,
    tame_retained = tame_retained,
    tame_retention_pct =
      100 * tame_retained / tame_total,
    
    native_total = native_total,
    native_retained = native_retained,
    native_retention_pct =
      100 * native_retained / native_total
  )
}

# run
retention_summary_m5 <- rbind(
  calc_retention(comparison_aci_m5,  "ACI"),
  calc_retention(comparison_pglc_m5, "pglc"),
  calc_retention(comparison_ngi_m5,  "ngi"),
  calc_retention(comparison_gpw_m5,  "gpw"),
  calc_retention(comparison_eccc_m5, "eccc")
)

retention_summary_m5


# cleaner output
retention_summary_m5 <- retention_summary_m5 |>
  transform(
    grass_retention_pct  = round(grass_retention_pct, 1),
    tame_retention_pct   = round(tame_retention_pct, 1),
    native_retention_pct = round(native_retention_pct, 1)
  )

retention_summary_m5

write.csv(
  retention_summary_m5,
  file.path(results_dir, "grassland_retention_summary_m5.csv"),
  row.names = FALSE
)
