###################################################################################
##                      ANALYSIS: CONSENSUS MAP REFINEMENT                       ##
###################################################################################

##----------------------- Alison Long - July 2026 -------------------------

# STEP 10: Use diagnostics results from consensus mapping analysis to refine
#          the final product. Create and test modifications to the map.

library(tidyverse)
library(terra)
library(sf)

results_dir <- "R:/CanadianGrasslands/WG_Classification/results/consensus_mapping"
frac_bin_file <- file.path(results_dir, "new_fusion_global_norm_modext_classdifficultynorm_31Aug2026/fraction_binary_support.tif")
frac_tri_file <- file.path(results_dir, "new_fusion_global_norm_modext_classdifficultynorm_31Aug2026/fraction_triclass_support.tif")

fraction_binary_support   <- rast(frac_bin_file)
fraction_triclass_support <- rast(frac_tri_file)

class(fraction_binary_support)
methods("global")

terra::global(fraction_binary_support, "mean", na.rm = TRUE)
terra::global(fraction_triclass_support, "mean", na.rm = TRUE)

global(fraction_binary_support, median, na.rm = TRUE)
global(fraction_triclass_support, median, na.rm = TRUE)

nt_stack_file <- file.path(results_dir, "new_fusion_global_norm_modext_26Aug2026/nt_stack.tif")
nt_stack <- rast(nt_stack_file)

# actual class prevalence across triclass stack
global_not_freq <- global(
  app(nt_stack == 0, mean, na.rm = TRUE),
  "mean",
  na.rm = TRUE
)[1,1]

global_tame_freq <- global(
  app(nt_stack == 1, mean, na.rm = TRUE),
  "mean",
  na.rm = TRUE
)[1,1]

global_native_freq <- global(
  app(nt_stack == 2, mean, na.rm = TRUE),
  "mean",
  na.rm = TRUE
)[1,1]

c(
  not = global_not_freq,
  tame = global_tame_freq,
  native = global_native_freq
)

# not       tame      native 
# 0.7425205 0.1073789 0.1501006

# prevalence before weighting:
mean(
  values(app(nt_stack == 0, mean, na.rm = TRUE)),
  na.rm = TRUE
)

mean(
  values(app(nt_stack == 1, mean, na.rm = TRUE)),
  na.rm = TRUE
)

mean(
  values(app(nt_stack == 2, mean, na.rm = TRUE)),
  na.rm = TRUE
)

#-------------------------------------------------
# EXPLORING 

results_dir <- "R:/CanadianGrasslands/WG_Classification/results/consensus_mapping"

hard_class_file     <- file.path(results_dir, "new_fusion_global_norm_modext_classdifficultynorm_31Aug2026/hard_class.tif")
missed_tame_file    <- file.path(results_dir, "new_fusion_global_norm_modext_classdifficultynorm_31Aug2026/missed_tame.tif")
tame_possible1_file <- file.path(results_dir, "new_fusion_global_norm_modext_classdifficultynorm_31Aug2026/tame_possible1.tif")
tame_possible2_file <- file.path(results_dir, "new_fusion_global_norm_modext_classdifficultynorm_31Aug2026/tame_possible2.tif")
entropy_norm_file   <- file.path(results_dir, "new_fusion_global_norm_modext_classdifficultynorm_31Aug2026/entropy_norm.tif")

p_not_file     <- file.path(results_dir, "new_fusion_global_norm_modext_classdifficultynorm_31Aug2026/p_not.tif")
p_tame_file    <- file.path(results_dir, "new_fusion_global_norm_modext_classdifficultynorm_31Aug2026/p_tame.tif")
p_native_file  <- file.path(results_dir, "new_fusion_global_norm_modext_classdifficultynorm_31Aug2026/p_native.tif")

hard_class     <- rast(hard_class_file)
missed_tame    <- rast(missed_tame_file)
tame_possible1 <- rast(tame_possible1_file)
tame_possible2 <- rast(tame_possible1_file)
entropy_norm   <- rast(entropy_norm_file)

p_not    <- rast(p_not_file)
p_tame   <- rast(p_tame_file)
p_native <- rast(p_native_file)

#-------- number of pixels ---------#
global(p_not, "sum", na.rm = TRUE)    # = 509,484,529 pixels
global(p_tame, "sum", na.rm = TRUE)   # =  73,220,128 pixels
global(p_native, "sum", na.rm = TRUE) # = 103,143,181 pixels

#-------- borderline grass ---------#
borderline_grass <- p_not > 0.4 & p_not < 0.6
plot(borderline_grass, main = "Borderline grass")

# number of pixels = 42,423,661
global(borderline_grass, "sum", na.rm = TRUE)

# fraction of all non-grass pixels = 0.08084522
not_wins <- hard_class == 1
fraction_borderline <- global(borderline_grass, "sum", na.rm = TRUE)[1,1] /
  global(not_wins,"sum", na.rm = TRUE)[1,1]
fraction_borderline

writeRaster(borderline_grass, file.path(results_dir, "borderline_grass.tif"),
            overwrite = TRUE, wopt = list(datatype = "INT1U",
                                          gdal = c("COMPRESS=ZSTD", "TILED=YES", "BIGTIFF=YES")))

#----- p_tame close second to p_not -----#
not_vs_tame_close <- (
  (p_not >= p_tame) &
    (p_not >= p_native) &
    (p_tame >= p_native) &
    ((p_not - p_tame) <= 0.1))

names(not_vs_tame_close) <- "not_vs_tame_close"

plot(not_vs_tame_close, main = "p_tame close second to p_not")

writeRaster(not_vs_tame_close, file.path(results_dir, "tame_close_to_not.tif"),
  overwrite = TRUE,wopt = list(datatype = "INT1U",
    gdal = c("COMPRESS=ZSTD", "TILED=YES", "BIGTIFF=YES")))

# number of pixels = 7,731,370
global(not_vs_tame_close, "sum", na.rm = TRUE)

# fraction of all non-grass pixels = 0.01473339
fraction_tame_close <- global(not_vs_tame_close, "sum", na.rm = TRUE)[1,1] /
  global(not_wins,"sum", na.rm = TRUE)[1,1]
fraction_tame_close

#----- p_tame close-ish second to p_not -----#
not_vs_tame_closeish <- (
  (p_not >= p_tame) &
    (p_not >= p_native) &
    (p_tame >= p_native) &
    ((p_not - p_tame) <= 0.2))

names(not_vs_tame_closeish) <- "not_vs_tame_closeish"

plot(not_vs_tame_closeish, main = "p_tame close-ish second to p_not")

writeRaster(not_vs_tame_closeish, file.path(results_dir, "tame_closeish_to_not.tif"),
            overwrite = TRUE,wopt = list(datatype = "INT1U",
                                         gdal = c("COMPRESS=ZSTD", "TILED=YES", "BIGTIFF=YES")))

# number of pixels = 14,938,657
global(not_vs_tame_closeish, "sum", na.rm = TRUE)

# fraction of all non-grass pixels = 0.02846805
fraction_tame_closeish <- global(not_vs_tame_closeish, "sum", na.rm = TRUE)[1,1] /
  global(not_wins,"sum", na.rm = TRUE)[1,1]
fraction_tame_closeish

#----- p_tame close second to p_native -----#
native_vs_tame_close <- (
  (p_native >= p_tame) &
    (p_native >= p_not) &
    ((p_native - p_tame) <= 0.1)
)

plot(native_vs_tame_close, main = "p_tame close second to p_native")

writeRaster(native_vs_tame_close, file.path(results_dir, "tame_close_to_native.tif"),
            overwrite = TRUE,wopt = list(datatype = "INT1U",
                                         gdal = c("COMPRESS=ZSTD", "TILED=YES", "BIGTIFF=YES")))

global(native_vs_tame_close, "sum", na.rm = TRUE)

#----- p_native close second to p_tame -----#
tame_vs_native_close <- (
  (p_tame >= p_native) &
    (p_tame >= p_not) &
    ((p_tame - p_native) <= 0.1)
)

plot(tame_vs_native_close, main = "p_tame close second to p_native")

writeRaster(tame_vs_native_close, file.path(results_dir, "native_close_to_tame.tif"),
            overwrite = TRUE,wopt = list(datatype = "INT1U",
                                         gdal = c("COMPRESS=ZSTD", "TILED=YES", "BIGTIFF=YES")))

global(tame_vs_native_close, "sum", na.rm = TRUE)

#----- non-grass pixels within 0.1 of becoming grass -----#
grass_runner_up <- app(c(p_tame, p_native), max, na.rm = TRUE)

names(grass_runner_up) <- "grass_runner_up"

not_vs_grass_close <- (
  (p_not >= grass_runner_up) &
    ((p_not - grass_runner_up) <= 0.1)
)

plot(not_vs_grass_close, main = "Non-grass pixels within 0.1 of becoming grass")

writeRaster(not_vs_grass_close, file.path(results_dir, "grass_close_to_not.tif"),
            overwrite = TRUE, wopt = list(datatype = "INT1U",
                                         gdal = c("COMPRESS=ZSTD", "TILED=YES", "BIGTIFF=YES")))

# number of pixels = 12,460,415
global(not_vs_grass_close, "sum", na.rm = TRUE)

# fraction of all non-grass pixels = 0.02374536
fraction_not_grass_close <-
  global(not_vs_grass_close, "sum", na.rm = TRUE)[1,1] /
  global(not_wins, "sum", na.rm = TRUE)[1,1]
fraction_not_grass_close

#------ missed tame -------#

# number of pixels = 37,170,678
global(missed_tame, "sum", na.rm = TRUE)

#------ tame possible 1 -------#

# number of pixels = 100,350,280
global(tame_possible1, "sum", na.rm = TRUE)

#------ tame possible 2 -------#

tame_possible2 <- p_tame >= 0.25 & entropy_norm >= 0.5

writeRaster(tame_possible2, file.path(results_dir, "tame_possible2.tif"),
            overwrite = TRUE, wopt = list(datatype = "INT1U",
                                          gdal = c("COMPRESS=ZSTD", "TILED=YES", "BIGTIFF=YES")))

# number of pixels = 106,984,471
global(tame_possible2, "sum", na.rm = TRUE)

##########################################################
##          HARD CLASSIFICATION MODIFICATIONS           ##
##########################################################

#------- Modification 1: mild non-grass penalty -------#

p_not_m1 <- p_not * 0.95

prob_sum_m1 <- p_not_m1 + p_tame + p_native

p_not_m1    <- p_not_m1 / prob_sum_m1
p_tame_m1   <- p_tame   / prob_sum_m1
p_native_m1 <- p_native / prob_sum_m1

hard_class_m1 <- which.max(c(p_not_m1, p_tame_m1, p_native_m1))

writeRaster(hard_class_m1, file.path(results_dir, "hard_class_m1_mild_nongrass_penalty.tif"),
            overwrite = TRUE, wopt = list(datatype = "INT1U",
                                          gdal = c("COMPRESS=ZSTD", "TILED=YES", "BIGTIFF=YES")))

# make sure probabilities still sum to 1
chk_m1 <- global(p_not_m1 + p_tame_m1 + p_native_m1, c("min","max","mean"), na.rm = TRUE)
print(chk_m1)

# measure how many pixels changed = 2,833,935
changed_m1 <- hard_class_m1 != hard_class
global(changed_m1, "sum", na.rm = TRUE)

# measure proportion of pixels that changed = 0.004132017
prop_changed_m1 <-
  global(changed_m1, "sum", na.rm = TRUE)[1,1] /
  global(!is.na(hard_class), "sum", na.rm = TRUE)[1,1]
prop_changed_m1

# build transition matrix (1 = non-grass, 2 = tame, 3 = native)
terra::crosstab(c(hard_class, hard_class_m1))
#                   m1 result
# original          1         2         3
# 1         521917702   1750501   1083434
# 2                 0  63179602         0
# 3                 0         0  97916599

# 1,750,501 non-grass became tame grass
# 1,083,434 non-grass became native grass

#------- Modification 2: moderate non-grass penalty -------#

p_not_m2 <- p_not * 0.90

prob_sum_m2 <- p_not_m2 + p_tame + p_native

p_not_m2    <- p_not_m2 / prob_sum_m2
p_tame_m2   <- p_tame   / prob_sum_m2
p_native_m2 <- p_native / prob_sum_m2

hard_class_m2 <- which.max(c(p_not_m2, p_tame_m2, p_native_m2))

writeRaster(hard_class_m2, file.path(results_dir, "hard_class_m2_moderate_nongrass_penalty.tif"),
            overwrite = TRUE, wopt = list(datatype = "INT1U",
                                          gdal = c("COMPRESS=ZSTD", "TILED=YES", "BIGTIFF=YES")))

# make sure probabilities still sum to 1
chk_m2 <- global(p_not_m2 + p_tame_m2 + p_native_m2, c("min","max","mean"), na.rm = TRUE)
print(chk_m2)

# measure how many pixels changed = 5,851,101
changed_m2 <- hard_class_m2 != hard_class
global(changed_m2, "sum", na.rm = TRUE)

# measure proportion of pixels that changed = 0.008531194
prop_changed_m2 <-
  global(changed_m2, "sum", na.rm = TRUE)[1,1] /
  global(!is.na(hard_class), "sum", na.rm = TRUE)[1,1]
prop_changed_m2

# build transition matrix
terra::crosstab(c(hard_class, hard_class_m2))
#                  m2 result
# original         1         2         3
#        1 518900536   3641126   2209975
#        2         0  63179602         0
#        3         0         0  97916599

# 3,641,126 non-grass became tame grass
# 2,209,975 non-grass became native grass

#------- Modification 3: strong non-grass penalty -------#

p_not_m3 <- p_not * 0.80

prob_sum_m3 <- p_not_m3 + p_tame + p_native

p_not_m3    <- p_not_m3 / prob_sum_m3
p_tame_m3   <- p_tame   / prob_sum_m3
p_native_m3 <- p_native / prob_sum_m3

hard_class_m3 <- which.max(c(p_not_m3, p_tame_m3, p_native_m3))

writeRaster(hard_class_m3, file.path(results_dir, "hard_class_m3_strong_nongrass_penalty.tif"),
            overwrite = TRUE, wopt = list(datatype = "INT1U",
                                          gdal = c("COMPRESS=ZSTD", "TILED=YES", "BIGTIFF=YES")))

# make sure probabilities still sum to 1
chk_m3 <- global(p_not_m3 + p_tame_m3 + p_native_m3, c("min","max","mean"), na.rm = TRUE)
print(chk_m3)

# measure how many pixels changed = 12,356,266
changed_m3 <- hard_class_m3 != hard_class
global(changed_m3, "sum", na.rm = TRUE)

# measure proportion of pixels that changed = 0.01801605
prop_changed_m3 <-
  global(changed_m3, "sum", na.rm = TRUE)[1,1] /
  global(!is.na(hard_class), "sum", na.rm = TRUE)[1,1]
prop_changed_m3

# build transition matrix
terra::crosstab(c(hard_class, hard_class_m3))
#                  m3 result
# original         1         2         3
# 1        512395371   7803428   4552838
# 2                0  63179602         0
# 3                0         0  97916599

# 7,803,428 non-grass became tame grass
# 4,552,838 non-grass became native grass

#------- Modification 4: boost all grass -------#

grass_boost <- 1.10

p_tame_m4 <- p_tame * grass_boost
p_native_m4 <- p_native * grass_boost

prob_sum_m4 <- p_not + p_tame_m4 + p_native_m4

p_not_m4    <- p_not / prob_sum_m4
p_tame_m4   <- p_tame_m4 / prob_sum_m4
p_native_m4 <- p_native_m4 / prob_sum_m4

hard_class_m4 <- which.max(c(p_not_m4, p_tame_m4, p_native_m4))

writeRaster(hard_class_m4, file.path(results_dir, "hard_class_m4_boost_all_grass.tif"),
            overwrite = TRUE, wopt = list(datatype = "INT1U",
                                          gdal = c("COMPRESS=ZSTD", "TILED=YES", "BIGTIFF=YES")))

# make sure probabilities still sum to 1
chk_m4 <- global(p_not_m4 + p_tame_m4 + p_native_m4, c("min","max","mean"), na.rm = TRUE)
print(chk_m4)

# measure how many pixels changed = 5,226,733
changed_m4 <- hard_class_m4 != hard_class
global(changed_m4, "sum", na.rm = TRUE)

# measure proportion of pixels that changed = 0.0076208347
prop_changed_m4 <-
  global(changed_m4, "sum", na.rm = TRUE)[1,1] /
  global(!is.na(hard_class), "sum", na.rm = TRUE)[1,1]
prop_changed_m4

# build transition matrix
terra::crosstab(c(hard_class, hard_class_m4))
#                  m4 result
# original         1         2         3
# 1        519524904   3276142   1950591
# 2                0  63179602         0
# 3                0         0  97916599

# 3,276,142 non-grass became tame grass
# 1,950,591 non-grass became native grass

#------- Modification 5: boost tame grass -------#

p_tame_m5 <- p_tame * 1.25

prob_sum_m5 <- p_not + p_tame_m5 + p_native

p_not_m5    <- p_not / prob_sum_m5
p_tame_m5   <- p_tame_m5 / prob_sum_m5
p_native_m5 <- p_native / prob_sum_m5

hard_class_m5 <- which.max(c(p_not_m5, p_tame_m5, p_native_m5))

writeRaster(hard_class_m5, file.path(results_dir, "hard_class_m5_boost_tame_grass.tif"),
            overwrite = TRUE, wopt = list(datatype = "INT1U",
                                          gdal = c("COMPRESS=ZSTD", "TILED=YES", "BIGTIFF=YES")))

# make sure probabilities still sum to 1
chk_m5 <- global(p_not_m5 + p_tame_m5 + p_native_m5, c("min","max","mean"), na.rm = TRUE)
print(chk_m5)

# measure how many pixels changed = 12,580,916
changed_m5 <- hard_class_m5 != hard_class
global(changed_m5, "sum", na.rm = TRUE)

# measure proportion of pixels that changed = 0.0183435965
prop_changed_m5 <-
  global(changed_m5, "sum", na.rm = TRUE)[1,1] /
  global(!is.na(hard_class), "sum", na.rm = TRUE)[1,1]
prop_changed_m5

# build transition matrix
terra::crosstab(c(hard_class, hard_class_m5))
#                  m5 result
# original        1         2         3
# 1       516599762   8151875         0
# 2               0  63179602         0
# 3               0   4429041  93487558

# 8,151,875 non-grass became tame grass
# 4,429,041 native became tame grass

# original probabilites in pixels that changed
global(mask(p_not, changed_m5), "mean", na.rm=TRUE)
global(mask(p_tame, changed_m5), "mean", na.rm=TRUE)
global(mask(p_native, changed_m5), "mean", na.rm=TRUE)

#------- Modification 6: ambiguous pixels to tame -------#
alt_tame <- (
  hard_class == 1 &
    p_tame >= 0.30 &
    entropy_norm >= 0.50
)

hard_class_alt <- classify(...)
