
###################################################################################
##                 SETUP FOR GRASSLAND CONSENSUS MAPPING ANALYSIS                ##
###################################################################################

##----------------------- Alison Long - June 2026 -------------------------##

# step 8b: create redundancy penalty raster to assist creation of consensus grassland map 

# penalty = low if highly redundant
# penalty = high if unique

#----- 0) Setup ------#

# Load packages
library(tidyverse)
library(terra)

# Directories and files
data_dir <- "R:/CanadianGrasslands/WG_Classification/results/accuracy_assessment/binary/ci_derived"
out_dir  <- "R:/CanadianGrasslands/WG_Classification/results/accuracy_assessment/binary/ci_derived"

pairwise_long <- file.path(data_dir, "pairwise_binary_spatial_agreement_01Jun2026.csv")

#----- 1) Symmetrize matrix -----#

# make symmetric table (instead of only triangle)
pairwise_full <- bind_rows(
  pairwise_long,
  pairwise_long %>%
    transmute(
      raster_i = raster_j,
      raster_j = raster_i,
      mcc = mcc
    )
)

#----- 2) Compute mean similarity per raster -----#
mean_similarity <- pairwise_full %>%
  group_by(raster_i) %>%
  summarise(
    mean_mcc = mean(mcc, na.rm = TRUE)
  ) %>%
  rename(raster = raster_i)

#----- 3) Convert similarity to penalty -----#
penalty_df <- mean_similarity %>%
  mutate(
    penalty = 1 / (1 + mean_mcc)
  )

#----- 4) Clamp values to prevent zero weights and overly extreme penalties -----#
penalty_df <- penalty_df %>%
  mutate(
    penalty = pmin(pmax(penalty, 0.2), 1)
  )

#----- 5) Convert to raster stack -----#
library(terra)

template <- rast("R:/CanadianGrasslands/WG_Classification/processed/aafc_aci/aligned/aafc_aci_merged_reproj_algn_binary_noshrub.tif")

penalty_rasters <- lapply(penalty_df$raster, function(r_name) {
  
  val <- penalty_df$penalty[penalty_df$raster == r_name]
  
  r <- template
  values(r) <- val
  
  names(r) <- r_name
  
  return(r)
})

penalty_bin_stack <- rast(penalty_rasters)

#----- 6) Save -----#
writeRaster(
  penalty_bin_stack,
  filename = file.path(out_dir, "penalty_bin.tif"),
  overwrite = TRUE,
  wopt = list(
    datatype="FLT4S",
    gdal=c("COMPRESS=ZSTD","TILED=YES","BIGTIFF=YES")
  )
)