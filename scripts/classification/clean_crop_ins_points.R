
# v2 of crop insurance-derived groundtruth points built with bad ecodistricts
# doing my best to clean them 

# ----------------------------
# 1) Setup
# ----------------------------

# load packages
library(readxl)
library(readr)
library(tidyverse)

# ----------------------------
# 2) File paths 
# ----------------------------

points_excel_path <- "R:/CanadianGrasslands/WG_Classification/processed/groundtruthed_data_points/ci_derived/point_raster_extract_v2_grass_v_nongrass.xlsx"
nt_points_excel_path <- "R:/CanadianGrasslands/WG_Classification/processed/groundtruthed_data_points/ci_derived/point_raster_extract_v2_native_v_tame.xlsx"
metadata_csv_path <- "R:/CanadianGrasslands/WG_Classification/github/tomorrows-prairies/processeddatasets/Grass/metadata_full_extent_binary_pli_sep_noplow.csv"

out_dir <- "R:/CanadianGrasslands/WG_Classification/processed/groundtruthed_data_points/ci_derived"

# ----------------------------
# 3) Configuration
# ----------------------------

# ----------------------------
# 4) Load points and normalize basic fields
# ----------------------------
points <- readxl::read_excel(
  points_excel_path,
  col_types = c("numeric", "text", "numeric", rep("guess", 4), rep("numeric", 26))
)

names(points)
str(points)
unique(points$LC_Class_N)
unique(points$Source_ID)

# Rename columns
points2 <- points %>%
  rename(
    Pt_ID_O    = Point_ID,
    Sorc_ID  = Source_ID,
    Yer_Cll  = Year_Colle,
    LC_Cl_N  = LC_Class_N,
    Ecds_ID = ecodist,
    Provinc  = PROV,
    Ecrgn_n  = ecoreg,
    aci2024  = aafc_aci_merged_reproj_algn_binary_noshrub,
    #aci2025  = aafc_aci_2025_merged_reproj_algn_binary_noshrub,
    luts10   = luts10_reproj_algn_binary_masked,
    luts30   = luts30_reproj_algn_binary_masked,
    cec      = nalcms_reproj_algn_binary_noshrub_masked,
    ngi      = cfga_ngi_merged_reproj_algn_binary,
    cgr      = cgr_reproj_algn_binary_masked, 
    eccc     = eccc_reclass_reproj_algn_binary_masked,
    esa      = esa_reproj_algn_binary_noshrub_masked,
    glc      = glc_reproj_algn_binary_noshrub_masked,
    cnvc     = cnvc_algn_binary_masked,
    lcc      = lcc_reproj_algn_binary_noshrub_masked,
    pudl     = pudl_reclass_reproj_algn_binary_noshrub, 
    pli_ap   = sk_pli_ap_reproj_algn_binary_noshrub,
    pli_cu   = sk_pli_cu_reproj_algn_binary_noshrub,
    pli_mg   = sk_pli_mg_reproj_algn_binary_noshrub,
    pli_mm   = sk_pli_mm_reproj_algn_binary_noshrub,
    lcr      = lcr_reproj_algn_binary_masked,
    gpw      = gpw_reproj_algn_binary_noshrub_masked,
    plvi     = plvi_reproj_algn_binary_gte2015,
    hfi      = hfi_reproj_algn_binary_gte2015,
    gcc      = gcc_reproj_algn_binary,
    dw       = dw_reproj_algn_masked_binary,
    esri     = esri_reproj_algn_masked_binary,
    mgi      = mb_mgi_reproj_algn_binary_noshrub,
    pglc     = pglc_reproj_algn_binary,
    plow     = plow_reproj_algn_intact_masked)

str(points2)

# ----------------------------
# 5) Clean up ecodistricts/ecoregions
# ----------------------------
# Get mapping (one row per ecodistrict)
ecodist_to_regions <- points2 %>%
  distinct(Ecds_ID, Ecrgn_n)

# Check whether an ecodistrict belongs to multiple regions/provinces
ecodist_to_regions %>%
  count(Ecds_ID) %>%
  filter(n > 1)

#  Ecds_ID     n
#1     618     2
#2     815     2
#3     828     2
#4     833     2
#5     838     2

points2 %>%
  distinct(Ecds_ID, Ecrgn_n) %>%
  filter(Ecds_ID %in%
           (points2 %>%
              distinct(Ecds_ID, Ecrgn_n) %>%
              count(Ecds_ID) %>%
              filter(n > 1) %>%
              pull(Ecds_ID)))

#      Ecds_ID Ecrgn_n               
#1     618 Western Alberta Upland
#2     833 Mixed Grassland       
#3     838 Cypress Upland        
#4     815 Mixed Grassland       
#5     828 Mixed Grassland       
#6     833 Cypress Upland        
#7     815 Cypress Upland        
#8     828 Moist Mixed Grassland 
#9     618 Peace Lowland         
#10    838 Mixed Grassland  

# Correct ecoregions are:
# 618 = 145 = Western Alberta Upland
# 815 = 159 = Mixed Grassland 
# 828 = 159 = Mixed Grassland 
# 833 = 159 = Mixed Grassland 
# 838 = 160 = Cypress Upland

# Make corrections
points2 <- points2 %>%
  mutate(Ecrgn_n = case_when(
    Ecds_ID == "618" ~ "Western Alberta Upland",
    Ecds_ID == "815" ~ "Mixed Grassland",
    Ecds_ID == "828" ~ "Mixed Grassland",
    Ecds_ID == "833" ~ "Mixed Grassland",
    Ecds_ID == "838" ~ "Cypress Upland",
    TRUE ~ Ecrgn_n
  ))

# Get mapping (one row per ecodistrict)
ecodist_to_regions <- points2 %>%
  distinct(Ecds_ID, Ecrgn_n)

# Check whether an ecodistrict belongs to multiple regions/provinces
ecodist_to_regions %>%
  count(Ecds_ID) %>%
  filter(n > 1)

# export
writexl::write_xlsx(points2, file.path(out_dir, "point_raster_extract_v2_grass_v_nongrass_clean.xlsx"))

#######################################################################################
#######################################################################################

# For native vs tame

# ----------------------------
# 4) Load points and normalize basic fields
# ----------------------------
ntpoints <- readxl::read_excel(
  nt_points_excel_path,
  col_types = c("numeric", "text", "numeric", "text", "numeric", rep("guess", 2), rep("numeric", 15))
)

names(ntpoints)
str(ntpoints)
unique(ntpoints$LC_Class_N)
unique(ntpoints$Source_ID)

# Rename columns
points3 <- ntpoints %>%
  rename(
    Pt_ID_O    = Point_ID,
    Sorc_ID  = Source_ID,
    Yer_Cll  = Year_Colle,
    LC_Cl_N  = LC_Class_N,
    Ecds_ID = ecodist,
    Provinc  = PROV,
    Ecrgn_n  = ecoreg,
    aci      = aafc_aci_merged_reproj_algn_natvtam_noshrub,
    ngi      = cfga_ngi_merged_reproj_algn_natvtam_good,
    cgr      = cgr_reproj_algn_natvtam_masked, 
    eccc     = eccc_reclass_reproj_algn_natvtam_masked,
    cnvc     = cnvc_algn_natvtam_masked,
    pudl     = pudl_reclass_reproj_algn_natvtam_noshrub, 
    pli_ap   = sk_pli_ap_reproj_algn_natvtam_noshrub,
    pli_cu   = sk_pli_cu_reproj_algn_natvtam_noshrub,
    pli_mg   = sk_pli_mg_reproj_algn_natvtam_noshrub,
    pli_mm   = sk_pli_mm_reproj_algn_natvtam_noshrub,
    gpw      = gpw_reproj_algn_natvtam_noshrub_masked,
    plvi     = plvi_reproj_algn_natvtam_gte2015,
    hfi      = hfi_reproj_algn_binary_gte2015,
    mgi      = mb_mgi_reproj_algn_natvtam_noshrub,
    pglc     = pglc_reproj_algn_natvtam)

str(points3)

# ----------------------------
# 5) Clean up ecodistricts/ecoregions
# ----------------------------
# Get mapping (one row per ecodistrict)
ecodist_to_regions <- points3 %>%
  distinct(Ecds_ID, Ecrgn_n)

# Check whether an ecodistrict belongs to multiple regions/provinces
ecodist_to_regions %>%
  count(Ecds_ID) %>%
  filter(n > 1)

#  Ecds_ID     n
#1     618     2
#2     815     2
#3     828     2
#4     833     2
#5     838     2

points3 %>%
  distinct(Ecds_ID, Ecrgn_n) %>%
  filter(Ecds_ID %in%
           (points3 %>%
              distinct(Ecds_ID, Ecrgn_n) %>%
              count(Ecds_ID) %>%
              filter(n > 1) %>%
              pull(Ecds_ID)))

#      Ecds_ID Ecrgn_n               
#1     618 Western Alberta Upland
#2     833 Mixed Grassland       
#3     838 Cypress Upland        
#4     815 Mixed Grassland       
#5     828 Mixed Grassland       
#6     833 Cypress Upland        
#7     815 Cypress Upland        
#8     828 Moist Mixed Grassland 
#9     618 Peace Lowland         
#10    838 Mixed Grassland  

# Correct ecoregions are:
# 618 = 145 = Western Alberta Upland
# 815 = 159 = Mixed Grassland 
# 828 = 159 = Mixed Grassland 
# 833 = 159 = Mixed Grassland 
# 838 = 160 = Cypress Upland

# Make corrections
points3 <- points3 %>%
  mutate(Ecrgn_n = case_when(
    Ecds_ID == "618" ~ "Western Alberta Upland",
    Ecds_ID == "815" ~ "Mixed Grassland",
    Ecds_ID == "828" ~ "Mixed Grassland",
    Ecds_ID == "833" ~ "Mixed Grassland",
    Ecds_ID == "838" ~ "Cypress Upland",
    TRUE ~ Ecrgn_n
  ))

# Get mapping (one row per ecodistrict)
ecodist_to_regions <- points3 %>%
  distinct(Ecds_ID, Ecrgn_n)

# Check whether an ecodistrict belongs to multiple regions/provinces
ecodist_to_regions %>%
  count(Ecds_ID) %>%
  filter(n > 1)

# export
writexl::write_xlsx(points3, file.path(out_dir, "point_raster_extract_v2_native_v_tame_clean.xlsx"))
