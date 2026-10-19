
# validate the binary and triclass spatial agreement analysis results (pilot)

# load packages
library(terra)
library(sf)
library(tidyverse)

# Set the temporary directory to a drive with sufficient space
terraOptions(tempdir = "R:/CanadianGrasslands/WG_classification/scratch")

####--- DIRECTORIES ---####

# pilot area
pa_dir <- "R:/CanadianGrasslands/processeddatasets/pilotarea"

# ecodistricts
ecodist_dir <- "R:/CanadianGrasslands/processeddatasets/boundaries/national-ecological-framework-for-canada/Ecodistricts"

# ccs
ccs_dir <- "R:/CanadianGrasslands/processeddatasets/boundaries/census-consolidated-subdivisions"

# processed data directory
proc_data_dir <- "R:/CanadianGrasslands/WG_Classification/processed"

# data directories
aci_dir <- file.path(proc_data_dir, "aafc_aci/aligned")
luts30_dir <- file.path(proc_data_dir, "aafc_luts_2020_30m_v5/masked")
luts10_dir <- file.path(proc_data_dir, "aafc_luts_2020_10m_v1/masked")
cec_dir  <- file.path(proc_data_dir, "cec_nalcms/masked")
ngi_dir  <- file.path(proc_data_dir, "cfga_ngi/aligned")
cgr_dir  <- file.path(proc_data_dir, "cgr_gam_adv/masked")
eccc_dir <- file.path(proc_data_dir, "eccc_grass/masked")
esa_dir  <- file.path(proc_data_dir, "esa_worldcover/masked")
glc_dir  <- file.path(proc_data_dir, "glc_fcs30d/masked")
cnvc_dir <- file.path(proc_data_dir, "natureserve_cnvc/aligned")
lcc_dir  <- file.path(proc_data_dir, "nrcan_lcc/masked")
pudl_dir <- file.path(proc_data_dir, "phjv_pudl/aligned")
lcr_dir  <- file.path(proc_data_dir, "statcan_lcr/masked")
gpw_dir  <- file.path(proc_data_dir, "wri_gpw_ann/masked")
plow_dir <- file.path(proc_data_dir, "wwf_plowprint/masked")
pli_dir  <- file.path(proc_data_dir, "sk_pli/aligned")
plvi_dir <- file.path(proc_data_dir, "ab_plvi/aligned")
hfi_dir  <- file.path(proc_data_dir, "abmi_hfi/aligned")
dw_dir   <- file.path(proc_data_dir, "dw_2025_grass/masked")
esri_dir <- file.path(proc_data_dir, "esri_rangelands/masked")
mgi_dir  <- file.path(proc_data_dir, "mb_mgi/aligned")
gcc_dir  <- file.path(proc_data_dir, "bc_gcc/aligned")
pglc_dir <- file.path(proc_data_dir, "aafc_pglc/aligned")

####--- FILES ---####

# Pilot area
pa_file <- file.path(pa_dir, "pilot_area_aligned_cnvc_Jan2026.tif")

# ecodistricts
ecodist_file <- file.path(ecodist_dir, "ecodistricts_clean_rast_Oct2026.tif")

# ccs
ccs_file <- file.path(ccs_dir, "ccs_aligned_class_sa_rast.tif")

# Study area
#sa_file <- file.path(sa_dir, "class_study_area_aligned_cnvc_Dec2025.tif")

# binary rasters
aci2024_bin_file  <- file.path(aci_dir, "aafc_aci_merged_reproj_algn_binary_noshrub.tif")
aci2025_bin_file  <- file.path(aci_dir, "aafc_aci_2025_merged_reproj_algn_binary_noshrub.tif")
luts30_bin_file <- file.path(luts30_dir, "luts30_reproj_algn_binary_masked.tif")
luts10_bin_file <- file.path(luts10_dir, "luts10_reproj_algn_binary_masked.tif")
cec_bin_file  <- file.path(cec_dir, "nalcms_reproj_algn_binary_noshrub_masked.tif")
ngi_bin_file  <- file.path(ngi_dir, "cfga_ngi_merged_reproj_algn_binary.tif")
cgr_bin_file  <- file.path(cgr_dir, "cgr_reproj_algn_binary_masked.tif")
eccc_bin_file <- file.path(eccc_dir, "eccc_reclass_reproj_algn_binary_masked.tif")
cnvc_bin_file <- file.path(cnvc_dir, "cnvc_algn_binary.tif")
lcc_bin_file  <- file.path(lcc_dir, "lcc_reproj_algn_binary_noshrub_masked.tif")
pudl_bin_file <- file.path(pudl_dir, "pudl_reclass_reproj_algn_binary_noshrub.tif")
lcr_bin_file  <- file.path(lcr_dir, "lcr_reproj_algn_binary_masked.tif")
gpw_bin_file  <- file.path(gpw_dir, "gpw_reproj_algn_binary_noshrub_masked.tif")
plow_bin_file <- file.path(plow_dir, "plow_reproj_algn_intact_masked.tif")
pli_bin_file  <- file.path(pli_dir, "sk_pli_merged_reproj_algn_binary_noshrub.tif")
esa_bin_file  <- file.path(esa_dir, "esa_reproj_algn_binary_noshrub_masked.tif")
glc_bin_file  <- file.path(glc_dir, "glc_reproj_algn_binary_noshrub_masked.tif")
plvi_bin_file <- file.path(plvi_dir, "plvi_reproj_algn_binary_gte2015.tif")
hfi_bin_file  <- file.path(hfi_dir, " hfi_reproj_algn_binary_gte2015.tif")
dw_bin_file   <- file.path(dw_dir, " dw_reproj_algn_masked_binary.tif")
esri_bin_file <- file.path(esri_dir, "esri_reproj_algn_masked_binary.tif")
mgi_bin_file  <- file.path(mgi_dir, "mb_mgi_reproj_algn_binary_noshrub.tif")
gcc_bin_file  <- file.path(gcc_dir, "gcc_reproj_algn_binary.tif")
pglc_bin_file <- file.path(pglc_dir, "pglc_reproj_algn_binary.tif")

# triclass rasters
aci2024_nt_file  <- file.path(aci_dir,  "aafc_aci_merged_reproj_algn_natvtam_noshrub.tif")
aci2025_nt_file  <- file.path(aci_dir,  "aafc_aci_2025_merged_reproj_algn_natvtam_noshrub.tif")
ngi_nt_file  <- file.path(ngi_dir,  "cfga_ngi_merged_reproj_algn_natvtam_good.tif")
cgr_nt_file  <- file.path(cgr_dir,   "cgr_reproj_algn_natvtam_masked.tif")
eccc_nt_file <- file.path(eccc_dir,  "eccc_reclass_reproj_algn_natvtam_masked.tif")
cnvc_nt_file <- file.path(cnvc_dir,  "cnvc_algn_natvtam.tif")
pudl_nt_file <- file.path(pudl_dir, "pudl_reclass_reproj_algn_natvtam_noshrub.tif")
pli_nt_file  <- file.path(pli_dir,  "sk_pli_merged_reproj_algn_natvtam_noshrub.tif")
gpw_nt_file  <- file.path(gpw_dir,   "gpw_reproj_algn_natvtam_noshrub_masked.tif")
plvi_nt_file <- file.path(plvi_dir, "plvi_reproj_algn_natvtam.tif")
hfi_nt_file  <- file.path(hfi_dir,  "hfi_reproj_algn_binary.tif") # tame only - same as binary
mgi_nt_file  <- file.path(mgi_dir, "mb_mgi_reproj_algn_natvtam_noshrub.tif")
pglc_nt_file <- file.path(pglc_dir, "pglc_reproj_algn_natvtam.tif")

################## BINARY CHECK - ALL ECODISTRICTS ########################

####--- LOAD RASTERS ---####

pilot <- rast(pa_file)

r1 <- rast(aci2024_bin_file)
r2 <- rast(aci2025_bin_file)

####--- CROP RASTERS TO PILOT ---####

r1 <- mask(crop(r1, pilot), pilot)
r2 <- mask(crop(r2, pilot), pilot)

plot(r1)
plot(r2)

####--- CALCULATE CONFUSION COUNTS ---####

# extract values
v1 <- values(r1, mat = FALSE)
v2 <- values(r2, mat = FALSE)

valid <- !is.na(v1) & !is.na(v2)

v1 <- v1[valid]
v2 <- v2[valid]

# convert to binary presence
b1 <- v1 > 0
b2 <- v2 > 0

# calculate confusion counts
n11 <- sum(b1 & b2)
n10 <- sum(b1 & !b2)
n01 <- sum(!b1 & b2)
n00 <- sum(!b1 & !b2)

c(
  n11 = n11,
  n10 = n10,
  n01 = n01,
  n00 = n00
)

# verify global counts
table(
  r1 = as.integer(b1),
  r2 = as.integer(b2)
)

#         r2
# r1       0       1
#  0  542177   95668
#  1   23937 1061802

####--- CALCULATE METRICS ---####

n11 <- as.numeric(n11)
n10 <- as.numeric(n10)
n01 <- as.numeric(n01)
n00 <- as.numeric(n00)

# Jaccard distance
d_jaccard <- 1 - n11 / (n11 + n10 + n01)
d_jaccard # = 0.1012395

# F1
precision <- n11 / (n11 + n10)
recall    <- n11 / (n11 + n01)

f1 <- 2 * precision * recall /
  (precision + recall)
f1 # = 0.9466813

# MCC
mcc <- (n11 * n00 - n10 * n01) /
  sqrt(
    (n11 + n10) *
      (n11 + n01) *
      (n00 + n10) *
      (n00 + n01)
  )
mcc # = 0.85119

################## BINARY CHECK - ONE ECODISTRICT ##################

####--- LOAD RASTERS ---####

pilot <- rast(pa_file)

r1 <- rast(aci2024_bin_file)
r2 <- rast(aci2025_bin_file)

####--- CROP RASTERS TO PILOT ---####

r1 <- mask(crop(r1, pilot), pilot)
r2 <- mask(crop(r2, pilot), pilot)

plot(r1)
plot(r2)

####--- CROP RASTERS TO ONE ECODISTRICT ---####

ecodist <- rast(ecodist_file)
ecodist

ecodist <- crop(ecodist, pilot)
ecodist <- mask(ecodist, pilot)

ccs1 <- ecodist == 832

r1_832 <- mask(r1, ccs1, maskvalues = 0)
r2_832 <- mask(r2, ccs1, maskvalues = 0)

plot(r1_832)
plot(r2_832)

####--- CALCULATE CONFUSION COUNTS ---####

# extract values
v1 <- values(r1_832, mat = FALSE)
v2 <- values(r2_832, mat = FALSE)

valid <- !is.na(v1) & !is.na(v2)

v1 <- v1[valid]
v2 <- v2[valid]

# convert to binary presence
b1 <- v1 > 0
b2 <- v2 > 0

# calculate confusion counts
n11 <- sum(b1 & b2)
n10 <- sum(b1 & !b2)
n01 <- sum(!b1 & b2)
n00 <- sum(!b1 & !b2)

c(
  n11 = n11,
  n10 = n10,
  n01 = n01,
  n00 = n00
)

# verify global counts
table(
  r1 = as.integer(b1),
  r2 = as.integer(b2)
)

#        r2
# r1      0      1
#  0 148656  72914
#  1  11099 436582

####--- CALCULATE METRICS ---####

n11 <- as.numeric(n11)
n10 <- as.numeric(n10)
n01 <- as.numeric(n01)
n00 <- as.numeric(n00)

# Jaccard distance
d_jaccard <- 1 - n11 / (n11 + n10 + n01)
d_jaccard # = 0.1613788

# F1
precision <- n11 / (n11 + n10)
recall    <- n11 / (n11 + n01)

f1 <- 2 * precision * recall /
  (precision + recall)
f1 # = 0.9122284

# MCC
mcc <- (n11 * n00 - n10 * n01) /
  sqrt(
    (n11 + n10) *
      (n11 + n01) *
      (n00 + n10) *
      (n00 + n01)
  )
mcc # = 0.7132822

################## BINARY CHECK - ONE CCS ##################

####--- LOAD RASTERS ---####

pilot <- rast(pa_file)

r1 <- rast(aci2024_bin_file)
r2 <- rast(aci2025_bin_file)

####--- CROP RASTERS TO PILOT ---####

r1 <- mask(crop(r1, pilot), pilot)
r2 <- mask(crop(r2, pilot), pilot)

plot(r1)
plot(r2)

####--- CROP RASTERS TO ONE CCS ---####

ccs <- rast(ccs_file)
ccs
plot(ccs)
levels(ccs)
file.info(ccs_file)$size

ccs_lookup <- levels(ccs)[[1]]

write_csv(ccs_lookup, file.path(ccs_dir, "ccs_aligned_class_sa_rast_lookup.csv"))

ccs <- crop(ccs, pilot)
ccs <- mask(ccs, pilot)

ccs1 <- ccs == 832

r1_832 <- mask(r1, ccs1, maskvalues = 0)
r2_832 <- mask(r2, ccs1, maskvalues = 0)

plot(r1_832)
plot(r2_832)

####--- CALCULATE CONFUSION COUNTS ---####

# extract values
v1 <- values(r1_832, mat = FALSE)
v2 <- values(r2_832, mat = FALSE)

valid <- !is.na(v1) & !is.na(v2)

v1 <- v1[valid]
v2 <- v2[valid]

# convert to binary presence
b1 <- v1 > 0
b2 <- v2 > 0

# calculate confusion counts
n11 <- sum(b1 & b2)
n10 <- sum(b1 & !b2)
n01 <- sum(!b1 & b2)
n00 <- sum(!b1 & !b2)

c(
  n11 = n11,
  n10 = n10,
  n01 = n01,
  n00 = n00
)

# verify global counts
table(
  r1 = as.integer(b1),
  r2 = as.integer(b2)
)


####--- CALCULATE METRICS ---####

n11 <- as.numeric(n11)
n10 <- as.numeric(n10)
n01 <- as.numeric(n01)
n00 <- as.numeric(n00)

# Jaccard distance
d_jaccard <- 1 - n11 / (n11 + n10 + n01)
d_jaccard # = 

# F1
precision <- n11 / (n11 + n10)
recall    <- n11 / (n11 + n01)

f1 <- 2 * precision * recall /
  (precision + recall)
f1 # = 

# MCC
mcc <- (n11 * n00 - n10 * n01) /
  sqrt(
    (n11 + n10) *
      (n11 + n01) *
      (n00 + n10) *
      (n00 + n01)
  )
mcc # = 

################## TRI-CLASS CHECK - ALL ECODISTRICTS ###################

####--- LOAD RASTERS ---####

pilot <- rast(pa_file)

r1 <- rast(aci2024_nt_file)
r2 <- rast(aci2025_nt_file)

####--- CROP RASTERS TO PILOT ---####

r1 <- mask(crop(r1, pilot), pilot)
r2 <- mask(crop(r2, pilot), pilot)

plot(r1)
plot(r2)

####--- CALCULATE CONFUSION COUNTS ---####

# extract values
v1 <- values(r1, mat = FALSE)
v2 <- values(r2, mat = FALSE)

valid <- !is.na(v1) & !is.na(v2)

v1 <- v1[valid]
v2 <- v2[valid]

# build the confusion matrix
classes <- c(0, 1, 2)

C <- table(
  factor(v1, levels = classes),
  factor(v2, levels = classes)
)

C

#        0      1      2
# 0 542177   9868  85800
# 1  14885  76087  60444
# 2   9052   1385 923886

# verify row totals
row0 <- sum(C[1,])
row1 <- sum(C[2,])
row2 <- sum(C[3,])

row0 # = 637845
row1 # = 151416
row2 # = 934323

# verify column totals
col0 <- sum(C[,1])
col1 <- sum(C[,2])
col2 <- sum(C[,3])

col0 # = 566114
col1 # = 87340
col2 # = 1070130

####--- CALCULATE METRICS ---####

# per-class Jaccard distance
TP0 <- C[1,1]
FN0 <- sum(C[1,]) - TP0
FP0 <- sum(C[,1]) - TP0
jaccard_class0 <- TP0 / (TP0 + FP0 + FN0)
jaccard_class0 # = 0.8192683

TP1 <- C[2,2]
FN1 <- sum(C[2,]) - TP1
FP1 <- sum(C[,2]) - TP1
jaccard_class1 <- TP1 / (TP1 + FP1 + FN1)
jaccard_class1 # = 0.4677412

TP2 <- C[3,3]
FN2 <- sum(C[3,]) - TP2
FP2 <- sum(C[,3]) - TP2
jaccard_class2 <- TP2 / (TP2 + FP2 + FN2)
jaccard_class2 # = 0.8550011

# per-class F1
precision0 <- TP0 / (TP0 + FP0)
recall0 <- TP0 / (TP0 + FN0)
f1_class0 <- 2 * precision0 * recall0 /
  (precision0 + recall0)
f1_class0 # = 0.9006569

precision1 <- TP1 / (TP1 + FP1)
recall1 <- TP1 / (TP1 + FN1)
f1_class1 <- 2 * precision1 * recall1 /
  (precision1 + recall1)
f1_class1 # = 0.637362

precision2 <- TP2 / (TP2 + FP2)
recall2 <- TP2 / (TP2 + FN2)
f1_class2 <- 2 * precision2 * recall2 /
  (precision2 + recall2)
f1_class2 # = 0.9218335

# macro F1
f1_macro <- mean(c(f1_class0, f1_class1, f1_class2), na.rm = TRUE)
f1_macro # = 0.8199508

# multiclass MCC
C <- matrix(as.numeric(C), nrow = 3, dimnames = dimnames(C))
t_sum <- sum(C)
sum_row_prod <- sum(rowSums(C) * colSums(C))
trace_C <- sum(diag(C))
numerator <- (trace_C * t_sum) - sum_row_prod
denom_row <- t_sum^2 - sum(colSums(C)^2)
denom_col <- t_sum^2 - sum(rowSums(C)^2)
mcc_multi <- numerator / sqrt(denom_row * denom_col)
mcc_multi # = 0.8123484

################## TRI-CLASS CHECK - ONE ECODISTRICT ###################

####--- LOAD RASTERS ---####

pilot <- rast(pa_file)

r1 <- rast(aci2024_nt_file)
r2 <- rast(aci2025_nt_file)

####--- CROP RASTERS TO PILOT ---####

r1 <- mask(crop(r1, pilot), pilot)
r2 <- mask(crop(r2, pilot), pilot)

plot(r1)
plot(r2)

####--- CROP RASTERS TO ONE ECODISTRICT ---####

ecodist <- rast(ecodist_file)
ecodist

ecodist <- crop(ecodist, pilot)
ecodist <- mask(ecodist, pilot)

ccs1 <- ecodist == 832

r1_832 <- mask(r1, ccs1, maskvalues = 0)
r2_832 <- mask(r2, ccs1, maskvalues = 0)

plot(r1_832)
plot(r2_832)

####--- CALCULATE CONFUSION COUNTS ---####

# extract values
v1 <- values(r1_832, mat = FALSE)
v2 <- values(r2_832, mat = FALSE)

valid <- !is.na(v1) & !is.na(v2)

v1 <- v1[valid]
v2 <- v2[valid]

# build the confusion matrix
classes <- c(0, 1, 2)

C <- table(
  factor(v1, levels = classes),
  factor(v2, levels = classes)
)

C

#        0      1      2
# 0 148656   7305  65609
# 1   6434  38017  23151
# 2   4665    675 374739

# verify row totals
row0 <- sum(C[1,])
row1 <- sum(C[2,])
row2 <- sum(C[3,])

row0 # = 221570
row1 # = 67602
row2 # = 380079

# verify column totals
col0 <- sum(C[,1])
col1 <- sum(C[,2])
col2 <- sum(C[,3])

col0 # = 159755
col1 # = 45997
col2 # = 463499

####--- CALCULATE METRICS ---####

# per-class Jaccard distance
TP0 <- C[1,1]
FN0 <- sum(C[1,]) - TP0
FP0 <- sum(C[,1]) - TP0
jaccard_class0 <- TP0 / (TP0 + FP0 + FN0)
jaccard_class0 # = 0.6389162

TP1 <- C[2,2]
FN1 <- sum(C[2,]) - TP1
FP1 <- sum(C[,2]) - TP1
jaccard_class1 <- TP1 / (TP1 + FP1 + FN1)
jaccard_class1 # = 0.5029901

TP2 <- C[3,3]
FN2 <- sum(C[3,]) - TP2
FP2 <- sum(C[,3]) - TP2
jaccard_class2 <- TP2 / (TP2 + FP2 + FN2)
jaccard_class2 # = 0.7992914

# per-class F1
precision0 <- TP0 / (TP0 + FP0)
recall0 <- TP0 / (TP0 + FN0)
f1_class0 <- 2 * precision0 * recall0 /
  (precision0 + recall0)
f1_class0 # = 0.7796814

precision1 <- TP1 / (TP1 + FP1)
recall1 <- TP1 / (TP1 + FN1)
f1_class1 <- 2 * precision1 * recall1 /
  (precision1 + recall1)
f1_class1 # = 0.6693193

precision2 <- TP2 / (TP2 + FP2)
recall2 <- TP2 / (TP2 + FN2)
f1_class2 <- 2 * precision2 * recall2 /
  (precision2 + recall2)
f1_class2 # = 0.8884513

# macro F1
f1_macro <- mean(c(f1_class0, f1_class1, f1_class2), na.rm = TRUE)
f1_macro # = 0.7791507

# multiclass MCC
C <- matrix(as.numeric(C), nrow = 3, dimnames = dimnames(C))
t_sum <- sum(C)
sum_row_prod <- sum(rowSums(C) * colSums(C))
trace_C <- sum(diag(C))
numerator <- (trace_C * t_sum) - sum_row_prod
denom_row <- t_sum^2 - sum(colSums(C)^2)
denom_col <- t_sum^2 - sum(rowSums(C)^2)
mcc_multi <- numerator / sqrt(denom_row * denom_col)
mcc_multi # = 0.710993
