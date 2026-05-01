
###################################################################################
##                    ANALYSIS: SPATIAL AGREEMENT SUMMARIES                      ##
###################################################################################

##----------------------- Alison Long - February 2026 -------------------------

# step 5: run descriptive data summaries to assess similarities and differences
#        among existing grassland datasets. This version handles all grassland
#        data products regardless of extent.

# load packages
library(terra)
library(ggplot2)
library(pheatmap)
library(vegan)
library(dplyr)
library(tidyr)
library(readr)
library(ggrepel)
library(ape)
library(sf)
library(viridis)
library(ggspatial)  # for scalebar / north arrow 
library(stringr)
library(patchwork)
library(lme4)
library(broom.mixed)

# Set the temporary directory to a drive with sufficient space
terraOptions(tempdir = "R:/CanadianGrasslands/WG_classification/scratch")

# ============================================================
# Pairwise similarity/difference for binary rasters (0/1) with
# varying extents/coverage. Comparisons are computed only on the
# intersection of valid (non-NA) pixels for each pair.
#
# Outputs:
# - Pairwise distance matrices: Jaccard, Sorensen/Dice, (1-Kappa)
# - Pairwise overlap matrices: n_overlap, frac_overlap_i, frac_overlap_j
# - PCoA plots + clustering + heatmaps
# - Optional union-grid coverage + consensus + uncertainty maps
#
# Notes:
# - Jaccard/Sorensen are presence-focused, avoiding dominance by shared zeros. [1](https://arxiv.org/html/2405.06559v1)
# - PCoA is appropriate for ordinating objects using arbitrary distances (e.g., Jaccard). [3](https://fundy.ext.unb.ca/metrics.htm)[4](https://researchportal.helsinki.fi/en/publications/measuring-habitat-fragmentation-an-evaluation-of-landscape-patter/)
# - Cohen's kappa measures chance-corrected agreement for categorical labels. [2](https://github.com/andland/logisticPCA)
# ============================================================

# ----------------------------
# DEFINE SETTINGS AND INPUTS 
# ----------------------------
#RASTER_DIR   <- "R:/CanadianGrasslands/WG_Classification/github/tomorrows-prairies/processeddatasets/Grass/full_binary"
#RASTER_GLOB  <- "\\.tif$"
META_CSV     <- "R:/CanadianGrasslands/WG_Classification/github/tomorrows-prairies/processeddatasets/Grass/metadata_full_extent_binary.csv"  # columns: raster,year,model,source (year optional)
OUT_DIR      <- "R:/CanadianGrasslands/WG_Classification/results/outputs_pairwise_partial_extent" 

####--- DIRECTORIES ---####

# study area
sa_dir <- "R:/CanadianGrasslands/processeddatasets/studyarea"

# grasslands
raw_grass_data_dir <- "R:/CanadianGrasslands/existingdatasets"
proc_grass_data_dir <- "R:/CanadianGrasslands/WG_Classification/processed"
#out_dir <- "R:/CanadianGrasslands/WG_Classification/results/accuracy_assessment"

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
cnvc_masked_dir   <- file.path(proc_grass_data_dir, "natureserve_cnvc/masked")
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

####--- FILES ---####

# Study area
#sa_file <- file.path(sa_dir, "class_study_area_aligned_cnvc_Dec2025.tif")

# grass raster files
aci_file    <- file.path(aci_aligned_dir, "aafc_aci_merged_reproj_algn_binary_noshrub.tif")
luts30_file <- file.path(luts30_masked_dir, "luts30_reproj_algn_binary_masked.tif")
luts10_file <- file.path(luts10_masked_dir, "luts10_reproj_algn_binary_masked.tif")
cec_file    <- file.path(cec_masked_dir, "nalcms_reproj_algn_binary_noshrub_masked.tif")
ngi_file    <- file.path(ngi_aligned_dir, "cfga_ngi_merged_reproj_algn_binary.tif")
cgr_file    <- file.path(cgr_masked_dir, "cgr_reproj_algn_binary_masked.tif")
eccc_file   <- file.path(eccc_masked_dir, "eccc_reclass_reproj_algn_binary_masked.tif")
esa_file    <- file.path(esa_masked_dir, "esa_reproj_algn_binary_noshrub_masked.tif")
glc_file    <- file.path(glc_masked_dir, "glc_reproj_algn_binary_noshrub_masked.tif")
cnvc_file   <- file.path(cnvc_masked_dir, "cnvc_algn_binary_masked.tif")
lcc_file    <- file.path(lcc_masked_dir, "lcc_reproj_algn_binary_noshrub_masked.tif")
pudl_file   <- file.path(pudl_aligned_dir, "pudl_reclass_reproj_algn_binary_noshrub.tif")
pli_file    <- file.path(pli_aligned_dir, "sk_pli_merge_reproj_algn_binary_noshrub.tif")
#pli_cu_file <- file.path(pli_aligned_dir, "sk_pli_cu_reproj_algn_binary_noshrub.tif")
#pli_mg_file <- file.path(pli_aligned_dir, "sk_pli_mg_reproj_algn_binary_noshrub.tif")
#pli_mm_file <- file.path(pli_aligned_dir, "sk_pli_mm_reproj_algn_binary_noshrub.tif")
lcr_file    <- file.path(lcr_masked_dir, "lcr_reproj_algn_binary_masked.tif")
gpw_file    <- file.path(gpw_masked_dir, "gpw_reproj_algn_binary_noshrub_masked.tif")
plow_file   <- file.path(plow_masked_dir, "plow_reproj_algn_intact_masked.tif")
plvi_file   <- file.path(plvi_aligned_dir, "plvi_reproj_algn_binary_gte2015.tif")
hfi_file    <- file.path(hfi_aligned_dir, "hfi_reproj_algn_binary_gte2015.tif")
gcc_file    <- file.path(gcc_aligned_dir, "gcc_reproj_algn_binary.tif")
dw_file     <- file.path(dw_masked_dir, "dw_reproj_algn_masked_binary.tif")
esri_file   <- file.path(esri_masked_dir, "esri_reproj_algn_masked_binary.tif")
mgi_file    <- file.path(mgi_aligned_dir, "mb_mgi_reproj_algn_binary_noshrub.tif")

# Sampling within each pairwise overlap to speed up large rasters:
# - NULL = use all overlapping pixels (may be slow/large)
# - integer = randomly sample up to this many overlapping pixels per pair
SAMPLE_N_PER_PAIR <- 100000

# If two rasters do not overlap at all, set their distance to this value:
NO_OVERLAP_DISTANCE <- 1

# If overlap exists but there are no presences in either raster within overlap
# (TP+FP+FN == 0), Jaccard is undefined; choose how to handle:
# - "zero" = treat as identical (distance 0)
# - "na"   = set NA
JACCARD_EMPTY_BEHAVIOR <- "zero"

# Optional: build union-grid coverage/consensus maps across all rasters
MAKE_UNION_CONSENSUS <- TRUE

SEED <- 42
set.seed(SEED)
#dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# -----------------------
# ECODISTRICT SUMMARY SETTINGS
# -----------------------
ECODISTRICT_SHP   <- "R:/CanadianGrasslands/processeddatasets/boundaries/ecodistricts_clean.shp"  
ECODISTRICT_LAYER <- NULL                        # "ecodistricts" for gpkg; NULL for shapefile
ECODISTRICT_ID    <- "ECODISTRIC"                    # unique ID field in ecodistrict polygons

# Agreement is only meaningful where >= this many rasters have non-NA values
MIN_RASTERS_FOR_AGREEMENT <- 2

# Per-ecodistrict sampling within each raster-pair overlap 
SAMPLE_N_PER_PAIR_PER_ECO <- NULL   # NULL = use all overlapping cells per pair per ecodistrict

# Minimum number of overlapping (non-NA in both) cells needed to compute metrics
MIN_OVERLAP_CELLS_PER_PAIR <- 50

# Output files
OUT_ECODISTRICT_CSV <- file.path(OUT_DIR, "ecodistrict_agreement_disagreement_summary.csv")
OUT_ECODISTRICT_GPKG <- file.path(OUT_DIR, "ecodistrict_agreement_disagreement_summary.gpkg")
WRITE_ECODISTRICT_GPKG <- TRUE

OUT_ECO_PAIRWISE_LONG  <- file.path(OUT_DIR, "ecodistrict_pairwise_metrics_long.csv")
OUT_ECO_PAIRWISE_SUM   <- file.path(OUT_DIR, "ecodistrict_pairwise_summary.csv")

# -----------------------
# READ METADATA
# -----------------------
# metadata.csv: raster,year,model,source
meta <- read_csv(META_CSV, show_col_types = FALSE) %>%
  mutate(
    model  = as.factor(model),
    source = as.factor(source),
    year = if ("year" %in% names(.)) as.POSIXct(year, tz = "UTC") else NA
  )

# -----------------------
# LOAD RASTERS
# -----------------------

####--- READ INTO TERRA ---####

files <- c(
  aci_file, luts10_file, luts30_file, cec_file, ngi_file, cgr_file, eccc_file,
  esa_file, glc_file, cnvc_file, lcc_file, pudl_file, pli_file,
  #pli_cu_file, pli_mg_file, pli_mm_file, 
  lcr_file, gpw_file, plow_file,
  plvi_file, hfi_file, gcc_file, dw_file, esri_file, mgi_file
)

keys <- c(
  "aci","luts10","luts30","cec","ngi","cgr","eccc", "esa","glc",
  "cnvc","lcc","pudl","pli",
  #"pli_cu","pli_mg","pli_mm",
  "lcr","gpw","plow","plvi","hfi","gcc","dw","esri","mgi"
)

names(files) <- keys
rlist <- lapply(files, rast)  # each is single-layer SpatRaster

#files <- list.files(RASTER_DIR, pattern = RASTER_GLOB, full.names = TRUE)
#if (length(files) == 0) stop("No rasters found in RASTER_DIR matching RASTER_GLOB.")

files_df <- tibble(
  file = files,
  raster = basename(files),
  raster_nox = tools::file_path_sans_ext(basename(files))
)

meta2 <- meta %>%
  left_join(files_df, by = "raster") %>%
  mutate(file = ifelse(is.na(file),
                       files_df$file[match(raster, files_df$raster_nox)],
                       file)) %>%
  mutate(file = as.character(file))

if (any(is.na(meta2$file))) {
  missing <- meta2$raster[is.na(meta2$file)]
  stop("Could not match these metadata rasters to files: ", paste(missing, collapse = ", "))
}

rlist <- lapply(meta2$file, terra::rast)
names(rlist) <- meta2$raster

# -----------------------
# BASIC CHECKS / STANDARDIZE VALUES
# -----------------------
# Convert to integer 0/1, keep NA as NA
rlist <- lapply(rlist, function(r) {
  r <- round(r)
  r[r < 0] <- 0
  r[r > 1] <- 1
  r
})

# Integrity table
safe_ext2 <- function(r) {
  e <- tryCatch(terra::ext(r), error = function(e) NULL)
  if (is.null(e)) return(c(xmin=NA_real_, xmax=NA_real_, ymin=NA_real_, ymax=NA_real_))
  c(xmin = terra::xmin(e), xmax = terra::xmax(e), ymin = terra::ymin(e), ymax = terra::ymax(e))
}

# Build a simple data.frame by iterating directly over rlist (no rowwise/which)
check2 <- lapply(seq_along(rlist), function(i) {
  r <- rlist[[i]]
  ex <- safe_ext2(r)
  data.frame(
    name    = names(rlist)[i],
    nlyr    = terra::nlyr(r),
    nrow    = terra::nrow(r),
    ncol    = terra::ncol(r),
    ncell   = terra::ncell(r),
    has_crs = !is.na(terra::crs(r)),
    xmin    = ex["xmin"],
    xmax    = ex["xmax"],
    ymin    = ex["ymin"],
    ymax    = ex["ymax"],
    stringsAsFactors = FALSE
  )
})
check2 <- dplyr::bind_rows(check2)
print(check2)

bad2 <- check2 %>% filter(is.na(xmin) | is.na(xmax) | is.na(ymin) | is.na(ymax) | ncell == 0)
if (nrow(bad2) > 0) {
  warning("Rasters with NA extent or zero cells: ", paste(bad2$name, collapse = ", "))
}

# look at one raster
r <- terra::rast("R:/CanadianGrasslands/WG_Classification/github/tomorrows-prairies/processeddatasets/Grass/full_binary/aci_merged_reproj_algn_binary_noshrub.tif")
terra::crs(r)      # shows a CRS 
terra::res(r)      # 
terra::ext(r)      # 
terra::geom(r)     # list of georeferencing; expected to give error

# ------------------------------------------------------------
# Helper: confusion counts & derived metrics on two vectors
# ------------------------------------------------------------
pair_metrics_from_vectors <- function(a, b, jaccard_empty = c("zero","na")) {
  jaccard_empty <- match.arg(jaccard_empty)
  
  # Expect a,b are 0/1 with no NA
  TP <- sum(a == 1 & b == 1)
  FP <- sum(a == 1 & b == 0)
  FN <- sum(a == 0 & b == 1)
  TN <- sum(a == 0 & b == 0)
  N  <- length(a)
  
  # Jaccard distance (presence-focused) [1](https://arxiv.org/html/2405.06559v1)
  denomJ <- TP + FP + FN
  if (denomJ == 0) {
    d_jac <- if (jaccard_empty == "zero") 0 else NA_real_
  } else {
    d_jac <- 1 - (TP / denomJ)
  }
  
  # Sørensen/Dice distance (presence-focused)
  denomS <- (2 * TP + FP + FN)
  if (denomS == 0) {
    d_sor <- 0
  } else {
    d_sor <- 1 - (2 * TP / denomS)
  }
  
  # Cohen's kappa (agreement-style) [2](https://github.com/andland/logisticPCA)
  po <- (TP + TN) / N
  p_yes_a <- (TP + FP) / N
  p_no_a  <- (FN + TN) / N
  p_yes_b <- (TP + FN) / N
  p_no_b  <- (FP + TN) / N
  pe <- p_yes_a * p_yes_b + p_no_a * p_no_b
  kappa <- if ((1 - pe) == 0) NA_real_ else (po - pe) / (1 - pe)
  
  list(
    TP = TP, FP = FP, FN = FN, TN = TN, N = N,
    d_jaccard = d_jac,
    d_sorensen = d_sor,
    kappa = kappa
  )
}

# ------------------------------------------------------------
# Helper: extract overlapping pixel vectors for a pair (with explicit terra calls)
# ------------------------------------------------------------
extract_pair_overlap <- function(rA, rB, sample_n = NULL) {
  # Safe extents
  eA <- tryCatch(terra::ext(rA), error = function(e) NULL)
  eB <- tryCatch(terra::ext(rB), error = function(e) NULL)
  if (is.null(eA) || is.null(eB)) return(list(ok = FALSE))
  
  # Quick bbox overlap test
  if (terra::xmin(eA) >= terra::xmax(eB) || terra::xmin(eB) >= terra::xmax(eA) ||
      terra::ymin(eA) >= terra::ymax(eB) || terra::ymin(eB) >= terra::ymax(eA)) {
    return(list(ok = FALSE))
  }
  
  # Use terra::intersect (NOT dplyr::intersect)
  e <- tryCatch(terra::intersect(eA, eB), error = function(e) NULL)
  if (is.null(e)) return(list(ok = FALSE))
  
  # Crop to overlap
  a <- tryCatch(terra::crop(rA, e, snap = "near"), error = function(e) NULL)
  b <- tryCatch(terra::crop(rB, e, snap = "near"), error = function(e) NULL)
  if (is.null(a) || is.null(b)) return(list(ok = FALSE))
  
  # Align grids if needed
  if (!terra::compareGeom(a, b, stopOnError = FALSE)) {
    b <- tryCatch(terra::resample(b, a, method = "near"), error = function(e) NULL)
    if (is.null(b)) return(list(ok = FALSE))
  }
  
  # Stack using c (not terra::c)
  s <- tryCatch(c(a, b), error = function(e) NULL)
  if (is.null(s)) return(list(ok = FALSE))
  
  # Sample or take all values
  if (!is.null(sample_n)) {
    df <- tryCatch(
      terra::spatSample(s, size = sample_n, method = "random",
                        na.rm = TRUE, as.df = TRUE, values = TRUE),
      error = function(e) NULL
    )
    if (is.null(df)) return(list(ok = FALSE))
    df <- df[stats::complete.cases(df), , drop = FALSE]
    if (nrow(df) == 0) return(list(ok = FALSE))
    va <- df[[1]]; vb <- df[[2]]
  } else {
    v <- tryCatch(terra::values(s), error = function(e) NULL)
    if (is.null(v)) return(list(ok = FALSE))
    v <- v[stats::complete.cases(v), , drop = FALSE]
    if (nrow(v) == 0) return(list(ok = FALSE))
    va <- v[,1]; vb <- v[,2]
  }
  
  list(ok = TRUE, va = va, vb = vb, n = length(va))
}



# ------------------------------------------------------------
# Pairwise computation across ALL rasters (even partial extents)
# ------------------------------------------------------------
n <- length(rlist)
ras_names <- names(rlist)

# Matrices
D_jac <- matrix(NA_real_, n, n, dimnames = list(ras_names, ras_names))
D_sor <- matrix(NA_real_, n, n, dimnames = list(ras_names, ras_names))
K_mat <- matrix(NA_real_, n, n, dimnames = list(ras_names, ras_names))

N_overlap <- matrix(0, n, n, dimnames = list(ras_names, ras_names))
Frac_overlap_i <- matrix(NA_real_, n, n, dimnames = list(ras_names, ras_names))
Frac_overlap_j <- matrix(NA_real_, n, n, dimnames = list(ras_names, ras_names))

# Precompute valid-cell counts per raster for overlap fractions
# (counts non-NA cells)
valid_counts <- sapply(rlist, function(r) global(!is.na(r), "sum", na.rm = TRUE)[1,1])

message("Computing pairwise metrics for ", n, " rasters...")

# Should all be TRUE and finite:
all_ok <- vapply(rlist, function(r) {
  e <- terra::ext(r)
  all(is.finite(c(terra::xmin(e), terra::xmax(e), terra::ymin(e), terra::ymax(e))))
}, logical(1))

stopifnot(all(all_ok))

# This should work and return a SpatExtent
terra::intersect(terra::ext(rlist[[1]]), terra::ext(rlist[[2]]))

# Initialize diagonals
for (i in seq_len(n)) {
  D_jac[i, i] <- 0
  D_sor[i, i] <- 0
  K_mat[i, i] <- 1
  N_overlap[i, i] <- valid_counts[i]
  Frac_overlap_i[i, i] <- 1
  Frac_overlap_j[i, i] <- 1
}

pairs <- utils::combn(n, 2)
message(sprintf("Computing pairwise metrics for %d rasters (%d pairs)...", n, ncol(pairs)))

for (k in seq_len(ncol(pairs))) {
  i <- pairs[1, k]
  j <- pairs[2, k]
  
  ov <- extract_pair_overlap(rlist[[i]], rlist[[j]], sample_n = SAMPLE_N_PER_PAIR)
  
  if (!isTRUE(ov$ok)) {
    D_jac[i, j] <- NO_OVERLAP_DISTANCE; D_jac[j, i] <- NO_OVERLAP_DISTANCE
    D_sor[i, j] <- NO_OVERLAP_DISTANCE; D_sor[j, i] <- NO_OVERLAP_DISTANCE
    K_mat[i, j]  <- NA_real_;          K_mat[j, i]  <- NA_real_
    N_overlap[i, j] <- 0;              N_overlap[j, i] <- 0
    Frac_overlap_i[i, j] <- 0;         Frac_overlap_i[j, i] <- 0
    Frac_overlap_j[i, j] <- 0;         Frac_overlap_j[j, i] <- 0
    next
  }
  
  m <- pair_metrics_from_vectors(ov$va, ov$vb, jaccard_empty = JACCARD_EMPTY_BEHAVIOR)
  
  D_jac[i, j] <- m$d_jaccard; D_jac[j, i] <- m$d_jaccard
  D_sor[i, j] <- m$d_sorensen; D_sor[j, i] <- m$d_sorensen
  K_mat[i, j]  <- m$kappa;     K_mat[j, i]  <- m$kappa
  
  N_overlap[i, j] <- m$N;      N_overlap[j, i] <- m$N
  
  # Fraction overlaps (explicit symmetry)
  Frac_overlap_i[i, j] <- m$N / valid_counts[i]
  Frac_overlap_i[j, i] <- m$N / valid_counts[j]
  Frac_overlap_j[i, j] <- m$N / valid_counts[j]
  Frac_overlap_j[j, i] <- m$N / valid_counts[i]
}

# check
stopifnot(
  is.matrix(D_jac), nrow(D_jac) == n, ncol(D_jac) == n,
  is.matrix(D_sor), nrow(D_sor) == n, ncol(D_sor) == n,
  is.matrix(K_mat), nrow(K_mat) == n, ncol(K_mat) == n,
  is.matrix(N_overlap), nrow(N_overlap) == n, ncol(N_overlap) == n,
  is.matrix(Frac_overlap_i), nrow(Frac_overlap_i) == n, ncol(Frac_overlap_i) == n,
  is.matrix(Frac_overlap_j), nrow(Frac_overlap_j) == n, ncol(Frac_overlap_j) == n
)

# Save matrices
saveRDS(D_jac, file.path(OUT_DIR, "pairwise_jaccard_distance_matrix.rds"))
saveRDS(D_sor, file.path(OUT_DIR, "pairwise_sorensen_distance_matrix.rds"))
saveRDS(K_mat, file.path(OUT_DIR, "pairwise_kappa_matrix.rds"))

#-------------------------------------------------------------------

# alternate method for pairwise computations :
# (pairwise loop above is very slow)

# Fast pairwise metrics for binary/categorical rasters using terra

# optional parallel
# library(future.apply); plan(multisession, workers = max(1, parallel::detectCores() - 1))

# --- helper for κ, Jaccard, Sørensen from confusion counts ---
.metrics_from_counts <- function(n00, n01, n10, n11, jaccard_empty = "zero") {
  # Presence is "1", absence "0"
  # n11: both 1 (TP), n10: r1=1,r2=0 (FP), n01: r1=0,r2=1 (FN), n00: both 0 (TN)
  N <- n00 + n01 + n10 + n11
  if (N == 0) {
    return(list(
      d_jaccard = if (identical(jaccard_empty, "one")) 1 else 0,  # or as you define
      d_sorensen = 0,
      kappa = NA_real_,
      N = 0
    ))
  }
  # Jaccard on presences only
  denom_j <- (n11 + n10 + n01)
  j <- if (denom_j == 0) {
    if (identical(jaccard_empty, "one")) 1 else 0
  } else n11 / denom_j
  d_jac <- 1 - j
  
  # Sørensen (Dice)
  denom_s <- (2 * n11 + n10 + n01)
  s <- if (denom_s == 0) 1 else (2 * n11) / denom_s
  d_sor <- 1 - s
  
  # Cohen's kappa
  po <- (n11 + n00) / N
  p1_yes <- (n11 + n10) / N
  p1_no  <- (n01 + n00) / N
  p2_yes <- (n11 + n01) / N
  p2_no  <- (n10 + n00) / N
  pe <- p1_yes * p2_yes + p1_no * p2_no
  kappa <- if (isTRUE(all.equal(1, pe))) NA_real_ else (po - pe) / (1 - pe)
  
  list(d_jaccard = d_jac, d_sorensen = d_sor, kappa = kappa, N = N)
}

# --- extract counts via crosstab (exact) or sampling (approx) ---
.fast_pair_counts <- function(r1, r2, sample_n = NULL, max_cells_for_exact = 5e7) {
  # Quick overlap test
  e1 <- terra::ext(r1); e2 <- terra::ext(r2)
  if (!terra::relate(e1, e2, "intersects")) return(list(ok = FALSE))
  
  e <- terra::intersect(e1, e2)
  if (is.null(e)) return(list(ok = FALSE))
  
  # Cheap crop (no values yet)
  r1o <- terra::crop(r1, e, values = FALSE, snap = "near")
  r2o <- terra::crop(r2, e, values = FALSE, snap = "near")
  
  # Align grids (assume binary/categorical so "near" is appropriate)
  if (!terra::compareGeom(r1o, r2o, stopOnError = FALSE)) {
    r2o <- terra::resample(r2o, r1o, method = "near")
  }
  
  s <- c(r1o, r2o)
  names(s) <- c("a", "b")
  
  # Mask to overlapping valid cells
  both_ok <- !is.na(s$a) & !is.na(s$b)
  s <- terra::mask(s, both_ok, maskvalues = FALSE)
  
  n_over <- terra::global(!is.na(s$a), "sum", na.rm = TRUE)[1, 1]
  if (is.na(n_over) || n_over == 0) return(list(ok = FALSE))
  
  # Choose exact vs sample
  do_exact <- is.null(sample_n) && (n_over <= max_cells_for_exact)
  
  if (do_exact) {
    # Exact counts with crosstab (fast for categorical)
    # Ensure values are 0/1 (if already binary, this is a no-op)
    # Comment out if you KNOW inputs are only 0/1.
    s$a <- classify(s$a, rbind(cbind(-Inf, 0, 0), cbind(0, Inf, 1)))
    s$b <- classify(s$b, rbind(cbind(-Inf, 0, 0), cbind(0, Inf, 1)))
    
    tab <- terra::crosstab(s, useNA = FALSE, long = FALSE)
    # Ensure 2x2 layout with zeros if absent
    getv <- function(tbl, r, c) if (r %in% rownames(tbl) && c %in% colnames(tbl)) tbl[r, c] else 0
    n00 <- getv(tab, "0", "0")
    n01 <- getv(tab, "0", "1")
    n10 <- getv(tab, "1", "0")
    n11 <- getv(tab, "1", "1")
    return(list(ok = TRUE, n00 = n00, n01 = n01, n10 = n10, n11 = n11))
  } else {
    # Sample-based counts
    df <- terra::spatSample(s, size = if (is.null(sample_n)) 1e6 else sample_n,
                            method = "random", na.rm = TRUE, values = TRUE)
    # Coerce to 0/1 if needed
    df$a <- as.integer(df$a > 0)
    df$b <- as.integer(df$b > 0)
    tab <- table(df$a, df$b)
    n00 <- ifelse("0" %in% rownames(tab) && "0" %in% colnames(tab), tab["0","0"], 0)
    n01 <- ifelse("0" %in% rownames(tab) && "1" %in% colnames(tab), tab["0","1"], 0)
    n10 <- ifelse("1" %in% rownames(tab) && "0" %in% colnames(tab), tab["1","0"], 0)
    n11 <- ifelse("1" %in% rownames(tab) && "1" %in% colnames(tab), tab["1","1"], 0)
    return(list(ok = TRUE, n00 = n00, n01 = n01, n10 = n10, n11 = n11))
  }
}

# --- main: compute pairwise matrices (fast + parallel) ---
compute_pairwise_metrics_fast <- function(rlist, sample_n = NULL,
                                          jaccard_empty = "zero",
                                          parallel = TRUE) {
  n <- length(rlist)
  ras_names <- names(rlist)
  D_jac <- matrix(NA_real_, n, n, dimnames = list(ras_names, ras_names))
  D_sor <- matrix(NA_real_, n, n, dimnames = list(ras_names, ras_names))
  K_mat <- matrix(NA_real_, n, n, dimnames = list(ras_names, ras_names))
  N_overlap <- matrix(0, n, n, dimnames = list(ras_names, ras_names))
  Frac_overlap_i <- matrix(NA_real_, n, n, dimnames = list(ras_names, ras_names))
  Frac_overlap_j <- matrix(NA_real_, n, n, dimnames = list(ras_names, ras_names))
  
  # valid non-NA counts per raster (same as your code)
  valid_counts <- sapply(rlist, function(r) terra::global(!is.na(r), "sum", na.rm = TRUE)[1,1])
  
  # diagonals
  diag(D_jac) <- 0; diag(D_sor) <- 0; diag(K_mat) <- 1
  diag(N_overlap) <- valid_counts
  diag(Frac_overlap_i) <- 1; diag(Frac_overlap_j) <- 1
  
  pairs <- utils::combn(n, 2)
  idx <- seq_len(ncol(pairs))
  
  .one <- function(k) {
    i <- pairs[1, k]; j <- pairs[2, k]
    cnt <- .fast_pair_counts(rlist[[i]], rlist[[j]], sample_n = sample_n)
    if (!isTRUE(cnt$ok)) {
      list(i=i, j=j,
           d_j=NO_OVERLAP_DISTANCE, d_s=NO_OVERLAP_DISTANCE, kap=NA_real_,
           N=0, fi=0, fj=0, ok=FALSE)
    } else {
      m <- .metrics_from_counts(cnt$n00, cnt$n01, cnt$n10, cnt$n11,
                                jaccard_empty = jaccard_empty)
      N <- m$N
      fi <- if (valid_counts[i] > 0) N / valid_counts[i] else 0
      fj <- if (valid_counts[j] > 0) N / valid_counts[j] else 0
      list(i=i, j=j, d_j=m$d_jaccard, d_s=m$d_sorensen, kap=m$kappa,
           N=N, fi=fi, fj=fj, ok=TRUE)
    }
  }
  
  res_list <- if (parallel) {
    # future.apply path (uncomment plan() before calling)
    # future_lapply(idx, .one)
    lapply(idx, .one)  # fall back to serial if you haven't set up plan()
  } else {
    lapply(idx, .one)
  }
  
  # Fill upper triangle then mirror
  for (res in res_list) {
    i <- res$i; j <- res$j
    D_jac[i, j] <- res$d_j; D_sor[i, j] <- res$d_s; K_mat[i, j] <- res$kap
    N_overlap[i, j] <- res$N
    Frac_overlap_i[i, j] <- res$fi
    Frac_overlap_j[i, j] <- res$fj
  }
  # Mirror
  D_jac[lower.tri(D_jac)] <- t(D_jac)[lower.tri(D_jac)]
  D_sor[lower.tri(D_sor)] <- t(D_sor)[lower.tri(D_sor)]
  K_mat[lower.tri(K_mat)] <- t(K_mat)[lower.tri(K_mat)]
  N_overlap[lower.tri(N_overlap)] <- t(N_overlap)[lower.tri(N_overlap)]
  Frac_overlap_i[lower.tri(Frac_overlap_i)] <- t(Frac_overlap_i)[lower.tri(Frac_overlap_i)]
  Frac_overlap_j[lower.tri(Frac_overlap_j)] <- t(Frac_overlap_j)[lower.tri(Frac_overlap_j)]
  
  list(
    D_jac = D_jac,
    D_sor = D_sor,
    K_mat = K_mat,
    N_overlap = N_overlap,
    Frac_overlap_i = Frac_overlap_i,
    Frac_overlap_j = Frac_overlap_j
  )
}


# ------------------------------------------------------------
# PCoA, clustering, and heatmaps
# ------------------------------------------------------------
# PCoA on Jaccard distances
pcoa_res <- cmdscale(D_jac, eig = TRUE, k = 2)

pcoa_df <- data.frame(
  x = pcoa_res$points[,1],
  y = pcoa_res$points[,2],
  raster = rownames(D_jac)
)

pcoa_df$label_prefix <- sub("_.*", "", pcoa_df$raster)

ggplot(pcoa_df, aes(x, y, label = label_prefix)) +
  geom_point(size = 3) +
  geom_label_repel(
    size = 3.2,
    fill = scales::alpha("white", 0.8),
    color = "black",
    label.size = 0.2,
    box.padding = 0.25,
    point.padding = 0.2,
    max.overlaps = Inf,
    seed = 42
  ) +
  theme_minimal() +
  labs(title = "PCoA based on Jaccard Distance")

ggsave(file.path(OUT_DIR, "pcoa_jaccard.png"), width = 7, height = 6)

# Variance percentages even for non‑Euclidean distances - robust PCoA with Cailliez correction
pcoa_ape <- ape::pcoa(as.dist(D_jac), correction = "cailliez")
pcoa_df <- data.frame(
  x = pcoa_ape$vectors[, 1],
  y = pcoa_ape$vectors[, 2],
  raster = rownames(D_jac),
  stringsAsFactors = FALSE
) |>
  mutate(
    label_prefix = ifelse(grepl("_", raster), sub("_.*", "", raster), raster),
    raster_base  = tools::file_path_sans_ext(raster)
  )
ax1_pct <- round(100 * pcoa_ape$values$Relative_eig[1], 1)
ax2_pct <- round(100 * pcoa_ape$values$Relative_eig[2], 1)

# PCoA plot with axis labels and variance explained
p <- ggplot(pcoa_df, aes(x, y, label = label_prefix)) + 
  geom_point(size = 3) +
  geom_label_repel(
    size = 3.2,
    fill = scales::alpha("white", 0.8),
    color = "black",
    label.size = 0.2,
    box.padding = 0.25,
    point.padding = 0.2,
    max.overlaps = Inf,
    seed = 42) +
  theme_minimal(base_size = 12) +
  labs(
    title = "PCoA based on Jaccard Distance",
    x = sprintf("Axis 1 (%.1f%%)", ax1_pct),
    y = sprintf("Axis 2 (%.1f%%)", ax2_pct))

print(p)

ggsave(file.path(OUT_DIR, "pcoa_jaccard_v2.png"), 
       plot = p, width = 7, height = 6, dpi = 300)

# Try to join by both with/without extension
meta_join <- meta %>%
  mutate(
    raster_base = tools::file_path_sans_ext(raster),
    raster_full = raster
  )

pcoa_df2 <- pcoa_df %>%
  # First try exact match on full name
  left_join(meta_join, by = c("raster" = "raster_full")) %>%
  # If some models are still NA, try matching by basename (sans extension)
  mutate(model = ifelse(is.na(model),
                        meta_join$model[match(raster_base, meta_join$raster_base)],
                        model)) %>%
  mutate(model = factor(model))

set.seed(42)  # reproducible label placement

# Ellipses need >= 3 points per group; drop groups with < 3
valid_groups <- names(which(table(pcoa_df2$model) >= 3))
pcoa_df2$ellipse_ok <- pcoa_df2$model %in% valid_groups

p <- ggplot(pcoa_df2, aes(x, y, color = model, label = label_prefix)) +  
  geom_point(size = 3) +
  # 95% ellipses for groups with enough points
  #stat_ellipse(
  #data = subset(pcoa_df2, ellipse_ok),
  #aes(group = model),
  #level = 0.95,
  #type = "norm",      # or "t" for small samples; "norm" is common
  #linewidth = 0.6,
  #alpha = 0.5
  #) +
  geom_label_repel(
    size = 3.2,
    fill = scales::alpha("white", 0.8),
    color = "black",
    label.size = 0.2,
    box.padding = 0.25,
    point.padding = 0.2,
    max.overlaps = Inf,
    seed = 42) +
  theme_minimal(base_size = 12) +
  labs(
    title = "PCoA based on Jaccard Distance",
    x = sprintf("Axis 1 (%.1f%%)", ax1_pct),
    y = sprintf("Axis 2 (%.1f%%)", ax2_pct),
    color = "Model"
  )

print(p)

ggsave(file.path(OUT_DIR, "pcoa_jaccard_v2.png"), #with_ellipses.png"),
       plot = p, width = 7, height = 6, dpi = 300)

# Ensure metadata aligns 1:1 with the distance matrix rows
df_perm <- data.frame(
  raster = rownames(D_jac),
  raster_base = tools::file_path_sans_ext(rownames(D_jac)),
  stringsAsFactors = FALSE
) %>%
  left_join(meta_join, by = c("raster" = "raster_full")) %>%
  mutate(model = ifelse(is.na(model),
                        meta_join$model[match(raster_base, meta_join$raster_base)],
                        model)) %>%
  mutate(model = factor(model))

# Sanity check
stopifnot(nrow(df_perm) == nrow(D_jac))

# Run PERMANOVA (999 permutations by default)
set.seed(42)
permanova <- vegan::adonis2(as.dist(D_jac) ~ model, data = df_perm, permutations = 999)

print(permanova)
# Prints ANOVA-like table with F, R2, and p-value for model.
# R² = proportion of distance variation explained by model
# Pr(>F) (permutation p-value) tests whether group centroids differ in multivariate space
# Significant result means models differ in their Jaccard‑space composition

### Hierarchical clustering
#hc <- hclust(as.dist(D_jac), method = "average")
#png(file.path(OUT_DIR, "clustering_jaccard.png"), width = 800, height = 600)
#plot(hc, main = "Hierarchical Clustering (Jaccard)")
#dev.off()

# Build prefixes from the distance matrix rownames
label_prefix <- ifelse(
  grepl("_", rownames(D_jac)),
  sub("_.*", "", rownames(D_jac)),
  rownames(D_jac)
)

# Hierarchical clustering 
hc <- hclust(as.dist(D_jac), method = "average")

# Replace dendrogram labels with prefixes (preserves clustering)
hc$labels <- label_prefix[match(hc$labels, rownames(D_jac))]

png(file.path(OUT_DIR, "clustering_jaccard.png"), width = 800, height = 600)
plot(hc, main = "Hierarchical Clustering (Jaccard)")
dev.off()

### Heatmap of the Jaccard distances
#pheatmap(
#  D_jac,
#  clustering_method = "average",
#  display_numbers = TRUE,
#  filename = file.path(OUT_DIR, "heatmap_jaccard.png")
#)

# Compute label prefixes from current matrix dimnames
label_prefix_row <- ifelse(
  grepl("_", rownames(D_jac)),
  sub("_.*", "", rownames(D_jac)),
  rownames(D_jac)
)
label_prefix_col <- ifelse(
  grepl("_", colnames(D_jac)),
  sub("_.*", "", colnames(D_jac)),
  colnames(D_jac)
)

# Optional: ensure uniqueness in case of duplicate prefixes
# (helps readability; pheatmap can show repeated labels but unique is nicer)
#label_prefix_row <- make.unique(label_prefix_row, sep = "_")
#label_prefix_col <- make.unique(label_prefix_col, sep = "_")

pheatmap(
  D_jac,
  clustering_method = "average",
  display_numbers = TRUE,               # keep to print values
  number_format = "%.2f",               # nicer formatting for distances
  labels_row = label_prefix_row,        # <-- short labels on rows
  labels_col = label_prefix_col,        # <-- short labels on cols
  filename = file.path(OUT_DIR, "heatmap_jaccard.png")
)

# checks
compareGeom(rlist[[1]], rlist[[2]])
ov <- extract_pair_overlap(rlist[[1]], rlist[[2]])
str(ov)
summary(as.vector(N_overlap))
lobal(rlist[[1]], "sum", na.rm = TRUE)
global(rlist[[2]], "sum", na.rm = TRUE)
range(D_jac)
do.call(rbind, lapply(rlist, function(r) as.vector(ext(r))))


a <- terra::crop(rlist[[1]], terra::ext(rlist[[1]]))
b <- terra::crop(rlist[[2]], terra::ext(rlist[[2]]))
c(a, b)


# ============================================================
# ECODISTRICT-LEVEL SUMMARY OF AGREEMENT/DISAGREEMENT
# ============================================================

message("Building ecodistrict-level agreement/disagreement summary...")

# ---- Read ecodistrict polygons as SpatVector
eco <- if (is.null(ECODISTRICT_LAYER)) {
  vect(ECODISTRICT_SHP)
} else {
  vect(ECODISTRICT_SHP, layer = ECODISTRICT_LAYER)
}

if (!(ECODISTRICT_ID %in% names(eco))) {
  stop("ECODISTRICT_ID field '", ECODISTRICT_ID, "' not found in ecodistrict data. Fields: ",
       paste(names(eco), collapse = ", "))
}

# ---- Choose target CRS + resolution from first raster
r0 <- rlist[[1]]
target_crs <- crs(r0)

# Ensure ecodistricts (SpatVector) use the same CRS as the target raster (SpatRaster)
target_crs <- terra::crs(r0, proj = TRUE)  # get CRS string from raster

if (!terra::same.crs(eco, r0)) {
  eco <- terra::project(eco, target_crs)
}

# ---- Reproject rasters to target CRS if needed (keeps categorical, nearest-neighbor)
rlist_proj <- lapply(rlist, function(r) {
  if (!same.crs(r, target_crs)) {
    project(r, target_crs, method = "near")
  } else {
    r
  }
})

# ---- Union extent across rasters (in target CRS)
ext_union <- ext(rlist_proj[[1]])
if (length(rlist_proj) > 1) {
  for (k in 2:length(rlist_proj)) ext_union <- union(ext_union, ext(rlist_proj[[k]]))
}


# Convert numeric extent to SpatExtent
ext_union_spat <- terra::ext(ext_union)

# Create template raster covering union extent
template <- terra::rast(
  extent    = ext_union_spat,
  resolution = terra::res(r0),
  crs        = target_crs
)

template


# ---- Resample each raster to template grid (nearest neighbor); keep NA where no data
# This can take time/memory if rasters are large.
r_aligned <- lapply(rlist_proj, function(r) resample(r, template, method = "near"))
r_stack <- rast(r_aligned)
names(r_stack) <- names(rlist)

# ---- Pixel-level summaries:
# n_valid = number of non-NA rasters per cell
# n1      = number of rasters == 1 per cell (ignoring NAs)
n_valid <- app(r_stack, fun = function(...) {
  x <- c(...)
  sum(!is.na(x))
})

n1 <- app(r_stack, fun = function(...) {
  x <- c(...)
  sum(x == 1, na.rm = TRUE)
})

# Avoid division by 0
consensus_p <- n1 / n_valid  # fraction of "presence" among available rasters

# Agreement/disagreement definitions (only where enough rasters have data)
enough <- n_valid >= MIN_RASTERS_FOR_AGREEMENT

agree_all <- enough & (n1 == 0 | n1 == n_valid)     # full agreement: all 0 or all 1
agree_presence <- enough & (n1 == n_valid)          # full agreement on presence (all 1)
agree_absence  <- enough & (n1 == 0)                # full agreement on absence (all 0)
disagree <- enough & !(agree_all)                   # mixed 0/1 among available rasters

# Optional: uncertainty score (0 = full agreement, max at p=0.5)
# This is a handy continuous "disagreement intensity" metric.
uncertainty <- enough * (4 * consensus_p * (1 - consensus_p))

# Make sure these are numeric (zonal stats need numeric rasters)
agree_all_r      <- as.int(agree_all)
agree_pres_r     <- as.int(agree_presence)
agree_abs_r      <- as.int(agree_absence)
disagree_r       <- as.int(disagree)
enough_r         <- as.int(enough)

# ---- Rasterize ecodistrict IDs to the template grid
# (each cell gets the ecodistrict id; NA outside polygons)
eco_id_r <- rasterize(eco, template, field = ECODISTRICT_ID)

# ---- Helper to compute zonal sums & means safely
zonal_sum <- function(x) zonal(x, eco_id_r, fun = "sum", na.rm = TRUE)
zonal_mean <- function(x) zonal(x, eco_id_r, fun = "mean", na.rm = TRUE)

# ---- Zonal aggregation:
# Count cells within ecodistricts with enough raster coverage
z_enough   <- zonal_sum(enough_r)
z_agree    <- zonal_sum(agree_all_r)
z_agree1   <- zonal_sum(agree_pres_r)
z_agree0   <- zonal_sum(agree_abs_r)
z_disagree <- zonal_sum(disagree_r)

# Mean consensus (presence fraction) and mean uncertainty, computed over valid-enough cells
# We mask to only 'enough' cells to avoid mixing in NA/low-coverage areas.
consensus_masked <- mask(consensus_p, enough, maskvalues = 0)
uncert_masked    <- mask(uncertainty, enough, maskvalues = 0)

z_cons_mean <- zonal_mean(consensus_masked)
z_unc_mean  <- zonal_mean(uncert_masked)

# ---- Combine into one summary table
# Each zonal() returns a data.frame with columns: zone, value
summary_df <- dplyr::tibble(
  !!ECODISTRICT_ID := z_enough[,1],
  cells_enough = z_enough[,2],
  cells_agree_all = z_agree[,2],
  cells_agree_presence = z_agree1[,2],
  cells_agree_absence  = z_agree0[,2],
  cells_disagree = z_disagree[,2]
) %>%
  left_join(dplyr::tibble(!!ECODISTRICT_ID := z_cons_mean[,1], consensus_presence_mean = z_cons_mean[,2]),
            by = ECODISTRICT_ID) %>%
  left_join(dplyr::tibble(!!ECODISTRICT_ID := z_unc_mean[,1], uncertainty_mean = z_unc_mean[,2]),
            by = ECODISTRICT_ID) %>%
  mutate(
    frac_agree_all = ifelse(cells_enough > 0, cells_agree_all / cells_enough, NA_real_),
    frac_agree_presence = ifelse(cells_enough > 0, cells_agree_presence / cells_enough, NA_real_),
    frac_agree_absence  = ifelse(cells_enough > 0, cells_agree_absence / cells_enough, NA_real_),
    frac_disagree = ifelse(cells_enough > 0, cells_disagree / cells_enough, NA_real_)
  )

# ---- Write CSV
readr::write_csv(summary_df, OUT_ECODISTRICT_CSV)
message("Wrote ecodistrict summary CSV: ", OUT_ECODISTRICT_CSV)

# ---- Optionally write a GPKG with attributes joined back to polygons
if (WRITE_ECODISTRICT_GPKG) {
  eco_out <- eco
  eco_out <- merge(eco_out, summary_df, by.x = ECODISTRICT_ID, by.y = ECODISTRICT_ID, all.x = TRUE)
  if (file.exists(OUT_ECODISTRICT_GPKG)) file.remove(OUT_ECODISTRICT_GPKG)
  writeVector(eco_out, OUT_ECODISTRICT_GPKG, filetype = "GPKG")
  message("Wrote ecodistrict summary GPKG: ", OUT_ECODISTRICT_GPKG)
}

# ============================================================
# ECODISTRICT-LEVEL PAIRWISE AGREEMENT/DISAGREEMENT 
# ============================================================

message("Computing ecodistrict-level pairwise agreement/disagreement...")

# ---- Read ecodistrict polygons
eco <- if (is.null(ECODISTRICT_LAYER)) {
  terra::vect(ECODISTRICT_SHP)
} else {
  terra::vect(ECODISTRICT_SHP, layer = ECODISTRICT_LAYER)
}

if (!(ECODISTRICT_ID %in% names(eco))) {
  stop("ECODISTRICT_ID field '", ECODISTRICT_ID, "' not found in ecodistrict data. Fields: ",
       paste(names(eco), collapse = ", "))
}

# ---- Choose target CRS + resolution from first raster
r0 <- rlist[[1]]
target_crs <- terra::crs(r0)

# Project ecodistricts to raster CRS if needed
if (!terra::same.crs(eco, r0)) {
  eco <- terra::project(eco, target_crs)
}

# ---- Project rasters to target CRS if needed (nearest for categorical 0/1)
rlist_proj <- lapply(rlist, function(r) {
  if (!terra::same.crs(r, target_crs)) terra::project(r, target_crs, method = "near") else r
})

# ---- Build union template so all rasters align to identical grid
ext_union <- terra::ext(rlist_proj[[1]])
if (length(rlist_proj) > 1) {
  for (k in 2:length(rlist_proj)) ext_union <- terra::union(ext_union, terra::ext(rlist_proj[[k]]))
}

template <- terra::rast(ext_union, resolution = terra::res(r0), crs = target_crs)

# ---- Resample all rasters to template grid (nearest neighbor)
r_aligned <- lapply(rlist_proj, function(r) terra::resample(r, template, method = "near"))
r_stack <- terra::rast(r_aligned)
names(r_stack) <- names(rlist)

ras_names <- terra::names(r_stack)
n_ras <- length(ras_names)

# ------------------------------------------------------------
# Helper: metrics for one ecodistrict (data.frame of raster cols)
# ------------------------------------------------------------
eco_pairwise_metrics <- function(df_zone, zone_id,
                                 sample_n = NULL,
                                 min_overlap = 1,
                                 jaccard_empty_behavior = "zero") {
  # df_zone: columns are rasters in ras_names (numeric 0/1/NA)
  out <- vector("list", length = (n_ras * (n_ras - 1)) / 2)
  idx_out <- 1
  
  # Precompute valid counts per raster within the ecodistrict
  valid_zone <- sapply(ras_names, function(nm) sum(!is.na(df_zone[[nm]])))
  
  for (i in seq_len(n_ras - 1)) {
    a0 <- df_zone[[ras_names[i]]]
    for (j in (i + 1):n_ras) {
      b0 <- df_zone[[ras_names[j]]]
      
      ok <- stats::complete.cases(a0, b0)
      N <- sum(ok)
      
      if (N < min_overlap) {
        out[[idx_out]] <- tibble::tibble(
          eco_id = zone_id,
          raster_i = ras_names[i],
          raster_j = ras_names[j],
          N = N,
          TP = NA_integer_, FP = NA_integer_, FN = NA_integer_, TN = NA_integer_,
          d_jaccard = if (N == 0) NO_OVERLAP_DISTANCE else NA_real_,
          d_sorensen = if (N == 0) NO_OVERLAP_DISTANCE else NA_real_,
          kappa = NA_real_,
          frac_overlap_i_zone = ifelse(valid_zone[i] > 0, N / valid_zone[i], NA_real_),
          frac_overlap_j_zone = ifelse(valid_zone[j] > 0, N / valid_zone[j], NA_real_)
        )
        idx_out <- idx_out + 1
        next
      }
      
      # Optional sampling within this zone & pair
      if (!is.null(sample_n) && N > sample_n) {
        take <- sample(which(ok), size = sample_n, replace = FALSE)
        a <- a0[take]; b <- b0[take]
      } else {
        a <- a0[ok]; b <- b0[ok]
      }
      
      m <- pair_metrics_from_vectors(a, b, jaccard_empty = jaccard_empty_behavior)
      
      out[[idx_out]] <- tibble::tibble(
        eco_id = zone_id,
        raster_i = ras_names[i],
        raster_j = ras_names[j],
        N = m$N,
        TP = m$TP, FP = m$FP, FN = m$FN, TN = m$TN,
        d_jaccard = m$d_jaccard,
        d_sorensen = m$d_sorensen,
        kappa = m$kappa,
        frac_overlap_i_zone = ifelse(valid_zone[i] > 0, m$N / valid_zone[i], NA_real_),
        frac_overlap_j_zone = ifelse(valid_zone[j] > 0, m$N / valid_zone[j], NA_real_)
      )
      idx_out <- idx_out + 1
    }
  }
  
  dplyr::bind_rows(out)
}

# ------------------------------------------------------------
# Extract raster values within ecodistricts
# NOTE: This can be large. If memory is tight, see the "low-memory" option below.
# ------------------------------------------------------------
# extract() returns a data.frame with polygon ID column (ID) and one column per raster layer
ex <- terra::extract(r_stack, eco, ID = TRUE, cells = FALSE, xy = FALSE)

# Rename the extract ID column to something stable
names(ex)[1] <- ".__poly_id__"

# Map extract polygon IDs back to ecodistrict IDs (attribute field)
# terra::extract ID is row index in 'eco' by default when ID=TRUE.
poly_to_eco <- tibble::tibble(
  .__poly_id__ = seq_len(nrow(eco)),
  eco_id = eco[[ECODISTRICT_ID]][,1]
)

ex <- dplyr::left_join(ex, poly_to_eco, by = ".__poly_id__")

# ------------------------------------------------------------
# Compute pairwise metrics per ecodistrict
# ------------------------------------------------------------
# Split by eco_id
ex_split <- split(ex[, ras_names, drop = FALSE], ex$eco_id)

message("Ecodistricts found: ", length(ex_split), " (computing pairwise metrics within each...)")

eco_long <- lapply(names(ex_split), function(zid) {
  dfz <- ex_split[[zid]]
  eco_pairwise_metrics(
    df_zone = dfz,
    zone_id = zid,
    sample_n = SAMPLE_N_PER_PAIR_PER_ECO,
    min_overlap = MIN_OVERLAP_CELLS_PER_PAIR,
    jaccard_empty_behavior = JACCARD_EMPTY_BEHAVIOR
  )
}) |> dplyr::bind_rows()

# Add convenient similarity versions (agreement-ish)
eco_long <- eco_long %>%
  dplyr::mutate(
    jaccard_similarity  = 1 - d_jaccard,
    sorensen_similarity = 1 - d_sorensen,
    kappa_disagreement  = 1 - kappa
  )

# Write long output
readr::write_csv(eco_long, OUT_ECO_PAIRWISE_LONG)
message("Wrote ecodistrict pairwise long table: ", OUT_ECO_PAIRWISE_LONG)

# ------------------------------------------------------------
# Ecodistrict summary across pairs (mean + N-weighted mean)
# ------------------------------------------------------------
eco_summary <- eco_long %>%
  dplyr::filter(!is.na(d_jaccard), !is.na(d_sorensen)) %>%
  dplyr::group_by(eco_id) %>%
  dplyr::summarise(
    n_pairs = dplyr::n(),
    # Simple means (unweighted)
    mean_d_jaccard = mean(d_jaccard, na.rm = TRUE),
    mean_d_sorensen = mean(d_sorensen, na.rm = TRUE),
    mean_kappa = mean(kappa, na.rm = TRUE),
    mean_jaccard_similarity = mean(jaccard_similarity, na.rm = TRUE),
    mean_sorensen_similarity = mean(sorensen_similarity, na.rm = TRUE),
    
    # Weighted by overlap N (more overlapping cells => more weight)
    wmean_d_jaccard = stats::weighted.mean(d_jaccard, w = pmax(N, 1), na.rm = TRUE),
    wmean_d_sorensen = stats::weighted.mean(d_sorensen, w = pmax(N, 1), na.rm = TRUE),
    wmean_kappa = stats::weighted.mean(kappa, w = pmax(N, 1), na.rm = TRUE),
    
    # Coverage / overlap diagnostics
    mean_N_overlap = mean(N, na.rm = TRUE),
    median_N_overlap = stats::median(N, na.rm = TRUE),
    .groups = "drop"
  )

readr::write_csv(eco_summary, OUT_ECO_PAIRWISE_SUM)
message("Wrote ecodistrict pairwise summary: ", OUT_ECO_PAIRWISE_SUM)

# ==============================================
# Map ecodistrict-level results
# ==============================================

# ----------------
# RESULTS AS INPUTS
# ----------------
# Point towards analysis results
OUT_ECODISTRICT_GPKG <- file.path(OUT_DIR, "ecodistrict_agreement_disagreement_summary.gpkg")
OUT_ECODISTRICT_CSV  <- file.path(OUT_DIR, "ecodistrict_agreement_disagreement_summary.csv")

# Pairwise summary CSV (mean/wmean distances per ecodistrict)
OUT_ECO_PAIRWISE_SUM <- file.path(OUT_DIR, "ecodistrict_pairwise_summary.csv")

# ----------------
# READ DATA
# ----------------
eco_sf <- NULL

if (file.exists(OUT_ECODISTRICT_GPKG)) {
  message("Reading GPKG with attributes: ", OUT_ECODISTRICT_GPKG)
  eco_sf <- sf::st_read(OUT_ECODISTRICT_GPKG, quiet = TRUE)
} else {
  message("GPKG not found; reading shapefile + joining CSV…")
  stopifnot(file.exists(ECODISTRICT_SHP), file.exists(OUT_ECODISTRICT_CSV))
  eco_sf <- sf::st_read(ECODISTRICT_SHP, quiet = TRUE)
  stopifnot(ECODISTRICT_ID %in% names(eco_sf))
  eco_tab <- readr::read_csv(OUT_ECODISTRICT_CSV, show_col_types = FALSE)
  # Join by ID field; the CSV column is also named ECODISTRICT_ID
  eco_sf <- eco_sf %>%
    left_join(eco_tab, by = setNames(ECODISTRICT_ID, ECODISTRICT_ID))
}

# If geometry is missing or invalid, try to fix:
if (is.null(sf::st_geometry(eco_sf))) {
  stop("No geometry found after read. Check inputs.")
}
eco_sf <- sf::st_make_valid(eco_sf)

# ----------------
# METRICS & LABELS
# ----------------
# Core metrics 
metric_labels <- c(
  frac_agree_all        = "Full Agreement (All 0 or All 1)",
  frac_agree_presence   = "Agreement on Presence (All 1)",
  frac_agree_absence    = "Agreement on Absence (All 0)",
  frac_disagree         = "Disagreement (Mixed 0/1)",
  consensus_presence_mean = "Mean Consensus (Presence Fraction)",
  uncertainty_mean      = "Uncertainty (0=Agree, Max at 0.5)"
)

# Metrics you want to map (these must exist in eco_sf)
metrics <- c(
  "frac_agree_all",
  "frac_disagree",
  "consensus_presence_mean",
  "uncertainty_mean"
)

# Ensure these exist in the attribute table; warn otherwise
missing_metrics <- setdiff(names(metric_labels), names(eco_sf))
if (length(missing_metrics) > 0) {
  warning("These expected metrics were not found in the ecodistrict data: ",
          paste(missing_metrics, collapse = ", "))
}

# Keep only ecodistricts that:
#  - have at least one valid raster cell (cells_enough > 0)
#  - AND have at least one non-NA metric value
eco_sf_data <- eco_sf %>%
  filter(cells_enough > 0) %>% 
  filter_at(vars(all_of(metrics)), any_vars(!is.na(.)))

# Confirm subset
message("Number of ecodistricts with data: ", nrow(eco_sf_data))

# ----------------
# COMMON THEMES
# ----------------
theme_map <- theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text  = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right"
  )

fill_scale_01 <- scale_fill_viridis(
  option = "C", direction = 1, na.value = "grey90",
  limits = c(0, 1), oob = scales::squish
)

# ----------------
# HELPER: single metric map
# ----------------
plot_metric_map <- function(eco_sf, metric, title = NULL, limits_01 = TRUE) {
  if (!metric %in% names(eco_sf)) {
    stop("Metric '", metric, "' not found in ecodistrict attributes.")
  }
  p <- ggplot(eco_sf_data) +
    geom_sf(aes(fill = .data[[metric]]), color = "grey30", size = 0.1) +
    (if (limits_01) fill_scale_01 else scale_fill_viridis(option = "C", na.value = "grey90")) +
    theme_map +
    labs(
      title = title %||% metric,
      fill  = if (metric %in% names(metric_labels)) metric_labels[[metric]] else metric
    ) #+
  #annotation_scale(location = "bl", width_hint = 0.2, text_cex = 0.7) +
  #annotation_north_arrow(location = "bl", which_north = "true",
  #                       pad_x = unit(0.5, "cm"), pad_y = unit(0.8, "cm"),
  #                       style = north_arrow_fancy_orienteering)
  
  return(p)
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

# ----------------
# EXAMPLE: single-metric maps
# ----------------
# 1) Full agreement
p1 <- plot_metric_map(eco_sf_data, "frac_agree_all", "Proportion of Ecodistrict with Full Agreement (All 0 or All 1)")
ggsave(file.path(OUT_DIR, "map_frac_agree_all.png"), p1, width = 8, height = 6, dpi = 300)

# 2) Disagreement
p2 <- plot_metric_map(eco_sf_data, "frac_disagree", "Proportion of Ecodistrict with Disagreement (Mixed 0/1)")
ggsave(file.path(OUT_DIR, "map_frac_disagree.png"), p2, width = 8, height = 6, dpi = 300)

# 3) Mean consensus (presence fraction)
p3 <- plot_metric_map(eco_sf_data, "consensus_presence_mean", "Ecodistricts — Mean Consensus (Presence Fraction)")
ggsave(file.path(OUT_DIR, "map_consensus_presence_mean.png"), p3, width = 8, height = 6, dpi = 300)

# 4) Uncertainty
p4 <- plot_metric_map(eco_sf_data, "uncertainty_mean", "Ecodistricts — Uncertainty")
ggsave(file.path(OUT_DIR, "map_uncertainty_mean.png"), p4, width = 8, height = 6, dpi = 300)

# ----------------
# FACETED SMALL-MULTIPLES
# ----------------
# Pivot selected metrics long and facet
metrics_to_facet <- intersect(names(metric_labels), names(eco_sf_data))  # only those present

eco_long <- eco_sf_data %>%
  sf::st_drop_geometry() %>%
  select(!!ECODISTRICT_ID, all_of(metrics_to_facet)) %>%
  pivot_longer(cols = -all_of(ECODISTRICT_ID), names_to = "metric", values_to = "value") %>%
  mutate(metric_lab = recode(metric, !!!metric_labels))

eco_sf_long <- eco_sf_data %>%
  dplyr::select(!!rlang::sym(ECODISTRICT_ID)) %>%  # geometry is preserved implicitly
  dplyr::inner_join(eco_long, by = ECODISTRICT_ID)

p_facets <- ggplot(eco_sf_long) +
  geom_sf(aes(fill = value), color = "grey30", size = 0.05) +
  fill_scale_01 +
  facet_wrap(~ metric_lab, ncol = 3) +
  theme_map +
  labs(title = "Ecodistrict-level Metrics", fill = "Value")

ggsave(file.path(OUT_DIR, "maps_ecodistrict_faceted.png"),
       p_facets, width = 12, height = 8, dpi = 300)

# ----------------
# OPTIONAL: Label ecodistrict IDs (centroid labels)
# ----------------
# (Uncomment if desired; can get cluttered for many polygons)
# eco_cent <- st_point_on_surface(eco_sf_data) %>% 
#   mutate(lbl = .data[[ECODISTRICT_ID]])
# p_labels <- p2 + 
#   geom_sf_text(data = eco_cent, aes(label = lbl), size = 2, color = "black")
# ggsave(file.path(OUT_DIR, "map_frac_disagree_with_labels.png"), p_labels, width = 8, height = 6, dpi = 300)

# ----------------
# OPTIONAL: Pairwise summary maps (mean/wmean distances, kappa)
# ----------------
if (file.exists(OUT_ECO_PAIRWISE_SUM)) {
  eco_pair <- read_csv(OUT_ECO_PAIRWISE_SUM, show_col_types = FALSE)
  # Expect columns like: eco_id, mean_d_jaccard, wmean_d_jaccard, mean_kappa, etc.
  id_col <- if (ECODISTRICT_ID %in% names(eco_pair)) ECODISTRICT_ID else "eco_id"
  if (!(id_col %in% names(eco_pair))) {
    warning("Could not find ecodistrict ID column in pairwise summary; skipping pairwise maps.")
  } else {
    eco_pair_sf <- eco_sf_data %>%
      left_join(eco_pair, by = setNames(id_col, ECODISTRICT_ID))
    
    # Map weighted mean Jaccard distance (0..1, higher = more dissimilar)
    if ("wmean_d_jaccard" %in% names(eco_pair_sf)) {
      p_jac <- ggplot(eco_pair_sf) +
        geom_sf(aes(fill = wmean_d_jaccard), color = "grey30", size = 0.1) +
        scale_fill_viridis(option = "C", na.value = "grey90", limits = c(0, 1)) +
        theme_map +
        labs(title = "Ecodistricts — Weighted Mean Jaccard Distance", fill = "WMean Jaccard")
      
      ggsave(file.path(OUT_DIR, "map_wmean_jaccard.png"), p_jac, width = 8, height = 6, dpi = 300)
    }
    
    # Map mean kappa (-1..1, higher = more agreement)
    if ("mean_kappa" %in% names(eco_pair_sf)) {
      p_kap <- ggplot(eco_pair_sf) +
        geom_sf(aes(fill = mean_kappa), color = "grey30", size = 0.1) +
        scale_fill_viridis(option = "C", na.value = "grey90") +
        theme_map +
        labs(title = "Ecodistricts — Mean Kappa", fill = "Mean Kappa")
      
      ggsave(file.path(OUT_DIR, "map_mean_kappa.png"), p_kap, width = 8, height = 6, dpi = 300)
    }
  }
}

message("Maps written to: ", OUT_DIR)

# ================================
# Report-style map of ecodistrict metrics (single page)
# ================================

# ----------------
# DEFINE INPUTS 
# ----------------
# Coverage threshold for showing as "valid" (greys out below this)
COVERAGE_THRESHOLD   <- 1  # ecodistricts with cells_enough < 1 get greyed out

# Output file
OUT_FIG              <- file.path(OUT_DIR, "report_ecodistrict_metrics.png")

# ----------------
# READ DATA
# ----------------
if (file.exists(OUT_ECODISTRICT_GPKG)) {
  message("Reading GPKG with attributes: ", OUT_ECODISTRICT_GPKG)
  eco_sf <- sf::st_read(OUT_ECODISTRICT_GPKG, quiet = TRUE)
} else {
  message("GPKG not found; reading shapefile + joining CSV…")
  stopifnot(file.exists(ECODISTRICT_SHP), file.exists(OUT_ECODISTRICT_CSV))
  eco_sf <- sf::st_read(ECODISTRICT_SHP, quiet = TRUE)
  stopifnot(ECODISTRICT_ID %in% names(eco_sf))
  eco_tab <- readr::read_csv(OUT_ECODISTRICT_CSV, show_col_types = FALSE)
  # Join on the ecodistrict ID
  eco_sf <- eco_sf %>%
    left_join(eco_tab, by = setNames(ECODISTRICT_ID, ECODISTRICT_ID))
}

# Validate geometry
if (is.null(sf::st_geometry(eco_sf))) stop("No geometry found. Check inputs.")
eco_sf <- sf::st_make_valid(eco_sf)

# ----------------
# METRICS & LABELS
# ----------------
metrics <- c(
  "frac_agree_all",
  "frac_disagree",
  "consensus_presence_mean",
  "uncertainty_mean"
)

metric_labels <- c(
  frac_agree_all           = "Full Agreement (All 0 or All 1)",
  frac_disagree            = "Disagreement (Mixed 0/1)",
  consensus_presence_mean  = "Mean Consensus (Presence Fraction)",
  uncertainty_mean         = "Uncertainty (0 = Full Agreement, Max at 0.5)"
)

missing <- setdiff(metrics, names(eco_sf))
if (length(missing)) {
  warning("Missing expected metrics in ecodistrict attributes: ",
          paste(missing, collapse = ", "))
}

# Grey out ecodistricts with insufficient coverage
if (!("cells_enough" %in% names(eco_sf))) {
  warning("'cells_enough' not found; cannot grey out uncovered areas. Proceeding without masking.")
  eco_sf$coverage_ok <- TRUE
} else {
  eco_sf <- eco_sf %>%
    mutate(coverage_ok = !is.na(cells_enough) & cells_enough >= COVERAGE_THRESHOLD)
}

# Helper: set values to NA where coverage is not OK (so they plot as grey)
for (m in intersect(metrics, names(eco_sf))) {
  eco_sf[[paste0(m, "_masked")]] <- ifelse(eco_sf$coverage_ok, eco_sf[[m]], NA_real_)
}

# ----------------
# THEMING & SCALES
# ----------------
theme_map <- theme_minimal(base_size = 11) +
  theme(
    panel.grid       = element_blank(),
    axis.text        = element_blank(),
    axis.title       = element_blank(),
    axis.ticks       = element_blank(),
    legend.position  = "right",
    legend.title     = element_text(size = 10),
    legend.text      = element_text(size = 9)
  )

# Fixed 0–1 scale for comparability across panels
fill_scale_01 <- scale_fill_viridis(
  option = "C",
  limits = c(0, 1),
  oob    = scales::squish,
  na.value = "grey92"
)

# ----------------
# MAP BUILDER
# ----------------
build_panel <- function(data, metric_masked, panel_title) {
  ggplot(data) +
    geom_sf(aes(fill = .data[[metric_masked]]), color = "grey35", linewidth = 0.1) +
    fill_scale_01 +
    theme_map +
    labs(
      title = panel_title,
      fill  = "Value (0–1)"
    ) #+
  # Optional map furniture (comment out if not desired)
  #annotation_scale(location = "bl", width_hint = 0.25, text_cex = 0.7) +
  #annotation_north_arrow(
  #  location = "bl", which_north = "true",
  #  pad_x = unit(0.5, "cm"), pad_y = unit(0.9, "cm"),
  #  style = north_arrow_fancy_orienteering
  #)
}

p1 <- build_panel(eco_sf, "frac_agree_all_masked",          metric_labels[["frac_agree_all"]])
p2 <- build_panel(eco_sf, "frac_disagree_masked",           metric_labels[["frac_disagree"]])
p3 <- build_panel(eco_sf, "consensus_presence_mean_masked", metric_labels[["consensus_presence_mean"]])
p4 <- build_panel(eco_sf, "uncertainty_mean_masked",        metric_labels[["uncertainty_mean"]])

# Arrange into a 2x2 with shared legend
report_title <- "Ecodistrict Metrics"
report_caption <- if ("cells_enough" %in% names(eco_sf)) {
  "Grey areas indicate ecodistricts without sufficient coverage (cells_enough < 1)."
} else {
  NULL
}

# Patchwork: collect legends so only one legend is shown
layout <- (p1 + p2) / (p3 + p4)
report_fig <- layout +
  plot_layout(guides = "collect") &
  theme(legend.position = "right") #+
#plot_annotation(title = report_title, caption = report_caption)

# Save
ggsave(OUT_FIG, plot = report_fig, width = 12, height = 9, dpi = 300)
message("Wrote report figure: ", OUT_FIG)

# Or can save as vector PDF 
# ggsave(sub("\\.png$", ".pdf", OUT_FIG), plot = report_fig, width = 12, height = 9)
