
###################################################################################
##                   ANALYSIS: SPATIAL AGREEMENT SUMMARIES                       ##
###################################################################################

##----------------------- Alison Long - January 2026 -------------------------

# step 5: run descriptive data summaries to assess similarities and differences
#        among existing grassland datasets. This version compares grassland datasets
#        containing more than one grassland class (native/intact vs. tame/disturbed).

# load packages
library(terra)
library(sf)
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(vegan)
library(dplyr)
library(tidyr)
library(readr)
library(ca)

# Set the temporary directory to a drive with sufficient space
terraOptions(tempdir = "R:/CanadianGrasslands/WG_classification/scratch")

# ============================================================
# 3-class categorical raster comparison (0/1/2) with partial extents
# Pairwise comparisons are computed ONLY on the intersection of
# valid (non-NA) pixels for each raster pair.
#
# Classes:
#   0 = not grass
#   1 = tame grass
#   2 = native grass
#
# Outputs:
# - Pairwise distance matrices:
#     * Simple matching distance (mismatch rate on categories)
#     * 1 - Cohen's kappa (chance-adjusted disagreement)
#     * Optional: one-vs-all Jaccard per class + macro average
# - Pairwise overlap matrices (counts & fractions)
# - PCoA plots (cmdscale) on chosen distance matrix
# - Heatmaps + clustering
# - Optional: Correspondence Analysis (CA) on raster x class table
# - Optional: union-grid coverage/consensus maps
#
# Notes:
# - PCoA supports ordination from arbitrary dissimilarity matrices. [3](https://fundy.ext.unb.ca/metrics.htm)[4](https://researchportal.helsinki.fi/en/publications/measuring-habitat-fragmentation-an-evaluation-of-landscape-patter/)
# - CA is a PCA-like method tailored for nominal/categorical tables. [5](https://ourcodingclub.github.io/tutorials/ordination/)[6](https://academic.oup.com/jpe/article/9/6/649/2623732)
# - Cohen’s kappa provides chance-corrected agreement for categorical labels. [1](https://github.com/andland/logisticPCA)
# - Jaccard is presence-focused; here used optionally one-vs-all per class. [2](https://arxiv.org/html/2405.06559v1)
# ============================================================

# -----------------------
# USER SETTINGS (EDIT)
# -----------------------
RASTER_DIR   <- "R:/CanadianGrasslands/WG_Classification/github/tomorrows-prairies/processeddatasets/Grass/pilot_binary"
RASTER_GLOB  <- "\\.tif$"
META_CSV     <- "R:/CanadianGrasslands/WG_Classification/github/tomorrows-prairies/processeddatasets/Grass/pilot_binary/metadata.csv"  # columns: raster,year,model,source (year optional)
OUT_DIR      <- "R:/CanadianGrasslands/WG_Classification/results/outputs_pairwise_partial_extent_pilot_natvtam"
PILOT_SHP    <- "R:/CanadianGrasslands/processeddatasets/pilotarea/pilot_area_aligned_cnvc_Jan2026.shp"

# Sampling within each pairwise overlap to speed up large rasters:
# NULL = use all overlapping pixels (can be slow/large)
SAMPLE_N_PER_PAIR <- NULL

# If two rasters do not overlap at all, set their distance to:
NO_OVERLAP_DISTANCE <- 1

# Optional Jaccard (one-vs-all) per class and macro-average?
COMPUTE_ONEVALL_JACCARD <- TRUE

# Which distance matrix to use for PCoA/cluster?
# Options: "match" or "kappa" or "jac_macro"
PCOA_DISTANCE <- "match"

# Optional: Correspondence Analysis on raster x class composition?
RUN_CA <- TRUE

# Optional: Build union-grid coverage/consensus maps across all rasters?
MAKE_UNION_CONSENSUS <- TRUE

# Optional: produce a detailed disagreement map for a chosen pair:
PAIR_MAP_A <- NULL  # e.g., "map1.tif" or "map1" (basename without extension)
PAIR_MAP_B <- NULL

SEED <- 42

set.seed(SEED)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# -----------------------
# READ METADATA
# -----------------------
meta <- read_csv(META_CSV, show_col_types = FALSE) %>%
  mutate(
    model  = as.factor(model),
    source = as.factor(source),
    year = if ("year" %in% names(.)) as.POSIXct(year, tz = "UTC") else NA
  )

# -----------------------
# LOAD RASTERS
# -----------------------
files <- list.files(RASTER_DIR, pattern = RASTER_GLOB, full.names = TRUE)
if (length(files) == 0) stop("No rasters found in RASTER_DIR matching RASTER_GLOB.")

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

rlist <- lapply(meta2$file, rast)
names(rlist) <- meta2$raster
n <- length(rlist)
ras_names <- names(rlist)

# -----------------------
# STANDARDIZE TO {0,1,2} + NA
# -----------------------
# Convert to integer; clamp to [0,2].
rlist <- lapply(rlist, function(r) {
  r <- round(r)
  r[r < 0] <- 0
  r[r > 2] <- 2
  r
})

# -----------------------
# PRECOMPUTE VALID PIXEL COUNTS PER RASTER
# -----------------------
valid_counts <- sapply(rlist, function(r) global(!is.na(r), "sum", na.rm = TRUE)[1,1])

# ============================================================
# Helper: extract overlapping pixel vectors for a pair
# ============================================================
extract_pair_overlap <- function(rA, rB, sample_n = NULL) {
  e <- intersect(ext(rA), ext(rB))
  if (is.null(e)) return(list(ok = FALSE))
  
  a <- crop(rA, e, snap = "near")
  b <- crop(rB, e, snap = "near")
  
  # If geometry differs, resample B to A grid in overlap using nearest neighbor
  if (!compareGeom(a, b, stopOnError = FALSE)) {
    b <- resample(b, a, method = "near")
  }
  
  s <- c(a, b)
  
  if (!is.null(sample_n)) {
    df <- spatSample(s, size = sample_n, method = "random",
                     na.rm = TRUE, as.df = TRUE, values = TRUE)
    df <- df[complete.cases(df), , drop = FALSE]
    if (nrow(df) == 0) return(list(ok = FALSE))
    va <- df[[1]]
    vb <- df[[2]]
  } else {
    v <- values(s)
    v <- v[complete.cases(v), , drop = FALSE]
    if (nrow(v) == 0) return(list(ok = FALSE))
    va <- v[,1]
    vb <- v[,2]
  }
  
  list(ok = TRUE, va = va, vb = vb, n = length(va))
}

# ============================================================
# Helper: multiclass metrics from vectors
# ============================================================
metrics_3class <- function(a, b, compute_jaccard = TRUE) {
  # a,b are {0,1,2} with no NA
  # Confusion table 3x3 (ensure all levels included)
  tab <- table(factor(a, levels = 0:2), factor(b, levels = 0:2))
  N <- sum(tab)
  agree <- sum(diag(tab))
  
  # Simple matching distance (mismatch rate on categories)
  d_match <- 1 - (agree / N)
  
  # Cohen's kappa (multiclass) computed from confusion table
  # kappa = (po - pe) / (1 - pe)
  po <- agree / N
  rowm <- rowSums(tab) / N
  colm <- colSums(tab) / N
  pe <- sum(rowm * colm)
  kappa <- if ((1 - pe) == 0) NA_real_ else (po - pe) / (1 - pe)
  d_kappa <- if (is.na(kappa)) NA_real_ else (1 - kappa)
  
  out <- list(
    N = N,
    confusion = tab,
    d_match = d_match,
    kappa = kappa,
    d_kappa = d_kappa
  )
  
  # Optional: one-vs-all Jaccard distances per class, plus macro-average
  # For class c: treat c as "presence" (1) and others as "absence" (0).
  if (compute_jaccard) {
    d_j <- numeric(3)
    names(d_j) <- c("notgrass(0)", "tame(1)", "native(2)")
    for (c in 0:2) {
      ac <- as.integer(a == c)
      bc <- as.integer(b == c)
      TP <- sum(ac == 1 & bc == 1)
      FP <- sum(ac == 1 & bc == 0)
      FN <- sum(ac == 0 & bc == 1)
      denom <- TP + FP + FN
      # If denom==0, neither map ever uses that class in overlap -> treat as identical for that class
      d_j[c+1] <- if (denom == 0) 0 else (1 - TP / denom)
    }
    out$d_jaccard_onevsall <- d_j
    out$d_jaccard_macro <- mean(d_j)
  }
  
  out
}

# ============================================================
# Pairwise matrices
# ============================================================
D_match <- matrix(NA_real_, n, n, dimnames = list(ras_names, ras_names))
D_kappa <- matrix(NA_real_, n, n, dimnames = list(ras_names, ras_names))
K_mat   <- matrix(NA_real_, n, n, dimnames = list(ras_names, ras_names))

D_jac0  <- matrix(NA_real_, n, n, dimnames = list(ras_names, ras_names))
D_jac1  <- matrix(NA_real_, n, n, dimnames = list(ras_names, ras_names))
D_jac2  <- matrix(NA_real_, n, n, dimnames = list(ras_names, ras_names))
D_jac_macro <- matrix(NA_real_, n, n, dimnames = list(ras_names, ras_names))

N_overlap <- matrix(0, n, n, dimnames = list(ras_names, ras_names))
Frac_overlap_i <- matrix(NA_real_, n, n, dimnames = list(ras_names, ras_names))
Frac_overlap_j <- matrix(NA_real_, n, n, dimnames = list(ras_names, ras_names))

message("Computing pairwise metrics for ", n, " rasters (3-class)...")

for (i in seq_len(n)) {
  # diagonal
  D_match[i,i] <- 0
  D_kappa[i,i] <- 0
  K_mat[i,i] <- 1
  N_overlap[i,i] <- valid_counts[i]
  Frac_overlap_i[i,i] <- 1
  Frac_overlap_j[i,i] <- 1
  
  if (COMPUTE_ONEVALL_JACCARD) {
    D_jac0[i,i] <- 0
    D_jac1[i,i] <- 0
    D_jac2[i,i] <- 0
    D_jac_macro[i,i] <- 0
  }
  
  for (j in (i+1):n) {
    ov <- extract_pair_overlap(rlist[[i]], rlist[[j]], sample_n = SAMPLE_N_PER_PAIR)
    
    if (!isTRUE(ov$ok)) {
      # No overlap or no valid pixels
      D_match[i,j] <- D_match[j,i] <- NO_OVERLAP_DISTANCE
      D_kappa[i,j] <- D_kappa[j,i] <- NA_real_
      K_mat[i,j]   <- K_mat[j,i]   <- NA_real_
      N_overlap[i,j] <- N_overlap[j,i] <- 0
      Frac_overlap_i[i,j] <- Frac_overlap_i[j,i] <- 0
      Frac_overlap_j[i,j] <- Frac_overlap_j[j,i] <- 0
      
      if (COMPUTE_ONEVALL_JACCARD) {
        D_jac0[i,j] <- D_jac0[j,i] <- NO_OVERLAP_DISTANCE
        D_jac1[i,j] <- D_jac1[j,i] <- NO_OVERLAP_DISTANCE
        D_jac2[i,j] <- D_jac2[j,i] <- NO_OVERLAP_DISTANCE
        D_jac_macro[i,j] <- D_jac_macro[j,i] <- NO_OVERLAP_DISTANCE
      }
      next
    }
    
    m <- metrics_3class(ov$va, ov$vb, compute_jaccard = COMPUTE_ONEVALL_JACCARD)
    
    D_match[i,j] <- D_match[j,i] <- m$d_match
    D_kappa[i,j] <- D_kappa[j,i] <- m$d_kappa
    K_mat[i,j]   <- K_mat[j,i]   <- m$kappa
    
    N_overlap[i,j] <- N_overlap[j,i] <- m$N
    Frac_overlap_i[i,j] <- Frac_overlap_j[j,i] <- m$N / valid_counts[i]
    Frac_overlap_j[i,j] <- Frac_overlap_i[j,i] <- m$N / valid_counts[j]
    
    if (COMPUTE_ONEVALL_JACCARD) {
      D_jac0[i,j] <- D_jac0[j,i] <- m$d_jaccard_onevsall[1]
      D_jac1[i,j] <- D_jac1[j,i] <- m$d_jaccard_onevsall[2]
      D_jac2[i,j] <- D_jac2[j,i] <- m$d_jaccard_onevsall[3]
      D_jac_macro[i,j] <- D_jac_macro[j,i] <- m$d_jaccard_macro
    }
  }
}

# ============================================================
# Save outputs
# ============================================================
saveRDS(D_match, file.path(OUT_DIR, "pairwise_distance_simple_matching.rds"))
saveRDS(D_kappa, file.path(OUT_DIR, "pairwise_distance_1_minus_kappa.rds"))
saveRDS(K_mat,   file.path(OUT_DIR, "pairwise_kappa.rds"))

saveRDS(N_overlap, file.path(OUT_DIR, "pairwise_overlap_n.rds"))
saveRDS(Frac_overlap_i, file.path(OUT_DIR, "pairwise_overlap_frac_of_i.rds"))
saveRDS(Frac_overlap_j, file.path(OUT_DIR, "pairwise_overlap_frac_of_j.rds"))

write_csv(as.data.frame(D_match), file.path(OUT_DIR, "pairwise_distance_simple_matching.csv"))
write_csv(as.data.frame(D_kappa), file.path(OUT_DIR, "pairwise_distance_1_minus_kappa.csv"))
write_csv(as.data.frame(K_mat),   file.path(OUT_DIR, "pairwise_kappa.csv"))