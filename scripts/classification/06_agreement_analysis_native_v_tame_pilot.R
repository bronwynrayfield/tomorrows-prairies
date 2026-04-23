
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
library(pheatmap)
library(vegan)
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
# USER SETTINGS 
# -----------------------
RASTER_DIR   <- "R:/CanadianGrasslands/WG_Classification/github/tomorrows-prairies/processeddatasets/Grass/pilot_nattam"
RASTER_GLOB  <- "\\.tif$"
META_CSV     <- "R:/CanadianGrasslands/WG_Classification/github/tomorrows-prairies/processeddatasets/Grass/pilot_nattam/metadata.csv"  # columns: raster,year,model,source (year optional)
OUT_DIR      <- "R:/CanadianGrasslands/WG_Classification/results/spatial_agreement_analysis/native_v_tame"
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
PAIR_MAP_A <- "cnvc_algn_crop_pa_natvtam"  # e.g., "map1.tif" or "map1" (basename without extension)
PAIR_MAP_B <- "ngi_merged_reproj_algn_crop_pa_natvtam"

SEED <- 42

set.seed(SEED)
#dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# -----------------------
# READ METADATA
# -----------------------
meta <- read_csv(META_CSV, show_col_types = FALSE) %>%
  mutate(
    raster = as.character(raster),
    model  = as.factor(model),
    source_simp = as.factor(source_simp),
    year = if ("year" %in% names(.)) {
      y_chr <- as.character(year)
      # If full dates appear, keep only the first 4 chars; blank -> NA
      y4    <- substr(y_chr, 1, 4)
      y_int <- suppressWarnings(as.integer(y4))
      ifelse(is.na(y_int) | y4 == "" | y4 == "NA", NA_integer_, y_int)
    } else {
      NA_integer_
    }
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

# Attach file paths to metadata (supporting raster name with/without extension)
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

message("Loading ", nrow(meta2), " rasters...")
rlist <- lapply(meta2$file, rast)
names(rlist) <- meta2$raster
n <- length(rlist)
ras_names <- names(rlist)

# -----------------------
# OPTIONAL: CLIP TO PILOT AREA
# -----------------------
if (!is.null(PILOT_SHP)) {
  message("Clipping all rasters to pilot area: ", PILOT_SHP)
  # use terra::vect for better integration
  pilot <- tryCatch(vect(PILOT_SHP), error = function(e) NULL)
  if (is.null(pilot)) {
    warning("Could not read PILOT_SHP; proceeding without clipping.")
  } else {
    # Reproject pilot to first raster CRS if needed
    r1_crs <- crs(rlist[[1]], proj=TRUE)
    if (!same.crs(pilot, rlist[[1]])) {
      pilot <- project(pilot, r1_crs)
    }
    # crop + mask each raster
    rlist <- lapply(rlist, function(r) {
      rr <- suppressWarnings(crop(r, pilot))
      rr <- suppressWarnings(mask(rr, pilot))
      rr
    })
  }
}

diag_empty <- tibble(
  raster = names(rlist),
  ncol = purrr::map_int(rlist, ~ suppressWarnings(tryCatch(ncol(.x), error = function(e) NA_integer_))),
  nrow = purrr::map_int(rlist, ~ suppressWarnings(tryCatch(nrow(.x), error = function(e) NA_integer_)))
) %>%
  mutate(is_empty = is.na(ncol) | is.na(nrow) | ncol == 0 | nrow == 0)

print(diag_empty %>% filter(is_empty))

# -----------------------
# CONFIRM SAME GRID AND EXTENT
# -----------------------
# ===========================
# Alignment audit for rlist
# ===========================
stopifnot(length(rlist) > 0)

# Helper readers that are NA-safe and type-stable
get_crs   <- function(r) tryCatch(crs(r, proj=TRUE), error=function(e) NA_character_)
get_nrow  <- function(r) tryCatch(as.integer(nrow(r)), error=function(e) NA_integer_)
get_ncol  <- function(r) tryCatch(as.integer(ncol(r)), error=function(e) NA_integer_)
get_xres  <- function(r) tryCatch(as.numeric(res(r)[1]), error=function(e) NA_real_)
get_yres  <- function(r) tryCatch(as.numeric(res(r)[2]), error=function(e) NA_real_)
get_xmin  <- function(r) tryCatch(as.numeric(xmin(ext(r))), error=function(e) NA_real_)
get_xmax  <- function(r) tryCatch(as.numeric(xmax(ext(r))), error=function(e) NA_real_)
get_ymin  <- function(r) tryCatch(as.numeric(ymin(ext(r))), error=function(e) NA_real_)
get_ymax  <- function(r) tryCatch(as.numeric(ymax(ext(r))), error=function(e) NA_real_)
get_ox    <- function(r) tryCatch(as.numeric(origin(r)[1]), error=function(e) NA_real_)
get_oy    <- function(r) tryCatch(as.numeric(origin(r)[2]), error=function(e) NA_real_)
get_nlyr  <- function(r) tryCatch(as.integer(nlyr(r)), error=function(e) NA_integer_)

# Ensure all rasters have non-NA names for reporting
if (is.null(names(rlist))) names(rlist) <- paste0("r", seq_along(rlist))
bad_nm <- which(is.na(names(rlist)) | !nzchar(names(rlist)))
if (length(bad_nm) > 0) {
  names(rlist)[bad_nm] <- tools::file_path_sans_ext(basename(meta2$file[bad_nm]))
}

# Build report
grid_report <- tibble(
  raster   = names(rlist),
  crs      = purrr::map_chr(rlist, get_crs),
  nrow     = purrr::map_int(rlist, get_nrow),
  ncol     = purrr::map_int(rlist, get_ncol),
  xres     = purrr::map_dbl(rlist, get_xres),
  yres     = purrr::map_dbl(rlist, get_yres),
  xmin     = purrr::map_dbl(rlist, get_xmin),
  xmax     = purrr::map_dbl(rlist, get_xmax),
  ymin     = purrr::map_dbl(rlist, get_ymin),
  ymax     = purrr::map_dbl(rlist, get_ymax),
  origin_x = purrr::map_dbl(rlist, get_ox),
  origin_y = purrr::map_dbl(rlist, get_oy),
  nlyr     = purrr::map_int(rlist, get_nlyr)
)

print(grid_report, n = nrow(grid_report))

# Reference (first raster)
r_ref <- rlist[[1]]
ref_vals <- list(
  crs   = get_crs(r_ref),
  nrow  = get_nrow(r_ref),
  ncol  = get_ncol(r_ref),
  xres  = get_xres(r_ref),
  yres  = get_yres(r_ref),
  xmin  = get_xmin(r_ref),
  xmax  = get_xmax(r_ref),
  ymin  = get_ymin(r_ref),
  ymax  = get_ymax(r_ref),
  ox    = get_ox(r_ref),
  oy    = get_oy(r_ref)
)

mismatch <- grid_report %>%
  mutate(
    crs_ok   = crs   == ref_vals$crs,
    nrow_ok  = nrow  == ref_vals$nrow,
    ncol_ok  = ncol  == ref_vals$ncol,
    xres_ok  = xres  == ref_vals$xres,
    yres_ok  = yres  == ref_vals$yres,
    xmin_ok  = xmin  == ref_vals$xmin,
    xmax_ok  = xmax  == ref_vals$xmax,
    ymin_ok  = ymin  == ref_vals$ymin,
    ymax_ok  = ymax  == ref_vals$ymax,
    ox_ok    = origin_x == ref_vals$ox,
    oy_ok    = origin_y == ref_vals$oy
  )

# Summaries
message("All CRS identical? ", all(mismatch$crs_ok,  na.rm = FALSE))
message("All nrow identical? ", all(mismatch$nrow_ok, na.rm = FALSE))
message("All ncol identical? ", all(mismatch$ncol_ok, na.rm = FALSE))
message("All res identical? x: ", all(mismatch$xres_ok, na.rm = FALSE),
        " | y: ", all(mismatch$yres_ok, na.rm = FALSE))
message("All extents identical? xmin: ", all(mismatch$xmin_ok, na.rm = FALSE),
        " xmax: ", all(mismatch$xmax_ok, na.rm = FALSE),
        " ymin: ", all(mismatch$ymin_ok, na.rm = FALSE),
        " ymax: ", all(mismatch$ymax_ok, na.rm = FALSE))
message("All origins identical? ox: ", all(mismatch$ox_ok, na.rm = FALSE),
        " oy: ", all(mismatch$oy_ok, na.rm = FALSE))

# Differences (if any)
differences <- mismatch %>%
  filter(!(crs_ok & nrow_ok & ncol_ok & xres_ok & yres_ok & xmin_ok & xmax_ok & ymin_ok & ymax_ok & ox_ok & oy_ok))

if (nrow(differences) == 0) {
  message("✅ All rasters share the same grid and extent.")
} else {
  message("⚠️ Some rasters differ from the reference grid:")
  print(differences %>%
          select(raster, ends_with("_ok"), nrow, ncol, xres, yres, xmin, xmax, ymin, ymax, origin_x, origin_y))
}

# Optional: strict pairwise geometry check
pair_bad <- list()
for (i in seq_along(rlist)) {
  for (j in seq.int(i+1, length(rlist))) {
    ok <- tryCatch(compareGeom(rlist[[i]], rlist[[j]], stopOnError = FALSE), error=function(e) FALSE)
    if (!isTRUE(ok)) {
      pair_bad[[length(pair_bad)+1]] <- c(names(rlist)[i], names(rlist)[j])
    }
  }
}
if (length(pair_bad) == 0) {
  message("✅ All pairwise compareGeom checks passed.")
} else {
  message("⚠️ These pairs fail compareGeom():")
  print(do.call(rbind, pair_bad))
}

# ---- Ensure non-NA, unique raster names ----
make_safe_names <- function(nms, files) {
  out <- nms
  # Fill missing/empty with basenames
  bad <- which(is.na(out) | !nzchar(out))
  if (length(bad) > 0) {
    fallback <- tools::file_path_sans_ext(basename(files[bad]))
    out[bad] <- fallback
  }
  # Enforce uniqueness (append index suffix only to duplicates)
  if (anyDuplicated(out)) {
    dup_seen <- list()
    for (k in seq_along(out)) {
      nm <- out[k]
      if (is.na(nm) || !nzchar(nm)) nm <- paste0("r", k)
      if (!is.null(dup_seen[[nm]])) {
        dup_seen[[nm]] <- dup_seen[[nm]] + 1
        out[k] <- paste0(nm, "_", dup_seen[[nm]])
      } else {
        dup_seen[[nm]] <- 1
        out[k] <- nm
      }
    }
  }
  out
}




# -----------------------
# STANDARDIZE TO {0,1,2} + NA
# -----------------------
rlist <- lapply(rlist, function(r) {
  r <- round(r)
  r[r < 0] <- 0
  r[r > 2] <- 2
  r
})

# ============================================================
# Helper: determine if a SpatRaster is empty (0 rows or 0 cols)
# ============================================================
is_empty_raster <- function(r) {
  if (is.null(r)) return(TRUE)
  nc <- suppressWarnings(try(ncol(r), silent = TRUE))
  nr <- suppressWarnings(try(nrow(r), silent = TRUE))
  if (inherits(nc, "try-error") || inherits(nr, "try-error")) return(TRUE)
  is.na(nc) || is.na(nr) || nc == 0 || nr == 0
}

# ----- DROP EMPTY RASTERS (0x0) BEFORE PAIRWISE COMPARISONS -----
empty_idx <- which(vapply(rlist, is_empty_raster, logical(1)))
if (length(empty_idx) > 0) {
  message("Dropping ", length(empty_idx), " empty raster(s) (likely outside clip): ",
          paste(names(rlist)[empty_idx], collapse = ", "))
  keep_idx <- setdiff(seq_along(rlist), empty_idx)
  rlist    <- rlist[keep_idx]
  meta2    <- meta2[keep_idx, , drop = FALSE]
}
# Refresh names and n after dropping empties
ras_names <- names(rlist)
n <- length(rlist)
if (n < 2) stop("Need at least two non-empty rasters to compare after clipping/standardization.")

# -----------------------
# PRECOMPUTE VALID PIXEL COUNTS PER RASTER
# -----------------------
valid_counts <- sapply(rlist, function(r) {
  v <- tryCatch(global(!is.na(r), "sum", na.rm = TRUE)[1,1], error = function(e) NA_real_)
  if (is.na(v)) 0 else v
})

# ============================================================
# Helper: extract overlapping pixel vectors for a pair
# ============================================================
extract_pair_overlap <- function(rA, rB, sample_n = NULL) {
  # If either raster is empty, bail early
  if (is_empty_raster(rA) || is_empty_raster(rB)) return(list(ok = FALSE))
  
  # ext() can error when rasters are malformed; wrap it
  e <- tryCatch({
    intersect(ext(rA), ext(rB))
  }, error = function(e) NULL)
  
  if (is.null(e)) return(list(ok = FALSE))
  
  a <- tryCatch(crop(rA, e, snap = "near"), error = function(e) NULL)
  b <- tryCatch(crop(rB, e, snap = "near"), error = function(e) NULL)
  if (is.null(a) || is.null(b) || is_empty_raster(a) || is_empty_raster(b)) {
    return(list(ok = FALSE))
  }
  
  if (!compareGeom(a, b, stopOnError = FALSE)) {
    b <- tryCatch(resample(b, a, method = "near"), error = function(e) NULL)
    if (is.null(b) || is_empty_raster(b)) return(list(ok = FALSE))
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
  tab <- table(factor(a, levels = 0:2), factor(b, levels = 0:2))
  N <- sum(tab)
  agree <- sum(diag(tab))
  
  # Simple matching distance (mismatch rate on categories)
  d_match <- 1 - (agree / N)
  
  # Cohen's kappa (multiclass)
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
  
  # Optional: one-vs-all Jaccard distances per class + macro average
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
      d_j[c+1] <- if (denom == 0) 0 else (1 - TP / denom)
    }
    out$d_jaccard_onevsall <- d_j
    out$d_jaccard_macro <- mean(d_j)
  }
  
  out
}

# =====================================================================
# (4) Pairwise geometry audit with safe names (replaces your old audit)
# =====================================================================
safe_nm <- function(idx) {
  nm <- names(rlist)[idx]
  if (is.na(nm) || !nzchar(nm)) paste0("r", idx) else nm
}

pair_bad <- list()
n <- length(rlist)
for (i in seq_len(n)) {
  if (i < n) {  # <-- critical guard: skip inner loop when i==n
    for (j in (i+1):n) {
      ok <- tryCatch(compareGeom(rlist[[i]], rlist[[j]], stopOnError = FALSE),
                     error = function(e) FALSE)
      if (!isTRUE(ok)) {
        pair_bad[[length(pair_bad)+1]] <- c(safe_nm(i), safe_nm(j))
      }
    }
  }
}
if (length(pair_bad) == 0) {
  message("✅ All pairwise compareGeom checks passed.")
} else {
  message("⚠️ These pairs fail compareGeom():")
  print(do.call(rbind, pair_bad))
}

# =====================================================================
# (5) Rebuild matrices with final names + SAFE name-based pairwise loop
# =====================================================================

# Reconfirm/repair names right before building matrices
make_safe_names <- function(nms, files) {
  out <- nms
  bad <- which(is.na(out) | !nzchar(out))
  if (length(bad) > 0) {
    fallback <- tools::file_path_sans_ext(basename(files[bad]))
    out[bad] <- fallback
  }
  # Enforce uniqueness (append stable suffix only to duplicates)
  if (anyDuplicated(out)) {
    dup_seen <- list()
    for (k in seq_along(out)) {
      nm <- out[k]
      if (is.na(nm) || !nzchar(nm)) nm <- paste0("r", k)
      if (!is.null(dup_seen[[nm]])) {
        dup_seen[[nm]] <- dup_seen[[nm]] + 1
        out[k] <- paste0(nm, "_", dup_seen[[nm]])
      } else {
        dup_seen[[nm]] <- 1
        out[k] <- nm
      }
    }
  }
  out
}

ras_names <- make_safe_names(names(rlist), meta2$file)
names(rlist) <- ras_names

# Recompute counts now that rlist is final
valid_counts <- sapply(rlist, function(r) {
  v <- tryCatch(global(!is.na(r), "sum", na.rm = TRUE)[1,1], error = function(e) NA_real_)
  if (is.na(v)) 0 else v
})

# Build matrices sized to the final n with matching dimnames
ras_names <- names(rlist)
n <- length(rlist)
stopifnot(n >= 2)

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
  ri <- ras_names[i]  # name of row i
  
  # Diagonal (name-based indexing)
  D_match[ri, ri] <- 0
  D_kappa[ri, ri] <- 0
  K_mat[ri, ri]   <- 1
  N_overlap[ri, ri] <- valid_counts[i]
  Frac_overlap_i[ri, ri] <- 1
  Frac_overlap_j[ri, ri] <- 1
  
  if (COMPUTE_ONEVALL_JACCARD) {
    D_jac0[ri, ri] <- 0
    D_jac1[ri, ri] <- 0
    D_jac2[ri, ri] <- 0
    D_jac_macro[ri, ri] <- 0
  }
  
  if (i < n) {  # <-- critical guard: skip inner loop when i==n
    for (j in (i+1):n) {
      rj <- ras_names[j]  # name of col j
      
      ov <- tryCatch(
        extract_pair_overlap(rlist[[i]], rlist[[j]], sample_n = SAMPLE_N_PER_PAIR),
        error = function(e) {
          message("extract_pair_overlap() error on pair: ", ri, " vs ", rj, " :: ", conditionMessage(e))
          structure(list(ok = FALSE), class = "list")
        }
      )
      
      if (!isTRUE(ov$ok)) {
        # No overlap or invalid => fill with NO_OVERLAP_DISTANCE / NA (separate assignments)
        D_match[ri, rj] <- NO_OVERLAP_DISTANCE
        D_match[rj, ri] <- NO_OVERLAP_DISTANCE
        D_kappa[ri, rj] <- NA_real_
        D_kappa[rj, ri] <- NA_real_
        K_mat[ri, rj]   <- NA_real_
        K_mat[rj, ri]   <- NA_real_
        N_overlap[ri, rj] <- 0
        N_overlap[rj, ri] <- 0
        Frac_overlap_i[ri, rj] <- 0
        Frac_overlap_i[rj, ri] <- 0
        Frac_overlap_j[ri, rj] <- 0
        Frac_overlap_j[rj, ri] <- 0
        if (COMPUTE_ONEVALL_JACCARD) {
          D_jac0[ri, rj] <- NO_OVERLAP_DISTANCE
          D_jac0[rj, ri] <- NO_OVERLAP_DISTANCE
          D_jac1[ri, rj] <- NO_OVERLAP_DISTANCE
          D_jac1[rj, ri] <- NO_OVERLAP_DISTANCE
          D_jac2[ri, rj] <- NO_OVERLAP_DISTANCE
          D_jac2[rj, ri] <- NO_OVERLAP_DISTANCE
          D_jac_macro[ri, rj] <- NO_OVERLAP_DISTANCE
          D_jac_macro[rj, ri] <- NO_OVERLAP_DISTANCE
        }
        next
      }
      
      m <- metrics_3class(ov$va, ov$vb, compute_jaccard = COMPUTE_ONEVALL_JACCARD)
      
      D_match[ri, rj] <- m$d_match
      D_match[rj, ri] <- m$d_match
      D_kappa[ri, rj] <- m$d_kappa
      D_kappa[rj, ri] <- m$d_kappa
      K_mat[ri, rj]   <- m$kappa
      K_mat[rj, ri]   <- m$kappa
      
      N_overlap[ri, rj] <- m$N
      N_overlap[rj, ri] <- m$N
      Frac_overlap_i[ri, rj] <- m$N / max(1, valid_counts[i])
      Frac_overlap_i[rj, ri] <- m$N / max(1, valid_counts[j])
      Frac_overlap_j[ri, rj] <- m$N / max(1, valid_counts[j])
      Frac_overlap_j[rj, ri] <- m$N / max(1, valid_counts[i])
      
      if (COMPUTE_ONEVALL_JACCARD) {
        D_jac0[ri, rj] <- m$d_jaccard_onevsall[1]
        D_jac0[rj, ri] <- m$d_jaccard_onevsall[1]
        D_jac1[ri, rj] <- m$d_jaccard_onevsall[2]
        D_jac1[rj, ri] <- m$d_jaccard_onevsall[2]
        D_jac2[ri, rj] <- m$d_jaccard_onevsall[3]
        D_jac2[rj, ri] <- m$d_jaccard_onevsall[3]
        D_jac_macro[ri, rj] <- m$d_jaccard_macro
        D_jac_macro[rj, ri] <- m$d_jaccard_macro
      }
    }
  }
}

# ============================================================
# Save core outputs (RDS + CSV)
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

write_csv(as.data.frame(N_overlap),     file.path(OUT_DIR, "pairwise_overlap_n.csv"))
write_csv(as.data.frame(Frac_overlap_i), file.path(OUT_DIR, "pairwise_overlap_frac_of_i.csv"))
write_csv(as.data.frame(Frac_overlap_j), file.path(OUT_DIR, "pairwise_overlap_frac_of_j.csv"))

if (COMPUTE_ONEVALL_JACCARD) {
  saveRDS(D_jac0,  file.path(OUT_DIR, "pairwise_distance_jaccard_notgrass0.rds"))
  saveRDS(D_jac1,  file.path(OUT_DIR, "pairwise_distance_jaccard_tame1.rds"))
  saveRDS(D_jac2,  file.path(OUT_DIR, "pairwise_distance_jaccard_native2.rds"))
  saveRDS(D_jac_macro, file.path(OUT_DIR, "pairwise_distance_jaccard_macro.rds"))
  
  write_csv(as.data.frame(D_jac0),  file.path(OUT_DIR, "pairwise_distance_jaccard_notgrass0.csv"))
  write_csv(as.data.frame(D_jac1),  file.path(OUT_DIR, "pairwise_distance_jaccard_tame1.csv"))
  write_csv(as.data.frame(D_jac2),  file.path(OUT_DIR, "pairwise_distance_jaccard_native2.csv"))
  write_csv(as.data.frame(D_jac_macro), file.path(OUT_DIR, "pairwise_distance_jaccard_macro.csv"))
}

# ============================================================
# Heatmaps (PNG) for selected matrices
# ============================================================
hm_out <- function(mat, fname, main=NULL) {
  png(file.path(OUT_DIR, fname), width = 1600, height = 1400, res = 200)
  pheatmap(mat, clustering_distance_rows = "euclidean",
           clustering_distance_cols = "euclidean",
           clustering_method = "average",
           main = main, display_numbers = FALSE, angle_col = 45)
  dev.off()
}

hm_out(D_match, "heatmap_simple_matching.png", "Simple Matching Distance (0..1, lower=better)")
if (any(!is.na(D_kappa))) {
  hm_out(D_kappa, "heatmap_1_minus_kappa.png", "1 - Cohen's kappa (lower=better, NA=no overlap)")
}
if (COMPUTE_ONEVALL_JACCARD) {
  hm_out(D_jac_macro, "heatmap_jaccard_macro.png", "Jaccard Macro (one-vs-all, lower=better)")
}
dev.cur()
dev.off()


png_path <- file.path(OUT_DIR, "heatmap_simple_matching.png")
print(file.info(png_path)$size)

# Open via default browser/viewer (avoids Explorer locks/preview)
browseURL(png_path)


summary(as.vector(D_match))
table(is.na(D_match))
unique_vals <- sort(unique(as.vector(D_match)))
print(unique_vals)


pheatmap::pheatmap(
  D_match,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "average",
  main = "Simple Matching Distance (0..1, lower=better)",
  display_numbers = FALSE,
  angle_col = 45,
  na_col = "grey90",
  color = colorRampPalette(c("#08306B", "#4292C6", "#FEF0D9", "#FB6A4A", "#99000D"))(200)
)

pheatmap::pheatmap(
  D_match,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "average",
  main = "Simple Matching Distance (0..1, lower=better)",
  display_numbers = FALSE,
  angle_col = 45,
  na_col = "grey90",
  color = colorRampPalette(c("#08306B", "#4292C6", "#FEF0D9", "#FB6A4A", "#99000D"))(200),
  filename = file.path(OUT_DIR, "heatmap_simple_matching.png"),
  width = 1600/200,  # inches = pixels/dpi
  height = 1400/200
)

# ============================================================
# PCoA (cmdscale) on selected distance
# ============================================================
get_dist_matrix <- function(which = PCOA_DISTANCE) {
  if (which == "match") return(D_match)
  if (which == "kappa") return(D_kappa)
  if (which == "jac_macro") return(D_jac_macro)
  stop("PCOA_DISTANCE must be 'match', 'kappa', or 'jac_macro'.")
}
D_use <- get_dist_matrix(PCOA_DISTANCE)

# Keep only rows/cols with finite (non-NA) distances
keep <- which(apply(D_use, 1, function(x) all(is.finite(x) | is.na(x))))
# Remove rows/cols that contain any NA outside diagonal
# (More conservative: keep only rasters with no NA in row/col)
nonna <- which(apply(D_use, 1, function(x) all(is.finite(x))))
use_idx <- intersect(keep, nonna)
D_sub <- D_use[use_idx, use_idx, drop=FALSE]

if (nrow(D_sub) >= 3) {
  # Force symmetry and zero diagonal, then dist object
  diag(D_sub) <- 0
  D_sub <- (D_sub + t(D_sub)) / 2
  d <- as.dist(D_sub)
  pc <- cmdscale(d, k = 2, eig = TRUE)
  pc_df <- tibble(
    raster = rownames(D_sub),
    Axis1 = pc$points[,1],
    Axis2 = pc$points[,2]
  ) %>%
    left_join(meta2 %>% select(raster, model, source_simp, year), by = "raster") %>%
    mutate(
      # Keep year as factor for distinct colors (or as.integer for a gradient)
      year = factor(year),
      label_short = sub("_.*$", "", raster)
    )
  gg <- ggplot(
    pc_df %>% mutate(label_short = sub("_.*$", "", raster)),
    aes(Axis1, Axis2, label = label_short, color = year, shape = source_simp)
  ) +
    
    geom_hline(yintercept = 0, linetype = 3, color = "grey60") +
    geom_vline(xintercept = 0, linetype = 3, color = "grey60") +
    geom_point(size = 3) +
    ggrepel::geom_text_repel(size = 3, max.overlaps = 20) +
    theme_minimal(base_size = 12) +
    labs(
      title = paste0("PCoA on '", PCOA_DISTANCE, "' distance"),
      subtitle = "Lower distances = more similar",
      color = "Year", shape = "Source"
    )
  ggsave(filename = file.path(OUT_DIR, paste0("PCoA_", PCOA_DISTANCE, "_year", ".png")), gg, width = 9, height = 7, dpi = 220)
} else {
  warning("PCoA skipped: not enough rasters with finite distances after NA filtering.")
}

pc_df %>% select(raster, year) %>% print(n=Inf)

unique(pc_df$year)
table(pc_df$year, useNA="ifany")


if (nrow(D_sub) >= 3) {
  diag(D_sub) <- 0
  D_sub <- (D_sub + t(D_sub)) / 2
  d <- as.dist(D_sub)
  pc <- cmdscale(d, k = 2, eig = TRUE)
  
  pc_df <- tibble(
    raster = rownames(D_sub),
    Axis1  = pc$points[,1],
    Axis2  = pc$points[,2]
  ) %>%
    left_join(meta2 %>% select(raster, model, source_simp, year), by = "raster") %>%
    mutate(
      year = factor(year),
      label_short = sub("_.*$", "", raster)
    )
  
  src_levels <- levels(factor(pc_df$source_simp))
  base_shapes <- c(16, 17, 15, 3, 7, 8)
  shape_vals <- setNames(base_shapes[(seq_along(src_levels)-1) %% length(base_shapes) + 1], src_levels)
  
  gg <- ggplot(pc_df, aes(x = Axis1, y = Axis2)) +
    geom_hline(yintercept = 0, linetype = 3, color = "grey60") +
    geom_vline(xintercept = 0, linetype = 3, color = "grey60") +
    geom_point(aes(color = year, shape = source_simp), size = 4, show.legend = TRUE, na.rm = TRUE) +
    ggrepel::geom_text_repel(
      data = pc_df,
      aes(x = Axis1, y = Axis2, label = label_short),
      inherit.aes = FALSE,
      size = 4,
      max.overlaps = 30,
      show.legend = FALSE,
      na.rm = TRUE
    ) +
    theme_minimal(base_size = 12) +
    scale_shape_manual(values = shape_vals, drop = FALSE) +
    scale_color_brewer(palette = "Dark2", na.value = "grey60") +
    labs(
      title = paste0("PCoA on '", PCOA_DISTANCE, "' distance"),
      subtitle = "Lower distances = more similar",
      color = "Year", shape = "Source"
    )
  
  ggsave(file.path(OUT_DIR, paste0("PCoA_", PCOA_DISTANCE, "_year.png")), gg, width = 9, height = 7, dpi = 220)
} else {
  warning("PCoA skipped: not enough rasters with finite distances after NA filtering.")
}



# PCoA that plots source and model:

D_use <- get_dist_matrix(PCOA_DISTANCE)

# Conservative keep: rows/cols with finite values
nonna <- which(apply(D_use, 1, function(x) all(is.finite(x))))
use_idx <- nonna
D_sub <- D_use[use_idx, use_idx, drop = FALSE]

if (nrow(D_sub) >= 3) {
  diag(D_sub) <- 0
  D_sub <- (D_sub + t(D_sub)) / 2
  d <- as.dist(D_sub)
  pc <- cmdscale(d, k = 2, eig = TRUE)
  
  pc_df <- tibble::tibble(
    raster = rownames(D_sub),
    Axis1 = pc$points[,1],
    Axis2 = pc$points[,2]
  ) %>%
    left_join(meta2 %>% select(raster, model, source_simp, year), by = "raster") %>%
    mutate(
      model  = as.factor(ifelse(is.na(model),  "UnknownModel",  as.character(model))),
      source_simp = as.factor(ifelse(is.na(source_simp), "UnknownSource", as.character(source_simp)))
    )
  
  # If >6 sources, set manual shapes (recycle) to avoid palette warning
  src_levels <- levels(pc_df$source_simp)
  base_shapes <- c(16, 17, 15, 3, 7, 8)
  shape_vals <- base_shapes[(seq_along(src_levels) - 1) %% length(base_shapes) + 1]
  
  gg <- ggplot(
    pc_df %>% mutate(label_short = sub("_.*$", "", raster)), 
    aes(Axis1, Axis2, label = label_short, color = model, shape = source_simp)) +
    geom_hline(yintercept = 0, linetype = 3, color = "grey60") +
    geom_vline(xintercept = 0, linetype = 3, color = "grey60") +
    geom_point(size = 3, na.rm = TRUE, show.legend = TRUE) +
    ggrepel::geom_text_repel(size = 3, max.overlaps = 30, na.rm = TRUE, show.legend = FALSE) +
    theme_minimal(base_size = 12) +
    scale_shape_manual(values = setNames(shape_vals, src_levels), drop = FALSE) +
    labs(
      title = paste0("PCoA on '", PCOA_DISTANCE, "' distance"),
      subtitle = "Lower distances = more similar",
      color = "Model", shape = "Source"
    )
  
  ggsave(filename = file.path(OUT_DIR, paste0("PCoA_", PCOA_DISTANCE, "_model", ".png")),
         plot = gg, width = 9, height = 7, dpi = 220)
} else {
  warning("PCoA skipped: not enough rasters with finite distances after NA filtering.")
}

# ============================================================
# Correspondence Analysis (CA) on raster x class table
# ============================================================
# -----------------------
# Build raster × class (0/1/2) count table robustly
# -----------------------
# This avoids freq() argument issues and counts exactly what we need.

class_counts <- matrix(0, nrow = length(rlist), ncol = 3,
                       dimnames = list(names(rlist), c("notgrass0", "tame1", "native2")))

for (i in seq_along(rlist)) {
  r <- rlist[[i]]
  # Make sure we work with a single-band, integer-like layer
  if (terra::nlyr(r) > 1) r <- r[[1]]
  
  # Count pixels equal to the class value; NA are ignored by na.rm=TRUE
  c0 <- tryCatch(as.numeric(terra::global(r == 0, "sum", na.rm = TRUE)[1,1]), error = function(e) NA_real_)
  c1 <- tryCatch(as.numeric(terra::global(r == 1, "sum", na.rm = TRUE)[1,1]), error = function(e) NA_real_)
  c2 <- tryCatch(as.numeric(terra::global(r == 2, "sum", na.rm = TRUE)[1,1]), error = function(e) NA_real_)
  
  class_counts[i, ] <- c(c0, c1, c2)
}

# Write the counts to CSV
readr::write_csv(
  tibble::tibble(
    raster = rownames(class_counts),
    notgrass0 = class_counts[, "notgrass0"],
    tame1     = class_counts[, "tame1"],
    native2   = class_counts[, "native2"]
  ),
  file.path(OUT_DIR, "raster_by_class_counts.csv")
)

# -----------------------
# Correspondence Analysis (CA) on the counts table 
# -----------------------
if (isTRUE(RUN_CA)) {
  message("Running Correspondence Analysis (CA) on class compositions (using robust counts)...")
  
  tot <- rowSums(class_counts, na.rm = TRUE)
  cc_sub <- class_counts[tot > 0, , drop = FALSE]
  
  if (nrow(cc_sub) >= 3) {
    
    ca_fit <- ca::ca(cc_sub)
    
    # Shorten raster labels to text before the first underscore
    short_labels <- sub("_.*$", "", rownames(cc_sub))
    
    # Extract coordinates
    row_coords <- as.data.frame(ca_fit$rowcoord)   # rasters
    col_coords <- as.data.frame(ca_fit$colcoord)   # classes
    
    colnames(row_coords) <- c("CA1", "CA2")
    colnames(col_coords) <- c("CA1", "CA2")
    
    row_coords$label <- short_labels
    row_coords$type  <- "Raster"
    col_coords$label <- rownames(col_coords)
    col_coords$type  <- "Class"
    
    # Optional scaling for classes (keep at 1 for a standard biplot)
    SCALE_FACTOR <- 1.0
    col_coords$CA1 <- col_coords$CA1 * SCALE_FACTOR
    col_coords$CA2 <- col_coords$CA2 * SCALE_FACTOR
    
    # ---- Compute padded limits from BOTH rasters and classes ----
    all_x <- c(row_coords$CA1, col_coords$CA1)
    all_y <- c(row_coords$CA2, col_coords$CA2)
    
    # If labels extend beyond points, add padding in data space
    PAD <- 0.12  # 12% padding on each side; adjust to 0.15–0.20 if labels still clip
    
    rng_x <- range(all_x, na.rm = TRUE)
    rng_y <- range(all_y, na.rm = TRUE)
    pad_x <- diff(rng_x) * PAD
    pad_y <- diff(rng_y) * PAD
    xlim  <- c(rng_x[1] - pad_x, rng_x[2] + pad_x)
    ylim  <- c(rng_y[1] - pad_y, rng_y[2] + pad_y)
    
    # ---- Plot with equal aspect ratio and unclipped text ----
    png(file.path(OUT_DIR, "CA_raster_by_class.png"),
        width = 1600, height = 1200, res = 200)
    op <- par(no.readonly = TRUE)
    on.exit(par(op), add = TRUE)
    
    # margins: leave room for labels
    par(mar = c(4.5, 4.5, 3.5, 2))
    par(xpd = NA)          # allow drawing outside plot region
    par(pty = "s")         # square plotting region
    
    plot(NA, NA,
         xlim = xlim, ylim = ylim,
         xlab = "CA1", ylab = "CA2",
         main = "Correspondence Analysis (rasters × classes)",
         asp = 1)          # keep axes on same scale
    
    abline(h = 0, v = 0, col = "grey70", lty = 3)
    
    # Raster points (black) + labels
    points(row_coords$CA1, row_coords$CA2,
           pch = 16, col = "black", cex = 1.2)
    text(row_coords$CA1, row_coords$CA2,
         labels = row_coords$label,
         pos = 3, cex = 0.95, col = "black", offset = 0.35)
    
    # Class points (red triangles) + labels
    points(col_coords$CA1, col_coords$CA2,
           pch = 17, col = "red3", cex = 1.6)
    text(col_coords$CA1, col_coords$CA2,
         labels = col_coords$label,
         pos = 3, cex = 1.1, col = "red3", offset = 0.35)
    
    legend("topright",
           legend = c("Rasters", "Classes"),
           pch = c(16, 17),
           col = c("black", "red3"),
           bty = "n", cex = 1.05)
    
    dev.off()
  } else {
    warning("CA skipped: too few rasters with non-zero class counts.")
  }
}

# check what values are present in rasters
for (nm in names(rlist)) {
  print(nm)
  print(freq(rlist[[nm]], digits = 0, value=TRUE, useNA="ifany"))
}

unique_vals <- lapply(rlist, function(r) sort(unique(values(r))))
unique_vals

# ============================================================
# Optional: Union-grid coverage & consensus maps
# ============================================================
if (isTRUE(MAKE_UNION_CONSENSUS)) {
  message("Building union-grid coverage & consensus maps (nearest resampling)...")
  # Use CRS + resolution of the first raster; union extent across rasters
  r1 <- rlist[[1]]
  r1_res <- res(r1)
  r1_crs <- crs(r1, proj=TRUE)
  
  # compute union extent across all rasters
  ex <- ext(r1)
  for (i in 2:n) {
    ei <- ext(rlist[[i]])
    ex <- ext(min(xmin(ex), xmin(ei)),
              max(xmax(ex), xmax(ei)),
              min(ymin(ex), ymin(ei)),
              max(ymax(ex), ymax(ei)))
  }
  
  template <- rast(ext = ex, resolution = r1_res, crs = r1_crs)
  
  # resample each raster to template (nearest)
  r_resamp <- lapply(rlist, function(r) {
    tryCatch(resample(r, template, method = "near"), error = function(e) NULL)
  })
  # Drop NULLs if any resampling failed
  keep_idx <- which(vapply(r_resamp, function(x) !is.null(x), logical(1)))
  if (length(keep_idx) < 1) {
    warning("All resamples failed; skipping union-grid products.")
  } else {
    r_resamp <- r_resamp[keep_idx]
    names(r_resamp) <- ras_names[keep_idx]
    s <- rast(r_resamp)  # SpatRaster with layers
    
    # Coverage = number of non-NA layers per cell
    coverage <- app(s, fun = function(v) sum(!is.na(v)))
    names(coverage) <- "coverage"
    
    # Per-class counts
    count_class <- function(val) app(s, fun = function(v) sum(v == val, na.rm=TRUE))
    cnt0 <- count_class(0); names(cnt0) <- "count_class0"
    cnt1 <- count_class(1); names(cnt1) <- "count_class1"
    cnt2 <- count_class(2); names(cnt2) <- "count_class2"
    
    # Consensus class: which class has the max count; NA where coverage==0
    # (Ties broken by the first max encountered)
    consensus <- app(c(cnt0, cnt1, cnt2), fun = function(v) {
      # v length 3
      if (sum(v) == 0) return(NA)
      which.max(v) - 1  # returns 0/1/2
    })
    names(consensus) <- "consensus_class"
    
    # Consensus fraction = max(count)/coverage; NA when coverage==0
    max_count <- app(c(cnt0, cnt1, cnt2), fun = function(v) max(v, na.rm=TRUE))
    names(max_count) <- "max_count"
    consensus_fraction <- max_count / coverage
    names(consensus_fraction) <- "consensus_fraction"
    consensus_fraction[coverage == 0] <- NA
    
    # Write outputs
    writeRaster(coverage,            filename = file.path(OUT_DIR, "union_coverage_count.tif"),      overwrite = TRUE)
    writeRaster(cnt0,                filename = file.path(OUT_DIR, "union_count_class0.tif"),         overwrite = TRUE)
    writeRaster(cnt1,                filename = file.path(OUT_DIR, "union_count_class1.tif"),         overwrite = TRUE)
    writeRaster(cnt2,                filename = file.path(OUT_DIR, "union_count_class2.tif"),         overwrite = TRUE)
    writeRaster(consensus,           filename = file.path(OUT_DIR, "union_consensus_class.tif"),      overwrite = TRUE)
    writeRaster(consensus_fraction,  filename = file.path(OUT_DIR, "union_consensus_fraction.tif"),   overwrite = TRUE)
  }
}

# ============================================================
# Optional: Pairwise disagreement map (0=match, 1=mismatch) on overlap
# ============================================================
make_pair_disagreement <- function(nameA, nameB) {
  # Try match by 'raster' (as in meta2$raster)
  a_idx <- which(ras_names %in% c(nameA, tools::file_path_sans_ext(nameA)))
  b_idx <- which(ras_names %in% c(nameB, tools::file_path_sans_ext(nameB)))
  if (length(a_idx) != 1 || length(b_idx) != 1) {
    # Fallback: try match by basename of file
    base_map <- setNames(seq_along(rlist), tools::file_path_sans_ext(basename(meta2$file)))
    if (length(a_idx) != 1) a_idx <- base_map[[tools::file_path_sans_ext(nameA)]]
    if (length(b_idx) != 1) b_idx <- base_map[[tools::file_path_sans_ext(nameB)]]
  }
  if (length(a_idx) != 1 || length(b_idx) != 1 || is.na(a_idx) || is.na(b_idx)) {
    warning("Could not resolve pair indices for disagreement map: ", nameA, " vs ", nameB)
    return(invisible(NULL))
  }
  rA <- rlist[[a_idx]]
  rB <- rlist[[b_idx]]
  e <- intersect(ext(rA), ext(rB))
  if (is.null(e)) {
    warning("No overlap for pair map: ", nameA, " vs ", nameB)
    return(invisible(NULL))
  }
  a <- crop(rA, e, snap="near")
  b <- crop(rB, e, snap="near")
  if (!compareGeom(a,b, stopOnError=FALSE)) {
    b <- resample(b, a, method="near")
  }
  # mismatch: NA where either NA, else 0 if equal, 1 if different
  # Use classify/ifel to stay lazy
  BOTH_OK <- !is.na(a) & !is.na(b)
  mismatch <- ifel(BOTH_OK, ifel(a == b, 0, 1), NA)
  names(mismatch) <- paste0("mismatch_", ras_names[a_idx], "_vs_", ras_names[b_idx])
  out_tif <- file.path(OUT_DIR, paste0(names(mismatch), ".tif"))
  writeRaster(mismatch, filename = out_tif, overwrite = TRUE)
  
  # quick plot summary as PNG
  vals <- values(mismatch)
  tab <- table(factor(vals, levels=c(0,1), labels=c("match","mismatch")), useNA="no")
  dfp <- as.data.frame(tab)
  png(file.path(OUT_DIR, paste0(names(mismatch), "_bar.png")), width=1000, height=800, res=150)
  print(
    ggplot(dfp, aes(x = Var1, y = Freq, fill = Var1)) +
      geom_col() +
      scale_fill_manual(values=c("#4daf4a","#e41a1c")) +
      labs(x=NULL, y="Pixel count", title = names(mismatch)) +
      theme_minimal(base_size = 12) +
      theme(legend.position = "none")
  )
  dev.off()
  invisible(out_tif)
}

if (!is.null(PAIR_MAP_A) && !is.null(PAIR_MAP_B)) {
  message("Producing disagreement map for pair: ", PAIR_MAP_A, " vs ", PAIR_MAP_B)
  try(make_pair_disagreement(PAIR_MAP_A, PAIR_MAP_B), silent = TRUE)
}

message("All done. Outputs in: ", OUT_DIR)
