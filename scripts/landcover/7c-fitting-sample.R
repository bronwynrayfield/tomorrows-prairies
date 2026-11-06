# Purpose:
#   Draws the fitting sample for every conversion model, once per 2010 origin
#   class, so the RF (8b) and GLM benchmark (8a) fit on identical rows and
#   ecodistrict folds. Model scripts only read these files and subset them.
#
#   One frame per origin, built from the multiclass response (7-response.R):
#     native : persistence + annual, settlement, tame positives
#     tame   : persistence + annual, settlement positives
#     annual : persistence + tame positives
#   Each binary model takes persistence plus one destination (option A: zeros
#   are persistence only; competing destinations are absent, not 0). The same
#   frame is the multiclass training set. Shared persistence rows make
#   coefficients and AUCs comparable across transitions from one origin.
#   Outcomes not listed (e.g. 274 ha tame -> other) are excluded.
#
#   Training population (unchanged from the RF script):
#     (1) 2010 class == origin
#     (2) outside CPCAD 2010 (protected pixels are structural zeros:
#         native->annual 0.0100% inside vs 0.276% outside)
#     (3) predictor-valid mask (4f)
#
#   Every class, persistence included, is capped at n_class, so each binary
#   frame is at most 1:1. A destination with fewer than n_class positives in
#   the population is taken whole, and its binary is imbalanced towards
#   persistence; model scripts downsample persistence if they need 1:1.
#
#   Two stages, two files:
#     train-points_<origin>.rds : x, y, lulc_2020 (factor), ecodistrict.
#                                 Drawn once; delete to redraw.
#     train-sample_<origin>.rds : points + predictors (transformed per 0b
#                                 predictor_transforms) + w. Re-extracted whenever
#                                 the points, 0b, the retained set, any predictor
#                                 raster or the 7a weights are newer than it, so a
#                                 predictor change never moves rows or folds.
#   Columns and attributes:
#     w         : case weight. Tame and annual persistence rows carry the 7a ecodistrict
#                 weight (negative_weights.csv), rescaled to mean 1; all others 1.
#     pop_n     : training-population pixel count per 2020 class (uncertain tame and annual
#                 pixels excluded); the per-class sampling fractions the
#                 multiclass intercept correction needs.
#     base_rate : named per destination, dest / (dest + persistence), from the
#                 hindcast map (uncertain tame and annual pixels count as persistence), so the
#                 prior correction targets the map rather than the training response.
#   Attributes are dropped by row subsetting: read them before subsetting.

source(file.path("scripts/landcover", "0a-setup.R"))
source(file.path("scripts/landcover", "0b-predictor-files.R"))

# First element of each is persistence; names become the lulc_2020 labels.
classes <- list(native = cls[c("native", "annual", "settlement", "tame")],
                tame   = cls[c("tame", "annual", "settlement")],
                annual = cls[c("annual", "tame", "settlement")])

n_class <- 150000   # per class. Native->annual has ~188k unprotected valid
# positives, the binding class; raising this past it
# leaves that binary short of 1:1.

# ---- inputs ----
sa_r <- rast(file.path(sa_filepath, sa_raster_filename))

in_sa         <- file.path("./processeddatasets", "predictors", "30m")
out_dir       <- "output/landcover/model"
retained_file <- file.path(in_sa, "retained-covariates.rds")
mask_file     <- file.path(in_sa, "predictor-valid-mask_study-area.tif")
resp_file     <- tagged_path("./processeddatasets", "response",
                             "response_multiclass.tif", "study-area")
lulc_file     <- tagged_path("./processeddatasets", "landcover", "lulc_2010.tif", "study-area")
pa_file       <- file.path("processeddatasets", "management", "Protected Areas Combined", "CPCAD_pa_2010.tif")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- predictors ----
final_set <- readRDS(retained_file)
stopifnot(all(final_set %in% names(predictor_files)))

rs <- rast(lapply(predictor_files[final_set], \(f) rast(f, win = ext(sa_r))))
names(rs) <- final_set

# newest upstream input that should trigger re-extraction
pred_mtime <- max(file.mtime(c(file.path("scripts/landcover", "0b-predictor-files.R"),
                               retained_file, predictor_files[final_set],
                               file.path(response_summary_filepath, "negative_weights.csv"))))
# ---- shared rasters ----
# Categories dropped so comparisons and freq() work on the integer codes.
resp <- crop(rast(resp_file), sa_r, snap = "out"); levels(resp) <- NULL
lulc <- rast(lulc_file); levels(lulc) <- NULL
msk  <- mask(crop(rast(mask_file), sa_r, snap = "out"), sa_r)
pa   <- mask(crop(rast(pa_file),   sa_r, snap = "out"), sa_r)

ecodistricts_sa <- vect(st_read(tagged_path("./processeddatasets", ecodistricts_filepath, ecodistricts_filename, "studyarea")))
ecodistrict_r  <- crop(rast(tagged_path("./processeddatasets", ecodistricts_filepath, paste0(tools::file_path_sans_ext(ecodistricts_filename), ".tif"), "studyarea")) , lulc)

map_file <- tagged_path("./processeddatasets", "landcover", "lulc_2020.tif", "study-area")
map20 <- crop(rast(map_file), sa_r, snap = "out"); levels(map20) <- NULL

# Positives are enumerated: native destinations are ~0.03% of valued cells,
# too sparse for rejection. Persistence is the majority, so rejection is cheap.
draw_all <- function(m, n) {
  i <- cells(m)
  xyFromCell(m, if (length(i) > n) sample(i, n) else i)
}
draw_reject <- function(m, n) {
  frac <- as.numeric(global(!is.na(m), "mean", na.rm = TRUE))
  p <- spatSample(m, ceiling(n * 3 / frac), method = "random", xy = TRUE,
                  values = TRUE, na.rm = FALSE)
  p <- as.matrix(p[!is.na(p[[3]]), c("x", "y")])
  p[sample(nrow(p), min(n, nrow(p))), , drop = FALSE]
}

for (origin in names(classes)) {
  cl <- classes[[origin]]
  points_file <- file.path(out_dir, paste0("train-points_", origin, ".rds"))
  sample_file <- file.path(out_dir, paste0("train-sample_", origin, ".rds"))
  cat("\n==== ", origin, " ====\n", sep = "")
  
  # ---- stage 1: points & response ----
  if (file.exists(points_file)) {
    pts <- readRDS(points_file)
  } else {
    resp_o <- mask(resp, subst(msk & pa == 0 & lulc == cl[[1]] & (resp %in% cl),
                               FALSE, NA))
    
    f_o   <- freq(resp_o)
    pop_n <- setNames(f_o$count[match(cl, f_o$value)], names(cl))
    stopifnot(!anyNA(pop_n))
    
    # Base rate from the hindcast map, not the training response: uncertain tame pixels
    # are NA in resp but persistence on the map, and the prior correction targets the map.
    f_map   <- freq(mask(map20, subst(msk & pa == 0 & lulc == cl[[1]] & (map20 %in% cl), FALSE, NA)))
    pop_map <- setNames(f_map$count[match(cl, f_map$value)], names(cl))
    base_rate_map <- pop_map[-1] / (pop_map[-1] + pop_map[[1]])
    
    set.seed(42)
    xy <- lapply(seq_along(cl), function(k) {
      m <- subst(resp_o == cl[[k]], FALSE, NA)
      p <- if (k == 1) draw_reject(m, n_class) else draw_all(m, n_class)
      data.frame(p, lulc_2020 = names(cl)[k])
    })
    pts <- do.call(rbind, xy)
    pts$lulc_2020 <- factor(pts$lulc_2020, levels = names(cl))
    pts$ecodistrict <- factor(terra::extract(ecodistrict_r, pts[, c("x", "y")])[, 2])
    pts <- pts[!is.na(pts$ecodistrict), ]
    pts$ecodistrict <- droplevels(pts$ecodistrict)
    
    attr(pts, "pop_n")     <- pop_n
    attr(pts, "base_rate") <- base_rate_map
    saveRDS(pts, points_file)
    rm(resp_o, xy)
  }
  
  # ---- stage 2: predictors ----
  if (!file.exists(sample_file) ||
      file.mtime(sample_file) < max(pred_mtime, file.mtime(points_file))) {
    d <- cbind(pts, apply_transforms(terra::extract(rs, pts[, c("x", "y")], ID = FALSE)))
    n_na <- sum(!complete.cases(d))
    d <- na.omit(d)
    
    # Tame and annual persistence rows carry the 7a ecodistrict weight (uncertain pixels are NA in
    # training but their 2010 class on the map, and concentrate where conversion is high); all other
    # rows 1. Rescaled to mean 1 over persistence rows so the 1:1 balance holds.
    d$w <- 1
    if (origin %in% c("tame", "annual")) {
      zw  <- read.csv(file.path(response_summary_filepath, "negative_weights.csv")) |> filter(origin == !!origin)
      neg <- d$lulc_2020 == origin
      d$w[neg] <- zw$w_neg[match(as.character(d$ecodistrict[neg]), as.character(zw$ecodistrict))]
      d$w[neg] <- d$w[neg] / mean(d$w[neg])
    }
    stopifnot(!anyNA(d$w))
    
    attr(d, "pop_n")     <- attr(pts, "pop_n")
    attr(d, "base_rate") <- attr(pts, "base_rate")
    saveRDS(d, sample_file)
    cat("extracted", length(final_set), "predictors;", n_na, "rows dropped for NA\n")
  } else {
    d <- readRDS(sample_file)
    cat("sample up to date\n")
  }
  
  cat("population:\n"); print(attr(d, "pop_n"))
  cat("base rate (%):\n"); print(round(100 * attr(d, "base_rate"), 4))
  cat("sampled:\n"); print(table(d$lulc_2020))
  cat("ecodistricts:", nlevels(d$ecodistrict), "\n")
  rm(pts, d)
}

# 
# ==== native ====
#   extracted 31 predictors; 0 rows dropped for NA
# population:
#   native     annual settlement       tame 
# 62153726     185947     228299     342095 
# base rate (%):
#   annual settlement       tame 
# 0.2983     0.3660     0.5474 
# sampled:
#   
#   native     annual settlement       tame 
# 149992     149997     149980     149980 
# ecodistricts: 126 
# 
# ==== tame ====
#   extracted 31 predictors; 0 rows dropped for NA
# population:
#   tame     annual settlement 
# 45962071    9627153     365125 
# base rate (%):
#   annual settlement 
# 14.3479     0.6313 
# sampled:
#   
#   tame     annual settlement 
# 149973     149972     149959 
# ecodistricts: 127 
#
# ==== annual ====
extracted 31 predictors; 0 rows dropped for NA
>   
  >   cat("population:\n"); print(attr(d, "pop_n"))
population:
  annual       tame settlement 
214188766    3777311     486179 
>   cat("base rate (%):\n"); print(round(100 * attr(d, "base_rate"), 4))
base rate (%):
  tame settlement 
1.5829     0.2066 
>   cat("sampled:\n"); print(table(d$lulc_2020))
sampled:
  
  annual       tame settlement 
149993     149968     149974 
>   cat("ecodistricts:", nlevels(d$ecodistrict), "\n")
ecodistricts: 127 
>   rm(pts, d)