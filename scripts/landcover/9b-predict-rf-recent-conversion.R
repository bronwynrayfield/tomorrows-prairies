# Purpose:
#   Applies the per-stratum conversion random forests (native, tame) to produce
#   P(grassland -> annual cropland over a 10-year window) at 30 m.
#
#   pred_period = 2010: predictors are the rasters the models were fitted on,
#     so the surface is an in-sample hindcast of 2010 -> 2020 conversion. Use it
#     to check spatial sense and calibration against the observed response.
#   pred_period = 2025: forward 10-year risk, assuming the 2010-2020
#     relationships hold. Requires period-specific predictor paths in
#     0b-predictor-files.R (Canfax net returns z-scored within 2020-2024;
#     climate from SSP2-4.5 rather than the 1991-2020 baseline).
#
#   Each stratum is predicted only over its own cells in the starting land
#   cover (lulc_<period>), not over the response's non-NA cells, so the same
#   masking works in 2025 when no response exists, and cells NA'd by a
#   response variant (e.g. excl_per native -> perennial) still get a value.
#   The two strata don't overlap, so they're combined with cover() for CLUES.
#
#   NB the ~6.8% of the study area lost to DSS NAs along SK waterways is blank.
#
# Outputs:
#   - pred_conversion_<stratum>_<period>_<res>m.tif : per-stratum probability
#   - pred_conversion_<period>_<res>m.tif           : strata combined, for CLUES
#   - pred_conversion_<period>_300m.tif             : display copy (30 m run)
# =============================================================================

source(file.path("scripts/landcover", "0a-setup.R"))
source(file.path("scripts/landcover", "0b-predictor-files.R"))

terraOptions(memfrac = 0.5)

# ---- inputs ----
# study area
sa_r <- rast(file.path(sa_filepath, sa_raster_filename))

pred_period <- 2010      # 2010 = fitting rasters; 2025 once 0b has its paths
pred_res    <- 300 #res_m     # set to 300 for a minutes-long preview first (not
# equivalent to the 30 m run: the forest is non-linear)
out_res     <- 300       # display resolution for the native run

# stratum codes in the lulc layers, from the shared cls_* definition in 0a-setup.R
strata <- c(native = cls_native, tame = cls_tame)

in_sa         <- file.path("./processeddatasets", "predictors", "30m")
model_dir     <- "output/landcover/model"
out_pred      <- file.path("./processeddatasets", "predictions", "recent-conversion")
retained_file <- file.path(in_sa, "retained-covariates.rds")
mask_file     <- file.path(in_sa, "predictor-valid-mask_study-area.tif")
lulc_file     <- tagged_path("./processeddatasets", "landcover", "lulc_2010.tif", "30m")
resp_file     <- file.path("./processeddatasets", "response",
                           "response_modelling-extent.tif")               
model_file    <- function(s) file.path(model_dir, sprintf("rf_conversion_%s.rds", s))

pf <- if (pred_period == 2010) predictor_files else
  stop("add period-specific predictor paths to 0b-predictor-files.R")

stopifnot(pred_res %% res_m == 0, out_res %% res_m == 0)
pred_fact <- pred_res / res_m
out_fact  <- out_res  / res_m

dir.create(out_pred, recursive = TRUE, showWarnings = FALSE)

gopt <- c("COMPRESS=DEFLATE", "ZLEVEL=9", "TILED=YES",
          "BLOCKXSIZE=256", "BLOCKYSIZE=256",
          "BIGTIFF=YES", "NUM_THREADS=ALL_CPUS", "PREDICTOR=3")

final_set <- readRDS(retained_file)
# temporary (delete later when rerun collinearity script)
dist_map <- c("d_elevator" = "dist_to_elevator",
              "d_road"     = "dist_to_road",
              "d_rail"     = "dist_to_rail",
              "d_popctr"   = "dist_to_popctr",
              "d_regctr"   = "dist_to_regctr")
i <- match(final_set, names(dist_map))
stopifnot(sum(!is.na(i)) == length(dist_map))
final_set[!is.na(i)] <- dist_map[i[!is.na(i)]]
stopifnot(all(final_set %in% names(predictor_files)))
# end of temp
stopifnot(length(final_set) == 31)   # conversion set; plowprint benchmark used 23

prior <- data.frame(stratum = s,
                    pi1 = <positives / cells in stratum, from the response within the valid mask>,
                    s1  = mean(d$y == 1))   # prevalence the trees were trained on
saveRDS(prior, file.path(model_dir, sprintf("prior_conversion_%s.rds", s)))


# ---- predictor stack ----
# Aggregate before cropping so a coarse run never materialises a 31-layer copy
# at 30 m. The stack is never masked or written; masking is done per stratum
# via an extra layer (below).
msk  <- rast(mask_file)
lulc <- crop(rast(lulc_file), msk, snap = "out")

rs <- rast(lapply(predictor_files[final_set], \(f) rast(f, win = ext(sa_r))))
names(rs) <- final_set
if (pred_res > res_m) rs <- aggregate(rs, fact = pred_fact, fun = "mean", na.rm = TRUE)
rs <- crop(rs, msk, snap = "out")

# ---- predict ----
# The stratum mask goes into the stack as an extra layer, NA outside the
# stratum. With na.rm = TRUE terra drops any row with an NA, so only the
# stratum's cells (88.9M native, 74.7M tame of the study area, fewer after the
# valid mask) reach ranger, instead of every valid cell in the study area,
# cropland included. ranger selects its own variables by name, so the extra
# column is ignored.

# Need to add this to fit script
# for (s in names(strata)) {
#   smsk <- subst(msk & lulc == strata[[s]], FALSE, NA)
#   pi1  <- global(mask(resp, smsk), "mean", na.rm = TRUE)[[1]]
#   d <- readRDS(file.path(model_dir, paste0("train-sample_conversion_", s, ".rds")))
#   s1 <- mean(d$converted == 1)
#   saveRDS(data.frame(stratum = s, pi1 = pi1, s1 = s1),
#           file.path(model_dir, sprintf("prior_conversion_%s.rds", s)))
# }
# end

pfun <- function(model, ..., w1, w0) {
  p <- predict(model, ..., num.threads = parallel::detectCores())$predictions[, "1"]
  p * w1 / (p * w1 + (1 - p) * w0)
}

pred_path <- function(s) file.path(out_pred, sprintf("pred_conversion_%s%d_%dm%s.tif",
                                                     if (is.null(s)) "" else paste0(s, "_"), pred_period, pred_res,
                                                     if (pred_res > res_m) "_preview" else ""))

# prior correction for balanced training (case-control): rescales odds by
# (pi1/s1) / ((1-pi1)/(1-s1)); exact only if sampling depended on class alone

fits   <- list()
priors <- list()
preds  <- list()
for (s in names(strata)) {
  fits[[s]] <- readRDS(model_file(s))
  stopifnot(setequal(fits[[s]]$forest$independent.variable.names, final_set))
  
  # pi1 = 2010-2020 population prevalence in the stratum, s1 = training
  # prevalence; both saved by the fit script so 2025 reuses the 2010 pi1
  priors[[s]] <- readRDS(file.path(model_dir, sprintf("prior_conversion_%s.rds", s)))
  w1 <- priors[[s]]$pi1 / priors[[s]]$s1
  w0 <- (1 - priors[[s]]$pi1) / (1 - priors[[s]]$s1)
  
  smsk <- subst(msk & lulc == strata[[s]], FALSE, NA)
  if (pred_res > res_m) smsk <- aggregate(smsk, fact = pred_fact, fun = "max", na.rm = TRUE)
  smsk <- resample(crop(smsk, rs, snap = "out"), rs, method = "near")
  names(smsk) <- "stratum_mask"
  
  t0 <- Sys.time()
  preds[[s]] <- predict(c(rs, smsk), fits[[s]], fun = pfun, na.rm = TRUE,
                        w1 = w1, w0 = w0,
                        filename = pred_path(s), overwrite = TRUE,
                        wopt = list(gdal = gopt, datatype = "FLT4S"))
  cat(s, "predicted in", format(Sys.time() - t0), "\n")
}



# fits  <- list()
# preds <- list()
# for (s in names(strata)) {
#   fits[[s]] <- readRDS(model_file(s))
#   stopifnot(setequal(fits[[s]]$forest$independent.variable.names, final_set))
#   
#   smsk <- subst(msk & lulc == strata[[s]], FALSE, NA)
#   if (pred_res > res_m) smsk <- aggregate(smsk, fact = pred_fact, fun = "max", na.rm = TRUE)
#   smsk <- resample(crop(smsk, rs, snap = "out"), rs, method = "near")
#   names(smsk) <- "stratum_mask"
#   
#   t0 <- Sys.time()
#   preds[[s]] <- predict(c(rs, smsk), fits[[s]], fun = pfun, na.rm = TRUE,
#                         filename = pred_path(s), overwrite = TRUE,
#                         wopt = list(gdal = gopt, datatype = "FLT4S"))
#   cat(s, "predicted in", format(Sys.time() - t0), "\n")
# }

pred <- cover(preds$native, preds$tame,
              filename = pred_path(NULL), overwrite = TRUE,
              wopt = list(gdal = gopt, datatype = "FLT4S"))

if (pred_res == res_m) {
  aggregate(pred, fact = out_fact, fun = "mean", na.rm = TRUE,
            filename = file.path(out_pred, sprintf("pred_conversion_%d_%dm.tif",
                                                   pred_period, out_res)),
            overwrite = TRUE, wopt = list(gdal = gopt, datatype = "FLT4S"))
}

# ---- QA ----
# Predicted mean should match the observed within-stratum rate (5-response.R
# baseline: native 0.0022, tame 0.174) ONLY if the training sample was drawn
# proportionally within stratum. If native was downsampled or case-weighted,
# expect the predicted mean to sit above the observed rate; the surface then
# needs a prior correction before it goes to CLUES.
# Accuracy at a 0.5 cutoff is omitted: at 0.22% positives it just reproduces
# the no-skill floor. Discrimination is judged by the LOEO AUC from the fit script.
cat("period / resolution:", pred_period, "/", pred_res, "m\n")
for (s in names(strata)) {
  pm <- global(preds[[s]], "mean", na.rm = TRUE)[[1]]
  cat(sprintf("%-6s predicted mean %.4f | OOB Brier %.4f", s, pm,
              fits[[s]]$prediction.error))
  if (pred_period == 2010 && pred_res == res_m) {
    resp <- crop(rast(resp_file), preds[[s]], snap = "out")
    om   <- global(mask(resp, preds[[s]]), "mean", na.rm = TRUE)[[1]]
    cat(sprintf(" | observed rate %.4f", om))
  }
  cat("\n")
}
plot(pred)

# For 2025: before trusting the surface, compare each predictor's range with
# its training range. The forest flatlines beyond the training envelope, so
# cells past it (Canfax: 62% raw -> 2.4% after z-scoring) get edge values.