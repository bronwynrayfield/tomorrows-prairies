# Purpose:
#   Applies the per-transition conversion RFs (8-rf-conversion.R) to produce
#   P(origin -> destination over a 10-year window), one surface per transition,
#   then mosaics the two origins into one surface per destination for CLUES.
#
#   pred_period = 2010: predictors are the rasters the models were fitted on,
#     so the surface is an in-sample hindcast of 2010 -> 2020 conversion. Use it
#     to check spatial sense and calibration, and to compare RF vs GLM.
#   pred_period = 2025: forward 10-year risk, assuming the 2010-2020
#     relationships hold. Requires period-specific predictor paths in
#     0b-predictor-files.R (Canfax net returns z-scored within 2020-2024;
#     climate from SSP2-4.5 rather than the 1991-2020 baseline).
#
#   Domain: every origin cell in the predictor-valid mask of the starting land
#   cover (lulc_<period>), protected cells included, so the same masking works
#   in 2025 when no response exists. The models were trained on unprotected
#   cells only, so protected values are what the model expects of that land if
#   it were unprotected; CLUES excludes protected land via its restriction
#   layer. QA is computed on the training population (unprotected), where the
#   observed rate is the base rate.
#
#   Prior correction: training is balanced, so odds are rescaled by
#   (br / s1) / ((1 - br) / (1 - s1)), with br (training-population base rate)
#   and s1 (sample prevalence) stored on each model by 8. Exact only if
#   sampling depended on class alone, which 7c's per-class draw satisfies.
#
# Outputs (in out_pred; _preview suffix when pred_res > res_m):
#   - pred_rf_conversion_<from>-<to>_<period>_<res>m.tif : per-transition probability
#   - pred_rf_<to>_<period>_<res>m.tif                    : origins combined, for CLUES
#                                                          (single-origin destinations:
#                                                          use the transition surface)
#   - pred_rf_<to>_<period>_<out_res>m.tif                : display copy (30 m run)

source(file.path("scripts/landcover", "0a-setup.R"))
source(file.path("scripts/landcover", "0b-predictor-files.R"))

terraOptions(memfrac = 0.5)

transitions <- data.frame(from = c("native", "tame",   "native",     "tame",       "native"),
                          to   = c("annual", "annual", "settlement", "settlement", "tame"))

# ---- inputs ----
sa_r <- rast(file.path(sa_filepath, sa_raster_filename))

pred_period <- 2010      # 2010 = fitting rasters; 2025 once 0b has its paths
pred_res    <- res_m       # 300 for a minutes-long preview first, res_m for the
# real run (not equivalent: the forest is non-linear)
out_res     <- 300       # display resolution for the 30 m run

in_sa         <- file.path("./processeddatasets", "predictors", "30m")
model_dir     <- "output/landcover/model"
out_pred      <- file.path("./processeddatasets", "predictions", "recent-conversion")
retained_file <- file.path(in_sa, "retained-covariates.rds")
mask_file     <- file.path(in_sa, "predictor-valid-mask_study-area.tif")
lulc_file     <- tagged_path("./processeddatasets", "landcover",
                             sprintf("lulc_%d.tif", pred_period), "study-area")
pa_file       <- file.path("processeddatasets", "management", "Protected Areas Combined", "CPCAD_pa_2010.tif")
resp_file     <- tagged_path("./processeddatasets", "response",
                             "response_multiclass.tif", "study-area")

pf <- if (pred_period == 2010) predictor_files else
  stop("add period-specific predictor paths to 0b-predictor-files.R")

stopifnot(pred_res %% res_m == 0, out_res %% res_m == 0)
pred_fact <- pred_res / res_m
out_fact  <- out_res  / res_m

dir.create(out_pred, recursive = TRUE, showWarnings = FALSE)

gopt <- c("COMPRESS=DEFLATE", "ZLEVEL=9", "TILED=YES",
          "BLOCKXSIZE=256", "BLOCKYSIZE=256",
          "BIGTIFF=YES", "NUM_THREADS=ALL_CPUS", "PREDICTOR=3")
wopt <- list(gdal = gopt, datatype = "FLT4S")

pred_path <- function(name, res = pred_res, ext = ".tif")
  file.path(out_pred, sprintf("pred_rf_%s_%d_%dm%s%s", name, pred_period, res,
                              if (pred_res > res_m) "_preview" else "", ext))

final_set <- readRDS(retained_file)
stopifnot(all(final_set %in% names(pf)))

# ---- predictor stack ----
# Aggregate before cropping so a coarse run never materialises a 31-layer copy
# at 30 m. The stack is never masked or written; masking is done per origin
# via an extra layer (below).
msk  <- rast(mask_file)
lulc <- crop(rast(lulc_file), msk, snap = "out")

rs <- rast(lapply(pf[final_set], \(f) rast(f, win = ext(msk))))
names(rs) <- final_set
if (pred_res > res_m) rs <- crop(aggregate(rs, fact = pred_fact, fun = "mean", na.rm = TRUE),
                                 msk, snap = "out")

# ---- origin masks ----
# One per origin, reused across its destinations. Added to the stack as an
# extra layer, NA outside the origin: with na.rm = TRUE terra drops any row
# with an NA, so only the origin's cells reach ranger. ranger selects its
# variables by name, so the extra column is ignored.
omsk <- lapply(c(native = cls_native, tame = cls_tame), function(v) {
  m <- subst(msk & lulc == v, FALSE, NA)
  if (pred_res > res_m) {
    m <- aggregate(m, fact = pred_fact, fun = "max", na.rm = TRUE)
    m <- resample(crop(m, rs, snap = "out"), rs, method = "near")
  }
  names(m) <- "origin_mask"
  m
})

pfun <- function(model, ..., w1, w0) {
  p <- predict(model, ..., num.threads = parallel::detectCores())$predictions[, "1"]
  p * w1 / (p * w1 + (1 - p) * w0)
}

# ---- predict ----
preds <- list()
priors <- list()
for (j in seq_len(nrow(transitions))) {
  from <- transitions$from[j]
  to   <- transitions$to[j]
  tag  <- sprintf("conversion_%s-%s", from, to)
  
  fit <- readRDS(file.path(model_dir, paste0("rf_", tag, ".rds")))
  stopifnot(setequal(fit$forest$independent.variable.names, final_set))
  priors[[tag]] <- c(br = fit$base_rate, s1 = fit$sample_prev)
  
  t0 <- Sys.time()
  preds[[tag]] <- predict(c(rs, omsk[[from]]), fit, fun = pfun, na.rm = TRUE,
                          w1 = fit$base_rate / fit$sample_prev,
                          w0 = (1 - fit$base_rate) / (1 - fit$sample_prev),
                          filename = pred_path(tag), overwrite = TRUE, wopt = wopt)
  cat(tag, "predicted in", format(Sys.time() - t0), "\n")
  write.table(data.frame(tag = tag, res = pred_res, minutes = as.numeric(Sys.time() - t0, units = "mins")),
              pred_path("timing", ext = ".csv"), sep = ",", row.names = FALSE,
              col.names = !file.exists(pred_path("timing", ext = ".csv")),
              append = file.exists(pred_path("timing", ext = ".csv")))
  rm(fit); gc()
}

# conversion_native-annual predicted in 1.740213 mins 
# conversion_tame-annual predicted in 4.198872 mins 

# ---- one surface per destination ----
# Origins are disjoint, so cover() is a mosaic, not an overlay.
for (dest in unique(transitions$to)) {
  p <- preds[transitions$to == dest]
  dest_r <- if (length(p) == 1) p[[1]] else
    cover(p[[1]], p[[2]], filename = pred_path(dest), overwrite = TRUE, wopt = wopt)
  if (pred_res == res_m)
    aggregate(dest_r, fact = out_fact, fun = "mean", na.rm = TRUE,
              filename = pred_path(dest, out_res), overwrite = TRUE, wopt = wopt)
}

# ---- QA ----
# A logistic model fitted by ML reproduces the sample prevalence exactly, so
# after the intercept shift the training-population mean should sit close to
# br at 30 m. It is not exact: the correction is non-linear in p, and the
# sample covers the population's predictor space unevenly.
if (pred_period == 2010) {
  resp <- crop(rast(resp_file), msk, snap = "out"); levels(resp) <- NULL   # integer codes for %in%
  pa   <- crop(rast(pa_file),   msk, snap = "out")
  # At 30 m resp and pa are already on the msk grid that rs is windowed to,
  # so resampling would only re-process every cell.
  if (pred_res > res_m) {
    resp <- aggregate(resp, fact = pred_fact, fun = "modal", na.rm = TRUE)
    pa   <- aggregate(pa,   fact = pred_fact, fun = "max",   na.rm = TRUE)
    resp <- resample(crop(resp, rs, snap = "out"), rs, method = "near")
    pa   <- resample(crop(pa,   rs, snap = "out"), rs, method = "near")
  }
  
  cat("\nperiod / resolution:", pred_period, "/", pred_res, "m\n")
  qa_file <- pred_path("qa", ext = ".csv")
  if (file.exists(qa_file)) file.remove(qa_file)
  for (j in seq_len(nrow(transitions))) {
    from <- transitions$from[j]
    to   <- transitions$to[j]
    tag  <- sprintf("conversion_%s-%s", from, to)
    tp   <- subst(pa == 0 & resp %in% cls[c(from, to)], FALSE, NA)
    obs  <- global(mask(resp == cls[[to]], mask(tp, preds[[tag]])), "mean", na.rm = TRUE)[[1]]
    pm_t <- global(mask(preds[[tag]], tp), "mean", na.rm = TRUE)[[1]]
    pm_a <- global(preds[[tag]], "mean", na.rm = TRUE)[[1]]
    cat(sprintf("%-28s predicted mean %.4f (training pop.) vs observed %.4f | %.4f whole domain\n",
                tag, pm_t, obs, pm_a))
    write.table(data.frame(tag = tag, res = pred_res, pred_mean_tp = pm_t, observed = obs,
                           base_rate = priors[[tag]][["br"]], pred_mean_all = pm_a),
                qa_file, sep = ",", row.names = FALSE,
                col.names = !file.exists(qa_file), append = file.exists(qa_file))
  }
}

# period / resolution: 2010 / 30 m
# conversion_native-annual     predicted mean 0.1015 (training pop.) vs base rate 0.0030 | 0.0982 whole domain
# conversion_tame-annual       predicted mean 0.1546 (training pop.) vs base rate 0.1968 | 0.1538 whole domain

# For 2025: before trusting the surface, compare each predictor's range with
# its training range. The forest flatlines beyond the training envelope, so
# cells past it (Canfax: 62% raw -> 2.4% after z-scoring) get edge values.