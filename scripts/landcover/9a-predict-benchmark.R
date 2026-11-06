# Purpose:
#   Projects the GLM benchmark (8a-glm-benchmark.R) for every transition, with
#   the same period, resolution, domain, origin masks and prior correction as
#   9b-predict-rf.R, so the two sets of surfaces can be compared cell for cell.
#   Kept separate from 9 so whichever model is dropped takes its script with it.
#
#   The model is rebuilt from glm_<tag>_coef.csv: coefficients are on the
#   standardised scale, so they are converted to raw units with the stored
#   centre and scale (beta_raw = b / scale; intercept absorbs -sum(b * centre /
#   scale)). No raster is standardised.
#
#   Prior correction: for a logistic model, balanced sampling shifts only the
#   intercept, by logit(s1) - logit(br). Adding logit(br) - logit(s1) back is
#   the exact counterpart of the odds rescaling in 9, with br and s1 from
#   glm_benchmark_summary.csv.
#
#   pred_period, domain and QA: see 9b-predict-rf.R. Unlike the forest, the GLM
#   extrapolates linearly on the logit scale beyond the training envelope.
#
# Outputs (in out_pred; _preview suffix when pred_res > res_m):
#   - pred_glm_conversion_<from>-<to>_<period>_<res>m.tif : per-transition probability
#   - pred_glm_<to>_<period>_<res>m.tif                    : origins combined
#   - pred_glm_<to>_<period>_<out_res>m.tif                : display copy (30 m run)

source(file.path("scripts/landcover", "0a-setup.R"))
source(file.path("scripts/landcover", "0b-predictor-files.R"))

terraOptions(memfrac = 0.5)

# Each transition: persistence + one destination from the 7c per-origin sample (option A).
transitions <- data.frame(from = c("native", "tame",   "native",     "tame",       "native", "annual", "annual"),
                          to   = c("annual", "annual", "settlement", "settlement", "tame",   "tame", "settlement"))

# ---- inputs ----
sa_r <- rast(file.path(sa_filepath, sa_raster_filename))

pred_period <- 2010
pred_res    <- res_m       # 300 for a preview, res_m for the real run
out_res     <- 300

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

# BR reminder to myself for later
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
  file.path(out_pred, sprintf("pred_glm_%s_%d_%dm%s%s", name, pred_period, res,
                              if (pred_res > res_m) "_preview" else "", ext))

final_set <- readRDS(retained_file)
stopifnot(all(final_set %in% names(pf)))

glm_sum <- read.csv(file.path(model_dir, "glm_benchmark_summary.csv"))

# ---- predictor stack and origin masks ----
msk  <- rast(mask_file)
lulc <- crop(rast(lulc_file), msk, snap = "out")

# win = reads a window without copying. A crop() on the 30 m stack would write
# a temporary 31-layer copy of the full extent; at 300 m the stack is small,
# so the crop is cheap there.
rs <- rast(lapply(pf[final_set], \(f) rast(f, win = ext(msk))))
names(rs) <- final_set
if (pred_res > res_m) rs <- crop(aggregate(rs, fact = pred_fact, fun = "mean", na.rm = TRUE),
                                 msk, snap = "out")

omsk <- lapply(c(native = cls_native, tame = cls_tame, annual = cls_annual), function(v) {
  m <- subst(msk & lulc == v, FALSE, NA)
  if (pred_res > res_m) {
    m <- aggregate(m, fact = pred_fact, fun = "max", na.rm = TRUE)
    m <- resample(crop(m, rs, snap = "out"), rs, method = "near")
  }
  names(m) <- "origin_mask"
  m
})

# terra passes each block as a data.frame; the origin_mask column is unused.
gfun <- function(model, data, ...)
  plogis(model$b0 + drop(as.matrix(apply_transforms(data)[names(model$beta)]) %*% model$beta))

# ---- predict ----
preds  <- list()
priors <- list()
for (j in seq_len(nrow(transitions))) {
  from <- transitions$from[j]
  to   <- transitions$to[j]
  tag  <- sprintf("conversion_%s-%s", from, to)
  
  co <- read.csv(file.path(model_dir, paste0("glm_", tag, "_coef.csv")))
  co <- co[match(c("(Intercept)", final_set), co$term), ]
  stopifnot(!anyNA(co$term))
  b    <- setNames(co$Estimate, co$term)
  beta <- b[final_set] / co$scale[-1]
  stopifnot(length(beta) == length(final_set), all(is.finite(beta)))
  g    <- glm_sum[glm_sum$transition == paste(from, to, sep = "-"), ]
  s1   <- g$n_pos / g$n
  b0   <- b[["(Intercept)"]] - sum(beta * co$center[-1]) +
    qlogis(g$base_rate) - qlogis(s1)
  priors[[tag]] <- c(br = g$base_rate, s1 = s1)
  
  t0 <- Sys.time()
  preds[[tag]] <- predict(c(rs, omsk[[from]]), list(b0 = b0, beta = beta),
                          fun = gfun, na.rm = TRUE,
                          filename = pred_path(tag), overwrite = TRUE, wopt = wopt)
  cat(tag, "predicted in", format(Sys.time() - t0), "\n")
  write.table(data.frame(tag = tag, res = pred_res, minutes = as.numeric(Sys.time() - t0, units = "mins")),
              pred_path("timing", ext = ".csv"), sep = ",", row.names = FALSE,
              col.names = !file.exists(pred_path("timing", ext = ".csv")),
              append = file.exists(pred_path("timing", ext = ".csv")))
}
# conversion_native-annual predicted in 30.34788 mins 
# conversion_tame-annual predicted in 30.36589 mins 
# conversion_native-settlement predicted in 30.94215 mins 
# conversion_tame-settlement predicted in 30.32889 mins 
# conversion_native-tame predicted in 31.48763 mins 
# conversion_annual-tame predicted in 33.5916 mins 
# conversion_annual-settlement predicted in 31.90372 mins

# ---- one surface per destination ----
for (dest in unique(transitions$to)) {
  p <- preds[transitions$to == dest]
  dest_r <- if (length(p) == 1) p[[1]] else
    writeRaster(Reduce(cover, p), pred_path(dest), overwrite = TRUE, wopt = wopt)
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
# conversion_native-annual     predicted mean 0.0034 (training pop.) vs observed 0.0030 | 0.0030 whole domain
# conversion_tame-annual       predicted mean 0.1473 (training pop.) vs observed 0.1732 | 0.1545 whole domain
# conversion_native-settlement predicted mean 0.0049 (training pop.) vs observed 0.0037 | 0.0044 whole domain
# conversion_tame-settlement   predicted mean 0.0078 (training pop.) vs observed 0.0079 | 0.0096 whole domain
# conversion_native-tame       predicted mean 0.0062 (training pop.) vs observed 0.0055 | 0.0060 whole domain
# conversion_annual-tame       predicted mean 0.0167 (training pop.) vs observed 0.0173 | 0.0215 whole domain
# conversion_annual-settlement predicted mean 0.0030 (training pop.) vs observed 0.0023 | 0.0039 whole domain

# ---- calibration-in-the-large ----
# After the prior correction, the population mean still sits above br: +3% tame->annual up to
# +45% annual->settlement (QA, 30 m run). One logit shift per transition sets the mean to br on
# the population br was computed on (7c: origin, unprotected, predictor-valid, map20 in {from, to},
# uncertain pixels included). Ranking and AUC are unchanged.
map20 <- crop(rast(tagged_path("./processeddatasets", "landcover", "lulc_2020.tif", "study-area")),
              msk, snap = "out"); levels(map20) <- NULL
pa30  <- crop(rast(pa_file), msk, snap = "out")   # pa above is resampled in a preview run

cal <- list()
set.seed(42)
for (j in seq_len(nrow(transitions))) {
  from <- transitions$from[j]
  to   <- transitions$to[j]
  tag  <- sprintf("conversion_%s-%s", from, to)
  pop  <- subst(msk & pa30 == 0 & lulc == cls[[from]] & map20 %in% cls[c(from, to)], FALSE, NA)
  xy   <- spatSample(pop, 1e6, na.rm = TRUE, xy = TRUE)[, c("x", "y")]
  p    <- na.omit(terra::extract(preds[[tag]], xy, ID = FALSE)[[1]])
  br   <- priors[[tag]][["br"]]
  delta <- uniroot(\(d) mean(plogis(qlogis(p) + d)) - br, c(-3, 1))$root
  cal[[tag]] <- app(preds[[tag]], \(x) plogis(qlogis(x) + delta),
                    filename = pred_path(paste0(tag, "_cal")), overwrite = TRUE, wopt = wopt)
  cat(sprintf("%-28s delta %.3f | mean %.4f -> %.4f (br %.4f, n %d)\n",
              tag, delta, mean(p), mean(plogis(qlogis(p) + delta)), br, length(p)))
}

for (dest in unique(transitions$to)) {
  p <- cal[transitions$to == dest]
  dest_r <- if (length(p) == 1) p[[1]] else
    writeRaster(Reduce(cover, p), pred_path(paste0(dest, "_cal")), overwrite = TRUE, wopt = wopt)
  if (pred_res == res_m)
    aggregate(dest_r, fact = out_fact, fun = "mean", na.rm = TRUE,
              filename = pred_path(paste0(dest, "_cal"), out_res), overwrite = TRUE, wopt = wopt)
}