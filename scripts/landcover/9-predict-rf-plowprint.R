# =============================================================================
# 9-predict-rf-plowprint.R
# =============================================================================
#
# Purpose:
#   Applies the random forest fitted in 8-fit-rf-plowprint.R across the study
#   area, producing a predicted probability surface (P of being in the
#   previously-or-currently-cultivated class).
#
#   Read this as a suitability-like surface, not a forecast: Plowprint is
#   cumulative, so a high value means "conditions resemble land that has been
#   cultivated", not "this will be converted". It exists to check that the
#   fitted model produces a spatially sensible surface, and to compare against
#   Olimb & Robinson's raster.
#
#   Predicted at 30 m over the same cells the model could train on (response
#   non-NA, predictor-valid), then aggregated to 300 m for display and for the
#   O&R comparison. NB the ~6.8% of the study area lost to DSS NAs along SK
#   waterways is blank here too.
#
# Outputs:
#   - pred_plowprint_30m.tif   : predicted probability, native resolution
#   - pred_plowprint_300m.tif  : mean of the 30 m surface, for display
# =============================================================================

source(file.path("scripts/landcover", "0a-setup.R"))
source(file.path("scripts/landcover", "0b-predictor-files.R"))

terraOptions(memfrac = 0.5)

# ---- inputs ----
# pred_res: resolution to predict at, in metres.
#   res_m (30) - predict at native resolution, then aggregate to out_res for
#                display. The real output; a multi-hour job.
#   300        - aggregate the predictors first, then predict. Minutes rather
#                than hours, for a quick look at whether the surface is
#                spatially sensible before committing to the full run.
#
# The two are NOT identical: averaging predictors and then predicting is not
# the same as predicting and then averaging, because the forest is non-linear.
# The coarse preview also smooths away fine structure the model responds to
# (terrain especially). Treat it as a preview, not a substitute.
#
# NB `res_m` is the native resolution from 0a-setup.R
pred_res      <- res_m   # run a preview at 300 before full resolution
out_res       <- 300     # display/comparison resolution for the native run

in_sa         <- file.path("./processeddatasets", "predictors", "30m")
model_dir     <- "output/landcover/model"
out_pred      <- file.path("./processeddatasets", "predictions", "plowprint")
retained_file <- file.path(in_sa, "retained-covariates.rds")
mask_file     <- file.path(in_sa, "predictor-valid-mask_study-area.tif")
resp_file     <- tagged_path("./processeddatasets", "response/response-plowprint",
                             "response_plowprint.tif", "30m", "255")
model_file    <- file.path(model_dir, "rf_plowprint.rds")

# aggregation factors derived from the native resolution in 0a-setup.R
stopifnot(pred_res %% res_m == 0, out_res %% res_m == 0)
pred_fact <- pred_res / res_m    # predictor aggregation before a coarse run
out_fact  <- out_res  / res_m    # output aggregation after a native run

dir.create(out_pred, recursive = TRUE, showWarnings = FALSE)

gopt <- c("COMPRESS=DEFLATE", "ZLEVEL=9", "TILED=YES",
          "BLOCKXSIZE=256", "BLOCKYSIZE=256",
          "BIGTIFF=YES", "NUM_THREADS=ALL_CPUS", "PREDICTOR=3")

fit       <- readRDS(model_file)
final_set <- readRDS(retained_file)

stopifnot(setequal(fit$forest$independent.variable.names, final_set))

# ---- masks ----
msk  <- rast(mask_file)
resp <- crop(rast(resp_file), msk, snap = "out")
pred_msk <- subst(msk & !is.na(resp), FALSE, NA)

# ---- predictor stack ----
# Order matters for I/O: aggregate BEFORE cropping, so a coarse run never
# materialises a 20-layer copy at 30 m. Masking is applied to the single-layer
# prediction at the end rather than to the stack, for the same reason -- one
# layer written instead of twenty.
rs <- rast(predictor_files[final_set])
names(rs) <- final_set   # rast() takes names from the files' band names

if (pred_res > res_m) {
  rs       <- aggregate(rs, fact = pred_fact, fun = "mean", na.rm = TRUE)
  pred_msk <- aggregate(pred_msk, fact = pred_fact, fun = "max", na.rm = TRUE)
}

rs       <- crop(rs, pred_msk, snap = "out")
pred_msk <- resample(pred_msk, rs, method = "near")

# ---- predict ----
# terra::predict feeds blocks of cells to the function, so the whole raster is
# never in memory. na.rm = TRUE skips incomplete rows rather than erroring;
# ranger's predict returns a matrix of class probabilities, and we keep the
# column for class 1.
pfun <- function(model, ...) {
  predict(model, ..., num.threads = parallel::detectCores())$predictions[, "1"]
}

pred_file <- file.path(out_pred,
                       sprintf("pred_plowprint_%dm%s.tif", pred_res,
                               if (pred_res > res_m) "_preview" else ""))

pred <- predict(rs, fit, fun = pfun, na.rm = TRUE)
pred <- mask(pred, pred_msk,
             filename = pred_file, overwrite = TRUE,
             wopt = list(gdal = gopt, datatype = "FLT4S"))

# ---- coarse version for display and comparison (30-m run only) ----
if (pred_res == res_m) {
  aggregate(pred, fact = out_fact, fun = "mean", na.rm = TRUE,
            filename = file.path(out_pred, sprintf("pred_plowprint_%dm.tif", out_res)),
            overwrite = TRUE, wopt = list(gdal = gopt, datatype = "FLT4S"))
}

# ---- QA ----
# Predicted mean should sit near the observed base rate if the model is
# calibrated (it should be: the training sample was drawn proportionally).
cat("resolution         :", pred_res, "m\n")
# resolution         : 300 m
print(global(pred, c("mean", "min", "max"), na.rm = TRUE))
# mean min max
# lyr1 0.6999891   0   1
plot(pred)
hist(values(pred), breaks = 50)

obs_msk <- subst(rast(mask_file) & !is.na(crop(rast(resp_file), rast(mask_file), snap = "out")), FALSE, NA)
obs <- mask(crop(rast(resp_file), obs_msk, snap = "out"), obs_msk)
cat("observed base rate :", round(global(obs, "mean", na.rm = TRUE)[[1]], 3), "\n")
# observed base rate : 0.713 