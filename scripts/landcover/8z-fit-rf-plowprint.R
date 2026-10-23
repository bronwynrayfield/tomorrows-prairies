# =============================================================================
# *Kept as a legacy - only one run once early on - but not updated
# Purpose:
#   Random forest on the 20 retained climate/biophysical predictors against the
#   Plowprint benchmark response from 7-response-plowprint.R. This is the
#   PIPELINE CHECK described there: it tests that predictor extraction, fold
#   assignment and spatial CV behave, before the same machinery is trusted on
#   the conversion-risk target.
#
#   Plowprint is cumulative ("ever cultivated"), so the response is close to a
#   soil-and-climate suitability surface by construction. A high AUC here means
#   the pipeline works, NOT that conversion is predictable. Do not read these
#   numbers across to the conversion model -- the positive class is a different
#   population and there is no conversion interval.
#
#   The training population comes from the response raster itself (LUTS
#   grassland, wetland and cropland, minus forest/water/settlement) intersected
#   with the predictor-valid mask. It is NOT the baseline-grassland population
#   used by the conversion model, so no grassland mask is applied here.
#
# Outputs:
#   - rf_plowprint.rds                 : fitted ranger model
#   - rf_plowprint_importance.csv      : permutation importance
#   - rf_plowprint_loeo.csv            : per-ecodistrict CV results
#   - train-sample_plowprint.rds       : the sampled training frame
# =============================================================================

source(file.path("scripts/landcover", "0a-setup.R"))
source(file.path("scripts/landcover", "0b-predictor-files.R"))

# ---- inputs ----
in_sa         <- file.path("./processeddatasets", "predictors", "30m")
out_dir       <- "output/landcover/model/plowprint"
retained_file <- file.path(in_sa, "retained-covariates.rds")
mask_file     <- file.path(in_sa, "predictor-valid-mask_study-area.tif")
resp_file     <- tagged_path("./processeddatasets", "response/response-plowprint",
                             "response_plowprint.tif", "30m", "255")
train_file    <- file.path(out_dir, "train-sample_plowprint.rds")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Sampling. The Plowprint classes are not rare (unlike the conversion
# response), so "proportional" is the default: it preserves the base rate,
# keeps predicted probabilities calibrated, and still yields plenty of both
# classes. "balanced" is kept for the conversion target, where the base rate
# makes proportional sampling impractical; note that it makes probabilities
# relative rather than absolute (AUC is rank-based, so it is unaffected).
sampling     <- "proportional"
n_total   <- 500000 
n_class      <- 150000    # per class if balanced
num_trees <- 500     
cv_trees  <- 250     
min_fold_n   <- 500       # skip folds smaller than this
min_fold_pos <- 25        # ...or with fewer positives than this

sa_r <- rast(file.path(sa_filepath, sa_raster_filename))

# ---- training population: response non-NA, predictor-valid ----
# Two separate constraints. The response raster is NA outside the LUTS
# grassland/wetland/cropland set, so it defines which pixels have a response at
# all; the valid mask defines which have a value for every predictor. The model
# trains on the intersection. NB the valid mask drops ~6.8% of the study area
# to DSS NAs concentrated along SK waterways (SKSID miscellaneous units), so
# those cells are absent here too.
#
# This population differs from the conversion model's, which is baseline
# grassland only: there, cropland pixels cannot convert and are ineligible,
# whereas here "ever cultivated" makes them legitimate positives.
resp <- rast(resp_file)
names(resp) <- "plowed"

msk  <- rast(mask_file)
resp <- crop(resp, msk, snap = "out")
train_msk  <- subst(msk & !is.na(resp), FALSE, NA)

final_set <- readRDS(retained_file)
rs <- rast(predictor_files[final_set])
names(rs) <- final_set   # rast() takes names from the files' band names

# ---- sample ----
# Manual oversampling: the training population is a modest fraction of the
# raster, and spatSample under-delivers with na.rm = TRUE at low valid
# fractions.
set.seed(42)

draw <- function(m, n) {
  p <- spatSample(m, n * 20, method = "random", xy = TRUE, values = TRUE,
                  na.rm = FALSE)
  p <- as.matrix(p[!is.na(p[[3]]), c("x", "y")])
  p[sample(nrow(p), min(n, nrow(p))), , drop = FALSE]
}

if (sampling == "balanced") {
  resp_d <- mask(resp, train_msk)
  pts <- rbind(draw(subst(resp_d == 1, FALSE, NA), n_class),
               draw(subst(resp_d == 0, FALSE, NA), n_class))
} else {
  pts <- draw(train_msk, n_total)
}

d <- terra::extract(rs, pts)
d$ID <- NULL
d$plowed <- terra::extract(resp, pts)[[1]]
d <- cbind(pts, d)
d <- na.omit(d)

# ---- ecodistrict folds ----
ecodistricts  <- st_read(tagged_path("./processeddatasets", ecodistricts_filepath,
                                     ecodistricts_filename, "modelling-extent"))
ecodistricts  <- st_transform(ecodistricts, crs(sa_r))
ecodistrict_r <- rasterize(vect(ecodistricts), sa_r, field = "ECODISTRIC")

d$ecodistrict <- factor(terra::extract(ecodistrict_r, d[, c("x", "y")])[, 2])
d <- d[!is.na(d$ecodistrict), ]
d$ecodistrict <- droplevels(d$ecodistrict)

d$plowed <- factor(d$plowed, levels = c(0, 1))

saveRDS(d, train_file)

cat("n =", nrow(d),
    "| positives =", sum(d$plowed == 1),
    sprintf("(%.3f)", mean(d$plowed == 1)),
    "| ecodistricts =", nlevels(d$ecodistrict), "\n")

# ---- fit on all data (OOB) ----
f <- as.formula(paste("plowed ~", paste(final_set, collapse = " + ")))

fit <- ranger(f, data = d, probability = TRUE, num.trees = num_trees,
              importance = "permutation", seed = 42,
              num.threads = parallel::detectCores())

oob_auc <- as.numeric(auc(roc(d$plowed, fit$predictions[, "1"], quiet = TRUE)))

saveRDS(fit, file.path(out_dir, "rf_plowprint.rds"))

imp <- data.frame(variable   = names(fit$variable.importance),
                  importance = unname(fit$variable.importance))
imp <- imp[order(-imp$importance), ]
write.csv(imp, file.path(out_dir, "rf_plowprint_importance.csv"), row.names = FALSE)

# ---- leave-one-ecodistrict-out CV ----
# OOB holds out rows at random, so climate varies within the held-out set and
# can discriminate. LOEO holds out a whole ecodistrict, within which climate is
# near-constant (ecodistrict R2 ~0.93 for prtot/tx_mean), so it cannot. The gap
# between the two is expected and mechanistic; report both.
folds <- levels(d$ecodistrict)

loeo <- lapply(folds, function(k) {
  te <- d[d$ecodistrict == k, ]
  if (nrow(te) < min_fold_n || sum(te$plowed == 1) < min_fold_pos ||
      length(unique(te$plowed)) < 2) {
    return(data.frame(ecodistrict = k, n = nrow(te),
                      n_pos = sum(te$plowed == 1), auc = NA_real_))
  }
  tr <- d[d$ecodistrict != k, ]
  m  <- ranger(f, data = tr, probability = TRUE, num.trees = cv_trees,
               seed = 42, num.threads = parallel::detectCores())
  p  <- predict(m, te)$predictions[, "1"]
  data.frame(ecodistrict = k, n = nrow(te), n_pos = sum(te$plowed == 1),
             auc = as.numeric(auc(roc(te$plowed, p, quiet = TRUE))))
})
loeo <- do.call(rbind, loeo)

write.csv(loeo, file.path(out_dir, "rf_plowprint_loeo.csv"), row.names = FALSE)

# Unweighted mean gives each ecodistrict one vote; the n-weighted mean is
# dominated by the largest districts. They answer different questions, so both
# are printed.
loeo_ok <- loeo[!is.na(loeo$auc), ]

cat("\nOOB AUC                 :", round(oob_auc, 3), "\n")
cat("LOEO AUC (mean of folds):", round(mean(loeo_ok$auc), 3),
    sprintf("(SD %.3f, %d/%d folds evaluated)\n",
            sd(loeo_ok$auc), nrow(loeo_ok), nrow(loeo)))
cat("LOEO AUC (n-weighted)   :",
    round(weighted.mean(loeo_ok$auc, loeo_ok$n), 3), "\n\n")
print(head(imp, 10), row.names = FALSE)