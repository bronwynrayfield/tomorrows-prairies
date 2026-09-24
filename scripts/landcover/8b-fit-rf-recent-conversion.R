# Purpose:
#   Random forest on the 31 retained predictors (biophysical, landscape context
#   and socio-economic) against the AAFC-ACI-derived
#   grassland-to-cropland conversion response from 7-response.R, fit separately
#   for each 2010 baseline class. Run after the plowprint PIPELINE CHECK.
#
#   Strata are fit separately because a pooled balanced draw would take its
#   positives 98.6% from tame (14,091,827 tame positives against 203,918
#   native), so native conversion would contribute almost nothing to where the
#   splits go. The two 2010 classes are disjoint, so per-stratum predictions
#   mosaic into the single to-cropland surface CLUES needs.
#
#   Training population is three constraints:
#     (1) 2010 baseline class == focal stratum. The response codes annual
#         cropland, other and settlement as 0, so without this the frame gains
#         ~480M ineligible negatives.
#     (2) Outside CPCAD 2010. Protected pixels are structural zeros: native
#         0.0100% inside vs 0.276% outside (28x), tame 3.50% vs 15.12% (4.3x).
#         Leaving them in lets the forest absorb protection into the covariates
#         protection correlates with (thin, steep, dry land), which deflates
#         predicted risk on unprotected marginal land. The mask costs 0.94% of
#         native positives and 0.80% of tame.
#     (3) Predictor-valid (DSS NAs, ~6.8% of the study area).
#
#   Sampling is balanced for both strata, so predicted probabilities are
#   RELATIVE, not absolute, and Brier is not comparable to the plowprint run.
#   The prior correction back to each stratum's own base rate (0.296% native,
#   15.14% tame) happens in the prediction script, before the mosaic. Without
#   it, a native and a tame pixel with equal predicted values read as equally
#   at risk when their true rates differ by ~51x.
#
#   Resuming: the training sample and completed LOEO folds are reused if their
#   files exist, so an overnight crash restarts where it stopped. To rerun a
#   stratum from scratch, delete its train-sample, _loeo.csv and _oof.csv.
#
# Outputs (stratum-tagged):
#   - train-sample_conversion_<stratum>.rds       : the sampled training frame
#   - rf_conversion_<stratum>.rds                 : fitted ranger model (full sample)
#   - rf_conversion_<stratum>_importance.csv      : permutation importance (full sample)
#   - rf_conversion_<stratum>_thin.csv            : OOB AUC per thinning replicate
#   - rf_conversion_<stratum>_thin_importance.csv : permutation importance, mean over replicates
#   - rf_conversion_<stratum>_correlogram.rds     : residual spline correlogram, thinned fit
#   - rf_conversion_<stratum>_loeo.csv            : per-ecodistrict CV results
#   - rf_conversion_<stratum>_oof.csv             : out-of-fold predictions
#   - rf_conversion_<stratum>_thin_oob.csv        : thinned refit predictions
#   - rf_conversion_<stratum>_calibration.csv     : binned calibration, thinned OOB and LOEO
#   - rf_conversion_<stratum>_calibration.png     : calibration plot

source(file.path("scripts/landcover", "0a-setup.R"))
source(file.path("scripts/landcover", "0b-predictor-files.R"))

strata <- c("native", "tame")            # run in order; native first as the priority

# ---- inputs ----
# study area
sa_r <- rast(file.path(sa_filepath, sa_raster_filename))

in_sa         <- file.path("./processeddatasets", "predictors", "30m")
out_dir       <- "output/landcover/model"
retained_file <- file.path(in_sa, "retained-covariates.rds")
mask_file     <- file.path(in_sa, "predictor-valid-mask_study-area.tif")
resp_file     <- tagged_path("./processeddatasets", "response",
                             "response.tif", "modelling-extent")
lulc_file     <- tagged_path("./processeddatasets", "landcover", "lulc_2010.tif", "30m")
pa_file       <- file.path("processeddatasets", "management","Protected Areas Combined", "CPCAD_pa_2010.tif")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Balanced for both strata. Native at 0.276% cannot be sampled proportionally,
# and the two strata must be fit identically for the mosaic to be coherent.
sampling     <- "balanced"
n_class      <- 150000    # per class. Native ceiling is ~188k positives after
# all three masks, so the native run sits near
# exhaustive; raising this will fail the stopifnot.
num_trees    <- 500
cv_trees     <- num_trees / 2   # CV only: 101 s per fold at 250 trees on native,
# ~2.5-3 h for the full LOEO. Reported model uses num_trees.

# A fold's AUC needs both classes: 752 is 98% positive (237 negatives of 10,625),
# 833 is 3.5%, so a size or positives-only floor passes folds with too few of one
# class to rank against. On native this excludes 22 of 130 folds; 31 at a floor
# of 50. Stricter floors can be applied afterwards from the _oof.csv; looser
# ones need a rerun, so this is the minimum reportable.
min_class_n  <- 25

# Thinning grid = the landscape-context focal window, read from the proportion
# layer paths so it cannot drift from 4d. One row per window means no two
# training rows share neighbourhood information.
focal_win_m  <- as.numeric(sub(".*proportion_(\\d+)m.*", "\\1", predictor_files[["crop"]]))
stopifnot(!is.na(focal_win_m))
n_thin_reps  <- 5         # each replicate redraws which row survives per cell

# spline.correlog is O(n^2): 300k rows asked for 670 Gb. 10k is ~0.8 Gb.
n_corr       <- 10000
corr_xmax    <- 100000    # m; native residual correlation is flat well before this

n_threads    <- parallel::detectCores(logical = FALSE)

# pROC's default direction = "auto" flips any fold scoring below 0.5 to above
# it, which biases fold AUCs upward. Fix the direction: controls 0, cases 1.
auc_fixed <- function(obs, p) {
  as.numeric(auc(roc(as.integer(as.character(obs)), p, levels = c(0, 1),
                     direction = "<", quiet = TRUE)))
}

append_csv <- function(x, file) {
  write.table(x, file, sep = ",", row.names = FALSE,
              col.names = !file.exists(file), append = file.exists(file))
}

tss_at <- function(obs, p, t) mean(p[obs == 1] >= t) + mean(p[obs == 0] < t) - 1

calib <- function(x) {
  x$bin <- cut(x$p, breaks = seq(0, 1, 0.1), include.lowest = TRUE)
  cbind(aggregate(cbind(p, obs) ~ bin, x, mean),
        n = aggregate(obs ~ bin, x, length)$obs)
}

# ---- predictors (shared across strata) ----
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

# 26 layers are on the modelling extent and 5 (the 4d focal proportions) on
# the study area, so a plain rast() fails on extent. win = reads a window
# without materialising a cropped copy. Belongs in 4f; remove once fixed there.
rs <- rast(lapply(predictor_files[final_set], \(f) rast(f, win = ext(sa_r))))
names(rs) <- final_set

f <- as.formula(paste("converted ~", paste(final_set, collapse = " + ")))

for (stratum in strata) {
  
  cls_focal  <- get(paste0("cls_", stratum))
  tag        <- paste0("conversion_", stratum)
  train_file <- file.path(out_dir, paste0("train-sample_", tag, ".rds"))
  out        <- function(suffix) file.path(out_dir, paste0("rf_", tag, suffix))
  cat("\n==== ", stratum, " ====\n", sep = "")
  
  # ---- sample (reused if present, so LOEO folds stay on the same rows) ----
  if (file.exists(train_file)) {
    d <- readRDS(train_file)
  } else {
    # Three separate constraints, see header. The response raster alone is NOT
    # sufficient: it is 0 everywhere outside the converted set, including on
    # classes that were never eligible to convert.
    resp <- crop(rast(resp_file), sa_r, snap = "out")
    names(resp) <- "converted"
    msk  <- mask(crop(rast(mask_file), sa_r, snap = "out"), sa_r)
    pa   <- mask(crop(rast(pa_file)  , sa_r, snap = "out"), sa_r)   # 1 = inside CPCAD 2010
    lulc <- mask(crop(rast(lulc_file), sa_r, snap = "out"), sa_r)
    train_msk <- subst(msk & !is.na(resp) & pa == 0 & lulc == cls_focal, FALSE, NA)
    resp_d <- mask(resp, train_msk)
    rm(resp, msk, pa, lulc, train_msk)
    
    # Positives are enumerated rather than sampled by rejection: native positives
    # outside CPCAD are 202,004 of ~685M valued cells (0.03%), so rejection would
    # need ~10,000x oversampling to fill n_class. Negatives are 10-11% of valued
    # cells in both strata, which rejection handles.
    set.seed(42)
    
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
    
    pos <- draw_all(subst(resp_d == 1, FALSE, NA), n_class)
    neg <- draw_reject(subst(resp_d == 0, FALSE, NA), n_class)
    stopifnot(nrow(pos) >= 0.95 * n_class, nrow(neg) >= 0.95 * n_class)
    pts <- rbind(pos, neg)
    
    d <- terra::extract(rs, pts)
    d$ID <- NULL
    d$converted <- terra::extract(resp_d, pts)[[1]]
    d <- cbind(pts, d)
    d <- na.omit(d)
    rm(pos, neg, pts)
    
    # ---- ecodistrict folds ----
    if (!exists("ecodistrict_r")) {
      ecodistricts_v <- vect(tagged_path("./processeddatasets", ecodistricts_filepath,
                                         ecodistricts_filename, "modelling-extent"))
      ecodistrict_r  <- rasterize(ecodistricts_v, sa_r, field = "ECODISTRIC")
    }
    d$ecodistrict <- factor(terra::extract(ecodistrict_r, d[, c("x", "y")])[, 2])
    d <- d[!is.na(d$ecodistrict), ]
    d$ecodistrict <- droplevels(d$ecodistrict)
    
    d$converted <- factor(d$converted, levels = c(0, 1))
    # conversion rate in the training population (all three masks); the prior
    # correction rescales predictions from sample prevalence to this
    attr(d, "base_rate") <- as.numeric(global(resp_d, "mean", na.rm = TRUE))
    saveRDS(d, train_file)
  }
  rm(resp_d)
  
  # Sample prevalence is set by design (1:1), not estimated. The prior
  # correction needs the realised value, since draw_all() can return fewer
  # rows than requested.
  stopifnot(!is.null(attr(d, "base_rate")))
  cat("n =", nrow(d),
      "| positives =", sum(d$converted == 1),
      sprintf("(%.4f sample prevalence)", mean(d$converted == 1)),
      sprintf("| base rate %.4f%%", 100 * attr(d, "base_rate")),
      "| ecodistricts =", nlevels(d$ecodistrict), "\n")
  
  # ---- full-sample fit: the reported model ----
  fit <- ranger(f, data = d, probability = TRUE, num.trees = num_trees,
                importance = "permutation", seed = 42, num.threads = n_threads)
  fit$base_rate   <- attr(d, "base_rate")
  fit$sample_prev <- mean(d$converted == 1)
  saveRDS(fit, out(".rds"))
  
  imp <- data.frame(variable   = names(fit$variable.importance),
                    importance = unname(fit$variable.importance))
  imp <- imp[order(-imp$importance), ]
  write.csv(imp, out("_importance.csv"), row.names = FALSE)
  
  # Naive OOB is kept for the record only. On native it was 0.995: 150k
  # positives from a ~202k pool at 30 m puts near-neighbour pixels on both
  # sides of each tree's bootstrap, so OOB rows are scored against near-copies.
  oob_naive <- auc_fixed(d$converted, fit$predictions[, "1"])
  
  # ---- spatially thinned OOB: the interpolation estimate ----
  # One row per focal window, then downsample the majority class to restore
  # 1:1. On native (seed 1): 24,536 thinned rows at prevalence 0.348, OOB 0.926;
  # rebalanced to 17,100 rows, OOB 0.923. The prevalence shift moved AUC by
  # 0.003, so the drop from 0.995 is pseudoreplication, not base rate.
  thin_r <- rast(ext(sa_r), resolution = focal_win_m, crs = crs(sa_r))
  cell   <- cellFromXY(thin_r, as.matrix(d[, c("x", "y")]))
  
  thin <- lapply(seq_len(n_thin_reps), function(r) {
    set.seed(r)
    i  <- sample(nrow(d))
    dt <- d[i[!duplicated(cell[i])], ]
    n_min <- min(table(dt$converted))
    keep  <- unlist(lapply(split(seq_len(nrow(dt)), dt$converted),
                           \(j) j[sample.int(length(j), n_min)]))
    db <- dt[keep, ]
    m  <- ranger(f, data = db, probability = TRUE, num.trees = num_trees,
                 importance = "permutation", seed = r, num.threads = n_threads)
    list(summary = data.frame(rep = r, n_thin = nrow(dt),
                              prev_thin = mean(dt$converted == 1),
                              n_bal = nrow(db),
                              oob_auc = auc_fixed(db$converted, m$predictions[, "1"])),
         oob = data.frame(rep = r, x = db$x, y = db$y, ecodistrict = db$ecodistrict,
                          obs = as.integer(db$converted == 1), p = m$predictions[, "1"]),
         imp = m$variable.importance[final_set],
         db = if (r == 1) db, m = if (r == 1) m)
  })
  
  thin_sum <- do.call(rbind, lapply(thin, `[[`, "summary"))
  write.csv(thin_sum, out("_thin.csv"), row.names = FALSE)
  
  thin_oob <- do.call(rbind, lapply(thin, `[[`, "oob"))
  write.csv(thin_oob, out("_thin_oob.csv"), row.names = FALSE)
  
  imp_thin <- data.frame(variable = final_set,
                         importance = rowMeans(sapply(thin, `[[`, "imp")))
  imp_thin <- imp_thin[order(-imp_thin$importance), ]
  write.csv(imp_thin, out("_thin_importance.csv"), row.names = FALSE)
  
  # ---- residual correlogram, thinned fit ----
  # On native the thinned-fit residuals cross zero at ~2-3 km and are flat to
  # 100 km, so the focal-window thinning distance is sufficient. Residuals from
  # the full-sample fit are flat at all lags because that fit absorbed the
  # short-range structure; they are not informative and are not computed here.
  # The thinned sample has few pairs below focal_win_m, so the short-lag end is
  # reconstructed from sparse data; the flat tail carries the conclusion.
  db1 <- thin[[1]]$db
  res <- as.integer(db1$converted == 1) - thin[[1]]$m$predictions[, "1"]
  set.seed(1)
  k   <- sample(nrow(db1), min(nrow(db1), n_corr))
  sc  <- ncf::spline.correlog(db1$x[k], db1$y[k], res[k], resamp = 0, xmax = corr_xmax)
  saveRDS(sc, out("_correlogram.rds"))
  rm(thin, db1, res)
  
  # ---- leave-one-ecodistrict-out CV: the transfer estimate ----
  # LOEO holds out a whole ecodistrict, within which climate is near-constant
  # (ecodistrict R2 ~0.93 for prtot/tx_mean), so it measures extrapolation to
  # unseen places. On native it gave 0.83 (88/130 folds, under an earlier
  # size-based gate), against 0.92 thinned OOB. With no residual dependence
  # beyond ~3 km, that gap is interpolation vs extrapolation, not leakage.
  # Report thinned OOB and LOEO; never report the naive OOB.
  #
  # Each fold appends its predictions to _oof.csv, then its summary row to
  # _loeo.csv. A fold counts as done only once its _loeo row exists, so on
  # resume any _oof rows from a fold that crashed mid-write are dropped.
  loeo_file <- out("_loeo.csv")
  oof_file  <- out("_oof.csv")
  
  # a _loeo.csv from the pre-oof version of this script would mark every fold done
  stopifnot(!file.exists(loeo_file) || "n_neg" %in% names(read.csv(loeo_file, nrows = 1)))
  
  done <- if (file.exists(loeo_file)) read.csv(loeo_file)$ecodistrict else character(0)
  if (file.exists(oof_file)) {
    o <- read.csv(oof_file)
    write.csv(o[o$ecodistrict %in% done, ], oof_file, row.names = FALSE)
    rm(o)
  }
  
  folds <- levels(d$ecodistrict)
  for (k in setdiff(folds, done)) {
    te    <- d[d$ecodistrict == k, ]
    n_pos <- sum(te$converted == 1)
    n_neg <- sum(te$converted == 0)
    row   <- data.frame(ecodistrict = k, n = nrow(te), n_pos = n_pos,
                        n_neg = n_neg, auc = NA_real_)
    if (min(n_pos, n_neg) >= min_class_n) {
      m <- ranger(f, data = d[d$ecodistrict != k, ], probability = TRUE,
                  num.trees = cv_trees, seed = 42, num.threads = n_threads,
                  importance = "none")
      p <- predict(m, te, num.threads = n_threads)$predictions[, "1"]
      append_csv(data.frame(ecodistrict = k,
                            obs = as.integer(te$converted == 1), p = p), oof_file)
      row$auc <- auc_fixed(te$converted, p)
    }
    append_csv(row, loeo_file)
  }
  
  # ---- summary ----
  loeo <- read.csv(loeo_file)
  oof  <- read.csv(oof_file)
  ok   <- !is.na(loeo$auc)
  
  # TSS from LOEO out-of-fold predictions. Max TSS = Youden index. With ~300k
  # rows, tuning one threshold on the scored predictions adds negligible
  # optimism, so the gap to TSS at 0.5 is a calibration signal: on native,
  # 0.64 at threshold 0.19 vs 0.45 at 0.5, because LOEO under-predicts in
  # held-out districts (shortfall grows ~0.11 per tenfold increase in the
  # positives a fold removes from training, holding prevalence fixed).
  r_oof <- roc(oof$obs, oof$p, levels = c(0, 1), direction = "<", quiet = TRUE)
  best  <- coords(r_oof, "best", best.method = "youden",
                  ret = c("threshold", "sensitivity", "specificity"))
  
  # Calibration and Brier from thinned OOB, not LOEO: LOEO removes a district's
  # positive analogues from training, so its probabilities are biased low in a
  # way the deployed model is not. On native, thinned OOB is near-diagonal with
  # mild under-confidence (max gap ~0.07); Brier 0.114, skill 0.54. Both at
  # sample prevalence, so not comparable to plowprint or to a landscape-scale
  # Brier after the prior correction.
  brier <- mean((thin_oob$p - thin_oob$obs)^2)
  cal_t <- calib(thin_oob)
  cal_l <- calib(oof)
  write.csv(rbind(cbind(source = "thinned_oob", cal_t), cbind(source = "loeo", cal_l)),
            out("_calibration.csv"), row.names = FALSE)
  
  png(out("_calibration.png"), width = 1600, height = 1600, res = 300)
  plot(cal_t$p, cal_t$obs, xlim = c(0, 1), ylim = c(0, 1), pch = 19, main = stratum,
       xlab = "Mean predicted", ylab = "Observed fraction converted")
  points(cal_l$p, cal_l$obs, pch = 1)
  abline(0, 1)
  legend("topleft", c("Thinned OOB", "LOEO"), pch = c(19, 1), bty = "n")
  dev.off()
  
  # Pooled AUC ranks all held-out predictions together, so between-fold offsets
  # in predicted probability count as discrimination. Fold prevalence runs from
  # 3.5% to 98% on native, so offsets are likely: read pooled against the fold
  # mean, not instead of it. A large gap points to calibration differences
  # between districts rather than better transfer.
  cat("\nOOB AUC, naive (do not report):", round(oob_naive, 3), "\n")
  cat("OOB AUC, thinned              :", round(mean(thin_sum$oob_auc), 3),
      sprintf("(range %.3f-%.3f over %d reps; %.0f rows at 1:1)\n",
              min(thin_sum$oob_auc), max(thin_sum$oob_auc), n_thin_reps,
              mean(thin_sum$n_bal)))
  cat("LOEO AUC, pooled              :", round(auc_fixed(oof$obs, oof$p), 3), "\n")
  cat("LOEO AUC, mean of folds       :", round(mean(loeo$auc[ok]), 3),
      sprintf("(SD %.3f, %d/%d folds evaluated)\n",
              sd(loeo$auc[ok]), sum(ok), nrow(loeo)))
  cat("LOEO AUC, n-weighted          :",
      round(weighted.mean(loeo$auc[ok], loeo$n[ok]), 3), "\n")
  cat("dropped folds hold            :",
      sprintf("%.1f%% of positives\n", 100 * sum(loeo$n_pos[!ok]) / sum(loeo$n_pos)))
  cat("dropped folds                 :", paste(loeo$ecodistrict[!ok], collapse = ", "), "\n\n")
  cat("TSS (LOEO), max               :", round(best$sensitivity + best$specificity - 1, 3),
      sprintf("(threshold %.3f)\n", best$threshold))
  cat("TSS (LOEO), at 0.5            :", round(tss_at(oof$obs, oof$p, 0.5), 3), "\n")
  cat("Brier (thinned OOB)           :", round(brier, 3),
      sprintf("(skill %.3f vs 0.25 no-skill at 1:1)\n", 1 - brier / 0.25))
  print(head(imp, 10), row.names = FALSE)
  cat("\nthinned importance:\n")
  print(head(imp_thin, 10), row.names = FALSE)
  
  rm(d, fit, imp, imp_thin, thin_sum, loeo, oof, cell, thin_r, thin_oob, r_oof, best, cal_t, cal_l, brier)
  gc()
}


# Plot variable importance
cat_of <- c(prtot = "Climatic", tx_mean = "Climatic", gddgrow_5 = "Climatic",
            oc = "Edaphic", bd = "Edaphic", clay = "Edaphic", ph = "Edaphic",
            drainage = "Edaphic", slope = "Edaphic", kp0 = "Edaphic",
            kp1500 = "Edaphic", ldepth = "Edaphic", psssm_mean = "Edaphic",
            psssm_sd = "Edaphic", tsm_planting = "Edaphic", tsm_growing = "Edaphic",
            tpi = "Topographic", tri = "Topographic", srr = "Topographic",
            northness = "Topographic", eastness = "Topographic",
            unprot_grs = "Landscape", tame = "Landscape", crop = "Landscape",
            dev = "Landscape", fors = "Landscape",
            net_return = "Socio-economic", dist_to_elevator = "Socio-economic",
            dist_to_road = "Socio-economic", dist_to_rail = "Socio-economic",
            dist_to_popctr = "Socio-economic", dist_to_regctr = "Socio-economic")

d$category <- factor(cat_of[as.character(d$variable)],
                     levels = c("Landscape", "Climatic", "Edaphic",
                                "Topographic", "Socio-economic"))

ggplot(d, aes(rel, variable, colour = stratum)) +
  geom_line(aes(group = variable), colour = "grey70") +
  geom_point(size = 2) +
  facet_grid(category ~ ., scales = "free_y", space = "free_y", switch = "y") +
  labs(x = "Relative permutation importance", y = NULL, colour = NULL)





top <- levels(d$variable)[levels(d$variable) %in%
                            head(n$variable[order(-n$importance)], 15)]
d15 <- d[d$variable %in% top, ]
d15$signed <- ifelse(d15$stratum == "Native", -d15$rel, d15$rel)

ggplot(d15, aes(signed, variable, fill = stratum)) +
  geom_col() +
  geom_vline(xintercept = 0, colour = "grey30") +
  scale_x_continuous(labels = abs) +
  labs(x = "Relative permutation importance", y = NULL, fill = NULL)
