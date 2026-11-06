# Purpose:
#   Random forest on the retained predictors (biophysical, landscape context
#   and socio-economic) for each conversion transition, fit on the shared
#   per-origin samples from 7c-fitting-sample.R (persistence plus one
#   destination, option A). Run after the plowprint PIPELINE CHECK. The GLM
#   benchmark (8a) uses the same rows and ecodistrict folds, so fold AUCs pair.
#
#   Origins are fit separately because a pooled balanced draw would take its
#   positives 98.1% from tame (9,627,153 tame->annual positives against
#   185,947 native->annual in the training population), so native conversion
#   would contribute almost nothing to where the splits go. The two 2010
#   classes are disjoint, so per-origin predictions mosaic into one surface
#   per destination.
#
#   Training population and sampling: see 7c-fitting-sample.R.
#
#   Sampling is balanced, so predicted probabilities are RELATIVE, not
#   absolute, and Brier is not comparable to the plowprint run. The prior
#   correction back to each transition's base rate (read from the 7c sample
#   attribute: native->annual 0.2983%, tame->annual 14.35%) happens in the
#   prediction script, before the mosaic. Without it, a native and a tame
#   pixel with equal predicted values read as equally at risk when their true
#   rates differ by ~48x.
#
#   Resuming: completed LOEO folds are reused if their files exist, so an
#   overnight crash restarts where it stopped. To rerun a transition from
#   scratch, delete its _loeo.csv and _oof.csv.
#
# Outputs (tag = conversion_<from>-<to>):
#   - rf_<tag>.rds                 : fitted ranger model (full sample)
#   - rf_<tag>_importance.csv      : permutation importance (full sample)
#   - rf_<tag>_thin.csv            : OOB AUC per thinning replicate
#   - rf_<tag>_thin_importance.csv : permutation importance, mean over replicates
#   - rf_<tag>_correlogram.rds     : residual spline correlogram, thinned fit
#   - rf_<tag>_loeo.csv            : per-ecodistrict CV results
#   - rf_<tag>_oof.csv             : out-of-fold predictions
#   - rf_<tag>_thin_oob.csv        : thinned refit predictions
#   - rf_<tag>_calibration.csv     : binned calibration, thinned OOB and LOEO
#   - rf_<tag>_calibration.png     : calibration plot

source(file.path("scripts/landcover", "0a-setup.R"))
source(file.path("scripts/landcover", "0b-predictor-files.R"))

# Each transition: persistence + one destination from the 7c per-origin sample (option A).
#transitions <- data.frame(from = c("native", "tame",   "native",     "tame",       "native", "annual", "annual"),
#                          to   = c("annual", "annual", "settlement", "settlement", "tame",   "tame", "settlement"))

transitions <- data.frame(from = c("annual", "annual"),
                          to   = c("tame", "settlement"))

# ---- inputs ----
sa_r <- rast(file.path(sa_filepath, sa_raster_filename))   # thinning grid extent

in_sa         <- file.path("./processeddatasets", "predictors", "30m")
out_dir       <- "output/landcover/model"
retained_file <- file.path(in_sa, "retained-covariates.rds")

num_trees    <- 500
cv_trees     <- num_trees / 2   # CV only: 101 s per fold at 250 trees on native
# (previous run), ~3.5 h for ~126 folds. Reported
# model uses num_trees.

# A fold's AUC needs both classes: in the previous run, 752 was 98% positive
# (237 negatives of 10,625) and 833 was 3.5%, so a size or positives-only floor
# passes folds with too few of one class to rank against. On native this
# excluded 22 of 130 folds (previous ecodistrict layer); 31 at a floor of 50.
# Stricter floors can be applied afterwards from the _oof.csv; looser ones need
# a rerun, so this is the minimum reportable. Same value as 8a.
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

# ---- predictors (shared across transitions) ----
final_set <- readRDS(retained_file)
stopifnot(all(final_set %in% names(predictor_files)))

f <- as.formula(paste("converted ~", paste(final_set, collapse = " + ")))

for (j in seq_len(nrow(transitions))) {
  from <- transitions$from[j]
  to   <- transitions$to[j]
  tag  <- sprintf("conversion_%s-%s", from, to)
  out  <- function(suffix) file.path(out_dir, paste0("rf_", tag, suffix))
  cat("\n==== ", from, " -> ", to, " ====\n", sep = "")
  
  # ---- sample (7c) ----
  # attributes are dropped by row subsetting, so read base_rate first
  d  <- readRDS(file.path(out_dir, paste0("train-sample_", from, ".rds")))
  br <- attr(d, "base_rate")[[to]]
  d  <- d[d$lulc_2020 %in% c(from, to), ]
  d$converted   <- factor(as.integer(d$lulc_2020 == to), levels = c(0, 1))
  d$ecodistrict <- droplevels(d$ecodistrict)
  
  # Sample prevalence is set by design (1:1), not estimated. The prior
  # correction needs the realised value, since 7c can return fewer rows than
  # requested.
  cat("n =", nrow(d),
      "| positives =", sum(d$converted == 1),
      sprintf("(%.4f sample prevalence)", mean(d$converted == 1)),
      sprintf("| base rate %.4f%%", 100 * br),
      "| ecodistricts =", nlevels(d$ecodistrict), "\n")
  
  # ---- full-sample fit: the reported model ----
  fit <- ranger(f, data = d, probability = TRUE, num.trees = num_trees,
                importance = "permutation", seed = 42, num.threads = n_threads, case.weights = d$w)
  fit$base_rate   <- br
  fit$sample_prev <- mean(d$converted == 1)
  saveRDS(fit, out(".rds"))
  
  imp <- data.frame(variable   = names(fit$variable.importance),
                    importance = unname(fit$variable.importance))
  imp <- imp[order(-imp$importance), ]
  write.csv(imp, out("_importance.csv"), row.names = FALSE)
  
  # Naive OOB is kept for the record only. On native it was 0.995 (previous
  # run): 150k positives from a ~186k pool at 30 m puts near-neighbour pixels
  # on both sides of each tree's bootstrap, so OOB rows are scored against
  # near-copies.
  oob_naive <- auc_fixed(d$converted, fit$predictions[, "1"])
  
  # ---- spatially thinned OOB: the interpolation estimate ----
  # One row per focal window, then downsample the majority class to restore
  # 1:1. On native (previous run, seed 1): 24,536 thinned rows at prevalence
  # 0.348, OOB 0.926; rebalanced to 17,100 rows, OOB 0.923. The prevalence
  # shift moved AUC by 0.003, so the drop from 0.995 is pseudoreplication,
  # not base rate.
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
                 importance = "permutation", seed = r, num.threads = n_threads, case.weights = db$w)
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
  # On native (previous run) the thinned-fit residuals cross zero at ~2-3 km
  # and are flat to 100 km, so the focal-window thinning distance is
  # sufficient. Residuals from the full-sample fit are flat at all lags because
  # that fit absorbed the short-range structure; they are not informative and
  # are not computed here. The thinned sample has few pairs below focal_win_m,
  # so the short-lag end is reconstructed from sparse data; the flat tail
  # carries the conclusion.
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
  # unseen places. On native (previous run) it gave 0.83 (88/130 folds, under
  # an earlier size-based gate), against 0.92 thinned OOB. With no residual
  # dependence beyond ~3 km, that gap is interpolation vs extrapolation, not
  # leakage. Report thinned OOB and LOEO; never report the naive OOB.
  #
  # Each fold appends its predictions to _oof.csv, then its summary row to
  # _loeo.csv. A fold counts as done only once its _loeo row exists, so on
  # resume any _oof rows from a fold that crashed mid-write are dropped.
  loeo_file <- out("_loeo.csv")
  oof_file  <- out("_oof.csv")
  
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
                  importance = "none", case.weights = d$w[d$ecodistrict != k])
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
  # optimism, so the gap to TSS at 0.5 is a calibration signal: on native
  # (previous run), 0.64 at threshold 0.19 vs 0.45 at 0.5, because LOEO
  # under-predicts in held-out districts (shortfall grows ~0.11 per tenfold
  # increase in the positives a fold removes from training, holding
  # prevalence fixed).
  r_oof <- roc(oof$obs, oof$p, levels = c(0, 1), direction = "<", quiet = TRUE)
  best  <- coords(r_oof, "best", best.method = "youden",
                  ret = c("threshold", "sensitivity", "specificity"))[1, ]
  
  # Calibration and Brier from thinned OOB, not LOEO: LOEO removes a district's
  # positive analogues from training, so its probabilities are biased low in a
  # way the deployed model is not. On native (previous run), thinned OOB is
  # near-diagonal with mild under-confidence (max gap ~0.07); Brier 0.114,
  # skill 0.54. Both at sample prevalence, so not comparable to plowprint or to
  # a landscape-scale Brier after the prior correction.
  brier <- mean((thin_oob$p - thin_oob$obs)^2)
  cal_t <- calib(thin_oob)
  cal_l <- calib(oof)
  write.csv(rbind(cbind(source = "thinned_oob", cal_t), cbind(source = "loeo", cal_l)),
            out("_calibration.csv"), row.names = FALSE)
  
  png(out("_calibration.png"), width = 1600, height = 1600, res = 300)
  plot(cal_t$p, cal_t$obs, xlim = c(0, 1), ylim = c(0, 1), pch = 19, main = tag,
       xlab = "Mean predicted", ylab = "Observed fraction converted")
  points(cal_l$p, cal_l$obs, pch = 1)
  abline(0, 1)
  legend("topleft", c("Thinned OOB", "LOEO"), pch = c(19, 1), bty = "n")
  dev.off()
  
  # Pooled AUC ranks all held-out predictions together, so between-fold offsets
  # in predicted probability count as discrimination.
  # cat("\nOOB AUC, naive (do not report):", round(oob_naive, 3), "\n")
  # cat("OOB AUC, thinned              :", round(mean(thin_sum$oob_auc), 3),
  #     sprintf("(range %.3f-%.3f over %d reps; %.0f rows at 1:1)\n",
  #             min(thin_sum$oob_auc), max(thin_sum$oob_auc), n_thin_reps,
  #             mean(thin_sum$n_bal)))
  # cat("LOEO AUC, pooled              :", round(auc_fixed(oof$obs, oof$p), 3), "\n")
  # cat("LOEO AUC, mean of folds       :", round(mean(loeo$auc[ok]), 3),
  #     sprintf("(SD %.3f, %d/%d folds evaluated)\n",
  #             sd(loeo$auc[ok]), sum(ok), nrow(loeo)))
  # cat("LOEO AUC, n-weighted          :",
  #     round(weighted.mean(loeo$auc[ok], loeo$n[ok]), 3), "\n")
  # cat("dropped folds hold            :",
  #     sprintf("%.1f%% of positives\n", 100 * sum(loeo$n_pos[!ok]) / sum(loeo$n_pos)))
  # cat("dropped folds                 :", paste(loeo$ecodistrict[!ok], collapse = ", "), "\n\n")
  # cat("TSS (LOEO), max               :", round(best$sensitivity + best$specificity - 1, 3),
  #     sprintf("(threshold %.3f)\n", best$threshold))
  # cat("TSS (LOEO), at 0.5            :", round(tss_at(oof$obs, oof$p, 0.5), 3), "\n")
  # cat("Brier (thinned OOB)           :", round(brier, 3),
  #     sprintf("(skill %.3f vs 0.25 no-skill at 1:1)\n", 1 - brier / 0.25))
  # print(head(imp, 10), row.names = FALSE)
  # cat("\nthinned importance:\n")
  # print(head(imp_thin, 10), row.names = FALSE)
  
  rm(d, fit, imp, imp_thin, thin_sum, loeo, oof, cell, thin_r, thin_oob, r_oof, best, cal_t, cal_l, brier)
  gc()
}


# Read back in the results and print summary
auc_fixed <- function(obs, p) as.numeric(auc(roc(as.integer(as.character(obs)), p,
                                                 levels = c(0, 1), direction = "<", quiet = TRUE)))
tss_at <- function(obs, p, t) mean(p[obs == 1] >= t) + mean(p[obs == 0] < t) - 1

for (tag in sprintf("conversion_%s-%s", transitions$from, transitions$to)) {
  out <- function(s) file.path(out_dir, paste0("rf_", tag, s))
  if (!file.exists(out("_loeo.csv"))) { cat("\n", tag, ": not run\n"); next }
  fit      <- readRDS(out(".rds"))
  thin_sum <- read.csv(out("_thin.csv"))
  thin_oob <- read.csv(out("_thin_oob.csv"))
  loeo     <- read.csv(out("_loeo.csv"))
  oof      <- read.csv(out("_oof.csv"))
  ok       <- !is.na(loeo$auc)
  obs_fit  <- fit$predictions[, "1"]   # OOB p, same row order as the sample
  best <- coords(roc(oof$obs, oof$p, levels = c(0, 1), direction = "<", quiet = TRUE),
                 "best", best.method = "youden",
                 ret = c("threshold", "sensitivity", "specificity"))[1, ]
  brier <- mean((thin_oob$p - thin_oob$obs)^2)
  
  cat("\n====", tag, "====\n")
  cat(sprintf("base rate %.4f%% | sample prevalence %.4f | n %d\n",
              100 * fit$base_rate, fit$sample_prev, fit$num.samples))
  cat(sprintf("OOB AUC thinned  %.3f (range %.3f-%.3f; %.0f rows at 1:1)\n",
              mean(thin_sum$oob_auc), min(thin_sum$oob_auc), max(thin_sum$oob_auc), mean(thin_sum$n_bal)))
  cat(sprintf("LOEO AUC pooled  %.3f | fold mean %.3f (SD %.3f, %d/%d) | n-weighted %.3f\n",
              auc_fixed(oof$obs, oof$p), mean(loeo$auc[ok]), sd(loeo$auc[ok]), sum(ok), nrow(loeo),
              weighted.mean(loeo$auc[ok], loeo$n[ok])))
  cat(sprintf("dropped folds hold %.1f%% of positives: %s\n",
              100 * sum(loeo$n_pos[!ok]) / sum(loeo$n_pos), paste(loeo$ecodistrict[!ok], collapse = ", ")))
  cat(sprintf("TSS max %.3f (thr %.3f) | at 0.5 %.3f | Brier %.3f (skill %.3f)\n",
              best$sensitivity + best$specificity - 1, best$threshold,
              tss_at(oof$obs, oof$p, 0.5), brier, 1 - brier / 0.25))
  print(head(read.csv(out("_importance.csv")), 10), row.names = FALSE)
  cat("thinned importance:\n")
  print(head(read.csv(out("_thin_importance.csv")), 10), row.names = FALSE)
}

# ==== conversion_native-annual ====
#   base rate 0.2983% | sample prevalence 0.5000 | n 299989
# OOB AUC thinned  0.925 (range 0.923-0.927; 17099 rows at 1:1)
# LOEO AUC pooled  0.906 | fold mean 0.829 (SD 0.092, 107/126) | n-weighted 0.839
# dropped folds hold 0.1% of positives: 728, 732, 733, 737, 741, 819, 824, 835, 836, 839, 840, 841, 843, 846, 847, 848, 851, 852, 854
# TSS max 0.644 (thr 0.198) | at 0.5 0.458 | Brier 0.112 (skill 0.551)
# variable importance
# crop 0.15754533
# unprot_grs 0.14375625
# prtot 0.13216479
# tx_mean 0.10961107
# dist_to_elevator 0.08612216
# psssm_sd 0.08430763
# gddgrow_5 0.07283189
# dist_to_popctr 0.06586325
# pop_density 0.06471381
# tame 0.06161139
# thinned importance:
#   variable importance
# crop 0.07930276
# unprot_grs 0.05096509
# prtot 0.03688907
# tx_mean 0.03657113
# psssm_sd 0.02785673
# gddgrow_5 0.01692109
# oc 0.01646072
# psssm_mean 0.01591482
# ph 0.01475893
# dist_to_elevator 0.01446740
# 
# ==== conversion_tame-annual ====
#   base rate 19.6799% | sample prevalence 0.5000 | n 299922
# OOB AUC thinned  0.763 (range 0.760-0.764; 33426 rows at 1:1)
# LOEO AUC pooled  0.746 | fold mean 0.724 (SD 0.069, 127/127) | n-weighted 0.728
# dropped folds hold 0.0% of positives: 
#   TSS max 0.364 (thr 0.443) | at 0.5 0.344 | Brier 0.199 (skill 0.202)
# variable importance
# crop 0.07241690
# bd 0.06058295
# unprot_grs 0.04433847
# dist_to_elevator 0.04292536
# tame 0.04246614
# prtot 0.03905890
# ph 0.03811073
# gddgrow_5 0.03743604
# oc 0.03584184
# tx_mean 0.03477050
# thinned importance:
#   variable  importance
# crop 0.030821389
# bd 0.019028626
# tame 0.014734651
# oc 0.014079078
# unprot_grs 0.012924555
# ph 0.009629520
# prtot 0.008556252
# dist_to_elevator 0.007810310
# gddgrow_5 0.007501365
# tx_mean 0.006878078
# 
# ==== conversion_tame-settlement ====
#   base rate 0.6729% | sample prevalence 0.5000 | n 299923
# OOB AUC thinned  0.891 (range 0.889-0.893; 18996 rows at 1:1)
# LOEO AUC pooled  0.898 | fold mean 0.858 (SD 0.041, 121/127) | n-weighted 0.855
# dropped folds hold 0.1% of positives: 777, 797, 799, 801, 802, 836
# TSS max 0.643 (thr 0.366) | at 0.5 0.591 | Brier 0.136 (skill 0.456)
# variable importance
# dist_to_road 0.11439187
# settl 0.10265488
# pop_density 0.08027642
# prtot 0.05868455
# unprot_grs 0.05506210
# dist_to_regctr 0.04708432
# tx_mean 0.04621235
# dist_to_elevator 0.04533073
# crop 0.04470260
# dist_to_popctr 0.04247546
# thinned importance:
#   variable  importance
# dist_to_road 0.075089379
# settl 0.040176583
# pop_density 0.025183229
# unprot_grs 0.020816121
# crop 0.017314595
# prtot 0.016842629
# tx_mean 0.013358687
# dist_to_elevator 0.009958840
# psssm_sd 0.008127812
# dist_to_popctr 0.007705789
# 
# ==== conversion_native-settlement ====
#   base rate 0.3660% | sample prevalence 0.5000 | n 299972
# OOB AUC thinned  0.919 (range 0.918-0.921; 17663 rows at 1:1)
# LOEO AUC pooled  0.940 | fold mean 0.885 (SD 0.053, 113/126) | n-weighted 0.891
# dropped folds hold 0.2% of positives: 728, 732, 737, 741, 835, 839, 840, 843, 847, 848, 851, 852, 854
# TSS max 0.740 (thr 0.300) | at 0.5 0.634 | Brier 0.117 (skill 0.532)
# variable importance
# settl 0.12733926
# unprot_grs 0.10152027
# pop_density 0.09674172
# dist_to_road 0.08952121
# dist_to_elevator 0.07386810
# tx_mean 0.06369607
# dist_to_rail 0.05952844
# prtot 0.05793263
# crop 0.05610987
# dist_to_popctr 0.05256734
# thinned importance:
#   variable importance
# dist_to_road 0.06124494
# settl 0.04821221
# unprot_grs 0.03310931
# crop 0.02800153
# pop_density 0.02273981
# prtot 0.01462050
# oc 0.01348728
# tx_mean 0.01346818
# dist_to_elevator 0.01247331
# dist_to_rail 0.01103389
# 
# ==== conversion_native-tame ====
#   base rate 0.5474% | sample prevalence 0.5000 | n 299972
# OOB AUC thinned  0.844 (range 0.844-0.845; 20936 rows at 1:1)
# LOEO AUC pooled  0.854 | fold mean 0.747 (SD 0.097, 108/126) | n-weighted 0.786
# dropped folds hold 0.0% of positives: 727, 728, 732, 733, 737, 741, 839, 840, 841, 843, 844, 846, 847, 848, 849, 851, 852, 854
# TSS max 0.562 (thr 0.319) | at 0.5 0.389 | Brier 0.162 (skill 0.353)
# variable importance
# crop 0.11722199
# unprot_grs 0.09327923
# prtot 0.09203712
# tame 0.08079609
# dist_to_elevator 0.08069250
# tx_mean 0.07310359
# psssm_mean 0.06169130
# gddgrow_5 0.06136884
# psssm_sd 0.06045722
# dist_to_popctr 0.05651305
# thinned importance:
#   variable  importance
# crop 0.046663678
# prtot 0.026395972
# unprot_grs 0.024602848
# psssm_sd 0.023263698
# tx_mean 0.018942766
# psssm_mean 0.018650858
# tame 0.015954345
# bd 0.014680278
# gddgrow_5 0.013138612
# net_return 0.009902713


# ---- importance heatmap, all transitions ----
# Thinned importance, not full-sample: ranger's permutation importance is OOB-based,
# so the full-sample version inherits the pseudoreplication that puts naive OOB
# AUC at 0.995 vs 0.925 thinned (native->annual).
# All thinned samples are rebalanced 1:1, so raw values share a scale across
# transitions; the pale tame->annual column reflects weaker skill (thinned OOB AUC
# 0.763 vs 0.84-0.93), not a units difference.
tags <- paste(transitions$from, transitions$to, sep = "-")
imp_all <- do.call(rbind, lapply(tags, \(t) cbind(transition = t,
                                                  read.csv(file.path(out_dir, paste0("rf_conversion_", t, "_thin_importance.csv"))))))
imp_all <- imp_all |>
  group_by(transition) |> mutate(rank = rank(-importance)) |> ungroup() |>
  mutate(variable   = reorder(variable, -rank, mean),   # best mean rank at top
         transition = factor(sub("-", " \u2192 ", transition), sub("-", " \u2192 ", tags)))

ggplot(imp_all, aes(transition, variable, fill = rank)) +
  geom_tile() + scale_fill_viridis_c(name = "Rank", direction =-1) +
  geom_text(aes(label = round(100 * importance, 1)), size = 2.2, colour = "white") +
  labs(x = NULL, y = NULL) + theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(out_dir, "rf_importance_heatmap.png"), width = 6, height = 8, dpi = 300)