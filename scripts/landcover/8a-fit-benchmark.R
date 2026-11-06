# Purpose:
#   Logistic GLM benchmark for each conversion transition, scored with the same
#   leave-one-ecodistrict-out (LOEO) CV as the RF in 8b-fit-rf-recent-conversion.R. 
#   Main effects only. A large RF-over-GLM gap means the forest is picking up
#   nonlinearity or interactions. A small gap, which is plausible for
#   settlement, supports using the simpler model.
#
#   Fits on the shared per-origin samples from 7c-fitting-sample.R, the same
#   rows and ecodistrict folds as 8b-fit-rf-recent-conversion.R, so fold AUCs pair
#   directly. Each transition takes persistence plus one destination
#   (option A).
#
#   Balanced sampling moves only the GLM intercept. AUC is unchanged, and since
#   no surface is predicted here, no prior correction is needed.
#
#   Predictors are standardised on the full sample. In a main-effects GLM this
#   is an affine reparameterisation, so fold predictions and AUCs equal those
#   of the unscaled fit. The only effect is that coefficients become comparable
#   across predictors.
#
#   Tame- and annual-origin persistence rows are weighted by ecodistrict (7c column w, from 7a):
#   uncertain tame and annual pixels are excluded from training, and the weights restore each
#   district's map persistence. Native weights are all 1. Evaluation (AUC, TSS,
#   Brier) is unweighted. pop_density, dist_to_road and settl are log-transformed
#   (0b predictor_transforms); base_rate is from the hindcast map (7c).
#
# Landscape proportions: other (withheld) and native (vifstep) are the reference, so each
# proportion's coefficient is its effect relative to native grassland + other land.
#
# Outputs (tag = conversion_<from>-<to>):
#   - glm_<tag>_coef.csv       : full-sample coefficients, standardised predictors
#   - glm_<tag>_loeo.csv       : per-ecodistrict CV results, with separation flag
#   - glm_<tag>_oof.csv        : out-of-fold predictions
#   - glm_benchmark_summary.csv: one row per transition

source(file.path("scripts/landcover", "0a-setup.R"))
source(file.path("scripts/landcover", "0b-predictor-files.R"))

# Each transition: persistence + one destination from the 7c per-origin sample (option A).
transitions <- data.frame(from = c("native", "tame",   "native",     "tame",       "native", "annual", "annual"),
                          to   = c("annual", "annual", "settlement", "settlement", "tame",   "tame", "settlement"))
# ---- inputs ----
in_sa         <- file.path("./processeddatasets", "predictors", "30m")
out_dir       <- "output/landcover/model"
retained_file <- file.path(in_sa, "retained-covariates.rds")
min_class_n <- 25       # same fold gate as the RF, so evaluated folds line up

auc_fixed <- function(obs, p) {
  as.numeric(auc(roc(as.integer(as.character(obs)), p, levels = c(0, 1),
                     direction = "<", quiet = TRUE)))
}
tss_at <- function(obs, p, t) mean(p[obs == 1] >= t) + mean(p[obs == 0] < t) - 1

# Separation (fitted p of 0/1) and non-convergence both raise warnings. A fold
# that triggers either still yields an AUC, because the ranking holds, but its
# coefficients are not interpretable. The flag is recorded rather than the fold
# being dropped.
fit_glm <- function(dat) {
  flag <- FALSE
  m <- withCallingHandlers(glm(f, family = binomial, data = dat, weights = w,
                               control = glm.control(maxit = 100)),
                           warning = function(w) {
                             if (grepl("numerically 0 or 1|did not converge", conditionMessage(w))) {
                               flag <<- TRUE; invokeRestart("muffleWarning")
                             }
                             if (grepl("non-integer #successes", conditionMessage(w))) invokeRestart("muffleWarning")
                           })
  list(m = m, flag = flag)
}

# ---- predictors (shared across transitions) ----
final_set <- readRDS(retained_file)

f <- as.formula(paste("converted ~", paste(final_set, collapse = " + ")))

summary_rows <- list()

for (j in seq_len(nrow(transitions))) {
  from <- transitions$from[j]
  to   <- transitions$to[j]
  tag  <- sprintf("conversion_%s-%s", from, to)
  out  <- function(suffix) file.path(out_dir, paste0("glm_", tag, suffix))
  cat("\n==== ", from, " -> ", to, " ====\n", sep = "")
  
  # ---- sample ----
  # attributes are dropped by row subsetting, so read base_rate first
  d  <- readRDS(file.path(out_dir, paste0("train-sample_", from, ".rds")))
  br <- attr(d, "base_rate")[[to]]
  d  <- d[d$lulc_2020 %in% c(from, to), ]
  d$converted   <- factor(as.integer(d$lulc_2020 == to), levels = c(0, 1))
  d$ecodistrict <- droplevels(d$ecodistrict)
  
  # keep the scaling so the coefficients can be applied to the predictor rasters after
  # apply_transforms() (9b); center/scale are on the transformed scale
  sc <- scale(d[final_set]); d[final_set] <- sc
  cat("n =", nrow(d), "| positives =", sum(d$converted == 1),
      sprintf("| base rate %.4f%%", 100 * br),
      "| ecodistricts =", nlevels(d$ecodistrict), "\n")
  
  # ---- full-sample fit: coefficients only ----
  g <- fit_glm(d)
  # In-sample AUC on the full sample. With 31 coefficients against ~300k rows the
  # GLM can't memorise near-neighbour pixels the way the RF's naive OOB did (0.995),
  # so this stands in for the RF's thinned OOB (0.92 native) as the interpolation estimate.
  auc_in <- auc_fixed(d$converted, fitted(g$m))
  co <- as.data.frame(summary(g$m)$coefficients)
  co <- cbind(term = rownames(co), co)
  co <- co[order(-abs(co$Estimate)), ]
  co$center <- attr(sc, "scaled:center")[co$term]   # NA for the intercept
  co$scale  <- attr(sc, "scaled:scale")[co$term]
  write.csv(co, out("_coef.csv"), row.names = FALSE)
  if (g$flag) cat("full-sample fit flagged (separation / non-convergence)\n")
  
  # ---- LOEO ----
  # At ~300k rows x 31 predictors each glm takes seconds, so the whole loop
  # finishes in minutes. Unlike the RF, there is no resume logic.
  folds <- lapply(levels(d$ecodistrict), function(k) {
    te    <- d[d$ecodistrict == k, ]
    n_pos <- sum(te$converted == 1)
    row   <- data.frame(ecodistrict = k, n = nrow(te), n_pos = n_pos,
                        n_neg = nrow(te) - n_pos, auc = NA_real_, flagged = NA)
    if (min(n_pos, nrow(te) - n_pos) < min_class_n) return(list(row = row))
    gk <- fit_glm(d[d$ecodistrict != k, ])
    p  <- predict(gk$m, te, type = "response")
    row$auc <- auc_fixed(te$converted, p)
    row$flagged <- gk$flag
    list(row = row, oof = data.frame(ecodistrict = k, obs = as.integer(te$converted == 1), p = p))
  })
  loeo <- do.call(rbind, lapply(folds, `[[`, "row"))
  oof  <- do.call(rbind, lapply(folds, `[[`, "oof"))
  write.csv(loeo, out("_loeo.csv"), row.names = FALSE)
  write.csv(oof,  out("_oof.csv"),  row.names = FALSE)
  ok <- !is.na(loeo$auc)
  
  # TSS from LOEO out-of-fold predictions, as in the RF. Max TSS = Youden index;
  # one threshold tuned on ~300k rows adds negligible optimism.
  r_oof <- roc(oof$obs, oof$p, levels = c(0, 1), direction = "<", quiet = TRUE)
  best  <- coords(r_oof, "best", best.method = "youden",
                  ret = c("threshold", "sensitivity", "specificity"))
  # Brier from the full-sample fit at 1:1. This is the GLM counterpart of the
  # RF's thinned-OOB Brier; both are at sample prevalence, so neither is a
  # landscape-scale Brier.
  brier <- mean((fitted(g$m) - (d$converted == 1))^2)
  
  s <- data.frame(transition = paste(from, to, sep = "-"), n = nrow(d),
                  n_pos = sum(d$converted == 1), base_rate = br,
                  folds_eval = sum(ok), folds_total = nrow(loeo),
                  folds_flagged = sum(loeo$flagged, na.rm = TRUE),
                  glm_auc_mean = mean(loeo$auc[ok]), glm_auc_sd = sd(loeo$auc[ok]),
                  glm_auc_pooled = auc_fixed(oof$obs, oof$p),
                  glm_auc_insample = auc_in,
                  tss_max = best$sensitivity + best$specificity - 1,
                  tss_thr = best$threshold,
                  tss_05  = tss_at(oof$obs, oof$p, 0.5),
                  sens = best$sensitivity, spec = best$specificity,
                  brier = brier, brier_skill = 1 - brier / 0.25)
  
  cat(sprintf("GLM LOEO AUC mean %.3f (SD %.3f, %d/%d folds, %d flagged) | pooled %.3f | TSS max %.3f (thr %.2f) | Brier %.3f\n",
              s$glm_auc_mean, s$glm_auc_sd, s$folds_eval, s$folds_total,
              s$folds_flagged, s$glm_auc_pooled,
              s$tss_max, s$tss_thr, s$brier))
  print(head(co[co$term != "(Intercept)", c("term", "Estimate", "Std. Error")], 10), row.names = FALSE)
  
  summary_rows[[j]] <- s
  rm(d, g, co, folds, loeo, oof)
  gc()
}

# ==== native -> annual ====
#   n = 299989 | positives = 149997 | base rate 0.2983% | ecodistricts = 126 
# GLM LOEO AUC mean 0.800 (SD 0.098, 107/126 folds, 0 flagged) | pooled 0.881 | TSS max 0.597 (thr 0.44) | Brier 0.130
# term   Estimate  Std. Error
# crop  1.2439632 0.013768209
# bd  0.7133475 0.010717721
# tame  0.5670269 0.005616704
# fors  0.4997861 0.009205833
# gddgrow_5 -0.4890221 0.012293690
# prtot  0.4490524 0.012850613
# oc -0.4131260 0.008377588
# unprot_grs -0.3452263 0.015597421
# dist_to_popctr -0.3381885 0.008904609
# tri -0.2210934 0.006482693
# 
# ==== tame -> annual ====
#   n = 299945 | positives = 149972 | base rate 14.3479% | ecodistricts = 127 
# GLM LOEO AUC mean 0.763 (SD 0.082, 126/127 folds, 0 flagged) | pooled 0.764 | TSS max 0.391 (thr 0.46) | Brier 0.194
# term   Estimate  Std. Error
# bd  1.2773353 0.009121634
# crop  0.8431120 0.011145907
# gddgrow_5 -0.5496440 0.010164857
# prtot  0.4373157 0.009202740
# oc -0.2606381 0.005830614
# clay -0.2101517 0.006500832
# tx_mean -0.1950483 0.009663392
# dist_to_elevator -0.1829999 0.006856695
# tri -0.1799535 0.005431133
# psssm_mean  0.1798342 0.006148070
# 
# ==== native -> settlement ====
#   n = 299972 | positives = 149980 | base rate 0.3660% | ecodistricts = 126 
# GLM LOEO AUC mean 0.862 (SD 0.055, 113/126 folds, 0 flagged) | pooled 0.924 | TSS max 0.704 (thr 0.47) | Brier 0.105
# term   Estimate  Std. Error
# settl  1.0972987 0.010806286
# oc -0.7594198 0.009151363
# dist_to_road -0.7007609 0.006918792
# unprot_grs -0.6740218 0.014424401
# gddgrow_5 -0.4649279 0.014623393
# pop_density  0.4505396 0.014764240
# tame  0.4232751 0.006298606
# bd  0.4044544 0.011995745
# prtot  0.3966890 0.012728716
# dist_to_elevator -0.3784598 0.010211420
# 
# ==== tame -> settlement ====
#   n = 299932 | positives = 149959 | base rate 0.6313% | ecodistricts = 127 
# GLM LOEO AUC mean 0.842 (SD 0.048, 121/127 folds, 0 flagged) | pooled 0.888 | TSS max 0.618 (thr 0.44) | Brier 0.133
# term   Estimate  Std. Error
# dist_to_road -1.1065632 0.007291123
# settl  1.0277969 0.011327808
# crop  0.6684018 0.013560503
# bd  0.5639196 0.010147569
# tx_mean -0.4435877 0.011427851
# prtot  0.4209517 0.011241333
# pop_density  0.4054190 0.012288626
# psssm_mean  0.3544201 0.008636928
# drainage -0.2798370 0.007239805
# dist_to_elevator  0.2715793 0.008607367
# 
# ==== native -> tame ====
#   n = 299972 | positives = 149980 | base rate 0.5474% | ecodistricts = 126 
# GLM LOEO AUC mean 0.709 (SD 0.117, 108/126 folds, 0 flagged) | pooled 0.836 | TSS max 0.528 (thr 0.43) | Brier 0.160
# term   Estimate  Std. Error
# tame  0.7103302 0.005477334
# gddgrow_5 -0.6319766 0.011399419
# crop  0.6251065 0.010028665
# unprot_grs -0.5188157 0.011057157
# oc -0.2855891 0.007864501
# tx_mean  0.2787518 0.010865011
# psssm_sd  0.2495440 0.007738356
# fors  0.1586160 0.007180456
# slope -0.1405979 0.005801928
# psssm_mean  0.1340427 0.007218435
# 
# ==== annual -> tame ====
#   n = 299961 | positives = 149968 | base rate 1.5829% | ecodistricts = 127 
# GLM LOEO AUC mean 0.809 (SD 0.067, 125/127 folds, 0 flagged) | pooled 0.828 | TSS max 0.499 (thr 0.46) | Brier 0.165
# term   Estimate  Std. Error
# crop -1.1201894 0.013940721
# bd -0.9581265 0.010067578
# oc  0.6883493 0.007179214
# prtot -0.6218839 0.009835278
# gddgrow_5  0.3830812 0.011165929
# tame  0.2473950 0.008560154
# tx_mean  0.2049246 0.010567603
# slope  0.1925180 0.006529282
# clay  0.1723989 0.007130728
# unprot_grs -0.1668303 0.014402818
#
# ==== annual -> settlement ====
#   n = 299967 | positives = 149974 | base rate 0.2066% | ecodistricts = 127 
# GLM LOEO AUC mean 0.857 (SD 0.059, 119/127 folds, 0 flagged) | pooled 0.885 | TSS max 0.613 (thr 0.46) | Brier 0.135
# term   Estimate  Std. Error
# dist_to_road -1.3216828 0.007707051
# settl  1.1342963 0.011502111
# oc  0.8199542 0.007971219
# tx_mean -0.7399509 0.012084995
# prtot -0.4916318 0.011277893
# gddgrow_5  0.3582387 0.014126459
# bd -0.3124840 0.011479226
# tame  0.2450997 0.011909143
# unprot_grs -0.2408642 0.016187549
# tpi  0.1795029 0.005720818

# ---- GLM coefficient heatmap, all transitions ----
# Standardised coefficients (log-odds per SD), signed. Conditional linear effects,
# not comparable to RF permutation importance: correlated predictors (e.g. the
# crop/tame/unprot_grs focal proportions) can split effects or flip sign.
# No SEs/stars: ~300k spatially autocorrelated rows make every term "significant".
tags <- paste(transitions$from, transitions$to, sep = "-")
coef_all <- do.call(rbind, lapply(tags, \(t) cbind(transition = t,
                                                   read.csv(file.path(out_dir, paste0("glm_conversion_", t, "_coef.csv")))))) |>
  filter(term != "(Intercept)") |>
  group_by(transition) |> mutate(rank = rank(-abs(Estimate))) |> ungroup() |>
  mutate(term       = reorder(term, -rank, mean),
         transition = factor(sub("-", " \u2192 ", transition), sub("-", " \u2192 ", tags)))

ggplot(coef_all, aes(transition, term, fill = Estimate)) +
  geom_tile() + geom_text(aes(label = round(Estimate, 2)), size = 2.2) +
  scale_fill_gradient2(name = "Coefficient", low = "#2166ac", high = "#b2182b") +
  labs(x = NULL, y = NULL) + theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(out_dir, "glm_coef_heatmap.png"), width = 6, height = 8, dpi = 300)



summary_tbl <- do.call(rbind, summary_rows)
write.csv(summary_tbl, file.path(out_dir, "glm_benchmark_summary.csv"), row.names = FALSE)
print(summary_tbl, row.names = FALSE)
# transition             n  n_pos   base_rate folds_eval folds_total folds_flagged glm_auc_mean glm_auc_sd glm_auc_pooled glm_auc_insample   tss_max   tss_thr    tss_05      sens      spec
# native-annual     299989 149997 0.002982804        107         126             0    0.7997535 0.09791644      0.8811522        0.8961714 0.5970444 0.4369881 0.5913294 0.8202624 0.7767820
# tame-annual       299945 149972 0.143478836        126         127             0    0.7626988 0.08196138      0.7643367        0.7735999 0.3914504 0.4581935 0.3896624 0.7506171 0.6408333
# native-settlement 299972 149980 0.003659692        113         126             0    0.8621418 0.05455598      0.9238700        0.9281379 0.7042187 0.4712492 0.7031759 0.8597819 0.8444368
# tame-settlement   299932 149959 0.006313105        121         127             0    0.8422000 0.04832796      0.8882023        0.8920468 0.6177374 0.4438936 0.6149238 0.8315118 0.7862256
# native-tame       299972 149980 0.005473886        108         126             0    0.7085296 0.11651521      0.8359198        0.8430036 0.5279253 0.4334053 0.5224151 0.8237370 0.7041883
# annual-tame       299961 149968 0.015828752        125         127             0    0.8090089 0.06695324      0.8280283        0.8371152 0.4993181 0.4647357 0.4967551 0.7593654 0.7399528
# annual-settlement 299967 149974 0.002065815        119         127             0    0.8571685 0.05907506      0.8846795        0.8886171 0.6127609 0.4550543 0.6092948 0.8108739 0.8018870
# brier brier_skill
# 0.1302141   0.4791434
# 0.1935816   0.2256737
# 0.1054411   0.5782356
# 0.1329270   0.4682920
# 0.1603523   0.3585906
# 0.1652626   0.3389498
# 0.1351368   0.4594528