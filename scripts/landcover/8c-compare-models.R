# Purpose:
#   Compares the RF (8-rf-conversion.R) against the GLM benchmark
#   (8b-glm-benchmark.R) for every transition where both have been run. Both
#   fit the same 7c rows and ecodistrict folds with the same fold gate, so
#   fold AUCs and out-of-fold predictions pair row for row.
#
#   Three questions:
#     transfer      : paired LOEO fold AUC (RF - GLM), and how many folds RF wins.
#                     Folds share ~99% of their training data, so wins are not
#                     independent: read them as consistency, not a test.
#     interpolation : RF thinned OOB vs GLM in-sample AUC. Different rows
#                     (thinned vs full sample), so the gap is approximate.
#     agreement     : rank correlation of the two models' out-of-fold
#                     predictions, and overlap of their top top_frac rows.
#                     A proxy for whether the choice would change CLUES
#                     allocation; computed on the balanced sample, not the
#                     landscape, so re-check on the mapped surfaces.
#
# Outputs:
#   - model_comparison_summary.csv : one row per transition
#   - model_comparison_folds.csv   : paired fold AUCs
#   - model_comparison_folds.png   : RF vs GLM fold AUC, one panel per transition

source(file.path("scripts/landcover", "0a-setup.R"))

out_dir  <- "output/landcover/model"
top_frac <- 0.10   # share of highest-ranked rows compared between models

auc_fixed <- function(obs, p) {
  as.numeric(auc(roc(obs, p, levels = c(0, 1), direction = "<", quiet = TRUE)))
}
tss_max <- function(obs, p) {
  b <- coords(roc(obs, p, levels = c(0, 1), direction = "<", quiet = TRUE), "best",
              best.method = "youden", ret = c("sensitivity", "specificity"))[1, ]
  b$sensitivity + b$specificity - 1
}

glm_sum <- read.csv(file.path(out_dir, "glm_benchmark_summary.csv"))

rows  <- list()
pairs <- list()

for (tr in glm_sum$transition) {
  tag  <- paste0("conversion_", tr)
  frf  <- function(s) file.path(out_dir, paste0("rf_",  tag, s))
  fglm <- function(s) file.path(out_dir, paste0("glm_", tag, s))
  if (!file.exists(frf("_loeo.csv"))) next
  
  # ---- transfer: paired folds ----
  lg <- read.csv(fglm("_loeo.csv"))
  lr <- read.csv(frf("_loeo.csv"))
  pr <- merge(lg[!is.na(lg$auc), c("ecodistrict", "n_pos", "auc")],
              lr[!is.na(lr$auc), c("ecodistrict", "auc")],
              by = "ecodistrict", suffixes = c("_glm", "_rf"))
  pr$diff <- pr$auc_rf - pr$auc_glm
  pairs[[tr]] <- cbind(transition = tr, pr)
  
  # ---- agreement: out-of-fold rows ----
  # Both scripts score the same 7c rows in the same order within each fold,
  # so a within-fold index aligns them. The obs check catches any drift.
  og <- read.csv(fglm("_oof.csv"))
  orf <- read.csv(frf("_oof.csv"))
  og$i  <- ave(seq_len(nrow(og)),  og$ecodistrict,  FUN = seq_along)
  orf$i <- ave(seq_len(nrow(orf)), orf$ecodistrict, FUN = seq_along)
  o <- merge(og, orf, by = c("ecodistrict", "i"), suffixes = c("_glm", "_rf"))
  stopifnot(all(o$obs_glm == o$obs_rf))
  
  top_g <- o$p_glm >= quantile(o$p_glm, 1 - top_frac)
  top_r <- o$p_rf  >= quantile(o$p_rf,  1 - top_frac)
  
  thin <- read.csv(frf("_thin.csv"))
  g    <- glm_sum[glm_sum$transition == tr, ]
  
  rows[[tr]] <- data.frame(
    transition    = tr,
    folds_paired  = nrow(pr),
    glm_auc_loeo  = mean(pr$auc_glm),
    rf_auc_loeo   = mean(pr$auc_rf),
    diff_mean     = mean(pr$diff),
    diff_median   = median(pr$diff),
    rf_wins       = sum(pr$diff > 0),
    glm_auc_interp = g$glm_auc_insample,
    rf_auc_interp  = mean(thin$oob_auc),
    glm_auc_pooled = auc_fixed(o$obs_glm, o$p_glm),
    rf_auc_pooled  = auc_fixed(o$obs_rf,  o$p_rf),
    glm_tss_max    = tss_max(o$obs_glm, o$p_glm),
    rf_tss_max     = tss_max(o$obs_rf,  o$p_rf),
    spearman       = cor(o$p_glm, o$p_rf, method = "spearman"),
    top_overlap    = sum(top_g & top_r) / sum(top_g),
    glm_top_prec   = mean(o$obs_glm[top_g]),
    rf_top_prec    = mean(o$obs_rf[top_r]))
  
  rm(lg, lr, og, orf, o)
}

summary_tbl <- do.call(rbind, rows)
folds_tbl   <- do.call(rbind, pairs)
write.csv(summary_tbl, file.path(out_dir, "model_comparison_summary.csv"), row.names = FALSE)
write.csv(folds_tbl,   file.path(out_dir, "model_comparison_folds.csv"),   row.names = FALSE)
print(summary_tbl, row.names = FALSE, digits = 3)
#        transition folds_paired glm_auc_loeo rf_auc_loeo diff_mean diff_median rf_wins glm_auc_interp rf_auc_interp glm_auc_pooled rf_auc_pooled glm_tss_max rf_tss_max spearman top_overlap glm_top_prec rf_top_prec
#     native-annual          107        0.800       0.829    0.0298      0.0178      70          0.895         0.925          0.880         0.906       0.596      0.644    0.883       0.470        0.942       0.989
#       tame-annual          127        0.701       0.724    0.0232      0.0207      97          0.720         0.763          0.713         0.746       0.309      0.364    0.811       0.509        0.778       0.828
# native-settlement          113        0.847       0.885    0.0374      0.0356     100          0.920         0.919          0.915         0.940       0.683      0.740    0.883       0.396        0.971       0.983
#   tame-settlement          121        0.826       0.858    0.0314      0.0286     106          0.883         0.891          0.878         0.898       0.595      0.643    0.882       0.288        0.963       0.950
#       native-tame          108        0.709       0.747    0.0380      0.0266      83          0.843         0.844          0.836         0.854       0.528      0.562    0.867       0.435        0.864       0.902

# Points above the line are folds the RF wins. Size by positives, so small
# folds with noisy AUCs read as such.
png(file.path(out_dir, "model_comparison_folds.png"),
    width = 1200 * length(pairs), height = 1300, res = 300)
par(mfrow = c(1, length(pairs)), mar = c(4, 4, 2, 1))
for (tr in names(pairs)) {
  p <- pairs[[tr]]
  plot(p$auc_glm, p$auc_rf, xlim = c(0.4, 1), ylim = c(0.4, 1), main = tr,
       xlab = "GLM fold AUC", ylab = "RF fold AUC",
       cex = 0.4 + 1.6 * sqrt(p$n_pos / max(p$n_pos)))
  abline(0, 1)
}
dev.off()