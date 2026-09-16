# Pairwise correlation + VIF screen to flag collinear predictors before the RF.
# RF tolerates collinearity for prediction, but it distorts variable importance
# and effect plots, so redundant predictors are screened out.
source(file.path("scripts/landcover", "0a-setup.R"))
source(file.path("scripts/landcover", "0b-predictor-files.R"))

# ---- inputs ----
out_dir  <- "output/landcover/collinearity"
out_sa   <- file.path("./processeddatasets", "predictors", paste0(res_m, "m"))
sample_n   <- 20000      # target complete (non-NA) sample size
cor_method <- "spearman" # "spearman" (monotonic) or "pearson" (linear)
cor_th     <- 0.7        # |r| at/above this flags a pair as collinear
vif_th     <- 5          # VIF at/above this flags multicollinearity (5 = stricter)

sa_r <- rast(file.path(sa_filepath, sa_raster_filename))
ecoregions <- st_read(tagged_path("./processeddatasets", ecoregions_filepath, ecoregions_filename, "modelling-extent"))
ecoregion_r <- rast(file.path("./processeddatasets", ecoregions_filepath, paste0(tools::file_path_sans_ext(ecoregions_filename), "_modelling-extent.tif")))
msk <- rast(file.path(out_sa, "predictor-valid-mask_study-area.tif"))

# baseline (2010) grassland: native + tame. Screening samples the training
# population, not the whole study area -- remaining grassland skews toward
# marginal land, so covariate relationships within it can differ from those
# across all land cover. 
grass <- rast(file.path("./processeddatasets", aafc_filepath, "binary",
                        "AAFC2010_ACI2009-2012_grass.tif"))
grass <- crop(grass, msk, snap = "out")
msk   <- subst(msk & (grass > 0), FALSE, NA)   # 0-cells must be NA for spatSample

# ---- load rasters ----
files <- predictor_files

missing <- files[!file.exists(files)]
if (length(missing)) {
  stop("missing predictor layers:\n  ", paste(names(missing), missing, sep = ": ", collapse = "\n  "))
}

sa_file <- function(nm) file.path(out_sa, paste0(nm, "_study-area.tif"))
r <- rast(sapply(names(files), sa_file))
names(r) <- names(files)

# ---- sample values ----
# Sample xy, not cell numbers: msk is on the study-area grid while the
# predictor files are on the modelling-extent grid, and cell indices don't
# transfer between grids. xy is grid-independent (and exact here, since the
# SA grid is a crop of the ME grid).
#
# Oversampling is done manually rather than via na.rm = TRUE: only ~15% of
# cells are eligible (grassland within the valid mask), and at that valid
# fraction spatSample's internal oversampling under-delivers and warns rather
# than iterating up to the requested n. So draw 20x with na.rm = FALSE, keep
# the rows where the mask value (column 3, after x and y) is non-NA, then
# randomly subsample to sample_n. The subsample must be random, not head(),
# because spatSample returns cells in raster order.
set.seed(42)
pts <- spatSample(msk, sample_n * 20, method = "random",
                  xy = TRUE, values = TRUE, na.rm = FALSE)
pts <- as.matrix(pts[!is.na(pts[[3]]), c("x", "y")])

vals <- cbind(pts, terra::extract(r, pts))
vals$ID <- NULL
vals <- na.omit(vals)
vals <- vals[sample(nrow(vals), min(sample_n, nrow(vals))), ]

saveRDS(vals, file.path(out_sa, "covariate-sample_training.rds"))

# closure sets: two exact linear dependencies in the predictor set.
#   texture   clay + sand + silt ~ 100
#   focal     native + tame + crop + dev + fors + other = 1 at every pixel
# Both are dropped a priori, not by threshold - with all members present the
# correlation matrix is singular and VIF is meaningless. One member of each set 
# is withheld as the reference: clay is retained over sand/silt as the more 
# interpretable; other is withheld per 4d, where it is computed but not 
# exported. The retained five proportions still sum to (1 - other) and remain 
# strongly collinear where other is small;vifstep arbitrates that.
#
# categorical, excluded from the screen but retained in the model:
#   grass_origin  binary native/tame. vif() needs numeric, and this is the
#     cohort split (0.22% vs 17.4% conversion), in by design rather than by
#     threshold. Expect it to lead importance and compress the rest.
#   ecoregion     25-level factor. Dummy-coding gives 25 columns summing to 1,
#     the same singularity as above, so VIF is undefined for it. Collinearity
#     with climate and soil is real but unquantifiable here (ecoregions are
#     delineated on those variables) - it nests above the LOEO folds, so the
#     with/without comparison is on LOEO AUC, not OOB.
drop_cols <- c("x", "y", "silt", "sand", "other", "grass_origin", "ecoregion")
vals_cov <- vals[, setdiff(names(vals), drop_cols)]


# ---- pairwise correlation ----
cmat <- cor(vals_cov, method = cor_method, use = "pairwise.complete.obs")
write.csv(round(cmat, 3), file.path(out_dir, "correlation_matrix.csv"))
pairs <- as.data.frame(as.table(cmat)) |>
  filter(Var1 != Var2, abs(Freq) >= cor_th) |>
  mutate(pair = paste(pmin(as.character(Var1), as.character(Var2)),
                      pmax(as.character(Var1), as.character(Var2)), sep = " - ")) |>
  distinct(pair, .keep_all = TRUE) |>
  transmute(pair, r = round(Freq, 3)) |>
  arrange(desc(abs(r)))

# ---- VIF screen ----
# vif_all: VIF for every variable at once, no selection. Each is regressed on
# all the others; VIF = 1/(1-R^2), so 10 means 90% of that variable is
# explained by the rest. Diagnostic only -- shows where the collinearity is.
vif_all  <- vif(vals_cov)

# vif_step: iterative removal. Computes all VIFs, drops the single worst
# variable, recomputes, repeats until every VIF < th. Because it drops one at
# a time, the surviving set depends on removal order - it gives *a* set below
# the threshold, not the only one.
#
# keep = protects variables from removal regardless of VIF; vifstep drops
# others to accommodate them.
#   climate (gddgrow_5, prtot, tx_mean): VIFs of 36/12/30, so vifstep would
#     strip the climate block. Retained because the CLUES scenario is SSP2-4.5
#     - without climate terms the model cannot respond to the scenario at all.
#   clay: the texture closure reference. sand and silt are already dropped
#     above, so if vifstep also removed clay the predictor set would carry no
#     texture information. keep pins the one member deliberately chosen.
#
# Input is vals_cov: all predictors except the two closure references
# (sand/silt, other). tn_mean and nrcdd are left in so vifstep decides them
# here rather than being excluded by hand upstream.
vif_step <- vifstep(vals_cov, th = vif_th, keep = c("gddgrow_5", "prtot", "tx_mean", "clay"))

# ---- correlation heatmap ----
hm <- as.data.frame(as.table(cmat)) |>
  ggplot(aes(Var1, Var2, fill = Freq)) +
  geom_tile() +
#  geom_text(aes(label = round(Freq, 2)), size = 3) +
  geom_text(aes(label = ifelse(abs(Freq) >= cor_th, round(Freq, 2), "")), size = 2.5) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-1, 1), name = cor_method) +
  labs(x = NULL, y = NULL, title = "Predictor correlation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#ggsave(file.path(out_dir, "correlation_heatmap.png"), hm, width = 7, height = 6, dpi = 150)
ggsave(file.path(out_dir, "correlation_heatmap.png"), hm, width = 11, height = 10, dpi = 150)

# ---- report ----
cat("Correlation method:", cor_method, "\n\n")
# Correlation method spearman 
cat("Pairs with |r| >=", cor_th, ":\n"); print(pairs, row.names = FALSE)
# Pairs with |r| >= 0.7 :
#   pair      r
# kp1500 - kp33  0.934
# tn_mean - tx_mean  0.828
# crop - native -0.796
# nrcdd - prtot -0.765
# ph - tx_mean  0.736
# prtot - tsm_growing  0.730
# bd - nrcdd  0.715
# bd - ph  0.709
cat("\nVIF (all variables):\n"); print(transform(vif_all, VIF = round(VIF, 2)))
# VIF (all variables):
#   Variables   VIF
# 1         prtot 11.61
# 2       tn_mean 12.02
# 3       tx_mean 29.79
# 4     gddgrow_5 35.62
# 5         nrcdd  8.29
# 6            oc  2.24
# 7            bd  4.52
# 8          clay  2.08
# 9            ph  3.62
# 10     drainage  3.13
# 11     watertbl  2.59
# 12          kp0  2.58
# 13         kp33 19.46
# 14       kp1500 17.86
# 15       ldepth  1.73
# 16   ldepth_max  1.42
# 17   psssm_mean  2.57
# 18     psssm_sd  2.73
# 19 tsm_planting  9.01
# 20  tsm_growing  6.79
# 21         elev 28.98
# 22          tri  1.89
# 23          tpi  1.20
# 24          srr  1.17
# 25    northness  1.02
# 26     eastness  1.02
# 27        slope  1.69
# 28   unprot_grs  2.56
# 29       native 25.88
# 30         tame  8.67
# 31         crop 11.03
# 32          dev  1.37
# 33         fors  3.58
# 34   net_return  2.88
# 35   d_elevator  3.01
# 36       d_road  1.25
# 37       d_rail  2.05
# 38     d_popctr  2.69
# 39     d_regctr  2.12

cat("\nvifstep (th =", vif_th, ") -- variables flagged for exclusion:\n")
print(vif_step@excluded)
# [1] "elev"         "native"       "kp33"         "tn_mean"      "nrcdd"        "tsm_planting"
cat("\nRetained set:\n"); print(vif_step@results$Variables)
# Retained set:
#   [1] "prtot"       "tx_mean"     "gddgrow_5"   "oc"          "bd"          "clay"        "ph"          "drainage"   
# [9] "watertbl"    "kp0"         "kp1500"      "ldepth"      "ldepth_max"  "psssm_mean"  "psssm_sd"    "tsm_growing"
# [17] "tri"         "tpi"         "srr"         "northness"   "eastness"    "slope"       "unprot_grs"  "tame"       
# [25] "crop"        "dev"         "fors"        "net_return"  "d_elevator"  "d_road"      "d_rail"      "d_popctr"   
# [33] "d_regctr"       

# vifstep output; manual drops are applied after the diagnostics below, which
# supply the evidence for them
vif_set <- vif_step@results$Variables
vif(vals_cov[, vif_set])
# Variables      VIF
# 1        prtot 6.060397
# 2      tx_mean 6.232994
# 3    gddgrow_5 5.365917
# 4           oc 2.124140
# 5           bd 4.109704
# 6         clay 1.901151
# 7           ph 3.289288
# 8     drainage 2.947956
# 9     watertbl 2.506332
# 10         kp0 1.716421
# 11      kp1500 2.134001
# 12      ldepth 1.585835
# 13  ldepth_max 1.366774
# 14  psssm_mean 2.156723
# 15    psssm_sd 2.289552
# 16 tsm_growing 2.675862
# 17         tri 1.791996
# 18         tpi 1.140516
# 19         srr 1.122572
# 20   northness 1.024043
# 21    eastness 1.018543
# 22       slope 1.688722
# 23  unprot_grs 2.354606
# 24        tame 1.523713
# 25        crop 2.353267
# 26         dev 1.166490
# 27        fors 2.176893
# 28  net_return 2.509276
# 29  d_elevator 2.950418
# 30      d_road 1.228350
# 31      d_rail 1.996731
# 32    d_popctr 2.646275
# 33    d_regctr 2.026562
length(vif_set)
# [1] 33

# ---- variance explained by ecoregion --------------------------------------
# Not a redundancy check -- ecoregion is included as a regional catch-all, not
# for independent information. But a covariate that is close to constant within
# ecoregions will split its importance with the ecoregion factor, the same way
# two correlated continuous predictors do. Worth knowing before reading
# variable importance.

ecoregion_r <- mask(crop(ecoregion_r, sa_r), sa_r)
names(ecoregion_r) <- "ecoregion"

d <- vals[, c("x", "y", vif_set)]

d$ecoregion <- factor(terra::extract(ecoregion_r, d[, c("x", "y")])[, 2])
d <- d[!is.na(d$ecoregion), ]
cat("n =", nrow(d), "across", nlevels(droplevels(d$ecoregion)), "ecoregions\n")
# n = 20000 across 9 ecoregions
d$ecoregion <- droplevels(d$ecoregion)

r2_eco <- sapply(setdiff(names(d), c("x", "y", "ecoregion")), function(v) {
  summary(lm(d[[v]] ~ d$ecoregion))$r.squared
})
print(round(sort(r2_eco, decreasing = TRUE), 3))
# prtot     tx_mean          ph  d_elevator          bd  psssm_mean   gddgrow_5 tsm_growing        fors    watertbl 
# 0.708       0.580       0.563       0.537       0.519       0.383       0.382       0.354       0.325       0.300 
# d_popctr    psssm_sd          oc  net_return    drainage  unprot_grs      d_rail        crop         tri        tame 
# 0.289       0.222       0.208       0.189       0.184       0.158       0.125       0.124       0.122       0.110 
# kp0        clay    d_regctr      kp1500       slope      d_road  ldepth_max      ldepth         dev         srr 
# 0.096       0.085       0.068       0.061       0.055       0.051       0.035       0.029       0.028       0.009 
# eastness         tpi   northness 
# 0.004       0.001       0.001 

ecodistricts <- st_read(tagged_path("./processeddatasets", ecodistricts_filepath,
                                    ecodistricts_filename, "modelling-extent"))
ecodistricts <- st_transform(ecodistricts, crs(sa_r))
ecodistrict_r <- rasterize(vect(ecodistricts), sa_r, field = "ECODISTRIC")
ecodistrict_r <- mask(ecodistrict_r, sa_r)

d$ecodistrict <- factor(terra::extract(ecodistrict_r, d[, c("x", "y")])[, 2])
d <- d[!is.na(d$ecodistrict), ]
d$ecodistrict <- droplevels(d$ecodistrict)

r2_dist <- sapply(setdiff(names(d), c("x", "y", "ecoregion", "ecodistrict")), function(v) {
  summary(lm(d[[v]] ~ d$ecodistrict))$r.squared
})
print(round(sort(r2_dist, decreasing = TRUE), 3))
# prtot     tx_mean   gddgrow_5  d_elevator    psssm_sd  psssm_mean          ph          bd  net_return tsm_growing 
# 0.933       0.927       0.884       0.873       0.872       0.828       0.783       0.761       0.756       0.752 
# d_regctr        clay          oc    d_popctr        fors      d_rail      kp1500    watertbl  unprot_grs    drainage 
# 0.710       0.690       0.650       0.647       0.493       0.474       0.452       0.448       0.418       0.414 
# crop         kp0  ldepth_max         tri        tame       slope         dev      d_road      ldepth         srr 
# 0.377       0.372       0.356       0.329       0.314       0.202       0.158       0.137       0.124       0.035 
# northness    eastness         tpi 
# 0.026       0.017       0.012 

# ---- variance explained by province ----------------------------------------
# DSS attributes come from independently-conducted provincial surveys with
# different depth conventions and mapping intensity, and the maps show visible
# steps at the AB/SK and SK/MB borders. A soil covariate that mostly encodes
# province acts as another regional stratifier (same mechanism as the OOB/LOEO
# gap) and will imprint provincial seams on prediction maps. Flag: high
# province R2 AND low within-province variation (sd ratio = sqrt(1 - R2), so
# per-province SDs below identify where the signal collapses).
prov <- geodata::gadm("CAN", level = 1, path = file.path(out_sa, "gadm"))
prov <- project(prov[prov$NAME_1 %in% c("Alberta", "Saskatchewan", "Manitoba"), ], crs(sa_r))

d$province <- factor(terra::extract(prov, d[, c("x", "y")])$NAME_1)
d <- d[!is.na(d$province), ]

r2_prov <- sapply(setdiff(names(d), c("x", "y", "ecoregion", "ecodistrict", "province")), function(v) {
  summary(lm(d[[v]] ~ d$province))$r.squared
})
print(round(sort(r2_prov, decreasing = TRUE), 3))
# net_return    watertbl    psssm_sd       prtot    drainage  ldepth_max        clay   gddgrow_5  unprot_grs     tx_mean 
# 0.541       0.354       0.312       0.306       0.255       0.239       0.201       0.195       0.195       0.133 
# fors      d_rail         kp0  psssm_mean  d_elevator      ldepth tsm_growing        tame        crop         tri 
# 0.129       0.123       0.114       0.113       0.082       0.063       0.061       0.044       0.042       0.042 
# kp1500      d_road         dev    d_popctr       slope          ph          oc    d_regctr          bd         srr 
# 0.041       0.034       0.034       0.033       0.024       0.023       0.022       0.020       0.013       0.002 
# tpi    eastness   northness 
# 0.002       0.000       0.000 
# per-province mean and SD for the soil variables: a province with SD near 0
# (e.g. a depth convention capping at 100 cm) is acting as a province code
soil_vars <- intersect(c("oc", "bd", "clay", "sand", "ph", "drainage", "watertbl",
                         "kp0", "kp33", "kp1500", "ldepth", "ldepth_max", "slope"),
                       names(d))
prov_stats <- do.call(rbind, lapply(soil_vars, function(v) {
  data.frame(var  = v,
             prov = levels(d$province),
             mean = round(tapply(d[[v]], d$province, mean), 2),
             sd   = round(tapply(d[[v]], d$province, sd), 2),
             row.names = NULL)
}))
print(prov_stats)
# var         prov   mean    sd
# 1          oc      Alberta  17.86  3.49
# 2          oc     Manitoba  19.96  2.57
# 3          oc Saskatchewan  18.15  2.54
# 4          bd      Alberta   1.32  0.07
# 5          bd     Manitoba   1.32  0.05
# 6          bd Saskatchewan   1.33  0.05
# 7        clay      Alberta  19.88  3.91
# 8        clay     Manitoba  25.14  5.36
# 9        clay Saskatchewan  23.04  2.68
# 10         ph      Alberta   7.09  0.38
# 11         ph     Manitoba   6.87  0.38
# 12         ph Saskatchewan   7.08  0.39
# 13   drainage      Alberta   3.18  0.59
# 14   drainage     Manitoba   4.05  1.11
# 15   drainage Saskatchewan   2.74  0.74
# 16   watertbl      Alberta   0.05  0.27
# 17   watertbl     Manitoba   1.10  1.20
# 18   watertbl Saskatchewan   0.05  0.26
# 19        kp0      Alberta  52.27  3.31
# 20        kp0     Manitoba  50.15  9.24
# 21        kp0 Saskatchewan  49.35  2.54
# 22     kp1500      Alberta  16.97  4.62
# 23     kp1500     Manitoba  12.95  4.73
# 24     kp1500 Saskatchewan  16.89  5.49
# 25     ldepth      Alberta  19.97 20.69
# 26     ldepth     Manitoba  16.66 13.23
# 27     ldepth Saskatchewan  12.17  2.25
# 28 ldepth_max      Alberta 112.05 14.96
# 29 ldepth_max     Manitoba 101.53  5.30
# 30 ldepth_max Saskatchewan 100.02  0.39
# 31      slope      Alberta   5.46  5.72
# 32      slope     Manitoba   4.35  9.33
# 33      slope Saskatchewan   6.43  4.97
write.csv(prov_stats, file.path(out_dir, "province_stats_soil.csv"), row.names = FALSE)

# ---- where are the DSS NAs? ------------------------------------------------
# The DSS layers are the binding NA constraint (~35.5M cells on kp0/kp1500).
# Dropping a layer only recovers study area if its NAs are attribute gaps in
# surveyed polygons; NAs where DSS has no polygon coverage at all are lost
# either way. Separate draw from d: that sample is all-valid by construction,
# so it contains no NA rows to count. Sampled over the whole study area
# (not the valid mask), NAs retained.
na_n <- 500000   # ~0.1% precision on a proportion; cheap relative to a full pass

set.seed(43)
na_pts <- spatSample(sa_r, na_n * 3, method = "random",
                     xy = TRUE, values = TRUE, na.rm = FALSE)
na_pts <- as.matrix(na_pts[!is.na(na_pts[[3]]), c("x", "y")])
na_pts <- na_pts[sample(nrow(na_pts), min(na_n, nrow(na_pts))), ]

na_d <- terra::extract(r, na_pts)
na_d$ID <- NULL
na_d$province <- factor(terra::extract(prov, na_pts)$NAME_1)
na_d <- na_d[!is.na(na_d$province), ]

dss_vars <- intersect(c("drainage", "watertbl", "kp0", "kp33", "kp1500",
                        "ldepth", "ldepth_max", "slope"), names(na_d))

# % NA per DSS layer within each province
na_by_prov <- sapply(dss_vars, function(v) {
  round(100 * tapply(is.na(na_d[[v]]), na_d$province, mean), 1)
})
cat("\n% NA by province:\n"); print(na_by_prov)
# % NA by province:
#   drainage watertbl kp0 kp33 kp1500 ldepth ldepth_max slope
# Alberta           2.2      1.1 2.2  2.3    2.2    0.1        0.1   0.1
# Manitoba          6.2      6.2 6.8  6.8    6.8    6.8        6.8   6.2
# Saskatchewan      1.7      0.0 9.8  9.8    9.8    9.8        9.8   1.7

# is a layer's NA set exclusive, or shared with the other DSS layers?
# exclusive NAs are what dropping that layer would actually recover
dss_na <- sapply(dss_vars, function(v) is.na(na_d[[v]]))
excl <- sapply(dss_vars, function(v) {
  others <- rowSums(dss_na[, setdiff(dss_vars, v), drop = FALSE]) > 0
  sum(dss_na[, v] & !others)
})
cat("\nDSS NA overlap (n =", nrow(na_d), "sampled):\n")
# DSS NA overlap (n = 499428 sampled):
print(data.frame(var        = dss_vars,
                 pct_na     = round(100 * colMeans(dss_na), 1),
                 pct_na_excl = round(100 * excl / nrow(na_d), 2),
                 row.names  = NULL))
# var pct_na pct_na_excl
# 1   drainage    2.5        0.00
# 2   watertbl    1.3        0.00
# 3        kp0    6.8        0.00
# 4       kp33    6.8        0.05
# 5     kp1500    6.8        0.00
# 6     ldepth    6.1        0.00
# 7 ldepth_max    6.1        0.00
# 8      slope    1.8        0.00
cat("\nany DSS layer NA:", round(100 * mean(rowSums(dss_na) > 0), 1), "%\n")
# any DSS layer NA: 6.8 %
write.csv(as.data.frame(na_by_prov), file.path(out_dir, "province_na_dss.csv"))

# ---- manual drops: jurisdictional artifacts --------------------------------
# Applied here, after the diagnostics above, because those tables are the
# evidence. DSS attributes come from independently-conducted provincial surveys
# with different depth conventions and mapping intensity, and the maps show
# steps at the AB/SK and SK/MB borders. Province R2 vs ecodistrict R2 separates
# real ecology from survey convention -- province boundaries have no ecological
# meaning, so a high ratio means the variable is partly encoding jurisdiction
# (values from the run of 2026-08; re-check against the printout above if the
# sample or source layers change):
#
#              province   ecodistrict   ratio
#   watertbl      0.289       0.371      0.78   <- dropped
#   ldepth_max    0.245       0.383      0.64   <- dropped
#   ldepth        0.062       0.121      0.51      kept, but see SD note below
#   drainage      0.185       0.377      0.49      kept (mechanistically relevant)
#   kp0           0.132       0.357      0.37      kept
#   kp1500        0.033       0.463      0.07      kept (tracks ecology, not province)
#   ph            0.017       0.790      0.02      kept (interpolated, province-blind)
#
# The interpolated national products (ph, bd, oc; all < 0.03) are essentially
# province-blind, which is what shows the DSS pattern is survey-derived rather
# than pedological -- real provincial soil differences would appear in both.
#
# ldepth_max: per-province SDs make the mechanism explicit --
#   AB 112.05 (SD 14.96) | MB 101.53 (SD 5.30) | SK 100.02 (SD 0.39)
# An SD of 0.39 cm province-wide is not a soil property. The Saskatchewan Soil
# Resource Database manual (SKSIDv4, AAFC Sask. Land Resource Unit, 2009,
# https://soilsofsask.ca/documents/sksid_usermanual.pdf) documents that SK
# profiles were typically examined to a 1 m "control section" -- the depth
# values are pinned at that protocol limit, not measured. Not in Olimb &
# Robinson's covariate set either, so dropping costs nothing relative to
# precedent.
#
# watertbl: AB 0.05 (SD 0.27) | SK 0.05 (SD 0.26) | MB 1.10 (SD 1.20).
# Near-constant outside MB, and collapsed to a single legend class in 6b, so
# it is dropped primarily on information grounds. (SK reports 0.0% NA here but
# 9.8% on the KP layers, suggesting zeros are being used in place of NA.)
#
# NB: NA sets across the DSS layers are fully nested (pct_na_excl = 0.00 for
# every layer), so these drops recover no additional study area.
manual_drop <- c("watertbl", "ldepth_max")
final_set <- setdiff(vif_set, manual_drop)

cat("\nvifstep retained :", length(vif_set), "\n")
# vifstep retained : 33 
cat("manual drops     :", paste(manual_drop, collapse = ", "), "\n")
# manual drops     : watertbl, ldepth_max
cat("final set        :", length(final_set), "\n")
# final set        : 31 
print(vif(vals_cov[, final_set]))
# Variables      VIF
# 1        prtot 6.015263
# 2      tx_mean 6.104209
# 3    gddgrow_5 5.542645
# 4           oc 2.174727
# 5           bd 4.385874
# 6         clay 1.944330
# 7           ph 3.343898
# 8     drainage 1.722802
# 9          kp0 1.635314
# 10      kp1500 1.837376
# 11      ldepth 1.532775
# 12  psssm_mean 2.184189
# 13    psssm_sd 2.270419
# 14 tsm_growing 2.741213
# 15         tri 1.809182
# 16         tpi 1.197176
# 17         srr 1.153202
# 18   northness 1.024074
# 19    eastness 1.033287
# 20       slope 1.665545
# 21  unprot_grs 2.347722
# 22        tame 1.550770
# 23        crop 2.329628
# 24         dev 1.207536
# 25        fors 2.047951
# 26  net_return 2.555316
# 27  d_elevator 2.979618
# 28      d_road 1.242466
# 29      d_rail 1.960998
# 30    d_popctr 2.530897
# 31    d_regctr 1.938696
saveRDS(final_set, file.path(out_sa, "retained-covariates.rds"))
