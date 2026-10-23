# =============================================================================
# Purpose:
# Pairwise correlation + VIF screen to flag collinear predictors before the RF.
# RF tolerates collinearity for prediction, but it distorts variable importance
# and effect plots, so redundant predictors are screened out.
# =============================================================================
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
ecoregions_sa  <- st_read(tagged_path("./processeddatasets", ecoregions_filepath, ecoregions_filename, "studyarea"))
ecoregions_r   <- rast(tagged_path("./processeddatasets", ecoregions_filepath, paste0(tools::file_path_sans_ext(ecoregions_filename), ".tif"), "studyarea"))
ecodistricts_r <- rast(tagged_path("./processeddatasets", ecodistricts_filepath, paste0(tools::file_path_sans_ext(ecodistricts_filename), ".tif"), "studyarea"))
msk <- rast(file.path(out_sa, "predictor-valid-mask_study-area.tif"))

# baseline (2010) grassland: native + tame. Screening samples the training
# population, not the whole study area - remaining grassland skews toward
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
  geom_text(aes(label = ifelse(abs(Freq) >= 0.65, round(Freq, 2), "")), size = 2.5) +
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
# kp1500 - kp33  0.933
# tn_mean - tx_mean  0.822
# crop - native -0.793
# nrcdd - prtot -0.768
# ph - tx_mean  0.732
# bd - nrcdd  0.718
# bd - ph  0.703
cat("\nVIF (all variables):\n"); print(transform(vif_all, VIF = round(VIF, 2)))
# VIF (all variables):
#   Variables   VIF
# 1             prtot  8.88
# 2           tn_mean 11.97
# 3           tx_mean 26.99
# 4         gddgrow_5 26.56
# 5             nrcdd  8.14
# 6                oc  2.29
# 7                bd  4.63
# 8              clay  2.06
# 9                ph  3.41
# 10         drainage  2.07
# 11              kp0  2.73
# 12             kp33 19.66
# 13           kp1500 17.08
# 14           ldepth  1.69
# 15       psssm_mean  2.38
# 16         psssm_sd  2.40
# 17             elev 22.43
# 18              tri  1.75
# 19              tpi  1.18
# 20              srr  1.17
# 21        northness  1.02
# 22         eastness  1.02
# 23            slope  1.73
# 24       unprot_grs  2.56
# 25             tame  8.97
# 26             crop 11.42
# 27            settl  1.91
# 28             fors  3.90
# 29           native 27.03
# 30       net_return  2.84
# 31 dist_to_elevator  3.02
# 32     dist_to_road  1.27
# 33     dist_to_rail  2.00
# 34   dist_to_popctr  2.69
# 35   dist_to_regctr  2.02
# 36      pop_density  1.52


cat("\nvifstep (th =", vif_th, ") -- variables flagged for exclusion:\n")
print(vif_step@excluded)
# [1] "native"  "elev"    "kp33"    "tn_mean" "nrcdd"  
cat("\nRetained set:\n"); print(vif_step@results$Variables)
# Retained set:
# [1] "prtot"            "tx_mean"          "gddgrow_5"        "oc"               "bd"               "clay"             "ph"               "drainage"         "kp0"              "kp1500"          
# [11] "ldepth"           "psssm_mean"       "psssm_sd"         "tri"              "tpi"              "srr"              "northness"        "eastness"         "slope"            "unprot_grs"      
# [21] "tame"             "crop"             "settl"            "fors"             "net_return"       "dist_to_elevator" "dist_to_road"     "dist_to_rail"     "dist_to_popctr"   "dist_to_regctr"  
# [31] "pop_density"

# vifstep output; manual drops are applied after the diagnostics below, which
# supply the evidence for them
vif_set <- vif_step@results$Variables
vif(vals_cov[, vif_set])
# Variables      VIF
# 1             prtot 5.582993
# 2           tx_mean 5.358291
# 3         gddgrow_5 4.601819
# 4                oc 2.136866
# 5                bd 4.394834
# 6              clay 1.896324
# 7                ph 3.146322
# 8          drainage 1.840846
# 9               kp0 1.706614
# 10           kp1500 1.912383
# 11           ldepth 1.561171
# 12       psssm_mean 2.138445
# 13         psssm_sd 2.279272
# 14              tri 1.719838
# 15              tpi 1.215276
# 16              srr 1.161556
# 17        northness 1.030112
# 18         eastness 1.020572
# 19            slope 1.664370
# 20       unprot_grs 2.298788
# 21             tame 1.543340
# 22             crop 2.252188
# 23            settl 1.563928
# 24             fors 2.079137
# 25       net_return 2.392814
# 26 dist_to_elevator 2.924034
# 27     dist_to_road 1.262195
# 28     dist_to_rail 1.980270
# 29   dist_to_popctr 2.548878
# 30   dist_to_regctr 1.995321
# 31      pop_density 1.369435
length(vif_set)
# [1] 31

final_set <- vif_set
saveRDS(final_set, file.path(out_sa, "retained-covariates.rds"))

# ---- variance explained by ecoregion --------------------------------------
# Caluculate between-ecoregion share of variance for each predictor
names(ecoregions_r) <- "ecoregion"

d <- vals[, c("x", "y", vif_set)]

d$ecoregion <- factor(terra::extract(ecoregions_r, d[, c("x", "y")])[, 2])
d <- d[!is.na(d$ecoregion), ]
cat("n =", nrow(d), "across", nlevels(droplevels(d$ecoregion)), "ecoregions\n")
# n = 19997 across 7 ecoregions
d$ecoregion <- droplevels(d$ecoregion)

r2_eco <- sapply(setdiff(names(d), c("x", "y", "ecoregion")), function(v) {
  summary(lm(d[[v]] ~ d$ecoregion))$r.squared
})
print(round(sort(r2_eco, decreasing = TRUE), 3))
# prtot          tx_mean               ph dist_to_elevator               bd       psssm_mean        gddgrow_5             fors   dist_to_popctr         psssm_sd               oc         drainage 
# 0.713            0.580            0.546            0.542            0.525            0.393            0.375            0.323            0.294            0.222            0.205            0.189 
# net_return       unprot_grs     dist_to_rail             crop              tri             tame              kp0             clay   dist_to_regctr           kp1500            slope     dist_to_road 
# 0.186            0.151            0.127            0.121            0.119            0.112            0.095            0.091            0.074            0.062            0.055            0.054 
# settl           ldepth      pop_density              srr         eastness        northness              tpi 
# 0.033            0.029            0.008            0.006            0.005            0.001            0.001 


# Caluculate between-ecodistrict share of variance for each predictor
names(ecodistricts_r) <- "ecodistrict"
d$ecodistrict <- factor(terra::extract(ecodistricts_r, d[, c("x", "y")])[, 2])
d <- d[!is.na(d$ecodistrict), ]
cat("n =", nrow(d), "across", nlevels(droplevels(d$ecodistrict)), "ecodistricts\n")
# n = 19997 across 127 ecodistricts
#d$ecodistrict <- droplevels(d$ecodistrict)

r2_dist <- sapply(setdiff(names(d), c("x", "y", "ecoregion", "ecodistrict")), function(v) {
  summary(lm(d[[v]] ~ d$ecodistrict))$r.squared
})
print(round(sort(r2_dist, decreasing = TRUE), 3))
# prtot          tx_mean        gddgrow_5 dist_to_elevator         psssm_sd       psssm_mean               ph               bd       net_return   dist_to_regctr             clay               oc 
# 0.935            0.928            0.889            0.875            0.875            0.827            0.787            0.769            0.760            0.717            0.702            0.662 
# dist_to_popctr    fors     dist_to_rail           kp1500         drainage       unprot_grs             crop              kp0              tri             tame            slope            settl 
# 0.652            0.507            0.476            0.461            0.436            0.409            0.369            0.366            0.326            0.324            0.205            0.182 
# dist_to_road    ldepth      pop_density              srr        northness         eastness              tpi 
# 0.143            0.119            0.044            0.033            0.029            0.020            0.008 

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
# net_return         psssm_sd            prtot         drainage        gddgrow_5             clay       unprot_grs             fors          tx_mean     dist_to_rail              kp0       psssm_mean 
# 0.540            0.314            0.304            0.254            0.205            0.205            0.188            0.135            0.129            0.122            0.114            0.103 
# dist_to_elevator ldepth            tame              tri           kp1500             crop     dist_to_road            settl   dist_to_popctr            slope               ph               oc 
# 0.086            0.063            0.047            0.042            0.041            0.040            0.040            0.034            0.032            0.025            0.024            0.022 
# dist_to_regctr      bd      pop_density              tpi              srr         eastness        northness 
# 0.022            0.016            0.003            0.001            0.001            0.000            0.000 

# This is legacy. Only kept to justify why I eliminated ldepth_max and watertbl
# # per-province mean and SD for the soil variables: a province with SD near 0
# # (e.g. a depth convention capping at 100 cm) is acting as a province code
# soil_vars <- intersect(c("oc", "bd", "clay", "sand", "ph", "drainage", "watertbl",
#                          "kp0", "kp33", "kp1500", "ldepth", "ldepth_max", "slope"),
#                        names(d))
# prov_stats <- do.call(rbind, lapply(soil_vars, function(v) {
#   data.frame(var  = v,
#              prov = levels(d$province),
#              mean = round(tapply(d[[v]], d$province, mean), 2),
#              sd   = round(tapply(d[[v]], d$province, sd), 2),
#              row.names = NULL)
# }))
# print(prov_stats)
# # var         prov   mean    sd
# # 1          oc      Alberta  17.86  3.49
# # 2          oc     Manitoba  19.96  2.57
# # 3          oc Saskatchewan  18.15  2.54
# # 4          bd      Alberta   1.32  0.07
# # 5          bd     Manitoba   1.32  0.05
# # 6          bd Saskatchewan   1.33  0.05
# # 7        clay      Alberta  19.88  3.91
# # 8        clay     Manitoba  25.14  5.36
# # 9        clay Saskatchewan  23.04  2.68
# # 10         ph      Alberta   7.09  0.38
# # 11         ph     Manitoba   6.87  0.38
# # 12         ph Saskatchewan   7.08  0.39
# # 13   drainage      Alberta   3.18  0.59
# # 14   drainage     Manitoba   4.05  1.11
# # 15   drainage Saskatchewan   2.74  0.74
# # 16   watertbl      Alberta   0.05  0.27
# # 17   watertbl     Manitoba   1.10  1.20
# # 18   watertbl Saskatchewan   0.05  0.26
# # 19        kp0      Alberta  52.27  3.31
# # 20        kp0     Manitoba  50.15  9.24
# # 21        kp0 Saskatchewan  49.35  2.54
# # 22     kp1500      Alberta  16.97  4.62
# # 23     kp1500     Manitoba  12.95  4.73
# # 24     kp1500 Saskatchewan  16.89  5.49
# # 25     ldepth      Alberta  19.97 20.69
# # 26     ldepth     Manitoba  16.66 13.23
# # 27     ldepth Saskatchewan  12.17  2.25
# # 28 ldepth_max      Alberta 112.05 14.96
# # 29 ldepth_max     Manitoba 101.53  5.30
# # 30 ldepth_max Saskatchewan 100.02  0.39
# # 31      slope      Alberta   5.46  5.72
# # 32      slope     Manitoba   4.35  9.33
# # 33      slope Saskatchewan   6.43  4.97
# write.csv(prov_stats, file.path(out_dir, "province_stats_soil.csv"), row.names = FALSE)

# ---- where are the DSS NAs? ----
# Haven't rerun this since updating predictors, come back to it when trying to do imputation
# # The DSS layers are the binding NA constraint (~35.5M cells on kp0/kp1500).
# # Dropping a layer only recovers study area if its NAs are attribute gaps in
# # surveyed polygons; NAs where DSS has no polygon coverage at all are lost
# # either way. Separate draw from d: that sample is all-valid by construction,
# # so it contains no NA rows to count. Sampled over the whole study area
# # (not the valid mask), NAs retained.
# na_n <- 500000   # ~0.1% precision on a proportion; cheap relative to a full pass
# 
# set.seed(43)
# na_pts <- spatSample(sa_r, na_n * 3, method = "random",
#                      xy = TRUE, values = TRUE, na.rm = FALSE)
# na_pts <- as.matrix(na_pts[!is.na(na_pts[[3]]), c("x", "y")])
# na_pts <- na_pts[sample(nrow(na_pts), min(na_n, nrow(na_pts))), ]
# 
# na_d <- terra::extract(r, na_pts)
# na_d$ID <- NULL
# na_d$province <- factor(terra::extract(prov, na_pts)$NAME_1)
# na_d <- na_d[!is.na(na_d$province), ]
# 
# dss_vars <- intersect(c("drainage", "ph", "kp0", "kp33", "kp1500",
#                         "ldepth", "slope"), names(na_d))
# 
# # % NA per DSS layer within each province
# na_by_prov <- sapply(dss_vars, function(v) {
#   round(100 * tapply(is.na(na_d[[v]]), na_d$province, mean), 1)
# })
# cat("\n% NA by province:\n"); print(na_by_prov)
# # % NA by province:
# #   drainage watertbl kp0 kp33 kp1500 ldepth ldepth_max slope
# # Alberta           2.2      1.1 2.2  2.3    2.2    0.1        0.1   0.1
# # Manitoba          6.2      6.2 6.8  6.8    6.8    6.8        6.8   6.2
# # Saskatchewan      1.7      0.0 9.8  9.8    9.8    9.8        9.8   1.7
# 
# # is a layer's NA set exclusive, or shared with the other DSS layers?
# # exclusive NAs are what dropping that layer would actually recover
# dss_na <- sapply(dss_vars, function(v) is.na(na_d[[v]]))
# excl <- sapply(dss_vars, function(v) {
#   others <- rowSums(dss_na[, setdiff(dss_vars, v), drop = FALSE]) > 0
#   sum(dss_na[, v] & !others)
# })
# cat("\nDSS NA overlap (n =", nrow(na_d), "sampled):\n")
# # DSS NA overlap (n = 499428 sampled):
# print(data.frame(var        = dss_vars,
#                  pct_na     = round(100 * colMeans(dss_na), 1),
#                  pct_na_excl = round(100 * excl / nrow(na_d), 2),
#                  row.names  = NULL))
# # var pct_na pct_na_excl
# # 1   drainage    2.5        0.00
# # 2   watertbl    1.3        0.00
# # 3        kp0    6.8        0.00
# # 4       kp33    6.8        0.05
# # 5     kp1500    6.8        0.00
# # 6     ldepth    6.1        0.00
# # 7 ldepth_max    6.1        0.00
# # 8      slope    1.8        0.00
# cat("\nany DSS layer NA:", round(100 * mean(rowSums(dss_na) > 0), 1), "%\n")
# # any DSS layer NA: 6.8 %
# write.csv(as.data.frame(na_by_prov), file.path(out_dir, "province_na_dss.csv"))
