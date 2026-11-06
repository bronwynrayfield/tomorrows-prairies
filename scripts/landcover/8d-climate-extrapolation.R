# =============================================================================
# Purpose:
#   Checks whether simulation-period (future) climate predictors fall outside
#   the range seen in RF fitting, for each conversion transition fit in 8a.
#
# Inputs:
#   - Fitting samples from 7c-fitting-sample.R (train-sample_<from>.rds)
#   - Retained covariate list (retained-covariates.rds)
#   - Future climate rasters (modelling-extent) from 3-crop-mask
#
# Outputs (tag = conversion_<from>-<to>, as in 8a):
#   - climate_extrapolation_summary.csv: training range, share of cells out of
#     range and max exceedance, per variable and transition
#   - climate_n_oor_<tag>.tif: number of climate variables out of range
#   - climate_mess_<tag>.tif: MESS (< 0 = novel, incl. novel combinations)
#
# Notes:
#   - Non-imputed fut layers only: imputed cells are filled from neighbours,
#     so they can't tell us anything about extrapolation
#   - Range comes from the rows each 8a model was fit on (persistence plus one
#     destination), not the hist raster or the full origin sample: both are
#     wider than what the RF saw and would undercount extrapolation
#   - Evaluated only on cells in the transition's origin class, since those
#     are the only cells that model scores
#   - ranger holds predictions flat beyond the training range, so exceedance
#     size matters more than whether it happens; exc is expressed as a fraction
#     of the training range (0 = inside)
# =============================================================================

source(file.path("scripts/landcover", "0a-setup.R"))
source(file.path("scripts/landcover", "0b-predictor-files.R"))
library(predicts)

# same transitions as 8a
transitions <- data.frame(from = c("native", "tame",   "native",     "tame",       "native"),
                          to   = c("annual", "annual", "settlement", "settlement", "tame"))

out_dir <- "output/landcover/model"
in_sa   <- file.path("./processeddatasets", "predictors", "30m")

#--- Load ----
# climate vars the models actually use (vifstep dropped tn_mean; nr_cdd dropped mechanistically)
final_set <- readRDS(file.path(in_sa, "retained-covariates.rds"))
vars <- intersect(c("gddgrow_5", "hddheat_18", "prtot", "tn_mean", "tx_mean", "nr_cdd"), final_set)

sa_r <- rast(file.path(sa_filepath, sa_raster_filename))
fut <- rast(lapply(vars, \(v) rast(tagged_path("./processeddatasets", climate_filepath,
                                               get(paste0(v, "_fut_filename")), "modelling-extent"))))
names(fut) <- vars
fut <- mask(crop(fut, sa_r), sa_r)
fact <- 30   # 30 m -> 900 m: 900x fewer cells, still ~7-10 blocks across each ~6-9 km climate cell
fut  <- aggregate(fut, fact, "mean", na.rm = TRUE)

rngs <- lapply(seq_len(nrow(transitions)), \(j) {
  smp <- readRDS(file.path(out_dir, paste0("train-sample_", transitions$from[j], ".rds")))
  t(sapply(smp[smp$lulc_2020 %in% c(transitions$from[j], transitions$to[j]), vars], range, na.rm = TRUE))
})
names(rngs) <- sprintf("%s-%s", transitions$from, transitions$to)
rngs

# LULC at simulation start
origin_r <- rast(tagged_path("./processeddatasets", grass_consensus_filepath, grass_consensus_filename, "study-area"))

# Narrowest training range per origin (largest min, smallest max across its transitions),
# so p_oor is an upper bound for every transition from that origin. Mins are identical
# within origin; maxes differ: native gddgrow_5 2022-2062 (39 GDD, ~3.5% of range),
# prtot 783-801 mm; tame prtot 788-796 mm.
res <- list()
for (o in unique(transitions$from)) {
  r   <- rngs[transitions$from == o]
  rng <- cbind(do.call(pmax, lapply(r, \(x) x[, 1])), do.call(pmin, lapply(r, \(x) x[, 2])))
  
  # tame includes native: native->tame moves cells into tame during the simulation
  cond <- if (o == "tame") origin_r %in% c(cls_native, cls_tame) else origin_r == get(paste0("cls_", o))
  w    <- aggregate(cond, fact, "sum", na.rm = TRUE)   # origin pixels per block
  fo   <- mask(fut, w, maskvalues = c(0, NA))
  
  below <- fo < rng[, 1]
  above <- fo > rng[, 2]
  oor   <- below | above
  exc   <- ((fo - rng[, 2]) * above + (rng[, 1] - fo) * below) / (rng[, 2] - rng[, 1])
  
  res[[o]] <- data.frame(
    origin    = o, var = vars,
    train_min = rng[, 1], train_max = rng[, 2],
    fut_min   = unlist(global(fo, "min", na.rm = TRUE)),
    fut_max   = unlist(global(fo, "max", na.rm = TRUE)),
    p_oor     = unlist(global(oor, "mean", weights = w, na.rm = TRUE)),   # share of origin pixels
    max_exc   = unlist(global(exc, "max", na.rm = TRUE))
  )
  
  writeRaster(sum(oor), file.path(out_dir, sprintf("climate_n_oor_%s.tif", o)),
             overwrite = TRUE, datatype = "INT1U", gdal = c("COMPRESS=DEFLATE"))

  # MESS reference = pooled training rows for this origin, so it's slightly less
  # conservative than p_oor; it's here for novel combinations, not range
  smp <- readRDS(file.path(out_dir, paste0("train-sample_", o, ".rds")))
  smp <- smp[smp$lulc_2020 %in% c(o, transitions$to[transitions$from == o]), vars]
  writeRaster(mess(fo, smp), file.path(out_dir, sprintf("climate_mess_%s.tif", o)),
             overwrite = TRUE, gdal = c("COMPRESS=DEFLATE"))
}

res <- do.call(rbind, res)
print(res)
#                  origin       var  train_min  train_max     fut_min   fut_max        p_oor   max_exc
# native.gddgrow_5 native gddgrow_5 929.690796 2022.38403 1088.797916 2329.6469 0.5485519952 0.2811977
# native.prtot     native     prtot 292.224487  782.96106  300.873352  843.6585 0.0001319397 0.1236863
# native.tx_mean   native   tx_mean   7.245274   13.70260    8.400085   14.8733 0.0398594760 0.1812986
# tame.gddgrow_5     tame gddgrow_5 977.213379 2063.91382 1088.797916 2329.6469 0.3964657518 0.2445320
# tame.prtot         tame     prtot 292.656006  787.98047  300.873352  843.6585 0.0001094988 0.1124071
# tame.tx_mean       tame   tx_mean   7.092031   13.67002    8.400085   14.8733 0.0290843213 0.1829253
write.csv(res, file.path(out_dir, "climate_extrapolation_summary.csv"), row.names = FALSE)

g     <- mask(fut[["gddgrow_5"]], w, maskvalues = c(0, NA))
gmax  <- min(sapply(rngs[transitions$from == "native"], \(x) x["gddgrow_5", 2]))   # narrowest native max, 2022
g_exc <- ifel(g > gmax, g - gmax, 0)   # GDD past the training max; 0 = within range
plot(g_exc)
writeRaster(g_exc, file.path(out_dir, "climate_gdd_exc_native.tif"), overwrite = TRUE)

# hist check: 0.009% of native pixels exceed the training GDD max (2022) under the
# 1991-2020 normal vs 55% under 2021-2050, so exceedance is climate shift, not sampling
gh <- rast(tagged_path("./processeddatasets", climate_filepath, gddgrow_5_hist_filename, "modelling-extent"))
gh <- aggregate(mask(crop(gh, sa_r), sa_r), fact, "mean", na.rm = TRUE)
global(gh > 2022.384, "mean", weights = mask(w, w, maskvalues = 0), na.rm = TRUE)