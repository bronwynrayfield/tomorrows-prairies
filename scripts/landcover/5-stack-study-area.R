# 5-stack-study-area.R
# tomorrows-prairies
#
# Assembles the biophysical predictor stack from the per-block outputs of
# 4b-4e, masks each layer to the study area, and writes it to a single
# directory with a _study-area suffix. Also computes the union of NA cells
# across the stack (excluding water) -- ranger drops any row with an NA in
# any predictor, so the union, not any single layer, sets the predictor
# eligibility envelope. The training population is this envelope intersected
# with baseline (2010) grassland, since only grassland cells can convert.

source(file.path("scripts/landcover", "0a-setup.R"))
source(file.path("scripts/landcover", "0b-predictor-files.R"))
terraOptions(memfrac = 0.5)

gopt_for <- function(datatype) {
  base <- c("COMPRESS=DEFLATE", "ZLEVEL=9", "TILED=YES",
            "BLOCKXSIZE=256", "BLOCKYSIZE=256",
            "BIGTIFF=YES", "NUM_THREADS=ALL_CPUS")
  c(base, if (grepl("^FLT", datatype)) "PREDICTOR=3" else "PREDICTOR=2")
}

out_sa <- file.path("./processeddatasets", "predictors", "30m")
sa_r <- rast(file.path(sa_filepath, sa_raster_filename))
water <- rast(file.path("./processeddatasets", aafc_filepath, "binary", "AAFC2010_water.tif"))

# ---- source layers -------------------------------------------------------
files <- predictor_files

missing <- files[!file.exists(files)]
if (length(missing)) {
  stop("missing predictor layers:\n  ", paste(names(missing), missing, sep = ": ", collapse = "\n  "))
}

sa_file <- function(nm) file.path(out_sa, paste0(nm, "_study-area.tif"))

# ---- mask to study area --------------------------------------------------
# skip only if the cached copy exists AND is newer than its source, so
# regenerated upstream layers invalidate their _study-area copies
for (nm in names(files)) {
  if (file.exists(sa_file(nm)) && file.mtime(sa_file(nm)) > file.mtime(files[[nm]])) {
    cat("skip", nm, "\n"); next
  }
  cat("masking", nm, format(Sys.time()), "\n")
  r <- rast(files[[nm]])
  names(r) <- nm
  mask(crop(r, sa_r), sa_r,
       filename = sa_file(nm), overwrite = TRUE,
       wopt = list(gdal = gopt_for(datatype(r))))
}

# ---- assemble and check --------------------------------------------------
stk <- rast(sapply(names(files), sa_file))
names(stk) <- names(files)

# geometry standardised upstream in 3-crop; but still double check
stopifnot(compareGeom(stk, sa_r, stopOnError = FALSE))

n_sa <- global(!is.na(sa_r), "sum", na.rm = TRUE)[[1]]

# per-layer NA counts: any layer well above the others is the binding one
na_by_layer <- global(is.na(stk), "sum", na.rm = TRUE)
na_by_layer$valid    <- ncell(stk) - na_by_layer$sum
na_by_layer$lost_vs_sa <- na_by_layer$valid - n_sa      # negative = cells lost

# union: cells valid in every layer, excluding water. This is the predictor
# eligibility envelope - not the training population, which is additionally
# restricted to baseline grassland (see 6a).
wtr <- crop(water, sa_r)
wtr <- subst(wtr, NA, 0) == 1          # NA outside water == not water
ok <- !any(is.na(stk)) & !wtr
writeRaster(ok, file.path(out_sa, "predictor-valid-mask_study-area.tif"),
            datatype = "INT1U", overwrite = TRUE, gdal = gopt_for("INT1U"))

n_ok <- global(ok, "sum", na.rm = TRUE)[[1]]

cat("\nstudy area cells          :", format(n_sa, big.mark = ","), "\n")
# study area cells          : 519,027,074 
cat("valid in every predictor  :", format(n_ok, big.mark = ","), "\n")
# valid in every predictor  : 474,953,890 
cat("lost to NA union + water  :", format(n_sa - n_ok, big.mark = ","),
    sprintf("(%.2f%%)\n", 100 * (n_sa - n_ok) / n_sa))
# lost to NA union + water  : 44,073,184 (8.49%)

print(na_by_layer[order(na_by_layer$lost_vs_sa), ])
# sum     valid lost_vs_sa
# kp33         472891728 483240294  -35786780
# kp0          472620372 483511650  -35515424
# kp1500       472620372 483511650  -35515424
# ldepth       468893249 487238773  -31788301
# ldepth_max   468893249 487238773  -31788301
# drainage     450685115 505446907  -13580167
# oc           448803917 507328105  -11698969
# bd           448803917 507328105  -11698969
# clay         448803917 507328105  -11698969
# sand         448803917 507328105  -11698969
# silt         448803917 507328105  -11698969
# ph           448803917 507328105  -11698969
# slope        446965867 509166155   -9860919
# watertbl     444251628 511880394   -7146680
# tsm_planting 439939324 516192698   -2834376
# tsm_growing  439939324 516192698   -2834376
# psssm_mean   438725748 517406274   -1620800
# psssm_sd     438725748 517406274   -1620800
# prtot        437469681 518662341    -364733
# tn_mean      437469681 518662341    -364733
# tx_mean      437469681 518662341    -364733
# gddgrow_5    437469681 518662341    -364733
# nrcdd        437469681 518662341    -364733
# elev         437104948 519027074          0
# tri          437104948 519027074          0
# tpi          437104948 519027074          0
# srr          437104948 519027074          0
# northness    437104948 519027074          0
# eastness     437104948 519027074          0

write.csv(data.frame(layer = rownames(na_by_layer), na_by_layer),
          file.path(out_sa, "predictor-na-summary.csv"), row.names = FALSE)
