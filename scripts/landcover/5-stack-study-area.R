# =============================================================================
# Purpose:
# Assemble the predictor stack, mask each layer to the study area, and write
# to disk with a _study-area suffix. 
#
# Also computes the union of NA cells across the stack (excluding water) 
# -ranger drops any row with an NA in any predictor, so the union, not any 
# single layer, sets the predictor eligibility envelope
# =============================================================================

source(file.path("scripts/landcover", "0a-setup.R"))
source(file.path("scripts/landcover", "0b-predictor-files.R"))
terraOptions(memfrac = 0.5)

gopt_for <- function(datatype) {
  base <- c("COMPRESS=DEFLATE", "ZLEVEL=9", "TILED=YES",
            "BLOCKXSIZE=256", "BLOCKYSIZE=256",
            "BIGTIFF=YES", "NUM_THREADS=ALL_CPUS")
  c(base, if (grepl("^FLT", datatype)) "PREDICTOR=3" else "PREDICTOR=2")
}

out_sa <- file.path("./processeddatasets", "predictors", paste0(res_m, "m"))
sa_r <- rast(file.path(sa_filepath, sa_raster_filename))
water <- rast(file.path("./processeddatasets", aafc_filepath, "binary", "AAFC2010_water.tif"))

# ---- source layers ----
files <- predictor_files

missing <- files[!file.exists(files)]
if (length(missing)) {
  stop("missing predictor layers:\n  ", paste(names(missing), missing, sep = ": ", collapse = "\n  "))
}

sa_file <- function(nm) file.path(out_sa, paste0(nm, "_study-area.tif"))

# ---- mask to study area ----
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

# ---- assemble and check ----
stk <- rast(sapply(names(files), sa_file))
names(stk) <- names(files)

# geometry standardised upstream in 3-crop; but still double check
stopifnot(compareGeom(stk, sa_r, stopOnError = FALSE))

# count of non-NA cells in study area
n_sa <- global(!is.na(sa_r), "sum", na.rm = TRUE)[[1]]

# per-layer NA counts: any layer well above the others is the binding one
na_by_layer <- global(is.na(stk), "sum", na.rm = TRUE)
na_by_layer$valid <- ncell(stk) - na_by_layer$sum
na_by_layer$lost_vs_sa <- na_by_layer$valid - n_sa      # negative = cells lost

# union: cells valid in every layer, excluding water. This is the predictor
# eligibility envelope - not the training population, which is additionally
# restricted to baseline grassland (see 6a).
wtr <- mask(crop(water, sa_r), sa_r)
wtr <- subst(wtr, NA, 0) == 1          # NA outside water == not water
ok <- !any(is.na(stk)) & !wtr
writeRaster(ok, file.path(out_sa, "predictor-valid-mask_study-area.tif"),
            datatype = "INT1U", overwrite = TRUE, gdal = gopt_for("INT1U"))

n_ok <- global(ok, "sum", na.rm = TRUE)[[1]]

cat("\nstudy area cells          :", format(n_sa, big.mark = ","), "\n")
# study area cells          : 519,027,074 
cat("valid in every predictor  :", format(n_ok, big.mark = ","), "\n")
# valid in every predictor  : 452,372,582 
cat("lost to NA union + water  :", format(n_sa - n_ok, big.mark = ","),
    sprintf("(%.2f%%)\n", 100 * (n_sa - n_ok) / n_sa))
# lost to NA union + water  : 66,654,492 (12.84%)

print(na_by_layer[order(na_by_layer$lost_vs_sa), ])
# layer       sum     valid lost_vs_sa
# 12             kp33 472877733 483254289  -35772785
# 11              kp0 472606377 483525645  -35501429
# 13           kp1500 472606377 483525645  -35501429
# 14           ldepth 468879128 487252894  -31774180
# 30       net_return 468059031 488072991  -30954083
# 36      pop_density 456764855 499367167  -19659907
# 10         drainage 450674814 505457208  -13569866
# 6                oc 448803917 507328105  -11698969
# 7                bd 448803917 507328105  -11698969
# 8              clay 448803917 507328105  -11698969
# 9                ph 448803917 507328105  -11698969
# 23            slope 446955600 509176422   -9850652
# 15       psssm_mean 438725748 517406274   -1620800
# 16         psssm_sd 438725748 517406274   -1620800
# 1             prtot 437469681 518662341    -364733
# 2           tn_mean 437469681 518662341    -364733
# 3           tx_mean 437469681 518662341    -364733
# 4         gddgrow_5 437469681 518662341    -364733
# 5             nrcdd 437469681 518662341    -364733
# 17             elev 437104948 519027074          0
# 18              tri 437104948 519027074          0
# 19              tpi 437104948 519027074          0
# 20              srr 437104948 519027074          0
# 21        northness 437104948 519027074          0
# 22         eastness 437104948 519027074          0
# 24       unprot_grs 437104948 519027074          0
# 25             tame 437104948 519027074          0
# 26             crop 437104948 519027074          0
# 27            settl 437104948 519027074          0
# 28             fors 437104948 519027074          0
# 29           native 437104948 519027074          0
# 31 dist_to_elevator 437104948 519027074          0
# 32     dist_to_road 437104948 519027074          0
# 33     dist_to_rail 437104948 519027074          0
# 34   dist_to_popctr 437104948 519027074          0
# 35   dist_to_regctr 437104948 519027074          0


write.csv(data.frame(layer = rownames(na_by_layer), na_by_layer),
          file.path(out_sa, "predictor-na-summary.csv"), row.names = FALSE)
