# =============================================================================
# Purpose: 
# Neighbourhood class proportions from lulc_2010 (3b) and the forest binary
# Baseline predictors: lulc_2020 is never used here.
#
# AWS instance - r6i.2xlarge (8 vCPU, 64 GB)
#
# Window: 3210 m (107 cells, halo 53), nearest odd-cell match to 2 miles,
#   the section-scale counterpart of the 810 m (quarter-section) terrain window.
#
# Classes: native, tame, crop, dev (settlement), fors (forest = lulc class 5 &
#   forest binary), other (class 5 minus forest). Wetland and barren not split
#   out (unreliable lineage; negligible area). The six sum to 1, so `other` is
#   computed for the closure check but withheld from 0b-predictor-files.R.
#   unprot (native + tame outside CPCAD protected areas) is a subset, not a
#   partition class - excluded from the closure check.
#
# Algorithm: separable box sums (107x1 then 1x107, ~53x fewer ops). NA -> 0 in
#   both numerator and valid-cell count, prop = num / den, which reproduces
#   focal(na.rm = TRUE) exactly (max|diff| = 0; test-separable-focal.R). den is
#   shared across classes. Tiled 4x4 because full-extent focal crashes on Windows;
#   halo must equal the focal radius.
#
# Edge: focal runs on the modelling extent, then is cropped and masked to sa_r.
#   me_r's buffer must be >= 53 cells around sa_r except along the US border,
#   where lulc is NA. There, windows within 1.6 km (53 rows) are truncated to
#   their Canadian portion (>= ~5,800 of 11,449 cells)
# =============================================================================

source(file.path("scripts/landcover", "0a-setup.R"))

# ---- Parameters ----
win_size_m <- 3210
# Other sizes considered but not run
# win_sizes_m <- c(
#   3200,
#   6000,
#   10000
# )
tile_nx    <- 4
tile_ny    <- 4

win_cells <- round(win_size_m / res_m)
if (win_cells %% 2 == 0) win_cells <- win_cells + 1
halo      <- (win_cells - 1) / 2
stopifnot(halo == 53)   # check - retile if the window changes

wr <- matrix(1, nrow = 1, ncol = win_cells)
wc <- matrix(1, nrow = win_cells, ncol = 1)

# ---- Paths ----
binary_dir <- file.path("./processeddatasets", aafc_filepath, "binary")
pred_dir   <- file.path("processeddatasets", "predictors", sprintf("%dm", res_m))
tile_dir   <- "D:/tile_scratch/landscape_context"
dir.create(file.path(".", pred_dir), recursive = TRUE, showWarnings = FALSE)
dir.create(tile_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Load ----
lulc_2010 <- rast(tagged_path("./processeddatasets", "landcover", "lulc_2010.tif", "modelling-extent"))
forest    <- rast(file.path(binary_dir, "AAFC2010_forest.tif"))
sa_r      <- rast(file.path(sa_filepath, sa_raster_filename))
pa_2010   <- rast(file.path("processeddatasets", "management", "Protected Areas Combined", "CPCAD_pa_2010.tif")) 
stopifnot(compareGeom(lulc_2010, forest, pa_2010, stopOnError = FALSE))

# Forest must sit entirely inside lulc class 5, or the six-way split is not a
# partition and the closure argument fails.
# leak <- global(forest == 1 & lulc_2010 != cls_other, "sum", na.rm = TRUE)[[1]]
# message(sprintf("forest outside lulc class 5: %s px", format(leak, big.mark = ",")))
# stopifnot(leak == 0)

# Class definitions as functions of the (extended) lulc and forest tiles.
# class_num <- list(
#   native = function(l, f) as.numeric(l == cls_native),
#   tame   = function(l, f) as.numeric(l == cls_tame),
#   crop   = function(l, f) as.numeric(l == cls_annual),
#   dev    = function(l, f) as.numeric(l == cls_settlement),
#   fors   = function(l, f) as.numeric(l == cls_other & f == 1),
#   other  = function(l, f) as.numeric(l == cls_other & f != 1)
# )

class_num <- list(
  # pa_2010 is 0/1 inside me_r, NA outside - test == 0 for unprotected.
  # NA only where lulc is also NA, so the ifel(is.na(tl), 0, ...) guard covers it.
  unprot = function(l, f, p) as.numeric(l %in% c(cls_native, cls_tame) & p == 0)
)

# Just run a subset
# class_num <- list(
#   tame   = function(l, f, p) as.numeric(l == cls_tame),
#   crop   = function(l, f, p) as.numeric(l == cls_annual),
#   unprot = function(l, f, p) as.numeric(l %in% c(cls_native, cls_tame) & p == 0)
# )

# ---- Tiled separable focal, all classes per tile ----
# Tiles cut on cell indices so they abut exactly. Each tile: crop to the
# clamped buffer, extend to the full (unclamped) buffer with NA fill, build
# den once and num per class, two box-sum passes each, crop to core, write.

box_sum <- function(x, f) {
  h <- focal(x, w = wr, fun = "sum", na.rm = FALSE,
             filename = sub("\\.tif$", "_h.tif", f), overwrite = TRUE,
             wopt = list(datatype = "INT2U"))
  focal(h, w = wc, fun = "sum", na.rm = FALSE, filename = f, overwrite = TRUE,
        wopt = list(datatype = "INT2U"))
}

e  <- ext(lulc_2010)
rs <- res(lulc_2010)[1]
cb <- round(seq(0, ncol(lulc_2010), length.out = tile_nx + 1))
rb <- round(seq(0, nrow(lulc_2010), length.out = tile_ny + 1))

core_files <- setNames(vector("list", length(class_num)), names(class_num))

for (i in seq_len(tile_nx)) for (j in seq_len(tile_ny)) {
  tid  <- sprintf("t%02d%02d", i, j)
  core <- ext(xmin(e) + cb[i]     * rs, xmin(e) + cb[i + 1] * rs,
              ymax(e) - rb[j + 1] * rs, ymax(e) - rb[j]     * rs)
  buf_full <- ext(xmin(core) - halo * rs, xmax(core) + halo * rs,
                  ymin(core) - halo * rs, ymax(core) + halo * rs)
  buf_clmp <- ext(max(xmin(buf_full), xmin(e)), min(xmax(buf_full), xmax(e)),
                  max(ymin(buf_full), ymin(e)), min(ymax(buf_full), ymax(e)))
  
  # NA fill beyond the physical edge; the num/den substitution turns it into 0.
  tl <- extend(crop(lulc_2010, buf_clmp), buf_full)
  tf <- extend(crop(forest,    buf_clmp), buf_full)
  tp <- extend(crop(pa_2010, buf_clmp), buf_full)
  
  den <- ifel(is.na(tl), 0, 1,
              filename = file.path(tile_dir, sprintf("den_%s.tif", tid)),
              overwrite = TRUE, wopt = list(datatype = "INT1U"))
  db  <- box_sum(den, file.path(tile_dir, sprintf("denbox_%s.tif", tid)))
  
  for (nm in names(class_num)) {
    num <- ifel(is.na(tl), 0, class_num[[nm]](tl, tf, tp),
                filename = file.path(tile_dir, sprintf("num_%s_%s.tif", nm, tid)),
                overwrite = TRUE, wopt = list(datatype = "INT1U"))
    nb  <- box_sum(num, file.path(tile_dir, sprintf("numbox_%s_%s.tif", nm, tid)))
    pc  <- crop(nb / db, core,
                filename = file.path(tile_dir, sprintf("core_%s_%s.tif", nm, tid)),
                overwrite = TRUE, wopt = list(datatype = "FLT4S"))
    core_files[[nm]] <- c(core_files[[nm]],
                          file.path(tile_dir, sprintf("core_%s_%s.tif", nm, tid)))
  }
  
  # Drop this tile's intermediates; keep only core_* files. Caps transient
  # disk at ~one tile plus the accumulating cores.
  unlink(file.path(tile_dir,
                   sprintf(c("den_%s.tif", "denbox_%s.tif", "denbox_%s_h.tif"), tid)))
  unlink(list.files(tile_dir, sprintf("^num.*_%s(_h)?\\.tif$", tid),
                    full.names = TRUE))
  
  message(sprintf("tile %s done (%d/%d)", tid, (i - 1) * tile_ny + j, tile_nx * tile_ny))
}

# ---- Virtual mosaic, crop to study area, write, dispose ----
# vrt() mosaics without writing (the cores abut exactly on the same grid, which
# the cell-index tiling guarantees). Each class's cores are deleted as soon as
# its output exists, so cores shrink as outputs grow.
props <- list()
for (nm in names(class_num)) {
  v   <- vrt(core_files[[nm]])
  out <- tagged_path("./processeddatasets", file.path(aafc_filepath, paste0("proportion_", sprintf("%dm", win_cells * res_m))),
                     sprintf("%s_prop_%dm.tif", nm, win_cells * res_m), "2010_modelling-extent")
  props[[nm]] <- mask(crop(v, sa_r), sa_r,
                      filename = out, overwrite = TRUE,
                      wopt = list(datatype = "FLT4S",
                                  gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3")))
  names(props[[nm]]) <- nm
  unlink(core_files[[nm]])
  message(sprintf("wrote %s", out))
}

# ---- QA/QC ----
p <- rast(props)
pth <- function(nm) tagged_path("./processeddatasets", file.path("predictors", sprintf("%dm", res_m)),
                                sprintf("%s_prop_%dm.tif", nm, win_cells * res_m), sprintf("%dm", res_m))
p <- rast(sapply(c("native", "tame", "crop", "dev", "fors", "other"), pth))


print(minmax(p))                                   # all layers in [0, 1]

s <- sum(p)                                        # closure: six sum to 1
print(minmax(s))
stopifnot(all(abs(minmax(s) - 1) < 1e-4))

print(global(is.na(p), "sum", na.rm = TRUE))       # identical NA footprint

# Withheld reference: if `other` is near zero over most of the ecozone the
# retained five are close to closed -- see header.
print(global(p$other, c("mean", "sd"), na.rm = TRUE))
hist(p$other, main = "other (withheld closure reference)", xlab = "proportion")

plot(p, nc = 3)

unlink(tile_dir, recursive = TRUE)
