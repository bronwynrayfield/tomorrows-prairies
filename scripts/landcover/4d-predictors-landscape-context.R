# =============================================================================
# 4d-predictors-landscape-context.R
# =============================================================================
#
# Focal proportion predictors: neighbourhood composition around each pixel as
# fractions of native grassland, tame grassland, annual cropland, settlement,
# forest and other. Reads lulc_2010 from 3b and never re-derives it; lulc_2020
# is NEVER used here -- these are baseline predictors.
#
# Window: 3210 m (107 cells), closest odd-cell approximation to 2 miles --
# section-grid aligned, the DLS logic of the 810 m terrain window one scale up.
#
# CLOSURE -- six computed, five used:
#   The six proportions sum to 1 at every valid pixel, so VIF is infinite as a
#   set. `other` is withheld as the closure reference and NOT passed to
#   0b-predictor-files.R (sand/silt treatment). It is still computed: it is the
#   exhaustiveness check, and a later reference-class change would otherwise
#   cost a full rerun. The retained five still sum to (1 - other); where
#   `other` is small they are near-closed and vifstep may flag them -- the QA
#   below prints its distribution so that is visible before screening.
#
# Class sourcing:
#   native, tame, annual, settlement from lulc_2010 (codes in 0a). forest is
#   the 3a LUTS binary and sits inside lulc class 5, so `other` = class 5 minus
#   forest. Wetland is deliberately NOT split out (mapping lineage unreliable;
#   the drainage signal is carried by dss_drainage). Barren is not split out
#   (negligible area; the watertbl failure mode).
#
# TODO -- unprotected grassland proportion: blocked on the CPCAD protection
#   layer; see placeholder at end. Not a closure member (subset of native+tame,
#   not a partition class).
#
# ALGORITHM -- separable box filter:
#   A uniform focal mean is separable: 107x1 then 1x107 box sums, ~53x fewer
#   operations than the 107x107 window. focal(na.rm = TRUE) does not compose
#   across two passes, so NA is removed first and restored after:
#     num  = class indicator with NA -> 0     box sum = class cells in window
#     den  = validity indicator, NA -> 0      box sum = valid cells in window
#     prop = num / den, masked back to the source
#   This reproduces na.rm = TRUE semantics EXACTLY (verified against direct
#   focal: max|diff| = 0, zero NA mismatches, interior and ecozone-boundary
#   crops -- test-separable-focal.R). Tiles are extended with NA fill before
#   the substitution, so the physical raster edge behaves like the internal NA
#   boundary: beyond-edge contributes 0 to num and den alike.
#
#   den depends only on the NA footprint, so it is computed once per tile and
#   shared across all six classes: 7 box-sum pairs per tile, not 12.
#
#   Still tiled (4x4, halo = 53): the Windows focal crash is about raster
#   size, not window size. The halo MUST equal the focal radius.
#
# DISK BUDGET (fits in 50 GB free on D:):
#   tile_dir MUST be on D: (never tempdir() -- that is C: on Windows).
#   Datatypes: num/den INT1U, box sums INT2U (max 107^2 = 11449), final ratio
#   FLT4S. Per-tile intermediates (~2.5 GB) are deleted at the end of each
#   tile. Core files accumulate to ~26.5 GB by the end of the tile loop --
#   this is the peak. The output loop uses vrt() (virtual mosaic, nothing
#   written) and deletes each class's cores as soon as its output exists, so
#   cores shrink as outputs grow. Outputs are DEFLATE-compressed. Peak ~30 GB;
#   50 GB free is enough with headroom. Before launching, confirm the dead
#   direct-run tiles and D:/terra_tmp_test are actually deleted.
#
# Edge handling: focal domain is the modelling extent; cropped and masked to
#   sa_r at the end, per convention. Confirm once that me_r's buffer is >= 53
#   cells on all four sides of sa_r, or study-area edge pixels silently get
#   truncated windows.
#
# Dependencies:  0a-setup.R (incl. cls_* codes, sa_filepath, sa_raster_filename),
#                3b-lulc-maps.R, 3a-crop-data.R (forest binary)
# Consumed by:   0b-predictor-files.R (five layers; `other` withheld)
# Verified by:   test-separable-focal.R
# =============================================================================

source(file.path("scripts/landcover", "0a-setup.R"))

# ---- Parameters ----
win_size_m <- 3200
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
stopifnot(halo == 53)   # tripwire: retile consciously if the window changes

wr <- matrix(1, nrow = 1, ncol = win_cells)
wc <- matrix(1, nrow = win_cells, ncol = 1)

# ---- Paths ----
binary_dir <- file.path("./processeddatasets", aafc_filepath, "binary")
pred_dir   <- file.path("processeddatasets", "predictors", sprintf("%dm", res_m))
tile_dir   <- "D:/tile_scratch/landscape_context"
dir.create(file.path(".", pred_dir), recursive = TRUE, showWarnings = FALSE)
dir.create(tile_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Load ----
lulc_2010 <- rast(tagged_path("./processeddatasets", "landcover", "lulc_2010.tif", "30m"))
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
  # pa_2010 is 0/1 inside me_r, NA outside -- test == 0 for unprotected.
  # NA only where lulc is also NA, so the ifel(is.na(tl), 0, ...) guard covers it.
  unprot = function(l, f, p) as.numeric(l %in% c(cls_native, cls_tame) & p == 0)
)

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
  out <- tagged_path("./processeddatasets", file.path("predictors", sprintf("%dm", res_m)),
                     sprintf("%s_prop_%dm.tif", nm, win_cells * res_m), sprintf("%dm", res_m))
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
