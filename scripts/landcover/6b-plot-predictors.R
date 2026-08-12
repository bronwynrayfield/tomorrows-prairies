# 6b-plot-predictors.R
# Paged multi-panel maps of the retained predictors (supplementary figure,
# styled after Fig. S2 in cobi.14346: grey panel strips, grey background outside
# the displayed population). Reads the retained set written by
# 6a-screen-collinearity-predictors.R so the figure stays in sync with the screen.
#
# Colour scales are continuous, clipped at clip_q and squished beyond, rather
# than quantile classes: with quantile breaks, colour distance does not
# correspond to value distance, so a near-constant variable looks as varied as a
# genuinely variable one.
source(file.path("scripts/landcover", "0a-setup.R"))
source(file.path("scripts/landcover", "0b-predictor-files.R"))

# ---- inputs ----
# population to display:
#   "valid"    - predictor-valid cells only (the supplementary figure)
#   "training" - baseline grassland within the valid mask (internal)
#   "sa"       - every study-area cell, valid or not. NAs render as the grey
#                panel background, so this version SHOWS the coverage gaps --
#                notably the SKSID miscellaneous units (hillwash, meadow,
#                runway, wetland) that carry no soil layer record and leave the
#                DSS covariates NA along SK waterways. For discussing
#                missingness, not for the supplement.
# Colour scale limits always come from cells that have data in the chosen
# population, so the legend describes the values actually shown.
population      <- "sa"

in_sa           <- file.path("./processeddatasets", "predictors", "30m")
out_dir         <- "output/landcover/figures"
plot_stack_file <- file.path(in_sa, "predictor-plot-stack_300m.tif")
cover_stack_file <- file.path(in_sa, "predictor-coverage-stack_300m.tif")
retained_file   <- file.path(in_sa, "retained-covariates.rds")
mask_file       <- file.path(in_sa, "predictor-valid-mask_study-area.tif")
sa_file         <- file.path(sa_filepath, sa_raster_filename)
grass_file      <- file.path("./processeddatasets", aafc_filepath, "binary",
                             "AAFC2010_ACI2009-2012_grass.tif")
break_n  <- 20000        # cells sampled at 30 m to derive scale limits
agg_fact <- 10           # 30 m -> 300 m plotting resolution
cover_th <- 0.5          # population = "sa" only: a 300 m cell is shown as
# missing unless at least this fraction of its 30 m
# cells have data
clip_q   <- c(.02, .98)  # percentile clip for colour scale limits; values
# outside are squished to the end colours so a few
# outliers don't wash out the map
ncol_pg  <- 3            # panels per row
nrow_pg  <- 4            # rows per page -> 12 panels/page
pg_w     <- 7.5          # page width (in)
pg_h     <- 10           # page height (in)

# any extra layers to include beyond the retained set (e.g. landcover /
# distance covariates); add name = path entries and a matching entry in
# predictor_labels in 0b-predictor-files.R
extra_files <- c(
  # crop_prop = tagged_path("./processeddatasets", ..., "modelling-extent"),
  # dist_crop = tagged_path("./processeddatasets", ..., "modelling-extent")
)

all_files <- c(predictor_files, extra_files)

# panel order for the figure: climate, then topography (slope last), then soil
# with the DSS-derived variables first. Kept here rather than reordering
# predictor_files, so the stack layer order used by the model is untouched.
# Names not in the retained set are dropped; any retained name missing from
# this vector is appended at the end rather than silently lost.
panel_order <- c(
  # climate
  "prtot", "tx_mean", "gddgrow_5",
  # topography
  "tri", "tpi", "srr", "northness", "eastness", "slope",
  # soil -- DSS first
  "drainage", "kp0", "kp33", "kp1500", "ldepth", "ldepth_max", "watertbl",
  # soil -- interpolated national products and moisture
  "oc", "bd", "clay", "sand", "ph",
  "psssm_mean", "psssm_sd", "tsm_planting", "tsm_growing"
)

final_set <- readRDS(retained_file)
final_set <- intersect(names(all_files), c(final_set, names(extra_files)))
final_set <- c(intersect(panel_order, final_set), setdiff(final_set, panel_order))

# ---- display mask ----
# For "sa" the mask is the study area itself, so cells that are NA in a
# predictor stay in the frame and render as background -- that is the point of
# this version. For the other two the mask is the predictor-valid intersection.
if (population == "sa") {
  msk <- rast(sa_file)
  msk <- !is.na(msk)
} else {
  msk <- rast(mask_file)
}
if (population == "training") {
  grass <- rast(grass_file)
  grass <- crop(grass, msk, snap = "out")        # grass is on the ME grid
  msk   <- msk & (subst(grass, NA, 0) > 0)
}
msk <- subst(msk, FALSE, NA)   # 0-cells must be NA, not FALSE, to exclude them

# ---- break sample (30 m, same population as the maps) ----
# Manual oversampling rather than na.rm = TRUE: the eligible fraction is low
# enough that spatSample's internal oversampling under-delivers and warns.
# Column 3 (after x and y) is the mask value. Subsample randomly, not head(),
# because spatSample returns cells in raster order.
set.seed(42)
pts <- spatSample(msk, break_n * 20, method = "random",
                  xy = TRUE, values = TRUE, na.rm = FALSE)
pts <- as.matrix(pts[!is.na(pts[[3]]), c("x", "y")])

rs <- rast(all_files[final_set])
names(rs) <- final_set   # rast() takes names from the files' internal band
# names, which are neither unique nor meaningful

smp <- terra::extract(rs, pts)
smp$ID <- NULL
# Drop only all-NA rows, not incomplete ones: for population = "sa" many cells
# are NA in the DSS variables but valid elsewhere, and requiring complete cases
# would compute every scale limit from the fully-surveyed subset. Limits are
# taken per variable with na.rm, so n varies by panel.
smp <- smp[rowSums(!is.na(smp)) > 0, , drop = FALSE]
smp <- smp[sample(nrow(smp), min(break_n, nrow(smp))), ]

stopifnot(nrow(smp) > 0, !length(setdiff(final_set, names(smp))))

# ---- plotting stack (aggregate once, cached) ----
# Built from rs, whose layer order is final_set, so the cache is only valid for
# that exact order. Rebuild whenever the count differs (e.g. a variable was
# dropped) -- and delete the cache manually if the ORDER changed but the count
# did not, since names are assigned positionally below.
cache_ok <- file.exists(plot_stack_file) &&
  nlyr(rast(plot_stack_file)) == length(final_set)

if (!cache_ok) {
  if (file.exists(plot_stack_file)) file.remove(plot_stack_file)
  aggregate(rs, fact = agg_fact, fun = "mean", na.rm = TRUE,
            filename = plot_stack_file,
            wopt = list(gdal = "COMPRESS=DEFLATE", datatype = "FLT4S"))
}
pr <- rast(plot_stack_file)
names(pr) <- final_set
pr <- pr[[final_set]]

# For "sa" the 300 m mean is computed with na.rm = TRUE, so a coarse cell keeps
# a value if ANY of its 100 fine cells is valid -- which would fill in exactly
# the narrow riparian gaps this version exists to show. Compute the valid
# fraction per coarse cell and blank cells below cover_th.
if (population == "sa") {
  if (!file.exists(cover_stack_file) ||
      nlyr(rast(cover_stack_file)) != length(final_set)) {
    if (file.exists(cover_stack_file)) file.remove(cover_stack_file)
    aggregate(!is.na(rs), fact = agg_fact, fun = "mean", na.rm = FALSE,
              filename = cover_stack_file,
              wopt = list(gdal = "COMPRESS=DEFLATE", datatype = "FLT4S"))
  }
  cov <- rast(cover_stack_file)
  names(cov) <- final_set
  pr <- mask(pr, cov[[final_set]] >= cover_th, maskvalues = c(0, NA))
}

# crop/mask to the display population (cache is at modelling extent; the two
# 300 m grids don't share an origin, so snap out then resample the mask on)
msk300 <- aggregate(msk, agg_fact, fun = "max", na.rm = TRUE)
pr     <- crop(pr, msk300, snap = "out")
msk300 <- resample(msk300, pr, method = "near")
pr     <- mask(pr, msk300)

# ---- panel builder ----
fmt <- function(x) formatC(signif(x, 3), format = "fg", big.mark = ",")

# manual scale limits where the percentile clip is unhelpful -- mainly signed
# variables that read better with symmetric limits around zero
limits_override <- list(
  # tpi       = c(-5, 5),
  # northness = c(-1, 1),
  # eastness  = c(-1, 1)
)

make_panel <- function(var) {
  lim <- limits_override[[var]]
  if (is.null(lim)) lim <- unname(quantile(smp[[var]], clip_q, na.rm = TRUE))
  if (length(lim) != 2 || !all(is.finite(lim)) || diff(lim) <= 0) {
    lim <- range(smp[[var]], na.rm = TRUE, finite = TRUE)       # near-constant
  }
  if (length(lim) != 2 || !all(is.finite(lim)) || diff(lim) <= 0) {
    lim <- if (length(lim) == 2 && is.finite(lim[1])) lim[1] + c(-1, 1) else c(0, 1)
  }
  
  rl <- pr[[var]]
  names(rl) <- predictor_labels[var]
  
  ggplot() +
    geom_spatraster(data = rl, maxcell = 1e6) +
    facet_wrap(~lyr) +
    scale_fill_viridis_c(
      limits = lim, oob = scales::squish, na.value = NA, name = NULL,
      breaks = lim, labels = fmt(lim),
      guide = guide_colourbar(barwidth = unit(42, "pt"), barheight = unit(5, "pt"),
                              direction = "horizontal", ticks = FALSE)
    ) +
    theme_void() +
    theme(
      strip.background  = element_rect(fill = "grey75", colour = "grey30"),
      strip.text        = element_text(size = 8, margin = margin(2, 2, 2, 2)),
      panel.background  = element_rect(fill = "grey92", colour = "grey30"),
      legend.position   = "inside",
      legend.position.inside = c(0.9, 0.97),
      legend.justification   = c(1, 1),
      legend.key.size   = unit(7, "pt"),
      legend.text       = element_text(size = 5),
      legend.margin     = margin(3, 8, 3, 8),
      legend.background = element_rect(fill = alpha("white", 0.7), colour = NA),
      plot.margin       = margin(1, 1, 1, 1)
    )
}

panels <- lapply(final_set, make_panel)
names(panels) <- final_set

# ---- assemble and write pages ----
pages <- split(final_set, ceiling(seq_along(final_set) / (ncol_pg * nrow_pg)))

for (i in seq_along(pages)) {
  pg <- wrap_plots(panels[pages[[i]]], ncol = ncol_pg, nrow = nrow_pg)
  ggsave(file.path(out_dir, sprintf("predictor-maps_%s_page%d.png", population, i)),
         pg, width = pg_w, height = pg_h, dpi = 300, bg = "white")
}