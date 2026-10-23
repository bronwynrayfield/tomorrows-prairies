# =============================================================================
# Purpose:
# Visualize predictors for QC/QC and communication
#
# Paged multi-panel maps of the retained predictors (supplementary figure,
# styled after Fig. S2 in cobi.14346: grey panel strips, grey background outside
# the displayed population). Reads the retained set written by
# 6a-screen-collinearity-predictors.R so the figure stays in sync with the screen.
#
# Colour scales are continuous, clipped at clip_q and squished beyond, rather
# than quantile classes: with quantile breaks, colour distance does not
# correspond to value distance, so a near-constant variable looks as varied as a
# genuinely variable one.
# =============================================================================
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
# panel groups = figure sections; order within each group is panel order
panel_groups <- list(
  climatic       = c("prtot", "tx_mean", "gddgrow_5", "nrcdd"),
  edaphic        = c("drainage", "ph", "kp0", "kp33", "kp1500", "ldepth",
                     "oc", "bd", "clay", "psssm_mean", "psssm_sd"),
  topographic    = c("elev", "slope", "tri", "tpi", "srr", "northness", "eastness"),
  landscape      = c("native", "tame", "unprot_grs", "crop", "settl", "fors"),
  socio_economic = c("net_return", "dist_to_elevator", "dist_to_road", "dist_to_rail", "dist_to_popctr", "dist_to_regctr", "pop_density")
)
panel_order <- unlist(panel_groups, use.names = FALSE)

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
# Sample points are drawn once and cached, so every variable is sampled at the same cells
# and a single variable can be re-extracted later. Columns are refreshed when missing or
# older than their input file.
pts_file <- file.path(in_sa, sprintf("break-points_%s.rds", population))
if (file.exists(pts_file)) {
  pts <- readRDS(pts_file)
} else {
  set.seed(42)
  pts <- spatSample(msk, break_n * 20, method = "random", xy = TRUE, values = TRUE, na.rm = FALSE)
  pts <- as.matrix(pts[!is.na(pts[[3]]), c("x", "y")])
  pts <- pts[sample(nrow(pts), min(nrow(pts), break_n)), ]
  saveRDS(pts, pts_file)
}

sa_path <- function(nm) file.path(in_sa, paste0(nm, "_study-area.tif"))
rs <- rast(sapply(final_set, sa_path))
names(rs) <- final_set

smp_file <- file.path(in_sa, sprintf("break-sample_%s.rds", population))
smp  <- if (file.exists(smp_file)) readRDS(smp_file) else data.frame(row.names = seq_len(nrow(pts)))
todo <- if (file.exists(smp_file)) final_set[!final_set %in% names(smp) |
                                               file.mtime(sapply(final_set, sa_path)) > file.mtime(smp_file)] else final_set
if (length(todo)) {
  smp[todo] <- terra::extract(rs[[todo]], pts)[todo]
  saveRDS(smp, smp_file)
}
smp <- smp[final_set]   # limits use na.rm per variable, so all-NA rows need no dropping

stopifnot(nrow(smp) > 0, !length(setdiff(final_set, names(smp))))

# ---- plotting stack (aggregate once, cached) ----
# Built from rs, whose layer order is final_set, so the cache is only valid for
# that exact order. Rebuild whenever the count differs (e.g. a variable was
# dropped) -- and delete the cache manually if the ORDER changed but the count
# did not, since names are assigned positionally below.
# one 300 m file per variable, rebuilt only if missing or older than its input
agg_dir <- file.path(in_sa, "plot-300m"); dir.create(agg_dir, showWarnings = FALSE)
pr_file <- file.path(agg_dir, sprintf("plot-stack_%s.tif", population))
d_vars <- grep("^dist_", final_set, value = TRUE)

if (!file.exists(pr_file)) {
  agg_path <- function(nm, kind = "mean") file.path(agg_dir, sprintf("%s_%s_300m.tif", nm, kind))
  stale <- function(out, nm) !file.exists(out) || file.mtime(out) < file.mtime(sa_path(nm))
  for (nm in final_set) if (stale(agg_path(nm), nm))
    aggregate(rs[[nm]], fact = agg_fact, fun = if (nm == "drainage") "modal" else "mean", na.rm = TRUE, filename = agg_path(nm),
              overwrite = TRUE, wopt = list(gdal = "COMPRESS=DEFLATE", datatype = "FLT4S"))
  pr <- rast(sapply(final_set, agg_path))
  names(pr) <- final_set
  
  # For "sa" the 300 m mean is computed with na.rm = TRUE, so a coarse cell keeps
  # a value if ANY of its 100 fine cells is valid -- which would fill in exactly
  # the narrow riparian gaps this version exists to show. Compute the valid
  # fraction per coarse cell and blank cells below cover_th.
  if (population == "sa") {
    for (nm in final_set) if (stale(agg_path(nm, "cover"), nm))
      aggregate(!is.na(rs[[nm]]), fact = agg_fact, fun = "mean", na.rm = FALSE, filename = agg_path(nm, "cover"),
                overwrite = TRUE, wopt = list(gdal = "COMPRESS=DEFLATE", datatype = "FLT4S"))
    cov <- rast(sapply(final_set, agg_path, kind = "cover"))
    names(cov) <- final_set
    pr <- mask(pr, cov[[final_set]] >= cover_th, maskvalues = c(0, NA))
  }
  
  # crop/mask to the display population (cache is at modelling extent; the two
  # 300 m grids don't share an origin, so snap out then resample the mask on)
  msk300 <- aggregate(msk, agg_fact, fun = "max", na.rm = TRUE)
  pr     <- crop(pr, msk300, snap = "out")
  msk300 <- resample(msk300, pr, method = "near")
  pr     <- mask(pr, msk300)
  
  # convert distance predictors from m to km for display
  pr[[d_vars]] <- pr[[d_vars]] / 1000

  # save  
  writeRaster(pr, pr_file, overwrite = TRUE)
}
pr <- rast(pr_file)
smp <- readRDS(smp_file)[final_set]; smp[d_vars] <- smp[d_vars] / 1000   # reload so reruns can't double-convert

# ---- panel builder ----
# decimals set by the span of the limits, so both ends share a precision (~2 sig figs on the range):
# pH 6.10–7.59 (span 1.49) -> 1 dp; 0.012–0.34 -> 2 dp; 1,500–4,500 -> 0 dp
fmt <- function(x) formatC(x, format = "f", digits = max(0, 1 - floor(log10(diff(range(x))))), big.mark = ",")

# manual scale limits where the percentile clip is unhelpful -- mainly signed
# variables that read better with symmetric limits around zero
limits_override <- list(
  drainage = c(1, 7)
  # tpi       = c(-5, 5),
  # northness = c(-1, 1),
  # eastness  = c(-1, 1)
)
labels_override <- list(
  drainage = c("VR", "R", "W", "MW", "I", "P", "VP")
  )

# scale limits shared by maps and histograms, so colour = same value in both figures
get_lim <- function(var) {
  lim <- limits_override[[var]]
  if (is.null(lim)) lim <- unname(quantile(smp[[var]], clip_q, na.rm = TRUE))
  if (length(lim) != 2 || !all(is.finite(lim)) || diff(lim) <= 0) {
    lim <- range(smp[[var]], na.rm = TRUE, finite = TRUE)       # near-constant
  }
  if (length(lim) != 2 || !all(is.finite(lim)) || diff(lim) <= 0) {
    lim <- if (length(lim) == 2 && is.finite(lim[1])) lim[1] + c(-1, 1) else c(0, 1)
  }
  lim
}

make_panel <- function(var) {
  lim <- get_lim(var)
  
  rl <- pr[[var]]
  names(rl) <- predictor_labels[var]
  
  ggplot() +
    geom_spatraster(data = rl, maxcell = 1e6) +
    facet_wrap(~lyr) +
    scale_fill_viridis_c(
      limits = lim, oob = scales::squish, na.value = NA, name = NULL,
      breaks = if (is.null(labels_override[[var]])) lim else seq(lim[1], lim[2]),
      labels = labels_override[[var]] %||% fmt(lim),
      guide = guide_colourbar(barwidth = unit(if (is.null(labels_override[[var]])) 42 else 70, "pt"), 
                              barheight = unit(5, "pt"),
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

# Distributions of the same 30 m sample used for the colour limits (n = break_n,
# minus NAs per variable). Bars are filled by value on the map's scale, so the
# squished tails beyond clip_q show as the end colours. Ordinal variables
# (labels_override) get one bar per class instead of binned counts.
log_vars <- c("pop_density")   # extreme right skew: 98th pct << max, linear bins collapse

make_hist <- function(var) {
  lim <- get_lim(var)
  d <- data.frame(x = smp[[var]][!is.na(smp[[var]])], lyr = predictor_labels[var])
  lab <- labels_override[[var]]
  p <- if (is.null(lab)) {
    ggplot(d, aes(x)) +
      geom_histogram(aes(fill = after_stat(if (var %in% log_vars) expm1(x) else x)), bins = 40) +
      scale_x_continuous(transform = if (var %in% log_vars) "log1p" else "identity",
                         breaks = if (var %in% log_vars) c(0, 10, 100, 1000, 10000) else waiver())
  } else {
    ggplot(d, aes(factor(x, levels = seq(lim[1], lim[2])), fill = x)) +
      geom_bar() + scale_x_discrete(labels = lab, drop = FALSE)
  }
  p +
    facet_wrap(~lyr) +
    scale_fill_viridis_c(limits = lim, oob = scales::squish, guide = "none") +
    theme_minimal(base_size = 5) +
    theme(
      strip.background = element_rect(fill = "grey75", colour = "grey30"),
      strip.text       = element_text(size = 8, margin = margin(2, 2, 2, 2)),
      panel.border     = element_rect(fill = NA, colour = "grey30"),
      panel.grid.minor = element_blank(),
      axis.title       = element_blank(),
      plot.margin      = margin(1, 1, 1, 1)
    )
}

hists <- lapply(final_set, make_hist)
names(hists) <- final_set

# ---- assemble and write pages, one figure per group ----
# retained names in no group go to "other" rather than being silently lost
groups <- lapply(panel_groups, intersect, final_set)
groups$other <- setdiff(final_set, panel_order)
groups <- groups[lengths(groups) > 0]

for (kind in c("maps", "hists")) {
  pl <- if (kind == "maps") panels else hists
  for (g in names(groups)) {
    pages <- split(groups[[g]], ceiling(seq_along(groups[[g]]) / (ncol_pg * nrow_pg)))
    for (i in seq_along(pages)) {
      nr <- ceiling(length(pages[[i]]) / ncol_pg)   # rows actually used; page height scales with it
      pg <- wrap_plots(pl[pages[[i]]], ncol = ncol_pg, nrow = nr)
      ggsave(file.path(out_dir, sprintf("predictor-%s_%s_%s_page%d.png", kind, population, g, i)),
             pg, width = pg_w, height = pg_h * nr / nrow_pg, dpi = 300, bg = "white")
    }
  }
}


# ---- summary table ----
# min/max/mean/n: exact, all 30 m cells in the display population (full pass, cached).
# P2/P50/P98: from smp (n = break_n); P2/P98 are the colour-scale limits in the figures.
stats_file <- file.path(in_sa, sprintf("predictor-stats_%s.rds", population))
if (!file.exists(stats_file)) {
  rs_pop <- if (population == "sa") rs else mask(rs, msk)
  saveRDS(global(rs_pop, c("min", "max", "mean", "notNA"), na.rm = TRUE), stats_file)
}
st <- readRDS(stats_file)
st[d_vars, c("min", "max", "mean")] <- st[d_vars, c("min", "max", "mean")] / 1000

q <- t(sapply(final_set, function(v) quantile(smp[[v]], c(.02, .5, .98), na.rm = TRUE)))
tab <- data.frame(st[final_set, c("min")], q, st[final_set, c("max", "mean")])
tab["drainage", "mean"] <- NA          # mean of ordinal class codes isn't meaningful
tab <- t(apply(tab, 1, function(x) fmt(x[!is.na(x)])[match(seq_along(x), which(!is.na(x)))]))

tab <- data.frame(Section   = rep(names(groups), lengths(groups))[match(final_set, unlist(groups))],
                  Predictor = predictor_labels[final_set],
                  N         = format(st[final_set, "notNA"], big.mark = ","),
                  tab, check.names = FALSE)
names(tab)[4:9] <- c("Min", "P2", "Median", "P98", "Max", "Mean")

ft <- flextable::flextable(tab) |> flextable::merge_v(j = "Section") |> flextable::autofit()
flextable::save_as_docx(ft, path = file.path(out_dir, sprintf("predictor-summary_%s.docx", population)))


