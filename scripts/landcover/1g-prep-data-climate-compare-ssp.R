# =============================================================================
# Purpose: Prep climate data and compare SSPs
# =============================================================================

library(terra)
dir  <- "D:/tomorrows-prairies/existingdatasets/Climate/CanDCS-M6/geotiff"
vars <- c("tx_mean","tn_mean","prtot","gddgrow_5","hddheat_18")
fut  <- "2021-2050"
base <- "1991-2020"        # historic baseline for context

f  <- function(v, s, h) rast(file.path(dir, sprintf("%s_%s_%s_p50.tif", v, s, h)))
gm <- function(r) global(r, "mean", na.rm = TRUE)[1, 1]

# scan all files: layers + value range
fs  <- list.files(dir, pattern = "\\.tif$", full.names = TRUE)
for (f in fs) {
  r  <- rast(f)
  mm <- minmax(r, compute = TRUE)
  cat(sprintf("%-48s nlyr=%d  range=[%.1f, %.1f]\n",
              basename(f), nlyr(r), mm[1], mm[2]))
}

# load one for a visual + metadata check
r <- rast(file.path(dir, "tx_mean_ssp245_2021-2050_p50.tif"))
print(r)                 # CRS, resolution, extent, nlyr
global(r, "notNA")       # non-NA cell count
plot(r)                  # eyeball

tab <- do.call(rbind, lapply(vars, function(v) {
  a   <- f(v, "ssp245", fut)            # SSP2-4.5 future
  b   <- f(v, "ssp585", fut)            # SSP5-8.5 future
  h   <- f(v, "ssp245", base)           # baseline
  d   <- b - a                          # inter-scenario difference
  chg <- a - h                          # SSP2-4.5 projected change
  data.frame(
    variable          = v,
    mean_ssp245       = round(gm(a), 2),
    mean_ssp585       = round(gm(b), 2),
    mean_diff         = round(gm(d), 3),
    max_abs_diff      = round(max(abs(minmax(d, compute = TRUE))), 3),
    proj_change_245   = round(gm(chg), 2)
  )
}))

print(tab, row.names = FALSE)
write.csv(tab, file.path(dir, "ssp_comparison_2021-2050.csv"), row.names = FALSE)