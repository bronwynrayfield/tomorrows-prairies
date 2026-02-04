library(terra)

# colours (same scheme you set earlier)
col_to <- c(
  crop    = "#A07040",  # brown
  builtup  = "#800080",  # purple
  past    = "#2E8B57"
)

# helper: build a categorical SpatRaster of destination group for a given decade
make_decade_layer <- function(r_class, r_year, start, end_exclusive) {
  # year mask for the decade (e.g., 1990–1999 with end_exclusive=2000)
  r_year_cal <- r_year + 1900
  decade_mask <- r_year_cal >= start & r_year_cal < end_exclusive
  
  # destination class code = to_code
  to_code <- r_class %% 100
  
  # map to_code -> destination group (cropland, builtup, past, else NA)
  # 3=crop, 6=builtup, 4=upast, 10=ipast
  to_group <- classify(
    to_code,
    rcl = matrix(
      c(
        3, 1,   # cropland → 1
        6, 2,   # builtup  → 2
        4, 3,   # upast    → 3 (past)
        10, 3   # ipast    → 3 (past)
      ),
      ncol = 2, byrow = TRUE
    ),
    others = NA
  )
  
  # keep only pixels that changed in the decade
  out <- mask(to_group, decade_mask, maskvalues = 0, updatevalue = NA)
  
  # attach labels so legend shows names
  levels(out) <- data.frame(ID = c(1,2,3),
                            to_group = c("cropland","builtup","past"))
  
  out
}

r_1990s <- make_decade_layer(grass_loss_class, grass_loss_year, 1990, 2000)
r_2000s <- make_decade_layer(grass_loss_class, grass_loss_year, 2000, 2010)
r_2010s <- make_decade_layer(grass_loss_class, grass_loss_year, 2010, 2020)


# ensure palette order matches level IDs: 1=cropland, 2=builtup, 3=past
pal <- col_to[c("crop","builtup","past")]

par(mfrow = c(1,1), mar = c(3,3,2,5))

plot(r_1990s, col = pal, axes = FALSE, plg=list(title="Converted to"),
     main = "Grass loss destinations (1990–2000)")
plot(r_2000s, col = pal, axes = FALSE, plg=list(title="Converted to"),
     main = "Grass loss destinations (2000–2010)")
plot(r_2010s, col = pal, axes = FALSE, plg=list(title="Converted to"),
     main = "Grass loss destinations (2010–2020)")


plot_decade_with_outline <- function(r_decade, title, pal, pilot_v,
                                     outline_col = "black", outline_lwd = 2) {
  plot(r_decade, col = pal[c("crop","builtup","past")],
       axes = FALSE, plg = list(title = "Converted to"),
       main = title)
  lines(pilot_v, col = outline_col, lwd = outline_lwd)
}

plot_decade_with_outline(r_1990s, "Grass loss destinations (1990–2000)", pal, pilot)

par(mfrow = c(1,1), mar = c(3,3,2,5))
plot_decade_with_outline(r_1990s, "Grass loss destinations (1990–2000)", pal, pilot)
plot_decade_with_outline(r_2000s, "Grass loss destinations (2000–2010)", pal, pilot)
plot_decade_with_outline(r_2010s, "Grass loss destinations (2010–2020)", pal, pilot)



bg <- mask(grass_loss_class, pilot)       # same extent/res as your rasters
bg[] <- 1                          # fill all pilot pixels with 1

# Define a helper to plot with a grey background
plot_decade_with_outline_bg <- function(r_decade, title, pal, pilot_v,
                                        outline_col = "black", outline_lwd = 2) {
  # Plot light grey background
  plot(bg, col = "grey90", legend = FALSE, axes = FALSE, main = title)
  # Overlay colored change pixels
  plot(r_decade, add = TRUE,
       col = pal[c("cropland","builtup","past")],
       plg = list(title = "Converted to"))
  # Outline the pilot area
  lines(pilot_v, col = outline_col, lwd = outline_lwd)
}

plot_decade_with_outline_bg(r_1990s, "Grass loss destinations (1990–2000)", pal, pilot)
plot_decade_with_outline_bg(r_2000s, "Grass loss destinations (2000–2010)", pal, pilot)
plot_decade_with_outline_bg(r_2010s, "Grass loss destinations (2010–2020)", pal, pilot)

