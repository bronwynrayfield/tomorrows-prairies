# install.packages(c("maptiles","terra","sf"))  # if needed
library(maptiles)
library(terra)
library(sf)

# 0) Prep objects ------------------------------------------------------------
# pilot: sf or SpatVector; r_class/r_year already in your equal-area CRS
pilot_v <- if (inherits(pilot, "SpatVector")) pilot else vect(pilot)
if (!identical(crs(pilot_v), crs(r_class))) pilot_v <- project(pilot_v, r_class)

# 1) Fetch satellite basemap over pilot extent (Web Mercator), then reproject --
# Provider options: "Esri.WorldImagery" (satellite), "CartoDB.PositronNoLabels" (clean grayscale streets)
bm_merc <- get_tiles(pilot_v, provider = "Esri.WorldImagery", zoom = 11)  # raise zoom for more detail
# Reproject to your raster CRS so it lines up:
bm <- project(bm_merc, r_class)

# 2) Convert RGB basemap to greyscale (desaturate) ---------------------------
# maptiles returns an RGB SpatRaster (3 bands); make 1-band gray
to_gray <- function(r) 0.299*r[[1]] + 0.587*r[[2]] + 0.114*r[[3]]
bm_gray <- to_gray(bm)

# 3) Your decade layer function (from earlier) -------------------------------
# pal must match classify IDs: 1=cropland, 2=builtup, 3=past
pal <- c(cropland = "#8B4513", builtup = "#800080", past = "#2E8B57")

plot_decade_with_basemap <- function(r_decade, title, pal, pilot_v) {
  # Basemap first (greyscale)
  plot(bm_gray, col = gray.colors(256), legend = FALSE, axes = FALSE, main = title)
  
  # Decade categories on top (no legend to avoid NA legend issues per decade)
  plot(r_decade, add = TRUE, col = pal[c("cropland","builtup","past")], legend = FALSE)
  
  # Pilot outline (two-pass for contrast)
  lines(pilot_v, col = "white", lwd = 3.5)
  lines(pilot_v, col = "black", lwd = 1.5)
  
  # Static legend (stable across decades)
  legend("right", inset = 0.02, title = "Converted to",
         fill = pal[c("cropland","builtup","past")],
         legend = c("cropland","builtup","past"),
         bty = "n", cex = 0.9)
}

# 4) Draw the three maps -----------------------------------------------------
par(mfrow = c(1,1), mar = c(3,3,2,5))
plot_decade_with_basemap(r_1990s, "Grass loss destinations (1990–2000)", pal, pilot_v)
plot_decade_with_basemap(r_2000s, "Grass loss destinations (2000–2010)", pal, pilot_v)
plot_decade_with_basemap(r_2010s, "Grass loss destinations (2010–2020)", pal, pilot_v)

# Optional: save to PNGs
png("map_grass_loss_1990_2000.png", width=1600, height=1200, res=200)
plot_decade_with_basemap(r_1990s, "Grass loss destinations (1990–2000)", pal, pilot_v)
dev.off()
png("map_grass_loss_2000_2010.png", width=1600, height=1200, res=200)
plot_decade_with_basemap(r_2000s, "Grass loss destinations (2000–2010)", pal, pilot_v)
dev.off()
png("map_grass_loss_2010_2020.png", width=1600, height=1200, res=200)
plot_decade_with_basemap(r_2010s, "Grass loss destinations (2010–2020)", pal, pilot_v)
dev.off()
