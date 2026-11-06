library(terra)

# 1) decode and mask to changed pixels only
to_code <- grass_loss_class %% 100
year_cal <- grass_loss_year + 1900
changed  <- year_cal >= 1985 & year_cal <= 2022
to_code_change <- mask(to_code, changed, maskvalues = 0, updatevalue = NA)

# 2) unique class codes in the order terra will plot (sorted, no NA)
vals <- sort(unique(values(to_code_change)))
vals <- vals[!is.na(vals)]

# 3) class labels (metadata)
class_lookup <- c(
  `1`="forest", `2`="shrub", `3`="crop", `4`="upast", `5`="grass",
  `6`="builtup", `7`="bare", `8`="wetland", `9`="water", `10`="ipast",
  `11`="forest/shrub", `12`="forest/wetland", `13`="shrub/forest",
  `14`="wetland/forest", `15`="upast/grass", `16`="grass/upast"
)
labs <- class_lookup[as.character(vals)]

# 4) define a **named palette by code**, then align it to 'vals'
pal_by_code <- c(
  `1`="green",  # forest
  `2`="#66A61E",  # shrub
  `3`="#A07040",  # crop  (brown)
  `4`="#2E8B57",  # upast
  `5`="#32CD32",  # grass
  `6`="#800080",  # builtup (purple)
  `7`="#808080",  # bare (grey)
  `8`="#1E90FF",  # wetland
  `9`="#4169E1",  # water
  `10`="#228B22", # ipast
  `11`="#6A3D9A",
  `12`="#8DA0CB",
  `13`="#B3B3B3",
  `14`="#FDBF6F",
  `15`="#66A61E",
  `16`="#FFFF99"
)
cols_vec <- unname(pal_by_code[as.character(vals)])  # <-- aligns colors to codes

# 5) light-grey background inside the raster extent
bg <- to_code; values(bg) <- 1

par(mar = c(3,3,2,6))

plot(bg, col = "grey90", legend = FALSE, axes = FALSE,
     main = "Grass loss — destination class (changed pixels only)")
plot(to_code_change, add = TRUE, col = cols_vec, legend = FALSE)

# add pilot boundary (convert sf → SpatVector for terra)
lines(vect(pilot), col = "white", lwd = 3)   # outer white stroke for contrast
lines(vect(pilot), col = "black", lwd = 1.2) # inner black outline

# add legend
legend("right", inset = 0.02, title = "Converted to",
       fill = cols_vec, legend = labs, bty = "n", cex = 0.9)




# change year
# 1. Decode the encoded years (year - 1900)
year_cal <- grass_loss_year + 1900
names(year_cal) <- "year"

# 2. Mask out no-change values (usually 0 → 1900)
year_cal[year_cal < 1985 | year_cal > 2022] <- NA

# 3. Choose a nice continuous color ramp for years
yr_range <- range(values(year_cal), na.rm = TRUE)
yr_colors <- hcl.colors(diff(yr_range) + 1, "Viridis")

# 4. Plot
par(mar = c(3, 3, 2, 6))
plot(year_cal,
     col = yr_colors,
     plg = list(title = "Year of change"),
     axes = FALSE,
     main = "Grass loss — year of most recent change")

# 5. Add pilot outline (two-pass for visibility)
lines(vect(pilot), col = "white", lwd = 3)
lines(vect(pilot), col = "black", lwd = 1.2)


