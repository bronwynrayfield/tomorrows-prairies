##############################################
# Prepare predictor rasters

#--- Setup ----
source(file.path("./scripts", "0-setup.R"))
source(file.path("./scripts/pilot", "functions.R"))

# Load files 
pilot <- st_read(file.path("./processeddatasets", pilot_filepath, "GNP_pilot_area_crop_CSD.shp"))
lulc_pilot_pad <- rast(file.path("./processeddatasets", pilot_filepath, "AAFC2020_40kmbuffer_noUS_pilot_pad.tif"))
gnp <- st_read(file.path("./processeddatasets", pilot_filepath, "gnp_park_pilot.shp"))

#--- Map AAFC 2020 LULC in pilot area ----
lulc <- crop(lulc_pilot_pad, pilot)
lulc <- mask(lulc, pilot)

# 1) Category table
ct <- cats(lulc)[[1]]  # columns: Value, LU2020, Definition, IPCC_Class

# 2) Value sets
v_grass <- ct %>% filter(str_detect(LU2020, "Grassland$")) %>% pull(Value)

v_annual <- ct %>% 
  filter(LU2020 %in% c("Annual Cropland", "Land Converted to Annual Cropland", "Cropland")) %>% 
  pull(Value)

# Settlement EXCLUDING Roads
v_settlement <- ct %>% 
  filter(IPCC_Class == "Settlement", LU2020 != "Roads") %>% 
  pull(Value)

# 3) Reclass table (default → Other = 4)
# Codes: 1=Grassland, 2=Annual cropland, 3=Settlement, 4=Other
rcl <- data.frame(from = ct$Value, to = 4L)
rcl$to[rcl$from %in% v_grass]      <- 1L
rcl$to[rcl$from %in% v_annual]     <- 2L
rcl$to[rcl$from %in% v_settlement] <- 3L

lulc_4 <- classify(lulc, rcl)

# 4) Attach readable labels
levels(lulc_4) <- data.frame(
  value = 1:4,
  Class = c("Grassland", "Annual cropland", "Settlement", "Other")
)

pal <- c(
  "Grassland"       = "#5AAE61",
  "Annual cropland" = "#FEA500",
  "Settlement"      = "#7B3294",
  "Other"           = "#BDBDBD"
)

lv   <- levels(lulc_4)[[1]]
cols <- pal[lv$Class]

# Save to disk 
png("./results/landcover/pilot/aafc2020_4classes.png", width = 2000, height = 1600, res = 300)

plot(lulc_4, col = cols, axes = FALSE, legend = FALSE)
lines(terra::vect(pilot), lwd = 0.8, col = "red") # Pilot outline
plot(st_geometry(gnp), add = TRUE, border = "black", lwd = 0.8)

# Label: use st_coordinates() to get centroid
gnp_cent <- sf::st_coordinates(sf::st_centroid(gnp))

text(gnp_cent[1], gnp_cent[2], 
     labels = "Grasslands National Park", 
     cex = 0.5, font = 1)

dev.off()

