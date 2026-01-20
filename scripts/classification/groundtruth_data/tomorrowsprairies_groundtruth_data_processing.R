
# Tomorrows Prairies table summaries - Groundtruth Data

# add in lib
library(dplyr)
library(sf)
library(tools)

# set dir & load datasets
setwd('C:/Users/samantha.yeo/OneDrive - The Nature Conservancy/Documents/Projects/NCS/Canada/TomorrowsPrairies/Data')
gdb_path     <- "./04_study_area/groundtruth_points_grass/Final_QC_SK_TP_202509/Final_QC_SK_TP_202509/Final_QC_SK_TP_202509.gdb"

tp_studyarea <- st_read("./04_study_area/TP_study_area_prairies_ecozone/TP_study_area_prairies_ecozone.shp")
gt_points    <- st_read("./04_study_area/groundtruth_points_grass/groundtruth_points_grass.shp")
ecoR         <- st_read("./Ecoregions/ecoregions.shp")

# load provinces
base <- "https://services.arcgis.com/wjcPoefzjpzCgffS/arcgis/rest/services/Provinces_and_Territories_of_Canada/FeatureServer/0/query"
query_url <- paste0(base,"?where=1%3D1", "&outFields=*", "&f=geojson") 
provinces <- read_sf(query_url)

# clean data
# Québec
unique(provinces$Name_EN)
provinces$Name_EN[provinces$Name_EN == "Québec"] <- "Quebec"
unique(provinces$Name_EN)

# class names
unique(gt_points$LC_Class_O)
gt_points$LC_Class_O <- toTitleCase(tolower(gt_points$LC_Class_O))
gt_points$LC_Class_O[gt_points$LC_Class_O == "Clover - Red Clover"] <- "Clover - Red"
gt_points$LC_Class_O[gt_points$LC_Class_O == "Grass"] <- "Grassland"
gt_points$LC_Class_O[gt_points$LC_Class_O == "Pasture - Unimproved"] <- "Unimproved Pasture"
unique(gt_points$LC_Class_O)

####################################
# intersect with study zone 
# transform CRS
if (st_crs(gt_points) != st_crs(tp_studyarea)) {  gt_points <- st_transform(gt_points, st_crs(tp_studyarea))}

# intersect
hits_mat <- st_intersects(gt_points, tp_studyarea, sparse = F)
inside <- rowSums(hits_mat) > 0
gt_points$tp_studyarea <- as.integer(inside)  # creates 1/0

# Count points
inside_n  <- sum(gt_points$tp_studyarea == 1, na.rm = TRUE) ; print(inside_n)
outside_n <- sum(gt_points$tp_studyarea == 0, na.rm = TRUE) ; print(outside_n)
paste0(round(inside_n/outside_n*100, 2), '% are inside the TP Study Area')

# verif
pts_out <- gt_points[gt_points$tp_studyarea == 0, ]
bb <- st_bbox(pts_out)
zoom_factor <- 0.8 # smaller = closer
x_range <- (bb$xmax - bb$xmin) * zoom_factor / 2
y_range <- (bb$ymax - bb$ymin) * zoom_factor / 2
xmid <- (bb$xmax + bb$xmin) / 2
ymid <- (bb$ymax + bb$ymin) / 2
  
plot(st_geometry(gt_points[gt_points$tp_studyarea == 0, ]), pch = 20, col = "#2FFF99", cex = 0.9, 
     main = "Ground-truth Points vs TP Study Area", add = F, 
     xlim = c(xmid - x_range, xmid + x_range),
     ylim = c(ymid - y_range, ymid + y_range))
plot(st_geometry(gt_points[gt_points$tp_studyarea == 1, ]), pch = 20, col = "#2C7BB6", cex = 0.9, add = T)
plot(st_geometry(tp_studyarea), col = NA, border = "gray40", lwd = 1.5, add = T)

# Add legend on the right side
legend("right",
  legend = c(
    paste0("Inside (", inside_n, ")"), 
    paste0("Outside (", outside_n, ")")),
  col = c("#2C7BB6", "#2FFF99"), pch = 20, pt.cex = 1.0, bty = "n", inset = c(0, 0))


####################################
# intersect with provinces
if (st_crs(gt_points) != st_crs(provinces)) {  gt_points <- st_transform(gt_points, st_crs(provinces))}
# if (st_crs(gt_points) != st_crs(provinces)) {  provinces <- st_transform(provinces, st_crs(gt_points))}

# intersect
gt_points <- st_join(
  gt_points,
  provinces[, c("Name_EN")],   # keep only the province name column + geometry
  join = st_intersects
)

names(gt_points)[names(gt_points) == "Name_EN"] <- "Province"
sum(is.na(gt_points$Province)) # 28 points are NA...

# verif
bb <- st_bbox(gt_points)
plot(st_geometry(provinces),
  col = NA, border = "gray50",
  xlim = c(bb["xmin"], bb["xmax"]),
  ylim = c(bb["ymin"], bb["ymax"]),
  main = "Groundtruth Points by Province")

# Unique colors for provinces
prov_names <- unique(gt_points$Province)
prov_colors <- rainbow(sum(!is.na(prov_names)))
names(prov_colors) <- prov_names[!is.na(prov_names)]

# Plot points colored by province
plot(st_geometry(gt_points[!is.na(gt_points$Province), ]),
  col = prov_colors[gt_points$Province[!is.na(gt_points$Province)]],
  pch = 20, cex = 0.8, add = TRUE)
# Add in the NA vals
plot(st_geometry(gt_points[is.na(gt_points$Province), ]),
  col = "black",
  pch = 4, cex = 1.0, lwd = 1.2, add = TRUE)

# # Legend
# legend("topleft", legend = names(prov_colors), col = prov_colors, pch = 20, cex = 0.8, bty = "n")
# legend("topleft", legend = c(names(prov_colors), "NA"), col = c(prov_colors, "black"), pch = c(rep(20, length(prov_colors)), 4), pt.cex = 0.9, bty = "n")


####################################
# intersect with ecoR
if (st_crs(gt_points) != st_crs(ecoR)) {  ecoR <- st_transform(ecoR, st_crs(gt_points))}

# intersect both ecoregion and region_nam to df
ecoR_sel <- ecoR[, c("ECOREGION", "REGION_NAM")]
gt_points <- st_join(gt_points,
  ecoR_sel,
  join = st_intersects)

# Rename to requested column name
names(gt_points)[names(gt_points) == "ECOREGION"] <- "Ecoregion_ID"
names(gt_points)[names(gt_points) == "REGION_NAM"] <- "Ecoregion_name"

# NA vals
unique(gt_points$Ecoregion_name)
sum(is.na(gt_points$Ecoregion_name)) # 6 NA vals

# Draw ecoregions
bb <- st_bbox(gt_points)
plot(st_geometry(ecoR),
     col = NA, border = "gray50",
     xlim = c(bb["xmin"], bb["xmax"]),
     ylim = c(bb["ymin"], bb["ymax"]),
     main = "Ground-truth Points by Ecoregion")

# Define a color per ecoregion
ecoreg_levels <- sort(unique(gt_points$Ecoregion_name))
ecoreg_palette <- setNames(rainbow(length(ecoreg_levels)), ecoreg_levels)

# Add points, colored by ecoregion
plot(st_geometry(gt_points),
     col = ecoreg_palette[gt_points$Ecoregion_name],
     pch = 20, cex = 0.8, add = TRUE)


#############################
# write shapefile
st_write(
  gt_points,
  "./Output_Data/gt_points_modified.shp",
  delete_dsn = T    # overwrite if it already exists
)

###############################
# summarise data into tables 

# summarise by - LC_Class_Cl, LC_Class_O,	Province
dat_summary <- gt_points %>%
  st_drop_geometry(dat_summary) %>%
  group_by(Province, LC_Class_N, LC_Class_O) %>% 
  summarise(count = n(), .groups = "drop")
dat_summary$perc = round(dat_summary$count/sum(dat_summary$count)*100, 2)

# summarise by - LC_Class_Cl,	Province
dat_summary <- gt_points %>%
  st_drop_geometry(dat_summary) %>%
  group_by(Province, LC_Class_N) %>% 
  summarise(count = n(), .groups = "drop")
dat_summary$perc = round(dat_summary$count/sum(dat_summary$count)*100, 2)

# summarise by - EcoR,	Province
dat_summary <- gt_points %>%
  st_drop_geometry(dat_summary) %>%
  group_by(Province, Ecoregion_name) %>% 
  summarise(count = n(), .groups = "drop")
dat_summary$perc = round(dat_summary$count/sum(dat_summary$count)*100, 2)



