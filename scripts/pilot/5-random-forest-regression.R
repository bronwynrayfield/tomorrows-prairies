##############################################
# Conversion risk mapping using random forest regression

#--- Setup ----
source(file.path("./scripts", "0-setup.R"))
source(file.path("./scripts/pilot", "functions.R"))

# Load files 
pilot <- st_read(file.path(pilot_filepath, "GNP_pilot_area_crop_CSD.shp"))
gnp <- st_read(file.path(pilot_filepath, "gnp_park_pilot.shp"))
crop_suit_filled <- rast(file.path(pilot_filepath, "cropland_suitability_2019_noNAs.tif"))
dist_to_road <- rast(file.path(pilot_filepath, "dist_to_road.tif"))
cropland_prop <- rast(file.path(pilot_filepath, "cropland_proportion_3-2.tif"))
settle_prop <- rast(file.path(pilot_filepath, "settle_proportion_3-2.tif"))
cropland <- rast(file.path(pilot_filepath, "cropland_binary.tif"))
settle <- rast(file.path(pilot_filepath, "settle_binary.tif"))

# For plotting later, get centre coordinates of gnp
gnp_cent <- sf::st_coordinates(sf::st_centroid(gnp))


#--- Cropland regression model ----
# stack rasters (predictors)
pred_stack <- c(crop_suit_filled, cropland_prop, dist_to_road)
names(pred_stack) <- c("crop_suit", "cropland_prop", "dist_to_road")

# sample training data from rasters, INCLUDING coordinates
set.seed(123)
train_df <- terra::spatSample(
  pred_stack,
  size   = 5000,
  method = "random",
  na.rm  = TRUE,
  as.df  = TRUE,
  xy     = TRUE   # <-- this adds "x" and "y" columns
)

# add response variable by extracting from 'cropland' at those coords
vals <- terra::extract(cropland, train_df[, c("x", "y")])
train_df$cropland <- vals$LU2020

# drop any rows with NA in predictors/response
train_df <- na.omit(train_df)


# 70% training, 30% testing
train_id <- sample(1:nrow(train_df), 0.7 * nrow(train_df))
train_data <- train_df[train_id, ]
test_data  <- train_df[-train_id, ]

train_data$cropland <- factor(train_data$cropland, levels = c(0, 1))
test_data$cropland  <- factor(test_data$cropland,  levels = c(0, 1))


# Fit ranger model
rf_model <- ranger(
  cropland ~ crop_suit + cropland_prop + dist_to_road,
  data = train_data,
  num.trees = 500,
  mtry = 2,                 
  min.node.size = 10,
  importance = "impurity",
  probability = TRUE,
  seed = 123
)

print(rf_model)
rf_model$variable.importance

# # model performance - needs updating
# # extract probability of cropland = 1
# test_pred <- predict(rf_model, data = test_data)$predictions[, 2]
# 
# # predicted class using 0.5 threshold
# pred_class <- ifelse(test_pred > 0.5, 1, 0)
# 
# # accuracy
# acc <- mean(pred_class == as.numeric(as.character(test_data$cropland)))
# cat("Overall accuracy:", round(acc, 3), "\n")
# 
# # confusion matrix
# confusionMatrix(
#   factor(pred_class, levels = c(0, 1)),
#   test_data$cropland,
#   positive = "1"
# )
# 
# roc_obj <- roc(test_data$cropland, test_pred)
# auc_val <- auc(roc_obj)
# plot(roc_obj, col = "steelblue", main = paste("ROC Curve (AUC =", round(auc_val, 3), ")"))
# 
# barplot(sort(rf_model$variable.importance, decreasing = TRUE),
#         col = "tan", main = "Variable Importance (Impurity)",
#         ylab = "Mean Decrease in Node Impurity")



#---- Map of predicted grass->crop conversion risk ----
rf_pred <- terra::predict(
  pred_stack, rf_model,
  fun = function(model, data) predict(model, data)$predictions[, 2],
  na.rm = TRUE
)

plot(rf_pred)

# remove protected areas and already converted areas from the map
# force cropland to zero
rf_adj <- ifel(cropland == 1, 0, rf_pred)

# force GNP inside area to zero
rf_adj <- ifel(rasterize(vect(gnp), rf_adj, field = 1) == 1, 0, rf_adj)

## start from adjusted prediction
# cropland == 1  → 0 risk
rf_adj <- ifel(cropland == 1, 0, rf_pred)

# GNP → 0 risk
gnp_r  <- rasterize(vect(gnp), rf_adj, field = 1)
rf_adj <- ifel(gnp_r == 1, 0, rf_adj)

## build a mask: 1 = cropland or GNP, NA elsewhere
zero_mask <- ifel(cropland == 1, 1, NA)        # cropland
zero_mask <- cover(zero_mask, ifel(gnp_r == 1, 1, NA))  # add GNP

## zero layer = only cropland + GNP
zero_layer <- rf_adj
zero_layer[is.na(zero_mask)] <- NA   # keep values only where zero_mask == 1

## risk layer = everything else inside the study area
risk_layer <- rf_adj
risk_layer[!is.na(zero_mask)] <- NA  # drop cropland + GNP
risk_layer[risk_layer == 0] <- NA    # in case any other true zeros exist

## Plot: risk (yellow→red) + zero risk (grey), outside area white
# Save to disk 
png("./results/landcover/pilot/grassland_conversion_risk.png", width = 2000, height = 1600, res = 300)

cols <- rev(hcl.colors(100, "YlOrRd"))

plot(risk_layer,
     col    = cols,
     axes   = FALSE,
     legend = TRUE,
     zlim   = c(0, 1)
)

# Overlay 0-risk cropland + GNP in grey
plot(zero_layer,
     col    = "#E6E6E6",
     add    = TRUE,
     legend = FALSE
)

plot(st_geometry(gnp), add = TRUE, border = "black", lwd = 0.8)
text(gnp_cent[1], gnp_cent[2], "Grasslands National Park", cex = 0.5)

dev.off()



#---- Maps showing top % of cells at risk ----
### BAU in pilot area
area_ha <- 3000

# cell area in m2 and ha (from risk raster)
cell_area_m2 <- prod(res(rf_adj))      # 30 * 30 = 900 if 30 m
cell_area_ha <- cell_area_m2 / 10000   # 0.09 ha if 30 m

# number of pixels corresponding to 500 ha
n_target <- round(area_ha / cell_area_ha)
n_target

# number of *valid* cells in the study area (inside mask)
total_cells <- global(!is.na(rf_adj), "sum", na.rm = TRUE)[1, "sum"]

prop_target <- n_target / total_cells   # fraction of the study area
prop_target

# number of *valid* cells in the study area (inside mask)
total_cells <- global(!is.na(rf_adj), "sum", na.rm = TRUE)[1, "sum"]

prop_target <- n_target / total_cells   # fraction of the study area
prop_target

# quantile where about n_target cells are above the threshold
vals <- values(risk_layer, na.rm = TRUE)
thr  <- quantile(vals, probs = 1 - prop_target, na.rm = TRUE)
thr

# build top-risk layer
top_layer <- risk_layer
top_layer[top_layer < thr] <- NA

## Plot: risk (yellow→red) + zero risk (grey), outside area white, top % pixels cyan
# save to disk 
png("./results/landcover/pilot/grassland_conversion_risk_top6.png", width = 2000, height = 1600, res = 300)

cols <- rev(hcl.colors(100, "YlOrRd"))   # same risk palette as before

# Base: risk
plot(risk_layer,
     col    = cols,
     axes   = FALSE,
     legend = TRUE,
     zlim   = c(0, 1),
)

# 0-risk pixels (cropland + GNP) in light grey
plot(zero_layer,
     col    = "#E6E6E6",   # your lighter grey
     add    = TRUE,
     legend = FALSE
)

# Top 500 ha of risk in teal
plot(top_layer,
     col    = "cyan",
     add    = TRUE,
     legend = FALSE
)

# Park outline + label
plot(st_geometry(gnp), add = TRUE, border = "black", lwd = 0.8)
text(gnp_cent[1], gnp_cent[2], "Grasslands National Park", cex = 0.5)

dev.off()


### Rapid Conversion Scenario (80th percentile of risk)
# Extract at-risk pixel values (risk_layer excludes cropland + GNP already)
vals <- values(risk_layer, na.rm = TRUE)

# Threshold at the 80th percentile (rapid conversion scenario)
thr_rapid <- quantile(vals, probs = 0.80, na.rm = TRUE)
thr_rapid  # optional: print threshold

# Build raster of rapid-conversion pixels
rapid_layer <- risk_layer
rapid_layer[rapid_layer < thr_rapid] <- NA   # keep only >= 80th percentile

# Compute area of rapid conversion (optional but useful)
n_rapid <- global(!is.na(rapid_layer), "sum", na.rm = TRUE)[1, "sum"]

cell_area_m2 <- prod(res(risk_layer))
cell_area_ha <- cell_area_m2 / 10000
area_rapid_ha <- n_rapid * cell_area_ha
area_rapid_ha  # optional: print area in hectares

## Plot: risk (yellow→red) + zero risk (grey), outside area white, top % pixels fuscia
# Save to disk 
png("./results/landcover/pilot/grassland_conversion_risk.png", width = 2000, height = 1600, res = 300)

cols <- rev(hcl.colors(100, "YlOrRd"))   # yellow → red palette

# Base risk map
plot(risk_layer,
     col    = cols,
     axes   = FALSE,
     legend = TRUE,
     zlim   = c(0, 1)
)

# 0-risk pixels (cropland + GNP) in light grey
plot(zero_layer,
     col    = "#E6E6E6",
     add    = TRUE,
     legend = FALSE
)

# Rapid conversion scenario in a highlight colour
plot(rapid_layer,
     col    = "#E600E6",    # choose your colour (dark teal/blue), or "#1ABC9C"
     add    = TRUE,
     legend = FALSE
)

# Park outline + label
plot(st_geometry(gnp), add = TRUE, border = "black", lwd = 0.8)
text(gnp_cent[1], gnp_cent[2], "Grasslands National Park", cex = 0.5)
dev.off()


#--- Settlement regression model ----
# stack rasters (predictors)
pred_stack <- c(crop_suit_filled, settle_prop, dist_to_road)
names(pred_stack) <- c("crop_suit", "settle_prop", "dist_to_road")

# sample training data from rasters, INCLUDING coordinates
set.seed(123)
train_df <- terra::spatSample(
  pred_stack,
  size   = 5000,
  method = "random",
  na.rm  = TRUE,
  as.df  = TRUE,
  xy     = TRUE   # <-- this adds "x" and "y" columns
)

# add response variable by extracting from 'cropland' at those coords
vals <- terra::extract(settle, train_df[, c("x", "y")])
train_df$settle <- vals$LU2020

# drop any rows with NA in predictors/response
train_df <- na.omit(train_df)


# 70% training, 30% testing
train_id <- sample(1:nrow(train_df), 0.7 * nrow(train_df))
train_data <- train_df[train_id, ]
test_data  <- train_df[-train_id, ]

train_data$settle <- factor(train_data$settle, levels = c(0, 1))
test_data$settle  <- factor(test_data$settle,  levels = c(0, 1))


# Fit ranger model
rf_model <- ranger(
  settle ~ crop_suit + settle_prop + dist_to_road,
  data = train_data,
  num.trees = 500,
  mtry = 2,                 
  min.node.size = 10,
  importance = "impurity",
  probability = TRUE,
  seed = 123
)

print(rf_model)
rf_model$variable.importance

# # model performance - needs updating
# # extract probability of cropland = 1
# test_pred <- predict(rf_model, data = test_data)$predictions[, 2]
# 
# # predicted class using 0.5 threshold
# pred_class <- ifelse(test_pred > 0.5, 1, 0)
# 
# # accuracy
# acc <- mean(pred_class == as.numeric(as.character(test_data$cropland)))
# cat("Overall accuracy:", round(acc, 3), "\n")
# 
# # confusion matrix
# confusionMatrix(
#   factor(pred_class, levels = c(0, 1)),
#   test_data$cropland,
#   positive = "1"
# )
# 
# roc_obj <- roc(test_data$cropland, test_pred)
# auc_val <- auc(roc_obj)
# plot(roc_obj, col = "steelblue", main = paste("ROC Curve (AUC =", round(auc_val, 3), ")"))
# 
# barplot(sort(rf_model$variable.importance, decreasing = TRUE),
#         col = "tan", main = "Variable Importance (Impurity)",
#         ylab = "Mean Decrease in Node Impurity")



#---- Map of predicted grass->crop conversion risk ----
rf_pred <- terra::predict(
  pred_stack, rf_model,
  fun = function(model, data) predict(model, data)$predictions[, 2],
  na.rm = TRUE
)

plot(rf_pred)

# remove protected areas and already converted areas from the map
# force cropland to zero
rf_adj <- ifel(settle == 1, 0, rf_pred)

# force GNP inside area to zero
rf_adj <- ifel(rasterize(vect(gnp), rf_adj, field = 1) == 1, 0, rf_adj)

## start from adjusted prediction
# cropland == 1  → 0 risk
rf_adj <- ifel(settle == 1, 0, rf_pred)

# GNP → 0 risk
gnp_r  <- rasterize(vect(gnp), rf_adj, field = 1)
rf_adj <- ifel(gnp_r == 1, 0, rf_adj)

## build a mask: 1 = cropland or GNP, NA elsewhere
zero_mask <- ifel(settle == 1, 1, NA)        # cropland
zero_mask <- cover(zero_mask, ifel(gnp_r == 1, 1, NA))  # add GNP

## zero layer = only cropland + GNP
zero_layer <- rf_adj
zero_layer[is.na(zero_mask)] <- NA   # keep values only where zero_mask == 1

## risk layer = everything else inside the study area
risk_layer <- rf_adj
risk_layer[!is.na(zero_mask)] <- NA  # drop cropland + GNP
risk_layer[risk_layer == 0] <- NA    # in case any other true zeros exist

## Plot: risk (yellow→red) + zero risk (grey), outside area white
# Save to disk 
png("./results/landcover/pilot/grassland_conversion_risk_settle.png", width = 2000, height = 1600, res = 300)

cols <- rev(hcl.colors(100, "YlOrRd"))

plot(risk_layer,
     col    = cols,
     axes   = FALSE,
     legend = TRUE,
     zlim   = c(0, 1)
)

# Overlay 0-risk cropland + GNP in grey
plot(zero_layer,
     col    = "#E6E6E6",
     add    = TRUE,
     legend = FALSE
)

plot(st_geometry(gnp), add = TRUE, border = "black", lwd = 0.8)
text(gnp_cent[1], gnp_cent[2], "Grasslands National Park", cex = 0.5)
plot(st_geometry(pilot), add = TRUE, border = "black", lwd = 0.8)

dev.off()



#---- Maps showing top % of cells at risk ----
### BAU in pilot area
area_ha <- 3000

# cell area in m2 and ha (from risk raster)
cell_area_m2 <- prod(res(rf_adj))      # 30 * 30 = 900 if 30 m
cell_area_ha <- cell_area_m2 / 10000   # 0.09 ha if 30 m

# number of pixels corresponding to 500 ha
n_target <- round(area_ha / cell_area_ha)
n_target

# number of *valid* cells in the study area (inside mask)
total_cells <- global(!is.na(rf_adj), "sum", na.rm = TRUE)[1, "sum"]

prop_target <- n_target / total_cells   # fraction of the study area
prop_target

# number of *valid* cells in the study area (inside mask)
total_cells <- global(!is.na(rf_adj), "sum", na.rm = TRUE)[1, "sum"]

prop_target <- n_target / total_cells   # fraction of the study area
prop_target

# quantile where about n_target cells are above the threshold
vals <- values(risk_layer, na.rm = TRUE)
thr  <- quantile(vals, probs = 1 - prop_target, na.rm = TRUE)
thr

# build top-risk layer
top_layer <- risk_layer
top_layer[top_layer < thr] <- NA

## Plot: risk (yellow→red) + zero risk (grey), outside area white, top % pixels cyan
# save to disk 
png("./results/landcover/pilot/grassland_conversion_risk_top6.png", width = 2000, height = 1600, res = 300)

cols <- rev(hcl.colors(100, "YlOrRd"))   # same risk palette as before

# Base: risk
plot(risk_layer,
     col    = cols,
     axes   = FALSE,
     legend = TRUE,
     zlim   = c(0, 1),
)

# 0-risk pixels (cropland + GNP) in light grey
plot(zero_layer,
     col    = "#E6E6E6",   # your lighter grey
     add    = TRUE,
     legend = FALSE
)

# Top 500 ha of risk in teal
plot(top_layer,
     col    = "cyan",
     add    = TRUE,
     legend = FALSE
)

# Park outline + label
plot(st_geometry(gnp), add = TRUE, border = "black", lwd = 0.8)
text(gnp_cent[1], gnp_cent[2], "Grasslands National Park", cex = 0.5)

dev.off()


### Rapid Conversion Scenario (80th percentile of risk)
# Extract at-risk pixel values (risk_layer excludes cropland + GNP already)
vals <- values(risk_layer, na.rm = TRUE)

# Threshold at the 80th percentile (rapid conversion scenario)
thr_rapid <- quantile(vals, probs = 0.80, na.rm = TRUE)
thr_rapid  # optional: print threshold

# Build raster of rapid-conversion pixels
rapid_layer <- risk_layer
rapid_layer[rapid_layer < thr_rapid] <- NA   # keep only >= 80th percentile

# Compute area of rapid conversion (optional but useful)
n_rapid <- global(!is.na(rapid_layer), "sum", na.rm = TRUE)[1, "sum"]

cell_area_m2 <- prod(res(risk_layer))
cell_area_ha <- cell_area_m2 / 10000
area_rapid_ha <- n_rapid * cell_area_ha
area_rapid_ha  # optional: print area in hectares

## Plot: risk (yellow→red) + zero risk (grey), outside area white, top % pixels fuscia
# Save to disk 
png("./results/landcover/pilot/grassland_conversion_risk.png", width = 2000, height = 1600, res = 300)

cols <- rev(hcl.colors(100, "YlOrRd"))   # yellow → red palette

# Base risk map
plot(risk_layer,
     col    = cols,
     axes   = FALSE,
     legend = TRUE,
     zlim   = c(0, 1)
)

# 0-risk pixels (cropland + GNP) in light grey
plot(zero_layer,
     col    = "#E6E6E6",
     add    = TRUE,
     legend = FALSE
)

# Rapid conversion scenario in a highlight colour
plot(rapid_layer,
     col    = "#E600E6",    # choose your colour (dark teal/blue), or "#1ABC9C"
     add    = TRUE,
     legend = FALSE
)

# Park outline + label
plot(st_geometry(gnp), add = TRUE, border = "black", lwd = 0.8)
text(gnp_cent[1], gnp_cent[2], "Grasslands National Park", cex = 0.5)
dev.off()
