library(terra)
library(dplyr)
library(ranger)

# stack your rasters
pred_stack <- c(crop_suit_filled, cropland_prop)
names(pred_stack) <- c("crop_suit", "cropland_prop")

# sample training data from rasters
set.seed(123)
train_df <- as.data.frame(terra::spatSample(pred_stack, size = 5000, method = "random", 
                                            na.rm = TRUE, as.df = TRUE))
# add response variable
train_df$cropland <- terra::extract(cropland_focal, 
                                    terra::vect(data.frame(x=train_df$x, y=train_df$y)))$cropland
train_df <- na.omit(train_df)


# stack all three rasters together
stack_all <- c(crop_suit_filled, cropland_prop, dist_to_road, cropland)
names(stack_all) <- c("crop_suit", "cropland_prop", "dist_to_road", "cropland")

# sample training data from rasters
set.seed(123)
train_df <- as.data.frame(terra::spatSample(stack_all, size = 5000, method = "random", 
                                            na.rm = TRUE, as.df = TRUE))

# 70% training, 30% testing
train_id <- sample(1:nrow(train_df), 0.7 * nrow(train_df))
train_data <- train_df[train_id, ]
test_data  <- train_df[-train_id, ]

train_data$cropland <- factor(train_data$cropland, levels = c(0, 1))
test_data$cropland  <- factor(test_data$cropland,  levels = c(0, 1))


# Fit the ranger model
rf_model <- ranger(
  cropland ~ crop_suit + cropland_prop + dist_to_road,
  data = train_data,
  num.trees = 500,
  mtry = 2,                  # 2 predictors
  min.node.size = 10,
  importance = "impurity",
  probability = TRUE,
  seed = 123
)

print(rf_model)
rf_model$variable.importance

# 1. Extract probability of cropland = 1
test_pred <- predict(rf_model, data = test_data)$predictions[, 2]

# 2. Predicted class using 0.5 threshold
pred_class <- ifelse(test_pred > 0.5, 1, 0)

# 3. Accuracy
acc <- mean(pred_class == as.numeric(as.character(test_data$cropland)))
cat("Overall accuracy:", round(acc, 3), "\n")

# 4. Confusion matrix
library(caret)
confusionMatrix(
  factor(pred_class, levels = c(0, 1)),
  test_data$cropland,
  positive = "1"
)

library(pROC)
roc_obj <- roc(test_data$cropland, test_pred)
auc_val <- auc(roc_obj)
plot(roc_obj, col = "steelblue", main = paste("ROC Curve (AUC =", round(auc_val, 3), ")"))

barplot(sort(rf_model$variable.importance, decreasing = TRUE),
        col = "tan", main = "Variable Importance (Impurity)",
        ylab = "Mean Decrease in Node Impurity")


predictors <- stack_all[[c("crop_suit", "cropland_prop", "dist_to_road")]]


rf_pred <- terra::predict(
  predictors, rf_model,
  fun = function(model, data) predict(model, data)$predictions[, 2],
  na.rm = TRUE
)

plot(rf_pred)#, main = "Predicted Probability of Cropland")

zoom(rf_pred)

# Choose an area of interest (extent in map units)
# You can also define this interactively:
zoom_area <- drawExtent()  # click+drag on plot to define

# Crop both rasters to that extent
rf_zoom  <- crop(rf_pred, zoom_area)
prop_zoom <- crop(cropland_prop, zoom_area)

# Compare visually
par(mfrow = c(1, 2))
plot(prop_zoom, main = "Cropland Proportion (local window)", col = terrain.colors(50))
plot(rf_zoom, main = "Predicted Probability of Cropland", col = viridis::viridis(50))

# Plot your raster
plot(cropland_prop)

# Interactively draw a rectangle (click + drag on the map)
zoom_area <- ext(draw(cropland_prop))


