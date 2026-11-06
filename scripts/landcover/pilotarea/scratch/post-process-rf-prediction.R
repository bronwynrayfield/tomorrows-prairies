library(terra)

# Mask out existing cropland
rf_pred_non_crop <- mask(rf_pred, cropland_focal, maskvalues = 1, updatevalue = NA)

# make sure it's a SpatRaster (not a dataset)
rf_pred_non_crop <- rast(rf_pred_non_crop)

# 99th percentile across all non-NA cells
threshold_99 <- quantile(values(rf_pred_non_crop), 0.99, na.rm = TRUE)

# top 1% mask
rf_top1 <- rf_pred_non_crop > threshold_99

# plot with cyan for top 1% pixels
plot(rf_pred_non_crop, col = hcl.colors(20, "YlOrRd", rev = TRUE),
     legend.args = list(text = "Probability", side = 4, line = 2.5, cex = 0.8))
plot(rf_top1, add = TRUE, col = c(NA, "cyan"), legend = FALSE)



