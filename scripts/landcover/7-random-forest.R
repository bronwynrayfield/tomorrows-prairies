# Setup
source(file.path("scripts/landcover", "0-setup.R"))
source(file.path("scripts/landcover/functions", "functions-general.R"))
source(file.path("scripts/landcover/functions", "reclass-aafc-luts.R"))

# ---- Resolution setting ----
target_res <- 30
agg_fact   <- target_res / 30


# Modelling extent
me_r <- rast(file.path(sa_filepath, me_raster_filename))
sa_r <- rast(file.path(sa_filepath, sa_raster_filename))
sa <- st_read(file.path(sa_filepath, sa_filename))

# Ecodistricts (for cropland suitability imputation)
ecoregions <- st_read(tagged_path("./processeddatasets", ecoregions_filepath,
                                    ecoregions_filename, "modelling-extent"))

ecoregion_zones_30 <- rasterize(ecoregions, me_r, field = "ECOREGION1")

# save to predictors folder
writeRaster(ecoregion_zones_30, 
            "D:/tomorrows-prairies/processeddatasets/predictors/30m/ecoregions.tif",
            overwrite = TRUE)


# Ecodistricts (for cropland suitability imputation)
ecodistricts <- st_read(tagged_path("./processeddatasets", ecodistricts_filepath,
                                    ecodistricts_filename, "modelling-extent"))

ecozones <- st_read(tagged_path("./processeddatasets", ecozones_filepath,
                                    ecozones_filename, "modelling-extent"))

# list all files in the 30m predictors folder
pred_files <- list.files(
  "D:/tomorrows-prairies/processeddatasets/predictors/30m",
  pattern = "\\.tif$",
  full.names = TRUE
)

# print to check what's there
print(pred_files)

# read them all in as a stack
pred_stack <- rast(pred_files)

# check names and properties
names(pred_stack)
nlyr(pred_stack)

names(pred_stack) <- c("cropland_prop", "dist_to_crop", "dist_to_road", "ecoregion", "crop_suit")

# verify
names(pred_stack)


# ---- Response variable ----
# 1 = annual cropland in 2020
# 0 = native or tame grassland in 2020

lulc2010 <- rast(tagged_path(file.path("./processeddatasets", aafc_filepath), "cropped", aafc2010_filename, "modelling-extent"))
lulc2015 <- rast(tagged_path(file.path("./processeddatasets", aafc_filepath), "cropped", aafc2015_filename, "modelling-extent"))
lulc2020 <- rast(tagged_path(file.path("./processeddatasets", aafc_filepath), "cropped", aafc2020_filename, "modelling-extent"))

# binary cropland layer for 2020
is_crop_2020 <- reclass_aafc(lulc2020, scheme = "binary", binary_class = "crop", check = TRUE) %>%
  crop(sa_r) %>% mask(sa_r)

# binary grassland layer for 2020 (native + tame)
# need a combined "grassland" class — this requires two separate calls and combining
is_native_2020 <- reclass_aafc(lulc2020, scheme = "binary", binary_class = "native", check = TRUE)
is_tame_2020   <- reclass_aafc(lulc2020, scheme = "binary", binary_class = "tame",   check = TRUE)

is_grass_2020 <- ifel(is_native_2020 == 1 | is_tame_2020 == 1, 1, 0) %>%
  crop(sa_r) %>% mask(sa_r)

# check
freq(is_crop_2020)
freq(is_grass_2020)


# Sampling
n_sample <- min(
  global(is_crop_2020  == 1, "sum", na.rm = TRUE)[1,1],
  global(is_grass_2020 == 1, "sum", na.rm = TRUE)[1,1],
  50000
)

cat("Sampling", n_sample, "of each class\n")

set.seed(123)

train_1 <- spatSample(ifel(is_crop_2020  == 1, 1, NA),
                      size = n_sample, method = "random",
                      na.rm = TRUE, as.df = TRUE, xy = TRUE)
train_1$response <- 1

train_0 <- spatSample(ifel(is_grass_2020 == 1, 1, NA),
                      size = n_sample, method = "random",
                      na.rm = TRUE, as.df = TRUE, xy = TRUE)
train_0$response <- 0

train_df <- bind_rows(train_1, train_0)


# extract covariate values at sample points
cov_vals <- terra::extract(pred_stack, train_df[, c("x", "y")])
train_df <- bind_cols(train_df, cov_vals %>% select(-ID))

# drop NAs
train_df <- na.omit(train_df)

# set factors
train_df$ecoregion <- as.factor(train_df$ecoregion)
train_df$response  <- as.factor(train_df$response)

# check balance and ecoregion distribution
table(train_df$response)
table(train_df$ecoregion)


# filter to ecoregions with at least 100 samples
ecoregions_keep <- names(table(train_df$ecoregion)[table(train_df$ecoregion) >= 100])

train_df <- train_df %>%
  filter(ecoregion %in% ecoregions_keep) %>%
  mutate(ecoregion = droplevels(ecoregion))

table(train_df$ecoregion)

# fit model
rf_suit <- ranger(
  response ~ crop_suit + cropland_prop + dist_to_road + ecoregion,
  data          = train_df,
  num.trees     = 500,
  mtry          = 2,
  min.node.size = 10,
  importance    = "impurity",
  probability   = TRUE,
  respect.unordered.factors = "partition",
  seed          = 123
)

print(rf_suit)
rf_suit$variable.importance


library(pROC)
library(caret)

# OOB predicted probabilities
oob_probs <- rf_suit$predictions[, 2]
actual    <- as.numeric(as.character(train_df$response))

# Brier score
mean((oob_probs - actual)^2)

# ROC and AUC
roc_obj <- roc(actual, oob_probs)
auc(roc_obj)
plot(roc_obj, col = "steelblue",
     main = paste("ROC Curve (AUC =", round(auc(roc_obj), 3), ")"))

# confusion matrix at 0.5 threshold
pred_class <- factor(ifelse(oob_probs >= 0.5, 1, 0), levels = c(0, 1))
actual_fac <- factor(actual, levels = c(0, 1))
confusionMatrix(pred_class, actual_fac, positive = "1")

# variable importance plot
barplot(sort(rf_suit$variable.importance, decreasing = TRUE),
        col  = "steelblue",
        main = "Variable importance — cropland suitability model",
        ylab = "Mean decrease in node impurity")
















crop_codes_2020  <- lut_2020$Value[lut_2020$LU2020 %in% 
                                     c("Annual Cropland", "Land Converted to Annual Cropland")]
grass_codes_2020 <- lut_2020$Value[lut_2020$LU2020 %in%
                                     c("Prairie Grassland", "Temperate Grassland",
                                       "Temperate Rainforest Grassland", "Boreal Grassland",
                                       "Montane Grassland", "Alpine Tundra Grassland",
                                       "Cropland", "Land Converted to Cropland")]

lulc2020_num <- as.numeric(lulc2020)

is_crop_2020  <- app(lulc2020_num, fun = function(x) as.integer(x %in% crop_codes_2020))
is_grass_2020 <- app(lulc2020_num, fun = function(x) as.integer(x %in% grass_codes_2020))

# mask both to sa_r
is_crop_2020  <- crop(is_crop_2020,  sa_r) %>% mask(sa_r)
is_grass_2020 <- crop(is_grass_2020, sa_r) %>% mask(sa_r)

# check counts
freq(is_crop_2020)
freq(is_grass_2020)

# ---- Sample training data ----
n_crop  <- global(is_crop_2020  == 1, "sum", na.rm = TRUE)[1,1]
n_grass <- global(is_grass_2020 == 1, "sum", na.rm = TRUE)[1,1]
n_sample <- min(n_crop, n_grass, 50000)

cat("Cropland cells:", n_crop, "\n")
cat("Grassland cells:", n_grass, "\n")
cat("Sampling:", n_sample, "of each\n")

set.seed(123)

train_1 <- spatSample(
  ifel(is_crop_2020 == 1, 1, NA),
  size   = n_sample,
  method = "random",
  na.rm  = TRUE,
  as.df  = TRUE,
  xy     = TRUE
)
train_1$response <- 1

train_0 <- spatSample(
  ifel(is_grass_2020 == 1, 1, NA),
  size   = n_sample,
  method = "random",
  na.rm  = TRUE,
  as.df  = TRUE,
  xy     = TRUE
)
train_0$response <- 0

# ---- Build training dataset ----
train_df <- bind_rows(train_1, train_0)

# stack covariates - now including distance to roads
pred_stack <- c(crop_suit_filled, cropland_prop, dist_to_road, ecoregion_zones)
names(pred_stack) <- c("crop_suit", "cropland_prop", "dist_to_road", "ecoregion")

# extract covariate values at sample points
cov_vals <- terra::extract(pred_stack, train_df[, c("x", "y")])
train_df <- bind_cols(train_df, cov_vals %>% select(-ID))

# drop NAs
train_df <- na.omit(train_df)

# set factors
train_df$ecoregion <- as.factor(train_df$ecoregion)
train_df$response  <- as.factor(train_df$response)

# check
table(train_df$response)
table(train_df$ecoregion)

# ---- Fit ranger model ----
rf_model_suit <- ranger(
  response ~ crop_suit + cropland_prop + dist_to_road + ecoregion,
  data          = train_df,
  num.trees     = 500,
  mtry          = 2,
  min.node.size = 10,
  importance    = "impurity",
  probability   = TRUE,
  respect.unordered.factors = "partition",
  seed          = 123
)

print(rf_model_suit)
rf_model_suit$variable.importance



# create ecodistrict zones at 30m
ecodistrict_zones <- rasterize(vect(ecodistricts), pred_stack[[1]], field = "ECODISTRIC")

# verify
freq(ecodistrict_zones) %>% nrow()

# extract ecodistrict ID at training points
ecodistrict_vals <- terra::extract(ecodistrict_zones, train_df[, c("x", "y")])
train_df$ecodistrict <- as.factor(ecodistrict_vals[, 2])

# check distribution
eco_dist_counts <- table(train_df$ecodistrict)
sum(eco_dist_counts >= 10)
sum(eco_dist_counts >= 20)





### to run

eco_dist_folds <- names(eco_dist_counts[eco_dist_counts >= 20])

cv_results_suit <- data.frame(
  ecodistrict = eco_dist_folds,
  auc         = NA,
  sensitivity = NA,
  specificity = NA
)

for (i in seq_along(eco_dist_folds)) {
  
  test_fold  <- train_df %>% filter(ecodistrict == eco_dist_folds[i])
  train_fold <- train_df %>% filter(ecodistrict != eco_dist_folds[i]) %>%
    mutate(ecoregion = droplevels(ecoregion))
  
  if (length(unique(test_fold$response)) < 2) next
  
  rf_cv <- ranger(
    response ~ crop_suit + cropland_prop + dist_to_road + ecoregion,
    data          = train_fold,
    num.trees     = 500,
    mtry          = 2,
    min.node.size = 10,
    probability   = TRUE,
    respect.unordered.factors = "partition",
    seed          = 123
  )
  
  preds  <- predict(rf_cv, data = test_fold)$predictions[, 2]
  actual <- as.numeric(as.character(test_fold$response))
  
  roc_cv <- roc(actual, preds, quiet = TRUE)
  cv_results_suit$auc[i]         <- auc(roc_cv)
  cv_results_suit$sensitivity[i] <- coords(roc_cv, "best")$sensitivity
  cv_results_suit$specificity[i] <- coords(roc_cv, "best")$specificity
}

summary(cv_results_suit$auc)
mean(cv_results_suit$auc, na.rm = TRUE)

# join ecoregion to cv results
cv_results_suit <- cv_results_suit %>%
  left_join(
    as.data.frame(ecodistricts) %>%
      select(ECODISTRIC, ECOREGION_) %>%
      mutate(ECODISTRIC = as.character(ECODISTRIC)),
    by = c("ecodistrict" = "ECODISTRIC")
  )

ecoregion_lookup <- as.data.frame(ecoregions) %>%
  select(ECOREGION_, ECOREGION1) %>%
  distinct()


# summarise by ecoregion
cv_results_suit %>%
  group_by(ECOREGION_) %>%
  summarise(
    n         = n(),
    mean_auc  = round(mean(auc, na.rm = TRUE), 3),
    min_auc   = round(min(auc, na.rm = TRUE), 3),
    mean_sens = round(mean(sensitivity, na.rm = TRUE), 3),
    mean_spec = round(mean(specificity, na.rm = TRUE), 3)
  ) %>%
  left_join(ecoregion_lookup, by = "ECOREGION_") %>%
  select(ECOREGION1, n, mean_auc, min_auc, mean_sens, mean_spec) %>%
  arrange(mean_auc)

cv_results_suit %>%
  group_by(ECOREGION_) %>%
  summarise(
    n          = n(),
    mean_auc   = round(mean(auc, na.rm = TRUE), 3),
    min_auc    = round(min(auc, na.rm = TRUE), 3),
    mean_sens  = round(mean(sensitivity, na.rm = TRUE), 3),
    mean_spec  = round(mean(specificity, na.rm = TRUE), 3)
  ) %>%
  left_join(ecoregion_lookup, by = "ECOREGION_") %>%
  select(ECOREGION1, n, mean_auc, min_auc, mean_sens, mean_spec) %>%
  arrange(mean_auc)


## end to run


# Prediction
# create 300m template
sa_r_300 <- aggregate(sa_r, fact = 10, fun = "modal")

# resample continuous predictors with bilinear
pred_stack_300 <- resample(pred_stack[[c("crop_suit", "cropland_prop", "dist_to_road")]], 
                           sa_r_300, method = "bilinear")

# resample categorical predictor with nearest neighbour
ecoregion_300 <- resample(pred_stack[["ecoregion"]], sa_r_300, method = "near")

# combine
pred_stack_300 <- c(pred_stack_300, ecoregion_300)

# verify
res(pred_stack_300)
names(pred_stack_300)


rf_suit_pred <- terra::predict(
  pred_stack_300, rf_suit,
  fun = function(model, data, ...) {
    library(ranger)
    predict(model, data)$predictions[, 2]
  },
  na.rm = TRUE
)

plot(rf_suit_pred, main = "Cropland suitability probability")

writeRaster(rf_suit_pred, 
            "D:/tomorrows-prairies/output/landcover/PCESC_cropland_suitability_prob_300m.tif",
            overwrite = TRUE)




# Plot
ecoregions_sf <- st_as_sf(ecoregions)
ecozone_sf    <- st_as_sf(ecozones)

# clip to study area
ecoregions_sf <- st_intersection(ecoregions_sf, sa)
ecozone_sf    <- st_intersection(ecozone_sf, sa)

# resample is_grass_2020 to 300m
is_grass_300 <- resample(is_grass_2020, sa_r_300, method = "near")

# load and resample protected areas to 300m
pa_hard_bin        <- rast(tagged_path("./processeddatasets", pa_filepath, pa_filename, "hard", "binary"))
pa_hard_bin_300 <- resample(pa_hard_bin, sa_r_300, method = "near")

# mask removes non-grassland and also protected areas
rf_suit_masked <- rf_suit_pred %>%
  mask(sa_r_300) %>%
  mask(is_grass_300, maskvalue = 0) %>%
  mask(pa_hard_bin_300, maskvalue = 1)

# quick check - should have fewer non-NA cells than before
global(!is.na(rf_suit_masked), "sum", na.rm = TRUE)


# plot
cols <- rev(hcl.colors(100, "YlOrRd"))

ggplot() +
  # grey background using study area polygon
  geom_sf(data = sa, fill = "grey85", colour = NA) +
  # prediction
  geom_spatraster(data = rf_suit_masked) +
  scale_fill_gradientn(colours = cols, na.value = NA,
                       limits = c(0, 1), name = "Probability") +
  # boundaries
  geom_sf(data = ecoregions_sf, fill = NA, colour = "black", linewidth = 0.3) +
  geom_sf(data = ecozone_sf,    fill = NA, colour = "black",  linewidth = 0.8) +
  # cities
  geom_sf(data = cities_sf, colour = "white", size = 2, shape = 21, fill = "white") +
  geom_sf_label(data = cities_sf, aes(label = city),
                size = 5, fill = "white", alpha = 0.7,
                label.size = 0.1, nudge_y = 30000) +
  theme_void() +
  theme(
    legend.title    = element_text(size = 16),
    legend.text     = element_text(size = 14),
    legend.key.size = unit(1.5, "cm"))

ggsave("PCESC_cropland_suitability_prob_masked1.png", 
       width = 20, height = 14, dpi = 300, bg = "white")



# hardcode from known results - clump 10787, candidate 4
c4_pt <- vect(data.frame(x = -103.32356, y = 51.47929),
              geom = c("x", "y"), crs = "EPSG:4326") %>%
  project(crs(pred_stack))

buf_field     <- buffer(c4_pt, width = 3000)
buf_landscape <- buffer(c4_pt, width = 15000)

# crop pred_stack to landscape buffer
pred_stack_c4 <- crop(pred_stack, buf_landscape) %>% mask(buf_landscape)

# predict at 30m for focal area
rf_suit_pred_c4 <- terra::predict(
  pred_stack_c4, rf_suit,
  fun = function(model, data, ...) {
    library(ranger)
    predict(model, data)$predictions[, 2]
  },
  na.rm = TRUE
)



# binary grassland layer for 2010
is_native_2010 <- reclass_aafc(lulc2010, scheme = "binary", binary_class = "native", check = FALSE)
is_tame_2015   <- reclass_aafc(lulc2015, scheme = "binary", binary_class = "tame",   check = FALSE)
is_grass_2010  <- ifel(is_native_2010 == 1 | is_tame_2015 == 1, 1, 0)

# crop to field buffer
is_grass_c4_field <- crop(is_grass_2010, buf_field) %>% mask(buf_field)
pa_c4_field       <- crop(pa_hard_bin,   buf_field) %>% mask(buf_field)

rf_suit_c4_field <- rf_suit_pred_c4 %>%
  crop(buf_field) %>%
  mask(buf_field) %>%
  mask(is_grass_c4_field, maskvalue = 0) %>%
  mask(pa_c4_field, maskvalue = 1)

# get field buffer bbox for plot limits
field_bbox <- st_bbox(buf_field_sf)

cols <- rev(hcl.colors(100, "YlOrRd"))

ggplot() +
  geom_sf(data = buf_field_sf, fill = "grey85", colour = NA) +
  geom_spatraster(data = rf_suit_c4_field) +
  scale_fill_gradientn(colours = cols, na.value = NA,
                       limits = c(0, 1), name = "Probability") +
  coord_sf(xlim = c(field_bbox["xmin"], field_bbox["xmax"]),
           ylim = c(field_bbox["ymin"], field_bbox["ymax"])) +
  theme_void() +
  theme(legend.position = "none")

ggsave("PCESC_cropland_suitability_c4_inset.png",
       width = 5, height = 5, dpi = 300, bg = "white")






# binary grassland from 2010 (native) and 2015 (tame)
# resample to 300m
lulc2015_300 <- resample(lulc2015, sa_r_300, method = "near")
lulc2010_300 <- resample(lulc2010, sa_r_300, method = "near")

# native grass from 2010
is_native_300 <- reclass_aafc(lulc2010_300, scheme = "binary", 
                              binary_class = "native", check = FALSE)

# tame grass from 2015
is_tame_300   <- reclass_aafc(lulc2015_300, scheme = "binary", 
                              binary_class = "tame", check = FALSE)

# combined grassland mask
is_grass_300 <- ifel(is_native_300 == 1 | is_tame_300 == 1, 1, 0) %>%
  crop(sa_r_300) %>%
  mask(sa_r_300)

# check
freq(is_grass_300)

# get top 5% threshold per ecoregion
# first mask prediction to study area and grassland
rf_suit_study <- rf_suit_pred %>%
  mask(sa_r_300) %>%
  mask(is_grass_300, maskvalue = 0) %>%
  mask(pa_hard_bin_300, maskvalue = 1)

# extract values with ecoregion zones
ecoregion_300 <- resample(pred_stack[["ecoregion"]], sa_r_300, method = "near")

# compute 95th percentile threshold per ecoregion
eco_vals <- data.frame(
  pred     = values(rf_suit_study,  na.rm = FALSE),
  ecoregion = values(ecoregion_300, na.rm = FALSE)
) %>%
  na.omit()

thresholds <- eco_vals %>%
  group_by(ecoregion) %>%
  summarise(threshold = quantile(lyr1, 0.95, na.rm = TRUE))

print(thresholds)

# build threshold raster
threshold_r <- classify(ecoregion_300,
                        as.matrix(thresholds[, c("ecoregion", "threshold")]))

# top 5% per ecoregion
top5_risk <- ifel(rf_suit_study >= threshold_r, 1, NA)

# check cell count
global(!is.na(top5_risk), "sum", na.rm = TRUE)

# plot
ggplot() +
  geom_sf(data = sa, fill = "grey85", colour = NA) +
  geom_spatraster(data = top5_risk) +
  scale_fill_gradientn(colours = "darkred", na.value = NA,
                       name = "") +
  geom_sf(data = ecoregions_sf, fill = NA, colour = "grey40", linewidth = 0.3) +
  geom_sf(data = ecozone_sf,    fill = NA, colour = "black",  linewidth = 0.8) +
  theme_void() +
  theme(legend.position = "none")





# global top 5% threshold
global_threshold <- quantile(values(rf_suit_masked, na.rm = TRUE), 0.975, na.rm = TRUE)
cat("Global 95th percentile threshold:", global_threshold, "\n")

# binary layer
top5_risk_global <- ifel(rf_suit_study >= global_threshold, 1, NA)

# check cell count - should be ~5% of non-NA cells
global(!is.na(top5_risk_global), "sum", na.rm = TRUE)
global(!is.na(rf_suit_study),    "sum", na.rm = TRUE)


library(ggnewscale)

top5_df <- as.data.frame(top5_risk_global, xy = TRUE) %>%
  filter(!is.na(lyr1))



# create a raster with value 1 for top 5%, NA elsewhere — already have this
# just plot it with a fixed colour using after_stat

top5_risk_factor <- as.factor(top5_risk_global)



ggplot() +
  geom_sf(data = sa, fill = "grey55", colour = NA) +
  geom_spatraster(data = rf_suit_masked, alpha=0.6) +
  scale_fill_gradientn(colours = cols, na.value = NA,
                       limits = c(0, 1), name = "Probability") +
  ggnewscale::new_scale_fill() +
  geom_spatraster(data = top5_risk_factor) +
  scale_fill_manual(values = c("1" = "cyan"), na.value = NA, guide = "none") +
  geom_sf(data = ecoregions_sf, fill = NA, colour = "black", linewidth = 0.3) +
  geom_sf(data = ecozone_sf,    fill = NA, colour = "black",  linewidth = 0.8) +
  geom_sf(data = cities_sf, colour = "white", size = 2, shape = 21, fill = "white") +
  geom_sf_label(data = cities_sf, aes(label = city),
                size = 2.5, fill = "white", alpha = 0.7,
                label.size = 0.1, nudge_y = 30000) +
  theme_void() +
  theme(
    legend.title    = element_text(size = 20),
    legend.text     = element_text(size = 18),
    legend.key.size = unit(2, "cm")
  )

ggsave("cropland_suitability_top5_overlay6.png",
       width = 20, height = 14, dpi = 300, bg = "white")


ggplot() +
  geom_sf(data = sa, fill = "grey55", colour = NA) +
  geom_spatraster(data = rf_suit_masked, alpha = 0.6) +
  scale_fill_gradientn(colours = cols, na.value = NA,
                       limits = c(0, 1), name = "Probability") +
  ggnewscale::new_scale_fill() +
  geom_spatraster(data = top5_risk_factor) +
  scale_fill_manual(values = c("1" = "cyan"), na.value = NA, guide = "none") +
  geom_sf(data = ecoregions_sf, fill = NA, colour = "black", linewidth = 0.3) +
  geom_sf(data = ecozone_sf,    fill = NA, colour = "black", linewidth = 0.8) +
  geom_sf(data = cities_sf, colour = "white", size = 5, shape = 21, fill = "white") +
  geom_sf_label(data = cities_sf, aes(label = city),
                size = 6, fill = "white", alpha = 0.7,
                label.size = 0.1, nudge_y = 30000) +
  theme_void() +
  theme(
    legend.title    = element_text(size = 20),
    legend.text     = element_text(size = 18),
    legend.key.size = unit(2, "cm")
  )

ggsave("cropland_suitability_top5_overlay6.png",
       width = 20, height = 14, dpi = 300, bg = "white")

















# build binary top 5% raster
# reclassify ecoregion zones to thresholds
threshold_r <- classify(ecoregion_300, 
                        as.matrix(thresholds[, c("ecoregion", "threshold")]))

# cells above threshold = 1, else 0
top5_risk <- ifel(rf_suit_study >= threshold_r, 1, NA)

# check
global(!is.na(top5_risk), "sum", na.rm = TRUE)

# plot
plot(top5_risk, col = "darkred", 
     main = "Top 5% at-risk grassland cells by ecoregion",
     legend = FALSE, axes = FALSE)
plot(st_geometry(ecoregions_sf), add = TRUE, border = "grey40", lwd = 0.3)
plot(st_geometry(ecozone_sf),    add = TRUE, border = "black",  lwd = 0.8)