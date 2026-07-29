# collinearity-check.R
# Pairwise correlation + VIF screen to flag collinear predictors before the RF.
# Reusable: point `files` at any set of ALIGNED rasters (same grid, extent, CRS).
# RF tolerates collinearity for prediction, but it distorts variable importance
# and effect plots, so redundant predictors are worth screening out.
library(terra)
library(usdm)
library(ggplot2)
library(dplyr)
library(sf)
source(file.path("scripts/landcover", "0-setup.R"))
sa <- st_read(file.path(sa_filepath, sa_filename))

# ---- inputs ----
# in_dir   <- "D:/tomorrows-prairies/existingdatasets/climate/CanDCS-M6/geotiff"
# pattern  <- "_ssp245_1991-2020"   # the climate covariates used for training
# out_dir  <- "D:/tomorrows-prairies/output/landcover/collinearity/climate"
in_dir   <- file.path("processeddatasets", dss_filepath)
pattern  <- "dss_"
out_dir  <- "output/landcover/collinearity/dss"
sample_n   <- 10000      # random non-NA cells to sample (use all if fewer exist)
cor_method <- "spearman" # "spearman" (monotonic) or "pearson" (linear)
cor_th     <- 0.7        # |r| at/above this flags a pair as collinear
vif_th     <- 5         # VIF at/above this flags multicollinearity (5 = stricter)

# ---- load raw rasters, crop + mask to study area ----
files <- list.files(in_dir, pattern = pattern, full.names = TRUE)
r     <- rast(files)
#names(r) <- sub("_ssp245.*$", "", basename(files))   # short variable names; edit per naming
names(r) <- sub("dss_", "", sub("\\.tif$", "", basename(files)))
sa_v  <- vect(st_transform(sa, crs(r)))
r     <- crop(r, sa_v) |> mask(sa_v)

# ---- sample values ----
vals <- na.omit(spatSample(r, size = sample_n, method = "random", na.rm = TRUE, as.df = TRUE))

# ---- pairwise correlation ----
cmat <- cor(vals, method = cor_method, use = "pairwise.complete.obs")
write.csv(round(cmat, 3), file.path(out_dir, "correlation_matrix.csv"))
pairs <- as.data.frame(as.table(cmat)) |>
  filter(Var1 != Var2, abs(Freq) >= cor_th) |>
  mutate(pair = paste(pmin(as.character(Var1), as.character(Var2)),
                      pmax(as.character(Var1), as.character(Var2)), sep = " - ")) |>
  distinct(pair, .keep_all = TRUE) |>
  transmute(pair, r = round(Freq, 3)) |>
  arrange(desc(abs(r)))

# ---- VIF screen ----
vif_all  <- vif(vals)
vif_step <- vifstep(vals, th = vif_th)

# ---- correlation heatmap ----
hm <- as.data.frame(as.table(cmat)) |>
  ggplot(aes(Var1, Var2, fill = Freq)) +
  geom_tile() +
  geom_text(aes(label = round(Freq, 2)), size = 3) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-1, 1), name = cor_method) +
  labs(x = NULL, y = NULL, title = "Predictor correlation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(out_dir, "correlation_heatmap.png"), hm, width = 7, height = 6, dpi = 150)

# ---- report ----
cat("Correlation method:", cor_method, "\n\n")
cat("Pairs with |r| >=", cor_th, ":\n"); print(pairs, row.names = FALSE)
cat("\nVIF (all variables):\n"); print(vif_all)
cat("\nvifstep (th =", vif_th, ") -- variables flagged for exclusion:\n")
print(vif_step@excluded)
cat("\nRetained set:\n"); print(vif_step@results$Variables)
vif(vals[, c("gddgrow_5", "prtot", "hddheat_18")])
vif(vals[, c("gddgrow_5", "prtot", "tx_mean")])