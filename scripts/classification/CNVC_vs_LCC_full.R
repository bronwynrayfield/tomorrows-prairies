
###############################################
## Spatial Comparison: CNVC vs LCC
## High-resolution PNG output version
## Author: Alison Long
## Date: March 2026
###############################################

library(terra)
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(scales)
library(readr)
library(stringr)
library(purrr)
library(gt)
library(glue)

############################################################
## USER PARAMETERS
############################################################

params <- list(
  raster_a = "R:/CanadianGrasslands/WG_Classification/processed/natureserve_cnvc/masked/cnvc_algn_native_masked_pe.tif",
  raster_b = "R:/CanadianGrasslands/WG_Classification/processed/nrcan_lcc/masked/lcc_reproj_algn_binary_noshrub_masked_pe.tif",
  raster_a_label = "CNVC",
  raster_b_label = "LCC",
  ecoreg = "R:/CanadianGrasslands/processeddatasets/boundaries/national-ecological-framework-for-canada/Ecoregions/ecoregions_clean_prairies.shp",
  class_key_csv = NULL,
  equal_area = FALSE,
  target_crs = "ESRI:102008",
  area_units = "ha",
  plot_fact = 4,
  top_changes = 12,
  author = "Alison Long",
  
  ## Output directory
  outdir = "R:/CanadianGrasslands/WG_Classification/results/CNVC_vs_LCC_for_DFL"
)

if (!dir.exists(params$outdir)) dir.create(params$outdir, recursive = TRUE)


############################################################
## Helper: Save ggplot to PNG (high resolution)
############################################################
save_png <- function(plot, filename, width=8, height=6, dpi=300) {
  ggsave(
    filename = file.path(params$outdir, filename),
    plot = plot,
    width = width, height = height, dpi = dpi
  )
}

############################################################
## Helper: Save a data frame as PNG using gt
############################################################
save_table_png <- function(df, filename, title = NULL) {
  tbl <- df |>
    gt() |>
    tab_header(title = title)
  
  gtsave(tbl, filename = file.path(params$outdir, filename), vwidth = 1200)
}

############################################################
## 1. LOAD & PREPARE RASTERS
############################################################

terraOptions(
  memfrac = 0.5,    # use only 50% of available RAM
  tempdir = tempdir()   # force disk-based temp rasters
)

A0 <- rast(params$raster_a)
B0 <- rast(params$raster_b)
names(A0) <- "A"
names(B0) <- "B"
plot(A0)
plot(B0)

## Reprojection / alignment logic
if (isTRUE(params$equal_area)) {
  target <- params$target_crs
  
  A1 <- if (!identical(crs(A0), target)) project(A0, target, method = "near") else A0
  B1 <- if (!identical(crs(B0), target)) project(B0, target, method = "near") else B0
  
  if (!compareGeom(A1, B1, stopOnError = FALSE)) {
    B1 <- resample(B1, A1, method = "near")
  }
  
} else {
  if (!compareGeom(A0, B0, stopOnError = FALSE)) {
    B1 <- project(B0, crs(A0), method = "near")
    B1 <- resample(B1, A0, method = "near")
  } else {
    B1 <- B0
  }
  A1 <- A0
}

## Overlapping valid pixels
both_valid <- !is.na(A1) & !is.na(B1)
if (global(both_valid, "sum", na.rm = TRUE)[[1]] == 0)
  stop("No overlapping valid cells.")

A <- mask(A1, both_valid)
B <- mask(B1, both_valid)

rm(A0, B0, A1, B1, both_valid)
gc()


############################################################
## 2. BINARY AGREEMENT METRICS
############################################################

rcl <- matrix(c(-Inf, 0.5, 0,
                0.5, Inf, 1),
              ncol = 3, byrow = TRUE)

A_bin <- classify(A, rcl)
B_bin <- classify(B, rcl)

valid <- !is.na(A_bin) | !is.na(B_bin)
A_mask <- mask(A_bin, valid, maskvalues = FALSE)
B_mask <- mask(B_bin, valid, maskvalues = FALSE)

freq(A_mask)
freq(B_mask)
freq(c(A_mask, B_mask))

freq(A_bin)
freq(B_bin)

freq(A)
freq(B)

xt <- crosstab(c(A_mask, B_mask), useNA = FALSE) |> as.data.frame()

xt$A <- as.numeric(as.character(xt$A))
xt$B <- as.numeric(as.character(xt$B))
colnames(xt) <- c("A", "B", "Freq")

## 2×2 confusion matrix
full <- expand.grid(A = c(0,1), B = c(0,1))
xt_full <- full |>
  left_join(xt, by = c("A","B")) |>
  mutate(Freq = ifelse(is.na(Freq), 0, Freq))

cm <- xt_full |>
  pivot_wider(names_from = B, values_from = Freq) |>
  arrange(A)

cm_mat <- as.matrix(cm[, c("0","1")])
rownames(cm_mat) <- c("0","1")

TN <- cm_mat[1,1]
FP <- cm_mat[1,2]
FN <- cm_mat[2,1]
TP <- cm_mat[2,2]
N  <- TN + FP + FN + TP

po <- (TP + TN) / N
p_yes_A <- (TP + FN) / N
p_yes_B <- (TP + FP) / N
p_no_A  <- (TN + FP) / N
p_no_B  <- (TN + FN) / N
pe <- p_yes_A*p_yes_B + p_no_A*p_no_B
kappa <- ifelse(abs(1 - pe) < 1e-12, NA, (po - pe) / (1 - pe))

precision <- TP / (TP + FP + 1e-12)
recall    <- TP / (TP + FN + 1e-12)
f1        <- 2 * precision * recall / (precision + recall + 1e-12)
jacc_binary <- TP / (TP + FP + FN + 1e-12)

metrics_tbl <- tibble(
  Metric = c(
    "TP","FP","FN","TN",
    "Overall Agreement","Kappa",
    "Precision","Recall","F1","Jaccard"
  ),
  Value = c(
    round(TP,0), 
    round(FP,0), 
    round(FN,0), 
    round(TN,0),
    round(po,2),
    round(kappa,2),
    round(precision,2),
    round(recall,2),
    round(f1,2),
    round(jacc_binary,2)
  )
) |>
  mutate(Value = format(Value, scientific = FALSE, trim = TRUE, big.mark = ","))

save_table_png(metrics_tbl, "01_binary_metrics_v7.png",
               "Binary Agreement Metrics")

#gc()


############################################################
## 3. AGREEMENT & CHANGE MAPS (PNG)
############################################################

agree <- A == B

change_num <- ifel(A == B, 0,
                   ifel(A == 1 & B == 0, 1, 2))

change_r <- as.factor(change_num)
levels(change_r) <- data.frame(
  ID = c(0,1,2),
  label = c("Agree",
            "CNVC=grass, LCC=not grass",
            "CNVC=not grass, LCC=grass")
)

## Sampling for plotting
#fact <- max(1, as.integer(params$plot_fact))
#pts <- spatSample(A, size = ncell(A)/fact^2,
#                  method = "regular", as.points = TRUE)
#coords <- crds(pts)

#A_vals  <- terra::extract(A, pts)[,2]
#B_vals  <- terra::extract(B, pts)[,2]
#AG_vals <- terra::extract(agree, pts)[,2]
#CH_vals <- terra::extract(change_r, pts)[,2]


# ---- Sampling factor ----
fact <- max(1, as.integer(params$plot_fact))

# ---- Regular spatial sampling for large rasters ----
# This picks one point per fact x fact block WITHOUT loading the full raster into memory.
target_n <- ceiling(ncell(A) / fact^2)

pts <- spatSample(
  A, 
  size   = target_n,
  method = "regular",
  as.points = TRUE
)

# Coordinates of sampled points
coords <- crds(pts)

# Convert XY coords to cell indices
cell_ids <- cellFromXY(A, coords)

# Extract values directly (FAST, LOW MEMORY)
A_vals  <- A[cell_ids]
B_vals  <- B[cell_ids]
AG_vals <- agree[cell_ids]
CH_vals <- change_r[cell_ids]

df_ag <- tibble(x=coords[,1], y=coords[,2],
                agree = factor(ifelse(AG_vals==1,"Agree","Disagree")))

df_ch <- tibble(x=coords[,1], y=coords[,2],
                change = factor(CH_vals,
                                levels=c("Agree",
                                         "CNVC=grass, LCC=not grass",
                                         "CNVC=not grass, LCC=grass")))

### Agreement map
p_ag <- ggplot(df_ag, aes(x, y, fill=agree)) +
  geom_raster() +
  scale_fill_manual(values=c("Agree"="#2ca25f","Disagree"="#de2d26")) +
  coord_equal() + theme_minimal() +
  labs(title="Agreement Map",
       subtitle=paste(params$raster_a_label,"vs",params$raster_b_label),
       x=NULL,y=NULL)

save_png(p_ag, "02_agreement_map.png")

### Change type map
p_ch <- ggplot(df_ch, aes(x, y, fill=change)) +
  geom_raster() +
  scale_fill_manual(values=c(
    "Agree"="#7f7f7f",
    "CNVC=grass, LCC=not grass"="#d95f02",
    "CNVC=not grass, LCC=grass"="#0096C7"
  )) +
  coord_equal() + theme_minimal() +
  labs(title="Map of Disagreement Type",
       subtitle=paste(params$raster_a_label,"vs",params$raster_b_label),
       fill="Type")

save_png(p_ch, "03_change_type_map.png")

gc()


############################################################
## 4. AREA SUMMARIES (PNG TABLES)
############################################################

cell_area_m2 <- global(cellSize(A), "mean", na.rm=TRUE)[[1]]
area_factor <- c(m=1, km=1e-6, ha=1e-4)[[params$area_units]]

area_by_class <- function(r) {
  f <- freq(r)
  tibble(class=f$value,
         area=f$count * cell_area_m2 * area_factor)
}

areaA <- area_by_class(A) |>
  mutate(label=ifelse(class==1,"Grass","Not grass"),
         map=params$raster_a_label)

areaB <- area_by_class(B) |>
  mutate(label=ifelse(class==1,"Grass","Not grass"),
         map=params$raster_b_label)

area_tbl <- bind_rows(areaA, areaB) |> arrange(label,map)

save_table_png(area_tbl, "04_area_by_class.png",
               glue("Area by Class ({params$area_units})"))

## Agreement area
f_change <- freq(change_r)
agreement_area_tbl <- tibble(
  category=f_change$value,
  area=f_change$count * cell_area_m2 * area_factor
)

save_table_png(agreement_area_tbl, "05_area_by_agreement_category.png",
               glue("Area by Agreement Category ({params$area_units})"))

gc()


############################################################
## 5. AGREEMENT BY ECOREGION
############################################################

ecoreg <- vect(params$ecoreg)
if (!identical(crs(ecoreg), crs(A))) {
  ecoreg <- project(ecoreg, crs(A))
}

stopifnot("ECOREGION_" %in% names(ecoreg))

zone_raster <- rasterize(ecoreg, A, field="ECOREGION_")

names(agree) <- "agree"
names(agree)

eco_agree <- zonal(agree, zone_raster, fun="sum", na.rm=TRUE)
eco_total <- zonal(!is.na(agree), zone_raster, fun="sum", na.rm=TRUE)

names(eco_agree)
names(eco_total)

eco_agree <- eco_agree |> rename(ECOREGION_ = ECOREGION_, agree_cells = agree)
eco_total <- eco_total |> rename(ECOREGION_ = ECOREGION_, total_cells = agree)

eco_df <- eco_total |>
  left_join(eco_agree, by="ECOREGION_") |>
  mutate(
    agree_cells = replace_na(agree_cells,0),
    disagree_cells = total_cells - agree_cells
  ) |>
  left_join(
    ecoreg |>
      as.data.frame() |>
      select(ECOREGION_, ECOREGION1),
    by="ECOREGION_"
  )

eco_summary <- eco_df |>
  group_by(ECOREGION_) |>
  summarise(
    ECOREGION1 = first(ECOREGION1),
    total_cells = sum(total_cells),
    agree_cells = sum(agree_cells),
    disagree_cells = sum(disagree_cells)
  ) |>
  mutate(
    cell_area_m2 = 900,
    total_ha = total_cells * 900 * 1e-4,
    agree_ha = agree_cells * 900 * 1e-4,
    disagree_ha = disagree_cells * 900 * 1e-4,
    agree_pct = percent(agree_cells / total_cells, 0.1)
  )

save_table_png(eco_summary, "06_ecoregion_agreement.png",
               "Agreement by Ecoregion")

gc()


############################################################
## DONE
############################################################

message("All outputs saved to: ", params$outdir)