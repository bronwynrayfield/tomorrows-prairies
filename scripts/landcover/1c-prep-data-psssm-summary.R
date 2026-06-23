# =============================================================================
# 1-prep-data-psssm-summary.R
# =============================================================================
#
# Purpose:
#   Downloads AAFC Percent Saturated Surface Soil Moisture (PSSSM) anomaly
#   GeoTIFFs and computes two covariate layers for use in cropland modelling:
#   the interannual mean and standard deviation of growing-season soil moisture
#   anomaly over 2011–2020.
#
# Inputs:
#   - AAFC PSSSM monthly anomaly GeoTIFFs downloaded from Canada's Open Data
#     portal (0.25 degree resolution, WGS84):
#       https://agriculture.canada.ca/atlas/data_donnees/psssm/
#       data_donnees/tif/anomaly/{year}/monthly/
#       SMDiffAvg_SMUDP2_{year}_Month_{month}_WGS84_P8.tif
#     Years: 2011–2020, Months: April (4) – August (8)
#
# Outputs:
#   - psssm_mean_2011_2020.tif : pixel-wise mean of annual growing-season
#       PSSM anomaly across 2011–2020
#   - psssm_sd_2011_2020.tif   : pixel-wise SD of annual growing-season
#       PSSM anomaly across 2011–2020 (interannual variability)
#   - Downloaded raw GeoTIFFs in soil/pssm/
#
# Notes:
#   - Annual growing-season value is the mean across April–August for each
#     year; SD is then computed across the 10 annual means (purely interannual)
#
# Dependencies:
#   - 0-setup.R (parameters, file paths, functions)
# =============================================================================

source(file.path("scripts/landcover", "0-setup.R"))

# --- Configuration ----
base_url  <- "https://agriculture.canada.ca/atlas/data_donnees/psssm/data_donnees/tif/anomaly"
years     <- 2011:2020
months    <- 4:8          # April to August
raw_dir   <- file.path("./existingdatasets", psssm_filepath)
                     
# Load modelling extent
me <- st_read(file.path(sa_filepath, me_filename))
# Crop extent: me with 1-degree buffer, in WGS84 to match downloaded rasters
# Get CRS from a single file
ref_file <- file.path("./existingdatasets", psssm_filepath,
                      sprintf("SMDiffAvg_SMUDP2_%d_Month_%d_WGS84_P8.tif", years[1], months[1]))
download.file(sprintf("%s/%d/monthly/%s", base_url, years[1], basename(ref_file)),
              ref_file, mode = "wb", quiet = TRUE)
pssm_crs <- crs(rast(ref_file))
me_ext  <- ext(vect(st_buffer(st_transform(me, pssm_crs), 1)))

# --- Download GeoTIFFs ----
for (yr in years) {
  for (mo in months) {
    fname <- sprintf("SMDiffAvg_SMUDP2_%d_Month_%d_WGS84_P8.tif", yr, mo)
    url   <- sprintf("%s/%d/monthly/%s", base_url, yr, fname)
    dest <- file.path(raw_dir, fname)
    download.file(url, dest, mode = "wb", quiet = TRUE)
    # Crop to buffered AOI immediately to reduce file size for all subsequent steps
    r <- crop(rast(dest), me_ext)
    writeRaster(r, dest, overwrite = TRUE)
  }
}

# --- Per-year mean across April–August ----
# Returns a single-layer SpatRaster: mean psssm anomaly for that year
annual_mean <- function(yr) {
  files <- sprintf(
    file.path(raw_dir, "SMDiffAvg_SMUDP2_%d_Month_%d_WGS84_P8.tif"),
    yr, months
  )
  stack <- rast(files)
  mean(stack, na.rm = TRUE)  # pixel-wise mean across months
}

annual_rasters <- lapply(years, annual_mean)

# --- Mean and SD across years ----
year_stack <- rast(annual_rasters)  # stack of 11 annual means
psssm_mean <- mean(year_stack, na.rm = TRUE)
psssm_sd   <- stdev(year_stack, na.rm = TRUE)

# --- Write output layers ----
writeRaster(psssm_mean, file.path(raw_dir, psssm_mean_filename), overwrite = TRUE)
writeRaster(psssm_sd, file.path(raw_dir, psssm_sd_filename), overwrite = TRUE)
