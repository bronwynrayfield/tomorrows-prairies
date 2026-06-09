# =============================================================================
# 1-4-prep-data-tsm-summary.R
# =============================================================================
#
# Purpose:
#   Downloads AAFC Total Soil Moisture (TSM) GeoTIFFs and computes two
#   covariate layers for use in cropland modelling: the 2020 growing-season
#   and planting-season mean root-zone soil moisture.
#
# Inputs:
#   - AAFC TSM weekly GeoTIFFs downloaded from Canada's Open Data portal
#     (interpolated map, native VSMB resolution):
#       https://agriculture.canada.ca/atlas/data_donnees/agClimate/
#       data_donnees/tif/vsmb/tsm/im/
#       tsm_im_{YYYYMMDD}.tif
#     Year: 2020
#     Planting season : weeks with dates in April–May   (months 4–5)
#     Growing season  : weeks with dates in June–August (months 6–8)
#
# Outputs:
#   - tsm_planting_2020.tif : pixel-wise mean of 2020 planting-season TSM
#   - tsm_growing_2020.tif  : pixel-wise mean of 2020 growing-season TSM
#   - Downloaded raw GeoTIFFs saved in existingdatasets/soil/tsm/
#
# Notes:
#   - The VSMB archive only extends to 2019 onward so a climatological
#     summary was not possible; 2020 is used as a snapshot temporally
#     matched to the 2020 Annual Crop Inventory response variable.
#   - TSM files are distributed weekly on irregular dates; the directory
#     listing is parsed to identify all available files, which are then
#     filtered to the target year and months.
#   - Files are cropped to the buffered modelling extent immediately after
#     download to reduce storage and processing overhead.
#
# Dependencies:
#   - 0-setup.R (parameters, file paths, functions)
# =============================================================================

source(file.path("scripts/landcover", "0-setup.R"))

# --- Configuration ----
base_url        <- "https://agriculture.canada.ca/atlas/data_donnees/agClimate/data_donnees/tif/vsmb/tsm/im/"
year <- 2020
planting_months <- 4:5   # April–May
growing_months  <- 6:8   # June–August
raw_dir         <- file.path("./existingdatasets", tsm_filepath)

dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

# --- Load modelling extent ----
me <- st_read(file.path(sa_filepath, me_filename))

# --- List available files and filter to 2020 and target months ----
listing_html <- readLines(base_url, warn = FALSE)

all_fnames <- unique(unlist(regmatches(
  listing_html,
  gregexpr("tsm_im_[0-9]{8}\\.tif", listing_html)
)))

file_dates    <- as.Date(regmatches(all_fnames, regexpr("[0-9]{8}", all_fnames)), format = "%Y%m%d")
target_fnames <- all_fnames[
  as.integer(format(file_dates, "%Y")) == year &
    as.integer(format(file_dates, "%m")) %in% c(planting_months, growing_months)
]

# --- Download and crop GeoTIFFs ----
for (fname in target_fnames) {
  dest <- file.path(raw_dir, fname)
  if (!file.exists(dest))
    download.file(sprintf("%s/%s", base_url, fname), dest, mode = "wb", quiet = TRUE)
}

ref_file   <- file.path(raw_dir, target_fnames[1])
me_ext     <- ext(vect(st_buffer(st_transform(me, crs(rast(ref_file))), 1)))
downloaded <- file.path(raw_dir, target_fnames)

for (f in downloaded) {
  writeRaster(crop(rast(f), me_ext), f, overwrite = TRUE)
}

# --- Seasonal means ----
dates           <- as.Date(regmatches(basename(downloaded), regexpr("[0-9]{8}", basename(downloaded))), format = "%Y%m%d")
planting_files  <- downloaded[as.integer(format(dates, "%m")) %in% planting_months]
growing_files   <- downloaded[as.integer(format(dates, "%m")) %in% growing_months]

tsm_planting <- mean(rast(planting_files), na.rm = TRUE)
tsm_growing  <- mean(rast(growing_files),  na.rm = TRUE)

# --- Write output layers ----
writeRaster(tsm_planting, file.path("./existingdatasets", tsm_filepath, tsm_planting_filename), overwrite = TRUE)
writeRaster(tsm_growing, file.path("./existingdatasets", tsm_filepath, tsm_growing_filename), overwrite = TRUE)
