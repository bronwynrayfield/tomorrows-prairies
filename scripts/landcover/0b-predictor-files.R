# 0b-predictor-files.R
# Named vector of predictor layer paths, shared by 6-screen-collinearity-predictors.R
# and 6b-plot-predictors.R. Sourced after 0-setup.R.

predictor_files <- c(
  # climate
  prtot      = tagged_path("./processeddatasets", climate_filepath, prtot_hist_filename, "modelling-extent"),
  tn_mean    = tagged_path("./processeddatasets", climate_filepath, tn_mean_hist_filename, "modelling-extent"),
  tx_mean    = tagged_path("./processeddatasets", climate_filepath, tx_mean_hist_filename, "modelling-extent"),
  gddgrow_5  = tagged_path("./processeddatasets", climate_filepath, gddgrow_5_hist_filename, "modelling-extent"),
  nrcdd      = tagged_path("./processeddatasets", climate_filepath, nr_cdd_hist_filename, "modelling-extent"),
  
  # soil
  oc           = tagged_path("./processeddatasets", soil_filepath, oc_filename, "modelling-extent"),
  bd           = tagged_path("./processeddatasets", soil_filepath, bd_filename, "modelling-extent"),
  clay         = tagged_path("./processeddatasets", soil_filepath, clay_filename, "modelling-extent"),
  sand         = tagged_path("./processeddatasets", soil_filepath, sand_filename, "modelling-extent"),
  silt         = tagged_path("./processeddatasets", soil_filepath, silt_filename, "modelling-extent"),
  ph           = tagged_path("./processeddatasets", soil_filepath, ph_filename, "modelling-extent"),
  drainage     = tagged_path("./processeddatasets", dss_filepath, dss_drainage_filename, "modelling-extent"),
  watertbl     = tagged_path("./processeddatasets", dss_filepath, dss_watertbl_filename, "modelling-extent"),
  kp0          = tagged_path("./processeddatasets", dss_filepath, dss_kp0_filename, "modelling-extent"),
  kp33         = tagged_path("./processeddatasets", dss_filepath, dss_kp33_filename, "modelling-extent"),
  kp1500       = tagged_path("./processeddatasets", dss_filepath, dss_kp1500_filename, "modelling-extent"),
  ldepth       = tagged_path("./processeddatasets", dss_filepath, dss_ldepth_filename, "modelling-extent"),
  ldepth_max   = tagged_path("./processeddatasets", dss_filepath, dss_ldepth_max_filename, "modelling-extent"),
  psssm_mean   = tagged_path("./processeddatasets", psssm_filepath, psssm_mean_filename, "modelling-extent"),
  psssm_sd     = tagged_path("./processeddatasets", psssm_filepath, psssm_sd_filename, "modelling-extent"),
  tsm_planting = tagged_path("./processeddatasets", tsm_filepath, tsm_planting_filename, "modelling-extent"),
  tsm_growing  = tagged_path("./processeddatasets", tsm_filepath, tsm_growing_filename, "modelling-extent"),
  
  # topography
  elev      = tagged_path("./processeddatasets", dem_filepath, dem_filename, "modelling-extent"),
  tri       = tagged_path("./processeddatasets", dem_filepath, tri_filename, "modelling-extent"),
  tpi       = tagged_path("./processeddatasets", dem_filepath, tpi_filename, "modelling-extent"),
  srr       = tagged_path("./processeddatasets", dem_filepath, srr_filename, "modelling-extent"),
  northness = tagged_path("./processeddatasets", dem_filepath, northness_filename, "modelling-extent"),
  eastness  = tagged_path("./processeddatasets", dem_filepath, eastness_filename, "modelling-extent"),
  slope     = tagged_path("./processeddatasets", dss_filepath, dss_slope_filename, "modelling-extent")
  
  # landscape context -- uncomment (and add a trailing comma above) once paths exist
  # grs        = tagged_path("./processeddatasets", ..., "modelling-extent"),
  # crop       = tagged_path("./processeddatasets", ..., "modelling-extent"),
  # dev        = tagged_path("./processeddatasets", ..., "modelling-extent"),
  # fors       = tagged_path("./processeddatasets", ..., "modelling-extent"),
  # wtr        = tagged_path("./processeddatasets", ..., "modelling-extent"),
  # bare       = tagged_path("./processeddatasets", ..., "modelling-extent"),
  # ease       = tagged_path("./processeddatasets", ..., "modelling-extent"),
  # past       = tagged_path("./processeddatasets", ..., "modelling-extent"),
  
  # socio-economic
  # net_return = tagged_path("./processeddatasets", ..., "modelling-extent"),
  # d_road     = tagged_path("./processeddatasets", ..., "modelling-extent"),
  # d_rail     = tagged_path("./processeddatasets", ..., "modelling-extent"),
  # d_irig     = tagged_path("./processeddatasets", ..., "modelling-extent"),
  # d_city     = tagged_path("./processeddatasets", ..., "modelling-extent"),
  # popdens    = tagged_path("./processeddatasets", ..., "modelling-extent")
)

# display names + units for figure panels -- verify units before submission
predictor_labels <- c(
  prtot        = "Total precipitation (mm)",
  tn_mean      = "Mean daily min. temperature (\u00B0C)",
  tx_mean      = "Mean daily max. temperature (\u00B0C)",
  gddgrow_5    = "Growing degree days > 5\u00B0C",
  nrcdd        = "Consecutive dry days (days)",
  oc           = "Soil organic carbon (%)",
  bd           = "Bulk density (g/cm\u00B3)",
  clay         = "Clay (%)",
  sand         = "Sand (%)",
  silt         = "Silt (%)",
  ph           = "Soil pH",
  drainage     = "Drainage class",
  watertbl     = "Water table depth (cm)",
  kp0          = "Water retention, 0 kPa (%)",
  kp33         = "Water retention, 33 kPa (%)",
  kp1500       = "Water retention, 1500 kPa (%)",
  ldepth       = "A horizon thickness (cm)",
  ldepth_max   = "Total soil depth (cm)",
  psssm_mean   = "Sat. surface soil moisture, mean (%)",
  psssm_sd     = "Sat. surface soil moisture, SD (%)",
  tsm_planting = "Total soil moisture, planting (mm)",
  tsm_growing  = "Total soil moisture, growing (mm)",
  elev         = "Elevation (m)",
  tri          = "Terrain ruggedness index (m)",
  tpi          = "Topographic position index (m)",
  srr          = "Surface relief ratio",
  northness    = "Northness",
  eastness     = "Eastness",
  slope        = "Slope (%)"
)