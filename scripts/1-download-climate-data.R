library(ncdf4)


# CanDCS-M6 is an upgraded, multivariate downscaled climate dataset that better preserves relationships between climate variables
# Maximum temperature
# Minimum temperature
# Total precipitation
# Growing degree days (5 degrees C)


# Open the OpenDAP NetCDF
url <- "https://pavics.ouranos.ca/twitcher/ows/proxy/thredds/dodsC/datasets/simulations/bias_adjusted/cmip6/climatedata_ca/CanDCS-M6/ensemble_members/30yAvg/annual_CanDCS-M6_climindices_ensemble_members_30yAvg.ncml"
nc <- nc_open(url)

names(nc$var)

# Pick a real variable
varname <- "ssp245_gddgrow_5"

# Get dimensions
print(nc$var[[varname]]$dim)

# Extract a small subset (e.g., lon=1, lat=1, time=1, member=1)
gdd_sample <- ncvar_get(nc, varid = varname, start = c(1, 1, 1, 1), count = c(1, 1, 1, 1))

print(gdd_sample)

# Clean up
nc_close(nc)