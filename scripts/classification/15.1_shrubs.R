
###################################################################################
##               ANALYSIS: ASSESSING SHRUB AGREEMENT AND AGREA                   ##
###################################################################################

##----------------------- Alison Long - Nov 2026 -------------------------

# STEP 15: Calculate pixelwise counts of shrub class and

# excluding datasets that have a shrub class but are non-independent:
#    - excluded CEC NALCMS (included NRCan LCC)
#    - excluded ACI 2024 (included ACI 2025)

# load packages
library(readxl)
library(terra)
library(sf)
library(tidyverse)

# if running on SageMaker
library(paws)

################## SETUP FOR RUNNING ON NIMBUS #################

# Set the temporary directory to a drive with sufficient space
terraOptions(tempdir = "R:/CanadianGrasslands/WG_classification/scratch")

####--- DIRECTORIES ---####

# pilot area
pa_dir <- "R:/CanadianGrasslands/processeddatasets/pilotarea"

# study area
sa_dir <- "R:/CanadianGrasslands/processeddatasets/studyarea"

# processed data
proc_data_dir <- "R:/CanadianGrasslands/WG_Classification/processed"

####--- Files ---####
# Pilot Area
pa_file <- file.path(pa_dir, "pilot_area_aligned_cnvc_Jan2026.tif")

# Study Area
sa_file <- file.path(sa_dir, "class_study_area_aligned_cnvc_Dec2025.tif")

# file paths and metadata
metadata_file <- file.path(proc_data_dir, "metadata_shrub_classes.xlsx")

################## SETUP FOR RUNNING ON SAGEMAKER #################

# temporary directory
terraOptions(tempdir = tempdir())

####--- DIRECTORIES ---####

# processed data
proc_data_dir <- "s3://sagemaker-gst-stage.sharing/alison-long/tomorrows-prairies/processed"

####--- Files ---####

# file paths and metadata
metadata_file <- file.path(proc_data_dir, "metadata/metadata_shrub_classes.xlsx")

save_object(
  object = "alison-long/tomorrows-prairies/processed/metadata/metadata_shrub_classes.xlsx",
  bucket = "sagemaker-gst-stage.sharing",
  file = "/tmp/metadata_shrub_classes.xlsx"
)

metadata <- read_excel("/tmp/metadata_shrub_classes.xlsx")

############# PRESENCE AND ABSENCE COUNTS ############

#----------- Load data
metadata <- read_excel(metadata_file)

# Check required columns exist
stopifnot(all(c("raster_name", "nimbus_file_path") %in% names(metadata)))

# Create named vector of raster paths
raster_paths <- setNames(
  metadata$nimbus_file_path,
  metadata$raster_name
)

# Load all rasters into a named list
rasters <- lapply(raster_paths, rast)

# List loaded raster names
names(rasters)

# create raster stack
raster_stack <- rast(metadata$nimbus_file_path)

# Apply layer names from metadata
names(raster_stack) <- metadata$raster_name

#----------- Check that all rasters have identical geometry
# Use the first raster as the reference
ref <- rasters[[1]]

alignment_check <- tibble(
  raster = names(rasters),
  aligned = sapply(
    rasters,
    function(x) compareGeom(
      ref,
      x,
      crs = TRUE,
      ext = TRUE,
      rowcol = TRUE,
      res = TRUE,
      stopOnError = FALSE
    )
  )
)

alignment_check

#------------- plot rasters
# reduce resolution for plotting
plot(raster_stack)

#plot(aggregate(raster_stack, fact = 200), nc = 3)

#------ confirm raster values
bind_rows(
  lapply(names(rasters), function(nm) {
    freq(rasters[[nm]]) |>
      as.data.frame() |>
      mutate(raster = nm)
  }))

#    layer value      count        raster
# 1      1     0 2662915935 aci2025_shrub
# 2      1     1   88136783 aci2025_shrub
# 3      1     0  692094934    pudl_shrub
# 4      1     1   20456239    pudl_shrub
# 5      1     0 4121891739     lcc_shrub
# 6      1     1  309979779     lcc_shrub
# 7      1     0 4168449603     esa_shrub
# 8      1     1  264149618     esa_shrub
# 9      1     0  260734390     pli_shrub
# 10     1     1    6369889     pli_shrub
# 11     1     0   78551449    plvi_shrub # changed due to gte2015?
# 12     1     1    2593708    plvi_shrub
# 13     1     0  120628498     mgi_shrub
# 14     1     1   14920433     mgi_shrub

#-------- pixelwise presence counts

# Count number of rasters with value == 1 at each pixel
presence_count <- app(
  raster_stack,
  fun = function(x) sum(x == 1, na.rm = TRUE))

presence_count2 <- app(
  raster_stack,
  fun = function(x) {
    if (all(is.na(x))) NA else sum(x == 1, na.rm = TRUE)})

# Plot
plot(presence_count,
     main = "Number of datasets indicating shrub presence")


# Count shrub absences (value = 0)
absence_count <- app(
  raster_stack,
  fun = function(x) sum(x == 0, na.rm = TRUE)
)

# Count datasets contributing information (0 or 1)
coverage_count <- app(
  raster_stack,
  fun = function(x) sum(!is.na(x))
)

absence_count <- app(
  raster_stack,
  fun = function(x) {
    if (all(is.na(x))) NA else sum(x == 0, na.rm = TRUE)
  }
)

coverage_count <- app(
  raster_stack,
  fun = function(x) {
    if (all(is.na(x))) NA else sum(!is.na(x))
  }
)

# check that all are equal
all.equal(
  values(coverage_count),
  values(presence_count) + values(absence_count))

#----- calculate proportions

presence_prop <- presence_count / coverage_count
absence_prop  <- absence_count / coverage_count


