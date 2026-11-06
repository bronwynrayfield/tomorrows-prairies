
###################################################################################
##          DATA SUMMARY: AREA PER REPORTING UNIT BY ORIGINAL CLASSES            ##
###################################################################################

##----------------------- Alison Long - October 2026 -------------------------

# STEP 5: Calculate pixel count and area in ha for each original class in each
#         dataset. script crops/masks all input rasters to study area and for each
#         reporting unit (ecodistrict or CCS), sums the number of pixels in the class, 
#         and converts to ha. the script then aggregates area to coarser reporting units. 

# metadata file has the following fields: dataset, class, file_s3

## ---- Setup ---- ##
# Temp dirs
tmp_root <- file.path(Sys.getenv("HOME"), "r_tmp")

#remove.packages("terra")
#install.packages("terra")

# Load packages
library(terra)
library(readr)
library(dplyr)
library(paws)
library(paws.storage)
library(ggplot2)

# Set terra options
terraOptions(
  tempdir = file.path(tmp_root, "terra"),
  memfrac = 0.6, # could go down to 0.4
  progress = 1
)

# Force a fresh run instead of silently re-using old scratch outputs
FORCE_REBUILD <- TRUE

## ---- Settings ---- ##
meta_file    <- "s3://sagemaker-gst-stage.sharing/alison-long/tomorrows-prairies/processed/?.csv"
ecodist_file <- "s3://sagemaker-gst-stage.sharing/alison-long/tomorrows-prairies/processed/boundaries/national-ecological-framework-for-canada/Ecodistricts/ecodistricts_clean_rast_Oct2026.tif"
ecoprov_lookup_file <-
ccs_file     <- "s3://sagemaker-gst-stage.sharing/alison-long/tomorrows-prairies/processed/boundaries/census-consolidated-subdivisions/ccs_aligned_class_sa_rast.tif"
ccs_lookup_file  <- "s3://sagemaker-gst-stage.sharing/alison-long/tomorrows-prairies/processed/boundaries/census-consolidated-subdivisions/ccs_aligned_class_sa_rast_lookup.csv"
prov_lookup_file <- 
pilot_file   <- "s3://sagemaker-gst-stage.sharing/alison-long/tomorrows-prairies/processed/pilot-area/pilot_area_aligned_cnvc_Jan2026.tif"
sa_file      <- "s3://sagemaker-gst-stage.sharing/alison-long/tomorrows-prairies/processed/study-area/temperate_scope_aligned_cnvc_Apr2026.tif"
out_dir      <- "s3://sagemaker-gst-stage.sharing/alison-long/tomorrows-prairies/results/spatial_agreement_analysis"
log_dir      <- "s3://sagemaker-gst-stage.sharing/alison-long/tomorrows-prairies/logs"
local_dir    <- file.path(tmp_root, "data")
dir.create(local_dir, showWarnings = FALSE, recursive = TRUE)

## ---- Message log ---- ##
log_file <- file.path(log_dir, paste0(
  "spatial_agreement_log_",
  format(Sys.time(), "%Y%m%d_%H%M%S"),
  ".log"
))

# Ensure local temp file first (upload to S3 at end)
log_file_local <- file.path(
  Sys.getenv("HOME"),
  paste0("spatial_agreement_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log")
)

log_con <- file(log_file_local, open = "wt")

sink(log_con, split = TRUE)          # captures console output

on.exit({
  sink()
  close(log_con)
}, add = TRUE)

log_msg <- function(..., level = "INFO") {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  msg <- paste0(
    "[", level, "] ",
    ts,
    " : ",
    paste(..., collapse = "")
  )
  
  message(msg)
  
  cat(msg, "\n", file = log_file_local, append = TRUE)
}

log_msg("---- SESSION INFO ----")

script_name <- tryCatch({
  cmd_args <- commandArgs(trailingOnly = FALSE)
  sub("^--file=", "", cmd_args[grep("^--file=", cmd_args)])
}, error = function(e) NA_character_)

log_msg("Script name: ", script_name)

session_lines <- capture.output(sessionInfo())

for (ln in session_lines) {
  log_msg(ln)
}

log_msg("R version: ", R.version.string)
log_msg("Platform: ", R.version$platform)
log_msg("Working directory: ", getwd())
log_msg("terra version: ", as.character(packageVersion("terra")))
log_msg("readr version: ", packageVersion("readr"))
log_msg("Temp root: ", tmp_root)
log_msg("GDAL_CACHEMAX: ", Sys.getenv("GDAL_CACHEMAX"))
log_msg("FORCE_REBUILD: ", FORCE_REBUILD)
log_msg("Sys.info: ",
        paste(names(Sys.info()), Sys.info(), collapse = "; "))
log_msg("terra options:")
capture.output(print(terraOptions())) |>
  lapply(log_msg)
gc_out <- capture.output(gc())
log_msg(paste(gc_out, collapse = " | "))

log_msg("---- END SESSION INFO ----")

log_msg("Starting spatial agreement analysis")

## ---- Reusable helper function for reading rasters ---- ##
read_s3_raster <- function(s3_path) {
  
  s3_uri <- sub("^s3://", "", s3_path)
  
  bucket <- sub("/.*", "", s3_uri)
  key <- sub("^[^/]*/", "", s3_uri)
  
  tmp_file <- tempfile(fileext = tools::file_ext(key))
  
  s3 <- paws.storage::s3()
  
  obj <- s3$get_object(
    Bucket = bucket,
    Key = key
  )
  
  writeBin(obj$Body, tmp_file)
  
  terra::rast(tmp_file)
}

## ---- Load metadata ---- ##
log_msg("Reading metadata from: ", meta_file)

s3_uri <- sub("^s3://", "", meta_file)

bucket <- sub("/.*", "", s3_uri)
key <- sub("^[^/]*/", "", s3_uri)

s3 <- paws.storage::s3()

obj <- s3$get_object(Bucket = bucket, Key = key)

meta <- readr::read_csv(rawConnection(obj$Body), show_col_types = FALSE) %>%
  rename(file = file_s3)

rfiles_s3 <- setNames(meta$file, meta$raster)

stopifnot(!anyDuplicated(meta$raster))

log_msg("Metadata rows loaded: ", nrow(meta))

## ---- Load pilot area ---- ##
log_msg("Loading pilot area")

# Read raster
pilot <- read_s3_raster(pilot_file)

log_msg("Pilot area loaded")

## ---- Load study area ---- ##
log_msg("Loading study area")

# Read raster
sa <- read_s3_raster(sa_file)

log_msg("Study area loaded")

## ---- Download rasters locally ---- ##
message("Downloading rasters locally...")
log_msg("Downloading rasters locally")

download_files <- function(files_named) {
  
  if (length(files_named) == 0)
    return(setNames(character(0), character(0)))
  
  sapply(names(files_named), function(nm) {
    
    log_msg("Processing raster: ", nm)
    log_msg("S3 path: ", files_named[[nm]])
    
    crop_file <- file.path(local_dir, paste0(nm, ".tif"))
    
    if (!file.exists(crop_file) || FORCE_REBUILD) {
      
      # download and read raster from S3
      r <- read_s3_raster(files_named[[nm]])
      
      # crop + mask to pilot
      r <- crop(r, pilot)
      r <- mask(r, pilot)
      
      writeRaster(
        r,
        crop_file,
        overwrite = TRUE,
        wopt = list(
          gdal = c(
            "COMPRESS=ZSTD",
            "TILED=YES",
            "BIGTIFF=YES"
          )
        )
      )
      
      rm(r)
      gc()
    }
    
    crop_file
    
  }, USE.NAMES = TRUE)
  
}

rfiles_local <- download_files(rfiles_s3)

log_msg("All rasters downloaded locally")

# Check file size to ensure correct download
bad_files <- names(rfiles_local)[file.size(rfiles_local) < 10000]
if (length(bad_files) > 0) {
  stop(paste("These rasters failed to download correctly:",
             paste(bad_files, collapse = ", ")))
}

## ---- Confirm raster alignment ---- ##
message("Checking alignment...")
log_msg("Checking raster alignment")

template <- rast(rfiles_local[[1]])

aligned_ok <- all(sapply(rfiles_local, function(f) {
  r <- rast(f)
  compareGeom(r, template, stopOnError = FALSE)
}))

if (!aligned_ok) {
  stop("Rasters are NOT perfectly aligned. You must align them.")
}

message("Rasters confirmed aligned — skipping resampling")
log_msg("Rasters confirmed aligned — skipping resampling")

## ---- Align rasters to a common grid (if needed) ---- ##
#message("Aligning rasters...")

#template <- rast(rfiles_local[[1]])

#aligned_files <- lapply(names(rfiles_local), function(nm) {
#  f <- rfiles_local[[nm]]
#  r <- rast(f)

#  if (!compareGeom(r, template, stopOnError = FALSE)) {
#    r <- resample(r, template, method = "near")
#  }

#  out <- file.path(local_dir, paste0(nm, "_aligned.tif"))
#  writeRaster(r, out, overwrite = TRUE)
#  out
#})

#names(aligned_files) <- names(rfiles_local)

## ---- Load ecodistrict raster ---- ##
log_msg("Loading pilot-clipped ecodistrict raster")

ecodist_r <- read_s3_raster(ecodist_file)

if (!compareGeom(template, pilot, stopOnError = FALSE)) {
  stop("Pilot raster is not aligned with template raster")
}

ecodist_r <- crop(ecodist_r, pilot)
ecodist_r <- mask(ecodist_r, pilot)

ecodist_rast <- as.int(ecodist_r)

# Ensure alignment
if (!compareGeom(ecodist_rast, template, stopOnError = FALSE)) {
  stop("Ecodistrict raster is not aligned with input rasters")
}

log_msg("Loaded ecodistrict raster")

## ---- Load census consolidated subdivision raster ---- ##
log_msg("Loading pilot-clipped census consolidated subdivision raster")

ccs_r <- read_s3_raster(ccs_file)

if (!compareGeom(template, pilot, stopOnError = FALSE)) {
  stop("Pilot raster is not aligned with template raster")
}

ccs_r <- crop(ccs_r, pilot)
ccs_r <- mask(ccs_r, pilot)

ccs_rast <- as.int(ccs_r)

# Ensure alignment
if (!compareGeom(ccs_rast, template, stopOnError = FALSE)) {
  stop("Census consolidated subdivision raster is not aligned with input rasters")
}

# load CCS lookup file
s3_uri <- sub("^s3://", "", ccs_lookup_file)

bucket <- sub("/.*", "", s3_uri)
key <- sub("^[^/]*/", "", s3_uri)

obj <- s3$get_object(Bucket = bucket, Key = key)

ccs_lookup <- readr::read_csv(rawConnection(obj$Body), show_col_types = FALSE)

stopifnot(all(c("ID", "CCSUID") %in% names(ccs_lookup)))

log_msg("Loaded census consolidated subdivision raster")

## ---- Build raster stack ---- ##
log_msg("Building raster stack")

r_stack <- rast(rfiles_local)
n_rasters <- nlyr(r_stack)
names(r_stack) <- names(rfiles_local)
ras_names <- names(r_stack)

log_msg("Raster stack built")

## ---- Initialize ecodistricts ---- ##
log_msg("Initializing ecodistricts")

# Get unique district IDs 
freq_table_eco <- freq(ecodist_rast)
district_vals <- as.integer(freq_table_eco$value)

if (anyDuplicated(district_vals)) {
  stop("Duplicate district IDs found")
}

n_districts <- length(district_vals)

# Map district IDs → index
district_index <- setNames(seq_along(district_vals), district_vals)

log_msg("Ecodistricts ready")

## ---- Initialize CCS ---- ##
log_msg("Initializing CCS")

# Get unique IDs 
freq_table_ccs <- freq(ccs_rast)
ccs_vals <- as.integer(freq_table_ccs$value)

missing_ids <- setdiff(ccs_vals, ccs_lookup$ID)

if (length(missing_ids) > 0) {
  stop("Lookup missing CCS IDs: ", paste(missing_ids, collapse = ", "))
}

if (anyDuplicated(ccs_vals)) {
  stop("Duplicate CCS IDs found")
}

n_ccs <- length(ccs_vals)

# Map CCS IDs → index
ccs_index <- setNames(seq_along(ccs_vals), ccs_vals)

log_msg("CCS ready")

## ---- Load lookup tables ---- ##
log_msg("Loading lookup tables")

read_s3_csv <- function(s3_path) {
  
  s3_uri <- sub("^s3://", "", s3_path)
  
  bucket <- sub("/.*", "", s3_uri)
  key <- sub("^[^/]*/", "", s3_uri)
  
  obj <- s3$get_object(
    Bucket = bucket,
    Key = key
  )
  
  read_csv(
    rawConnection(obj$Body),
    show_col_types = FALSE
  )
}

ecoprov_lookup <- read_s3_csv(ecoprov_lookup_file)
prov_lookup <- read_s3_csv(prov_lookup_file)

log_msg(
  "Loaded ",
  nrow(ecoprov_lookup),
  " ecodistrict lookup rows and ",
  nrow(prov_lookup),
  " province lookup rows"
)

## ---- Area constants ---- ##
cell_area_ha <- prod(res(template)) / 10000

log_msg(
  "Cell area (ha): ",
  round(cell_area_ha, 4)
)

## ---- Function: summarize one raster by reporting unit ---- ##
summarize_by_unit <- function(class_raster,
                              unit_raster,
                              unit_field) {
  
  x <- c(unit_raster, class_raster)
  
  names(x) <- c("unit_id", "class_value")
  
  out <- as.data.frame(
    terra::freq(
      x,
      digits = 0,
      value = TRUE
    )
  )
  
  names(out) <- c(
    unit_field,
    "class_value",
    "pixel_count"
  )
  
  out$area_ha <- out$pixel_count * cell_area_ha
  
  out
}

## ---- Summaries by ecodistrict ---- ##
log_msg("Calculating ecodistrict summaries")

eco_results <- vector(
  "list",
  length(ras_names)
)

for (i in seq_along(ras_names)) {
  
  nm <- ras_names[i]
  
  log_msg(
    "Processing ecodistrict summary: ",
    nm,
    " (",
    i,
    "/",
    length(ras_names),
    ")"
  )
  
  r <- rast(rfiles_local[[nm]])
  
  tmp <- summarize_by_unit(
    class_raster = r,
    unit_raster = ecodist_rast,
    unit_field = "ECODISTRICT"
  )
  
  tmp$raster <- nm
  
  eco_results[[i]] <- tmp
  
  rm(r, tmp)
  gc()
}

eco_area <- bind_rows(eco_results)

log_msg(
  "Ecodistrict summary rows: ",
  nrow(eco_area)
)

## ---- Summaries by CCS ---- ##
log_msg("Calculating CCS summaries")

ccs_results <- vector(
  "list",
  length(ras_names)
)

for (i in seq_along(ras_names)) {
  
  nm <- ras_names[i]
  
  log_msg(
    "Processing CCS summary: ",
    nm,
    " (",
    i,
    "/",
    length(ras_names),
    ")"
  )
  
  r <- rast(rfiles_local[[nm]])
  
  tmp <- summarize_by_unit(
    class_raster = r,
    unit_raster = ccs_rast,
    unit_field = "ID"
  )
  
  tmp$raster <- nm
  
  ccs_results[[i]] <- tmp
  
  rm(r, tmp)
  gc()
}

ccs_area <- bind_rows(ccs_results)

ccs_area <- ccs_area %>%
  left_join(
    ccs_lookup %>%
      select(ID, CCSUID),
    by = "ID"
  )

log_msg(
  "CCS summary rows: ",
  nrow(ccs_area)
)

## ---- Join raster metadata ---- ##
log_msg("Joining raster metadata")

meta_join <- meta %>%
  select(
    raster,
    dataset,
    class
  )

eco_area <- eco_area %>%
  left_join(meta_join, by = "raster")

ccs_area <- ccs_area %>%
  left_join(meta_join, by = "raster")

## ---- Join ecodistrict hierarchy ---- ##
eco_area <- eco_area %>%
  left_join(
    ecoprov_lookup,
    by = "ECODISTRICT"
  )

## ---- Ecoregion summaries ---- ##
ecoregion_area <- eco_area %>%
  group_by(
    raster,
    dataset,
    class,
    class_value,
    ECOREGION
  ) %>%
  summarise(
    pixel_count = sum(pixel_count),
    area_ha = sum(area_ha),
    .groups = "drop"
  )

## ---- Ecoprovince summaries ---- ##
ecoprov_area <- eco_area %>%
  group_by(
    raster,
    dataset,
    class,
    class_value,
    ECOPROVINCE
  ) %>%
  summarise(
    pixel_count = sum(pixel_count),
    area_ha = sum(area_ha),
    .groups = "drop"
  )

## ---- Ecozone summaries ---- ##
ecozone_area <- eco_area %>%
  group_by(
    raster,
    dataset,
    class,
    class_value,
    ECOZONE
  ) %>%
  summarise(
    pixel_count = sum(pixel_count),
    area_ha = sum(area_ha),
    .groups = "drop"
  )

## ---- Province summaries ---- ##
ccs_area <- ccs_area %>%
  left_join(
    prov_lookup,
    by = "CCSUID"
  )

province_area <- ccs_area %>%
  group_by(
    raster,
    dataset,
    class,
    class_value,
    PROVINCE
  ) %>%
  summarise(
    pixel_count = sum(pixel_count),
    area_ha = sum(area_ha),
    .groups = "drop"
  )

## ---- Write outputs ---- ##
log_msg("Writing output CSV files")

eco_out_file <- file.path( local_dir, "area_by_ecodistrict.csv")
ccs_out_file <- file.path( local_dir, "area_by_ccs.csv")
ecoregion_out_file <- file.path( local_dir, "area_by_ecoregion.csv")
ecoprov_out_file <- file.path( local_dir, "area_by_ecoprovince.csv")
ecozone_out_file <- file.path( local_dir, "area_by_ecozone.csv")
province_out_file <- file.path( local_dir, "area_by_province.csv")

write_csv(eco_area, eco_out_file)
write_csv(ccs_area, ccs_out_file)
write_csv(ecoregion_area, ecoregion_out_file)
write_csv(ecoprov_area, ecoprov_out_file)
write_csv(ecozone_area, ecozone_out_file)
write_csv(province_area, province_out_file)

log_msg("Finished writing CSV outputs")

## ---- Upload outputs to s3 ---- ##

# define function
upload_file_s3 <- function(local_file, s3_uri) {
  
  s3_uri2 <- sub("^s3://", "", s3_uri)
  
  bucket <- sub("/.*", "", s3_uri2)
  key <- sub("^[^/]*/", "", s3_uri2)
  
  s3$put_object(
    Bucket = bucket,
    Key = key,
    Body = file(local_file, "rb")
  )
}

upload_file_s3( eco_out_file, file.path(out_dir, "area_by_ecodistrict.csv"))
upload_file_s3(ccs_out_file, file.path(out_dir, "area_by_ccs.csv"))
upload_file_s3( ecoregion_out_file, file.path(out_dir, "area_by_ecoregion.csv"))
upload_file_s3(ecoprov_out_file, file.path(out_dir, "area_by_ecoprovince.csv"))
upload_file_s3(ecozone_out_file, file.path(out_dir, "area_by_ecozone.csv"))
upload_file_s3(province_out_file, file.path(out_dir, "area_by_province.csv"))

log_msg("All output files uploaded")

upload_file_s3(log_file_local, file.path(out_dir, basename(log_file_local)))

message("Log written to: ", log_file)


