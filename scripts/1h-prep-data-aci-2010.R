# =============================================================================
# 1h-prep-data-aci-2010.R
# =============================================================================
#
# Purpose:
#   Downloads AAFC Annual Crop Inventory (ACI) 2010 GeoTIFFs for ON, MB, SK,
#   merges them into a single raster covering the modelling extent, and
#   reprojects/crops to match the predictor stack grid. Used to distinguish
#   annual cropland from tame grassland within the LUTS "Cropland" class in
#   2010, since LUTS does not separate these two classes in that year.
#
# Inputs:
#   - AAFC ACI 2010 pre-packaged GeoTIFF zips (one per province) from Canada's
#     Open Data portal:
#       https://agriculture.canada.ca/atlas/data_donnees/annualCropInventory/
#       data_donnees/tif/2010/aci_2010_{prov}_v2.zip
#     Provinces: on, mb, sk
#
# Outputs:
#   - aci_2010_{prov}_v2.tif : unzipped raw GeoTIFFs, cropped to a buffered
#       modelling extent, in existingdatasets/aci/2010/
#   - aci_2010_merged.tif    : merged ON+MB+SK raster, native CRS/resolution
#   - aci_2010_30m.tif       : merged raster reprojected (nearest neighbour)
#       and cropped to the modelling extent, aligned to the 30m predictor
#       stack grid
#
# Notes:
#   - ACI is a classified (categorical) raster, so reprojection uses
#     method = "near" throughout to avoid inventing intermediate class codes
#   - Each tile is cropped to a buffered AOI immediately after download, to
#     cut file size before merge/project (same pattern as PSSSM prep)
#   - merge() is used rather than mosaic() since the three province tiles are
#     adjacent and non-overlapping
#   - Unzipped filenames are assumed to match the zip name (e.g.
#     aci_2010_on_v2.tif) -- check this on first run and adjust tif_name
#     below if the archive structure differs
#
# Dependencies:
#   - 0-setup.R (parameters, file paths, functions)
# =============================================================================
source(file.path("scripts/landcover", "0-setup.R"))

# --- Configuration ----
base_url <- "https://agriculture.canada.ca/atlas/data_donnees/annualCropInventory/data_donnees/tif/2010"
provs <- c("ab", "mb", "sk")

raw_dir <- capitalize_first(file.path("./existingdatasets", aci_filepath, "raw"))
cropped_dir <- capitalize_first(file.path("./existingdatasets", aci_filepath, "cropped"))
aci_nodata <- 0  # 0 = real background/nodata in these tiles

# Load modelling extent (used to crop tiles immediately after download)
me <- st_read(file.path(sa_filepath, me_filename))

# --- Download, unzip, crop tiles ----
for (prov in provs) {
  zip_name <- sprintf("aci_2010_%s_v2.zip", prov)
  zip_dest <- file.path(raw_dir, zip_name)
  download.file(file.path(base_url, zip_name), zip_dest, mode = "wb", quiet = TRUE)
  unzip(zip_dest, exdir = raw_dir)
  
  tif_name <- sprintf("aci_2010_%s_v2.tif", prov)
  tif_path <- file.path(raw_dir, tif_name)
}

# --- Crop each tile to the modelling extent, with NoData preserved ----
tif_raw   <- file.path(raw_dir, sprintf("aci_2010_%s_v2.tif", provs))
tif_crop  <- file.path(cropped_dir, sprintf("aci_2010_%s_v2_crop.tif", provs))

for (i in seq_along(provs)) {
  r <- rast(tif_raw[i])
  coltab(r) <- NULL                   # strip palette -> plain classified raster
  r <- subst(r, aci_nodata, NA)   # 0 is background, not real classes

  me_native <- st_transform(me, crs(r))
  r_crop <- crop(r, vect(st_buffer(me_native, 5000)))
  plot(r_crop)
  writeRaster(r_crop, tif_crop[i], overwrite = TRUE,
              datatype = "INT1U", NAflag = 0,
              gdal = "COMPRESS=DEFLATE")
}

# --- Merge province tiles ----
tiles <- lapply(tif_crop, rast)
aci_2010 <- do.call(merge, c(tiles, list(method = "near")))
plot(aci_2010)
writeRaster(aci_2010, file.path("./existingdatasets", aci_filepath, aci_filename),
            datatype = "INT1U", NAflag = 0,
            gdal = "COMPRESS=DEFLATE")

