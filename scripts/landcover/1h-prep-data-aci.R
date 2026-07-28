# =============================================================================
# 1h-prep-data-aci.R
# =============================================================================
#
# Purpose:
#   Downloads AAFC Annual Crop Inventory (ACI) GeoTIFFs for AB, SK and MB across
#   a configurable range of years, merges the province tiles into a single
#   raster per year, and reprojects/crops each to the predictor stack grid.
#   Used to distinguish annual cropland from tame grassland (ACI class 122)
#   within the LUTS "Cropland" class, and to characterise the multi-year
#   tame-grass / annual-crop rotation cycle.
#
# Inputs:
#   - AAFC ACI pre-packaged GeoTIFF zips, one per province per year:
#       https://agriculture.canada.ca/atlas/data_donnees/annualCropInventory/
#       data_donnees/tif/{year}/{zip}
#     The zip filename pattern is NOT stable across years (some carry a _v2
#     suffix), so candidate patterns are probed by HTTP HEAD before download.
#   - Official ACI class legend CSV (used for the legend drift check)
#   - Modelling extent polygon (me) and raster (me_r)
#
# Outputs:
#   - existingdatasets/aci/raw/{year}/       : zips and unzipped provincial tifs
#   - existingdatasets/aci/cropped/{year}/   : tiles cropped to buffered ME
#   - existingdatasets/aci/aci_{year}_merged.tif  : merged AB+SK+MB, native grid
#   - existingdatasets/aci/30m/aci_{year}_30m_modelling-extent.tif : aligned to me_r
#   - existingdatasets/aci/300m/aci_{year}_300m_modelling-extent.tif : aligned to 300 m me_r
#   - output/.../aci_url_index.csv          : resolved URLs and file sizes
#   - output/.../aci_class_frequencies.csv  : class counts by year (drift check)
#
# Notes:
#   - ACI is categorical: method = "near" everywhere. Never bilinear.
#   - 2009 and 2010 are 56 m; 2011 onward are 30 m. Projecting the early years
#     to the 30 m grid replicates each 56 m value across ~3-4 cells.
#   - The 300 m outputs are near-sampled from native resolution onto a fixed
#     grid, so a given 300 m cell samples the same ground location every year.
#     That makes them a valid 1-in-100 systematic sample for multi-year
#     trajectory work. Do NOT use modal aggregation for that purpose: modal
#     changes what a cell means between a 56 m year and a 30 m year.
#   - Class codes are not guaranteed stable across the series. Early
#     Prairie-only years may lean on generic classes (120 Agriculture,
#     121 Cropland) where later years resolve wheat/canola separately. Check
#     aci_class_frequencies.csv before applying a fixed annual-crop vector.
#
#
# Class code stability (checked 2009-2012, modelling extent):
#   - Codes 110 (Grassland) and 122 (Pasture/forages) are stable and large
#     throughout (~0.9-1.1M cells at 300 m), so the tame/native distinction is
#     supported across the series.
#   - Wheat is coded 140 in 2009, 2011, 2012 but 146 (Spring wheat) in 2010.
#     Winter wheat (145) appears only in 2009-2010. Generic 132 (Cereals)
#     appears from 2011. Derive annual-crop codes from the legend CSV rather
#     than hard-coding: a fixed vector keyed on 140 loses ~1M wheat cells in 2010.
#   - 120 (Agriculture, undifferentiated) is negligible (<100 cells/year), so
#     the ambiguity between annual and perennial crop in that class is not a
#     practical concern here.
#   - 130 (Too wet to be seeded) is a wet-year artifact: ~467k cells in 2010 and
#     ~577k in 2011 (prairie flood years) vs ~15-55k otherwise. Treat 130 as
#     MISSING in any persistence test, not as annual crop. Within the tame pool
#     it affects only 3.9% (2010) and 7.5% (2011) of pixels, since
#     "too wet to be seeded" is a cropland concept and rarely lands on forage.
#   - Fewer than 100 cells per year (<0.01% of the tame pool) are NA in any
#     single year; these are scattered, not a shared coverage seam.
#
# Dependencies:
#   - 0-setup.R (parameters, file paths, functions)
#   - curl (HEAD probe only)
# =============================================================================

source(file.path("scripts/landcover", "0-setup.R"))

# ---- Configuration ------------------------------------------------------
years <- 2009:2012
provs <- c("ab", "sk", "mb")

base_url <- "https://agriculture.canada.ca/atlas/data_donnees/annualCropInventory/data_donnees/tif"
legend_url <- paste0("https://agriculture.canada.ca/atlas/data_donnees/annualCropInventory/",
                     "supportdocument_documentdesupport/",
                     "aci_crop_classifications_iac_classifications_des_cultures.csv")

# Candidate zip name patterns, tried in order. Add to this if a year 404s.
zip_patterns <- c("aci_%d_%s_v4.zip", "aci_%d_%s_v2.zip", "aci_%d_%s.zip", "ACI_%d_%s.zip")

aci_nodata <- 0      # 0 = real background/nodata in these tiles
buffer_m   <- 5000   # buffer on modelling extent before cropping tiles
res_coarse <- 300    # coarse grid for trajectory diagnostics

raw_dir     <- file.path("./existingdatasets", aci_filepath, "raw")
cropped_dir <- file.path("./existingdatasets", aci_filepath, "cropped")
out_30_dir  <- file.path("./processeddatasets", aci_filepath, "30m")
out_300_dir <- file.path("./processeddatasets", aci_filepath, "300m")
summary_dir <- "./output/landcover/aci"

for (d in c(raw_dir, cropped_dir, out_30_dir, out_300_dir, summary_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# ---- Load modelling extent ----------------------------------------------
me   <- st_read(file.path(sa_filepath, me_filename))
me_r <- rast(file.path(sa_filepath, me_raster_filename))

me_r_coarse <- aggregate(me_r, fact = res_coarse / res, fun = "modal")

# ---- Resolve zip URLs before downloading anything -----------------------
head_probe <- function(u) {
  h <- curl::new_handle(nobody = TRUE, followlocation = TRUE)
  r <- try(curl::curl_fetch_memory(u, handle = h), silent = TRUE)
  if (inherits(r, "try-error") || r$status_code != 200) return(NA_real_)
  as.numeric(curl::parse_headers_list(r$headers)[["content-length"]] %||% NA)
}

url_index <- do.call(rbind, lapply(years, function(yr) {
  do.call(rbind, lapply(provs, function(pv) {
    for (p in zip_patterns) {
      zip_name <- sprintf(p, yr, pv)
      u <- file.path(base_url, yr, zip_name)
      sz <- head_probe(u)
      if (!is.na(sz)) {
        return(data.frame(year = yr, prov = pv, zip = zip_name,
                          url = u, size_mb = round(sz / 1e6, 1)))
      }
    }
    data.frame(year = yr, prov = pv, zip = NA, url = NA, size_mb = NA)
  }))
}))

write.csv(url_index, file.path(summary_dir, "aci_url_index.csv"), row.names = FALSE)
print(url_index)

# ---- Download, unzip, crop ----------------------------------------------
for (i in seq_len(nrow(url_index))) {
  yr <- url_index$year[i]
  pv <- url_index$prov[i]
  
  yr_raw     <- file.path(raw_dir, yr)
  yr_cropped <- file.path(cropped_dir, yr)
  dir.create(yr_raw,     recursive = TRUE, showWarnings = FALSE)
  dir.create(yr_cropped, recursive = TRUE, showWarnings = FALSE)
  
  crop_path <- file.path(yr_cropped, sprintf("aci_%d_%s_crop.tif", yr, pv))
  if (file.exists(crop_path)) next
  
  zip_dest <- file.path(yr_raw, url_index$zip[i])
  if (!file.exists(zip_dest)) {
    download.file(url_index$url[i], zip_dest, mode = "wb", quiet = TRUE)
  }
  
  # Read the tif name out of the archive rather than assuming it matches
  contents <- unzip(zip_dest, list = TRUE)$Name
  tif_name <- grep("\\.tif$", contents, value = TRUE, ignore.case = TRUE)[1]
  stopifnot(!is.na(tif_name))
  unzip(zip_dest, files = tif_name, exdir = yr_raw, junkpaths = TRUE)
  
  r <- rast(file.path(yr_raw, basename(tif_name)))
  coltab(r) <- NULL                       # strip palette -> plain classified raster
  
  me_native <- st_transform(me, crs(r))
  r <- crop(r, vect(st_buffer(me_native, buffer_m)))   # crop before subst
  r <- subst(r, aci_nodata, NA)
  
  writeRaster(r, crop_path, overwrite = TRUE,
              datatype = "INT1U", NAflag = 0,
              gdal = "COMPRESS=DEFLATE")
}

# ---- Merge provinces and align to predictor grid ------------------------
for (yr in years) {
  merged_path <- file.path("./existingdatasets", aci_filepath, sprintf("aci_%d_merged.tif", yr))
  path_30     <- file.path(out_30_dir,  sprintf("aci_%d_30m.tif", yr))
  path_300    <- file.path(out_300_dir, sprintf("aci_%d_%dm.tif", yr, res_coarse))
  
  tif_crop <- file.path(cropped_dir, yr, sprintf("aci_%d_%s_crop.tif", yr, provs))
   
  tiles <- sprc(lapply(tif_crop, rast))
  aci_merged <- merge(tiles, method = "near")
   
  writeRaster(aci_merged, merged_path, overwrite = TRUE,
              datatype = "INT1U", NAflag = 0, gdal = "COMPRESS=DEFLATE")
  
  # Project straight from native resolution onto each target grid
  aci_30_proj <- project(aci_merged, me_r, method = "near")
  aci_30 <- mask(crop(aci_30_proj, me_r), me_r)
  write_rast_tagged(aci_30, "./processeddatasets", file.path(aci_filepath, "30m"),  sprintf("aci_%d_30m", yr), "modelling-extent", datatype = "INT1U", gdal = "COMPRESS=DEFLATE")
  
  aci_300_proj <- project(aci_merged, me_r_coarse, method = "near")
  aci_300 <- mask(crop(aci_300_proj, me_r_coarse), me_r_coarse)
  write_rast_tagged(aci_300, "./processeddatasets", file.path(aci_filepath, "300m"),  sprintf("aci_%d_30m", yr), "modelling-extent", datatype = "INT1U", gdal = "COMPRESS=DEFLATE")
  
}

# ---- Legend drift check -------------------------------------------------
download.file(legend_url, file.path(summary_dir, "aci_class_legend.csv"),
              mode = "wb", quiet = TRUE)
legend <- read.csv(file.path(summary_dir, "aci_class_legend.csv"),
                   fileEncoding = "CP1252", check.names = FALSE)
legend <- legend[, c("Code", "Label")]
names(legend) <- c("code", "label")

class_freq <- do.call(rbind, lapply(years, function(yr) {
  f <- freq(rast(file.path(out_300_dir, sprintf("aci_%d_%dm.tif", yr, res_coarse))))
  data.frame(year = yr, code = f$value, count = f$count)
}))

class_wide <- class_freq |>
  tidyr::pivot_wider(names_from = year, values_from = count, values_fill = 0) |>
  arrange(code)

write.csv(class_wide, file.path(summary_dir, "aci_class_frequencies.csv"), row.names = FALSE)
print(class_wide, n = 100)

# Codes that appear in some years but not others -- the drift to watch for
present <- class_wide[, as.character(years)] > 0
class_wide$n_years <- rowSums(present)
print(class_wide[class_wide$n_years < length(years), c("code", "n_years")], n = 100)
