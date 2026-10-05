
library(terra)

# convert input rasters to COGs
files <- list.files("rasters/", pattern = "\\.tif$", full.names = TRUE)

for (f in files) {
  
  out <- file.path("cogs", paste0(tools::file_path_sans_ext(basename(f)), "_cog.tif"))
  
  r <- rast(f)
  
  writeRaster(r, out, filetype = "COG",
    gdal = c("COMPRESS=LZW", "BLOCKSIZE=256"), overwrite = TRUE)
  
  message("Finished: ", basename(f))
}

# convert outut rasters to COGs
files <- list.files("rasters/", pattern = "\\.tif$", full.names = TRUE)

for (f in files) {
  
  out <- file.path("cogs", paste0(tools::file_path_sans_ext(basename(f)), "_cog.tif"))
  
  r <- rast(f)
  
  writeRaster(r, out, filetype = "COG",
              gdal = c("COMPRESS=LZW", "BLOCKSIZE=256"), overwrite = TRUE)
  
  message("Finished: ", basename(f))
}
