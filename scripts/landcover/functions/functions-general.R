# Functions to write spatial data to disk with helpful tags (e.g. cropped)
# Generic tagging
collapse_tags <- function(...) {
  tags <- c(...)
  tags <- tags[nzchar(tags)]
  if (length(tags) == 0) NULL else paste(tags, collapse = "_")
}

# tagged_path <- function(root, rel_dir, filename, ...) {
#   stem <- tools::file_path_sans_ext(filename)
#   ext  <- tools::file_ext(filename)
#   tag  <- collapse_tags(...)
#   name <- if (is.null(tag)) stem else sprintf("%s_%s", stem, tag)
#   if (identical(ext, "")) ext <- "tif"
#   file.path(root, rel_dir, sprintf("%s.%s", name, ext))
# }

tagged_path <- function(root, rel_dir, filename, ..., driver = "ESRI Shapefile") {
  stem <- tools::file_path_sans_ext(filename)
  ext  <- tools::file_ext(filename)
  tag  <- collapse_tags(...)
  name <- if (is.null(tag)) stem else sprintf("%s_%s", stem, tag)
  
  if (driver == "GPKG") {
    return(file.path(root, rel_dir, sprintf("%s.gpkg", name)))
  }
  
  if (identical(ext, "")) ext <- "tif"
  file.path(root, rel_dir, sprintf("%s.%s", name, ext))
}

# Write vector data
# write_sf_tagged <- function(x, root, rel_dir, filename, ..., quiet = TRUE) {
#   stopifnot(inherits(x, "sf"))
#   stem <- tools::file_path_sans_ext(filename)
#   tag  <- collapse_tags(...)
#   layer <- if (is.null(tag)) stem else sprintf("%s_%s", stem, tag)
#   dsn <- file.path(root, rel_dir)
#   dir.create(dsn, recursive = TRUE, showWarnings = FALSE)
#   sf::st_write(
#     x,
#     dsn = dsn,
#     layer = layer,
#     driver = "ESRI Shapefile",
#     delete_layer = TRUE,
#     quiet = quiet
#   )
# }
write_sf_tagged <- function(x, root, rel_dir, filename, ..., 
                            driver = "ESRI Shapefile", quiet = TRUE) {
  stopifnot(inherits(x, "sf"))
  stem <- tools::file_path_sans_ext(filename)
  tag  <- collapse_tags(...)
  layer <- if (is.null(tag)) stem else sprintf("%s_%s", stem, tag)
  
  ext <- switch(driver,
                "GPKG"           = ".gpkg",
                "ESRI Shapefile" = ".shp",
                ".gpkg"
  )
  
  dsn <- switch(driver,
                "GPKG"           = file.path(root, rel_dir, paste0(layer, ext)),
                "ESRI Shapefile" = file.path(root, rel_dir)
  )
  
  dir.create(
    ifelse(driver == "GPKG", 
           file.path(root, rel_dir), 
           dsn),
    recursive = TRUE, showWarnings = FALSE
  )
  
  sf::st_write(
    x,
    dsn    = dsn,
    layer  = layer,
    driver = driver,
    delete_layer = TRUE,
    quiet  = quiet
  )
}

# Write raster data
write_rast_tagged <- function(r, root, rel_dir, filename, ..., overwrite = TRUE,
                              datatype = "FLT4S",
                              gdal = c("COMPRESS=LZW", "TILED=YES")) {
  stopifnot(inherits(r, "SpatRaster"))
  out <- tagged_path(root, rel_dir, filename, ...)
  dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
  
  terra::writeRaster(
    r, out,
    overwrite = overwrite,
    wopt = list(datatype = datatype, gdal = gdal)
  )
  
  invisible(out)
}

#Capitalize first letter of a string
capitalize_first <- function(x) {
  paste0(
    toupper(substr(x, 1, 1)),
    substr(x, 2, nchar(x))
  )
}
