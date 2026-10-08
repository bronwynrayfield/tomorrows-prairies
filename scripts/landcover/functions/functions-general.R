# Functions to write spatial data to disk with helpful tags (e.g. cropped)


# ---- Generic tagging ----
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


# ---- LULC change ----
# Ecodistrict -> ecoregion lookup
eco_lut <- function(ecoregions_sa, ecodistricts_sa) {
  eco_names <- as.data.frame(ecoregions_sa) |>
    distinct(ecoregion = ECOREGION_, ecoregion_name = ECOREGION1)
  as.data.frame(ecodistricts_sa) |>
    distinct(stratum_id = as.numeric(ECODISTRIC), ecoregion = ECOREGION_) |>
    left_join(eco_names, by = "ecoregion")
}

# ---- area of every (ecodistrict, from, to) combination ----
# packs three values per pixel (ecodistrict, 2010 class, 2020 class) into one integer, counts how often each integer occurs, then unpacks them
tabulate_transitions <- function(r_from, r_to, r_sid, codes = cls) {
  if (is.list(codes)) codes <- setNames(unlist(codes, use.names = FALSE), rep(names(codes), lengths(codes)))
  k <- 10^ceiling(log10(max(codes, minmax(r_from), minmax(r_to), na.rm = TRUE) + 1))   # 10 for single-digit codes; 100 for LUTS (max 16)
  r_sid <- ifel(is.na(r_sid), 0, r_sid)    # pixels outside ecodistricts kept as stratum 0, not dropped
  # combined code stays exact in FLT4S while < 16,777,216: k = 10 -> sid < 167,772; k = 100 -> sid < 1,677
  f  <- freq((r_sid * k + r_from) * k + r_to, bylayer = FALSE)
  v  <- f$value
  fc <- (v %/% k) %% k
  tc <- v %% k
  data.frame(
    stratum_id = v %/% k^2,
    from_code  = fc, to_code = tc,
    from       = names(codes)[match(fc, codes)],   # NA here = code not in lookup (caught by coverage_qc)
    to         = names(codes)[match(tc, codes)],
    pixels     = f$count,
    area_ha    = f$count * prod(res(r_from)) / 1e4
  )
}

# LULC builder
# Order matters: tame is tested before annual, so `annual` may be the whole
# LUTS ag domain (2010) and resolves to ag minus tame.
# LUTS binaries are mutually exclusive.
build_lulc <- function(native, tame, annual, settlement) {
  r <- ifel(native == 1, cls_native,
            ifel(tame, cls_tame,
                 ifel(annual, cls_annual,
                      ifel(settlement == 1, cls_settlement, cls_other))))
  levels(r) <- lulc_levels
  r
}

# Scoring for one transition: rows for the origin class only;
# 1 = focal, 0 = persistence, NA = competing. classify() sends all other codes to NA.
score_response <- function(from_cls, to_cls) {
  d <- data.frame(from = from_cls, to = lulc_levels[[1]])
  d$score <- ifelse(d$to == to_cls, 1, ifelse(d$to == from_cls, 0, NA))
  d
}

# One pass: pixel counts by stratum x from x to (classes 1-5)
tabulate_lulc <- function(l10, l20, strata) {
  f <- freq(strata * 100 + l10 * 10 + l20)
  data.frame(stratum_id = f$value %/% 100, from = f$value %/% 10 %% 10,
             to = f$value %% 10, n = f$count)
}

# Counting ACI years as evidence of tame
count_name <- function(code, yrs) sprintf("aci_%d_count_%s_%dm", code,
                                          paste(range(yrs), collapse = "-"), res_m)

# Apply tame evidence threshold
tame_ev <- function(counts, codes, min)
  subst(sum(counts[[paste0("aci", strsplit(codes, "\\+")[[1]])]]), NA, 0) >= min

# Summarises a tabulated, scored transition table (grouped beforehand):
# n_pos / n_neg = pixels scored 1 / 0; n_competing = pixels set to NA from the
# same origin; pos_ha = positives in ha; rate = n_pos / (n_pos + n_neg);
# competing_per_pos = n_competing / n_pos.
sum_scores <- function(d) summarise(d,
                                    n_pos       = sum(n[score %in% 1]),
                                    n_neg       = sum(n[score %in% 0]),
                                    n_competing = sum(n[is.na(score)]),
                                    pos_ha      = n_pos * cell_ha,
                                    rate        = n_pos / (n_pos + n_neg),
                                    competing_per_pos = n_competing / n_pos,
                                    .groups = "drop")
