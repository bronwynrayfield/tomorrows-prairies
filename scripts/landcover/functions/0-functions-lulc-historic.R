# suppressPackageStartupMessages({
#   library(terra)
#   library(dplyr)
#   library(sf)
#   library(ggplot2)
# })
# 
# # ---- lookup (same as before) ----
# .default_class_lookup <- function() {
#   c(
#     `1`  = "forest",
#     `2`  = "shrub",
#     `3`  = "crop",
#     `4`  = "upast",
#     `5`  = "grass",
#     `6`  = "builtup",
#     `7`  = "bare",
#     `8`  = "wetland",
#     `9`  = "water",
#     `10` = "ipast",
#     `11` = "forest/shrub",
#     `12` = "forest/wetland",
#     `13` = "shrub/forest",
#     `14` = "wetland/forest",
#     `15` = "upast/grass",
#     `16` = "grass/upast"
#   )
# }
# 
# # ---- full-run summariser (tiling + points/table) ----
# # NOTE: this is the version that worked conceptually for you:
# # - rasterize polygon IDs to raster grid using a temp field
# # - build "to_recent" = destination code only for years in window
# # - loop over extent in tiles; in each tile: as.points(na.rm=TRUE) + table(stratum_id,to_code)
# summarize_recent_loss_by_dest_full <- function(
#     r_class,
#     r_year,
#     ecodistricts,      # SpatVector (preferred) or sf
#     eco_id_field,
#     year_min,
#     year_max,
#     tile_size = 200000,
#     tempdir = "D:/terra-temp",
#     class_lookup = NULL
# ) {
#   if (is.null(class_lookup)) class_lookup <- .default_class_lookup()
#   
#   stopifnot(inherits(r_class, "SpatRaster"), inherits(r_year, "SpatRaster"))
#   if (!terra::compareGeom(r_class, r_year, stopOnError = FALSE)) {
#     stop("class and year rasters do not share the same geometry (CRS/res/extent).")
#   }
#   
#   dir.create(tempdir, recursive = TRUE, showWarnings = FALSE)
#   terra::terraOptions(todisk = TRUE, memfrac = 0.6, tempdir = tempdir)
#   
#   v_eco <- if (inherits(ecodistricts, "SpatVector")) ecodistricts else terra::vect(ecodistricts)
#   if (!(eco_id_field %in% names(v_eco))) stop("eco_id_field not found in ecodistricts attributes.")
#   
#   # pixel area
#   resm <- terra::res(r_class)
#   pix_area_ha <- (resm[1] * resm[2]) / 10000
#   
#   # year window is encoded as (year - 1900)
#   lo <- year_min - 1900
#   hi <- year_max - 1900
#   
#   # destination code
#   to_code <- r_class %% 100
#   
#   # to_recent raster on disk
#   to_recent <- terra::ifel(
#     (r_year >= lo) & (r_year <= hi),
#     to_code,
#     NA,
#     filename  = file.path(tempdir, sprintf("to_recent_%d_%d_%s.tif", year_min, year_max, sample(1e9,1))),
#     overwrite = TRUE
#   )
#   
#   # rasterize ecodistrict IDs (force correct attribute via temp field)
#   v_eco$.eco_id_tmp <- v_eco[[eco_id_field]]
#   r_sid <- terra::rasterize(
#     v_eco, r_class, field = ".eco_id_tmp",
#     filename  = file.path(tempdir, sprintf("eco_id_%s.tif", sample(1e9,1))),
#     overwrite = TRUE
#   )
#   
#   # tiling grid
#   e <- terra::ext(r_class)
#   xs <- seq(e$xmin, e$xmax, by = tile_size)
#   ys <- seq(e$ymin, e$ymax, by = tile_size)
#   
#   acc <- new.env(parent = emptyenv())
#   `%||%` <- function(a, b) if (is.null(a)) b else a
#   
#   add_counts <- function(stratum_ids, to_codes) {
#     key <- paste(stratum_ids, to_codes, sep = "_")
#     tab <- table(key)
#     nms <- names(tab)
#     for (k in seq_along(tab)) {
#       nm <- nms[k]
#       acc[[nm]] <- (acc[[nm]] %||% 0L) + as.integer(tab[[k]])
#     }
#   }
#   
#   for (i in seq_along(xs)) {
#     x0 <- xs[i]
#     x1 <- min(x0 + tile_size, e$xmax)
#     
#     for (j in seq_along(ys)) {
#       y0 <- ys[j]
#       y1 <- min(y0 + tile_size, e$ymax)
#       
#       et <- terra::ext(x0, x1, y0, y1)
#       
#       sid_t <- terra::crop(r_sid, et)
#       rec_t <- terra::crop(to_recent, et)
#       
#       n_recent <- terra::global(!is.na(rec_t), "sum", na.rm = TRUE)[1, 1]
#       if (is.na(n_recent) || n_recent == 0) next
#       
#       s <- terra::mask(c(sid_t, rec_t), rec_t)
#       names(s) <- c("stratum_id", "to_code")
#       
#       pts <- terra::as.data.frame(terra::as.points(s, na.rm = TRUE))
#       if (!nrow(pts)) next
#       
#       add_counts(pts$stratum_id, pts$to_code)
#     }
#   }
#   
#   keys <- ls(acc)
#   if (!length(keys)) return(data.frame())
#   
#   counts <- vapply(keys, function(k) acc[[k]], integer(1))
#   parts <- do.call(rbind, strsplit(keys, "_", fixed = TRUE))
#   
#   out <- data.frame(
#     stratum_id = as.integer(parts[, 1]),
#     to_code    = as.integer(parts[, 2]),
#     pixels     = as.numeric(counts)
#   )
#   
#   out$area_ha <- out$pixels * pix_area_ha
#   out$to_class <- unname(class_lookup[as.character(out$to_code)])
#   out$to_class[is.na(out$to_class)] <- paste0("code_", out$to_code)
#   
#   out[order(out$stratum_id, out$to_code), ]
# }
# 
# # ---- QA/QC helpers ----
# qc_loss_table <- function(df, name = "loss_df") {
#   stopifnot(all(c("stratum_id", "to_code", "area_ha") %in% names(df)))
#   
#   message("\n[QC] ", name)
#   message("  rows: ", nrow(df))
#   message("  unique ecodistricts: ", n_distinct(df$stratum_id))
#   message("  total area (ha): ", round(sum(df$area_ha, na.rm = TRUE), 3))
#   
#   message("  top destinations by area (ha):")
#   print(
#     df %>%
#       group_by(to_code, to_class) %>%
#       summarise(area_ha = sum(area_ha, na.rm = TRUE), .groups = "drop") %>%
#       arrange(desc(area_ha)) %>%
#       head(10)
#   )
# }
# 
# summarise_metrics <- function(df, prefix) {
#   # prefix chosen to keep shapefile field names <=10 chars
#   df %>%
#     group_by(stratum_id) %>%
#     summarise(
#       "{prefix}_bu" := sum(area_ha[to_code == 6], na.rm = TRUE),
#       "{prefix}_cr" := sum(area_ha[to_code == 3], na.rm = TRUE),
#       "{prefix}_bc" := sum(area_ha[to_code %in% c(3,6)], na.rm = TRUE),
#       .groups = "drop"
#     ) %>%
#     as.data.frame()
# }
# 
# qc_gpu_equals_parts <- function(m_gpu, m_g, m_pu) {
#   x <- m_gpu %>%
#     full_join(m_g,  by = "stratum_id") %>%
#     full_join(m_pu, by = "stratum_id")
#   
#   for (nm in setdiff(names(x), "stratum_id")) x[[nm]][is.na(x[[nm]])] <- 0
#   
#   dif_bu <- max(abs(x$gpu_bu - (x$g_bu + x$pu_bu)))
#   dif_cr <- max(abs(x$gpu_cr - (x$g_cr + x$pu_cr)))
#   dif_bc <- max(abs(x$gpu_bc - (x$g_bc + x$pu_bc)))
#   
#   message("\n[QC] gpu == g + pu max abs diffs (ha):")
#   message("  built-up: ", dif_bu)
#   message("  crop:     ", dif_cr)
#   message("  bu+crop:  ", dif_bc)
# }
# 
# 
# # map_metric_linear_grey0 <- function(eco_sf, field, title, subtitle = "Grey = zero", out_png = NULL) {
# #   stopifnot(field %in% names(eco_sf))
# #   
# #   tmp <- eco_sf
# #   tmp[[paste0(field, "_plot")]] <- ifelse(tmp[[field]] > 0, tmp[[field]], NA)
# #   
# #   p <- ggplot(tmp) +
# #     geom_sf(aes(fill = .data[[paste0(field, "_plot")]]), color = NA) +
# #     scale_fill_viridis_c(
# #       option = "C",
# #       na.value = "grey85",
# #       name = paste0(field, " (ha)")
# #     ) +
# #     labs(title = title, subtitle = subtitle) +
# #     theme_minimal()
# #   
# #   if (!is.null(out_png)) {
# #     ggsave(out_png, p, width = 9, height = 6, dpi = 300)
# #   }
# #   
# #   p
# # }
# # 
# # 
# # map_metric_log1p_grey0 <- function(eco_sf, field, title, subtitle = "Grey = zero", out_png = NULL) {
# #   stopifnot(field %in% names(eco_sf))
# #   tmp <- eco_sf
# #   tmp[[paste0(field, "_log1p")]] <- ifelse(tmp[[field]] > 0, log1p(tmp[[field]]), NA)
# #   
# #   p <- ggplot(tmp) +
# #     geom_sf(aes(fill = .data[[paste0(field, "_log1p")]]), color = NA) +
# #     scale_fill_viridis_c(option = "C", na.value = "grey85", name = paste0(field, "\nlog1p(ha)")) +
# #     labs(title = title, subtitle = subtitle) +
# #     theme_minimal()
# #   
# #   if (!is.null(out_png)) ggsave(out_png, p, width = 9, height = 6, dpi = 300)
# #   p
# # }
# 
# 
# map_metric_grey0 <- function(
#     eco_sf,
#     field,
#     title,
#     scale = c("linear", "log1p"),
#     subtitle = "Grey = zero",
#     prov = NULL,
#     cities = NULL,
#     city_labels = NULL,
#     out_png = NULL
# ) {
#   scale <- match.arg(scale)
#   stopifnot(field %in% names(eco_sf))
#   
#   tmp <- eco_sf
#   
#   plot_field <- if (scale == "linear") {
#     tmp[[paste0(field, "_plot")]] <- ifelse(tmp[[field]] > 0, tmp[[field]], NA)
#     paste0(field, "_plot")
#   } else {
#     tmp[[paste0(field, "_plot")]] <- ifelse(tmp[[field]] > 0, log1p(tmp[[field]]), NA)
#     paste0(field, "_plot")
#   }
#   
#   fill_name <- if (scale == "linear") {
#     paste0(field, " (ha)")
#   } else {
#     paste0(field, "\nlog1p(ha)")
#   }
#   
#   p <- ggplot(tmp) +
#     geom_sf(aes(fill = .data[[plot_field]]), color = NA)
#   
#   # province boundaries
#   if (!is.null(prov)) {
#     p <- p + geom_sf(data = prov, fill = NA, linewidth = 0.3)
#   }
#   
#   # city points
#   if (!is.null(cities)) {
#     p <- p + geom_sf(data = cities, color = "black", size = 1)
#   }
#   
#   # city labels (Fix B: deterministic)
#   if (!is.null(city_labels)) {
#     p <- p + geom_sf_text(
#       data = city_labels,
#       aes(label = NAME),
#       size = 3,
#       check_overlap = FALSE
#     )
#   }
#   
#   p <- p +
#     scale_fill_viridis_c(
#       option = "C",
#       na.value = "grey85",
#       name = fill_name
#     ) +
#     labs(title = title, subtitle = subtitle) +
#     coord_sf(
#       xlim = c(bb["xmin"], bb["xmax"]),
#       ylim = c(bb["ymin"], bb["ymax"]),
#       expand = FALSE
#     ) +
#     theme_minimal() +
#     theme(axis.title = element_blank())
#   
#   if (!is.null(out_png)) {
#     dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)
#     ggplot2::ggsave(filename = out_png, plot = p, width = 9, height = 6, dpi = 300)
#     message("Saved: ", normalizePath(out_png, winslash = "/", mustWork = FALSE))
#   }
#   
#   p
# }
# 



# Class codes from cls_* in 0a-setup.R; the names become the labels in every output.
cls <- c(other = cls_other, native = cls_native, tame = cls_tame,
         annual = cls_annual, settlement = cls_settlement)

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


