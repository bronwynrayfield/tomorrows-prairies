#' Reclassify AAFC LUTS raster to named groups
#'
#' @param lulc        A terra SpatRaster — raw AAFC LUTS raster (factor or numeric)
#' @param scheme      One of "4class", "5class", or "binary"
#' @param binary_class For scheme = "binary": which group to code as 1.
#'                    One of "native", "crop", "tame", "settlement", "other"
#' @param check       If TRUE, prints a frequency table with group labels for QA/QC
#'
#' @return A reclassified SpatRaster
#'
#' @note RAT differences between years:
#'   2010 only:   85  = Newly-Detected Road <10 years         -> other
#'   2015 & 2020: 52  = Annual Cropland                       -> crop
#'                56  = Land Converted to Annual Cropland     -> crop
#'                92  = Snow and Ice                          -> other
#'   All other codes are shared between 2010 and 2015.
#'
#' @examples
#' reclass_aafc(lulc2010, scheme = "4class", check = TRUE)
#' reclass_aafc(lulc2010, scheme = "binary", binary_class = "native", check = TRUE)
#' reclass_aafc(lulc2020, scheme = "binary", binary_class = "crop",   check = TRUE)

reclass_aafc <- function(lulc, scheme = "5class", binary_class = NULL, check = FALSE) {
  
  # ---- Define groups ----
  # native     : Prairie/Temperate/Boreal/Montane/Alpine grassland (63-68)
  # tame       : Cropland (51), Land Converted to Cropland (55)
  # crop       : Annual Cropland (52), Land Converted to Annual Cropland (56)
  #              NOTE: annual cropland classes only present in 2020
  # settlement : Settlement classes (21,22,24,28,29,81,82,84,88,89)
  # forest     : Forest (41,42,43,44,47,48,49)
  # water      : Water (31)
  # other      : Roads (25),
  #              Wetland (73,74,75,76), Other Land (91), Snow & Ice (92),
  #              Newly-Detected Road <10 years (85) [2010 only]
  # ---- Define groups ----
  groups <- list(
    native     = c(63, 64, 65, 66, 67, 68),
    tame       = c(51, 55),
    crop       = c(52, 56),
    settlement = c(21, 22, 24, 28, 29, 81, 82, 84, 88, 89),
    forest     = c(41, 42, 43, 44, 47, 48, 49),
    water      = c(31),
    other      = c(25, 73, 74, 75, 76, 85, 91, 92)
  )
  
  # ---- Define output codes and labels per scheme ----
  if (scheme == "5class") {
    # Collapse forest, water, other into single "other" group for this scheme
    groups$other <- c(groups$forest, groups$water, groups$other)
    groups$forest <- NULL
    groups$water  <- NULL
    to     <- c(native = 1L, tame = 2L, crop = 3L, settlement = 4L, other = 5L)
    labels <- c("1" = "Native", "2" = "Tame", "3" = "Crop",
                "4" = "Settlement", "5" = "Other")
    
  } else if (scheme == "4class") {
    groups$crop  <- c(groups$crop, groups$tame)
    groups$tame  <- NULL
    groups$other <- c(groups$forest, groups$water, groups$other)
    groups$forest <- NULL
    groups$water  <- NULL
    to     <- c(native = 1L, crop = 3L, settlement = 4L, other = 5L)
    labels <- c("1" = "Native", "3" = "Crop (incl. tame)",
                "4" = "Settlement", "5" = "Other")
    
  } else if (scheme == "binary") {
    if (is.null(binary_class) || !binary_class %in% names(groups)) {
      stop("For scheme = 'binary', supply binary_class: one of ",
           paste(names(groups), collapse = ", "))
    }
    to     <- setNames(ifelse(names(groups) == binary_class, 1L, 0L), names(groups))
    labels <- c("0" = paste0("Not ", binary_class), "1" = binary_class)
    
  } else {
    stop("scheme must be one of: '4class', '5class', 'binary'")
  }
  
  # ---- Build from/to matrix ----
  rcl <- data.frame(
    from  = unlist(groups, use.names = FALSE),
    group = rep(names(groups), lengths(groups))
  )
  rcl$to  <- unname(to[rcl$group])
  rcl_mat <- as.matrix(rcl[, c("from", "to")])
  
  # ---- Handle factor raster ----
  if (is.factor(lulc)) {
    lut         <- levels(lulc)[[1]]
    id_col      <- names(lut)[1]
    lulc        <- as.numeric(lulc)
    idx_to_code <- setNames(lut[[id_col]], seq_len(nrow(lut)))
    lulc        <- classify(lulc, cbind(as.integer(names(idx_to_code)),
                                        as.integer(idx_to_code)))
  }
  
  # ---- Validate all raster codes are covered ----
  raster_codes <- unique(values(lulc, na.rm = TRUE))
  all_defined  <- unlist(groups, use.names = FALSE)
  uncovered    <- setdiff(raster_codes, all_defined)
  if (length(uncovered) > 0) {
    warning("The following codes are not covered by any group and will be NA: ",
            paste(uncovered, collapse = ", "))
  }
  
  # ---- Reclassify ----
  out <- terra::classify(lulc, rcl_mat, others = NA)
  
  # ---- QA/QC check ----
  if (check) {
    f       <- as.data.frame(freq(out))
    f$label <- labels[as.character(f$value)]
    f       <- f[, c("value", "label", "count")]
    f$pct   <- round(100 * f$count / sum(f$count, na.rm = TRUE), 2)
    cat("\n--- reclass_aafc QA/QC: scheme =", scheme)
    if (scheme == "binary") cat(", binary_class =", binary_class)
    cat(" ---\n")
    print(f, row.names = FALSE)
    cat("Total cells:", sum(f$count, na.rm = TRUE), "\n\n")
  }
  
  return(out)
}