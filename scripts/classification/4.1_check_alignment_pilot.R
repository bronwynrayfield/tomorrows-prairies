
###################################################################################
##            VALIDATION: CHECK ALIGNMENT OF RASTERS IN PILOT AREA               ##
###################################################################################

##----------------------- Alison Long - February 2026 -------------------------

# step 4: check alignment of rasters cropped to pilot area in step 3
#         to prep for pilot analysis

# load packages
library(terra)
library(sf)
library(tidyverse)

# Set the temporary directory to a drive with sufficient space
terraOptions(tempdir = "R:/CanadianGrasslands/WG_classification/scratch")

####--- DIRECTORIES ---####

# pilot area
pa_dir <- "R:/CanadianGrasslands/processeddatasets/pilotarea"

# processed data directory
proc_data_dir <- "R:/CanadianGrasslands/WG_Classification/processed"

# pilot area data directories
aci_pilot_dir <- file.path(proc_data_dir, "aafc_aci/pilot")
luts_pilot_dir <- file.path(proc_data_dir, "aafc_luts_2020_30m_v5/pilot")
cec_pilot_dir <- file.path(proc_data_dir, "cec_nalcms/pilot")
ngi_pilot_dir <- file.path(proc_data_dir, "cfga_ngi/pilot")
cgr_pilot_dir <- file.path(proc_data_dir, "cgr_gam_adv/pilot")
eccc_pilot_dir <- file.path(proc_data_dir, "eccc_grass/pilot")
esa_pilot_dir <- file.path(proc_data_dir, "esa_worldcover/pilot")
glc_pilot_dir <- file.path(proc_data_dir, "glc_fcs30d/pilot")
cnvc_pilot_dir <- file.path(proc_data_dir, "natureserve_cnvc/pilot")
lcc_pilot_dir <- file.path(proc_data_dir, "nrcan_lcc/pilot")
pudl_pilot_dir <- file.path(proc_data_dir, "phjv_pudl/pilot")
lcr_pilot_dir <- file.path(proc_data_dir, "statcan_lcr/pilot")
gpw_pilot_dir <- file.path(proc_data_dir, "wri_gpw_ann/pilot")
plow_pilot_dir <- file.path(proc_data_dir, "wwf_plowprint/pilot")

####--- FILES ---####

# Pilot area
pa_file <- file.path(pa_dir, "GNP_pilot_area_crop_CSD.shp")

# Study area
#sa_file <- file.path(sa_dir, "class_study_area_aligned_cnvc_Dec2025.tif")

# rasters
aci_file <- file.path(aci_pilot_dir, "aci_merged_reproj_algn_crop_pa_binary_noshrub.tif")
luts_file <- file.path(luts_pilot_dir, "luts_reproj_algn_crop_pa_binary.tif")
cec_file <- file.path(cec_pilot_dir, "nalcms_reproj_algn_crop_pa_binary_noshrub.tif")
ngi_file <- file.path(ngi_pilot_dir, "ngi_merged_reproj_algn_crop_pa_binary.tif")
cgr_file <- file.path(cgr_pilot_dir, "cgr_merged_reproj_algn_crop_pa_binary.tif")
eccc_file <- file.path(eccc_pilot_dir, "eccc_reproj_algn_crop_pa_binary.tif")
cnvc_file <- file.path(cnvc_pilot_dir, "cnvc_crop_pa_binary.tif")
lcc_file <- file.path(lcc_pilot_dir, "lcc_reproj_algn_crop_pa_binary_noshrub.tif")
pudl_file <- file.path(pudl_pilot_dir, "pudl_reproj_algn_crop_pa_binary.tif")
lcr_file <- file.path(lcr_pilot_dir, "lcr_reproj_algn_crop_pa_binary.tif")
gpw_file <- file.path(gpw_pilot_dir, "gpw_reproj_algn_crop_pa_binary_noshrub.tif")
plow_file <- file.path(plow_pilot_dir, "plow_reproj_algn_crop_pa_binary.tif")

####--- LOAD RASTERS ---####
aci_pa <- rast(aci_file)
luts_pa <- rast(luts_file)
cec_pa <- rast(cec_file)
ngi_pa <- rast(ngi_file)
cgr_pa <- rast(cgr_file)
eccc_pa <- rast(eccc_file)
cnvc_pa <- rast(cnvc_file)
lcc_pa <- rast(lcc_file)
pudl_pa <- rast(pudl_file)
lcr_pa <- rast(lcr_file)
gpw_pa <- rast(gpw_file)
plow_pa <- rast(plow_file)

# create a list of rasters
rlist <- list(
  aci   = aci_pa,
  luts  = luts_pa,
  cec   = cec_pa,
  ngi   = ngi_pa,
  cgr   = cgr_pa,
  eccc  = eccc_pa,
  cnvc  = cnvc_pa,
  lcc   = lcc_pa,
  pudl  = pudl_pa,
  lcr   = lcr_pa,
  gpw   = gpw_pa,
  plow  = plow_pa
)

####--- CHECK ALIGNMENT ---####

# Determine if a SpatRaster layer is categorical (factor) using levels()
is_categorical <- function(r) {
  # TRUE if any layer has a non-null levels table
  any(!vapply(levels(r), is.null, logical(1)))
}

# Summarize a raster
summarize_raster <- function(r, name = NULL) {
  name <- if (is.null(name)) deparse(substitute(r)) else name
  has_lvls <- any(!vapply(levels(r), is.null, logical(1)))
  list(
    name       = name,
    crs        = as.character(crs(r)),
    nlyr       = nlyr(r),
    dim        = paste(nrow(r), "rows x", ncol(r), "cols"),
    res        = res(r),
    origin     = origin(r),
    extent     = as.vector(ext(r)),   # xmin, xmax, ymin, ymax
    ncell      = ncell(r),
    has_levels = has_lvls,
    is_factor  = has_lvls,            # alias for readability
    datatype   = terra::datatype(r),
    NAflag     = tryCatch(NAflag(r), error = function(e) NA) # may be NA/NULL depending on file
  )
}

# Pairwise geometry comparison using terra::compareGeom
pairwise_compare <- function(rasters, names = NULL) {
  n <- length(rasters)
  if (n < 2) return(invisible(NULL))
  if (is.null(names)) names <- paste0("r", seq_len(n))
  out <- list()
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      same_crs  <- (crs(rasters[[i]]) == crs(rasters[[j]]))
      same_res  <- all(res(rasters[[i]]) == res(rasters[[j]]))
      same_org  <- all(origin(rasters[[i]]) == origin(rasters[[j]]))
      cg        <- try(compareGeom(
        rasters[[i]], rasters[[j]],
        crs = TRUE, rowcol = TRUE, res = TRUE,
        stopOnError = FALSE),
        silent = TRUE)
      out[[paste(names[i], names[j], sep = " vs ")]] <- list(
        same_crs    = same_crs,
        same_res    = same_res,
        same_origin = same_org,
        compareGeom = cg
      )
    }
  }
  out
}

# Full alignment report
check_alignment <- function(rlist, rnames = NULL, print_report = TRUE) {
  if (is.null(rnames)) rnames <- names(rlist)
  if (is.null(rnames) || any(!nzchar(rnames))) {
    rnames <- paste0("r", seq_along(rlist))
  }
  
  per_raster <- lapply(seq_along(rlist), function(i) summarize_raster(rlist[[i]], rnames[i]))
  pw <- pairwise_compare(rlist, rnames)
  
  if (print_report) {
    cat("\n=== Per-raster summary ===\n")
    for (i in seq_along(per_raster)) {
      s <- per_raster[[i]]
      cat(sprintf(
        "\n[%s]\n CRS: %s\n Layers: %d\n Dim: %s\n Res: %s\n Origin: %s\n Extent: c(%.6f, %.6f, %.6f, %.6f)\n Cells: %d\n Levels? %s  Factor? %s\n Data type: %s\n NAflag: %s\n",
        s$name, s$crs, s$nlyr, s$dim, paste(s$res, collapse=", "),
        paste(s$origin, collapse=", "), s$extent[1], s$extent[2], s$extent[3], s$extent[4],
        s$ncell, s$has_levels, s$is_factor, s$datatype, as.character(s$NAflag)
      ))
    }
    if (length(pw)) {
      cat("\n=== Pairwise alignment (compareGeom) ===\n")
      for (nm in names(pw)) {
        entry <- pw[[nm]]
        cat(sprintf("\n%s:\n  same_crs: %s\n  same_res: %s\n  same_origin: %s\n",
                    nm, entry$same_crs, entry$same_res, entry$same_origin))
        cat("  compareGeom result: ")
        if (inherits(entry$compareGeom, "try-error")) {
          cat("error running compareGeom\n")
        } else {
          print(entry$compareGeom)
        }
      }
    }
  }
  
  invisible(list(per_raster = per_raster, pairwise = pw))
}


# Example usage:
rep <- check_alignment(rlist)
