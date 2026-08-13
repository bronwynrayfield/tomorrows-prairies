# =============================================================================
# 7-response-plowprint.R
# =============================================================================
#
# Purpose:
#   Builds a Plowprint-style benchmark response for comparison against
#   Olimb & Robinson (2019). This is a PIPELINE CHECK, not a model that feeds
#   CLUES: it tests whether the predictor stack, extraction and spatial CV
#   machinery reproduce the kind of result they obtained, before those tools
#   are trusted on the harder conversion-risk target.
#
#   Response:  1 = previously or currently cultivated (Plowprint 1/2/3)
#              0 = intact habitat (Plowprint 0, restricted -- see below)
#              NA = outside the training set
#
#   This is NOT a change response. It has no 2010 baseline, no grass_origin,
#   and no conversion interval. Because Plowprint is cumulative
#   ("ever cultivated"), it is close to a soil-and-climate suitability surface
#   by construction, which is why biophysical predictors do well on it. Expect
#   a HIGH AUC and read it as "the pipeline works", not as evidence about the
#   conversion model. It says nothing about conversion timing or the marginal
#   pixel.
#
# Inputs:
#   - WWF Plowprint 2024 (aligned to the modelling extent by an earlier script)
#   - AAFC LUTS 2020 (cropped to modelling extent)
#   - Study area raster (sa_r) : Prairie ecozone
#
# Outputs:
#   - response_plowprint_30m.tif   : 1 / 0 / NA
#   - response_plowprint_300m.tif  : near-sampled, for quick diagnostics
#   - plowprint_response_summary.csv
#
# -----------------------------------------------------------------------------
# CLASS SCHEME (determined empirically -- the raster carries no RAT or coltab)
# -----------------------------------------------------------------------------
#   Verified by crosstab against LUTS 2020 within the study area at 300 m:
#     0        not plowed. Annual Cropland is only 0.3% of this class, so it
#              essentially never misses a real crop pixel.
#     1, 2     small classes (~1.4M px at 30 m combined), 61-65% of each sits
#              in LUTS perennial Cropland. Possibly conversion-period bins or a
#              hay-specific flag; not documented. Folded into the positive
#              class -- too small to move the model either way.
#     3        main plowed class. 97% of LUTS Annual Cropland is class 3.
#     40, 60   developed / roads.            } excluded: non-habitat, and
#     100      water and roads mask.         } calling them intact is indefensible
#
# -----------------------------------------------------------------------------
# TRAINING SET DEFINITION
# -----------------------------------------------------------------------------
#   Plowprint 0 means "not plowed", NOT "intact grassland". Within the Prairies
#   Ecozone it breaks down as:
#     native grassland     810,729  52.0%   keep
#     perennial cropland   345,117  22.2%   keep -- see note below
#     forest               217,164  13.9%   DROP
#     wetland              131,652   8.5%   keep
#     settlement + roads    24,945   1.6%   DROP
#     water                 18,219   1.2%   DROP
#     Other Land             5,883   0.4%   DROP
#     annual cropland        3,984   0.3%   (negligible)
#
#   Forest, water, settlement, roads and Other Land are dropped: they are not
#   candidates for cropping, and leaving them in would let biophysical
#   predictors separate forest climate from cropland climate -- an easy problem
#   that inflates AUC without testing anything. O&R almost certainly made these
#   exclusions even though the paper does not list them.
#
#   PERENNIAL CROPLAND IS KEPT AS A NEGATIVE, with reservations. Plowprint
#   codes it 0 (not plowed), but it is tame forage -- cultivated land under
#   Plowprint's own "previously or currently cultivated" definition. The
#   product is internally inconsistent here, presumably because ACI never
#   flagged those pixels as an annual crop class in any year. This is the same
#   tame-grass ambiguity that drove the ACI work for the main response, showing
#   up in WWF's product rather than ours. It is left as published rather than
#   hand-edited, since the point of this run is comparability.
#
#   Consequence: Plowprint's positive class and this project's tame stratum are
#   NOT the same population. Do not read across between the two models.
#
#   Wetland is kept as intact habitat, which is defensible but INCONSISTENT
#   with the main response, where wetland conversion is out of scope entirely.
#   Acceptable for a pipeline check; note it in the write-up.
#
# Notes:
#   - Plowprint is categorical: method = "near" everywhere, never modal.
#   - Vintage: the 2024 release. Olimb & Robinson used the 2016 version.
#     Plowprint is cumulative, so 2024 = 2016 plus eight more years of
#     accumulation; the positive class here is larger than theirs. The spatial
#     pattern is dominated by the same historical footprint, so this is
#     acceptable for a pipeline check -- but state the vintage difference
#     rather than implying replication.
#   - WWF paused the 2025 Plowprint analysis because of changes in the
#     underlying government data and is reviewing alternative datasets, so the
#     2024 release is currently the last one.
#   - Study area differs too: O&R worked the US Great Plains, this is the
#     Canadian Prairies Ecozone. The question is "does my pipeline reproduce
#     their kind of result", not "do I get their numbers".
#
# Dependencies:
#   - 0-setup.R (parameters, file paths, functions)
# =============================================================================

source(file.path("scripts/landcover", "0-setup.R"))

# ---- Parameters ----
pp_plowed   <- c(1, 2, 3)
pp_intact   <- 0
pp_nonhab   <- c(40, 60, 100)

summary_dir <- "./output/landcover/conversion-response"
out_dir     <- "response/response-plowprint"

dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path("./processeddatasets", out_dir), recursive = TRUE, showWarnings = FALSE)

# LUTS 2020 classes retained in the analysis domain
domain_labels <- c(
  "Prairie Grassland", "Boreal Grassland", "Montane Grassland",
  "Temperate Grassland", "Alpine Tundra Grassland",
  "Prairie Wetland", "Boreal Wetland", "Temperate Wetland",
  "Cropland", "Land Converted to Cropland",
  "Annual Cropland", "Land Converted to Annual Cropland"
)

# ---- Load ----
sa_r <- rast(file.path(sa_filepath, sa_raster_filename))

plowprint <- rast(tagged_path("./processeddatasets", plowprint_filepath,
                              plowprint_filename, "modelling-extent"))

lulc2020 <- rast(tagged_path(file.path("./processeddatasets", aafc_filepath), "cropped",
                             aafc2020_filename, "modelling-extent"))

stopifnot(compareGeom(plowprint, lulc2020, stopOnError = FALSE))

# ---- Restrict to the study area (Prairie ecozone) ----
pp_sa   <- mask(crop(plowprint, sa_r), sa_r)
luts_sa <- mask(crop(lulc2020,  sa_r), sa_r)

# ---- Domain: LUTS grassland, wetland and cropland only ----
rat <- levels(luts_sa)[[1]]
lab_col <- names(rat)[2]                      # 2nd column holds the labels
domain_codes <- rat[[1]][rat[[lab_col]] %in% domain_labels]

stopifnot(length(domain_codes) > 0)
setdiff(domain_labels, rat[[lab_col]])        # any label that didn't match

in_domain <- ifel(luts_sa %in% domain_codes, 1, NA)

# ---- Response ----
pp_binary <- ifel(pp_sa %in% pp_plowed, 1,
                  ifel(pp_sa == pp_intact, 0, NA))    # 40/60/100 -> NA

response_pp <- mask(pp_binary, in_domain)
names(response_pp) <- "plowed"

write_rast_tagged(response_pp, "./processeddatasets", out_dir,
                  "response_plowprint.tif", "30m",
                  datatype = "INT1U", NAflag = 255, overwrite = TRUE)

# 300 m version for quick diagnostics. Near-sampled onto a fixed grid so each
# cell is a point sample of the same ground location, consistent with the ACI
# 300 m layers -- a 1-in-100 systematic sample, not a block summary.
sa_r_300 <- project(sa_r, aggregate(sa_r, fact = 10), method = "near")
response_pp_300 <- project(response_pp, sa_r_300, method = "near")

write_rast_tagged(response_pp_300, "./processeddatasets", out_dir,
                  "response_plowprint.tif", "300m",
                  datatype = "INT1U", NAflag = 255, overwrite = TRUE)

# ---- QA/QC ----
f <- freq(response_pp)
print(f)

n_pos <- f$count[f$value == 1]
n_neg <- f$count[f$value == 0]

cat("base rate: ", round(n_pos / (n_pos + n_neg), 4), "\n", sep = "")

# Composition of each class, to document what the model is actually separating
comp <- function(v, nm) {
  z <- freq(mask(luts_sa, ifel(response_pp == v, 1, NA)))
  data.frame(class = nm, luts = as.character(z$value), count = z$count)
}
composition <- rbind(comp(1, "plowed"), comp(0, "intact"))
composition$share <- with(composition, ave(count, class, FUN = \(x) x / sum(x)))
print(composition[order(composition$class, -composition$count), ])

summary_pp <- data.frame(
  parameter = c("plowprint_positive_codes", "plowprint_intact_code",
                "plowprint_excluded_codes", "domain_source", "domain_classes",
                "study_area", "n_positive", "n_negative", "base_rate",
                "vintage", "date_run", "script"),
  value = c(
    paste(pp_plowed, collapse = ","),
    pp_intact,
    paste(pp_nonhab, collapse = ","),
    "AAFC LUTS 2020",
    paste(domain_labels, collapse = "; "),
    "Prairie Ecozone (sa_r)",
    n_pos, n_neg, round(n_pos / (n_pos + n_neg), 4),
    "WWF Plowprint 2024 (O&R used 2016)",
    as.character(Sys.Date()),
    "5-response-plowprint.R"
  )
)
write.csv(summary_pp, file.path(summary_dir, "plowprint_response_summary.csv"),
          row.names = FALSE)
write.csv(composition, file.path(summary_dir, "plowprint_response_composition.csv"),
          row.names = FALSE)