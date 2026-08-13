# =============================================================================
# 5-response-sensitivity.R
# =============================================================================
#
# Purpose:
#   Builds alternative versions of the grass-to-crop conversion response
#   variable under different definitional choices, so the sensitivity of the
#   Random Forest model to those choices can be tested.
#
#   This script produces RESPONSE LAYERS ONLY. It does not fit models, and it
#   does not touch ./processeddatasets/response/. All outputs go to
#   ./processeddatasets/response-sensitivity/.
#
# Variants:
#   min2      Tame = LUTS 2010 cropland AND ACI 122 in >= 2 of 2009-2012.
#             Outcome = LUTS 2020 annual cropland. Perennial cropland scored 0.
#             This was the original main definition; retained as a control and
#             as the "score perennial as 0" alternative to excl_per.
#   min3      As min2 but requiring >= 3 tame years. Tests the threshold, which
#             is a judgment call rather than a break found in the data (the
#             evidence distribution is monotone increasing: 13.4 / 18.5 / 27.4 /
#             40.7% at n_tame = 1/2/3/4 within LUTS cropland). Drops the tame
#             stratum by 33%; the removed n=2 group converts at 28.3% against
#             12.0% for the retained pixels, i.e. it carries the rotational-hay
#             signature. This is the largest and most informative perturbation.
#   aci110    As min2, additionally admitting LUTS 2010 cropland pixels coded
#             ACI 110 in all four years as tame. Tests whether the excluded
#             110-in-cropland pool belongs in the baseline. Adds only 2.8% to
#             the tame stratum, so read a null ranking result cautiously: the
#             perturbation may simply be too small to move anything.
#   broadout  As min2 but counting LUTS 2020 "Cropland" and "Land Converted to
#             Cropland" (perennial forage) as conversion. DIAGNOSTIC ONLY --
#             tame pixels were already LUTS cropland in 2010, so a tame pixel
#             that stayed tame is scored 1 and the tame rate hits 99.5%. The
#             variant is uninformative as a model target; it exists to
#             demonstrate that the perennial classes mean "converted" for
#             native and "unchanged" for tame, hence excl_per's origin
#             conditioning. The native half is the meaningful part
#             (203,918 -> 625,847, a 3.1x increase).
#   excl_per  As min2, but native pixels whose 2020 class is perennial cropland
#             are set to NA rather than 0. Those ~418k pixels are native sod
#             that was broken and seeded to forage: labelling them
#             "did not convert" is false, and at 2:1 against the ~204k true
#             native positives it is material contamination of the negative
#             class. NA says "we do not model this transition"; 0 says "this
#             land resisted conversion". Only the first is true. Tame pixels
#             are untouched -- perennial cropland in 2020 genuinely means a
#             tame pixel did not convert to annual crop.
#
# Inputs:
#   - AAFC LUTS 2010 binaries : native, tame&crop (cropped to modelling extent)
#   - AAFC LUTS 2020          : binary annual crop, plus the factor raster
#   - AAFC ACI 2009-2012, 30 m, aligned to the modelling extent
#   - Study area raster (sa_r), ecoregions
#
# Outputs (./processeddatasets/response-sensitivity/):
#   - response_{id}_30m.tif          : binary response per variant (NA = excluded)
#   - grassland_origin_{id}_30m.tif  : 1 = native (LUTS), 2 = tame
#   (./output/landcover/conversion-response/):
#   - sensitivity_summary.csv        : counts and rates by variant and origin
#   - sensitivity_agreement.csv      : pixel agreement of each variant vs min2
#   - sensitivity_by_ecoregion.csv   : per-variant ecoregion rates
#
# Notes:
#   - ACI 130 (too wet to be seeded) is treated as UNOBSERVED in the evidence
#     count, not as evidence against tame. It is a wet-year artifact (~467k
#     cells in 2010, ~577k in 2011 at 300 m) and lands overwhelmingly on
#     cropland, affecting only 3.9% / 7.5% of the tame pool in those years.
#   - The baseline window (2009-2012) is kept DISJOINT from the outcome window
#     (LUTS 2020) in every variant, so baseline membership is never conditioned
#     on the outcome. Requiring persistence THROUGH the window would delete
#     early converters, which are the most informative positives.
#   - Consequence: the tame stratum's effective conversion interval is
#     2012-2020, while the native stratum's is 2010-2020.
#   - Rate denominators count pixels present in the response (non-NA), not all
#     grassland pixels. These differ for excl_per and must not be confused.
#   - The comparison that matters for CLUES is whether the RANKING of the
#     predicted surface shifts, not whether AUC changes. Fit each variant, then
#     compare prediction surfaces by Spearman correlation and top-decile
#     agreement. That happens downstream, not here.
#   - Runtime is several hours at 30 m for all five variants.
#
# Dependencies:
#   - 0-setup.R (parameters, file paths, functions)
# =============================================================================

source(file.path("scripts/landcover", "0-setup.R"))

# ---- Parameters & paths ----
aci_years   <- 2009:2012
tame_code   <- 122
native_code <- 110   # ACI native grassland; used only by the aci110 variant
unobs_code  <- 130   # too wet to be seeded: unobserved, not evidence against

summary_dir <- "./output/landcover/conversion-response"
out_dir     <- "response/response-sensitivity"

dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path("./processeddatasets", out_dir), recursive = TRUE, showWarnings = FALSE)

# ---- Variant definitions ----
variants <- tibble::tribble(
  ~id,         ~tame_min_years, ~admit_aci110, ~outcome,
  "min2",      2,               FALSE,         "annual",
  "min3",      3,               FALSE,         "annual",
  "aci110",    2,               TRUE,          "annual",
  "broadout",  2,               FALSE,         "broad",
  "excl_per",  2,               FALSE,         "annual_excl_perennial"
)

# ---- Load files ----
sa_r              <- rast(file.path(sa_filepath, sa_raster_filename))
native2010        <- rast(file.path("./processeddatasets", aafc_filepath, "binary", "AAFC2010_native.tif"))
luts_tamecrop2010 <- rast(file.path("./processeddatasets", aafc_filepath, "binary", "AAFC2010_tame&crop.tif"))
crop2020          <- rast(file.path("./processeddatasets", aafc_filepath, "binary", "AAFC2020_crop.tif"))

lulc2020 <- rast(tagged_path(file.path("./processeddatasets", aafc_filepath), "cropped",
                             aafc2020_filename, "modelling-extent"))

aci <- rast(file.path("./processeddatasets", file.path(aci_filepath, "30m"),
                      sprintf("aci_%d_30m_modelling-extent.tif", aci_years)))
names(aci) <- paste0("y", aci_years)

ecoregions <- st_read(tagged_path("./processeddatasets", ecoregions_filepath,
                                  ecoregions_filename, "modelling-extent"))

stopifnot(compareGeom(aci, native2010, stopOnError = FALSE))
stopifnot(compareGeom(crop2020, native2010, stopOnError = FALSE))

# ---- Evidence and outcome layers (computed once, reused across variants) ----
tame_evidence   <- sum(aci == tame_code,   na.rm = TRUE)
native_evidence <- sum(aci == native_code, na.rm = TRUE)

# LUTS 2020 separates annual from perennial cropland; LUTS 2010 does not. That
# asymmetry is why ACI is needed on the baseline side but not on the outcome.
perennial2020  <- ifel(lulc2020 == "Cropland" |
                         lulc2020 == "Land Converted to Cropland", 1, 0)
crop2020_broad <- ifel(crop2020 == 1 | perennial2020 == 1, 1, 0)

# ---- Ecoregion zones (Prairie ecozone only) ----
ecoregions_r <- rasterize(vect(ecoregions |> filter(ECOZONE_ID == 10)),
                          sa_r, field = "ECOREGION1")

zsum <- function(r, nm) setNames(zonal(r, ecoregions_r, fun = "sum", na.rm = TRUE),
                                 c("ecoregion", nm))

# ---- Build each variant ----
summary_rows   <- list()
agreement_rows <- list()
ecoregion_rows <- list()

for (i in seq_len(nrow(variants))) {
  v  <- variants[i, ]
  id <- v$id
  message("variant: ", id)
  
  # ---- Baseline ----
  tame <- ifel(luts_tamecrop2010 == 1 & tame_evidence >= v$tame_min_years, 1, 0)
  
  if (v$admit_aci110) {
    tame <- ifel(luts_tamecrop2010 == 1 & native_evidence == length(aci_years), 1, tame)
  }
  
  grass  <- ifel(native2010 == 1 | tame == 1, 1, 0)
  origin <- ifel(native2010 == 1, 1, ifel(tame == 1, 2, NA))
  names(grass)  <- "grass2010"
  names(origin) <- "grass_origin"
  
  # ---- Outcome ----
  # excl_per is origin-conditioned: perennial cropland in 2020 means CONVERTED
  # for a native pixel but UNCHANGED for a tame one, so only native pixels are
  # mislabelled by scoring 0 and only those are dropped.
  #
  # The !is.na(origin) guard is REQUIRED: terra propagates NA through `&`
  # rather than short-circuiting, so `origin == 1` returns NA (not FALSE)
  # outside grassland, and without the guard the NA propagates into `outcome`
  # and then into `response`, blanking ~98.5k cells that should be 0.
  outcome <- switch(v$outcome,
                    annual                = crop2020,
                    broad                 = crop2020_broad,
                    annual_excl_perennial = ifel(!is.na(origin) & origin == 1 & perennial2020 == 1,
                                                 NA, crop2020),
                    stop("unknown outcome: ", v$outcome)
  )
  
  response <- ifel(grass == 1 & outcome == 1, 1, 0)
  names(response) <- "conversion"
  
  write_rast_tagged(response, "./processeddatasets", out_dir,
                    sprintf("response_%s.tif", id), "30m",
                    datatype = "INT1U", NAflag = 255, overwrite = TRUE)
  write_rast_tagged(origin, "./processeddatasets", out_dir,
                    sprintf("grassland_origin_%s.tif", id), "30m",
                    datatype = "INT1U", NAflag = 255, overwrite = TRUE)
  
  # ---- Summary by origin ----
  # Denominator counts pixels present in the response, not all grassland: these
  # differ for excl_per, where ~418k native pixels are NA.
  s <- zonal(response, origin, fun = "sum", na.rm = TRUE)
  names(s) <- c("origin", "n_conversions")
  s$n_grassland <- zonal(!is.na(response), origin, fun = "sum", na.rm = TRUE)[, 2]
  s <- s |>
    mutate(variant         = id,
           tame_min_years  = v$tame_min_years,
           admit_aci110    = v$admit_aci110,
           outcome         = v$outcome,
           origin          = c("native (LUTS)", "tame (ACI)")[origin],
           conversion_rate = n_conversions / n_grassland) |>
    select(variant, tame_min_years, admit_aci110, outcome,
           origin, n_grassland, n_conversions, conversion_rate)
  
  summary_rows[[id]] <- s
  print(s)
  
  # ---- Agreement with min2 ----
  if (id == "min2") {
    grass_ref    <- grass
    response_ref <- response
  } else {
    stopifnot(exists("grass_ref"))   # min2 must run first
    agreement_rows[[id]] <- data.frame(
      variant          = id,
      grass_added      = as.numeric(global(grass == 1 & grass_ref == 0, "sum", na.rm = TRUE)),
      grass_removed    = as.numeric(global(grass == 0 & grass_ref == 1, "sum", na.rm = TRUE)),
      response_added   = as.numeric(global(response == 1 & response_ref == 0, "sum", na.rm = TRUE)),
      response_removed = as.numeric(global(response == 0 & response_ref == 1, "sum", na.rm = TRUE)),
      response_na      = as.numeric(global(is.na(response) & !is.na(response_ref), "sum", na.rm = TRUE))
    )
  }
  
  # ---- Ecoregion rates ----
  response_p <- mask(crop(response, sa_r), sa_r)
  grass_p    <- mask(crop(grass,    sa_r), sa_r)
  origin_p   <- mask(crop(origin,   sa_r), sa_r)
  # Must be restricted to grassland: response is 0 (not NA) outside the
  # baseline, so !is.na(response_p) alone counts the entire ecoregion.
  inresp_p   <- !is.na(response_p) & grass_p == 1
  
  e <- zsum(response_p, "n_conversions") |>
    left_join(zsum(inresp_p, "n_grassland"), by = "ecoregion") |>
    left_join(zsum(ifel(origin_p == 2, response_p, NA), "n_conv_tame"),  by = "ecoregion") |>
    left_join(zsum(ifel(origin_p == 2, inresp_p,   NA), "n_grass_tame"), by = "ecoregion") |>
    mutate(variant        = id,
           n_conv_native  = n_conversions - n_conv_tame,
           n_grass_native = n_grassland   - n_grass_tame,
           rate_all       = n_conversions / n_grassland,
           rate_tame      = n_conv_tame   / n_grass_tame,
           rate_native    = n_conv_native / n_grass_native)
  
  stopifnot(all(e$n_conv_native >= 0), all(e$n_grass_native >= 0))
  ecoregion_rows[[id]] <- e
}

# ---- Write comparison tables ----
sensitivity_summary <- bind_rows(summary_rows)
write.csv(sensitivity_summary, file.path(summary_dir, "sensitivity_summary.csv"), row.names = FALSE)
print(sensitivity_summary)

sensitivity_agreement <- bind_rows(agreement_rows)
write.csv(sensitivity_agreement, file.path(summary_dir, "sensitivity_agreement.csv"), row.names = FALSE)
print(sensitivity_agreement)

sensitivity_by_ecoregion <- bind_rows(ecoregion_rows) |>
  select(variant, ecoregion, n_grassland, n_conversions,
         n_grass_tame, n_conv_tame, n_grass_native, n_conv_native,
         rate_all, rate_tame, rate_native) |>
  arrange(variant, desc(rate_tame))
write.csv(sensitivity_by_ecoregion, file.path(summary_dir, "sensitivity_by_ecoregion.csv"), row.names = FALSE)

# ---- Summary values (30 m, modelling extent) ----
#   variant   origin  n_grassland  n_conversions     rate
#   min2      native   92,380,345        203,918  0.00221
#   min2      tame     80,956,335     14,091,827  0.17407
#   min3      native   92,380,345        203,918  0.00221
#   min3      tame     53,987,560      6,459,966  0.11966
#   aci110    native   92,380,345        203,918  0.00221
#   aci110    tame     83,190,958     14,267,741  0.17151
#   broadout  native   92,380,345        625,847  0.00677
#   broadout  tame     80,956,335     80,520,685  0.99462   <- diagnostic only
#   excl_per  native   91,958,416        203,918  0.00222
#   excl_per  tame     80,956,335     14,091,827  0.17407
#
# excl_per drops 421,929 native pixels relative to min2 (not the 418,353 in the
# LUTS crosstab -- that figure is Prairie Grassland only, and native also
# includes Boreal, Montane, Temperate and Alpine grassland). Verified against
# min2: NA footprints differ by exactly 421,929 and no shared pixel differs in
# value.