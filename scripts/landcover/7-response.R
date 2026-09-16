# Purpose:
#   Derives the binary response variable for the Random Forest conversion risk
#   model: grassland pixels in 2010 that converted to ANNUAL cropland by 2020.
#   Also summarises conversion by ecoregion and produces a map.
#
#   Grassland in 2010 is the union of two strata:
#     (a) NATIVE -- native grassland in AAFC LUTS 2010.
#     (b) TAME   -- AAFC LUTS 2010 Cropland that is coded Pasture/forages
#                   (ACI class 122) in at least 2 of 2009-2012.
#   LUTS 2010 has no tame grassland class, which is the entire reason ACI is
#   in this pipeline. LUTS 2020 DOES separate annual from perennial cropland,
#   so no ACI is needed on the outcome side.
#
# Inputs:
#   - AAFC LUTS 2010 binaries : AAFC2010_native.tif, AAFC2010_tame&crop.tif
#   - AAFC LUTS 2020 binary   : AAFC2020_crop.tif (annual cropland)
#   - AAFC ACI 2009-2012, 30 m, aligned to the modelling extent
#   - Study area raster (sa_r) and shapefile (sa); ecozones and ecoregions
#
# Outputs:
#   - AAFC2010_ACI2009-2012_grass.tif  : 1 = in the 2010 grassland baseline
#   - grassland_origin_2010_30m.tif    : 1 = native (LUTS), 2 = tame (ACI)
#   - response_30m.tif                 : 1 = converted, 0 = not converted
#   - conversion_by_origin.csv, conversion_by_ecoregion.csv
#   - response_definition.csv          : parameter provenance for this run
#   - grassland_conversion_to_cropland_map.png
#
# -----------------------------------------------------------------------------
# TAME GRASSLAND DEFINITION
# -----------------------------------------------------------------------------
#   A 2010 pixel is TAME if it is LUTS 2010 Cropland AND coded ACI 122 in at
#   least `tame_min_years` (= 2) of 2009-2012. ACI 130 (too wet to be seeded)
#   is treated as UNOBSERVED, not as evidence against.
#
#   Why a multi-year rule:
#     A single 122 year does not distinguish permanent pasture from the forage
#     phase of an annual rotation. Single-year 122 pixels convert with near
#     certainty and would form a large, low-information positive class that
#     teaches the model something close to a tautology.
#
#   Why the windows are disjoint:
#     The baseline window (2009-2012) never overlaps the outcome (LUTS 2020),
#     so baseline membership is not conditioned on the outcome. Note the rule
#     counts EVIDENCE, not persistence: a pixel tame in 2009-2010 and broken in
#     2011 still qualifies at 2 tame years and is retained as a positive.
#     Requiring persistence THROUGH the window would delete early converters,
#     which are the most informative positives.
#
#   Why 130 is unobserved:
#     130 is a wet-year artifact (~467k cells in 2010, ~577k in 2011 at 300 m
#     vs ~15-55k otherwise -- the prairie flood years). "Too wet to be seeded"
#     is a cropland concept and rarely lands on forage: within the tame pool it
#     affects only 3.9% (2010) and 7.5% (2011) of pixels.
#
#   Why >= 2 and not 3:
#     The evidence distribution within LUTS cropland is monotone increasing
#     (13.4 / 18.5 / 27.4 / 40.7% at n_tame = 1/2/3/4), so there is NO natural
#     break. The cut is a judgment call, not a discovered boundary. >= 2 drops
#     the 13.4% with the weakest evidence; >= 3 would drop 31.9%. The removed
#     n=2 group converts at 28.3% against 12.0% for the retained pixels, i.e.
#     it carries the rotational-hay signature -- so >= 3 is a large
#     perturbation and is run as variant "min3" in 5-response-sensitivity.R.
#     Model-ranking stability across the two must be reported.
#
#   Why ACI 110 is NOT used:
#     LUTS defines native; ACI defines tame. Each product is used where it is
#     authoritative. ACI 110 within LUTS native grassland is redundant (those
#     pixels are already in the baseline). ACI 110 within LUTS CROPLAND
#     (153,796 cells at 300 m) was tested three ways and rejected: the
#     persistence distribution is DECREASING (38/30/17/15%, the mirror of the
#     tame pool), and even the fully persistent n=4 subgroup converts at 8.4%
#     against 0.22% for the native stratum -- 39x higher. They do not behave
#     like native grassland, and they would have no coherent grass_origin
#     (not native, since LUTS disagrees; not tame, since ACI 122 disagrees).
#     Tested as variant "aci110" in the sensitivity script.
#
#   Validation:
#     ACI 122 and 110 land in near-opposite places in LUTS (78% vs 8% cropland;
#     8% vs 68% native), so two independently produced products agree on the
#     native/tame axis. Of ACI 122 pixels, the remaining 14% fall in wetland,
#     forest, settlement or water -- ACI over-calling 122 on spectrally
#     hay-like cover (roads alone are ~19k cells at 300 m, a 56 m mixed-pixel
#     artifact). Excluded by design.
#
# -----------------------------------------------------------------------------
# OUTCOME DEFINITION AND SCOPE
# -----------------------------------------------------------------------------
#   Conversion = LUTS 2020 "Annual Cropland" or "Land Converted to Annual
#   Cropland". LUTS routes all change through the "Land Converted to" classes,
#   so Prairie Grassland -> Annual Cropland is zero by construction.
#
#   EXCLUDED BY DESIGN: native grassland -> "Land Converted to Cropland"
#   (perennial forage), 418,353 pixels, roughly 2x the annual signal. This is
#   real, irreversible grassland loss, but it is a transition CLUES does not
#   allocate (CLUES moves land into annual cropland or developed only).
#   Including it would train the model on a process the simulation cannot act
#   on. Tested as variant "broadout".
#
#   KNOWN ISSUE: those 418k pixels are currently scored 0, i.e. "did not
#   convert", which is false -- the sod was broken. Against ~204k true native
#   positives that is 2:1 contamination of the native negative class. Variant
#   "excl_per" sets them to NA instead (421,929 pixels once Boreal, Montane,
#   Temperate and Alpine grassland are included), which asserts nothing false
#   while keeping the CLUES alignment. excl_per is a candidate to REPLACE this
#   definition as the main response; see 5-response-sensitivity.R.
#
#   The exclusion must be origin-conditioned. Perennial cropland in 2020 means
#   CONVERTED for a native pixel but UNCHANGED for a tame one (tame pixels were
#   already LUTS cropland in 2010). Applying it uniformly would strip ~660k
#   legitimate tame negatives; scoring perennial as conversion uniformly drives
#   the tame rate to 99.5% and destroys the response.
#
# -----------------------------------------------------------------------------
# RESULTS AND CONSEQUENCES
# -----------------------------------------------------------------------------
#   Baseline    173,336,680 px    Conversions  14,295,745
#     native     92,380,345 px       203,918   (0.22%)
#     tame       80,956,335 px    14,091,827   (17.4%)
#
#   - 77x base-rate gap between strata. Tame is 98.6% of all conversions, and
#     97-99.7% within every ecoregion. A map of this response is effectively a
#     map of tame grass breaking, NOT of native prairie loss -- state this
#     explicitly for conservation audiences, and report the by-origin table
#     rather than the pooled one.
#   - The gap makes pooled model evaluation unreliable: a pooled AUC can look
#     healthy while the model fits only the tame process. Score the
#     leave-one-ecodistrict-out folds SEPARATELY by grass_origin before the
#     CLUES handoff. If native collapses toward 0.5, stratum-specific models
#     plus isotonic calibration on the out-of-fold predictions are needed to
#     put both strata on a common scale (CLUES allocates by ranking).
#   - The native stratum reproduces the earlier native-only response exactly
#     (108,679 / 43,730 / 35,095 / 9,468 / 2,431 / 1,423 / 107 by ecoregion),
#     so the tame stratum is purely additive and nothing was disturbed.
#   - Tame conversion rates vary only ~3x across ecoregions (10.0-30.4%) while
#     counts vary ~130x: the count map shows where tame grass IS, the rate map
#     shows where it is under pressure. Aspen Parkland is the native outlier at
#     1.20%, 5x any other ecoregion.
#   - grass_origin is a fixed 2010 attribute and the surface is defined only on
#     the 2010 grassland domain. If CLUES needs suitability outside that domain,
#     or can move a pixel into grassland and later ask its conversion risk, an
#     origin rule for those pixels is required. Unresolved.
#
# Notes:
#   - ACI 2009 and 2010 are 56 m; 2011 onward are 30 m. The early years are
#     nearest-neighbour projected to the 30 m grid, so they carry 56 m
#     blockiness. The LUTS intersection removes most of the consequence.
#   - ACI grass/forage classes are not equally separable across the prairies:
#     Moist-Mixed and Mixed are distinguishable, Fescue is not. Response-side
#     measurement error therefore varies across LOEO folds, and the Fescue fold
#     is noisier for reasons unrelated to model transferability.
#   - The ecoregion summary covers the Prairie ecozone only (ECOZONE_ID 10);
#     ~22% of the response falls outside it, within the modelling extent.
#
# Dependencies:
#   - 0-setup.R (parameters, file paths, functions)
#   - 3-crop-data.R outputs (LUTS 2010/2020 binaries)
#   - 1h-prep-data-aci.R outputs (ACI 2009-2012, 30 m, aligned)
# Related:
#   - 7-response-sensitivity.R : variants min2 / min3 / aci110 / broadout / excl_per
# =============================================================================

source(file.path("scripts/landcover", "0-setup.R"))

# ---- Parameters & paths ----
aci_years <- 2009:2012
aci_tame_codes <- c(122) # not using 110 intentionally
tame_min_years <- 2
unobs_code     <- 130   # too wet to be seeded: treated as unobserved, not as evidence against

summary_dir <- "./output/landcover/conversion-response"

# ---- Load files ----
sa_r <- rast(file.path(sa_filepath, sa_raster_filename))
sa   <- st_read(file.path(sa_filepath, sa_filename))

native2010 <- rast(file.path("./processeddatasets", aafc_filepath, "binary", "AAFC2010_native.tif"))
luts_tamecrop2010 <- rast(file.path("./processeddatasets", aafc_filepath, "binary", "AAFC2010_tame&crop.tif"))
crop2020 <- rast(file.path("./processeddatasets", aafc_filepath, "binary", "AAFC2020_crop.tif"))

aci <- rast(file.path("./processeddatasets", file.path(aci_filepath, "30m"), sprintf("aci_%d_30m_modelling-extent.tif", aci_years)))
names(aci) <- paste0("y", aci_years)

ecozones   <- st_read(tagged_path("./processeddatasets", ecozones_filepath,   ecozones_filename,   "modelling-extent"))
ecoregions <- st_read(tagged_path("./processeddatasets", ecoregions_filepath, ecoregions_filename, "modelling-extent"))

# ---- ACI tame grass mask ----
tame_evidence <- sum(aci == aci_tame_codes, na.rm = TRUE)

# ---- Grassland baseline ~2010 ----
# LUTS cropland that ACI calls tame grass
tame2010 <- ifel(luts_tamecrop2010 == 1 & tame_evidence >= tame_min_years, 1, 0)

grass2010 <- ifel(native2010 == 1 | tame2010 == 1, 1, 0)
names(grass2010) <- "grass2010"

# Origin of each grassland pixel: 1 = native (LUTS), 2 = tame (ACI)
grass_origin <- ifel(native2010 == 1, 1, ifel(tame2010 == 1, 2, NA))
names(grass_origin) <- "grass_origin"

# ---- Response variable ----
response <- ifel(grass2010 == 1 & crop2020 == 1, 1, 0)
names(response) <- "conversion"

# ---- Save to disk ----
writeRaster(grass2010, file.path("./processeddatasets", aafc_filepath, "binary", "AAFC2010_ACI2009-2012_grass.tif"), datatype = "INT1U", overwrite = TRUE)
write_rast_tagged(grass_origin, "./processeddatasets", "predictors/30m", "grassland_origin_2010_modelling-extent.tif", "30m", datatype = "INT1U")
write_rast_tagged(response,     "./processeddatasets", "response", "response.tif", "30m", datatype = "INT1U")
#grass2010 <- rast(file.path("./processeddatasets", aafc_filepath, "binary", "AAFC2010_ACI2009-2012_grass.tif"))
#grass_origin <- rast(file.path("./processeddatasets", "predictors/30m", "grassland_origin_2010_30m.tif"))
#response <- rast(file.path("./processeddatasets", "response", "response_30m.tif"))
  
# ---- QA/QC ----
freq(grass2010)
#   layer value     count
# 1     1     0 511878693
# 2     1     1 173336680
freq(response)
#   layer value     count
# 1     1     0 671018281
# 2     1     1  14295745

# Conversions split by grassland origin
conversion_by_origin <- zonal(response, grass_origin, fun = "sum", na.rm = TRUE)
names(conversion_by_origin) <- c("origin", "n_conversions")
conversion_by_origin$n_grassland <- freq(grass_origin)$count
conversion_by_origin <- conversion_by_origin |>
  mutate(origin = c("native (LUTS)", "tame (ACI)")[origin],
         conversion_rate = n_conversions / n_grassland)

print(conversion_by_origin)
#         origin n_conversions n_grassland conversion_rate
# 1 native (LUTS)        203918    92380345     0.002207374
# 2    tame (ACI)      14091827    80956335     0.174067008
write.csv(conversion_by_origin,
          file.path(summary_dir, "conversion_by_origin.csv"), row.names = FALSE)

# ---- Record definition provenance ----
response_definition <- data.frame(
  parameter = c("aci_years", "tame_code", "tame_min_years", "unobserved_code",
                "native_source", "tame_source", "outcome_source",
                "outcome_classes", "baseline_window", "outcome_year", "resolution_m",
                "n_grassland", "n_conversions", "date_run", "script"),
  value = c(
    paste(range(aci_years), collapse = "-"),
    aci_tame_codes,
    tame_min_years,
    unobs_code,
    "AAFC LUTS 2010 native grassland classes",
    "LUTS 2010 Cropland AND ACI 122 in >= tame_min_years of aci_years",
    "AAFC LUTS 2020",
    "Annual Cropland; Land Converted to Annual Cropland",
    paste(range(aci_years), collapse = "-"),
    2020,
    res,
    freq(grass2010)$count[freq(grass2010)$value == 1],
    freq(response)$count[freq(response)$value == 1],
    as.character(Sys.Date()),
    "5-response.R"
  ),
  note = c(
    "ACI years used to establish tame grass evidence",
    "ACI class code for Pasture/forages",
    "Minimum years coded tame_code to qualify; sensitivity: 3",
    "ACI class treated as unobserved in the persistence count",
    "LUTS defines native; ACI not used (see 110 rationale in header)",
    "Baseline window disjoint from outcome window to avoid conditioning on outcome",
    "No ACI needed on outcome side: LUTS 2020 separates annual from perennial cropland",
    "Excludes 'Land Converted to Cropland' (perennial): CLUES allocates annual cropland only",
    "Tame stratum effective interval 2012-2020; native 2010-2020",
    "",
    "Native predictor/response grid",
    "Pixels in 2010 grassland baseline",
    "Pixels converted to annual cropland by 2020",
    "", ""
  )
)

write.csv(response_definition,
          file.path(summary_dir, "response_definition.csv"), row.names = FALSE)


# ---- Summarise conversion by ecoregion within the study area ----
# Filter to Prairie ecozone only (ECOZONE_ID = 10)
response_prairies     <- mask(crop(response,     sa_r), sa_r)
grass2010_prairies    <- mask(crop(grass2010,    sa_r), sa_r)
grass_origin_prairies <- mask(crop(grass_origin, sa_r), sa_r)

ecoregions_prairie <- ecoregions |> filter(ECOZONE_ID == 10)
ecoregions_r <- rasterize(vect(ecoregions_prairie), sa_r, field = "ECOREGION1")

zsum <- function(r, nm) {
  z <- zonal(r, ecoregions_r, fun = "sum", na.rm = TRUE)
  setNames(z, c("ecoregion", nm))
}

conversion_by_ecoregion <-
  zsum(response_prairies, "n_conversions") |>
  left_join(zsum(grass2010_prairies, "n_grassland"), by = "ecoregion") |>
  left_join(zsum(ifel(grass_origin_prairies == 2, response_prairies, NA), "n_conv_tame"), by = "ecoregion") |>
  left_join(zsum(ifel(grass_origin_prairies == 2, 1, NA), "n_grass_tame"), by = "ecoregion") |>
  mutate(n_conv_native  = n_conversions - n_conv_tame,
         n_grass_native = n_grassland   - n_grass_tame,
         rate_all       = n_conversions  / n_grassland,
         rate_tame      = n_conv_tame    / n_grass_tame,
         rate_native    = n_conv_native  / n_grass_native) |>
  arrange(desc(rate_tame))

print(conversion_by_ecoregion)
#                    ecoregion n_conversions n_grassland n_conv_tame n_grass_tame n_conv_native
# 1 Southwest Manitoba Uplands         44188      169864       44081       145064           107
# 2        Lake Manitoba Plain        468999     2081927      467576      1655401          1423
# 3           Fescue Grassland        289782     5160246      280314      1375715          9468
# 4             Aspen Parkland       5771402    38253498     5662723     29196594        108679
# 5      Moist Mixed Grassland       2041724    28580398     2006629     12632114         35095
# 6            Mixed Grassland       2512602    70751576     2468872     17192996         43730
# 7             Cypress Upland         92743     7023559       90312       902569          2431
# n_grass_native   rate_all rate_tame  rate_native
# 1          24800 0.26013752 0.3038728 0.0043145161
# 2         426526 0.22527159 0.2824548 0.0033362562
# 3        3784531 0.05615663 0.2037588 0.0025017631
# 4        9056904 0.15087253 0.1939515 0.0119995751
# 5       15948284 0.07143791 0.1588514 0.0022005502
# 6       53558580 0.03551302 0.1435975 0.0008164892
# 7        6120990 0.01320456 0.1000610 0.0003971580
write.csv(conversion_by_ecoregion,
          file.path(summary_dir, "conversion_by_ecoregion.csv"), row.names = FALSE)

# Sanity checks
summary(conversion_by_ecoregion$n_conversions)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 44188  191263  468999 1603063 2277163 5771402 
sum(conversion_by_ecoregion$n_conversions < 100, na.rm = TRUE)   # expect 0
sum(conversion_by_ecoregion$n_conversions == 0,  na.rm = TRUE)   # expect 0


# ---- Map ----
# Aggregate to 3x target resolution for plotting efficiency
response_plot <- aggregate(response_prairies, fact = 20, fun = "max")
response_df   <- as.data.frame(response_plot, xy = TRUE) |>
  filter(conversion == 1)

ecoregions_sf <- st_as_sf(ecoregions) |> st_intersection(sa)
ecozone_sf    <- st_as_sf(ecozones)   |> st_intersection(sa)

ecoregions_label <- ecoregions_sf |>
  left_join(conversion_by_ecoregion, by = c("ECOREGION1" = "ecoregion")) |>
  filter(n_conversions >= 1000)

map <- ggplot() +
  geom_tile(data = response_df, aes(x = x, y = y),
            fill = "darkred", alpha = 0.8) +
  geom_sf(data = ecoregions_sf, fill = NA, colour = "grey40", linewidth = 0.3) +
  geom_sf(data = ecozone_sf,    fill = NA, colour = "black",  linewidth = 0.7) +
  geom_sf_label(data = ecoregions_label, aes(label = ECOREGION1),
                size = 2.5, fill = "white", alpha = 0.7, label.size = 0.1) +
  theme_void()

ggsave(file.path("./output/landcover/conversion-response", "grassland_conversion_to_cropland_map.png"),
       plot = map, width = 10, height = 7, dpi = 300, bg = "white")
