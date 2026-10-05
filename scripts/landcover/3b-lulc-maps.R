# Purpose:
#   Builds the two LULC state maps for model training: LULC 2010
#   and LULC 2020. This is the script that decides
#   what counts as native grassland, tame grassland, annual cropland or
#   settlement. Predictors (4d) and the response (7) read these maps.
#
#   Combining AAFC LUTS and AAFC: LUTS bounds the domains, ACI splits one of them.
#     - LUTS defines native grassland, the agricultural domain, settlement and
#       everything else, at both dates.
#     - ACI does one job: split the LUTS agricultural domain into tame
#       grassland and annual cropland, using classes 110 and 122 only
#   LUTS 2010 has no tame class, which is why ACI is used. LUTS 2020 does 
#   split annual from perennial, but its split is not used here, instead ACI
#   splits the union of the two, so both maps are built by the same rule. LUTS
#   2020's split is a response variant in 7.
#
# Class scheme (identical in both years, so the maps cross-tabulate):
#   1 = native grassland   2 = tame grassland   3 = annual cropland
#   4 = settlement         5 = other            NA = outside modelling extent
#
#     conversion = lulc_2010 %in% c(1,2) & lulc_2020 == 3
#     excluded   = (lulc_2010 == 1 & lulc_2020 == 2) |
#                  (lulc_2010 %in% c(1,2) & lulc_2020 == 4)   -> NA
#
# Tame definition:
#   A pixel is TAME if it is inside the LUTS agricultural domain for its year
#   AND coded 110 (grass/prairie) or 122 (pasture/forage) by ACI in >= tame_min_years of that year's window.
#   The rule counts EVIDENCE, not persistence, so a pixel tame in 2009-2010 and
#   broken in 2011 still qualifies and is retained - early converters are the
#   most informative positives. Pixels below the threshold fall through to
#   annual cropland, which is the intended rotational-hay treatment.
#   The windows (2009-2012, 2019-2022) are disjoint so baseline membership is
#   never conditioned on the outcome.
#
# Domain:
#   luts_ag_2010 and luts_ag_2020 are from LUTS. They must not reference any
#   ACI-derived object, and their totals cannot vary with aci_tame_codes:
#     luts_ag_2010 = 386,424,415 px    luts_ag_2020 = 387,177,855 px
#   Both go to lulc_definition.csv. `other` is LUTS-defined and must also match
#   across variants (193,153,740 px at 2020).
#
# Notes:
#   - Class codes derived from 0a-setup.R
#   - ACI 2009-2010 are 56 m, nearest-neighbour projected to 30 m; the t1
#     window is native 30 m. The two windows differ in effective grain.
#   - AAFC re-releases individual province-years. aci_url_index.csv from 1h
#     pins which file each year came from; commit it alongside the scripts.
#
# Inputs:
#   - LUTS 2010 binaries: native, tame&crop, settlement
#   - LUTS 2020 binaries: native, annual cropland, perennial cropland, settlement
#   - ACI 2009-2012 and 2019-2022, 30 m, aligned to the modelling extent
# Outputs:
#   - lulc_2010_30m.tif, lulc_2020_30m.tif
#   - aci_tame_evidence_{2010,2020}_30m.tif
#   - grassland_origin_2010_30m.tif, AAFC2010_ACI2009-2012_grass.tif
#   - lulc_definition.csv
# Dependencies:  0a-setup.R, 3a-crop-data.R, 1h-prep-data-aci.R
# Used by:   4d-predictors-landscape-context.R (lulc_2010 only), 7-response.R
# =============================================================================

source(file.path("scripts/landcover", "0a-setup.R"))

# ---- Parameters ----
tame_min_years <- 2        # shared by both timesteps

# ACI evidence windows
aci_years_2010 <- 2009:2012  # brackets 2010
aci_years_2020 <- 2019:2022  # brackets 2020; see header, window choice is open

# ACI codes counted as tame evidence
#aci_tame_codes <- c(110, 122)
aci_tame_codes <- 122

# ---- Paths ----
binary_dir  <- file.path("./processeddatasets", aafc_filepath, "binary")
aci_dir <- file.path("./processeddatasets", file.path(aci_filepath, "30m"))
summary_dir <- "./output/landcover/lulc"
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

f_2010_native    <- "AAFC2010_native_studyarea.tif"
f_2010_tamecrop  <- "AAFC2010_tame&crop_studyarea.tif"
f_2010_settlement <- "AAFC2010_settlement_studyarea.tif"
f_2020_native    <- "AAFC2020_native_studyarea.tif"
f_2020_annual    <- "AAFC2020_crop_studyarea.tif"        # Annual Cropland + Land Converted to Annual Cropland
f_2020_tame <- "AAFC2020_tame_studyarea.tif"   # Cropland + Land Converted to Cropland
f_2020_settlement <- "AAFC2020_settlement_studyarea.tif"

# ---- Load ----
# Boundaries (study area)
ecoregions_sa  <- st_read(tagged_path("./processeddatasets", ecoregions_filepath, ecoregions_filename, "studyarea"))
ecoregions_r   <- rast(tagged_path("./processeddatasets", ecoregions_filepath, paste0(tools::file_path_sans_ext(ecoregions_filename), ".tif"), "studyarea"))
ecodistricts_sa  <- st_read(tagged_path("./processeddatasets", ecodistricts_filepath, ecodistricts_filename, "studyarea"))
ecodistricts_r <- rast(tagged_path("./processeddatasets", ecodistricts_filepath, paste0(tools::file_path_sans_ext(ecodistricts_filename), ".tif"), "studyarea"))

# LULC (study area)
lulc2025       <- rast(tagged_path("./processeddatasets", grass_consensus_filepath, grass_consensus_filename, "studyarea"))
luts2010       <- rast(tagged_path(file.path("./processeddatasets", aafc_filepath), "cropped", aafc2010_filename, "studyarea"))
luts2020       <- rast(tagged_path(file.path("./processeddatasets", aafc_filepath), "cropped", aafc2020_filename, "studyarea"))
native2010     <- rast(file.path(binary_dir, f_2010_native))
tamecrop2010   <- rast(file.path(binary_dir, f_2010_tamecrop))
settlement2010 <- rast(file.path(binary_dir, f_2010_settlement))
# define ag in 2010 as undifferentiated Cropland (tamecrop)
luts_ag_2010   <- tamecrop2010 == 1

native2020     <- rast(file.path(binary_dir, f_2020_native))
annual2020     <- rast(file.path(binary_dir, f_2020_annual))
tame2020_luts  <- rast(file.path(binary_dir, f_2020_tame))
settlement2020 <- rast(file.path(binary_dir, f_2020_settlement))
# define ag in 2020 as union of annual cropland and tame grass
luts_ag_2020   <- annual2020 == 1 | tame2020_luts == 1

# ACI
aci_2010 <- rast(file.path(aci_dir, sprintf("aci_%d_30m_modelling-extent.tif", aci_years_2010)))
names(aci_2010) <- paste0("y", aci_years_2010)

aci_2020 <- rast(file.path(aci_dir, sprintf("aci_%d_30m_modelling-extent.tif", aci_years_2020)))
names(aci_2020) <- paste0("y", aci_years_2020)


#--- Check grassland extents and frequencies in raw LUTS ----
# LUTS 2010
freq(luts2010)
# layer                                                     value     count
# 1      1                                                Settlement   1935475
# 2      1                               High Reflectance Settlement   1057398
# 3      1                                         Settlement Forest    225493
# 4      1                                                     Roads  11006795
# 5      1                                      Vegetated Settlement   3526845
# 6      1                          Very High Reflectance Settlement    464057
# 7      1                                                     Water  19659907
# 8      1                                                    Forest  16445377
# 9      1                                            Forest Wetland   8311269
# 10     1               Forest Regenerating after Harvest <20 years    929666
# 11     1       Forest Wetland Regenerating after Harvest <20 years     92943
# 12     1             Forest Regenerating after Harvest 20-29 years      9194
# 13     1     Forest Wetland Regenerating after Harvest 20-29 years      5721
# 14     1                  Forest Regenerating after Fire <20 years     21503
# 15     1                                                  Cropland 332667091
# 16     1                                Land Converted to Cropland    513266
# 17     1                                         Prairie Grassland  88484274
# 18     1                                          Boreal Grassland     57683
# 19     1                                         Montane Grassland    382687
# 20     1                                           Prairie Wetland  29095288
# 21     1                                            Boreal Wetland    462382
# 22     1                       Newly-Detected Settlement <10 years    769767
# 23     1      Newly-Detected High Reflectance Settlement <10 years    264193
# 24     1                Newly-Detected Settlement Forest <10 years     56106
# 25     1                             Newly-Detected Road <10 years     75838
# 26     1             Newly-Detected Vegetated Settlement <10 years   1152973
# 27     1 Newly-Detected Very High Reflectance Settlement <10 years    176131
# 28     1                                                Other Land    715305

freq(native2010)
# layer value     count
# 1     1     0 429639983
# 2     1     1  88924644
# There are 8,003,218 ha of native grassland in the LUTS 2010 in the study area raster
freq(tamecrop2010)
# layer value     count
# 1     1     0 185384270
# 2     1     1 333180357
freq(settlement2010)
# layer value     count
# 1     1     0 508936189
# 2     1     1   9628438

# LUTS 2020
freq(luts2020)
# layer                                                     value     count
# 1      1                                                Settlement   2686590
# 2      1                               High Reflectance Settlement   1342799
# 3      1                                         Settlement Forest    285393
# 4      1                                                     Roads  11082633
# 5      1                                      Vegetated Settlement   4667407
# 6      1                          Very High Reflectance Settlement    646276
# 7      1                                                     Water  19665425
# 8      1                                                    Forest  16191316
# 9      1                                            Forest Wetland   8218524
# 10     1               Forest Regenerating after Harvest <20 years    976149
# 11     1       Forest Wetland Regenerating after Harvest <20 years    113924
# 12     1             Forest Regenerating after Harvest 20-29 years     61790
# 13     1     Forest Wetland Regenerating after Harvest 20-29 years     10755
# 14     1                  Forest Regenerating after Fire <20 years     16207
# 15     1                                                  Cropland  82662765
# 16     1                                           Annual Cropland 249537392
# 17     1                                Land Converted to Cropland    890596
# 18     1                         Land Converted to Annual Cropland    787669
# 19     1                                         Prairie Grassland  87566282
# 20     1                                          Boreal Grassland     56957
# 21     1                                         Montane Grassland    382363
# 22     1                                           Prairie Wetland  28181710
# 23     1                                            Boreal Wetland    459784
# 24     1                       Newly-Detected Settlement <10 years     58471
# 25     1      Newly-Detected High Reflectance Settlement <10 years    213190
# 26     1                Newly-Detected Settlement Forest <10 years     26761
# 27     1             Newly-Detected Vegetated Settlement <10 years    971718
# 28     1 Newly-Detected Very High Reflectance Settlement <10 years     89004
# 29     1                                                Other Land    714777

freq(native2020)
# layer value     count
# 1     1     0 430559025
# 2     1     1  88005602
# There are 7,920,504 ha of native grassland in the LUTS 2020 in the study area raster
# = 82,714 ha less than in 2010
freq(tame2020_luts)
# layer value     count
# 1     1     0 435011266
# 2     1     1  83553361
# There are 7,519,802 ha of tame grassland in the LUTS 2020 in the study area raster
freq(annual2020)
# layer value     count
# 1     1     0 268239566
# 2     1     1 250325061
freq(settlement2020)
# layer value     count
# 1     1     0 507577018
# 2     1     1  10987609


# LULC 2025
freq(lulc2025)
# layer value     count
# 1     1     1 100226817
# 2     1     2  36417806
# 3     1     3 300084392
# 4     1     4  16023110
# 5     1     5  65820330

# There are 9,020,414 ha of native grassland in the Grassland Consensus 2025 map in the study area raster
# There are 3,277,603 ha of tame grassland in the Grassland Consensus 2025 map in the study area raster


# Summarize transitions from native grassland using raw LUTS
# ignore transitions from tame/cropland because LUTS 2010 doesn't differentiate
# prep ecodistric - ecoregion lookup
eco_names <- as.data.frame(ecoregions_sa) %>%
  distinct(ecoregion = ECOREGION_, ecoregion_name = ECOREGION1)

lut <- as.data.frame(ecodistricts_sa) %>%
  distinct(stratum_id = as.numeric(ECODISTRIC), ecoregion = ECOREGION_) %>%
  left_join(eco_names, by = "ecoregion")

# apply transitions function
tr_luts <- tabulate_transitions(luts2010, luts2020, ecodistricts_r, codes = groups)
# add ecoregion
tr_luts <- left_join(tr_luts, lut, by = "stratum_id")
write.csv(tr, file.path("./processeddatasets", aafc_filepath, "transitions_luts_2010_2020_by_ecodistrict.csv"), row.names = FALSE)

# summarize native per ecoregion
native_by_ecoregion <- tr_luts %>%
  filter(!is.na(ecoregion)) %>%                       # drops the 61 edge pixels (stratum 0)
  group_by(ecoregion, ecoregion_name) %>%
  summarise(
    native_2010_ha = sum(area_ha[from == "native"]),
    native_2020_ha = sum(area_ha[to == "native"]),    # 2010 minus losses plus any gains into native
    to_crop_ha     = sum(area_ha[from == "native" & to == "crop"]),
    to_settle_ha   = sum(area_ha[from == "native" & to == "settlement"]),
    to_tame_ha     = sum(area_ha[from == "native" & to == "tame"]),
    .groups = "drop") %>%
  mutate(across(starts_with("to_"), ~ 100 * .x / native_2010_ha, .names = "{.col}_pct")) %>%
  arrange(desc(to_crop_ha))
# # A tibble: 7 × 10
# ecoregion ecoregion_name             native_2010_ha native_2020_ha to_crop_ha to_settle_ha to_tame_ha to_crop_ha_pct to_settle_ha_pct to_tame_ha_pct
# <int> <chr>                               <dbl>          <dbl>      <dbl>        <dbl>      <dbl>          <dbl>            <dbl>          <dbl>
# 1       156 Aspen Parkland                    815121.        784053.    9781.         9453.     11835.          1.20             1.16            1.45 
# 2       159 Mixed Grassland                  4820272.       4800505.    3936.         5356.     10476.          0.0816           0.111           0.217
# 3       157 Moist Mixed Grassland            1435344.       1410543.    3159.        10844.     10799.          0.220            0.755           0.752
# 4       158 Fescue Grassland                  340608.        336180.     852.         1399.      2177.          0.250            0.411           0.639
# 5       160 Cypress Upland                    550889.        548789.     219.          410.      1471.          0.0397           0.0744          0.267
# 6       162 Lake Manitoba Plain                38387.         37921.     128.          250.        88.5         0.334            0.652           0.230
# 7       163 Southwest Manitoba Uplands          2232           2158.       9.63         11.1       53.0         0.431            0.496           2.38 



# ---- ACI tame evidence ----
# Count of years in each window coded tame. Written to disk because these are
# the expensive steps and are reused by the sensitivity variants.

evidence_name <- function(yrs) sprintf("aci_tame_evidence_%s_%dm",
                                       paste(range(yrs), collapse = "-"), res_m)

tame_evidence_2010 <- app(aci_2010 %in% aci_tame_codes, fun = "sum", na.rm = TRUE,
                          filename = tagged_path("./processeddatasets", file.path(aci_filepath, "30m"),
                                                 evidence_name(aci_years_2010)),
                          overwrite = TRUE, wopt = list(datatype = "INT1U"))
names(tame_evidence_2010) <- "n_tame_years_2010"

tame_evidence_2020 <- app(aci_2020 %in% aci_tame_codes, fun = "sum", na.rm = TRUE,
                          filename = tagged_path("./processeddatasets", file.path(aci_filepath, "30m"),
                                                 evidence_name(aci_years_2020)),
                          overwrite = TRUE, wopt = list(datatype = "INT1U"))
names(tame_evidence_2020) <- "n_tame_years_2020"

# ---- LULC 2010 ----
# LUTS bounds the agricultural domain; ACI only splits it into tame vs annual.
# Pixels in the domain with forage evidence below the threshold fall through to
# annual cropland, which is the intended rotational-hay treatment.
tame2010     <- luts_ag_2010 & tame_evidence_2010 >= tame_min_years

lulc_2010 <- ifel(native2010 == 1, cls_native,
             ifel(tame2010, cls_tame,
             ifel(luts_ag_2010, cls_annual,
             ifel(settlement2010 == 1, cls_settlement, cls_other)))
)
names(lulc_2010) <- "lulc_2010"
levels(lulc_2010) <- lulc_levels # from 0a-setup.R

# ---- LULC 2020 ----
# Same construction. LUTS's own annual/tame split is
# deliberately not used, so both maps are built by the same rule.
tame2020     <- luts_ag_2020 & tame_evidence_2020 >= tame_min_years

lulc_2020 <- ifel(native2020 == 1, cls_native,
             ifel(tame2020, cls_tame,
             ifel(luts_ag_2020, cls_annual,
             ifel(settlement2020 == 1, cls_settlement, cls_other)))
)
names(lulc_2020) <- "lulc_2020"
levels(lulc_2020) <- levels(lulc_2010)

# ---- Derived layers ----
# grass_origin and grass2010 are both functions of lulc_2010
grass_origin <- ifel(lulc_2010 == cls_native, 1,
                     ifel(lulc_2010 == cls_tame, 2, NA))
names(grass_origin) <- "grass_origin"

grass2010 <- mask(lulc_2010 %in% c(cls_native, cls_tame), lulc_2010)
names(grass2010) <- "grass2010"

# ---- Write ----
write_rast_tagged(lulc_2010, "./processeddatasets", "landcover", "lulc_2010.tif", "30m", datatype = "INT1U")
write_rast_tagged(lulc_2020, "./processeddatasets", "landcover", "lulc_2020.tif", "30m_122", datatype = "INT1U")
write_rast_tagged(grass_origin, "./processeddatasets", "predictors/30m", "grassland_origin_2010.tif", "30m", datatype = "INT1U")
writeRaster(grass2010, file.path(binary_dir, "AAFC2010_ACI2009-2012_grass.tif"), datatype = "INT1U", overwrite = TRUE)

# lulc_2010    <- rast(tagged_path("./processeddatasets", "landcover", "lulc_2010.tif", "30m"))
# lulc_2020    <- rast(tagged_path("./processeddatasets", "landcover", "lulc_2020.tif", "30m"))
# grass_origin <- rast(tagged_path("./processeddatasets", "predictors/30m", "grassland_origin_2010.tif", "30m"))

# ---- QA/QC ----
freq(lulc_2010)
freq(lulc_2020)

# Expected against the previous 7-response.R run:
#   native  92,380,345    tame  80,956,335    grassland total 173,336,680
# Native should reproduce exactly (LUTS alone defines it). Tame may move
# slightly if the 2010 ACI files were regenerated by the updated 1h.
freq(grass2010)

# Full transition matrix. Script 7 reduces this to a binary; every response and
# sensitivity variant should be expressible as a selection on it.
tx <- crosstab(c(lulc_2010, lulc_2020), useNA = FALSE)
cell_ha <- (res_m^2) / 1e4
tx_ha <- tx * cell_ha
print(round(tx_ha))
write.csv(as.data.frame(tx), file.path(summary_dir, "lulc_transitions_2010_2020.csv"),
          row.names = FALSE)

# The three exclusion pools, read off tx. native->tame was ~418k under the
# LUTS-only definition. If either settlement pool is large relative to ~14M
# conversions, LUTS Settlement is absorbing more grassland than expected and
# the farmyard question gets urgent.
#   native -> tame        (tame forage)
#   native -> settlement
#   tame   -> settlement

# Farm infrastructure gap: LUTS settlement land that ACI still calls grass or
# forage. Split by 2010 class to separate the coverage hole (settlement in both
# years) from grassland lost to settlement.
aci_grass_2020 <- app(aci_2020 %in% c(110, 122), fun = "sum", na.rm = TRUE) >= tame_min_years
crosstab(c(lulc_2010, ifel(lulc_2020 == cls_settlement & aci_grass_2020, 1, 0)),
         useNA = FALSE)

# ---- Record definition provenance ----
lulc_definition <- data.frame(
  parameter = c("aci_years_2010", "aci_years_2020", "aci_tame_codes",
                "tame_min_years",
                "native_source", "settlement_source", "luts_domain_2010",
                "luts_domain_2020", "tame_rule", "annual_rule", "resolution_m",
                "n_native_2010", "n_tame_2010", "n_annual_2010", "n_settlement_2010",
                "n_native_2020", "n_tame_2020", "n_annual_2020", "n_settlement_2020",
                "date_run", "script"),
  value = c(
    paste(range(aci_years_2010), collapse = "-"),
    paste(range(aci_years_2020), collapse = "-"),
    paste(aci_tame_codes, collapse = ";"),
    tame_min_years,
    "AAFC LUTS native grassland classes, both years",
    "AAFC LUTS settlement/settlement classes, both years",
    "LUTS 2010 Cropland (annual + tame, undivided)",
    "LUTS 2020 Annual Cropland UNION tame Cropland",
    sprintf("LUTS domain AND tame-coded in >= %d years of window", tame_min_years),
    "LUTS domain AND NOT tame (residual within the domain)",
    res_m,
    freq(lulc_2010)$count[freq(lulc_2010)$value == cls_native],
    freq(lulc_2010)$count[freq(lulc_2010)$value == cls_tame],
    freq(lulc_2010)$count[freq(lulc_2010)$value == cls_annual],
    freq(lulc_2010)$count[freq(lulc_2010)$value == cls_settlement],
    freq(lulc_2020)$count[freq(lulc_2020)$value == cls_native],
    freq(lulc_2020)$count[freq(lulc_2020)$value == cls_tame],
    freq(lulc_2020)$count[freq(lulc_2020)$value == cls_annual],
    freq(lulc_2020)$count[freq(lulc_2020)$value == cls_settlement],
    as.character(Sys.Date()),
    "3b-lulc-maps.R"
  ),
  note = c(
    "Baseline evidence window; brackets 2010",
    "Outcome evidence window; disjoint from 2010 by construction",
    "ACI Pasture/forages only",
    "Adds Grassland as a drift correction for the 122 decline; see header",
    "Shared by both timesteps; sensitivity: 3",
    "LUTS is authoritative on native at both dates",
    "Includes farmyards; see KNOWN GAP in header",
    "The class ACI splits at 2010",
    "The class ACI splits at 2020; LUTS's own split is not used",
    "Applied identically at 2010 and 2020",
    "No ACI crop code is ever tested",
    "Native predictor/response grid",
    "", "", "", "", "", "", "", "",
    "", ""
  )
)

write.csv(lulc_definition, file.path(summary_dir, "lulc_definition.csv"), row.names = FALSE)