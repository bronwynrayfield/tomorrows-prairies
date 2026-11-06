# =============================================================================
# Purpose:
#   Builds the LULC 2010 and 2020 state maps (classes from 0a-setup.R) on the
#   modelling extent, masks them to the study area for the response, and writes
#   the ACI count rasters that every map variant thresholds.
#   No transition is scored here; that is 7a and 7b-response-sensitivity.
#
# LUTS defines native, the agricultural domain, settlement and other at both
# dates. ACI is used to split the 2010 agricultural domain into tame vs
# annual, because LUTS 2010 has no tame class. At 2020 LUTS's own
# annual/perennial split is used, except that a switch between tame and annual
# counts only if ACI shows a clean switch (clean0_k2; see LULC 2020 below).
#
# Notes:
#   - Tame 2010 counts evidence, not persistence: a pixel tame in 2009-2010 and
#     broken in 2011 still qualifies. Early converters are the most informative
#     positives.
#   - ACI 2009-2010 are 56 m resampled to 30 m; 2019-2022 are native 30 m.
#   - LUTS checks use studyarea, the build uses modelling-extent
# 
# Outputs:
#   - lulc_{2010,2020}.tif, tag _modelling-extent  
#   - lulc_{2010,2020}.tif,                        
#   - aci_{110,122}_count_{2009-2012,2019-2022}_30m   
#   - grassland_origin_2010_30m.tif, AAFC2010_ACI2009-2012_grass.tif 
#   - lulc_transitions_2010_2020.csv, lulc_definition.csv           
# =============================================================================

source(file.path("scripts/landcover", "0a-setup.R"))

# ---- Parameters ----
# See 0a-setup.R

# ---- Load ----
# Boundaries (study area)
sa_r <- rast(file.path(sa_filepath, sa_raster_filename))
ecoregions_sa  <- st_read(tagged_path("./processeddatasets", ecoregions_filepath, ecoregions_filename, "studyarea"))
ecoregions_r   <- rast(tagged_path("./processeddatasets", ecoregions_filepath, paste0(tools::file_path_sans_ext(ecoregions_filename), ".tif"), "studyarea"))
ecodistricts_sa  <- st_read(tagged_path("./processeddatasets", ecodistricts_filepath, ecodistricts_filename, "studyarea"))
ecodistricts_r <- rast(tagged_path("./processeddatasets", ecodistricts_filepath, paste0(tools::file_path_sans_ext(ecodistricts_filename), ".tif"), "studyarea"))

# LULC (study area) - used by the LUTS checks only
lulc2025       <- rast(tagged_path("./processeddatasets", grass_consensus_filepath, grass_consensus_filename, "study-area"))
luts2010       <- rast(tagged_path(file.path("./processeddatasets", aafc_filepath), "cropped", aafc2010_filename, "study-area"))
luts2020       <- rast(tagged_path(file.path("./processeddatasets", aafc_filepath), "cropped", aafc2020_filename, "study-area"))

# LUTS binaries (study area)
sa_name <- function(f) sub("\\.tif$", "_studyarea.tif", f)
native2010     <- rast(file.path("./processeddatasets", binary_filepath, sa_name(luts_filenames["native2010"])))
tamecrop2010   <- rast(file.path("./processeddatasets", binary_filepath, sa_name(luts_filenames["tamecrop2010"])))
settlement2010 <- rast(file.path("./processeddatasets", binary_filepath, sa_name(luts_filenames["settlement2010"])))
native2020     <- rast(file.path("./processeddatasets", binary_filepath, sa_name(luts_filenames["native2020"])))
annual2020     <- rast(file.path("./processeddatasets", binary_filepath, sa_name(luts_filenames["annual2020"])))
tame2020_luts  <- rast(file.path("./processeddatasets", binary_filepath, sa_name(luts_filenames["tame2020"])))
settlement2020 <- rast(file.path("./processeddatasets", binary_filepath, sa_name(luts_filenames["settlement2020"])))

# LUTS binaries (modelling extent) - used by the build
luts_me <- lapply(luts_filenames, \(f) rast(file.path("./processeddatasets", binary_filepath, f)))

# Define total cropland in 2010 (tame and annual crop) 
luts_ag_2010 <- luts_me$tamecrop2010 == 1

# ACI
aci_2010 <- rast(file.path("./processeddatasets", aci_filepath, "30m", sprintf("aci_%d_30m_modelling-extent.tif", aci_years_2010)))
names(aci_2010) <- paste0("y", aci_years_2010)

aci_2020 <- rast(file.path("./processeddatasets", aci_filepath, "30m", sprintf("aci_%d_30m_modelling-extent.tif", aci_years_2020)))
names(aci_2020) <- paste0("y", aci_years_2020)


# ---- Check grassland extents and frequencies in raw LUTS ----
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

# 2010: LUTS has no tame/annual split, so 4class (crop incl. tame); 2020: 5class
l10 <- reclass_aafc(luts2010, groups = groups, scheme = "4class")
l20 <- reclass_aafc(luts2020, groups = groups, scheme = "5class")
levels(l10) <- setNames(data.frame(c(1, 3, 4, 5),
                                   c("native grassland", "cropland (tame + annual)", "settlement", "other")),
                        names(lulc_levels))
levels(l20) <- lulc_levels
tx_luts5 <- crosstab(c(l10, l20), useNA = FALSE) * cell_ha
print(tx_luts5)
#                          class.1
# class                    annual cropland native grassland       other  settlement tame grassland
# cropland (tame + annual)     22458626.64             0.00      509.49    87445.53     7439650.47
# native grassland                18085.14       7920504.18        0.00    27725.94       36902.70
# other                           52543.71             0.00  7711859.97     7153.92       43249.32
# settlement                          0.00             0.00        0.00   866559.42           0.00


# LULC 2025
freq(lulc2025)
# layer value     count
# 1     1     1  95773537
# 2     1     2  57628305
# 3     1     3 291638177
# 4     1     4  14422618
# 5     1     5  59280901

# There are 8,619,618 ha of native grassland in the Grassland Consensus 2025 map in the study area raster
# There are 5,186,547 ha of tame grassland in the Grassland Consensus 2025 map in the study area raster


# Summarize transitions from native grassland using raw LUTS
# ignore transitions from tame/cropland because LUTS 2010 doesn't differentiate
# prep ecodistric - ecoregion lookup
eco_names <- as.data.frame(ecoregions_sa) %>%
  distinct(ecoregion = ECOREGION_, ecoregion_name = ECOREGION1)

lut <- as.data.frame(ecodistricts_sa) %>%
  distinct(stratum_id = as.numeric(ECODISTRIC), ecoregion = ECOREGION_) %>%
  left_join(eco_names, by = "ecoregion")

# apply transitions function
# Crop to LUTS grid: extra extent all NA; 356,902 ecodistrict cells (~321 km², US border) have no LUTS data and drop out of tabulation
ecodistricts_r <- crop(ecodistricts_r, luts2010)
tr_luts <- tabulate_transitions(luts2010, luts2020, ecodistricts_r, codes = groups)
# add ecoregion
tr_luts <- left_join(tr_luts, lut, by = "stratum_id")
write.csv(tr_luts, file.path("./processeddatasets", aafc_filepath, "transitions_luts_2010_2020_by_ecodistrict.csv"), row.names = FALSE)

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
native_by_ecoregion
# # A tibble: 7 × 10
#   ecoregion ecoregion_name             native_2010_ha native_2020_ha to_crop_ha to_settle_ha to_tame_ha to_crop_ha_pct to_settle_ha_pct to_tame_ha_pct
#       <int> <chr>                               <dbl>          <dbl>      <dbl>        <dbl>      <dbl>          <dbl>            <dbl>          <dbl>
# 1       156 Aspen Parkland                    815122.        784053.    9781.         9453.     11835.          1.20             1.16            1.45 
# 2       159 Mixed Grassland                  4820271.       4800503.    3936.         5356.     10476.          0.0816           0.111           0.217
# 3       157 Moist Mixed Grassland            1435346.       1410546.    3159.        10844.     10799.          0.220            0.755           0.752
# 4       158 Fescue Grassland                  340608.        336180.     852.         1399.      2177.          0.250            0.411           0.639
# 5       160 Cypress Upland                    550890.        548790.     219.          410.      1471.          0.0397           0.0744          0.267
# 6       162 Lake Manitoba Plain                38387.         37921.     128.          250.        88.5         0.334            0.652           0.230
# 7       163 Southwest Manitoba Uplands          2232           2158.       9.63         11.1       53.0         0.431            0.496           2.38 



# ---- ACI tame evidence ----
# Years coded 110 and years coded 122, per window, written to disk because
# every response variant for sensitivity is a threshold on them. 
# Each pixel-year carries one ACI code, so count110 + count122 = years coded 110 or 122.

aci_counts <- function(aci, yrs) {
  r <- rast(lapply(aci_count_codes, \(cd)
                   app(aci == cd, fun = "sum", na.rm = TRUE,
                       filename = tagged_path("./processeddatasets", file.path(aci_filepath, "30m"),
                                              count_name(cd, yrs)),
                       overwrite = TRUE, wopt = list(datatype = "INT1U"))))
  names(r) <- paste0("aci", aci_count_codes)
  r
}

# Generate counts for both ACI grass classes and years (some used later in response variants)
count_2010 <- aci_counts(aci_2010, aci_years_2010)
count_2020 <- aci_counts(aci_2020, aci_years_2020)


# ---- LULC 2010 (modelling extent) ----
# LUTS bounds the agricultural domain; ACI only splits it into tame vs annual.
# Pixels in the domain with forage evidence below the threshold are classified as annual cropland
tame2010 <- luts_ag_2010 & tame_ev(count_2010, aci_tame_code, tame_min_years)

lulc_2010 <- build_lulc(luts_me$native2010, tame2010, luts_ag_2010, luts_me$settlement2010)
names(lulc_2010) <- "lulc_2010"


# ---- LULC 2020 (modelling extent) ----
# LUTS's own split: perennial cropland = tame, annual cropland = annual.
lulc_2020_luts <- build_lulc(luts_me$native2020, luts_me$tame2020, luts_me$annual2020, luts_me$settlement2020)
names(lulc_2020_luts) <- "lulc_2020"

# Hindcast map (clean0_k2): LUTS, except a tame<->annual switch counts only on a clean
# ACI switch (>= tame_min_years forage years at one date, none at the other). Chosen
# from 7a rules against the Census of Agriculture (7d): right sign in 24/24 CDs, r = 0.52;
# net tame -7.4%, annual +2.4% (study area). Native outcomes stay LUTS.
ev10 <- subst(count_2010$aci110 + count_2010$aci122, NA, 0)
ev20 <- subst(count_2020$aci110 + count_2020$aci122, NA, 0)
lulc_2020 <- ifel(lulc_2010 == cls_tame & lulc_2020_luts == cls_annual & ev20 > 0, cls_tame, lulc_2020_luts)
lulc_2020 <- ifel(lulc_2010 == cls_annual & lulc_2020_luts == cls_tame & !(ev10 == 0 & ev20 >= tame_min_years),
                  cls_annual, lulc_2020)
levels(lulc_2020) <- levels(lulc_2020_luts)
names(lulc_2020) <- "lulc_2020"


# ---- Derived layers (study area) ----
lulc_2010_sa <- mask(crop(lulc_2010, sa_r), sa_r)
lulc_2020_sa <- mask(crop(lulc_2020, sa_r), sa_r)
lulc_2020_luts_sa <- mask(crop(lulc_2020_luts, sa_r), sa_r)

# grass_origin and grass2010 are both functions of lulc_2010
grass_origin <- ifel(lulc_2010_sa == cls_native, 1,
                     ifel(lulc_2010_sa == cls_tame, 2, NA))
names(grass_origin) <- "grass_origin"

grass2010 <- mask(lulc_2010_sa %in% c(cls_native, cls_tame), lulc_2010_sa)
names(grass2010) <- "grass2010"

# ---- Write ----
write_rast_tagged(lulc_2010, "./processeddatasets", "landcover", "lulc_2010.tif", "modelling-extent", datatype = "INT1U")
write_rast_tagged(lulc_2020, "./processeddatasets", "landcover", "lulc_2020.tif", "modelling-extent", datatype = "INT1U")
write_rast_tagged(lulc_2010_sa, "./processeddatasets", "landcover", "lulc_2010.tif", "study-area", datatype = "INT1U")
write_rast_tagged(lulc_2020_sa, "./processeddatasets", "landcover", "lulc_2020.tif", "study-area", datatype = "INT1U")
write_rast_tagged(grass_origin, "./processeddatasets", file.path(predictors_filepath, "30m"), "grassland_origin_2010.tif", "study-area", datatype = "INT1U")
writeRaster(grass2010, file.path("./processeddatasets", binary_filepath, "AAFC2010_ACI2009-2012_grass.tif"), datatype = "INT1U", overwrite = TRUE)
write_rast_tagged(lulc_2020_luts, "./processeddatasets", "landcover", "lulc_2020_luts.tif", "modelling-extent", datatype = "INT1U")
write_rast_tagged(lulc_2020_luts_sa, "./processeddatasets", "landcover", "lulc_2020_luts.tif", "study-area", datatype = "INT1U")

# ---- QA/QC (study area)----
f10 <- freq(lulc_2010_sa)
f10
# layer            value     count
# 1     1 native grassland  88924644
# 2     1   tame grassland  74688926
# 3     1  annual cropland 258491431
# 4     1       settlement   9628438
# 5     1            other  86831188
f20 <- freq(lulc_2020_luts_sa)
f20
# layer            value     count
# 1     1 native grassland  88005602
# 2     1   tame grassland  83553361
# 3     1  annual cropland 250325061
# 4     1       settlement  10987609
# 5     1            other  85692994
f20h <- freq(lulc_2020_sa)        # hindcast map
f20h
# layer            value     count
# 1     1 native grassland  88005602
# 2     1   tame grassland  69126999
# 3     1  annual cropland 264751423
# 4     1       settlement  10987609
# 5     1            other  85692994

#n_cls <- function(f, cls) f$count[f$value == cls]
n_cls <- function(f, cls) f$count[f$value == lulc_levels[[2]][lulc_levels[[1]] == cls]]
n_bin <- function(r) global(r, "sum", na.rm = TRUE)[[1]]


# LUTS binaries are mutually exclusive, so every LUTS-defined class must
# reproduce the study-area binaries exactly. This also ties the
# modelling-extent binaries (build) to the study-area binaries (checks).
stopifnot(
  "native 2010"     = n_cls(f10, cls_native)     == n_bin(native2010),
  "settlement 2010" = n_cls(f10, cls_settlement) == n_bin(settlement2010),
  "ag domain 2010"  = n_cls(f10, cls_tame) + n_cls(f10, cls_annual) == n_bin(tamecrop2010),
  "native 2020"     = n_cls(f20, cls_native)     == n_bin(native2020),
  "tame 2020"       = n_cls(f20, cls_tame)       == n_bin(tame2020_luts),
  "annual 2020"     = n_cls(f20, cls_annual)     == n_bin(annual2020),
  "settlement 2020" = n_cls(f20, cls_settlement) == n_bin(settlement2020),
  "ag domain 2020 (hindcast)" = n_cls(f20h, cls_tame) + n_cls(f20h, cls_annual) == n_bin(tame2020_luts) + n_bin(annual2020),
  "native 2020 (hindcast)"    = n_cls(f20h, cls_native) == n_bin(native2020),
  "settlement 2020 (hindcast)" = n_cls(f20h, cls_settlement) == n_bin(settlement2020)
)

# Cheeck against sum of native and tame in 2010
freq(grass2010)
# layer value     count
# 1     1     0 354951057
# 2     1     1 163613570

# Full transition matrix
tx <- crosstab(c(lulc_2010_sa, lulc_2020_sa), useNA = FALSE)
tx_ha <- tx * cell_ha
print(round(tx_ha))
#                    lulc_2020
# lulc_2010        annual cropland native grassland    other settlement tame grassland
# annual cropland         22831805                0      235      51431         380758
# native grassland           18085          7920504        0      27726          36903
# other                      52544                0  7711860       7154          43249
# settlement                     0                0        0     866559              0
# tame grassland            925194                0      274      36015        5760520
write.csv(as.data.frame(tx), file.path("./processeddatasets", "landcover", "lulc_transitions_2010_2020_study-area.csv"),
          row.names = FALSE)
