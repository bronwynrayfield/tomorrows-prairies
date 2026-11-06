# =============================================================================
# Purpose:
#   Sensitivity of the five responses to the tame definition. Each variant is
#   a (2010 rule, 2020 rule) pair thresholding the ACI count rasters from 3b.
#   Reference is run through the same code and must reproduce 7 exactly.
#   Variants and rationale: tame-definition-variants.xlsx
#
#   2010 variants move pixels between tame and annual inside the fixed LUTS ag
#   domain, so they change the tame and annual focal proportions. Tame
#   proportion is recomputed per 2010 variant on the modelling extent;
#   annual = ag - tame, with ag taken from the reference 4d layers.
#   2020-only variants leave every 2010 predictor unchanged.
#
#   Only responses a variant changes are written: a 2010 change affects the
#   tame-origin responses; a 2020 change affects all but native -> settlement
#   (native persistence and settlement are LUTS-defined at both dates).
#
#   Reference: 2010 tame = ACI 110+122 >= 2 years (2009-2012) within LUTS Cropland;
#   2020 tame = LUTS perennial. Variants: aci122 (the previous reference, so the
#   code-set switch is itself reported), min3 (threshold), aci_both (ACI 110+122
#   at 2020 too: the 2020 source), min3_both(k = 3 at both dates). Decision and evidence: "Tame grassland
#   definition: decision record". aci_out, aci_out_drift and aci_out_min3 were
#   dropped: they mix code sets between dates.
# =============================================================================

source(file.path("scripts/landcover", "0a-setup.R"))

prop_dir     <- file.path("./processeddatasets", "predictors/30m")
f_prop_tame   <- "SET_TO_4D_TAME_PROPORTION.tif"     # reference 4d layers, modelling extent
f_prop_annual <- "SET_TO_4D_ANNUAL_PROPORTION.tif"

# ---- Variants ----
# codes: ACI count layers summed as tame evidence. codes_2020 = NA: LUTS split.
variants <- tibble::tribble(
  ~variant,    ~codes_2010,   ~min_2010,          ~codes_2020,   ~min_2020,
  "reference", aci_tame_code, tame_min_years,     NA,            NA,
  "aci122",    "122",         tame_min_years,     NA,            NA,
  "min3",      aci_tame_code, tame_min_years + 1, NA,            NA,
  "aci_both",  aci_tame_code, tame_min_years,     aci_tame_code, tame_min_years,
  "min3_both", aci_tame_code, tame_min_years + 1, aci_tame_code, tame_min_years + 1
) |>
  mutate(chg_2010 = codes_2010 != codes_2010[1] | min_2010 != min_2010[1],
         chg_2020 = !is.na(codes_2020))

# ---- Load ----
sa_r           <- rast(file.path(sa_filepath, sa_raster_filename))
ecodistricts_r <- crop(rast(tagged_path("./processeddatasets", ecodistricts_filepath, paste0(tools::file_path_sans_ext(ecodistricts_filename), ".tif"), "studyarea")), sa_r)

lulc_2010_sa <- rast(tagged_path("./processeddatasets", "landcover", "lulc_2010.tif", "study-area"))
lulc_2020_sa <- rast(tagged_path("./processeddatasets", "landcover", "lulc_2020.tif", "study-area"))

load_counts <- function(yrs) {
  r <- rast(sapply(aci_count_codes, \(cd)
                   tagged_path("./processeddatasets", file.path(aci_filepath, "30m"), count_name(cd, yrs))))
  names(r) <- paste0("aci", aci_count_codes)
  r
}
count_2010_me <- load_counts(aci_years_2010)          # 2010 builds run on the modelling extent (focal)
count_2020_sa <- mask(crop(load_counts(aci_years_2020),sa_r), sa_r)   # 2020 builds only feed responses

# 2010 on the modelling extent: variant maps feed the tame focal proportion,
# which needs the buffer so windows at the study-area edge see real neighbours.
me_2010 <- lapply(luts_filenames[c("native2010", "tamecrop2010", "settlement2010")],
                  \(f) rast(file.path("./processeddatasets", binary_filepath, f)))
# 2020 on the study area only: variant maps feed the responses alone (all
# training predictors are 2010-state)
sa_name <- function(f) sub("\\.tif$", "_studyarea.tif", f)
sa_2020 <- lapply(luts_filenames[c("native2020", "annual2020", "tame2020", "settlement2020")],
                  \(f) rast(file.path("./processeddatasets", binary_filepath, sa_name(f))))

# LUTS ag domains, fixed across variants: each variant only splits them into tame vs annual.
# 2010: undifferentiated Cropland (modelling extent, for the focal).
luts_ag_2010_me <- me_2010$tamecrop2010 == 1
# 2020: annual U perennial. Under an ACI split a LUTS-annual pixel can be tame
# and a LUTS-perennial pixel can be annual, so the full union is the domain.
luts_ag_2020_sa <- sa_2020$annual2020 == 1 | sa_2020$tame2020 == 1
# TBD Focal ag proportion (tame + annual) from 4d's reference layers. The ag domain
# doesn't vary, so neither does this; variants recompute tame and get annual = ag - tame.
# ag_prop_me <- rast(file.path(prop_dir, f_prop_tame)) + rast(file.path(prop_dir, f_prop_annual))

# ---- Relabelling check (decides aci_out_drift) ----
# Stable tame: reference tame in 2010 that LUTS still calls perennial in 2020
# (51,770,979 px). How does ACI code it in 2019-2022?
# value: 0 = neither, 1 = 122 only, 10 = 110 only, 11 = both (each >= min_ref years)
# Result: 122 only 47.6%, both 6.4%, 110 only 21.7%, neither 24.3%.
# 122 alone recognises 54.0%; the summed aci_out_drift rule (110 + 122) 78.4%.
# 21.7% of stable forage is coded 110 rather than 122 in 2019-2022 (relabelling
# or forage falling out of management; either way still grass), so aci_out_drift is kept. The 21.6% that
# fails even the summed rule is ACI/LUTS disagreement no code choice fixes,
# which is why LUTS is the 2020 reference.
min_ref <- variants$min_2010[1]
# Population fixed to the 122 rule: the check asks whether 122 forage is
# relabelled 110, independent of the reference code set.
persist <- mask(crop(luts_ag_2010_me & tame_ev(count_2010_me, "122", min_ref), sa_r), sa_r) &
  lulc_2020_sa == cls_tame
relabel <- freq(ifel(persist, (count_2020_sa$aci110 >= min_ref) * 10 +
                       (count_2020_sa$aci122 >= min_ref), NA))
relabel$share <- relabel$count / sum(relabel$count)
print(relabel)
# layer value    count      share
# 1     1     0 12567432 0.24275052
# 2     1     1 24618535 0.47552771
# 3     1    10 11254132 0.21738302
# 4     1    11  3330880 0.06433875
freq(ifel(persist, (count_2020_sa$aci110 + count_2020_sa$aci122) >= min_ref, NA))
# layer value    count                    
# 1     1     0 11184973
# 2     1     1 40586006

# ---- Threshold calibration against LUTS 2020 ----
# Within the 2020 LUTS ag domain, which ACI evidence threshold best reproduces
# LUTS perennial? The best k is the evidence-based threshold for 2010, so both
# dates measure the same kind of land. Tested for 122 alone (the 2010 rule) and
# 122 only:  kappa 0.580 / 0.484 / 0.364 / 0.208 at k = 1-4; omission 57% at k = 2.
# 110 + 122: kappa 0.704 / 0.691 / 0.629 / 0.490 at k = 1-4; omission 32% at k = 2.
# k >= 3 is worse under both code sets; 122 alone misses most of LUTS perennial at 2020.
calib_one <- function(codes) {
   tab <- freq(ifel(luts_ag_2020_sa,
                   subst(sum(count_2020_sa[[paste0("aci", strsplit(codes, "\\+")[[1]])]]), NA, 0) * 10 +
                     (sa_2020$tame2020 == 1), NA)) |>
    mutate(n_years = value %/% 10, luts_tame = value %% 10)
  bind_rows(lapply(seq_along(aci_years_2020), \(k) summarise(tab, codes = codes, k = k,
                                                             tp = sum(count[n_years >= k & luts_tame == 1]),
                                                             fp = sum(count[n_years >= k & luts_tame == 0]),
                                                             fn = sum(count[n_years <  k & luts_tame == 1]),
                                                             tn = sum(count[n_years <  k & luts_tame == 0]))))
}

calib_k <- bind_rows(lapply(c("122", "110+122"), calib_one)) |>
  mutate(n = tp + fp + fn + tn,
         accuracy   = (tp + tn) / n,
         omission   = fn / (tp + fn),
         commission = fp / (tp + fp),
         pe    = ((tp + fp) * (tp + fn) + (fn + tn) * (fp + tn)) / n^2,
         kappa = (accuracy - pe) / (1 - pe))
print(calib_k)
#     codes k       tp       fp       fn        tn         n  accuracy  omission commission        pe     kappa
# 1     122 1 48684646 13236353 34868715 237088708 333878422 0.8559204 0.4173227 0.21376194 0.6571123 0.5798053
# 2     122 2 36105694  6730190 47447667 243594871 333878422 0.8377318 0.5678726 0.15711570 0.6856646 0.4837737
# 3     122 3 25012996  3732663 58540365 246592398 333878422 0.8134859 0.7006345 0.12985136 0.7067442 0.3639883
# 4     122 4 13120473  1613122 70432888 248711939 333878422 0.7842148 0.8429689 0.10948597 0.7277070 0.2075259
# 5 110+122 1 64478199 17819548 19075162 232505513 333878422 0.8894966 0.2282992 0.21652534 0.6266277 0.7040395
# 6 110+122 2 56579175  8997223 26974186 241327838 333878422 0.8922620 0.3228378 0.13720215 0.6516437 0.6907246
# 7 110+122 3 47910521  5140920 35642840 245184141 333878422 0.8778485 0.4265877 0.09690444 0.6703816 0.6294153
# 8 110+122 4 34012910  2340473 49540451 247984588 333878422 0.8446113 0.5929199 0.06438116 0.6953627 0.4899223
write.csv(calib_k, file.path(response_summary_filepath, "sensitivity", "threshold_calibration_2020.csv"), row.names = FALSE)


# ---- Variant loop ----
tabs <- list()
prev <- read.csv(file.path(response_summary_filepath, "sensitivity", "transitions_by_variant.csv"))
todo <- filter(variants, !variant %in% unique(prev$variant))
for (i in seq_len(nrow(todo))) {
  v <- todo[i, ]
  
  # 2010: rebuild on the modelling extent, recompute tame focal, mask to study area
  if (v$chg_2010) {
    tame_me <- luts_ag_2010_me & tame_ev(count_2010_me, v$codes_2010, v$min_2010)
    l10_me  <- build_lulc(me_2010$native2010, tame_me, luts_ag_2010_me, me_2010$settlement2010)
    #tame_prop <- focal_prop(l10_me == cls_tame)
    #props <- c(tame_prop, ag_prop_me - tame_prop)
    #names(props) <- c("tame_prop", "annual_prop")
    #write_rast_tagged(to_sa(props), "./processeddatasets", "predictors/30m", "focal_tame_annual_2010.tif",
    #                  paste0("study-area_", v$variant))
    l10 <- mask(crop(l10_me, sa_r), sa_r)
  } else {
    l10 <- lulc_2010_sa
  }
  
  # 2020: study area only
  if (v$chg_2020) {
    tame20 <- luts_ag_2020_sa & tame_ev(count_2020_sa, v$codes_2020, v$min_2020)
    l20 <- build_lulc(sa_2020$native2020, tame20, luts_ag_2020_sa, sa_2020$settlement2020)
  } else {
    l20 <- lulc_2020_sa
  }
  
  # Responses this variant changes
  affected <- responses |>
    filter((v$chg_2010 & from == cls_tame) |
             (v$chg_2020 & !(from == cls_native & to == cls_settlement)))
  tcode <- l10 * 10 + l20
  for (r in affected$response) {
    s <- filter(scores, response == r)
    out <- classify(tcode, cbind(s$from * 10 + s$to, s$score), others = NA)
    names(out) <- r
    write_rast_tagged(out, "./processeddatasets", "response", paste0("response_", r, ".tif"),
                      paste0("study-area_", v$variant), datatype = "INT1U")
  }
  
  tabs[[v$variant]] <- tabulate_lulc(l10, l20, ecodistricts_r) |> mutate(variant = v$variant)
}

tab <- bind_rows(prev, bind_rows(tabs))
write.csv(tab, file.path(response_summary_filepath, "sensitivity", "transitions_by_variant.csv"), row.names = FALSE)

# ---- Extent ----
extent <- bind_rows(
  tab |> group_by(variant, class = from) |> summarise(n = sum(n), .groups = "drop") |> mutate(year = 2010),
  tab |> group_by(variant, class = to)   |> summarise(n = sum(n), .groups = "drop") |> mutate(year = 2020)) |>
  mutate(ha = n * cell_ha)
write.csv(extent, file.path(response_summary_filepath, "sensitivity", "extent_by_variant.csv"), row.names = FALSE)

# LUTS-defined totals (native, settlement, other, ag domain) must not move
inv <- extent |>
  mutate(class = ifelse(class %in% c(cls_tame, cls_annual), "ag", as.character(class))) |>
  group_by(variant, year, class) |> summarise(n = sum(n), .groups = "drop") |>
  group_by(year, class) |> summarise(same = n_distinct(n) == 1, .groups = "drop")
stopifnot("LUTS-defined totals differ across variants" = all(inv$same))

# ---- By response ----
by_response <- tab |>
  inner_join(scores, by = c("from", "to"), relationship = "many-to-many") |>
  group_by(variant, response) |>
  sum_scores()

ref_variant <- variants$variant[1]
by_response <- by_response |>
  left_join(by_response |> filter(variant == ref_variant) |>
              select(response, n_pos_ref = n_pos, rate_ref = rate), by = "response") |>
  mutate(d_pos = n_pos - n_pos_ref, d_rate = rate - rate_ref)
print(by_response, n = Inf)
# A tibble: 20 × 12
# variant   response             n_pos    n_neg n_competing   pos_ha    rate competing_per_pos n_pos_ref rate_ref    d_pos    d_rate
# <chr>     <chr>                <dbl>    <dbl>       <dbl>    <dbl>   <dbl>             <dbl>     <dbl>    <dbl>    <dbl>     <dbl>
# 1 aci122    native_annual       200946 88005588      718096   18085. 0.00228            3.57      200946  0.00228        0  0       
# 2 aci122    native_settlement   308066 88005588      610976   27726. 0.00349            1.98      308066  0.00349        0  0       
# 3 aci122    native_tame         410030 88005588      509012   36903. 0.00464            1.24      410030  0.00464        0  0       
# 4 aci122    tame_annual       11022980 51770972      318806  992068. 0.176              0.0289  14225857  0.192   -3202877 -0.0160  
# 5 aci122    tame_settlement     316133 51770972    11025653   28452. 0.00607           34.9       400162  0.00662   -84029 -0.000549
# 6 aci_both  native_annual       325105 88005588      593937   29259. 0.00368            1.83      200946  0.00228   124159  0.00140 
# 7 aci_both  native_settlement   308066 88005588      610976   27726. 0.00349            1.98      308066  0.00349        0  0       
# 8 aci_both  native_tame         285871 88005588      633171   25728. 0.00324            2.21      410030  0.00464  -124159 -0.00140 
# 9 aci_both  tame_annual       25133347 49152333      403211 2262001. 0.338              0.0160  14225857  0.192   10907490  0.147   
# 10 aci_both  tame_settlement     400162 49152333    25136396   36015. 0.00808           62.8       400162  0.00662        0  0.00146 
# 11 min3      native_annual       200946 88005588      718096   18085. 0.00228            3.57      200946  0.00228        0  0       
# 12 min3      native_settlement   308066 88005588      610976   27726. 0.00349            1.98      308066  0.00349        0  0       
# 13 min3      native_tame         410030 88005588      509012   36903. 0.00464            1.24      410030  0.00464        0  0       
# 14 min3      tame_annual        6634370 47750646      228502  597093. 0.122              0.0344  14225857  0.192   -7591487 -0.0695  
# 15 min3      tame_settlement     226598 47750646     6636274   20394. 0.00472           29.3       400162  0.00662  -173564 -0.00190 
# 16 reference native_annual       200946 88005588      718096   18085. 0.00228            3.57      200946  0.00228        0  0       
# 17 reference native_settlement   308066 88005588      610976   27726. 0.00349            1.98      308066  0.00349        0  0       
# 18 reference native_tame         410030 88005588      509012   36903. 0.00464            1.24      410030  0.00464        0  0       
# 19 reference tame_annual       14225857 60059823      403211 1280327. 0.192              0.0283  14225857  0.192          0  0       
# 20 reference tame_settlement     400162 60059823    14228906   36015. 0.00662           35.6       400162  0.00662        0  0   
write.csv(by_response, file.path(response_summary_filepath, "sensitivity", "conversion_by_response_by_variant.csv"), row.names = FALSE)

# Reference must reproduce response created in 7-response.R
reference <- read.csv(file.path(response_summary_filepath, "conversion_by_response.csv")) |>
  select(response, ref_pos = n_pos, ref_neg = n_neg) |>
  left_join(filter(by_response, variant == ref_variant), by = "response")
stopifnot("reference does not reproduce 7" =
            all(reference$n_pos == reference$ref_pos & reference$n_neg == reference$ref_neg))