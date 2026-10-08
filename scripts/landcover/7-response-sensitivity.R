# Purpose:
#   Sensitivity of the five responses to the tame definition. Each variant is
#   a (2010 rule, 2020 rule) pair thresholding the ACI count rasters from 3b.
#   Canonical is run through the same code and must reproduce 7 exactly.
#   Variants and rationale: tame-definition-variants.xlsx
#
#   2010 variants move pixels between tame and annual inside the fixed LUTS ag
#   domain, so they change the tame and annual focal proportions. Tame
#   proportion is recomputed per 2010 variant on the modelling extent;
#   annual = ag - tame, with ag taken from the canonical 4d layers.
#   2020-only variants leave every 2010 predictor unchanged.
#
#   Only responses a variant changes are written: a 2010 change affects the
#   tame-origin responses; a 2020 change affects all but native -> settlement
#   (native persistence and settlement are LUTS-defined at both dates).
# =============================================================================

source(file.path("scripts/landcover", "0a-setup.R"))

prop_dir     <- file.path("./processeddatasets", "predictors/30m")
f_prop_tame   <- "SET_TO_4D_TAME_PROPORTION.tif"     # canonical 4d layers, modelling extent
f_prop_annual <- "SET_TO_4D_ANNUAL_PROPORTION.tif"

# ---- Variants ----
# codes: ACI count layers summed as tame evidence. codes_2020 = NA: LUTS split.
variants <- tibble::tribble(
  ~variant,        ~codes_2010, ~min_2010, ~codes_2020, ~min_2020,
  "reference",     aci_tame_code,       tame_min_years,         NA,          NA,
  "min3",          "122",       3,         NA,          NA,
#  "aci110",        "110+122",   2,         NA,          NA,
  "aci_out",       "122",       2,         "122",       2,
#  "aci_out_drift", "122",       2,         "110+122",   2,
  "aci_out_min3",  "122",       3,         "122",       3,
#  "aci110_both",   "110+122",   2,         "110+122",   2
) |>
  mutate(chg_2010 = codes_2010 != codes_2010[1] | min_2010 != min_2010[1],
         chg_2020 = !is.na(codes_2020))

# ---- Load ----
sa_r           <- rast(file.path(sa_filepath, sa_raster_filename))
ecodistricts_r <- rast(tagged_path("./processeddatasets", ecodistricts_filepath, paste0(tools::file_path_sans_ext(ecodistricts_filename), ".tif"), "studyarea"))

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
# TBD Focal ag proportion (tame + annual) from 4d's canonical layers. The ag domain
# doesn't vary, so neither does this; variants recompute tame and get annual = ag - tame.
ag_prop_me <- rast(file.path(prop_dir, f_prop_tame)) + rast(file.path(prop_dir, f_prop_annual))

# ---- Relabelling check (decides aci_out_drift) ----
# Canonical tame that LUTS still calls perennial in 2020: how does ACI code it
# in 2019-2022? value: 0 = neither, 1 = 122 only, 10 = 110 only, 11 = both.
# A large 110-only share means 122 was relabelled 110: keep aci_out_drift.
# Otherwise drop it and rely on aci110_both.
# min_ref <- variants$min_2010[1]
# persist <- lulc_2010_sa == cls_tame & lulc_2020_sa == cls_tame
# relabel <- freq(ifel(persist, (count_2020_sa$aci110 >= min_ref) * 10 +
#                        (count_2020_sa$aci122 >= min_ref), NA))
# relabel$share <- relabel$count / sum(relabel$count)
# print(relabel)

# ---- Variant loop ----
tabs <- list()
for (i in seq_len(nrow(variants))) {
  v <- variants[i, ]
  
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

tab <- bind_rows(tabs)
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
write.csv(by_response, file.path(response_summary_filepath, "sensitivity", "conversion_by_response_by_variant.csv"), row.names = FALSE)

# Canonical must reproduce 7
canon <- read.csv(file.path(response_summary_filepath, "conversion_by_response.csv")) |>
  select(response, ref_pos = n_pos, ref_neg = n_neg) |>
  left_join(filter(by_response, variant == ref_variant), by = "response")
stopifnot("canonical does not reproduce 7" =
            all(canon$n_pos == canon$ref_pos & canon$n_neg == canon$ref_neg))