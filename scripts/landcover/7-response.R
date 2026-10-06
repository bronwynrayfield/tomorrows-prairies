# Purpose:
#   Derives five binary response rasters from the 3b LULC maps (study area),
#   one per transition CLUES needs a suitability layer for:
#     native -> annual, tame -> annual, native -> settlement,
#     tame -> settlement, native -> tame
#   Summarises each by response and by ecoregion, and maps them.
#
# Scoring (option A): 1 = focal transition, 0 = stayed in the origin class,
#   NA = any other transition from the same origin, and all other origins.
#   Each model therefore contrasts one destination with persistence, so every
#   layer shares the same reference and competition between destinations is
#   left to the CLUES allocation rather than baked into the layers.
#   The choice is material: competing transitions outnumber positives 3.6x
#   (native->annual), 2.0x (native->settlement), 1.2x (native->tame) and
#   34.9x (tame->settlement); only tame->annual is insensitive (0.03x).
#
# Study area, 2010-2020 (ha, from the 3b transition matrix):
#   native -> annual      18,085   0.23%     tame -> annual    992,068  17.6%
#   native -> settlement  27,726   0.35%     tame -> settlement 28,452   0.61%
#   native -> tame        36,903   0.46%
#   Native loss totals 82,714 ha (= LUTS native decline); 45% of it is to tame.
#
# Known issues:
#   - Tame 2010 is ACI-defined and tame 2020 is LUTS-defined, so part of
#     tame -> annual is disagreement between products rather than conversion.
#     Variant aci_out sizes this.
#   - The tame -> annual vs native -> annual base-rate gap is ~77x. Stitching
#     origin-specific layers into one CLUES suitability layer requires
#     calibrated probabilities (e.g. isotonic on out-of-fold predictions).
#   - Native in 2010 is the residue that was not worth ploughing, so native
#     layers rank conversion within a truncated suitability range.
#   - ACI forage classes separate poorly in Fescue Grassland, so tame
#     measurement error varies across LOEO folds.
# =============================================================================

source(file.path("scripts/landcover", "0a-setup.R"))

summary_dir <- "./output/landcover/conversion-response"

# ---- Responses ----
responses <- tibble::tribble(
  ~response,           ~from,      ~to,
  "native_annual",     cls_native, cls_annual,
  "tame_annual",       cls_tame,   cls_annual,
  "native_settlement", cls_native, cls_settlement,
  "tame_settlement",   cls_tame,   cls_settlement,
  "native_tame",       cls_native, cls_tame
)

scores <- bind_rows(Map(\(r, f, t) mutate(score_response(f, t), response = r),
                        responses$response, responses$from, responses$to))

# ---- Load ----
sa              <- st_read(file.path(sa_filepath, sa_filename))
ecoregions_sa   <- st_read(tagged_path("./processeddatasets", ecoregions_filepath, ecoregions_filename, "studyarea"))
ecodistricts_sa <- st_read(tagged_path("./processeddatasets", ecodistricts_filepath, ecodistricts_filename, "studyarea"))
ecodistricts_r  <- rast(tagged_path("./processeddatasets", ecodistricts_filepath, paste0(tools::file_path_sans_ext(ecodistricts_filename), ".tif"), "studyarea"))

lulc_2010_sa <- rast(tagged_path("./processeddatasets", "landcover", "lulc_2010.tif", "study-area"))
lulc_2020_sa <- rast(tagged_path("./processeddatasets", "landcover", "lulc_2020.tif", "study-area"))

# ---- Response rasters ----
tcode <- lulc_2010_sa * 10 + lulc_2020_sa   # from * 10 + to

resp <- list()
for (r in responses$response) {
  s <- filter(scores, response == r)
  resp[[r]] <- classify(tcode, cbind(s$from * 10 + s$to, s$score), others = NA)
  names(resp[[r]]) <- r
  write_rast_tagged(resp[[r]], "./processeddatasets", "response", paste0("response_", r, ".tif"),
                    "study-area", datatype = "INT1U")
}

# ---- Summaries ----
cell_ha <- (res_m^2) / 1e4

summ <- tabulate_lulc(lulc_2010_sa, lulc_2020_sa, ecodistricts_r) |>
  left_join(eco_lut(ecoregions_sa, ecodistricts_sa), by = "stratum_id") |>
  inner_join(scores, by = c("from", "to"), relationship = "many-to-many")

sum_scores <- function(d) summarise(d,
                                    n_pos       = sum(n[score %in% 1]),
                                    n_neg       = sum(n[score %in% 0]),
                                    n_competing = sum(n[is.na(score)]),
                                    pos_ha      = n_pos * cell_ha,
                                    rate        = n_pos / (n_pos + n_neg),
                                    competing_per_pos = n_competing / n_pos,
                                    .groups = "drop")

by_response <- summ |> group_by(response) |> sum_scores()
print(by_response)
# A tibble: 5 × 7
# response             n_pos    n_neg n_competing  pos_ha    rate competing_per_pos
# <chr>                <dbl>    <dbl>       <dbl>   <dbl>   <dbl>             <dbl>
# 1 native_annual       200946 88005588      718096  18085. 0.00228            3.57  
# 2 native_settlement   308066 88005588      610976  27726. 0.00349            1.98  
# 3 native_tame         410030 88005588      509012  36903. 0.00464            1.24  
# 4 tame_annual       11022980 51770972      318806 992068. 0.176              0.0289
# 5 tame_settlement     316133 51770972    11025653  28452. 0.00607            34.9 
write.csv(by_response, file.path(summary_dir, "conversion_by_response.csv"), row.names = FALSE)

by_ecoregion <- summ |>
  filter(!is.na(ecoregion)) |>
  group_by(response, ecoregion, ecoregion_name) |>
  sum_scores() |>
  arrange(response, desc(rate))
print(by_ecoregion, n = Inf)
# A tibble: 35 × 9
# response          ecoregion ecoregion_name               n_pos    n_neg n_competing    pos_ha     rate competing_per_pos
# <chr>                 <int> <chr>                        <dbl>    <dbl>       <dbl>     <dbl>    <dbl>             <dbl>
# 1 native_annual           156 Aspen Parkland              108679  8711699      236525   9781.   0.0123             2.18   
# 2 native_annual           163 Southwest Manitoba Uplands     107    23981         712      9.63 0.00444            6.65   
# 3 native_annual           162 Lake Manitoba Plain           1423   421341        3762    128.   0.00337            2.64   
# 4 native_annual           158 Fescue Grassland              9468  3735329       39734    852.   0.00253            4.20   
# 5 native_annual           157 Moist Mixed Grassland        35095 15672705      240471   3159.   0.00223            6.85   
# 6 native_annual           159 Mixed Grassland              43730 53338941      175909   3936.   0.000819           4.02   
# 7 native_annual           160 Cypress Upland                2431  6097661       20898    219.   0.000399           8.60   
# 8 native_settlement       156 Aspen Parkland              105029  8711699      240175   9453.   0.0119             2.29   
# 9 native_settlement       157 Moist Mixed Grassland       120484 15672705      155082  10844.   0.00763            1.29   
# 10 native_settlement       162 Lake Manitoba Plain           2779   421341        2406    250.   0.00655            0.866  
# 11 native_settlement       163 Southwest Manitoba Uplands     123    23981         696     11.1  0.00510            5.66   
# 12 native_settlement       158 Fescue Grassland             15549  3735329       33653   1399.   0.00415            2.16   
# 13 native_settlement       159 Mixed Grassland              59511 53338941      160128   5356.   0.00111            2.69   
# 14 native_settlement       160 Cypress Upland                4553  6097661       18776    410.   0.000746           4.12   
# 15 native_tame             163 Southwest Manitoba Uplands     589    23981         230     53.0  0.0240             0.390  
# 16 native_tame             156 Aspen Parkland              131496  8711699      213708  11835.   0.0149             1.63   
# 17 native_tame             157 Moist Mixed Grassland       119987 15672705      155579  10799.   0.00760            1.30   
# 18 native_tame             158 Fescue Grassland             24185  3735329       25017   2177.   0.00643            1.03   
# 19 native_tame             160 Cypress Upland               16345  6097661        6984   1471.   0.00267            0.427  
# 20 native_tame             162 Lake Manitoba Plain            983   421341        4202     88.5  0.00233            4.27   
# 21 native_tame             159 Mixed Grassland             116398 53338941      103241  10476.   0.00218            0.887  
# 22 tame_annual             163 Southwest Manitoba Uplands   44081   100340         643   3967.   0.305              0.0146 
# 23 tame_annual             162 Lake Manitoba Plain         467576  1160554       27271  42082.   0.287              0.0583 
# 24 tame_annual             158 Fescue Grassland            280314  1091603        3798  25228.   0.204              0.0135 
# 25 tame_annual             156 Aspen Parkland             5662723 23333565      200274 509645.   0.195              0.0354 
# 26 tame_annual             157 Moist Mixed Grassland      2006629 10576066       49419 180597.   0.159              0.0246 
# 27 tame_annual             159 Mixed Grassland            2468872 14687413       36711 222198.   0.144              0.0149 
# 28 tame_annual             160 Cypress Upland               90312   811661         596   8128.   0.100              0.00660
# 29 tame_settlement         162 Lake Manitoba Plain          27227  1160554      467620   2450.   0.0229            17.2    
# 30 tame_settlement         156 Aspen Parkland              197925 23333565     5665072  17813.   0.00841           28.6    
# 31 tame_settlement         163 Southwest Manitoba Uplands     643   100340       44081     57.9  0.00637           68.6    
# 32 tame_settlement         157 Moist Mixed Grassland        49229 10576066     2006819   4431.   0.00463           40.8    
# 33 tame_settlement         158 Fescue Grassland              3795  1091603      280317    342.   0.00346           73.9    
# 34 tame_settlement         159 Mixed Grassland              36624 14687413     2468959   3296.   0.00249           67.4    
# 35 tame_settlement         160 Cypress Upland                 596   811661       90312     53.6  0.000734         152.  
write.csv(by_ecoregion, file.path(summary_dir, "conversion_by_ecoregion.csv"), row.names = FALSE)

# ---- Map ----
# Aggregate to 600 m (fact = 20) for plotting
response_df <- bind_rows(lapply(responses$response, \(r)
                                as.data.frame(aggregate(resp[[r]], fact = 20, fun = "max", na.rm = TRUE), xy = TRUE) |>
                                  setNames(c("x", "y", "value")) |>
                                  filter(value == 1) |>
                                  mutate(response = r)))

map <- ggplot() +
  geom_tile(data = response_df, aes(x = x, y = y), fill = "darkred", alpha = 0.8) +
  geom_sf(data = ecoregions_sa, fill = NA, colour = "grey40", linewidth = 0.3) +
  geom_sf(data = sa, fill = NA, colour = "black", linewidth = 0.7) +
  facet_wrap(~ response) +
  theme_void()

ggsave(file.path(summary_dir, "response_maps.png"),
       plot = map, width = 14, height = 9, dpi = 300, bg = "white")