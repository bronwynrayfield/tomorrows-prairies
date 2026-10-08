# Purpose:
#   Derives the response rasters from the 3b LULC maps (study area):
#   - one multiclass response: the 2020 class of every pixel that was native or
#     tame in 2010 (origin = grass_origin from 3b). This is the training target
#     for the per-origin models;
#   - five binary option-A responses, one per transition, kept as benchmarks
#     and transition summaries:
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
#   native -> annual      18,085   0.23%     
#   native -> settlement  27,726   0.35%
#   native -> tame        36,903   0.46%
#   tame   -> annual   1,280,327  19.2%
#   tame   -> settlement  36,015   0.66%
#
# Known issues:
#   - Tame 2010 is ACI-defined and tame 2020 is LUTS-defined, so part of
#     tame -> annual is disagreement between products rather than conversion.
#     Variant aci_out sizes this.
# =============================================================================

source(file.path("scripts/landcover", "0a-setup.R"))

# ---- Load ----
sa              <- st_read(file.path(sa_filepath, sa_filename))
ecoregions_sa   <- st_read(tagged_path("./processeddatasets", ecoregions_filepath, ecoregions_filename, "studyarea"))
ecodistricts_sa <- st_read(tagged_path("./processeddatasets", ecodistricts_filepath, ecodistricts_filename, "studyarea"))
ecodistricts_r  <- rast(tagged_path("./processeddatasets", ecodistricts_filepath, paste0(tools::file_path_sans_ext(ecodistricts_filename), ".tif"), "studyarea"))

lulc_2010_sa <- rast(tagged_path("./processeddatasets", "landcover", "lulc_2010.tif", "study-area"))
lulc_2020_sa <- rast(tagged_path("./processeddatasets", "landcover", "lulc_2020.tif", "study-area"))

# ---- Binary responses ----
tcode <- lulc_2010_sa * 10 + lulc_2020_sa   # from * 10 + to

resp <- list()
for (r in responses$response) {
  s <- filter(scores, response == r)
  resp[[r]] <- classify(tcode, cbind(s$from * 10 + s$to, s$score), others = NA)
  names(resp[[r]]) <- r
  write_rast_tagged(resp[[r]], "./processeddatasets", "response", paste0("response_", r, ".tif"),
                    "study-area", datatype = "INT1U")
}

# ---- Multiclass response ----
# 2020 class for every pixel native or tame in 2010; NA elsewhere.
resp_mc <- mask(lulc_2020_sa, lulc_2010_sa %in% c(cls_native, cls_tame), maskvalues = FALSE)
names(resp_mc) <- "lulc_2020"
write_rast_tagged(resp_mc, "./processeddatasets", "response", "response_multiclass.tif",
                  "study-area", datatype = "INT1U")

# ---- Summaries ----
tab <- tabulate_lulc(lulc_2010_sa, lulc_2020_sa, ecodistricts_r) |>
  left_join(eco_lut(ecoregions_sa, ecodistricts_sa), by = "stratum_id")

summ <- tab |> inner_join(scores, by = c("from", "to"), relationship = "many-to-many")

by_response <- summ |> group_by(response) |> sum_scores()
print(by_response)
# A tibble: 5 × 7
# response             n_pos    n_neg n_competing   pos_ha    rate competing_per_pos
# <chr>                <dbl>    <dbl>       <dbl>    <dbl>   <dbl>             <dbl>
# 1 native_annual       200946 88005588      718096   18085. 0.00228            3.57  
# 2 native_settlement   308066 88005588      610976   27726. 0.00349            1.98  
# 3 native_tame         410030 88005588      509012   36903. 0.00464            1.24  
# 4 tame_annual       14225857 60059823      403211 1280327. 0.192              0.0283
# 5 tame_settlement     400162 60059823    14228906   36015. 0.00662           35.6  
write.csv(by_response, file.path(response_summary_filepath, "conversion_by_response.csv"), row.names = FALSE)

by_ecoregion <- summ |>
  filter(!is.na(ecoregion)) |>
  group_by(response, ecoregion, ecoregion_name) |>
  sum_scores() |>
  arrange(response, desc(rate))
print(by_ecoregion, n = Inf)
# A tibble: 35 × 9
# response          ecoregion ecoregion_name               n_pos    n_neg n_competing    pos_ha     rate competing_per_pos
# <chr>                 <int> <chr>                        <dbl>    <dbl>       <dbl>     <dbl>    <dbl>             <dbl>
#   1 native_annual           156 Aspen Parkland              108679  8711699      236525   9781.   0.0123             2.18   
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
# 22 tame_annual             163 Southwest Manitoba Uplands   78837   189785        1369   7095.   0.293              0.0174 
# 23 tame_annual             162 Lake Manitoba Plain         934863  3447673       67946  84138.   0.213              0.0727 
# 24 tame_annual             158 Fescue Grassland            368515  1435702        6005  33166.   0.204              0.0163 
# 25 tame_annual             156 Aspen Parkland             6751834 26762684      233420 607665.   0.201              0.0346 
# 26 tame_annual             157 Moist Mixed Grassland      2752391 11448853       54544 247715.   0.194              0.0198 
# 27 tame_annual             159 Mixed Grassland            3195913 15846758       39045 287632.   0.168              0.0122 
# 28 tame_annual             160 Cypress Upland              140652   916668         764  12659.   0.133              0.00543
# 29 tame_settlement         162 Lake Manitoba Plain          67873  3447673      934936   6109.   0.0193            13.8    
# 30 tame_settlement         156 Aspen Parkland              230746 26762684     6754508  20767.   0.00855           29.3    
# 31 tame_settlement         163 Southwest Manitoba Uplands    1369   189785       78837    123.   0.00716           57.6    
# 32 tame_settlement         157 Moist Mixed Grassland        54348 11448853     2752587   4891.   0.00472           50.6    
# 33 tame_settlement         158 Fescue Grassland              6002  1435702      368518    540.   0.00416           61.4    
# 34 tame_settlement         159 Mixed Grassland              38943 15846758     3196015   3505.   0.00245           82.1    
# 35 tame_settlement         160 Cypress Upland                 763   916668      140653     68.7  0.000832         184.
write.csv(by_ecoregion, file.path(response_summary_filepath, "conversion_by_ecoregion.csv"), row.names = FALSE)

# Multiclass outcome amounts by origin (one row of the transition matrix each)
multiclass_by_origin <- tab |>
  filter(from %in% c(cls_native, cls_tame)) |>
  group_by(from, to) |> summarise(n = sum(n), .groups = "drop") |>
  group_by(from) |> mutate(share = n / sum(n), ha = n * cell_ha) |> ungroup() |>
  mutate(origin  = lulc_levels[[2]][match(from, lulc_levels[[1]])],
         outcome = lulc_levels[[2]][match(to,   lulc_levels[[1]])])
print(multiclass_by_origin)
# A tibble: 8 × 7
# from    to        n     share       ha origin           outcome         
# <dbl> <dbl>    <dbl>     <dbl>    <dbl> <chr>            <chr>           
# 1     1     1 88005588 0.990     7920503. native grassland native grassland
# 2     1     2   410030 0.00461     36903. native grassland tame grassland  
# 3     1     3   200946 0.00226     18085. native grassland annual cropland 
# 4     1     4   308066 0.00346     27726. native grassland settlement      
# 5     2     2 60059823 0.804     5405384. tame grassland   tame grassland  
# 6     2     3 14225857 0.190     1280327. tame grassland   annual cropland 
# 7     2     4   400162 0.00536     36015. tame grassland   settlement      
# 8     2     5     3049 0.0000408     274. tame grassland   other   
write.csv(multiclass_by_origin, file.path(response_summary_filepath, "multiclass_by_origin.csv"), row.names = FALSE)

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

ggsave(file.path(response_summary_filepath, "response_maps.png"),
       plot = map, width = 14, height = 9, dpi = 300, bg = "white")