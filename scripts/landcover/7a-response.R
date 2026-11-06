# =============================================================================
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
# Study area, 2010-2020 (ha, from the 3b hindcast transition matrix):
#   native -> annual      18,085   0.23%
#   native -> settlement  27,726   0.35%
#   native -> tame        36,903   0.46%
#   tame   -> annual     925,194  13.8% of 2010 tame (16.6% of agreed pixels)
#   tame   -> settlement  36,015   0.54% of 2010 tame (0.77% of agreed pixels)
#
# Tame-origin training excludes pixels where ACI and LUTS disagree on whether
# forage remains in 2020 (12,517,439 px, 16.8% of 2010 tame): LUTS annual but
# ACI forage >= 1 yr (3,945,923), or LUTS perennial but no ACI forage (8,571,516).
# The hindcast map keeps them tame. Their share rises with conversion across
# ecodistricts (Spearman 0.72), so agreed negatives are weighted by ecodistrict
# (tame_negative_weights.csv) to restore each district's map rate.
#
# Known issues:
#   - Excluded pixels may differ from agreed ones within a district (e.g.
#     rotational fields); the weights correct the district mix, not that.
#   - The clean-switch rule removes proportionally more tame -> annual in the dry
#     south-west (Fescue -47%, Cypress -46%, Mixed -32% vs ~25% elsewhere),
#     consistent with the 7d census map, where the south shows least loss.
# =============================================================================

source(file.path("scripts/landcover", "0a-setup.R"))

# ---- Load ----
sa <- st_read(file.path(sa_filepath, sa_filename))

lulc_2010_sa <- rast(tagged_path("./processeddatasets", "landcover", "lulc_2010.tif", "study-area"))
lulc_2020_sa <- rast(tagged_path("./processeddatasets", "landcover", "lulc_2020.tif", "study-area"))
lulc_2020_luts_sa <- rast(tagged_path("./processeddatasets", "landcover", "lulc_2020_luts.tif", "study-area"))

ecoregions_sa   <- st_read(tagged_path("./processeddatasets", ecoregions_filepath, ecoregions_filename, "studyarea"))
ecodistricts_sa <- st_read(tagged_path("./processeddatasets", ecodistricts_filepath, ecodistricts_filename, "studyarea"))
ecodistricts_r  <- crop(rast(tagged_path("./processeddatasets", ecodistricts_filepath, paste0(tools::file_path_sans_ext(ecodistricts_filename), ".tif"), "studyarea")) , lulc_2010_sa)

# 2020 forage evidence (years coded 110 or 122, 2019-2022), from the 3b count rasters
count_2020 <- rast(sapply(aci_count_codes, \(cd)
                          tagged_path("./processeddatasets", file.path(aci_filepath, "30m"), count_name(cd, aci_years_2020))))
ev20 <- subst(crop(sum(count_2020), lulc_2010_sa), NA, 0)
# 2010 forage evidence
count_2010 <- rast(sapply(aci_count_codes, \(cd)
                          tagged_path("./processeddatasets", file.path(aci_filepath, "30m"), count_name(cd, aci_years_2010))))
ev10 <- subst(crop(sum(count_2010), lulc_2010_sa), NA, 0)

# Pixels where ACI and LUTS disagree on the 2020 outcome -> NA in training; the
# hindcast map keeps their 2010 class. Restricted to the 2020 ag domain.
#   tame origin:   LUTS annual vs no ACI forage (ev20 == 0); 12,517,439 px (16.8%):
#                  3,945,923 LUTS annual with forage, 8,571,516 LUTS perennial without
#   annual origin: LUTS perennial vs clean ACI entry (ev10 == 0, ev20 >= tame_min_years);
#                  23,048,388 px (8.9%): 18,372,285 LUTS perennial not clean,
#                  4,676,103 clean ACI entry LUTS calls annual. Positives 4,230,640.
ag20 <- lulc_2020_luts_sa %in% c(cls_tame, cls_annual)
uncertain <- ag20 & (
  (lulc_2010_sa == cls_tame   & (lulc_2020_luts_sa == cls_annual) != (ev20 == 0)) |
    (lulc_2010_sa == cls_annual & (lulc_2020_luts_sa == cls_tame) != (ev10 == 0 & ev20 >= tame_min_years)))
global(uncertain & lulc_2010_sa == cls_tame, "sum", na.rm = TRUE)   # must be 12,517,439
#           sum                        
# lulc_2020 12517439
freq(ifel(lulc_2010_sa == cls_annual & ag20,
          (lulc_2020_luts_sa == cls_tame) * 10 + (ev10 == 0 & ev20 >= tame_min_years), NA))
# layer value     count                   
# 1     1     0 230638336
# 2     1     1   4676103
# 3     1    10  18372285
# 4     1    11   4230640


# ---- Binary responses ----
tcode <- lulc_2010_sa * 10 + lulc_2020_sa   # from * 10 + to

resp <- list()
for (r in responses$response) {
  s <- filter(scores, response == r)
  resp[[r]] <- mask(classify(tcode, cbind(s$from * 10 + s$to, s$score), others = NA), uncertain, maskvalues = TRUE)
  names(resp[[r]]) <- r
  write_rast_tagged(resp[[r]], "./processeddatasets", "response", paste0("response_", r, ".tif"),
                    "study-area", datatype = "INT1U")
}

# ---- Multiclass response ----
# 2020 class for every pixel native or tame in 2010; NA elsewhere.
resp_mc <- mask(lulc_2020_sa, lulc_2010_sa %in% c(cls_native, cls_tame, cls_annual), maskvalues = FALSE)
resp_mc <- mask(resp_mc, uncertain, maskvalues = TRUE)
names(resp_mc) <- "lulc_2020"
write_rast_tagged(resp_mc, "./processeddatasets", "response", "response_multiclass.tif",
                  "study-area", datatype = "INT1U")

# ---- Summaries ----
tab <- tabulate_lulc(lulc_2010_sa, lulc_2020_sa, ecodistricts_r) |>
  left_join(eco_lut(ecoregions_sa, ecodistricts_sa), by = "stratum_id")

# Full 2010 x 2020 transition matrix (ha), all origins
lab <- \(x) lulc_levels[[2]][match(x, lulc_levels[[1]])]
trans_all <- tab |> group_by(from, to) |> summarise(ha = sum(n) * cell_ha, .groups = "drop") |>
  mutate(from = lab(from), to = lab(to))
xtabs(ha ~ from + to, trans_all) |> round() |> print()
#                  to
# from             annual cropland native grassland    other settlement tame grassland
# annual cropland         22829032                0      235      51423         380668
# native grassland           18084          7920150        0      27723          36898
# other                      52535                0  7709044       7152          43242
# settlement                     0                0        0     866402              0
# tame grassland            925012                0      274      36004        5759393


tab_resp <- tabulate_lulc(lulc_2010_sa, ifel(uncertain, 9, lulc_2020_sa), ecodistricts_r) |>
  left_join(eco_lut(ecoregions_sa, ecodistricts_sa), by = "stratum_id")
summ <- tab_resp |> inner_join(scores, by = c("from", "to"), relationship = "many-to-many")

by_response <- summ |> group_by(response) |> sum_scores()
print(by_response)
# # A tibble: 6 × 7
# response             n_pos     n_neg n_competing  pos_ha    rate competing_per_pos
# <chr>                <dbl>     <dbl>       <dbl>   <dbl>   <dbl>             <dbl>
#   1 annual_tame        4229644 230612690      573978 380668. 0.0180             0.136 
# 2 native_annual       200933  88001672      718011  18084. 0.00228            3.57  
# 3 native_settlement   308028  88001672      610916  27723. 0.00349            1.98  
# 4 native_tame         409983  88001672      508961  36898. 0.00464            1.24  
# 5 tame_annual       10277907  51478541      403093 925012. 0.166              0.0392
# 6 tame_settlement     400044  51478541    10280956  36004. 0.00771           25.7   
write.csv(by_response, file.path(response_summary_filepath, "conversion_by_response.csv"), row.names = FALSE)

by_ecoregion <- summ |>
  filter(!is.na(ecoregion)) |>
  group_by(response, ecoregion, ecoregion_name) |>
  sum_scores() |>
  arrange(response, desc(rate))
print(by_ecoregion, n = Inf)
# # A tibble: 42 × 9
# response          ecoregion ecoregion_name               n_pos    n_neg n_competing    pos_ha     rate competing_per_pos
# <chr>                 <int> <chr>                        <dbl>    <dbl>       <dbl>     <dbl>    <dbl>             <dbl>
# 1 annual_tame             160 Cypress Upland               50317   925880         332   4529.   0.0515             0.00660
# 2 annual_tame             158 Fescue Grassland            203724  6560419       43813  18335.   0.0301             0.215  
# 3 annual_tame             156 Aspen Parkland             2042340 91411521      352301 183811.   0.0219             0.172  
# 4 annual_tame             159 Mixed Grassland             909885 56877978       47397  81890.   0.0157             0.0521 
# 5 annual_tame             157 Moist Mixed Grassland       873320 58661631       89653  78599.   0.0147             0.103  
# 6 annual_tame             163 Southwest Manitoba Uplands   14747  1070346        1864   1327.   0.0136             0.126  
# 7 annual_tame             162 Lake Manitoba Plain         135311 15104915       38618  12178.   0.00888            0.285  
# 8 native_annual           156 Aspen Parkland              108679  8711703      236525   9781.   0.0123             2.18   
# 9 native_annual           163 Southwest Manitoba Uplands     107    23981         712      9.63 0.00444            6.65   
# 10 native_annual           162 Lake Manitoba Plain           1423   421341        3762    128.   0.00337            2.64   
# 11 native_annual           158 Fescue Grassland              9468  3735328       39734    852.   0.00253            4.20   
# 12 native_annual           157 Moist Mixed Grassland        35095 15672728      240471   3159.   0.00223            6.85   
# 13 native_annual           159 Mixed Grassland              43730 53338924      175909   3936.   0.000819           4.02   
# 14 native_annual           160 Cypress Upland                2431  6097667       20898    219.   0.000399           8.60   
# 15 native_settlement       156 Aspen Parkland              105029  8711703      240175   9453.   0.0119             2.29   
# 16 native_settlement       157 Moist Mixed Grassland       120484 15672728      155082  10844.   0.00763            1.29   
# 17 native_settlement       162 Lake Manitoba Plain           2779   421341        2406    250.   0.00655            0.866  
# 18 native_settlement       163 Southwest Manitoba Uplands     123    23981         696     11.1  0.00510            5.66   
# 19 native_settlement       158 Fescue Grassland             15549  3735328       33653   1399.   0.00415            2.16   
# 20 native_settlement       159 Mixed Grassland              59511 53338924      160128   5356.   0.00111            2.69   
# 21 native_settlement       160 Cypress Upland                4553  6097667       18776    410.   0.000746           4.12   
# 22 native_tame             163 Southwest Manitoba Uplands     589    23981         230     53.0  0.0240             0.390  
# 23 native_tame             156 Aspen Parkland              131496  8711703      213708  11835.   0.0149             1.63   
# 24 native_tame             157 Moist Mixed Grassland       119987 15672728      155579  10799.   0.00760            1.30   
# 25 native_tame             158 Fescue Grassland             24185  3735328       25017   2177.   0.00643            1.03   
# 26 native_tame             160 Cypress Upland               16345  6097667        6984   1471.   0.00267            0.427  
# 27 native_tame             162 Lake Manitoba Plain            983   421341        4202     88.5  0.00233            4.27   
# 28 native_tame             159 Mixed Grassland             116398 53338924      103241  10476.   0.00218            0.887  
# 29 tame_annual             163 Southwest Manitoba Uplands   58536   149061        1369   5268.   0.282              0.0234 
# 30 tame_annual             162 Lake Manitoba Plain         681366  2792706       67946  61323.   0.196              0.0997 
# 31 tame_annual             156 Aspen Parkland             5053050 22542712      233420 454774.   0.183              0.0462 
# 32 tame_annual             157 Moist Mixed Grassland      2053614  9879888       54544 184825.   0.172              0.0266 
# 33 tame_annual             158 Fescue Grassland            193739  1226730        6005  17437.   0.136              0.0310 
# 34 tame_annual             159 Mixed Grassland            2162185 14059308       39045 194597.   0.133              0.0181 
# 35 tame_annual             160 Cypress Upland               75417   828136         764   6788.   0.0835             0.0101 
# 36 tame_settlement         162 Lake Manitoba Plain          67873  2792706      681439   6109.   0.0237            10.0    
# 37 tame_settlement         156 Aspen Parkland              230746 22542712     5055724  20767.   0.0101            21.9    
# 38 tame_settlement         163 Southwest Manitoba Uplands    1369   149061       58536    123.   0.00910           42.8    
# 39 tame_settlement         157 Moist Mixed Grassland        54348  9879888     2053810   4891.   0.00547           37.8    
# 40 tame_settlement         158 Fescue Grassland              6002  1226730      193742    540.   0.00487           32.3    
# 41 tame_settlement         159 Mixed Grassland              38943 14059308     2162287   3505.   0.00276           55.5    
# 42 tame_settlement         160 Cypress Upland                 763   828136       75418     68.7  0.000920          98.8 
write.csv(by_ecoregion, file.path(response_summary_filepath, "conversion_by_ecoregion.csv"), row.names = FALSE)


# Ecoregion conversion table: % of 2010 origin area, not `rate` (which excludes competing
# transitions from its denominator, e.g. tame_settlement 0.66% vs 0.54% of tame area)
conv_table <- by_ecoregion |>
  separate(response, c("origin", "dest"), sep = "_") |>
  mutate(area_ha = (n_pos + n_neg + n_competing) * cell_ha, pct = 100 * pos_ha / area_ha) |>
  select(origin, ecoregion_name, area_ha, dest, pos_ha, pct) |>
  pivot_wider(names_from = dest, values_from = c(pos_ha, pct))



# Multiclass outcome amounts by origin (one row of the transition matrix each)
multiclass_by_origin <- tab_resp |>
  filter(from %in% c(cls_native, cls_tame, cls_annual)) |>
  group_by(from, to) |> summarise(n = sum(n), .groups = "drop") |>
  group_by(from) |> mutate(share = n / sum(n), ha = n * cell_ha) |> ungroup() |>
  mutate(origin  = lulc_levels[[2]][match(from, lulc_levels[[1]])],
         outcome = lulc_levels[[2]][match(to,   lulc_levels[[1]])])
print(multiclass_by_origin)
# # A tibble: 14 × 7
# from    to         n     share        ha origin           outcome         
# <dbl> <dbl>     <dbl>     <dbl>     <dbl> <chr>            <chr>           
# 1     1     1  88001672 0.990      7920150. native grassland native grassland
# 2     1     2    409983 0.00461      36898. native grassland tame grassland  
# 3     1     3    200933 0.00226      18084. native grassland annual cropland 
# 4     1     4    308028 0.00346      27723. native grassland settlement      
# 5     2     2  51478541 0.689      4633069. tame grassland   tame grassland  
# 6     2     3  10277907 0.138       925012. tame grassland   annual cropland 
# 7     2     4    400044 0.00536      36004. tame grassland   settlement      
# 8     2     5      3049 0.0000408      274. tame grassland   other           
# 9     2     9  12514716 0.168      1126324. tame grassland   NA              
# 10     3     2   4229644 0.0164      380668. annual cropland  tame grassland  
# 11     3     3 230612690 0.892     20755142. annual cropland  annual cropland 
# 12     3     4    571366 0.00221      51423. annual cropland  settlement      
# 13     3     5      2612 0.0000101      235. annual cropland  other           
# 14     3     9  23043225 0.0892     2073890. annual cropland  NA  
write.csv(multiclass_by_origin, file.path(response_summary_filepath, "multiclass_by_origin.csv"), row.names = FALSE)


# Negative weights by ecodistrict and origin: uncertain pixels keep their 2010 class on the
# map but are NA in training; persist / (persist - unc) restores each district's map rate.
#   tame:   rho 0.72, w 1.07-2.47 (median 1.25)
#   annual: rho 0.85, w 1.03-3.00 (median 1.12)
# largest annual weight: 846 (8,400 ha annual, 10.8% to tame, w 3.00); 750 (95,600 ha, w 1.65)
dest <- c(tame = cls_annual, annual = cls_tame)
zw <- bind_rows(lapply(names(dest), \(o) {
  orig <- lulc_2010_sa == cls[[o]]
  z <- zonal(c(orig, orig & lulc_2020_sa == dest[[o]], orig & lulc_2020_sa == cls[[o]], orig & uncertain),
             ecodistricts_r, "sum", na.rm = TRUE)
  names(z) <- c("ecodistrict", "orig", "conv", "persist", "unc")
  mutate(z, origin = o)
})) |>
  mutate(na_share = unc / orig, w_neg = persist / (persist - unc))
zw |> group_by(origin) |>
  summarise(rho = cor(na_share, conv / orig, method = "spearman"),
            w_min = min(w_neg), w_med = median(w_neg), w_max = max(w_neg))
write.csv(zw, file.path(response_summary_filepath, "negative_weights.csv"), row.names = FALSE)
# # A tibble: 2 × 5
# origin   rho w_min w_med w_max
# <chr>  <dbl> <dbl> <dbl> <dbl>
# 1 annual 0.847  1.03  1.12  3.00
# 2 tame   0.724  1.07  1.25  2.47

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