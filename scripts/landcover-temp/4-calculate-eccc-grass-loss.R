#--- Setup ---
source(file.path("scripts", "0-setup.R"))

#--- Directories ---
# ECCC grassland loss directory
data_root <- file.path("existingdatasets/landcover/MostRecentChangeByClass/2022")

#--- Load prairie study area ---
sa <- st_read(ecozones_file)
ecoregions <- st_read(ecoregions_file)
template <- grass_loss  # or whichever raster you want as template
sa_proj <- sa |> st_make_valid() |> st_transform(crs(template))
sa_v <- vect(sa_proj)                 # preserves attributes, incl. ECOZONE
eco_id <- sa_proj$ECOZONE[1]         # tag results with this ID

grass_loss <- rast(file.path(data_root, "loss/c05_grass", "recentchange_2022_loss_c05_grass_b01_class.tif"))
upast_loss <- rast(file.path(data_root, "loss/c04_upast", "recentchange_2022_loss_c04_upast_b01_class.tif"))
ipast_loss <- rast(file.path(data_root, "loss/c10_ipast", "recentchange_2022_loss_c10_ipast_b01_class.tif"))

crop_gain <- rast(file.path(data_root, "gain/c03_crop", "recentchange_2022_gain_c03_crop_b01_class.tif"))
builtup_gain <- rast(file.path(data_root, "gain/c06_builtup", "recentchange_2022_gain_c06_builtup_b01_class.tif"))
upast_gain <- rast(file.path(data_root, "gain/c04_upast", "recentchange_2022_gain_c04_upast_b01_class.tif"))
ipast_gain <- rast(file.path(data_root, "gain/c10_ipast", "recentchange_2022_gain_c10_ipast_b01_class.tif"))

# --- Align helper (same as before) ---
align_to <- function(x, template) {
  if (!compareGeom(x, template, stopOnError = FALSE)) {
    resample(x, template, method = "near")
  } else x
}

loss_list <- list(
  c05_grass = align_to(grass_loss, template),
  c04_upast = align_to(upast_loss, template),
  c10_ipast = align_to(ipast_loss, template)
)

gain_list <- list(
  c03_crop    = align_to(crop_gain, template),
  c06_builtup = align_to(builtup_gain, template),
  c04_upast   = align_to(upast_gain, template),
  c10_ipast   = align_to(ipast_gain, template)
)

# --- Transition area for one from→to pair (same logic) ---
transition_area <- function(from_r, to_r, sa_v, change_val = 1) {
  r <- (from_r == change_val) & (to_r == change_val)
  r_sa <- mask(crop(r, sa_v), sa_v)
  cs <- cellSize(r_sa, unit = "m")
  area_m2 <- global(cs * (r_sa == 1), fun = "sum", na.rm = TRUE)[1, 1]
  if (is.na(area_m2)) 0 else area_m2
}

# --- Compute all transitions and include ECOZONE tag ---
lookup <- c(
  c05_grass = "Grass",
  c04_upast = "Tame pasture",
  c10_ipast = "Improved pasture",
  c03_crop = "Annual cropland",
  c06_builtup = "Built-up"
)

res <- imap_dfr(loss_list, function(from_r, from_name) {
  imap_dfr(gain_list, function(to_r, to_name) {
    tibble(
      ECOZONE = eco_id,
      from = lookup[[from_name]],
      to   = lookup[[to_name]],
      area_ha = transition_area(from_r, to_r, sa_v, change_val = 1) / 1e4
    )
  })
}) |>
  arrange(desc(area_ha))

print(res, n = Inf)









# --- 1) Precompute an SA mask on the template grid ---------------------------
template <- grass_loss  # or whichever you chose as the template
sa_proj  <- sa |> sf::st_make_valid() |> sf::st_transform(crs(template))
sa_v     <- vect(sa_proj)

# Binary mask: 1 inside SA, NA outside (cheap to reuse)
sa_mask  <- mask(template*0 + 1, sa_v)  # same extent/resolution as template

# Constant cell area (m²) in projected equal-area CRS
cell_area_m2 <- prod(res(template))

# --- 2) Faster transition-area function --------------------------------------
transition_area_fast <- function(from_r, to_r, sa_mask, change_val = 1) {
  # Align (nearest) if needed
  if (!compareGeom(from_r, sa_mask, stopOnError = FALSE)) {
    from_r <- resample(from_r, sa_mask, method = "near")
  }
  if (!compareGeom(to_r, sa_mask, stopOnError = FALSE)) {
    to_r <- resample(to_r, sa_mask, method = "near")
  }
  
  # Crop both to SA extent to reduce work
  from_c <- crop(from_r, sa_mask)
  to_c   <- crop(to_r, sa_mask)
  mask_c <- crop(sa_mask, from_c)
  
  # Logical intersection within SA (TRUE = 1, FALSE/NA = 0)
  # Use terra's arithmetic to keep it fast
  inter <- (from_c == change_val) & (to_c == change_val) & (mask_c == 1)
  
  # Count TRUE cells and multiply by constant cell area
  n_cells <- global(inter, "sum", na.rm = TRUE)[1, 1]
  if (is.na(n_cells)) n_cells <- 0
  area_m2 <- as.numeric(n_cells) * cell_area_m2
  area_m2
}

# --- 3) Use it in loop ---------------------------------------------------
# (Assumes loss_list and gain_list are already defined/aligned as before)
res <- purrr::imap_dfr(loss_list, function(from_r, from_name) {
  purrr::imap_dfr(gain_list, function(to_r, to_name) {
    message("Computing ", from_name, " → ", to_name)
    a_m2 <- transition_area_fast(from_r, to_r, sa_mask, change_val = 1)
    tibble::tibble(
      ECOZONE = sa_proj$ECOZONE[1],
      from_class = from_name,
      to_class   = to_name,
      area_ha    = a_m2 / 1e4
    )
  })
})






cell_area_m2 <- prod(res(template))        # constant in equal-area CRS
eco_id <- sa_proj$ECOZONE[1]

# --- Class ID mapping (from metadata) -----------------------------------
# we only need the classes present in lists
class_id <- c(
  c03_crop    = 3,
  c04_upast   = 4,
  c05_grass   = 5,
  c06_builtup = 6,
  c10_ipast   = 10
)

pretty_name <- c(
  c03_crop    = "Annual cropland",
  c04_upast   = "Tame pasture",
  c05_grass   = "Grassland",
  c06_builtup = "Built-up",
  c10_ipast   = "Improved pasture"
)

# --- rasters (already created above) ------------------------------------
loss_list <- list(
  c05_grass = grass_loss,
  c04_upast = upast_loss,
  c10_ipast = ipast_loss
)

gain_names <- c("c03_crop","c06_builtup","c04_upast","c10_ipast")  # "to" classes

# --- Helper: frequency of codes within SA for one "from" raster --------------
freq_in_sa <- function(r, sa_mask) {
  # align to mask if needed
  if (!compareGeom(r, sa_mask, stopOnError = FALSE)) {
    r <- resample(r, sa_mask, method = "near")
  }
  r_sa <- mask(r, sa_mask)   # NA outside SA
  
  f <- terra::freq(r_sa, digits = 0, value = TRUE)
  if (is.null(f) || nrow(f) == 0) {
    return(tibble::tibble(value = integer(), count = integer()))
  }
  tibble::tibble(value = as.integer(f[, "value"]),
                 count = as.numeric(f[, "count"]))
}


# --- Main: compute areas by looking up codes from freq table ------------------
res <- imap_dfr(loss_list, function(from_r, from_name) {
  message("Indexing codes for from = ", from_name)
  f_tbl <- freq_in_sa(from_r, sa_mask)
  
  from_id <- class_id[[from_name]]
  to_ids  <- class_id[gain_names]
  codes   <- from_id*100 + to_ids
  
  # look up counts for each target code
  counts <- tibble(code = as.integer(codes),
                   to_name = names(to_ids)) |>
    left_join(f_tbl, by = c("code" = "value")) |>
    mutate(count = replace_na(count, 0),
           area_ha = count * cell_area_m2 / 1e4,
           from_name = from_name)
  
  counts |>
    transmute(
      ECOZONE = eco_id,
      from = pretty_name[from_name],
      to   = pretty_name[to_name],
      area_ha
    )
})

# Order nicely
res <- res |> arrange(desc(area_ha))
print(res, n = Inf)

# Optional rollups
totals_from <- res |> group_by(ECOZONE, from) |> summarise(area_ha = sum(area_ha), .groups = "drop")
totals_to   <- res |> group_by(ECOZONE, to)   |> summarise(area_ha = sum(area_ha), .groups = "drop")
totals_from; totals_to


f_grass <- freq_in_sa(grass_loss, sa_mask)
head(f_grass[order(f_grass$value), ], 15)  # shows codes like 501, 503, 506, etc.
# Grass -> Crop should be code 503
f_grass |> filter(value == 503)





transition_area_code <- function(from_r, from_id, to_id, sa_mask) {
  code <- from_id * 100 + to_id
  
  # align
  if (!compareGeom(from_r, sa_mask, stopOnError = FALSE)) {
    from_r <- resample(from_r, sa_mask, method = "near")
  }
  
  # crop + mask
  r_sa <- mask(crop(from_r, sa_mask), sa_mask)
  
  # count cells equal to code
  hits <- (r_sa == code)
  n_cells <- global(hits, "sum", na.rm = TRUE)[1, 1]
  if (is.na(n_cells)) n_cells <- 0
  
  # cell area (constant in equal-area CRS)
  cell_area_m2 <- prod(res(sa_mask))
  n_cells * cell_area_m2 / 1e4   # hectares
}

from_id <- 5   # grass
to_id   <- 3   # crop

transition_area_code(grass_loss, from_id, to_id, sa_mask)






loss_list <- list(
  c05_grass = grass_loss,
  c04_upast = upast_loss,
  c10_ipast = ipast_loss
)
gain_ids <- c(c03_crop = 3, c06_builtup = 6, c04_upast = 4, c10_ipast = 10)

res <- purrr::imap_dfr(loss_list, function(from_r, from_name) {
  from_id <- class_id[[from_name]]
  purrr::imap_dfr(gain_ids, function(to_id, to_name) {
    area_ha <- transition_area_code(from_r, from_id, to_id, sa_mask)
    tibble::tibble(
      ECOZONE = eco_id,
      from    = pretty_name[[from_name]],
      to      = pretty_name[[to_name]],
      area_ha = area_ha
    )
  })
})


library(dplyr)

# res_annual has: ECOZONE, from, to, area_ha, area_km2_total, area_km2_peryr

years <- 2022 - 1984  # = 38

res_annual <- res %>%
  dplyr::mutate(
    area_km2_total = area_ha / 100,              # total km² over 1984–2022
    area_km2_peryr = area_km2_total / years      # annualized km²/yr
  ) %>%
  arrange(desc(area_km2_peryr))

print(res_annual, n = Inf)

library(dplyr)

# 1) Annual roll-up by FROM class (total loss per year from each source class)
annual_by_from <- res_annual %>%
  group_by(ECOZONE, from) %>%
  summarise(
    km2_peryr = sum(area_km2_peryr),
    total_km2 = sum(area_km2_total),
    total_ha  = sum(area_ha),
    .groups = "drop"
  ) %>%
  arrange(desc(km2_peryr))

# 2) Annual roll-up by TO class (total gain per year into each destination class)
annual_by_to <- res_annual %>%
  group_by(ECOZONE, to) %>%
  summarise(
    km2_peryr = sum(area_km2_peryr),
    total_km2 = sum(area_km2_total),
    total_ha  = sum(area_ha),
    .groups = "drop"
  ) %>%
  arrange(desc(km2_peryr))

# 3) Net annual change by class (gains – losses), km²/yr
net_by_class <- full_join(
  annual_by_to %>% rename(class = to, gains_km2_peryr = km2_peryr),
  annual_by_from %>% rename(class = from, losses_km2_peryr = km2_peryr),
  by = c("ECOZONE", "class")
) %>%
  mutate(
    gains_km2_peryr  = coalesce(gains_km2_peryr, 0),
    losses_km2_peryr = coalesce(losses_km2_peryr, 0),
    net_km2_peryr    = gains_km2_peryr - losses_km2_peryr
  ) %>%
  select(ECOZONE, class, gains_km2_peryr, losses_km2_peryr, net_km2_peryr) %>%
  arrange(desc(net_km2_peryr))

# 4) Nice one-liners if you want whole-SA totals (ignore ECOZONE)
annual_by_from_all <- res_annual %>%
  group_by(from) %>%
  summarise(km2_peryr = sum(area_km2_peryr), .groups = "drop") %>%
  arrange(desc(km2_peryr))

annual_by_to_all <- res_annual %>%
  group_by(to) %>%
  summarise(km2_peryr = sum(area_km2_peryr), .groups = "drop") %>%
  arrange(desc(km2_peryr))

net_by_class_all <- full_join(
  annual_by_to_all %>% rename(class = to, gains_km2_peryr = km2_peryr),
  annual_by_from_all %>% rename(class = from, losses_km2_peryr = km2_peryr),
  by = "class"
) %>%
  mutate(
    gains_km2_peryr  = coalesce(gains_km2_peryr, 0),
    losses_km2_peryr = coalesce(losses_km2_peryr, 0),
    net_km2_peryr    = gains_km2_peryr - losses_km2_peryr
  ) %>%
  arrange(desc(net_km2_peryr))

# 5) (Optional) export to CSV
# write.csv(annual_by_from,     "annual_by_from.csv", row.names = FALSE)
# write.csv(annual_by_to,       "annual_by_to.csv",   row.names = FALSE)
# write.csv(net_by_class,       "net_by_class.csv",   row.names = FALSE)
# write.csv(annual_by_from_all, "annual_by_from_all.csv", row.names = FALSE)
# write.csv(annual_by_to_all,   "annual_by_to_all.csv",   row.names = FALSE)
# write.csv(net_by_class_all,   "net_by_class_all.csv",   row.names = FALSE)









