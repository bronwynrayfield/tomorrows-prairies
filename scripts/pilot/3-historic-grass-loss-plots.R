##############################################
# Assess historic grass, ipast, and upast loss (tabular)
# Plot annual and decadal results

#--- Setup ----
source(file.path("./scripts", "0-setup.R"))
source(file.path("./scripts/pilot", "functions.R"))

# Load files 
pilot <- st_read(file.path("./processeddatasets", pilot_filepath, pilot_s_filename))
grass_loss_class <- rast(file.path("./processeddatasets", pilot_filepath, grass_loss_class_filename))
grass_loss_year <- rast(file.path("./processeddatasets", pilot_filepath, grass_loss_year_filename))  
ipast_loss_class <- rast(file.path("./processeddatasets", pilot_filepath, ipast_loss_class_filename))
ipast_loss_year <- rast(file.path("./processeddatasets", pilot_filepath, ipast_loss_year_filename))
upast_loss_class <- rast(file.path("./processeddatasets", pilot_filepath, upast_loss_class_filename))
upast_loss_year <- rast(file.path("./processeddatasets", pilot_filepath, upast_loss_year_filename))


#--- Calculate historic loss rates ----
# grass loss
annual_grass_loss <- summarize_loss(
    grass_loss_class,
    grass_loss_year,
    out_units = "ha"
)

# ipast loss
annual_ipast_loss <- summarize_loss(
  ipast_loss_class,
  ipast_loss_year,
  out_units = "ha"
)

# upast loss
annual_upast_loss <- summarize_loss(
  upast_loss_class,
  upast_loss_year,
  out_units = "ha"
)



# Regroup and summarize annual data
annual_loss_lumped <- bind_rows(
  annual_grass_loss,
  annual_ipast_loss,
  annual_upast_loss
) %>%
  mutate(
    from_prim = str_remove(from_class, "/.*"),
    to_prim   = str_remove(to_class,   "/.*"),
    from_group = case_when(
      from_prim %in% c("upast","ipast") ~ "pasture",
      TRUE                               ~ from_prim
    ),
    to_group = case_when(
      to_prim %in% c("upast","ipast") ~ "pasture",
      TRUE                            ~ to_prim
    )
  ) %>%
  group_by(year, from_group, to_group) %>%
  summarise(pixels = sum(pixels), area_ha = sum(area_ha), .groups = "drop")



dplyr::count(annual_loss_lumped, from_group, to_group, sort = TRUE)
head(annual_loss_lumped)

#--- Annual plots ----
# ensure complete grid
to_levels <- c("crop","builtup","pasture","shrub","bare","wetland","water","grass")

plot_df <- annual_loss_lumped %>%
  mutate(year = as.integer(year)) %>%
  tidyr::complete(
    year = seq.int(min(year, na.rm=TRUE), max(year, na.rm=TRUE)),
    from_group, to_group,
    fill = list(pixels = 0L, area_ha = 0)
  ) %>%
  mutate(
    to_group   = factor(to_group, levels = to_levels),
    from_group = factor(from_group, levels = c("grass","pasture"))
  )  %>%
  # drop self-transitions (compare as character to avoid factor-level issues)
  filter(as.character(from_group) != as.character(to_group)) %>%
  filter(year >= 1990, year <= 2020)

# color palette for destination classes (tweak as you like)
to_cols <- c(
  crop    = "#A07040",  # brown
  builtup = "#800080",  # purple
  pasture = "#2E8B57",
  shrub   = "#A09010",
  bare    = "#808080",  # grey
  wetland = "#1E90FF",
  water   = "#4169E1",
  grass   = "#66A61E"
  )
to_cols <- to_cols[names(to_cols) %in% levels(plot_df$to_group)]

# Plot with all destination classes
p_annual <- ggplot(plot_df, aes(x = year, y = area_ha, colour = to_group)) +
  geom_line(linewidth = 1) +
  facet_wrap(~ from_group, scales = "fixed") +
  scale_colour_manual(values = to_cols, name = "Converted to") +
  labs(
    x = "Year",
    y = "Area lost (ha)"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top",
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6)
        )

p_annual

ggsave(
  "./results/landcover/pilot/annual_loss_all_1990_2020.png",
  plot = p_annual,
  width = 8, height = 5, units = "in",
  dpi = 300,        
  bg = "white"      
)

# Plot with filtered destination classes
plot_df_filtered <- plot_df %>%
  filter(!to_group %in% c("water", "wetland", "bare"))

# re-subset your color vector to match what's left
to_cols_filt <- to_cols[names(to_cols) %in% levels(plot_df_filtered$to_group)]


p_annual_filtered <- ggplot(plot_df_filtered, aes(x = year, y = area_ha, colour = to_group)) +
  geom_line(linewidth = 1) +
  facet_wrap(~ from_group, scales = "fixed") +
  scale_colour_manual(values = to_cols_filt, name = "Converted to") +
  labs(
    x = "Year",
    y = "Area lost (ha)"
    ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top",
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6)
        )

p_annual_filtered

ggsave(
  "./results/landcover/pilot/annual_loss_filtered_1990_2020.png",
  plot = p_annual_filtered,
  width = 8, height = 5, units = "in",
  dpi = 300,        
  bg = "white"      
)



#--- Decadal plots ----
decadal <- annual_loss_lumped %>%
  filter(year >= 1990, year < 2020,
         from_group %in% c("grass","pasture")) %>%
  mutate(period = cut(year,
                      breaks = c(1990, 2000, 2010, 2020),
                      right = FALSE,
                      labels = c("1990–2000","2000–2010","2010–2020"))) %>%
  filter(!is.na(period)) %>%
  group_by(period, from_group, to_group) %>%
  summarise(pixels = sum(pixels), area_ha = sum(area_ha), .groups = "drop") %>%
  arrange(period, from_group, desc(area_ha)) %>%
  # drop self-transitions (compare as character to avoid factor-level issues)
  filter(as.character(from_group) != as.character(to_group))


# Plot with all destination classes
p_decadal <- ggplot(decadal, aes(x = period, y = area_ha, fill = to_group)) +
  geom_bar(stat = "identity", colour = "grey20") +
  facet_wrap(~ from_group) +
  scale_fill_manual(values = to_cols, name = "Converted to") +
  labs(x = "Decade", y = "Total loss (ha)") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top",
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6)
  )

p_decadal

ggsave(
  "./results/landcover/pilot/decadal_loss_all_1990_2020.png",
  plot = p_decadal,
  width = 8, height = 5, units = "in",
  dpi = 300,        
  bg = "white"      
)

# Plot with filtered destination classes
decadal_filtered <- decadal %>%
  filter(!to_group %in% c("bare", "wetland", "water"))

# make sure colors only include the remaining to_classes
to_cols_dec <- to_cols[names(to_cols) %in% unique(decadal_filtered$to_group)]

p_decadal_filtered <- ggplot(decadal_filtered, aes(x = period, y = area_ha, fill = to_group)) +
  geom_bar(stat = "identity", colour = "grey20") +
  facet_wrap(~ from_group) +
  scale_fill_manual(values = to_cols_dec, name = "Converted to") +
  labs(x = "Decade", y = "Total loss (ha)") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top",
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6)
  )

p_decadal_filtered

ggsave(
  "./results/landcover/pilot/decadal_loss_filtered_1990_2020.png",
  plot = p_decadal_filtered,
  width = 8, height = 5, units = "in",
  dpi = 300,        
  bg = "white"      
)


# Plot with 2 y-axes
pilot_area_ha <- sum(as.numeric(st_area(pilot))) / 10000
scaling_factor <- 100 / pilot_area_ha


p_decadal_filtered_dual <- ggplot(decadal_filtered, aes(x = period, y = area_ha, fill = to_group)) +
  geom_bar(stat = "identity", colour = "grey20") +
  facet_wrap(~ from_group) +
  scale_fill_manual(values = to_cols_dec, name = "Converted to") +
  scale_y_continuous(
    name = "Total loss (ha)",
    sec.axis = sec_axis(~ . * scaling_factor,
                        name = "Total loss (% of pilot area)")
  ) +
  labs(x = "Decade") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top",
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6)
  )

p_decadal_filtered_dual

ggsave(
  "./results/landcover/pilot/decadal_loss_filtered_dual_1990_2020.png",
  plot = p_decadal_filtered_dual,
  width = 8, height = 5, units = "in",
  dpi = 300,        
  bg = "white"      
)


