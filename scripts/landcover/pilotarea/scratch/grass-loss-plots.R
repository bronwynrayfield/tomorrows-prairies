annual_loss <- summarize_grass_loss(
    grass_loss_class,
    grass_loss_year,
    out_units = "ha"
)


annual_loss_lumped <- annual_loss %>%
  # keep only the part before "/" (primary class)
  mutate(
    from_class_primary = str_remove(from_class, "/.*"),
    to_class_primary   = str_remove(to_class, "/.*")
  ) %>%
  # lump pasture subclasses (after stripping)
  mutate(
    from_group = case_when(
      from_class_primary %in% c("upast", "ipast") ~ "past",
      TRUE ~ from_class_primary
    ),
    to_group = case_when(
      to_class_primary %in% c("upast", "ipast") ~ "past",
      TRUE ~ to_class_primary
    )
  ) %>%
  group_by(year, from_group, to_group) %>%
  summarise(
    pixels  = sum(pixels),
    area_ha = sum(area_ha),
    .groups = "drop"
  )

# peek
dplyr::count(annual_loss_lumped, from_group, to_group, sort = TRUE)
head(annual_loss_lumped)


library(dplyr)
library(tidyr)
library(ggplot2)

# define the full year range (integer, no NAs)
all_years <- seq.int(
  min(annual_loss_lumped$year, na.rm = TRUE),
  max(annual_loss_lumped$year, na.rm = TRUE),
  by = 1
)

plot_df <- annual_loss_lumped %>%
  mutate(year = as.integer(year)) %>%
  complete(
    year = all_years,
    from_group,
    to_group,
    fill = list(pixels = 0L, area_ha = 0)
  )

# optional: nicer labels / ordering
to_levels <- c("crop", "builtup", "past", "forest", "shrub", "bare", "wetland", "water", "other")
plot_df <- plot_df %>%
  mutate(
    to_group = factor(to_group, levels = intersect(to_levels, unique(to_group))),
    from_group = factor(from_group)  # keep all
  )

# color palette for destination classes (tweak as you like)
to_cols <- c(
  crop    = "#A07040",  # brown
  builtup = "#800080",  # purple
  past    = "#2E8B57",
  forest  = "#228B22",
  shrub   = "#66A61E",
  bare    = "#808080",  # grey
  wetland = "#1E90FF",
  water   = "#4169E1",
  other   = "grey50"
)
to_cols <- to_cols[names(to_cols) %in% levels(plot_df$to_group)]


p_annual <- ggplot(plot_df, aes(x = year, y = area_ha, colour = to_group)) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = to_cols, name = "Converted to") +
  labs(
    x = "Year",
    y = "Area lost (ha)",
    title = "Annual grassland loss by destination class"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank()
  )

p_annual





# filter out unwanted destination classes
plot_df_filtered <- plot_df %>%
  filter(!to_group %in% c("water", "wetland", "bare"))

# re-subset your color vector to match what's left
to_cols_filt <- to_cols[names(to_cols) %in% levels(plot_df_filtered$to_group)]

# plot (same style)
p_annual_filtered <- ggplot(plot_df_filtered,
                            aes(x = year, y = area_ha, colour = to_group)) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = to_cols_filt, name = "Converted to") +
  labs(
    x = "Year",
    y = "Area lost (ha)",
    title = "Annual grassland loss by destination class"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank()
  )

p_annual_filtered


decadal <- annual_loss_lumped %>%
  # keep only the years of interest and (optionally) only grass origin
  filter(year >= 1990, year < 2020, from_group == "grass") %>%
  # left-closed, right-open bins: [1990,2000), [2000,2010), [2010,2020)
  mutate(period = cut(
    year,
    breaks = c(1990, 2000, 2010, 2020),
    right = FALSE,
    labels = c("1990–2000", "2000–2010", "2010–2020")
  )) %>%
  filter(!is.na(period)) %>%
  group_by(period, to_group) %>%
  summarise(
    pixels  = sum(pixels),
    area_ha = sum(area_ha),
    .groups = "drop"
  ) %>%
  arrange(period, desc(area_ha))

decadal

to_cols_dec <- to_cols[names(to_cols) %in% unique(decadal$to_group)]

p_decadal <- ggplot(decadal, aes(x = period, y = area_ha, fill = to_group)) +
  geom_bar(stat = "identity", position = "stack", colour = "grey20") +
  scale_fill_manual(values = to_cols_dec, name = "Converted to") +
  labs(
    x = "Decade",
    y = "Total grassland loss (ha)",
    title = "Total grassland loss by destination class and decade"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank()
  )

p_decadal

# filter out unwanted destination classes
decadal_filtered <- decadal %>%
  filter(!to_group %in% c("bare", "wetland", "water"))

# make sure colors only include the remaining to_classes
to_cols_dec <- to_cols[names(to_cols) %in% unique(decadal_filtered$to_group)]

# plot
p_decadal_filtered <- ggplot(decadal_filtered,
                             aes(x = period, y = area_ha, fill = to_group)) +
  geom_bar(stat = "identity", position = "stack", colour = "grey20") +
  scale_fill_manual(values = to_cols_dec, name = "Converted to") +
  labs(
    x = "Decade",
    y = "Total grassland loss (ha)",
    title = "Total grassland loss by destination class and decade"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank()
  )

p_decadal_filtered




pilot_area_ha <- sum(as.numeric(st_area(pilot))) / 10000
decadal_pct <- decadal_filtered %>%
  mutate(area_pct = 100 * area_ha / pilot_area_ha)



scaling_factor <- 100 / pilot_area_ha

to_cols_dec <- to_cols[names(to_cols) %in% unique(decadal_filtered$to_group)]

p_decadal_dual <- ggplot(decadal_filtered,
                         aes(x = period, y = area_ha, fill = to_group)) +
  geom_bar(stat = "identity", position = "stack", colour = "grey20") +
  scale_fill_manual(values = to_cols_dec, name = "Converted to") +
  scale_y_continuous(
    name = "Total grassland loss (ha)",
    sec.axis = sec_axis(~ . * scaling_factor,
                        name = "Grassland loss (% of pilot area)")
  ) +
  labs(
    x = "Decade",
    title = "Decadal grassland loss by destination class"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank()
  )

p_decadal_dual

