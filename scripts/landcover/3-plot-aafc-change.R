library(tidyverse)

#--- Setup ---
source(file.path("scripts", "0-setup.R"))

# === Set parameters ===
file_00_20 <- file.path(out_lulc_historic_dir, "aafc-lulc-change-2000-2020-km2.csv")
file_10_20 <- file.path(out_lulc_historic_dir, "aafc-lulc-change-2010-2020-km2.csv")

# === Load data ===
change_00_20 <- read_csv(file_00_20) %>%
  mutate(period = "2000–2020")
  
change_10_20 <- read_csv(file_10_20) %>%
  mutate(period = "2010–2020")

# === Combine datasets ===
change_all <- bind_rows(change_00_20, change_10_20) 


####################
# Plot LULC class summaries
###################

# === Summarize annual km² of change *from* each class ===
loss_summary <- change_all %>%
  group_by(period, From_Class) %>%
  summarise(
    total_area_annual = sum(area_annual, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(!From_Class %in% c("Other Land", "Water", "Settlement", "Wetland"))

# === Plot paired bar chart ===
ggplot(loss_summary, aes(x = reorder(From_Class, -total_area_annual), y = total_area_annual, fill = period)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  labs(
    title = "Annual LULC Loss by Class",
    x = "LULC Class (source of change)",
    y = "Annual Change (km²)",
    fill = "Period"
  ) +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



gain_summary <- change_all %>%
  filter(From_Class != To_Class) %>%  # exclude unchanged
  filter(!To_Class %in% c("Other Land", "Water", "Wetland")) %>%
  group_by(period, To_Class) %>%
  summarise(total_area_annual = sum(area_annual, na.rm = TRUE), .groups = "drop")

ggplot(gain_summary, aes(x = reorder(To_Class, -total_area_annual), y = total_area_annual, fill = period)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  labs(
    title = "Annual LULC Gain by Class",
    x = "LULC Class (gain destination)",
    y = "Annual Change (km²)",
    fill = "Period"
  ) +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))





gain_df <- change_all %>%
  filter(From_Class != To_Class) %>%
  filter(!To_Class %in% c("Other Land", "Water", "Wetland")) %>%
  group_by(period, Class = To_Class) %>%
  summarise(area_annual = sum(area_annual), .groups = "drop") %>%
  mutate(type = "Gain")

loss_df <- change_all %>%
  filter(From_Class != To_Class) %>%
  filter(!From_Class %in% c("Other Land", "Water", "Wetland")) %>%
  group_by(period, Class = From_Class) %>%
  summarise(area_annual = -sum(area_annual), .groups = "drop") %>%  # NOTE: make losses negative
  mutate(type = "Loss")

combined_df <- bind_rows(gain_df, loss_df)

net_cropland <- bind_rows(gain_df, loss_df) %>%
  filter(Class == "Cropland") %>%
  group_by(period, Class) %>%
  summarise(area_annual = sum(area_annual), .groups = "drop") %>%
  mutate(type = "Net")

combined_all <- bind_rows(gain_df, loss_df, net_cropland)

combined_all$period <- factor(combined_all$period, levels = c("2010–2020", "2000–2020"))

ggplot(combined_all, aes(x = reorder(Class, area_annual), y = area_annual, fill = period)) +
  geom_col(data = filter(combined_all, type != "Net"),
           position = position_dodge(width = 0.7), width = 0.6) +
  geom_col(data = filter(combined_all, type == "Net"),
           position = position_dodge(width = 0.7), width = 0.2, colour = "black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = "Annual LULC Change by Class (Net Shown for Cropland)",
    x = "LULC Class",
    y = "Annual Change (km²)",
    fill = "Period"
  ) +
  coord_flip() +
  theme_minimal(base_size = 13) +
  scale_fill_manual(
    values = c("2000–2020" = "#1f77b4", "2010–2020" = "#ff7f0e"),
    breaks = c("2000–2020", "2010–2020"),  # this sets the legend order
    name = "Period"
  )




















####################
# Plot transition types
###################

transition_summary <- change_all %>%
  filter(From_Class != To_Class) %>%  # remove no-change
  filter(!From_Class %in% c("Other Land", "Water", "Settlement", "Wetland")) %>%
  filter(!To_Class   %in% c("Other Land", "Water", "Wetland")) %>%
  mutate(Transition = paste(From_Class, "→", To_Class)) %>%
  group_by(period, Transition) %>%
  summarise(
    total_km2_annual = sum(km2_annual, na.rm = TRUE),
    .groups = "drop"
  )


ggplot(transition_summary, aes(x = reorder(Transition, -total_km2_annual), y = total_km2_annual, fill = period)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  labs(
    title = "Annual LULC Change by Transition Type",
    x = "Transition",
    y = "Annual Change (km²)",
    fill = "Period"
  ) +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))





