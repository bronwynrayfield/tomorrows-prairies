
####################################################################################
##         PLOTTING GRASSLAND AREA ESTIMATES PER DATASET (PRAIRIES ECOZONE)       ##
####################################################################################

##----------------------- Alison Long - April 2026 -------------------------

# Using Sam Robinson's calculations for multi-year datasets
# See 

# ----------------------------
# 1) Load
# ----------------------------

# load packages
library(readxl)
library(tidyverse)
library(ggplot2)
library(ggrepel)

# set working directory
setwd("R:/CanadianGrasslands/WG_Classification/results")

# load data
df <- read_excel("R:/CanadianGrasslands/WG_Classification/processed/grass_area_pe.xlsx")

# examine data
head(df)
str(df)

# ----------------------------
# 2) Setup
# ----------------------------

# convert to long format
df_long <- df %>%
  select(-Other) %>%
  pivot_longer(
    cols = c(TotGrass, TameGrass, NatGrass),
    names_to = "GrassType",
    values_to = "Area"
  ) %>%
  filter(!is.na(Area))

# reorder GrassType levels
df_long <- df_long %>%
  mutate(
    GrassType = factor(
      GrassType,
      levels = c("TotGrass", "NatGrass", "TameGrass")
    )
  )

# create label data frame (last year per dataset and grass type)
label_df <- df_long %>%
  group_by(dataset, GrassType) %>%
  filter(year == min(year)) %>%
  ungroup()

# define labels for facets
facet_labels <- c(
  TotGrass  = "Total grassland",
  NatGrass  = "Natural grassland",
  TameGrass = "Tame grassland"
)

# ----------------------------
# 3) Plot
# ----------------------------

p <- ggplot(df_long, aes(x = year, y = Area/1000000, color = dataset, group = interaction(dataset, GrassType))) +
  geom_line() +
  geom_point(aes(shape = Provincial), size = 2.5) +
  #geom_text(data = label_df, aes(label = dataset), hjust = -0.1, size = 3, show.legend = FALSE) +
  geom_text_repel(data = label_df, aes(label = dataset), size = 3, show.legend = FALSE, direction = "y") +
  scale_shape_manual(
    values = c(`FALSE` = 16, `TRUE` = 17),
    name = "Provincial") +
  labs(x = "Year", y = "Grass area (millions of hectares)", color = "Dataset") +
  scale_y_continuous(labels = scales::label_number()) +
  theme_minimal() +
  theme(legend.position = "right",
        strip.text = element_text(size = 11),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)) +
  facet_wrap(~ GrassType, labeller = labeller(GrassType = facet_labels)) #, scales = "free_y")

# ----------------------------
# 3) Export to PNG
# ----------------------------
ggsave(
  filename = "grass_area_by_dataset.png",
  plot = p,
  width = 14,
  height = 6,
  units = "in",
  dpi = 300
)
