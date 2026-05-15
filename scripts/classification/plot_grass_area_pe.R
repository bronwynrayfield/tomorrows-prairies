
#######################################################################################
##              PLOT TOTAL GRASSLAND AREA BY DATASET            ##
#######################################################################################

##----------------------- Alison Long - April 2026 -------------------------

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)  # for label placement
library(ggthemes) # for colorblind palette

data_dir <- "R:/CanadianGrasslands/WG_Classification/results"

df <- read_excel(file.path(data_dir,"grass_area_pe_noprov.xlsx")) # without provincial datasets
#df <- read_excel(file.path(data_dir,"grass_area_pe.xlsx"))       # with provincial datasets

str(df)

# ------------------------ data wrangling --------------------------

# filter out Plowprint
df <- df %>%
  filter(dataset != "WWF")

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

# create a label data frame
label_df <- df_long %>%
  group_by(dataset, GrassType) %>%
  filter(year == min(year)) %>%
  ungroup()

# rename panels
facet_labels <- c(
  TotGrass  = "Total grassland",
  NatGrass  = "Natural grassland",
  TameGrass = "Tame grassland"
)

# rename datasets for legend
dataset_labels <- c(
  "ACI" = "AAFC Annual Crop Inventory (ACI)",
  "CEC" = "CEC NALCMS Landcover Map of North America (CEC)",
  "NGI" = "CFGA National Grassland Inventory (NGI)",
  "CGR" = "Central Grasslands Roadmap (CGR)",
  "CNVC" = "NatureServe Canadian National Vegetation Classification (CNVC)",
  "DW" = "WRI/Google Dynamic World (DW)",
  "ECCC" = "ECCC Canadian Prairies Grassland Change Product Summary (ECCC)",
  "ESA" = "ESA WorldCover (ESA)",
  "ESRI" = "ESRI LULC Time Series (ESRI)",
  "GLC" = "GLC_FCS30D (GLC)",
  "GPW" = "WRI Global Pasture Watch (GPW)",
  "LCC" = "NRCan Land Cover of Canada (LCC)",
  "PHJV" = "PHJV Potentially Undisturbed Lands (PHJV)",
  "SDLU" = "AAFC Semi-Decadal Land Use (SDLU)",
  "LCR" = "StatCan Land Cover Register (LCR)"
  #"WWF" = "WWF Plowprint"
)

# load colorblind-friendly palette
cb_pal <- ggthemes::colorblind_pal()(8)

# define shapes
default_shapes <- c(16, 16, 16, 16, 16, 16, 16, 16)
extra_shapes <- c(15, 17, 1, 2, 23, 0, 18)
all_shapes <- c(default_shapes, extra_shapes)

#-------------------------- plot -------------------------------
p <- ggplot(df_long,
            aes(x = year, y = Area/1000000, color = dataset, group = interaction(dataset, GrassType))) +
  geom_line() +
  #geom_point(aes(shape = Provincial), size = 2.5) +
  geom_point(aes(shape = dataset), size = 2.5) +
  geom_text_repel(data = label_df, aes(label = dataset), size = 3, show.legend = FALSE) + #direction = "y") +
  #geom_text_repel(data = label_df, aes(label = dataset), #nudge_y = 0.05,      
    #direction = "y", 
    #box.padding = 0.4, point.padding = 0.3, show.legend = FALSE) +
  scale_color_manual(name = "Dataset",
    labels = function(x) stringr::str_wrap(dataset_labels[x], width = 25),
    values = scales::hue_pal()(length(dataset_labels))) +
  #scale_color_manual(name   = "Dataset",
    #labels = function(x) stringr::str_wrap(dataset_labels[x], width = 25),
    #values = rep_len(cb_pal, length(dataset_labels))) +
  #scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), name = "Provincial") +
  scale_shape_manual(name   = "Dataset",
    values = all_shapes[seq_along(dataset_labels)]) +
  scale_y_continuous(labels = scales::label_number()) +
  labs(x = "Year", y = "Grassland area (millions of hectares)", color = "Dataset") +
  theme_minimal() +
  theme(strip.text = element_text(size = 14, face = "bold"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        panel.spacing = unit(1, "lines"),
        legend.text  = element_text(size = 6), # make legend text smaller
        legend.key.height = unit(8, "mm"), # add vertical space between legend items
        legend.position = "right") +
  facet_wrap(~ GrassType, labeller = labeller(GrassType = facet_labels))

p

#----------------------- export ----------------------------
ggsave(
  filename = file.path(data_dir, "grass_area_by_dataset_noprov_v10.png"),
  plot = p,
  width = 14,
  height = 6,
  units = "in",
  dpi = 300
)
