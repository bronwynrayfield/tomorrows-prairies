# Map of prior-corrected conversion risk (native + tame combined), coloured by
# percentile rank: bottom 80% as a gradient, then stepped top 20/10/5% classes.
#
# Percentiles rather than raw probability: native risk sits around 0.3% and
# tame around 15%, so a linear scale shows tame and almost nothing else. Ranks
# are comparable across strata only because both are prior-corrected.
# NB ranks here are of 300 m cell means, not 30 m cells.

source(file.path("scripts/landcover", "0a-setup.R"))
library(ggplot2); library(tidyterra)

pred_period <- 2010
out_pred    <- file.path("./processeddatasets", "predictions", "recent-conversion")
pred_file   <- file.path(out_pred, sprintf("pred_conversion_%d_300m_preview.tif", pred_period))
# after the 30 m run: sprintf("pred_conversion_%d_300m.tif", pred_period)

sa <- st_read(file.path(sa_filepath, sa_filename))

cuts <- c(0.5, 0.80, 0.90, 0.95)   # top 50 / 20 / 10 / 5 %

pred <- rast(pred_file)

# percentile rank of each cell within the surface (0-1)
v    <- values(pred, mat = FALSE)
rank <- setValues(pred, ecdf(v[!is.na(v)])(v))
cat("risk at top 50/20/10/5% thresholds:", signif(quantile(v, cuts, na.rm = TRUE), 3), "\n")

prov <- geodata::gadm("CAN", level = 1, path = tempdir()) |>
  subset(NAME_1 %in% c("Alberta", "Saskatchewan", "Manitoba"), NSE=TRUE) |>
  project(crs(pred)) |> crop(vect(sa))

cities <- vect(data.frame(name = c("Calgary", "Saskatoon", "Regina", "Winnipeg"),
                          lon  = c(-114.07, -106.67, -104.62, -97.14),
                          lat  = c(  51.04,   52.13,   50.45,  49.90)),
               geom = c("lon", "lat"), crs = "EPSG:4326") |>
  project(crs(pred))

# YlOrRd: gradient to the 80th percentile, hard steps above (paired stops)
p <- ggplot() +
  geom_spatvector(data = prov, fill = "grey95", colour = "grey60", linewidth = 0.3) +
  geom_spatraster(data = rank) +
  scale_fill_gradientn(
    colours = c("#ffffcc", "#fed976", "#fd8d3c", "#fd8d3c",
                "#e31a1c", "#e31a1c", "#800026", "#800026"),
    values  = c(0, cuts[1], cuts[1], cuts[2], cuts[2], cuts[3], cuts[3], 1),
    limits  = c(0, 1), breaks = cuts,
    labels  = c("Top 50%", "Top 20%", "Top 10%", "Top 5%"),
    na.value = NA, name = "Conversion risk") +
  geom_spatvector(data = cities, size = 1.5) +
  geom_spatvector_text(data = cities, aes(label = name),
                       size = 3.2, nudge_y = 25000, fontface = "bold") +
  guides(fill = guide_colourbar(barheight = unit(3.5, "cm"), barwidth = unit(0.4, "cm"))) +
#guides(fill = guide_colourbar(barheight = unit(6, "cm"))) +
  theme_void() +
  theme(legend.position = c(0.92, 0.5),legend.title = element_text(size = 9), legend.text = element_text(size = 5.5))


q  <- quantile(v, cuts, na.rm = TRUE)
df <- data.frame(risk = v[!is.na(v) & v > 0])   # log axis drops exact zeros
df$cls <- cut(df$risk, c(-Inf, q, Inf), labels = c("rest", "50", "20", "10", "5"))
lr <- log10(range(df$risk))
reg <- seq(lr[1], lr[2], length.out = 61)
br  <- sort(c(reg[sapply(reg, \(b) min(abs(b - log10(q)))) > 0.02], log10(q)))

h <- ggplot(df, aes(log10(risk), fill = cls)) +
  geom_histogram(breaks = br) +
  scale_x_continuous(breaks = -4:0,
                     labels = \(x) paste0(sapply(100 * 10^x, format, scientific = FALSE), "%")) +
#  scale_x_log10(labels = scales::label_percent(accuracy = 0.1)) +
  scale_fill_manual(values = c("#ffffcc", "#fed976", "#fd8d3c", "#e31a1c", "#800026"), guide = "none") +
  labs(x = "Conversion risk (10-yr)", y = NULL) +
  theme_minimal(base_size = 8) +
  theme(axis.text.y = element_blank(), panel.grid = element_blank(),
        plot.background = element_rect(fill = "white", colour = "grey80"))

p <- p + inset_element(h, left = 0.68, bottom = 0.66, right = 0.99, top = 0.99)


ggsave(file.path(out_pred, sprintf("map_conversion_%d.png", pred_period)),
       p, width = 10, height = 6, dpi = 300, bg = "white")
