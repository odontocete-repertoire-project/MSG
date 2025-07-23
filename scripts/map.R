
library(ggplot2)
library(sf)
library(ggrepel)
library(ggspatial)
library(oce)
library(ocedata)
library(remotes)
library(dplyr)

remotes::install_github("ropensci/rnaturalearthhires")

library(rnaturalearth)


# High-res world map
world_hi <- ne_countries(scale = "large", returnclass = "sf")

# Shift latitude bounds slightly north to include more of Central America
region_bounds <- c(xmin = -85, xmax = -45, ymin = -15, ymax = 20)

# Locations of interest
locations <- data.frame(
  name = c("Brazil", "Panama", "Costa Rica"),
  lon = c(-48.5, -82.24, -82.65),
  lat = c(-0.72, 9.34, 9.63)
)

locations_sf <- st_as_sf(locations, coords = c("lon", "lat"), crs = 4326)

# Plot
region_map <- ggplot(data = world_hi) +
  geom_sf(fill = "gray95", color = "gray60") +
  geom_sf(data = locations_sf, color = "black", size = 2.5) +
  geom_text_repel(data = locations, aes(x = lon, y = lat, label = name), size = 4) + 
  labs(x = "Longitude", y = "Latitude") +
  coord_sf(xlim = c(-85, -45), ylim = c(-5, 15), expand = FALSE) +
  scale_x_continuous(breaks = seq(-85, -45, by = 5), name = "Longitude") +
  scale_y_continuous(breaks = seq(-2, 12, by = 5), name = "Latitude") +
  annotation_scale(location = "bl", width_hint = 0.4) +
  annotation_north_arrow(location = "tr", which_north = "true", 
                         style = north_arrow_fancy_orienteering) +
#  labs(title = "Regional Overview: Para Coast to Costa Rica") +
  theme_minimal(base_size=14) +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.ticks.length = unit(0.2, "cm"),
    plot.margin = margin(10, 10, 10, 10)
  )

print(region_map)
### GANDOCA MAP


library(osmdata)
library(readr)
library(sf)
library(dplyr)

# Define bounding box for Gandoca
bbox <- c(
  xmin = -82.806738,
  ymin = 9.535939,
  xmax = -82.514686,
  ymax = 9.740446
)

coast_query <- opq(bbox = bbox) |>
  add_osm_feature(key = 'natural', value = 'coastline')

coast_sf <- osmdata_sf(coast_query)$osm_lines

# Load the CSV (adjust the path if needed)
setwd("/Volumes/cure_bio/MaiaProjects/manzanilloanalysis/scripts")
sightings <- read_csv("sightings.csv")

# Remove rows with missing lat/lon
sightings_clean <- sightings %>%
  filter(!is.na(Latitude), !is.na(Longitude))

# Then convert to sf
sightings_sf <- st_as_sf(sightings_clean, coords = c("Longitude", "Latitude"), crs = 4326)


# Define bounding box for Gandoca
bbox <- c(
  xmin = -82.561289,
  ymin = 9.573194,
  xmax = -82.514686,
  ymax = 9.740446
)



library(RColorBrewer)

# Define your palettes
pinks <- c("#ffe8ee", "#ffd6e0", "#feb4c7", "#fe7799", "#fe0849", "#710a25")
blues <- RColorBrewer::brewer.pal(6, "Blues")
purples <- RColorBrewer::brewer.pal(6, "Purples")

# Custom color mapping
custom_palette <- c(
  "MSG" = purples[4],             # Medium purple
  "Guiana" = pinks[4],  # Medium pink
  "Bottlenose" = blues[5]  # Darker blue
)



# Define your desired order of species (bottom to top)
species_order <- c("Guiana", "Bottlenose", "MSG")

# Convert species to factor with that order
sightings_sf <- sightings_sf %>%
  mutate(Species = factor(Species, levels = species_order))

buffer <- 0.02

xlim_zoom <- c(-82.76126403 - buffer, -82.550188 + buffer)
ylim_zoom <- c(9.583419 - buffer, 9.667076999 + buffer)


gandoca <- ggplot() +
  geom_sf(data = coast_sf, color = "black", fill = "grey80") +
  geom_sf(data = sightings_sf %>% filter(Species == "Guiana"), 
          aes(color = Species), size = 3, alpha = 0.8) +
  geom_sf(data = sightings_sf %>% filter(Species == "Bottlenose"), 
          aes(color = Species), size = 3, alpha = 0.8) +
  geom_sf(data = sightings_sf %>% filter(Species == "MSG"), 
          aes(color = Species), size = 3, alpha = 0.8) +
  coord_sf(xlim = xlim_zoom,
           ylim = ylim_zoom,
           expand = FALSE) +
  scale_color_manual(values = custom_palette, na.value = "gray50") +
  annotation_scale(location = "bl", width_hint = 0.3)  +
  labs(
    x = "Longitude",
    y = "Latitude",
    color = "Species") +
  theme_minimal(base_size = 14)
  

print(gandoca)




library(cowplot)

combo <- plot_grid(
  region_map,
  gandoca,
  ncol = 1,
  labels = c("a", "b"),
  label_size = 14,       
  label_fontface = "bold", 
  align = "v",            
  rel_heights = c(1.2, 1)
)

ggsave('map.jpg', combo, dpi=300)

