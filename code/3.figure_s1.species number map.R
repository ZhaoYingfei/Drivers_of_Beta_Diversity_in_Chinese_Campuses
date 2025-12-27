# Cleaning memory
cat("\014")
rm(list = ls())
gc()

library(sf)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(purrr)
library(cowplot)

setwd("C:/Users/qianh/Desktop/R/cam.pla.div")
getwd()

# Read data files
campus_checklist <- read.csv("./data/Campus_checklist_status.csv") # Campus plant species checklist
driving_factors <- read.csv("./data/Driving_factors.csv") # Campus driving factors data
china_map_province <- st_read("./data/China_map/bou2_4p.shp") # China provincial boundaries shapefile
china_map_southsea <- st_read("./data/China_map/south_sea.shp") # South China Sea boundaries shapefile

# Set coordinate reference system for province map
china_map_province <- st_set_crs(china_map_province, 4326)

# Data preprocessing: Convert variables to factors
campus_checklist$taxon_name <- as.factor(campus_checklist$taxon_name)
campus_checklist$status <- as.factor(campus_checklist$status)

# Calculate species counts per campus and status
species_counts <- campus_checklist %>%
  group_by(univ.links.uni.abbrev02, status) %>%
  summarise(species_count = n_distinct(taxon_name)) %>%
  ungroup()

# Calculate native, invasive and alien species counts
native_count <- species_counts %>%
  filter(status == "native") %>%
  select(univ.links.uni.abbrev02, native_count = species_count)

invasive_count <- species_counts %>%
  filter(status == "invasive") %>%
  select(univ.links.uni.abbrev02, invasive_count = species_count)

alien_count <- species_counts %>%
  filter(status == "alien") %>%
  select(univ.links.uni.abbrev02, non_invasive_count = species_count)

# Combine all species count columns
species_counts_combined <- native_count %>%
  left_join(invasive_count, by = "univ.links.uni.abbrev02") %>%
  left_join(alien_count, by = "univ.links.uni.abbrev02") %>%
  mutate(
    total_count = native_count + invasive_count + non_invasive_count
  )

# Merge species counts with driving factors data
campus_sf <- left_join(driving_factors, species_counts_combined, by = c("campus02" = "univ.links.uni.abbrev02"))
campus_points_sf <- st_as_sf(campus_sf, coords = c("lng", "lat"), crs = 4326)

write.csv(campus_sf, "./results/campus_plant_num.csv", row.names = FALSE)

set.seed(123)

# Define visualization parameters for size legend
size_breaks <- c(20, 40, 60, 80, 100, 200, 500, 800, 1200, 1500)
size_labels <- c("20", "40", "60", "80", "100", "200", "500", "800", "1200", "1500")
size_range  <- c(1, 6) 
# Note: Actual data range is 1-1424 species

# Create standalone legend for size scale
legend_plot_size <- ggplot(campus_points_sf, aes(x = 1, y = 1, size = total_count)) +
  geom_sf(color = "#4D3BA3", alpha = 0.8) +
  scale_size_continuous(
    range  = size_range,
    breaks = size_breaks,
    labels = size_labels,
    limits = c(0, 1500),
    guide  = guide_legend(
    title.position = "top",
    override.aes = list(color = "#4D3BA3", alpha = 0.8),
    trans="sqrt"
    )
  ) +
  theme_minimal(base_family = "serif") +
  labs(size = "Number of species") +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 16, family = "serif"),
    legend.text  = element_text(size = 12, family = "serif"),
    panel.grid   = element_blank(),
    axis.title   = element_blank(),
    axis.text    = element_blank(),
    axis.ticks   = element_blank()
  ) + 
    guides(size = guide_legend(nrow = 1, byrow = FALSE))

legend_size_grob <- ggpubr::get_legend(legend_plot_size)

# Function to create map with cropped view and inset
plot_map_crop <- function(count_var, plot_title) {
  # Main map plot
  main_plot <- ggplot() +
    geom_sf(data = china_map_province, fill = "white", color = "black", linewidth = 0.2) +
    geom_sf(
      data = campus_points_sf,
      aes_string(size = count_var),
      color = "#4D3BA3",
      shape = 16,
      alpha = 0.7,
      show.legend = FALSE
    ) +
    geom_sf(data = china_map_southsea, color = "black", size = 0.4, linetype = "dashed") +
    scale_size_continuous(
      range  = size_range,
      breaks = size_breaks,
      labels = size_labels,
      limits = c(0, 1500),
      trans = "sqrt"
    ) +
    coord_sf(
      ylim = c(
        st_bbox(china_map_province)["ymin"] + 0.25 * (st_bbox(china_map_province)["ymax"] - st_bbox(china_map_province)["ymin"]),
        st_bbox(china_map_province)["ymax"]
      )
    ) +
    theme_minimal(base_family = "serif") +
    theme(
      legend.position = "none",
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(size = 20, family = "serif", hjust = 0.5)
    ) +
    labs(title = plot_title, x = NULL, y = NULL)

  # Inset map for detailed view
  inset_plot <- ggplot() +
    geom_sf(data = china_map_province, fill = "white", color = "black", linewidth = 0.2) +
    geom_sf(
      data = campus_points_sf,
      aes_string(size = count_var),
      color = "#4D3BA3",
      show.legend = FALSE,
      shape = 16,
      alpha = 0.7
    ) +
    geom_sf(data = china_map_southsea, color = "black", size = 0.4, linetype = "dashed") +
    scale_size_continuous(
      range  = size_range,
      breaks = size_breaks,
      labels = size_labels,
      limits = c(0, 1500),
      trans = "sqrt"
    ) +
    coord_sf(
      xlim = c(105, 125),
      ylim = c(2, 20)
    ) +
    theme_void() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4)
    )
  
  ggdraw() +
    draw_plot(main_plot, 0, 0, 1, 1) +
    draw_plot(inset_plot, 0.75, 0.02, 0.24, 0.2)
}

# Define variables and corresponding plot titles
count_vars <- c(
  "total_count",
  "native_count", 
  "non_invasive_count",
  "invasive_count" 
)
plot_titles <- c("(A) All species", 
               "(B) Native species", 
                 "(C) Non-invasive species",
                 "(D) Invasive species"
                 )

# Generate all four maps using purrr::map2
plots <- purrr::map2(count_vars, plot_titles, plot_map_crop)

# Combine the four maps into a 2x2 grid
final_map <- cowplot::plot_grid(plotlist = plots, nrow = 2, align = "h")

# Combine maps with the size legend
combined_plot <- plot_grid(
  final_map,
  legend_size_grob,
  ncol = 1,        
  rel_heights = c(1, 0.15) 
)

# Export the final combined plot as high-resolution PNG
ggexport(combined_plot, filename = "./results/species_number_map.png",
         width = 3600,
         height = 3000,
         pointsize = 12,
         res = 300)

