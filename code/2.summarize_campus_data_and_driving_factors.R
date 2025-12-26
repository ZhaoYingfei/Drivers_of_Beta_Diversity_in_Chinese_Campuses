# ====================================================
# ====Summarize campus data and environmental data====
# ====================================================
cat("\014")  # Clear console
rm(list = ls())  # Clear environment

# Load required packages
library(raster)   # for handling raster data
library(sf)       # for reading and manipulating shapefiles
library(dplyr)    # for data manipulation
library(readxl)   # for reading Excel files

setwd("C:/Users/qianh/Desktop/R/cam.pla.div")
getwd()  # Confirm working directory

# ------------------------------------------
# Read and prepare the shapefile of Chinese cities
# ------------------------------------------
city.map <- sf::st_read(dsn = "./data/City_map/CN-shi-A.shp") %>%
  filter(!is.na(CityNameC)) %>%
  st_transform(4326)  # Convert projection to WGS84 (lat/lon)

# Merge polygons of the same city to avoid duplicates
city.map <- city.map %>%
  group_by(CityNameC) %>%
  summarise(SHI_D_ID = first(SHI_D_ID), .groups = "drop")

nrow(city.map)  # Number of unique cities

# ------------------------------------------
# Load national-scale environmental raster layers
# ------------------------------------------
dem <- raster("./data/Data_extracted/dem20240826.tif")  # Elevation (DEM)
mat <- raster("./data/Data_extracted/CHELSA_bio1_1981-2010_V.2.1.tif")  # Mean annual temperature
map <- raster("./data/Data_extracted/CHELSA_bio12_1981-2010_V.2.1.tif")  # Mean annual precipitation

# Check coordinate reference systems
crs(city.map)  # Should be "+proj=longlat +datum=WGS84 ..."
crs(dem)
crs(mat)
crs(map)

# ------------------------------------------
# Loop over each city to extract raster values
# ------------------------------------------
results <- list()  # Initialize result list

for (i in 1:nrow(city.map)) {
  city_i <- city.map[i, ]
  city_name <- city_i$CityNameC
  message(paste0("Processing city: ", city_name, " (", i, "/", nrow(city.map), ")"))
  
  # Convert sf object to Spatial object for raster masking
  city_sp <- as(city_i, "Spatial")
  
  # Crop and mask rasters by city polygon
  dem_mask <- mask(crop(dem, city_sp), city_sp)
  map_mask <- mask(crop(map, city_sp), city_sp)
  mat_mask <- mask(crop(mat, city_sp), city_sp)
  
  # Calculate mean values for each raster
  dem_mean <- cellStats(dem_mask, stat = "mean", na.rm = TRUE)
  map_mean <- cellStats(map_mask, stat = "mean", na.rm = TRUE)
  mat_mean <- cellStats(mat_mask, stat = "mean", na.rm = TRUE)
  
  # Store results in list
  results[[i]] <- data.frame(
    i = i,
    city_name = city_name,
    SHI_D_ID = city_i$SHI_D_ID,
    dem_mean = dem_mean,
    map_mean = map_mean,
    mat_mean = mat_mean
  )
}

# Combine all results into a single data frame
results_df <- bind_rows(results)

# ------------------------------------------
# Merge with campus data
# ------------------------------------------
campus_data <- read_excel("./data/Raw_data/campus_data_all.xlsx", sheet = "campus_data")
# wetland_data <- read_csv("./data/Wetland_data.csv",locale = locale(encoding = "GBK"))
# Join campus-level data with extracted environmental variables
campus_merged1 <- left_join(campus_data, results_df, by = "city_name")
# campus_merged2 <- campus_merged1 %>% left_join(wetland_data %>% select(city_name, wet_area, wet_num), by = "city_name")

# Remove universities 'fjsfu', 'gxsfu', 'znlykju', and 'zust' from the driving factors dataset
# Fujian Normal University, 
# Guangxi Normal University, 
# Central South University of Forestry and Technology, 
# and Zhejiang University of Science and Technology, do not have invasive plants
campus_merged3 <- campus_merged1 %>%
  filter(!campus02 %in% c("fjsfu", "gxsfu", "znlykju", "zust"))

# Export final result
write.csv(campus_merged3, "./data/dem_climate_mean.csv", row.names = FALSE)

# ====================================================
# ===============Compile GDP data=====================
# ====================================================
cat("\014")  # Clear console
rm(list = ls())  # Clear environment

# Load required packages
library(tidyverse)   
library(sf)       
library(terra)    
library(exactextractr)  
library(readr)  
library(readxl)

setwd("C:/Users/qianh/Desktop/R/cam.pla.div")
getwd()  

# Define paths to ArcInfo Grid raster directories (containing .adf files)
gdp2000_dir <- "./data/GDP_data/gdp2000"
gdp2005_dir <- "./data/GDP_data/gdp2005"
gdp2010_dir <- "./data/GDP_data/gdp2010"
gdp2015_dir <- "./data/GDP_data/gdp2015"
gdp2020_dir <- "./data/GDP_data/gdp2020"

# Read the list of cities to be processed
cities_need <- readr::read_csv("./data/dem_climate_mean.csv", show_col_types = FALSE) |>
  dplyr::select(city_name) |>
  distinct() |>
  mutate(city_name = trimws(city_name))

stopifnot(nrow(cities_need) > 0)

# Read shapefile of Chinese cities
city.map <- sf::st_read(dsn = "./data/City_map/CN-shi-A.shp") %>%
  dplyr::filter(!is.na(CityNameC)) |>
  st_transform(4326) |>
  group_by(CityNameC) |>
  summarise(SHI_D_ID = first(SHI_D_ID), .groups = "drop")

# Keep only the required cities and validate geometry
city_poly <- city.map |>
  dplyr::rename(city_name = CityNameC) |>
  dplyr::semi_join(cities_need, by = "city_name") |>
  sf::st_make_valid()

# Read raster data for five years
r2000 <- terra::rast(gdp2000_dir)
r2005 <- terra::rast(gdp2005_dir)
r2010 <- terra::rast(gdp2010_dir)
r2015 <- terra::rast(gdp2015_dir)
r2020 <- terra::rast(gdp2020_dir)

# Reproject all rasters to match 2020 grid
r2000 <- terra::project(r2000, r2020)
r2005 <- terra::project(r2005, r2020)
r2010 <- terra::project(r2010, r2020)
r2015 <- terra::project(r2015, r2020)

# Transform city polygons to the CRS of the rasters
city_poly_r <- sf::st_transform(city_poly, terra::crs(r2020) |> as.character())

# Extract mean GDP values for each polygon from rasters
extract_mean <- function(r, poly) exactextractr::exact_extract(r, poly, "mean")

df_gdp <- city_poly_r |>
  mutate(
    gdp_2000 = extract_mean(r2000, city_poly_r),
    gdp_2005 = extract_mean(r2005, city_poly_r),
    gdp_2010 = extract_mean(r2010, city_poly_r),
    gdp_2015 = extract_mean(r2015, city_poly_r),
    gdp_2020 = extract_mean(r2020, city_poly_r)
  ) |>
  st_drop_geometry() |>
  dplyr::select(city_name, gdp_2000, gdp_2005, gdp_2010, gdp_2015, gdp_2020)

# Merge extracted GDP with city list
df_gdp <- cities_need |>
  left_join(df_gdp, by = "city_name")

# Read CPI data
cpi_data <- read_excel("./data/GDP_data/CPI_datatable_2023.xlsx", sheet = "datatable")

# Keep average CPI (AVE) for the required years
cpi_avg <- cpi_data %>% 
  dplyr::select(YEAR, AVE) %>% 
  filter(YEAR %in% c(2000, 2005, 2010, 2015, 2020)) %>%
  mutate(AVE = as.numeric(AVE)) 

# Convert GDP data to long format
df_gdp_long <- df_gdp %>%
  pivot_longer(cols = starts_with("gdp_"), 
               names_to = "year", 
               names_prefix = "gdp_", 
               values_to = "gdp_value") 

# Join CPI data by year
df_gdp_long <- df_gdp_long %>%
  left_join(cpi_avg, by = c("year" = "YEAR"))

# Adjust GDP based on 2000 CPI
df_gdp_long <- df_gdp_long %>%
  mutate(gdp_adjusted = gdp_value / AVE * cpi_avg$AVE[cpi_avg$YEAR == 2000])

# Compute mean adjusted GDP for each city in five years
df_gdp_long_avg <- df_gdp_long %>%
  group_by(city_name) %>%
  mutate(gdp_mean = mean(gdp_adjusted, na.rm = TRUE)) %>%
  ungroup()
# Keep only city_name and gdp_mean, remove duplicates
df_gdp_mean <- df_gdp_long_avg %>%
  dplyr::select(city_name, gdp_mean) %>%  
  distinct(city_name, .keep_all = TRUE) 

# Save result to CSV
write.csv(df_gdp_mean, "./data/gdp_mean.csv", row.names = FALSE)

# ====================================================
# ==========compile population density data===========
# ====================================================
cat("\014")  # Clear console
rm(list = ls())  # Clear environment

# Load required packages
library(tidyverse)   
library(sf)       
library(terra)    
library(exactextractr)  
library(readr)  
library(readxl)

setwd("C:/Users/qianh/Desktop/R/cam.pla.div")
getwd()  

# Define paths to ArcInfo Grid raster directories (containing .adf files)
pd2000_dir <- "./data/pd_data/pd2000"
pd2005_dir <- "./data/pd_data/pd2005"
pd2010_dir <- "./data/pd_data/pd2010"
pd2015_dir <- "./data/pd_data/pd2015"
pd2020_dir <- "./data/pd_data/pd2020"

# Read the list of cities to be processed
cities_need <- readr::read_csv("./data/dem_climate_mean.csv", show_col_types = FALSE) |>
  dplyr::select(city_name) |>
  distinct() |>
  mutate(city_name = trimws(city_name))

stopifnot(nrow(cities_need) > 0)

# Read shapefile of Chinese cities
city.map <- sf::st_read(dsn = "./data/City_map/CN-shi-A.shp") %>%
  dplyr::filter(!is.na(CityNameC)) |>
  st_transform(4326) |>
  group_by(CityNameC) |>
  summarise(SHI_D_ID = first(SHI_D_ID), .groups = "drop")

# Keep only the required cities and validate geometry
city_poly <- city.map |>
  dplyr::rename(city_name = CityNameC) |>
  dplyr::semi_join(cities_need, by = "city_name") |>
  sf::st_make_valid()

# Read raster data for five years
r2000 <- terra::rast(pd2000_dir)
r2005 <- terra::rast(pd2005_dir)
r2010 <- terra::rast(pd2010_dir)
r2015 <- terra::rast(pd2015_dir)
r2020 <- terra::rast(pd2020_dir)

# Reproject all rasters to match 2020 grid
r2000 <- terra::project(r2000, r2020)
r2005 <- terra::project(r2005, r2020)
r2010 <- terra::project(r2010, r2020)
r2015 <- terra::project(r2015, r2020)

# Transform city polygons to the CRS of the rasters
city_poly_r <- sf::st_transform(city_poly, terra::crs(r2020) |> as.character())

# Extract mean population density values for each polygon from rasters
extract_mean <- function(r, poly) exactextractr::exact_extract(r, poly, "mean")

df_pd <- city_poly_r |>
  mutate(
    pd_2000 = extract_mean(r2000, city_poly_r),
    pd_2005 = extract_mean(r2005, city_poly_r),
    pd_2010 = extract_mean(r2010, city_poly_r),
    pd_2015 = extract_mean(r2015, city_poly_r),
    pd_2020 = extract_mean(r2020, city_poly_r)
  ) |>
  st_drop_geometry() |>
  dplyr::select(city_name, pd_2000, pd_2005, pd_2010, pd_2015, pd_2020)

# Merge extracted population density with city list
df_pd <- cities_need |>
  left_join(df_pd, by = "city_name")

# Compute mean population density for each city in five years and keep only city_name and gdp_mean
df_pd_mean <- df_pd %>%
  rowwise() %>%
  mutate(pd_mean = mean(c_across(starts_with("pd_")), na.rm = TRUE)) %>%
  ungroup()%>%
  dplyr::select(city_name, pd_mean)

# Save result to CSV
write_csv(df_pd_mean, "./data/pd_mean.csv")

# ====================================================
# ============compile urbanization level==============
# ====================================================
cat("\014")  # Clear console
rm(list = ls())  # Clear environment

# Load required packages
library(tidyverse)   
library(sf)       
library(terra)    
library(readr)  
library(readxl)
library(purrr)

setwd("C:/Users/qianh/Desktop/R/cam.pla.div")
getwd()  

# Read the list of cities to be processed
cities_need <- readr::read_csv("./data/dem_climate_mean.csv", show_col_types = FALSE) |>
  dplyr::select(city_name) |>
  distinct() |>
  mutate(city_name = trimws(city_name))

# Read shapefile of Chinese cities
city.map <- sf::st_read(dsn = "./data/City_map/CN-shi-A.shp") %>%
  dplyr::filter(!is.na(CityNameC)) |>
  st_transform(4326) |>
  group_by(CityNameC) |>
  summarise(SHI_D_ID = first(SHI_D_ID), .groups = "drop")

# Keep only the required cities and validate geometry
city_poly <- city.map |>
  dplyr::rename(city_name = CityNameC) |>
  dplyr::semi_join(cities_need, by = "city_name") |>
  sf::st_make_valid()

# Read raster data for each years
ul2000 <- sf::st_read("./data/ul_data/ul2000/GUB_Global_2000.shp")
ul2005 <- sf::st_read("./data/ul_data/ul2005/GUB_Global_2005.shp")
ul2010 <- sf::st_read("./data/ul_data/ul2010/GUB_Global_2010.shp")
ul2015 <- sf::st_read("./data/ul_data/ul2015/GUB_Global_2015.shp")
ul2018 <- sf::st_read("./data/ul_data/ul2018/GUB_Global_2018.shp")

# Disable s2 (GEOS intersection is faster in projected coordinates), 
# use equal-area projection
sf::sf_use_s2(FALSE)

# Use global equal-area projection (World Cylindrical Equal Area)
target_crs <- 6933

# Project to equal-area coordinate system and calculate total area
city_poly_eq <- city_poly |>
  st_make_valid() |>
  st_transform(target_crs) |>
  mutate(total_area = as.numeric(st_area(geometry))) |>  
  dplyr::select(city_name, total_area, geometry)

# Define a function to intersect with UL layer for a given year and calculate urban fraction
urban_fraction_by_year <- function(ul, city_eq, year) {
  ul_eq <- ul |>
    st_make_valid() |>
    st_transform(target_crs)

  # Single intersection (cities × urbanization level layer for that year)
  inter <- st_intersection(city_eq, ul_eq)

  # Aggregate intersected area by city name
  res <- inter |>
    mutate(a = as.numeric(st_area(geometry))) |>
    st_drop_geometry() |>
    group_by(city_name) |>
    summarise(urban_area = sum(a, na.rm = TRUE), .groups = "drop") |>
    right_join(st_drop_geometry(city_eq), by = "city_name") |>
    mutate(!!paste0("urbanization_", year) := coalesce(urban_area, 0) / total_area) |>
    dplyr::select(city_name, !!sym(paste0("urbanization_", year)))

  res
}

# Calculate for five years (one intersection per year), then combine results
yrs <- c(2000, 2005, 2010, 2015, 2018)
ul_list <- list(`2000` = ul2000, `2005` = ul2005, `2010` = ul2010, `2015` = ul2015, `2018` = ul2018)

results <- imap(ul_list, ~ urban_fraction_by_year(.x, city_poly_eq, .y))  # .y是年份字符串
df_urbanization <- reduce(results, left_join, by = "city_name")
# Check for monotonic increase
res <- df_urbanization %>%
  mutate(
    v2000_2005 = urbanization_2000 <= urbanization_2005,
    v2005_2010 = urbanization_2005 <= urbanization_2010,
    v2010_2015 = urbanization_2010 <= urbanization_2015,
    v2015_2018 = urbanization_2015 <= urbanization_2018
  ) %>%
  mutate(across(starts_with("v"), ~coalesce(.x, TRUE))) %>%
  mutate(any_violation = if_any(starts_with("v"), ~ .x == FALSE))

res %>%
  summarise(
    n_rows = n(),
    n_violations = sum(any_violation),
    pct_violations = mean(any_violation) * 100
  )
# List which cities have violations for which year pairs
violations_long <- res %>%
  filter(any_violation) %>%
  pivot_longer(starts_with("v"),
               names_to = "pair", values_to = "ok") %>%
  filter(!ok) %>%
  mutate(pair = dplyr::recode(pair,
    v2000_2005 = "2000 > 2005",
    v2005_2010 = "2005 > 2010",
    v2010_2015 = "2010 > 2015",
    v2015_2018 = "2015 > 2018"
  )) %>%
  dplyr::select(city_name, pair)

# View details
violations_long
# # A tibble: 3 × 2
#   city_name pair       
#   <chr>     <chr>      
# 1 扬州市    2015 > 2018
# 2 淮安市    2015 > 2018
# 3 盐城市    2015 > 2018

# Compute mean population density for each city in five years and keep only city_name and gdp_mean
df_urbanization_mean <- df_urbanization %>%
  rowwise() %>%
  mutate(ul_mean = mean(c_across(starts_with("urbanization_")), na.rm = TRUE)) %>%
  ungroup()%>%
  dplyr::select(city_name, ul_mean)

# Save the result to CSV
write_csv(df_urbanization_mean, "./data/ul_mean.csv")

#——————————————Aggregate all driving factors——————————————
cat("\014")  # Clear console
rm(list = ls())  # Clear environment

library(readr)
library(dplyr)

setwd("C:/Users/qianh/Desktop/R/cam.pla.div")
getwd()  

dem_climate <- read_csv("./data/dem_climate_mean.csv")
gdp <- read_csv("./data/gdp_mean.csv")
pd <- read_csv("./data/pd_mean.csv")
ul <- read_csv("./data/ul_mean.csv")

# Merge step by step by city_name
Driving_factors <- dem_climate %>%
  left_join(gdp, by = "city_name") %>%
  left_join(pd,  by = "city_name") %>%
  left_join(ul,  by = "city_name")

# Save as new CSV
write_csv(Driving_factors, "./data/Driving_factors.csv")
