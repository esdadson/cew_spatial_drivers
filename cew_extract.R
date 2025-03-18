# using renv package to update dependencies
# do not need to run everytime
# can load the packages using library() function (e.g. library(terra))
# renv::install(packages = c("landscapemetrics", "terra"), repos = 'https://rspatial.r-universe.dev')
# renv::snapshot()
# renv::update()
library(terra)
library(landscapemetrics)
library(tidyverse)
library(dplyr)
library(readr)

# load trapping data into the workspace
data_dir <- "/Volumes/cmjone25/Data/Original/pest-occurrence/corn_earworm/"
output_dir <- "/Volumes/cmjone25/Data/Raster/USA/pops_casestudies/corn_earworm/outputs"
trap_data <- read.csv(paste(data_dir, 
                            "trap_network_2020_21_22_23_24_counts.csv", 
                            sep = ""))

trap_data$CEW_count <- as.numeric(trap_data$CEW_count)
trap_data <- trap_data %>%
  select(trapID, year, county, lat, long, CEW_count) %>%
  group_by(trapID, year) %>%
  summarize(cum_CEW_count = sum(CEW_count, na.rm = T),
            lat = lat[1],
            long = long[2],
            county = county[1])

# View first 6 rows
head(trap_data)

# What is trap_data
class(trap_data)

# create vector of trap data and use plot to see distribution
trap_vect <- terra::vect(trap_data, 
                         geom = c("long", "lat"), 
                         crs = "EPSG:4326")

library(tidyverse)
plot(trap_vect)
maps::map.text('county', 'north carolina', add = TRUE)
maps::map('county', 'north carolina', fill = FALSE, add = TRUE)

# read cropland data layer (CDL) as raster files
cdl_path <- "/Volumes/cmjone25/Data/Raster/USA/landcover/"
cdl_data <- rast(paste0(cdl_path, "2020_30m_cdls/2020_30m_cdls.tif"))

# the CDL data is large, we need to crop the data to our study area, but first 
# we need to ensure the data are in the same coordinate reference system (CRS). 
# This is how we are able to preserve the integrity of operations on spatial data. 
# We first insure the spatial data are speaking in the same language (CRS), then 
# we can perform operations.
e <- terra::ext(trap_vect)
e[1] <- -78.5
e[2] <- -77
e[3] <- 35
e[4] <- 37
empty_vect <- vect(e, crs = "EPSG:4326")
empty_vect <- project(empty_vect, crs(cdl_data))
trap_vect <- project(trap_vect, crs(cdl_data))
cdl_data_cropped <- crop(cdl_data, empty_vect)
plot(cdl_data_cropped)
plot(trap_vect, add = T)

# list environmental variables in your workspace
ls()
# remove unneeded large files from your workspace
rm("cdl_data_cropped")
rm("cdl_data")
rm("e")

# example of a for loop extracting data from each year's cropland data layer
for (year in seq(2020, 2023)) {
  cdl_data <- rast(paste0(cdl_path, year, "_30m_cdls/", year, "_30m_cdls.tif"))
  cdl_data <- crop(cdl_data, empty_vect)
  trap_year <- trap_vect[trap_vect$year == year,]
  trap_buffer <- buffer(trap_year, 500)
  trap_extract <- terra::extract(cdl_data, trap_buffer)
  trap_year$percent_corn <- table(trap_extract)[,c("Corn")]/nrow(trap_extract)
  trap_year$percent_cotton <- table(trap_extract)[,c("Cotton")]/nrow(trap_extract)
  trap_year$percent_peanuts <- table(trap_extract)[,c("Peanuts")]/nrow(trap_extract)
  trap_year$percent_soybeans <- table(trap_extract)[,c("Soybeans")]/nrow(trap_extract)
  trap_year$percent_sorghum <- table(trap_extract)[,c("Sorghum")]/nrow(trap_extract)
  trap_year$percent_sweetcorn <- table(trap_extract)[,c("Sweet Corn")]/nrow(trap_extract)
  trap_year$percent_tobacco <- table(trap_extract)[,c("Tobacco")]/nrow(trap_extract)
  trap_year$egg_area <- rowSums(cbind(table(trap_extract)[,c("Corn")],
                                  table(trap_extract)[,c("Cotton")],
                                  table(trap_extract)[,c("Peanuts")],
                                  table(trap_extract)[,c("Soybeans")],
                                  table(trap_extract)[,c("Sorghum")],
                                  table(trap_extract)[,c("Sweet Corn")],
                                  table(trap_extract)[,c("Tobacco")]))
  assign(paste0("trap_", year), trap_year)
}

trap_vect <- rbind(trap_2020, trap_2021, trap_2022, trap_2023)
trap_vect <- project(trap_vect, "EPSG:4326")
trap_data <- cbind(as.data.frame(trap_vect), long = crds(trap_vect)[,1],
                   lat = crds(trap_vect)[,2])
write.csv(trap_data, file = paste0(output_dir, "/traps_extracted_data.csv"))

# example extracting landscape information
trap_point <- trap_data[1,]
trap_geo <- terra::vect(trap_point, geom = c("long", "lat"), crs = "EPSG:4326")
trap_geo <- terra::project(trap_geo, crs(cdl_2020))
trap_buffer <- terra::buffer(trap_geo, 1000)
trap_extract <- terra::extract(cdl_2020, trap_buffer)
percent_corn_1k <- table(trap_extract)[,c("Corn")]/nrow(trap_extract)
percent_primary_1k <- sum(table(trap_extract)[,c("Corn","Cotton")])/nrow(trap_extract)

# New nested for loop for year & buffer combined

buffer_lengths <- c(200, 400, 600, 800, 1000, 1200)

all_data <- list()

for (year in 2020:2023) {
  for (buffer in buffer_lengths) {
    
    cdl_data <- rast(paste0(cdl_path, year, "_30m_cdls/", year, "_30m_cdls.tif"))
    cdl_data <- crop(cdl_data, empty_vect)
    
    trap_year <- trap_vect[trap_vect$year == year,]
    
    trap_buffer <- buffer(trap_year, width = buffer)
    
    trap_extract <- terra::extract(cdl_data, trap_buffer)
    
    crop_list <- c("Corn", "Cotton", "Peanuts", "Soybeans", "Sorghum", "Sweet Corn", "Tobacco")
    
    for (crop in crop_list) {
      if (crop %in% colnames(table(trap_extract))) {
        trap_year[[paste0("percent_", tolower(gsub(" ", "", crop)))]] <- table(trap_extract)[,crop] / nrow(trap_extract)
      } else {
        trap_year[[paste0("percent_", tolower(gsub(" ", "", crop)))]] <- 0
      }
    }
    
    trap_year$egg_area <- sum(table(trap_extract)[,c("Corn")],
                              table(trap_extract)[,c("Cotton")],
                              table(trap_extract)[,c("Peanuts")],
                              table(trap_extract)[,c("Soybeans")],
                              table(trap_extract)[,c("Sorghum")],
                              table(trap_extract)[,c("Sweet Corn")],
                              table(trap_extract)[,c("Tobacco")])
    trap_year$year <- year
    trap_year$buffer <- buffer
    
    trap_year <- project(trap_year, "EPSG:4326")
    trap_data <- cbind(
      as.data.frame(trap_year),
      long = crds(trap_year)[,1],
      lat = crds(trap_year)[,2]
    )
    all_data[[paste0(year, "_", buffer)]] <- trap_data
  }
}
combined_data <- do.call(rbind, all_data)
write_csv(combined_data, file = paste0(output_dir, "/traps_extracted_data.csv"))
