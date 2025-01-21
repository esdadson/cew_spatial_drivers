# using renv package to update dependencies
# do not need to run everytime
# can load the packages using library() function (e.g. library(terra))
# renv::install(packages = c("landscapemetrics", "terra"), repos = 'https://rspatial.r-universe.dev')
# renv::snapshot()
# renv::update()
library(terra)
library(landscapemetrics)

# load trapping data into the workspace
data_dir <- "/Volumes/cmjone25/Data/Original/pest-occurrence/corn_earworm/"
output_dir <- "/Volumes/cmjone25/Data/Raster/USA/pops_casestudies/corn_earworm/outputs"
trap_data <- read.csv(paste(data_dir, 
                            "trap_network_2020_21_22_23_24_counts.csv", 
                            sep = ""))

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
  assign(paste0("trap_", year), trap_year)
}

trap_vect <- rbind(trap_2020, trap_2021, trap_2022, trap_2023)
trap_vect <- project(trap_vect, "EPSG:4326")
trap_data <- cbind(as.data.frame(trap_vect), long = crds(trap_vect)[,1],
                   lat = crds(trap_vect)[,2])
write_csv(trap_data, file = paste0(output_dir, "/traps_extracted_data.csv"))

# example extracting landscape information
trap_point <- trap_data[1,]
trap_geo <- terra::vect(trap_point, geom = c("long", "lat"), crs = "EPSG:4326")
trap_geo <- terra::project(trap_geo, crs(cdl_2020))
trap_buffer <- terra::buffer(trap_geo, 1000)
trap_extract <- terra::extract(cdl_2020, trap_buffer)
percent_corn_1k <- table(trap_extract)[,c("Corn")]/nrow(trap_extract)
percent_primary_1k <- sum(table(trap_extract)[,c("Corn","Cotton")])/nrow(trap_extract)
