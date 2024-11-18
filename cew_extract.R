install.packages('terra', repos='https://rspatial.r-universe.dev')
install.packages("landscapemetrics")
library(landscapemetrics)
library(terra)

data_dir <- "/Volumes/cmjone25/Data/Original/pest-occurrence/corn_earworm/"
trap_data <- read.csv(paste(data_dir, 
                            "trap_network_2020_21_22_23_24_counts.csv", 
                            sep = ""))
head(trap_data)

# create vector of trap data and use plot to see distribution
trap_vect <- vect(trap_data, geom = c("long", "lat"), crs = "EPSG:4326")
plot(trap_vect)

# read cropland data layer (CDL) as raster files
cdl_path <- "/Volumes/cmjone25/Data/Raster/USA/landcover/"
cdl_2020 <- rast(paste0(cdl_path, "2020_30m_cdls/2020_30m_cdls.tif"))

# the CDL data is large, we need to crop the data to our study area, but first 
# we need to ensure the data are in the same coordinate reference system (CRS). 
# This is how we are able to preserve the integrity of operations on spatial data. 
# We first insure the spatial data are speaking in the same language (CRS), then 
# we can perform operations.
trap_vect <- project(trap_vect, crs(cdl_2020))
cdl_2020_cropped <- crop(cdl_2020, trap_vect)
plot(cdl_2020_cropped)
plot(trap_vect, add = T)

> extract_metrics <- sample_lsm(cdl_2020,
                                traps_2020,
                                size = c(1000),
                                metric = 'area', 
                                level = "class",
                                plot_id = traps_2020$trapID,
                                shape = 'circle')
