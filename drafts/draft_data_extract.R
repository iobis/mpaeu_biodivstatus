
library(terra)
library(mgcv)
library(xgboost)
library(ggplot2)
fs::dir_create("proc-data")
# PUT THE GRIDS ON THAT FOLDER THAT WAS CREATED

europe <- vect("data/mpa_europe_starea_v2.shp") # The path to the study area
env_layers <- list.files("../mpaeu_sdm/data/env/current", full.names = T) # The path to the folder
# with the current period layers
env_layers <- env_layers[grepl("depthsurf", env_layers)]
env_layers <- env_layers[grepl("_mean", env_layers)]
env_layers <- env_layers[grepl("chl|no3|o2|par|siconc|so|tas|thetao|sws", env_layers)]

# The path to the folder with the terrain variables
terrain_layers <- list.files("../mpaeu_sdm/data/env/terrain", full.names = T)
terrain_layers <- terrain_layers[grepl("bathymetry_mean|distcoast|rugosity|wavefetch", terrain_layers)]
terrain_layers <- terrain_layers[!grepl("aux", terrain_layers)]

env_layers <- rast(c(env_layers, terrain_layers))

env_layers <- mask(crop(env_layers, europe), europe)

# Load grid
grid_6 <- vect("proc-data/grid_6.gpkg") #  Comment if you don't want grid_6
grid_5 <- vect("proc-data/grid_5.gpkg")

# Load richness
target_group <- "Aves"
# Ensure that all the richness files are on the folder "data/richness"
richness <- rast(paste0("data/richness/richness_", target_group, "_rf_202411.tif"))

# Start with non-grided version
richness_masked <- mask(richness, europe)
richness_masked <- crop(richness_masked, europe)
names(richness_masked) <- "richness"

# Extract by grid
richness_6 <- terra::extract(richness_masked, grid_6, fun = max, na.rm = T) #  Comment if you don't want grid_6
richness_5 <- terra::extract(richness_masked, grid_5, fun = max, na.rm = T)

env_6 <- terra::extract(env_layers, grid_6, fun = mean, na.rm = T) #  Comment if you don't want grid_6
env_5 <- terra::extract(env_layers, grid_5, fun = mean, na.rm = T)

grid_6 <- cbind(cbind(grid_6, data.frame(richness = richness_6[,-1])), env_6[,-1]) #  Comment if you don't want grid_6
grid_5 <- cbind(cbind(grid_5, data.frame(richness = richness_5[,-1])), env_5[,-1])

centroids_6 <- terra::centroids(grid_6) #  Comment if you don't want grid_6
centroids_5 <- terra::centroids(grid_5)

points_6 <- as.data.frame(centroids_6, geom = "XY") #  Comment if you don't want grid_6
points_5 <- as.data.frame(centroids_5, geom = "XY")

points_6 <- points_6[!is.na(points_6$richness),] #  Comment if you don't want grid_6
points_5 <- points_5[!is.na(points_5$richness),]

valid_data <- as.data.frame(c(richness_masked, env_layers), cell = T, xy = T)
valid_data <- valid_data[!is.na(valid_data$richness),]

varnames <- colnames(valid_data)[5:ncol(valid_data)]

valid_data[,varnames] <- apply(valid_data[,varnames], 2, round, digits = 2)
points_6[,varnames] <- apply(points_6[,varnames], 2, round, digits = 2)
points_5[,varnames] <- apply(points_5[,varnames], 2, round, digits = 2)

write.csv(valid_data, paste0("proc-data/", target_group, "native_data_pts.csv"), row.names = F)
write.csv(points_6, paste0("proc-data/", target_group, "grid6_data_pts.csv"), row.names = F)
write.csv(points_5, paste0("proc-data/", target_group, "grid5_data_pts.csv"), row.names = F)