################## MPA Europe - Biodiversity status exploring ##################
#################################### WP5 #######################################
# November of 2024
# Author: Silas Principe, Anna Adamo
# Contact: s.principe@unesco.org
#
# Plot richness maps
# Richness maps are available at the S3 bucket at https://mpaeu-dist.s3.amazonaws.com/results/diversity/
# Study area is available at https://github.com/iobis/mpaeu_studyarea

library(terra)
library(ggplot2)
library(sf)
library(ggpubr)
sf::sf_use_s2(FALSE)
fs::dir_create("figures")

europe <- vect("data/mpa_europe_starea_v2.shp")

wrld <- rnaturalearth::ne_countries(returnclass = "sf", scale = "small")
wrld <- st_crop(wrld, europe)

europe_buff <- sf::st_as_sf(buffer(europe, width = 1000))

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

names(env_layers) <- c(
    "Chlorophyll", "Nox", "O2", "PAR", "Sea ice", "Salinity", "Current speed", "Air temperature",
    "Sea surface temperature", "Bathymetry", "Distance to coast", "Rugosity", "Wave fetch"
)

env_layers_df <- as.data.frame(env_layers, xy = T)

env_layers_df <- tidyr::pivot_longer(env_layers_df, 3:15, names_to = "variable", values_to = "value")

plot_list <- lapply(seq_len(nlyr(env_layers)), function(x) NULL)

for (i in seq_along(unique(env_layers_df$variable))) {
    tv <- unique(env_layers_df$variable)[i]

    plot_list[[i]] <- ggplot() +
        geom_sf(data = wrld, fill = "grey90", color = "grey90") +
        geom_sf(data = europe_buff, fill = "grey90", color = "grey90") +
        geom_raster(data = env_layers_df[env_layers_df$variable == tv,], aes(x = x, y = y, fill = value)) +
        scale_fill_viridis_c() +
        #scale_fill_distiller(palette = "Blues", direction = 1) +
        theme_light() +
        xlab(NULL) + ylab(NULL) +
        ggtitle(tv) +
        theme(panel.border = element_blank(), legend.title = element_blank())
}

# Add plots of richness on the two remaining spots
richness <- rast("data/richness/richness_full_rf_202411.tif")
richness <- mask(crop(richness, europe), europe)
richness_df <- as.data.frame(richness, xy = T)
colnames(richness_df)[3] <- "Richness"

ric_raw <- read.csv("proc-data/raw_grid6_data_pts_metrics.csv")
ric_raw <- ric_raw[,c("sp", "h3_6")]
ric_raw <- ric_raw[!is.na(ric_raw$sp),]
names(ric_raw)[1] <- "Richness"
ric_raw_pol <- h3jsr::cell_to_polygon(ric_raw$h3_6, simple = F)
names(ric_raw_pol)[1] <- "h3_6"
ric_raw_pol <- dplyr::left_join(ric_raw_pol, ric_raw)

plot_list[[14]] <- ggplot() +
        geom_sf(data = wrld, fill = "grey90", color = "grey90") +
        geom_sf(data = europe_buff, fill = "grey90", color = "grey90") +
        geom_raster(data = richness_df, aes(x = x, y = y, fill = Richness)) +
        scale_fill_viridis_c() +
        #scale_fill_distiller(palette = "Blues", direction = 1) +
        theme_light() +
        xlab(NULL) + ylab(NULL) +
        ggtitle("Richness (modeled)") +
        theme(panel.border = element_blank(), legend.title = element_blank())

plot_list[[15]] <- ggplot() +
        geom_sf(data = wrld, fill = "grey90", color = "grey90") +
        geom_sf(data = europe_buff, fill = "grey90", color = "grey90") +
        geom_sf(data = ric_raw_pol, aes(fill = Richness), color = "#ffffff00") +
        scale_fill_viridis_c() +
        #scale_fill_distiller(palette = "Blues", direction = 1) +
        theme_light() +
        xlab(NULL) + ylab(NULL) +
        ggtitle("Richness (raw)") +
        theme(panel.border = element_blank(), legend.title = element_blank())


ggarrange(plotlist = plot_list, ncol = 3, nrow = 5)

ggsave(filename = file.path("figures", paste0("env-variables.jpg")),
    quality = 100, width = 10, height = 15)