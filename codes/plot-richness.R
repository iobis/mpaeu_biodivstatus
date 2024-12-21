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
library(patchwork)
sf::sf_use_s2(FALSE)
fs::dir_create("figures")

europe <- vect("data/mpa_europe_starea_v2.shp")

ric_preds <- list.files("data/richness", full.names = T)
ric_preds <- ric_preds[grepl("_rf_", ric_preds)]

titles <- gsub("data/richness/richness_", "", ric_preds)
titles <- gsub("_rf_202411.tif", "", titles)
titles <- gsub("full", "All groups", titles)

wrld <- rnaturalearth::ne_countries(returnclass = "sf", scale = "small")
wrld <- st_crop(wrld, europe)

europe_buff <- sf::st_as_sf(buffer(europe, width = 1000))

for (i in seq_along(ric_preds)) {

    r <- rast(ric_preds[i])
    r <- crop(mask(r, europe), europe)

    r <- as.data.frame(r, xy = T)
    colnames(r)[3] <- "Richness"

    p1 <- ggplot() +
        geom_sf(data = wrld, fill = "grey90", color = "grey90") +
        geom_sf(data = europe_buff, fill = "grey90", color = "grey90") +
        geom_raster(data = r, aes(x = x, y = y, fill = Richness)) +
        scale_fill_distiller(palette = "Blues", direction = 1) +
        theme_light() +
        xlab(NULL) + ylab(NULL) +
        ggtitle(titles[i]) +
        theme(panel.border = element_blank())
        #coord_sf(crs = "EPSG:3035")

    ggsave(filename = file.path("figures", paste0("richness-", titles[i], ".jpg")),
    plot = p1, quality = 100, width = 11, height = 7)
}

# For the full richness map, produce also a detailed view

full_richness <- rast("data/richness/richness_full_rf_202411.tif")
full_richness <- crop(mask(full_richness, europe), europe)

# From marineregions.org, pre-processed on QGIS to remove areas out of our region
iho <- vect("data/World_Seas_IHO_v3.shp")

full_richness_df <- as.data.frame(full_richness, xy = T)
colnames(full_richness_df)[3] <- "Richness"

p1 <- ggplot() +
    geom_sf(data = wrld, fill = "grey90", color = "grey90") +
    geom_sf(data = europe_buff, fill = "grey90", color = "grey90") +
    geom_raster(data = full_richness_df, aes(x = x, y = y, fill = Richness)) +
    scale_fill_distiller(palette = "Blues", direction = 1) +
    theme_light() +
    xlab(NULL) + ylab(NULL) +
    ggtitle("All groups") +
    theme(panel.border = element_blank())

# Baltic
baltic <- crop(full_richness, iho[iho$mpa_region == "baltic",])
baltic <- as.data.frame(baltic, xy = T)
colnames(baltic)[3] <- "Richness"

# Black sea
black_sea <- crop(full_richness, iho[iho$mpa_region == "black_sea",])
black_sea <- as.data.frame(black_sea, xy = T)
colnames(black_sea)[3] <- "Richness"

# Med
med <- crop(full_richness, iho[iho$mpa_region == "mediterranean",])
med <- as.data.frame(med, xy = T)
colnames(med)[3] <- "Richness"

get_lims <- function(what = "x", where = "baltic") {
    sl <- iho[iho$mpa_region == where,]
    sl <- ext(sl)
    if (what == "x") {
        unname(sl[1:2])
    } else {
        unname(sl[3:4])
    }
}

p2 <- ggplot() +
    geom_sf(data = wrld, fill = "grey90", color = "grey90") +
    geom_sf(data = europe_buff, fill = "grey90", color = "grey90") +
    geom_raster(data = baltic, aes(x = x, y = y, fill = Richness)) +
    scale_fill_distiller(palette = "Blues", direction = 1) +
    theme_light() +
    xlab(NULL) + ylab(NULL) +
    ggtitle("Baltic Sea") +
    coord_sf(xlim = get_lims("x", "baltic"), ylim = get_lims("y", "baltic")) +
    theme(panel.border = element_blank());p2

p3 <- ggplot() +
    geom_sf(data = wrld, fill = "grey90", color = "grey90") +
    geom_sf(data = europe_buff, fill = "grey90", color = "grey90") +
    geom_raster(data = black_sea, aes(x = x, y = y, fill = Richness)) +
    scale_fill_distiller(palette = "Blues", direction = 1) +
    theme_light() +
    xlab(NULL) + ylab(NULL) +
    ggtitle("Black Sea") +
    coord_sf(xlim = get_lims("x", "black_sea"), ylim = get_lims("y", "black_sea")) +
    theme(panel.border = element_blank());p3

p4 <- ggplot() +
    geom_sf(data = wrld, fill = "grey90", color = "grey90") +
    geom_sf(data = europe_buff, fill = "grey90", color = "grey90") +
    geom_raster(data = med, aes(x = x, y = y, fill = Richness)) +
    scale_fill_distiller(palette = "Blues", direction = 1) +
    theme_light() +
    xlab(NULL) + ylab(NULL) +
    ggtitle("Mediterranean Sea") +
    coord_sf(xlim = get_lims("x", "mediterranean"), ylim = get_lims("y", "mediterranean")) +
    theme(panel.border = element_blank());p4

(p1 & theme(legend.position = "bottom")) + ((p2/p3/p4) & theme(
    legend.position = "none", axis.text = element_blank(),
    axis.ticks = element_blank()
))

ggsave(filename = file.path("figures", paste0("richness-detailed.jpg")),
 quality = 100, width = 9, height = 8)
