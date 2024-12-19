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
