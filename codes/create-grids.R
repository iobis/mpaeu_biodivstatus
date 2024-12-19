################## MPA Europe - Biodiversity status exploring ##################
############################ WP3 - OBIS contribution ###########################
# November of 2024
# Author: Silas Principe
# Contact: s.principe@unesco.org
#
# Create analysis grid

library(sf)
library(dplyr)
library(h3jsr)
library(progress)
sf_use_s2(FALSE)

study_area <- read_sf("data/mpa_europe_starea_v2.shp")
study_area <- st_buffer(study_area, 0.2)
plot(study_area)

grid_area <- h3jsr::polygon_to_cells(study_area, res = 7)
grid_area_7 <- grid_area[[1]]

grid_area_6 <- h3jsr::get_parent(grid_area_7, res = 6)
grid_area_5 <- h3jsr::get_parent(grid_area_7, res = 5)

grid_area_6 <- unique(grid_area_6)
grid_area_5 <- unique(grid_area_5)

saveRDS(
    list(
        res_7 = grid_area_7,
        res_6 = grid_area_6,
        res_5 = grid_area_5
    ),
    file = "proc-data/grids_h3.rds"
)

to_polygon <- function(cells) {
    batches <- split(cells, ceiling(seq_along(cells) / 1000))

    result <- lapply(seq_along(batches), function(x) NULL)

    p <- progress_bar$new(total = length(batches))

    for (k in seq_along(batches)) {
        p$tick()

        result[[k]] <- sf::st_sf(h3jsr::cell_to_polygon(batches[[k]]))
        st_geometry(result[[k]]) <- "geometry"
    }

    result <- bind_rows(result)
}

grid_area_7_p <- to_polygon(grid_area_7)
grid_area_6_p <- to_polygon(grid_area_6)
grid_area_5_p <- to_polygon(grid_area_5)

sf::st_write(grid_area_7_p, dsn = "proc-data/grid_7.gpkg")
sf::st_write(grid_area_6_p, dsn = "proc-data/grid_6.gpkg")
sf::st_write(grid_area_5_p, dsn = "proc-data/grid_5.gpkg")

