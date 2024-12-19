################## MPA Europe - Biodiversity status exploring ##################
#################################### WP5 #######################################
# November of 2024
# Author: Silas Principe, Anna Adamo
# Contact: s.principe@unesco.org
#
# Process groups to add region and depth information

library(terra)
library(dplyr)

to_process <- list.files("proc-data-grp", full.names = T)

europe <- vect("data/mpa_europe_starea_v2.shp")
# From marineregions.org, pre-processed on QGIS to remove areas out of our region
iho <- vect("data/World_Seas_IHO_v3.shp")
# Correct north part
top_ext <- ext(c(ext(europe)[1:2]), 72, 90)
top_ext <- vect(top_ext)
top_ext$mpa_region <- "ne_atlantic"
iho <- rbind(iho, top_ext)
iho <- rasterize(iho, rast(res = 0.05, crs = "EPSG:4326"), field = "mpa_region")

iho <- mask(crop(iho, europe), europe)

# Correct those that were not assigned to any region based on nearest point
get_nearest_region <- function(x, y) {
    cell <- cellFromXY(iho, cbind(x, y))
    look <- matrix(c(rep(1, 12), 0, rep(1, 12)), nrow = 5, ncol = 5, byrow = T)
    adj <- adjacent(iho, cell, directions = look)
    vals <- iho[as.vector(adj)]
    vals <- as.character(unique(vals[,1]))
    if (all(is.na(vals))) {
        return(NA)
    } else if (length(na.omit(vals)) == 1) {
        return(na.omit(vals))
    } else {
        freq <- table(iho[as.vector(adj)][,1])
        return(names(which.max(freq)))
    }
}

for (k in seq_along(to_process)) {

    message("Processing ", k)

    dat <- read.csv(to_process[k])

    regions <- terra::extract(iho, dat[,c("lon", "lat")])

    dat$mpa_region <- regions[,2]

    if (any(is.na(dat$mpa_region))) {
        dat$mpa_region[is.na(dat$mpa_region)] <- unlist(lapply(
            which(is.na(dat$mpa_region)), function(k) {
                get_nearest_region(dat$lon[k], dat$lat[k])
            }
        ))
    }

    dat$depth <- ifelse(dat$btm < -200, "deep", "shallow")

    write.csv(dat, to_process[k], row.names = F)
    rm(dat)
}
