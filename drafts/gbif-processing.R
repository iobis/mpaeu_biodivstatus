################## MPA Europe - Biodiversity status exploring ##################
############################ WP3 - OBIS contribution ###########################
# November of 2024
# Author: Silas Principe
# Contact: s.principe@unesco.org
#
# Load and pre-process GBIF data

# GBIF data is organized by order, so we take advantage of that to be able to
# process all the data in batches.

library(DBI)
library(duckdb)
library(glue)
library(sf)
library(h3jsr)
sf_use_s2(FALSE)

study_area <- read_sf("data/mpa_europe_starea_v2.shp")
study_area <- st_buffer(study_area, 0.2)
plot(study_area)

path_gbif <- "../mpaeu_sdm/data/raw/gbif_20240726/"
wkt <- st_as_text(st_geometry(study_area))
wkt <- st_as_text(st_as_sfc(st_bbox(study_area)))

fs::dir_create("proc-data/gbif/")

orders <- list.files(path_gbif)

counter <- 0
for (od in orders) {
    counter <- counter + 1
    cat("Processing order", counter, "out of", length(orders), "\n")

    con <- dbConnect(duckdb(), , dbdir = "my-db.duckdb", read_only = FALSE)
    dbSendQuery(con, "install httpfs; load httpfs;")
    dbSendQuery(con, "install spatial; load spatial;")

    dbSendQuery(con, glue("
  create table geom_table as select
    *, ST_Point(decimallatitude, decimallongitude) AS geometry
  from read_parquet('{path_gbif}/{od}/*')
"))

    dbGetQuery(con, "select * from geom_table limit 1;")

    species <- dbGetQuery(con, glue("
  select *
  from geom_table
  where ST_Intersects(geometry, ST_GeomFromText('{wkt}'))
"))

    dbDisconnect(con)
    fs::file_delete("my-db.duckdb")

    if (nrow(species) < 1) next

    species$h3_5 <- species$h3_6 <- species$h3_7 <- NA

    batches <- split(seq_len(nrow(species)), ceiling(seq_len(nrow(species)) / 1000))

    i <- 0
    pb <- progress::progress_bar$new(total = length(batches))
    while (i < length(batches)) {
        pb$tick()
        i <- i + 1
        h3c <- suppressMessages(
            h3jsr::point_to_cell(
                species[batches[[i]], c("decimallongitude", "decimallatitude")],
                res = 7
            )
        )
        species[batches[[i]], "h3_7"] <- h3c
        species[batches[[i]], "h3_6"] <- get_parent(h3c, res = 6)
        species[batches[[i]], "h3_5"] <- get_parent(h3c, res = 5)
    }

    species$geometry <- NULL

    species <- st_as_sf(species,
        coords = c("decimallongitude", "decimallatitude"),
        crs = 4326, remove = FALSE
    )

    suppressWarnings(sfarrow::st_write_parquet(species, glue("proc-data/gbif/{od}.parquet")))
}

fs::file_delete(list.files(pattern = "my-db"))
