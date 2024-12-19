################## MPA Europe - Biodiversity status exploring ##################
############################ WP3 - OBIS contribution ###########################
# November of 2024
# Author: Silas Principe
# Contact: s.principe@unesco.org
#
# Load and pre-process OBIS data

library(DBI)
library(duckdb)
library(glue)
library(sf)
library(h3jsr)
sf_use_s2(FALSE)

study_area <- read_sf("data/mpa_europe_starea_v2.shp")
study_area <- st_buffer(study_area, 0.2)
plot(study_area)

path_obisfe <- "../mpaeu_sdm/data/raw/obis_20240625.parquet"
wkt <- st_as_text(st_geometry(study_area))
wkt <- st_as_text(st_as_sfc(st_bbox(study_area)))

sel_columns <- c(
    "id", "dataset_id",  "decimalLongitude", "decimalLatitude", "date_start", "date_mid", "date_end", "date_year",
    "scientificName", "originalScientificName", "minimumDepthInMeters", "maximumDepthInMeters", "depth",
    "coordinateUncertaintyInMeters", "absence", "taxonRank", "AphiaID", "redlist_category",
    "kingdom", "phylum", "class", '"order"', "family", "genus", "species", "institutionID",
    "collectionID", "datasetName", "basisOfRecord", "materialSampleID", "occurrenceID",
    "individualCount", "organismQuantity", "organismQuantityType",
    "occurrenceStatus", "eventID", "parentEventID", "samplingProtocol",
    "sampleSizeValue", "sampleSizeUnit", "samplingEffort", "eventDate",
    "locationID", "higherGeographyID", "higherGeography", "taxonID"
)
sel_columns <- paste0(sel_columns, collapse = ", ")

con <- dbConnect(duckdb(), dbdir = "my-db.duckdb", read_only = FALSE)
dbSendQuery(con, "install httpfs; load httpfs;")
dbSendQuery(con, "install spatial; load spatial;")

dbSendQuery(con, glue("
  create table geom_table as select
    {sel_columns}, ST_Point(decimalLatitude, decimalLongitude) AS geometry
  from read_parquet('{path_obisfe}')
"))

dbGetQuery(con, "select * from geom_table limit 1;")

species <- dbGetQuery(con, glue("
  select *
  from geom_table
  where ST_Intersects(geometry, ST_GeomFromText('{wkt}'))
"))

head(species)
nrow(species)

dbDisconnect(con)
fs::file_delete("my-db.duckdb")

species$h3_5 <- species$h3_6 <- species$h3_7 <- NA

batches <- split(seq_len(nrow(species)), ceiling(seq_len(nrow(species)) / 1000))

i = 0
pb <- progress::progress_bar$new(total = length(batches))
while (i < length(batches)) {
    pb$tick()
    i <- i + 1
    h3c <- suppressMessages(
        h3jsr::point_to_cell(
            species[batches[[i]],c("decimalLongitude", "decimalLatitude")], res = 7
        )
    )
    species[batches[[i]], "h3_7"] <- h3c
    species[batches[[i]], "h3_6"] <- get_parent(h3c, res = 6)
    species[batches[[i]], "h3_5"] <- get_parent(h3c, res = 5)
}

species$geometry <- NULL

species <- st_as_sf(species, coords = c("decimalLongitude", "decimalLatitude"),
                    crs = 4326, remove = FALSE)

head(species)

sfarrow::st_write_parquet(species, "proc-data/obis.parquet")
