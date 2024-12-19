################## MPA Europe - Biodiversity status exploring ##################
############################ WP3 - OBIS contribution ###########################
# November of 2024
# Author: Silas Principe
# Contact: s.principe@unesco.org
#
# Remove duplicates across OBIS/GBIF
# It uses the same base as the obissdm::outqc_dup_check
# See https://github.com/iobis/mpaeu_msdm

library(arrow)
library(dplyr)
library(DBI)
library(duckdb)
library(glue)
library(future)
library(furrr)
set.seed(2023)
outfolder <- "proc-data/undup"
fs::dir_create(outfolder)

species_list <- read.csv("data/all_splist_20240724.csv")
proc_path <- "../mpaeu_sdm/data/species"
proc_files <- list.files(proc_path)
proc_files <- as.numeric(gsub("key=", "", gsub("\\.parquet", "", proc_files)))
species_list <- species_list[species_list$AphiaID %in% proc_files,]

path_obis <- "../mpaeu_sdm/data/raw/obis_20240625.parquet"
path_gbif <- "../mpaeu_sdm/data/raw/gbif_20240726/"

undup_records <- function(i, species_list, path_obis, path_gbif, outfolder, proc_path, verbose = FALSE) {

    if (file.exists(file.path(outfolder, paste0("key=", species_list$taxonID[i], ".parquet")))) {
        return(species_list$taxonID[i])
    }

    base <- terra::rast("../mpaeu_sdm/data/env/current/thetao_baseline_depthsurf_mean.tif")

    sp_obis_id <- as.integer(species_list$taxonID[i])
    sp_gbif_id <- as.integer(species_list$gbif_speciesKey[i])

    sp_obis_id <- ifelse(is.na(sp_obis_id), 0, sp_obis_id)
    sp_gbif_id <- ifelse(is.na(sp_gbif_id), 0, sp_gbif_id)

    con <- dbConnect(duckdb())

    obis_sp <- dbGetQuery(con, glue(
        "select AphiaID, scientificName, decimalLongitude, decimalLatitude, date_mid, date_year, dataset_id, datasetName, occurrenceID, eventID, basisOfRecord
        from read_parquet('{path_obis}')
        where AphiaID = {sp_obis_id}
        "
    ))

    gbif_sp <- dbGetQuery(con, glue(
        "select specieskey, species, decimallongitude, decimallatitude, eventdate, day, month, year, datasetkey, basisofrecord
        from read_parquet('{path_gbif}/*/*')
        where specieskey = {sp_gbif_id}
        "
    ))

    to_remove <- c("FossilSpecimen", "FOSSIL_SPECIMEN", "LivingSpecimen", 
                "LIVING_SPECIMEN")

    obis_sp <- obis_sp %>%
        filter(!basisOfRecord %in% to_remove) %>%
        filter(!is.na(date_year)) %>%
        filter(!is.na(decimalLongitude) & !is.na(decimalLongitude)) %>%
        filter(date_year >= 1950)

    obis_sp$month <- lubridate::month(as.POSIXct(obis_sp$date_mid / 1000, origin = "1970-01-01"))

    gbif_sp <- gbif_sp %>%
        filter(!basisofrecord %in% to_remove) %>%
        filter(!is.na(year)) %>%
        filter(!is.na(decimallongitude) & !is.na(decimallongitude)) %>%
        filter(year >= 1950) %>%
        rename(scientificName = species, decimalLongitude = decimallongitude,
            decimalLatitude = decimallatitude, date_year = year)

    to_h3 <- function(dataset, resolution = 8, verbose = FALSE) {
        batches <- split(seq_len(nrow(dataset)), ceiling(seq_len(nrow(dataset))/100))
        if (verbose) {
            pb <- progress::progress_bar$new(total = length(batches))
        }
        h3_vec <- rep(NA, nrow(dataset))
        for (k in seq_len(length(batches))) {
            if (verbose) pb$tick()
            h3_vec[batches[[k]]] <- suppressMessages(
                h3jsr::point_to_cell(dataset[batches[[k]] ,c("decimalLongitude", "decimalLatitude")], res = resolution)
            )
        }
        return(h3_vec)
    }
    
    obis_sp$h3_8 <- to_h3(obis_sp, verbose = verbose)
    gbif_sp$h3_8 <- to_h3(gbif_sp, verbose = verbose)

    if (nrow(gbif_sp) > 0 & nrow(obis_sp) > 0) {
        all_data <- bind_rows(obis_sp, gbif_sp)
    } else if (nrow(gbif_sp) > 0) {
        all_data <- gbif_sp
    } else {
        all_data <- obis_sp
    }

    undup_data <- all_data %>%
        mutate(cell = factor(paste(h3_8, date_year, month, sep = "_"))) %>%
        distinct(cell, .keep_all = T) %>%
        select(-cell)

    proc_file <- arrow::read_parquet(file.path(proc_path, paste0("key=", species_list$taxonID[i], ".parquet")))
    proc_file$h3_6 <- to_h3(proc_file, resolution = 6)

    undup_data$h3_6 <- to_h3(undup_data, resolution = 6)

    undup_data <- undup_data %>%
        filter(h3_6 %in% proc_file$h3_6)

    undup_data$thetao <- terra::extract(base, undup_data[,c("decimalLongitude", "decimalLatitude")], ID = F)[,1]

    undup_data <- undup_data %>%
        filter(!is.na(thetao)) %>%
        select(-thetao)

    write_parquet(undup_data, file.path(outfolder, paste0("key=", species_list$taxonID[i], ".parquet")))

    DBI::dbDisconnect(con)

    return(species_list$taxonID[i])
}

plan(multisession, workers = 6)
result <- future_map(seq_len(nrow(species_list)),
                     species_list, path_obis, path_gbif,
                     outfolder, proc_path,
                     .f = undup_records,
                     .progress = TRUE, .options = furrr_options(seed = T))