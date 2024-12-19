################## MPA Europe - Biodiversity status exploring ##################
############################ WP3 - OBIS contribution ###########################
# November of 2024
# Author: Silas Principe
# Contact: s.principe@unesco.org
#
# Richness analysis based on raw data

library(arrow)
library(dplyr)
library(terra)
library(iNEXT)
library(h3jsr)
library(gsl)
library(data.table)

# Load processed data
undup <- "proc-data/undup"
ds <- open_dataset(undup)

# Add raw resolution cell
base <- terra::rast("../mpaeu_sdm/data/env/current/thetao_baseline_depthsurf_mean.tif")

ds_cell <- ds %>%
  map_batches(function(batch) {
    batch %>%
      as.data.frame() %>%
      mutate(raw_cell = as.integer(terra::cellFromXY(base, cbind(decimalLongitude, decimalLatitude)))) %>%
      as_record_batch()
  }) %>%
  collect()

# # Add h3_5
# ds_cell <- ds_cell %>%
#   map_batches(function(batch) {
#     batch %>%
#       as.data.frame() %>%
#       mutate(h3_5 = get_parent(h3_6, res = 5)) %>%
#       as_record_batch()
#   })

# ES50 - raw resolution (cell)
es50_raw_t1 <- ds_cell %>%
    group_by(raw_cell, scientificName) %>%
    summarize(ni = n()) %>%
    collect()

es50_raw_t2 <- es50_raw_t1 %>%
    ungroup() %>%
    group_by(raw_cell) %>%
    mutate(n = sum(ni))

esn <- 50
es50_raw_t3 <- es50_raw_t2 %>%
    group_by(raw_cell, scientificName) %>%
    mutate(
      hi = -(ni/n*log(ni/n)),
      si = (ni/n)^2,
      qi = ni/n,
      esi = case_when(
        n-ni >= esn ~ 1-exp(lngamma(n-ni+1)+lngamma(n-esn+1)-lngamma(n-ni-esn+1)-lngamma(n+1)),
        n >= esn ~ 1
      )
    )

es50_raw_t4 <- es50_raw_t3 %>%
    group_by(raw_cell) %>%
    summarize(
      n = sum(ni),
      sp = n(),
      shannon = sum(hi),
      simpson = sum(si),
      maxp = max(qi),
      es = sum(esi)
    )

raw_metrics_final <- es50_raw_t4 %>%
    mutate(
      hill_1 = exp(shannon),
      hill_2 = 1/simpson,
      hill_inf = 1/maxp
    )

# Create rasters
es_50 <- rast(base)
hill_1 <- rast(base)

europe <- vect("data/mpa_europe_starea_v2.shp")

es_50 <- mask(crop(es_50, europe), europe)
hill_1 <- mask(crop(hill_1, europe), europe)

plot(es_50)
plot(hill_1)

raw_data <- read.csv("proc-data/full_data_pts.csv")
colnames(raw_metrics_final)[1] <- "cell"
raw_data <- left_join(raw_data, raw_metrics_final)

write.csv(raw_data, "proc-data/raw_full_data_pts_metrics.csv", row.names = F)




# By grid 6 - ES50
es50_h6_t1 <- ds_cell %>%
    group_by(h3_6, scientificName) %>%
    summarize(ni = n()) %>%
    collect()

es50_h6_t2 <- es50_h6_t1 %>%
    ungroup() %>%
    group_by(h3_6) %>%
    mutate(n = sum(ni))

esn <- 50
es50_h6_t3 <- es50_h6_t2 %>%
    group_by(h3_6, scientificName) %>%
    mutate(
      hi = -(ni/n*log(ni/n)),
      si = (ni/n)^2,
      qi = ni/n,
      esi = case_when(
        n-ni >= esn ~ 1-exp(lngamma(n-ni+1)+lngamma(n-esn+1)-lngamma(n-ni-esn+1)-lngamma(n+1)),
        n >= esn ~ 1
      )
    )

es50_h6_t4 <- es50_h6_t3 %>%
    group_by(h3_6) %>%
    summarize(
      n = sum(ni),
      sp = n(),
      shannon = sum(hi),
      simpson = sum(si),
      maxp = max(qi),
      es = sum(esi)
    )

h6_metrics_final <- es50_h6_t4 %>%
    mutate(
      hill_1 = exp(shannon),
      hill_2 = 1/simpson,
      hill_inf = 1/maxp
    )

grid6_data <- read.csv("proc-data/grid6_data_pts.csv")

batches <- split(seq_len(nrow(grid6_data)), ceiling(seq_len(nrow(grid6_data))/1000))
grid6_data$h3_6 <- NA

for (i in seq_along(batches)) {
    grid6_data$h3_6[batches[[i]]] <- suppressMessages(h3jsr::point_to_cell(grid6_data[batches[[i]], c("x", "y")], res = 6))
}

grid6_data <- left_join(grid6_data, h6_metrics_final)

write.csv(grid6_data, "proc-data/raw_grid6_data_pts_metrics.csv", row.names = F)

un_h3 <- unique(es50_h6_t1$h3_6)
pb <- progress::progress_bar$new(total = length(un_h3))

es50_h6_t1_dt <- data.table(es50_h6_t1)

inext_output <- data.frame(
  h3_6 = un_h3,
  observed_q0 = NA,
  observed_q1 = NA,
  observed_q2 = NA,
  asymptote_rich = NA,
  asymptote_shan = NA,
  asymptote_simp = NA
)

for (i in seq_along(un_h3)) {
  pb$tick()

  in_data <- es50_h6_t1_dt[h3_6 == un_h3[i]]
  in_data <- in_data[!is.na(in_data$scientificName),]
  in_data <- data.frame(counts = in_data$ni, row.names = in_data$scientificName)

  in_result <- iNEXT(in_data, q = c(0,1,2), datatype = "abundance")

  obs <- in_result$iNextEst[[1]][in_result$iNextEst[[1]]$Method == "Observed",]

  inext_output$observed_q0[i] <- obs$qD[obs$Order.q == 0]
  inext_output$observed_q1[i] <- obs$qD[obs$Order.q == 1]
  inext_output$observed_q2[i] <- obs$qD[obs$Order.q == 2]

  inext_output$asymptote_rich[i] <- in_result$AsyEst[1]
  inext_output$asymptote_shan[i] <- in_result$AsyEst[2]
  inext_output$asymptote_shan[i] <- in_result$AsyEst[3]
}
