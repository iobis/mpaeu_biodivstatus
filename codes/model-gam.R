################## MPA Europe - Biodiversity status exploring ##################
#################################### WP5 #######################################
# November of 2024
# Author: Silas Principe, Anna Adamo
# Contact: s.principe@unesco.org
#
# Uses GAMs to explore the relationship between
# richness (expressed as number of species) from modeled data to 
# environmental variables

library(mgcv)
library(ggplot2)
library(plotly)
library(RColorBrewer)
library(akima)
library(dplyr)
library(gridExtra)
library(terra)
source("functions/response-curves.R")
set.seed(2023)

# Load full data
valid_data <- read.csv("proc-data/full_data_pts.csv")
valid_data$ID <- seq_len(nrow(valid_data))

# Fill wavefetch NAs with 0
valid_data$wavefetch[is.na(valid_data$wavefetch)] <- 0

colnames(valid_data)[grepl("^x|^y", colnames(valid_data))] <- c("lon", "lat")
colnames(valid_data)[seq_len(ncol(valid_data) - 1)] <- gsub(
    "_mean", "",
    tolower(colnames(valid_data))[seq_len(ncol(valid_data) - 1)]
)

# Get a stratified sample to fit the model
stratified_sample <- valid_data %>%
  group_by(richness) %>%
  slice_sample(n = 5, replace = FALSE) %>%
  ungroup()

evaluate_sample <- valid_data %>%
    filter(!ID %in% stratified_sample$ID) %>%
    group_by(richness) %>%
    slice_sample(n = 5, replace = FALSE) %>%
    ungroup()

# List variables of interest
variables <- c(
    # Position
    "lon", "lat",
    # terrain
    "rugosity", "coastdist", "bathymetry",
    # environment
    "thetao", "so", "siconc", "par"
)

variables_alt <- c(variables, "wavefetch")


# Model

form <- as.formula("richness ~ lon + lat + rugosity + coastdist + bathymetry + s(thetao, bs = 'cr') + s(so, bs = 'cr') + s(siconc, bs = 'cr') + s(par, bs = 'cr')")
m1_gam <- gam(form, family = nb(), data = stratified_sample)
summary(m1_gam)

# Get univariate partial response curves
get_resp_curves(stratified_sample, m1_gam, variables)

# Get bivariate for the some variables
best_comb <- combn(c("thetao", "so", "lat", "coastdist"), 2, simplify = TRUE)
plot_list <- list()
for (i in seq_len(ncol(best_comb))) {
    plot_list[[i]] <- resp_3d(stratified_sample, m1_gam,
        best_comb[1, i], best_comb[2, i], variables,
        dynamic = FALSE
    )
}
grid.arrange(grobs = plot_list, ncol = 3, nrow = 2)

# Can also get dynamic
resp_3d(stratified_sample, m1_gam,
    "coastdist", "thetao", variables,
    dynamic = TRUE
)


# Do the same, now with wavefetch
form <- as.formula("richness ~ lon + lat + rugosity + coastdist + bathymetry + s(thetao, bs = 'cr') + s(so, bs = 'cr') + s(siconc, bs = 'cr') + s(par, bs = 'cr') + s(wavefetch, bs = 'cr')")
m2_gam <- gam(form, family = nb(), data = stratified_sample)
summary(m2_gam)

# Get univariate partial response curves
get_resp_curves(stratified_sample, m2_gam, variables_alt)

# Get bivariate for the some variables
best_comb <- combn(c("thetao", "so", "lat", "coastdist"), 2, simplify = TRUE)
plot_list <- list()
for (i in seq_len(ncol(best_comb))) {
    plot_list[[i]] <- resp_3d(stratified_sample, m2_gam,
        best_comb[1, i], best_comb[2, i], variables_alt,
        dynamic = FALSE
    )
}
grid.arrange(grobs = plot_list, ncol = 3, nrow = 2)

# Can also get dynamic
resp_3d(stratified_sample, m2_gam,
    "coastdist", "thetao", variables_alt,
    dynamic = TRUE
)


# See how models perform with the validation dataset
pred_m1 <- predict(m1_gam, evaluate_sample, type = "response")
pred_m2 <- predict(m2_gam, evaluate_sample, type = "response")

sqrt(mean((evaluate_sample$richness - pred_m1)^2))
sqrt(mean((evaluate_sample$richness - pred_m2)^2))

AIC(m1_gam)
AIC(m2_gam)
anova(m1_gam, m2_gam)


# We can also fit a model similar to kriging, incorporating the spatial part
# using a Gaussian process basis which allows latitude and longitude to interact
# We do without wavefetch
form <- as.formula("richness ~ s(lon, lat, bs = 'gp', k = 100, m = 2) + rugosity + coastdist + bathymetry + s(thetao, bs = 'cr') + s(so, bs = 'cr') + s(siconc, bs = 'cr') + s(par, bs = 'cr')")
m3_gam <- gam(form, family = nb(), data = stratified_sample)
summary(m3_gam)
AIC(m3_gam)

pred_m3 <- predict(m3_gam, evaluate_sample, type = "response")

sqrt(mean((evaluate_sample$richness - pred_m3)^2))

# Get univariate partial response curves
get_resp_curves(stratified_sample, m3_gam, variables)

# Get bivariate for the some variables
best_comb <- combn(c("thetao", "so", "lat", "coastdist"), 2, simplify = TRUE)
plot_list <- list()
for (i in seq_len(ncol(best_comb))) {
    plot_list[[i]] <- resp_3d(stratified_sample, m1_gam,
        best_comb[1, i], best_comb[2, i], variables,
        dynamic = FALSE
    )
}
grid.arrange(grobs = plot_list, ncol = 3, nrow = 2)

# Can also get dynamic
resp_3d(stratified_sample, m3_gam,
    "coastdist", "thetao", variables,
    dynamic = TRUE
)

# Predict to the dataset
fit_pred <- predict(m3_gam, stratified_sample, type = "response")
fit_sf <- sf::st_as_sf(cbind(pred_richness = fit_pred, stratified_sample[,c("lon", "lat")]),
                       crs = "EPSG:4326", coords = c("lon", "lat"))

ggplot(fit_sf) +
    geom_sf(aes(color = pred_richness)) +
    scale_color_distiller(palette = "YlGnBu") +
    coord_sf() +
    theme_light()

# Predict to a grid
new_grid <- sf::st_make_grid(fit_sf, cellsize = 0.5, square = F)
new_grid <- sf::st_as_sf(new_grid)
new_grid$sf_ID <- seq_len(nrow(new_grid))
points_sf <- sf::st_as_sf(stratified_sample, coords = c("lon", "lat"), crs = "EPSG:4326")
joined_sf <- sf::st_join(points_sf, sf::st_as_sf(new_grid))
joined_sf <- joined_sf %>%
  group_by(sf_ID) %>% 
  summarise(across(2:15, ~mean(.x, na.rm = T)))

joined_sf_coords <- sf::st_centroid(joined_sf) %>%
    cbind(., sf::st_coordinates(.))
colnames(joined_sf_coords)[16:17] <- c("lon", "lat")

joined_sf$pred_richness <- predict(m3_gam, joined_sf_coords, type = "response")

new_grid <- sf::st_as_sf(new_grid)
new_grid$pred_richness <- NA
new_grid$pred_richness[joined_sf$sf_ID] <- joined_sf$pred_richness

ggplot(new_grid) +
    geom_sf(aes(fill = pred_richness), color = NA) +
    scale_fill_distiller(palette = "YlGnBu", na.value = NA) +
    coord_sf() +
    theme_light()

# Compare fits using a Bland-Altman plot
ba_plot_stats <- blandr::blandr.statistics(evaluate_sample$richness, pred_m3)
blandr::blandr.plot.ggplot(ba_plot_stats)


# With poisson instead of bn
form <- as.formula("richness ~ s(lon, lat, bs = 'gp', k = 100, m = 2) + rugosity + coastdist + bathymetry + s(thetao, bs = 'cr') + s(so, bs = 'cr') + s(siconc, bs = 'cr') + s(par, bs = 'cr')")
stratified_sample$richness <- as.integer(stratified_sample$richness)
m4_gam <- gam(form, family = Tweedie(p=1.1), data = stratified_sample)
summary(m4_gam)
AIC(m4_gam)

pred_m4 <- predict(m4_gam, evaluate_sample, type = "response")

sqrt(mean((evaluate_sample$richness - pred_m4)^2))

# Get univariate partial response curves
get_resp_curves(stratified_sample, m4_gam, variables)

# Get bivariate for the some variables
best_comb <- combn(c("thetao", "so", "lat", "coastdist"), 2, simplify = TRUE)
plot_list <- list()
for (i in seq_len(ncol(best_comb))) {
    plot_list[[i]] <- resp_3d(stratified_sample, m1_gam,
        best_comb[1, i], best_comb[2, i], variables,
        dynamic = FALSE
    )
}
grid.arrange(grobs = plot_list, ncol = 3, nrow = 2)

# Can also get dynamic
resp_3d(stratified_sample, m4_gam,
    "so", "thetao", variables,
    dynamic = TRUE
)

