################## MPA Europe - Biodiversity status exploring ##################
#################################### WP5 #######################################
# November of 2024
# Author: Silas Principe, Anna Adamo
# Contact: s.principe@unesco.org
#
# Uses Boosted Regression Trees to explore the relationship between
# richness (expressed as number of species) from modeled data to 
# environmental variables

library(dismo)
library(ggplot2)
library(plotly)
library(RColorBrewer)
library(akima)
library(dplyr)
library(gridExtra)
sf::sf_use_s2(FALSE)
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
# We will tune a boosted regression tree
# As we are predicting richness (count) we will use a poisson family
# To avoid overfitting and overly complex relationships we will define a priori
# some monotonic responses (either positive or negative), as follows:

#0 is random, 1 is positive and -1 is negative
# must be the same order as the variables
mon_responses <- c( 
    lon = 0, lat = 0,
    rugosity = 0, coastdist = -1, bathymetry = -1,
    thetao = 0, so = 0, siconc = -1, par = 1
)

mon_responses_alt <- c(mon_responses, wavefetch = -1)

m1_brt <- dismo::gbm.step(
    data = as.data.frame(stratified_sample),
    gbm.x = variables,
    gbm.y = "richness", tree.complexity = 3,
    var.monotone = mon_responses,
    family = "poisson", n.folds = 5,
    max.trees = 15000, learning.rate = 0.01, n.trees = 100
)
summary(m1_brt)

m1_brt$contributions %>%
    arrange(rel.inf) %>%
    mutate(var = factor(var, levels = var)) %>%
    ggplot() +
    geom_bar(aes(y = var, x = rel.inf), stat = "identity", fill = "#156ea5") +
    xlab("Contribution") + ylab("Variable") +
    scale_x_continuous(limits = c(0, 100), expand = c(0, 1)) +
    theme_light() + 
    theme(panel.grid = element_blank(), panel.grid.major.x = element_line(color = "gray80"),
    panel.border = element_blank())


# Get univariate partial response curves
resp_curves <- get_resp_curves(stratified_sample, m1_brt, variables)

# Save this plot
resp_curves$data$variable <- as.factor(resp_curves$data$variable)
levels(resp_curves$data$variable) <- c(
    "Depth", "Distance to Coast", "Latitude", "Longitude", "PAR",
    "Rugosity", "Sea ice", "Salinity", "Sea surface temperature"
)
resp_curves +
    theme(strip.background = element_rect(fill = "#ffffff", color = "black"),
          strip.text = element_text(color = "black", face = "bold"))
ggsave("figures/brt-respcurves.jpg", quality = 100)

# Get bivariate for the 4 most important variables
best_vars <- m1_brt$contributions[order(m1_brt$contributions$rel.inf, decreasing = T),][1:4,"var"]
best_comb <- combn(best_vars, 2, simplify = TRUE)
plot_list <- list()
for (i in seq_len(ncol(best_comb))) {
    plot_list[[i]] <- resp_3d(stratified_sample, m1_brt,
        best_comb[1, i], best_comb[2, i], variables,
        dynamic = FALSE
    )
}
grid.arrange(grobs = plot_list, ncol = 3, nrow = 2)

# Can also get dynamic
resp_3d(stratified_sample, m1_brt,
    "coastdist", "thetao", variables,
    dynamic = TRUE
)


# Do the same, now with wavefetch
m2_brt <- dismo::gbm.step(
    data = as.data.frame(stratified_sample),
    gbm.x = variables_alt,
    gbm.y = "richness", tree.complexity = 3,
    var.monotone = mon_responses_alt,
    family = "poisson", n.folds = 5,
    max.trees = 15000, learning.rate = 0.01, n.trees = 100
)
summary(m2_brt)

m2_brt$contributions %>%
    arrange(rel.inf) %>%
    mutate(var = factor(var, levels = var)) %>%
    ggplot() +
    geom_bar(aes(y = var, x = rel.inf), stat = "identity", fill = "#156ea5") +
    xlab("Contribution") + ylab("Variable") +
    scale_x_continuous(limits = c(0, 100), expand = c(0, 1)) +
    theme_light() + 
    theme(panel.grid = element_blank(), panel.grid.major.x = element_line(color = "gray80"),
    panel.border = element_blank())


# Get univariate partial response curves
get_resp_curves(stratified_sample, m2_brt, variables_alt)

# Get bivariate for the 4 most important variables
best_vars <- m2_brt$contributions[order(m2_brt$contributions$rel.inf, decreasing = T),][1:4,"var"]
best_comb <- combn(best_vars, 2, simplify = TRUE)
plot_list <- list()
for (i in seq_len(ncol(best_comb))) {
    plot_list[[i]] <- resp_3d(stratified_sample, m2_brt,
        best_comb[1, i], best_comb[2, i], variables_alt,
        dynamic = FALSE
    )
}
grid.arrange(grobs = plot_list, ncol = 3, nrow = 2)

# Can also get dynamic
resp_3d(stratified_sample, m2_brt,
    "coastdist", "thetao", variables_alt,
    dynamic = TRUE
)


# See how models perform with the validation dataset
pred_m1 <- predict(m1_brt, evaluate_sample, type = "response")
pred_m2 <- predict(m2_brt, evaluate_sample, type = "response")

sqrt(mean((evaluate_sample$richness - pred_m1)^2))
sqrt(mean((evaluate_sample$richness - pred_m2)^2))


# Because is a spatial problem, we can also do spatial cross validation to
# fit the model and latter get the cross-validation metrics
library(blockCV)

to_get_block <- terra::vect(
    stratified_sample[,c("richness", variables_alt)],
    crs = "EPSG:4326"
)
blocks_size <- blockCV::cv_spatial_autocor(x = to_get_block, column = variables[!variables %in% c("lon", "lat")])

blocks <- blockCV::cv_spatial(to_get_block, size = blocks_size$range)

# Will do only without wavefetch
m3_brt <- dismo::gbm.step(
    data = as.data.frame(stratified_sample),
    gbm.x = variables,
    gbm.y = "richness", tree.complexity = 3,
    var.monotone = mon_responses,
    family = "poisson", n.folds = 5,
    max.trees = 15000, learning.rate = 0.01, n.trees = 100,
    fold.vector = blocks$folds_ids, keep.fold.fit = T
)
summary(m3_brt)

rmse <- rep(NA, 5)
for (k in 1:5) {
    fit <- m3_brt$fold.fit[blocks$folds_ids == k]
    true <- stratified_sample$richness[blocks$folds_ids == k]
    rmse[k] <- sqrt(mean((true - exp(fit))^2))
}
mean(rmse)

m3_brt$contributions %>%
    arrange(rel.inf) %>%
    mutate(var = factor(var, levels = var)) %>%
    ggplot() +
    geom_bar(aes(y = var, x = rel.inf), stat = "identity", fill = "#156ea5") +
    xlab("Contribution") + ylab("Variable") +
    scale_x_continuous(limits = c(0, 100), expand = c(0, 1)) +
    theme_light() + 
    theme(panel.grid = element_blank(), panel.grid.major.x = element_line(color = "gray80"),
    panel.border = element_blank())


# Get univariate partial response curves
resp_curves <- get_resp_curves(stratified_sample, m3_brt, variables)

# Get bivariate for the 4 most important variables
best_vars <- m3_brt$contributions[order(m3_brt$contributions$rel.inf, decreasing = T),][1:4,"var"]
best_comb <- combn(best_vars, 2, simplify = TRUE)
plot_list <- list()
for (i in seq_len(ncol(best_comb))) {
    plot_list[[i]] <- resp_3d(stratified_sample, m3_brt,
        best_comb[1, i], best_comb[2, i], variables,
        dynamic = FALSE
    )
}
grid.arrange(grobs = plot_list, ncol = 3, nrow = 2)

# Can also get dynamic
resp_3d(stratified_sample, m3_brt,
    "coastdist", "thetao", variables,
    dynamic = TRUE
)

pred_m3 <- predict(m3_brt, evaluate_sample, type = "response")
sqrt(mean((evaluate_sample$richness - pred_m3)^2))

# Model 3 with spatial cross validation does much better than the others.

# Compare fits using a Bland-Altman plot
ba_plot_stats <- blandr::blandr.statistics(evaluate_sample$richness, pred_m3)
blandr::blandr.plot.ggplot(ba_plot_stats)
