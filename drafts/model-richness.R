################## MPA Europe - Biodiversity status exploring ##################
############################ WP3 - OBIS contribution ###########################
# November of 2024
# Author: Silas Principe
# Contact: s.principe@unesco.org
#
# Richness analysis based on modeled data

library(terra)
library(mgcv)
library(xgboost)
library(ggplot2)

europe <- vect("data/mpa_europe_starea_v2.shp")
iho <- vect("data/World_Seas_IHO_v3.shp")
# Correct north part
top_ext <- ext(c(ext(europe)[1:2]), 72, 90)
top_ext <- vect(top_ext)
top_ext$mpa_region <- "ne_atlantic"
iho <- rbind(iho, top_ext)
iho <- rasterize(iho, rast(res = 0.05, crs = "EPSG:4326"), field = "mpa_region")

env_layers <- list.files("../mpaeu_sdm/data/env/current", full.names = T)
env_layers <- env_layers[grepl("depthsurf", env_layers)]
env_layers <- env_layers[grepl("_mean", env_layers)]
env_layers <- env_layers[grepl("chl|no3|o2|par|siconc|so|tas|thetao|sws", env_layers)]

terrain_layers <- list.files("../mpaeu_sdm/data/env/terrain", full.names = T)
terrain_layers <- terrain_layers[grepl("bathymetry_mean|distcoast|rugosity|wavefetch", terrain_layers)]
terrain_layers <- terrain_layers[!grepl("aux", terrain_layers)]

env_layers <- c(rast(c(env_layers, terrain_layers)), iho)

env_layers <- mask(crop(env_layers, europe), europe)

# Load grid
grid_6 <- vect("proc-data/grid_6.gpkg")
grid_5 <- vect("proc-data/grid_5.gpkg")

# Load richness
richness <- rast("data/richness/richness_full_rf_202411.tif")

# Start with non-grided version
richness_masked <- mask(richness, europe)
richness_masked <- crop(richness_masked, europe)
names(richness_masked) <- "richness"

# Extract by grid
richness_6 <- terra::extract(richness_masked, grid_6, fun = max, na.rm = T)
richness_5 <- terra::extract(richness_masked, grid_5, fun = max, na.rm = T)

env_6 <- terra::extract(subset(env_layers, "mpa_region", negate = T), grid_6, fun = mean, na.rm = T)
env_5 <- terra::extract(subset(env_layers, "mpa_region", negate = T), grid_5, fun = mean, na.rm = T)

region_6 <- terra::extract(subset(env_layers, "mpa_region"), grid_6, fun = table, na.rm = T)
nams_region <- names(region_6[,-1])
region_6 <- apply(region_6[,-1], 1, which.max)
region_6 <- nams_region[region_6]

region_5 <- terra::extract(subset(env_layers, "mpa_region"), grid_5, fun = table, na.rm = T)
nams_region <- names(region_5[,-1])
region_5 <- apply(region_5[,-1], 1, which.max)
region_5 <- nams_region[region_5]

grid_6 <- cbind(cbind(grid_6, data.frame(richness = richness_6[,-1])), cbind(env_6[,-1], mpa_region = region_6))
grid_5 <- cbind(cbind(grid_5, data.frame(richness = richness_5[,-1])), cbind(env_5[,-1], mpa_region = region_5))

centroids_6 <- terra::centroids(grid_6)
centroids_5 <- terra::centroids(grid_5)

points_6 <- as.data.frame(centroids_6, geom = "XY")
points_5 <- as.data.frame(centroids_5, geom = "XY")

points_6 <- points_6[!is.na(points_6$richness),]
points_5 <- points_5[!is.na(points_5$richness),]

valid_data <- as.data.frame(c(richness_masked, env_layers), cell = T, xy = T)
valid_data <- valid_data[!is.na(valid_data$richness),]

varnames <- colnames(valid_data)[5:(ncol(valid_data)-1)]

valid_data[,varnames] <- apply(valid_data[,varnames], 2, round, digits = 2)
points_6[,varnames] <- apply(points_6[,varnames], 2, round, digits = 2)
points_5[,varnames] <- apply(points_5[,varnames], 2, round, digits = 2)

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

valid_data$mpa_region[is.na(valid_data$mpa_region)] <- unlist(lapply(
    which(is.na(valid_data$mpa_region)), function(k) {
        get_nearest_region(valid_data$x[k], valid_data$y[k])
    }
))

summary(valid_data)
# There are still 7 NA points for MPA region. Let's get the nearest valid value
not_valid <- valid_data[is.na(valid_data$mpa_region), c("x", "y")]
near_pts <- nearest(
    vect(not_valid, crs = "EPSG:4326", geom = c("x", "y")),
    vect(valid_data[!is.na(valid_data$mpa_region), c("x", "y")],
        crs = "EPSG:4326", geom = c("x", "y")
    )
)
valid_data$mpa_region[is.na(valid_data$mpa_region)] <- valid_data$mpa_region[!is.na(valid_data$mpa_region)][near_pts$to_id]
table(valid_data$mpa_region)

points_6$mpa_region[is.na(points_6$mpa_region)] <- unlist(lapply(
    which(is.na(points_6$mpa_region)), function(k) {
        get_nearest_region(points_6$x[k], points_6$y[k])
    }
))
summary(points_6)
table(points_6$mpa_region)

points_5$mpa_region[is.na(points_5$mpa_region)] <- unlist(lapply(
    which(is.na(points_5$mpa_region)), function(k) {
        get_nearest_region(points_5$x[k], points_5$y[k])
    }
))
summary(points_5)
table(points_5$mpa_region)


# Add info of shallow
# > 200m deep | <= 200 is shallow
valid_data$depth <- ifelse(valid_data$bathymetry_mean < -200, "deep", "shallow")
points_6$depth <- ifelse(points_6$bathymetry_mean < -200, "deep", "shallow")
points_5$depth <- ifelse(points_5$bathymetry_mean < -200, "deep", "shallow")

write.csv(valid_data, "proc-data/full_data_pts.csv", row.names = F)
write.csv(points_6, "proc-data/grid6_data_pts.csv", row.names = F)
write.csv(points_5, "proc-data/grid5_data_pts.csv", row.names = F)

# Verify correlation
usdm::vifstep(env_layers, th = 5)

env_layers <- subset(env_layers, "tas_mean", negate = T)
usdm::vifstep(env_layers, th = 5)

env_layers <- subset(env_layers, c("no3_mean"), negate = T)
usdm::vifstep(env_layers, th = 5)

env_layers <- subset(env_layers, c("o2_mean"), negate = T)
usdm::vifstep(env_layers, th = 5)
    
library(mapgl)
maplibre(bounds = sf::st_as_sf(grid_5)) |> 
  add_fill_layer(id = "nc_data",
                 source = sf::st_as_sf(grid_5),
                 fill_color = interpolate(
                column = "richness",
                values = c(0, 3286),
                stops = c("lightblue", "darkblue"),
                na_color = "lightgrey"
                ),
                 fill_opacity = 1)

valid_data <- as.data.frame(c(richness_masked, env_layers), cell = T, xy = T)
valid_data <- valid_data[!is.na(valid_data$richness),]








##### MODEL PART
library(dplyr)
library(mgcv)
valid_data <- read.csv("proc-data/full_data_pts.csv")

stratified_sample <- valid_data %>%
  group_by(richness) %>%
  slice_sample(n = 10, replace = FALSE) %>%
  ungroup()

# Try gam model
#variables <- colnames(stratified_sample)[grep("chl_mean", colnames(stratified_sample)):ncol(stratified_sample)]
variables <- c("bathymetry_mean", "thetao_mean", "siconc_mean", "PAR_mean_mean", "o2_mean", "y")

form <- paste(paste0("s(", variables,")"), collapse = "+")
form <- as.formula(paste("richness ~", form))

m1 <- gam(form, family = poisson(), data = stratified_sample) #Tweedie(p=2)
summary(m1)
plot(m1)
gam.check(m1)

# m2
form <- as.formula("richness ~ s(thetao_mean) + s(siconc_mean) + s(PAR_mean_mean) + bathymetry_mean + y")
m2 <- gam(form, family = poisson(), data = stratified_sample) #Tweedie(p=2)
summary(m2)
plot(m2)
gam.check(m2)

# m3
m3 <- gam(richness ~ s(x, y, bs = 'gp', k = 100, m = 2) + s(thetao_mean) + so_mean + bathymetry_mean,
          family = poisson(), data = stratified_sample)
summary(m3)
plot(m3)
gam.check(m3)

AIC(m1);AIC(m2);AIC(m3)

# m4
m4 <- gamm(
  richness ~ s(thetao_mean) + so_mean + bathymetry_mean,
  data = stratified_sample, family = poisson(),
  correlation = corSpatial(form = ~ x + y, type = 'gaussian'),
  niterPQL = 5
)
plot(m4$gam)
summary(m4$gam)

#### BRT
m5 <- dismo::gbm.step(
    data = as.data.frame(stratified_sample),
    gbm.x = c("y", "x", "wavefetch", "rugosity", "coastdist", "bathymetry_mean",
    "thetao_mean", "so_mean", "siconc_mean", "PAR_mean_mean"),
    gbm.y = "richness",
    family = "poisson", max.trees = 15000, learning.rate = 0.01
)
summary(m5)


m6 <- dismo::gbm.step(
    data = as.data.frame(stratified_sample),
    gbm.x = c("y", "wavefetch", "rugosity", "coastdist", "bathymetry_mean",
    "thetao_mean", "so_mean", "siconc_mean", "PAR_mean_mean"),
    gbm.y = "richness", tree.complexity = 3,
    var.monotone = c(0, 0, 0, -1, -1, 0, 0, -1, 1),
    family = "poisson", max.trees = 15000, learning.rate = 0.01
)
summary(m6)




get_resp_curves(stratified_sample, m6, c("y", "wavefetch", "rugosity", "coastdist", "bathymetry_mean",
    "thetao_mean", "so_mean", "siconc_mean", "PAR_mean_mean"))


# Try BRT
library(dismo)

m2 <- gbm.step(data = as.data.frame(stratified_sample),
               gbm.x = variables, gbm.y = "richness", family = "poisson")


ranges <- apply(valid_data[,variables], 2, range)
means <- apply(valid_data[,variables], 2, mean)

resp_curves <- lapply(variables, function(x) NULL)

for (i in seq_along(variables)) {
    to_pred <- data.frame(matrix(ncol = length(variables), nrow = 100))
    to_pred[,] <- rep(means, each = 100)
    to_pred[,i] <- seq(from = ranges[1,variables[i]], to = ranges[2,variables[i]], length.out = 100)
    colnames(to_pred) <- variables

    pred <- predict(m2, to_pred, type = "response")

    resp_curves[[i]] <- data.frame(value = to_pred[,i], prediction = pred, variable = variables[i])
}

resp_curves <- do.call("rbind", resp_curves)

ggplot(resp_curves) +
    geom_line(aes(x = value, y = prediction)) +
    facet_wrap(~variable, scales = "free") +
    theme_light()






# Try xgboost model
xg_vars <- obissdm::sdm_options("xgboost")
xg_grid <- expand.grid(xg_vars)

dat <- stratified_sample[,variables]
p <- stratified_sample[["richness"]]
cvfold <- sample(1:5, size = nrow(dat), replace = T)

# Tune
tune_xgboost <- function(dat, p, grid, cvfold) {
    tune_results <- list()
    for (k in 1:nrow(grid)) {
        cat("Tunning option", k, "out of", nrow(grid), "\n")
        cv_res <- rep(NA, length(unique(cvfold)))
        for (cv in unique(cvfold)) {
            cat("Fold", cv, "\n")
            fit_dat <- dat[cvfold != cv,]
            fit_p <- p[cvfold != cv]
            pred_dat <- dat[cvfold != cv,]
            pred_p <- p[cvfold != cv]

            fm <- xgboost::xgboost(data = as.matrix(fit_dat), label = fit_p, 
                max_depth = grid$depth[k], gamma = grid$gamma[k], 
                #scale_pos_weight = grid$scale_pos_weight[k],
                #eval_metric = "poisson-nloglik",
                nrounds = grid$rounds[k], 
                learning_rate = grid$shrinkage[k], verbose = 0, nthread = 1, 
                objective = "reg:tweedie")
                #objective = "count:poisson")
            pred <- predict(fm, as.matrix(pred_dat), type = "response")
            pred <- exp(pred)

            cv_res[cv] <- sqrt(mean((pred_p - pred)^2))
        }
        tune_results[[k]] <- cv_res
    }
    return(tune_results)
}

xg_tune <- tune_xgboost(dat, p, xg_grid, cvfold)

m2 <- xgboost::xgboost(data = as.matrix(dat), label = p, 
        max_depth = best_tune$depth, gamma = best_tune$gamma, 
        scale_pos_weight = best_tune$scale_pos_weight, nrounds = best_tune$rounds, 
        learning_rate = best_tune$shrinkage,
        verbose = 0, nthread = 1, 
        objective = "binary:logistic")

variables <- colnames(stratified_sample)[grep("chl_mean", colnames(stratified_sample)):ncol(stratified_sample)]

form <- paste(paste0("s(", variables,")"), collapse = "+")
form <- as.formula(paste("richness ~", form))

m1 <- gam(form, family = poisson(), data = stratified_sample)

dtrain <- xgb.DMatrix(data = as.matrix(stratified_sample[,variables]), label = as.integer(stratified_sample$richness))

# Train the model
xgb_model <- xgb.train(
  data = dtrain,
  #objective = "count:poisson",
  objective = "reg:tweedie",
  #eval_metric = "poisson-nloglik",  # Evaluation metric for Poisson models
  max_depth = 3,
  eta = 0.1,
  nrounds = 50
)
y_pred_log <- predict(xgb_model, newdata = as.matrix(stratified_sample[,variables]))

# Convert predictions to response scale
y_pred_response <- exp(y_pred_log)

# View predictions
head(y_pred_response)



# Get partial response curves
ranges <- apply(valid_data[,variables], 2, range)
means <- apply(valid_data[,variables], 2, mean)

resp_curves <- lapply(variables, function(x) NULL)

for (i in seq_along(variables)) {
    to_pred <- data.frame(matrix(ncol = length(variables), nrow = 100))
    to_pred[,] <- rep(means, each = 100)
    to_pred[,i] <- seq(from = ranges[1,variables[i]], to = ranges[2,variables[i]], length.out = 100)
    colnames(to_pred) <- variables

    pred <- predict(m1, to_pred, type = "response")

    resp_curves[[i]] <- data.frame(value = to_pred[,i], prediction = pred, variable = variables[i])
}

resp_curves <- do.call("rbind", resp_curves)

ggplot(resp_curves) +
    geom_line(aes(x = value, y = prediction)) +
    facet_wrap(~variable, scales = "free") +
    theme_light()

# Extract information
grid_richness <- terra::extract(richness, grid_6, fun = max, na.rm = T)

grid_richness <- grid_richness[,2]

grid_6$richness <- grid_richness

mapview::mapview(sf::st_as_sf(grid_6))




library(arrow)
library(dplyr)

r <- read_parquet("~/Research/obis/dataviz/edna-dashboard/data/species_sites_risk_full.parquet")

at_risk <- function(x) {ifelse(x > 0, 1, 0)}

proc <- r %>%
    mutate(across(5:13, ~at_risk(.x)))

proc <- proc %>%
    pivot_longer(5:13, names_to = "scenario", values_to = "risk")

proc2 <- proc %>%
    group_by(higherGeography, scenario) %>%
    distinct(scientificName, .keep_all = T) %>%
    summarise(in_risk = sum(risk, na.rm = T), total = n())

proc2 <- proc2 %>%
    filter(!grepl("dec50", scenario)) %>%
    mutate(perc_risk = (in_risk * 100) / total)

proc2$scenario <- gsub("baseline", "current", proc2$scenario)

write.csv(proc2, "risk_species_perc_bysite.csv", row.names = F)

proc3 <- proc %>%
    group_by(scenario) %>%
    distinct(scientificName, .keep_all = T) %>%
    summarise(in_risk = sum(risk, na.rm = T), total = n())

proc3 <- proc3 %>%
    filter(!grepl("dec50", scenario)) %>%
    mutate(perc_risk = (in_risk * 100) / total)

proc3$scenario <- gsub("baseline", "current", proc3$scenario)

write.csv(proc3, "risk_species_perc_total.csv", row.names = F)
