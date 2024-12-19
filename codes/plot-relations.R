################## MPA Europe - Biodiversity status exploring ##################
#################################### WP5 #######################################
# November of 2024
# Author: Silas Principe, Anna Adamo
# Contact: s.principe@unesco.org
#
# Plot relationship between richness and environmental variables

# Set the folder where files are
folder <- "proc-data"

full_data <- read.csv(file.path(folder, "full_data_pts.csv"))
grid6_data <- read.csv(file.path(folder, "grid6_data_pts.csv"))
grid5_data <- read.csv(file.path(folder, "grid5_data_pts.csv"))

grid6_data <- grid6_data[!is.na(grid6_data$richness),]
grid5_data <- grid5_data[!is.na(grid5_data$richness),]

plot_relationships <- function(dataset, variables = "all", n_highlight = 500) {

    require(ggplot2)
    require(dplyr)
    require(tidyr)

    vars <- colnames(dataset)
    vars <- vars[!vars %in% c("cell", "x", "X", "y", "Y", "richness", "ID")]

    if (variables != "all") {
        vars <- vars[vars %in% variables]
    }

    ed_ds <- dataset %>%
        select(all_of(c("richness", vars))) %>%
        pivot_longer(2:ncol(.), names_to = "variable", values_to = "value") 

    sample_high <- function(x, n) {
        total <- length(x)
        ret_vec <- rep(NA, total)
        if (n > total) {
            n <- total
        }
        to_do <- sample(1:total, n)
        ret_vec[to_do] <- "color"
        if (any(is.na(ret_vec))) {
            ret_vec[is.na(ret_vec)] <- "no_color"
        }
        return(ret_vec)
    }

    ed_ds <- ed_ds %>%
        group_by(variable) %>%
        mutate(highlight = sample_high(value, n = n_highlight))

    if (sum(ed_ds$highlight == "no_color") > 2000) {
        to_get <- sample(1:nrow(ed_ds[ed_ds$highlight == "no_color",]), 2000)
    }

    ggplot(ed_ds) +
        geom_point(data = ed_ds[ed_ds$highlight == "no_color",][to_get,], aes(x = value, y = richness), color = "grey70", alpha = 0.2) +
        geom_point(data = ed_ds[ed_ds$highlight == "color",], aes(x = value, y = richness), color = "#0661ab", alpha = 0.7) +
        facet_wrap(~variable, scales = "free_x") +
        theme_light()


}

plot_relationships(full_data)
plot_relationships(grid6_data)
plot_relationships(grid5_data)
