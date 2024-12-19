# Get and plot response curves

# 3D bivariate plot
resp_3d <- function(valid_data, model,
                    target_var1, target_var2, variables, dynamic = TRUE, zlim = NULL) {

    variables <- variables[!variables %in% c(target_var1, target_var2)]

    eval(parse(text = paste0(
        "xyz <- expand.grid(",
        target_var1, " = seq(min(valid_data$", target_var1, "), max(valid_data$", target_var1, "), length=50),",
        target_var2, " = seq(min(valid_data$", target_var2, "), max(valid_data$", target_var2, "), length=50),",
        paste(paste0(
            variables, "=", "mean(valid_data$", variables, ", na.rm = T)"
        ), collapse = ","),
        ")"
    )))

    xyz$z <- predict(model, newdata=xyz, type='response')

    if (!dynamic) {
        require(lattice)
        if (is.null(zlim)) {
           z <- xyz$z
           zlim = if (is.factor(z)) levels(z) else range(z, finite = TRUE)
        }
        cls <- colorRampPalette(rev(RColorBrewer::brewer.pal(11, 'RdYlBu')))(100)
        wireframe(as.formula(paste0("z ~ ", target_var1, " + ", target_var2)), data = xyz, 
            zlab = list("Richness", rot=90), zlim = zlim,
                drape = TRUE, col.regions = cls, scales = list(arrows = FALSE))
    } else {
        require(plotly)

        x_vals <- seq(min(xyz[[target_var1]], na.rm = T), max(xyz[[target_var1]], na.rm = T), length.out = 100)
        y_vals <- seq(min(xyz[[target_var2]], na.rm = T), max(xyz[[target_var2]], na.rm = T), length.out = 100)

        interpolation <- suppressWarnings(
            akima::interp(x = xyz[[target_var1]], y = xyz[[target_var2]], z = xyz$z,
                            xo = x_vals, yo = y_vals, duplicate = "strip")
        )

        z_vals <- interpolation$z
        cls <- colorRampPalette(rev(RColorBrewer::brewer.pal(11, 'RdYlBu')))(100)

        plot_ly(x = x_vals, y = y_vals, z = z_vals, type = "surface", colorscale = cls) %>%
        layout(scene = list(
            xaxis = list(title = target_var1),
            yaxis = list(title = target_var2),
            zaxis = list(title = 'Richness')
        ))
    }
}

# Get partial response curves
get_resp_curves <- function(valid_data, model, variables) {

    require(ggplot2)

    ranges <- apply(valid_data[,variables], 2, range, na.rm = T)
    means <- apply(valid_data[,variables], 2, mean, na.rm = T)

    resp_curves <- lapply(variables, function(x) NULL)

    for (i in seq_along(variables)) {
        to_pred <- data.frame(matrix(ncol = length(variables), nrow = 100))
        to_pred[,] <- rep(means, each = 100)
        to_pred[,i] <- seq(from = ranges[1,variables[i]], to = ranges[2,variables[i]], length.out = 100)
        colnames(to_pred) <- variables

        pred <- predict(model, to_pred, type = "response")

        resp_curves[[i]] <- data.frame(value = to_pred[,i], prediction = pred, variable = variables[i])
    }

    resp_curves <- do.call("rbind", resp_curves)

    ggplot(resp_curves) +
        geom_line(aes(x = value, y = prediction)) +
        facet_wrap(~variable, scales = "free") +
        theme_light()
}