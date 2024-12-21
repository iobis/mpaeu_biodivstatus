# PCA correlations




 data(iris)
 res.pca <- prcomp(iris[, -5],  scale = TRUE)
 # Extract the results for individuals
 ind <- get_pca_ind(res.pca)
 print(ind)





# Load necessary libraries
library(ggplot2)
library(dplyr)
library(mgcv)
library(factoextra)

ric <- read.csv("~/Downloads/full_data_plus_pts_AMA.csv")

ric_sel <- ric %>%
    select(-richness, -wf, -mpa_region, -zone)

pca <- prcomp(ric_sel,
             center = TRUE,
             scale = TRUE) 

pca$scale
print(pca)
summary(pca)

var <- factoextra::get_pca_ind(pca)
var$contrib

factoextra::fviz_pca_biplot(pca,
                geom.ind = "point",
                pointshape = 21, pointsize = 2.5, 
                fill.ind = data$mpa_region,  col.ind = "#FFFFFF", 
                palette = "simpson",
                alpha.var = 1, 
                col.var = "contrib", 
                gradient.cols = c("#999999","#0B0B0B","#9E1221"), 
                arrowsize = 1.0, labelsize = 6,
                legend.title=list(fill = "marine region", color = "contribution"), 
                title = "Principal Component Analysis", subtitle = "Full data by marine region",
                xlab="PC1 (39.5%)", ylab="PC2 (21.3%)", #check value in summary
                ) +
  theme(text = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14))

data <- ric %>% select(richness, zone, mpa_region) %>%
  mutate(PC1 = pca$x[, 1],
         PC2 = pca$x[, 2])

# Correlate PCA components with richness
cor(data$richness, data$PC1)
cor(data$richness, data$PC2)

# By zone
tapply(data[,c("richness", "PC1")], data$zone, function(x){
    cor(x[,1], x[,2])
})
tapply(data[,c("richness", "PC2")], data$zone, function(x){
    cor(x[,1], x[,2])
})

# By region
tapply(data[,c("richness", "PC1")], data$mpa_region, function(x){
    cor(x[,1], x[,2])
})
tapply(data[,c("richness", "PC2")], data$mpa_region, function(x){
    cor(x[,1], x[,2])
})


m2 <- gam(richness ~ s(PC1, bs = "cr") + s(PC2, bs = "cr"), data = data %>% slice_sample(n = 1000))
summary(m2)

ric_pc1 <- predict(m2, data.frame(
    PC1 = seq(min(data$PC1), max(data$PC1), length.out = 100),
    PC2 = mean(data$PC2)
))

plot(ric_pc1 ~ seq(min(data$PC1), max(data$PC1), length.out = 100))

ric_pc2 <- predict(m2, data.frame(
    PC2 = seq(min(data$PC2), max(data$PC2), length.out = 100),
    PC1 = mean(data$PC1)
))

plot(ric_pc2 ~ seq(min(data$PC2), max(data$PC2), length.out = 100))



var$contrib[order(var$contrib[,1], decreasing = T),]
var$contrib[order(var$contrib[,2], decreasing = T),]

fviz_contrib(pca, choice = "var", axes = 1, top = 10)
fviz_contrib(pca, choice = "var", axes = 2, top = 10)
