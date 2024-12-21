
library(gplots)
library(ggplot2) 
library(gridExtra)
library(factoextra)
library(dplyr)
library(tidyr)
library(ggpubr)
library(corrplot)
library(heatmaply)

####PREPARATION DATASET####
head(richness)

richness_new <- richness 
colnames(richness_new) <- richness[1, ]
richness_new 
head(richness_new)
richness_new <- richness_new[-1, ] 
richness_new 
head(richness_new)

is.numeric(richness_new)
print(typeof(richness_new))
sapply(richness_new, class)

richness_new$richnessrf2 <- as.numeric(as.character(richness_new$richnessrf2))
richness_new$wf <- as.numeric(as.character(richness_new$wf))
richness_new$slp <- as.numeric(as.character(richness_new$slp))
richness_new$tri <- as.numeric(as.character(richness_new$tri))
richness_new$rgs <- as.numeric(as.character(richness_new$rgs))
richness_new$dtc <- as.numeric(as.character(richness_new$dtc))
richness_new$btm <- as.numeric(as.character(richness_new$btm))
richness_new$sst <- as.numeric(as.character(richness_new$sst))
richness_new$tas <- as.numeric(as.character(richness_new$tas))
richness_new$sws <- as.numeric(as.character(richness_new$sws))
richness_new$sln <- as.numeric(as.character(richness_new$sln))
richness_new$sic <- as.numeric(as.character(richness_new$sic))
richness_new$phy <- as.numeric(as.character(richness_new$phy))
richness_new$ph <- as.numeric(as.character(richness_new$ph))
richness_new$par <- as.numeric(as.character(richness_new$par))
richness_new$oxy <- as.numeric(as.character(richness_new$oxy))
richness_new$ntr <- as.numeric(as.character(richness_new$ntr))
richness_new$kdp <- as.numeric(as.character(richness_new$kdp))
richness_new$chl <- as.numeric(as.character(richness_new$chl))
richness_new$lon <- as.numeric(as.character(richness_new$lon))
richness_new$lat <- as.numeric(as.character(richness_new$lat))
sapply(richness_new, class)

###CORRELOGRAM AND DENDROGRAM####
richness_new <- richness_new[,-1] 
head(richness_new)
richness_matrix <- cor(richness_new)
head(round(richness_matrix,2))

dist <- dist(richness_matrix[ , c(1:20)] , diag=TRUE)
hc <- hclust(dist)
plot(hc)

col_matrix <- colorRampPalette(c("#4477AA","#77AADD","#FFFFFF", "#EE9988", "#BB4444")) (3)
heatmap(x=richness_matrix, col=col_matrix,  symm=TRUE, cexRow = 1.5, cexCol = 1.5)                              
legend(x="topleft", legend=c("low", "medium", "high"), fill=col_matrix, lty=1:2, cex=1, box.lty=0, bg=NA) #no border and backgroudn legendbox 

corrplot(richness_matrix, type="full", order="hclust")
corrplot(richness_matrix, type="lower", order="hclust")

cor.mtest <- function(richness_matrix, ...) {
  richness_matrix <- as.matrix(richness_matrix)
  n <- ncol(richness_matrix)
  p.richness_matrix<- matrix(NA, n, n)
  diag(p.richness_matrix) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(richness_matrix[, i], richness_matrix[, j], ...)
      p.richness_matrix[i, j] <- p.richness_matrix[j, i] <- tmp$p.value
    }
  }
  colnames(p.richness_matrix) <- rownames(p.richness_matrix) <- colnames(richness_matrix)
  p.richness_matrix
}
p.richness_matrix <- cor.mtest(richness_new)
head(p.richness_matrix[, 1:21])

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(richness_matrix, method="color", col=col(200),  
         type="lower", order="hclust", 
         addCoef.col = "black", 
         number.cex=0.70, 
         tl.cex = 1.2, 
         tl.col="black", tl.srt=45, 
                  p.richness_matrix = p.richness_matrix, sig.level = 0.01, insig = "blank", 
         diag=FALSE)

####SCATTERPLOT####
head(richness_new)
p1 <- ggplot(richness_new, aes(x=sst, y=hill_q0)) + 
  geom_point(size=3, color="#339933", alpha = 0.7) +
  xlab("Sea surface temperature") + ylab("Hill q=0")  + 
  ggtitle("") +
  theme(text = element_text(size=15))
p2 <- ggplot(richness_new, aes(x=sst, y=es_50)) + 
  geom_point(size=3, color="#339933", alpha = 0.7) +
  xlab("Sea surface temperature") + ylab("ES50")  + 
  ggtitle("") +
  theme(text = element_text(size=15))

p3 <- ggplot(richness_new, aes(x=sic, y=hill_q0)) + 
  geom_point(size=3, color="#339933", alpha = 0.7) +
  xlab("Sea ice concentration") + ylab("Hill q=0")  + 
  ggtitle("") +
  theme(text = element_text(size=15))
p4 <- ggplot(richness_new, aes(x=sic, y=es_50)) + 
  geom_point(size=3, color="#339933", alpha = 0.7) +
  xlab("Sea ice concentration") + ylab("ES50")  + 
  ggtitle("") +
  theme(text = element_text(size=15))

p5 <- ggplot(richness_new, aes(x=par, y=hill_q0)) + 
  geom_point(size=3, color="#339933", alpha = 0.7) +
  xlab("Photosynthetically active radiation") + ylab("Hill q=0")  + 
  ggtitle("") +
  theme(text = element_text(size=15))
p6 <- ggplot(richness_new, aes(x=par, y=es_50)) + 
  geom_point(size=3, color="#339933", alpha = 0.7) +
  xlab("Photosynthetically active radiation") + ylab("ES50")  + 
  ggtitle("") +
  theme(text = element_text(size=15))

grid.arrange(p1, p2, p3, p4, p5, p6, ncol = 2, nrow = 3) #grid.arrange(p1, p2,...  = 2, nrow = 3)

grid.arrange(p1, p2, ncol = 2, nrow = 1) #grid.arrange(p1, p2,...  = 2, nrow = 3)
grid.arrange(p1, p2, ncol = 2, nrow = 1) #grid.arrange(p1, p2,...  = 2, nrow = 3)

#by depth zone
cols <- c("#0073C2FF", "#EFC000FF")
p1 <- ggplot(richness_new, aes(x=rgs, y=es_50, color = depth)) + 
  geom_point(size=3,alpha = 0.7) +
  scale_colour_manual(values = cols) +
  xlab("Rugosity") + ylab("Species richness (ES50)")  + 
  ggtitle("") + #+ theme(legend.position="none") #remove legend for the final combined grid plot
  theme(text = element_text(size=15))

p2 <- ggplot(richness_new, aes(x=rgs, y=hill_q0, color = depth)) + 
  geom_point(size=3,alpha = 0.7) +
  scale_colour_manual(values = cols) +
  xlab("Rugosity") + ylab("Species richness (Hill q=0)")  + 
  ggtitle("")+
  theme(text = element_text(size=15))

grid.arrange(p1, p2, ncol = 2, nrow = 1) #grid.arrange(p1, p2,...  = 2, nrow = 3)



####BOXPLOT####
head(richness_new)
legend_title <- "depth zone"
bxp1 <- ggplot(richness_new, aes(x=btm, y=mpa_region, fill=zone)) + 
  geom_boxplot(outlier.shape = NA) +
    labs(title="",x="bathymetry", y = "marine region") +
  scale_fill_manual(legend_title, values=c("#0073C2FF", "#EFC000FF")) + 
  coord_flip() + 
  theme(text = element_text(size=15))
bxp2 <- ggplot(richness_new, aes(x=dtc, y=mpa_region, fill=zone)) + 
  geom_boxplot(outlier.shape = NA) +
  labs(title="",x="distance to coast", y = "marine region") +
  scale_fill_manual(legend_title, values=c("#0073C2FF", "#EFC000FF")) + 
  coord_flip() + 
  theme(text = element_text(size=15))
bxp3 <- ggplot(richness_new, aes(x=lat, y=mpa_region, fill=zone)) + 
  geom_boxplot(outlier.shape = NA) +
  labs(title="",x="latitude", y = "marine region") +
  scale_fill_manual(legend_title, values=c("#0073C2FF", "#EFC000FF")) + 
  coord_flip() + 
  theme(text = element_text(size=15))
bxp4 <- ggplot(richness_new, aes(x=lon, y=mpa_region, fill=zone)) + 
  geom_boxplot(outlier.shape = NA) +
  labs(title="",x="longitude", y = "marine region") +
  scale_fill_manual(legend_title, values=c("#0073C2FF", "#EFC000FF")) + 
  coord_flip() + 
  theme(text = element_text(size=15))
bxp5 <- ggplot(richness_new, aes(x=par, y=mpa_region, fill=zone)) + 
  geom_boxplot(outlier.shape = NA) +
  labs(title="",x="photosynthetically active radiation", y = "marine region") +
  scale_fill_manual(legend_title, values=c("#0073C2FF", "#EFC000FF")) + 
  coord_flip() + 
  theme(text = element_text(size=15))
bxp6 <- ggplot(richness_new, aes(x=rgs, y=mpa_region, fill=zone)) + 
  geom_boxplot(outlier.shape = NA) +
  labs(title="",x="rugosity", y = "marine region") +
  scale_fill_manual(legend_title, values=c("#0073C2FF", "#EFC000FF")) + 
  coord_flip() + 
  xlim(0,500) + 
  theme(text = element_text(size=15))
bxp7 <- ggplot(richness_new, aes(x=sic, y=mpa_region, fill=zone)) + 
  geom_boxplot(outlier.shape = NA) +
  labs(title="",x="sea ice concentration", y = "marine region") +
  scale_fill_manual(legend_title, values=c("#0073C2FF", "#EFC000FF")) + 
  coord_flip() + 
  theme(text = element_text(size=15))
bxp8 <- ggplot(richness_new, aes(x=sln, y=mpa_region, fill=zone)) + 
  geom_boxplot(outlier.shape = NA) +
  labs(title="",x="salinity", y = "marine region") +
  scale_fill_manual(legend_title, values=c("#0073C2FF", "#EFC000FF")) + 
  coord_flip() +
  theme(text = element_text(size=15))
bxp9 <- ggplot(richness_new, aes(x=sst, y=mpa_region, fill=zone)) + 
  geom_boxplot(outlier.shape = NA) +
  labs(title="",x="sea surface temperature", y = "marine region") +
  scale_fill_manual(legend_title, values=c("#0073C2FF", "#EFC000FF")) + 
  coord_flip() + 
  theme(text = element_text(size=15))
bxp10 <- ggplot(richness_new, aes(x=wf, y=mpa_region, fill=zone)) + 
  geom_boxplot(outlier.shape = NA) +
  labs(title="",x="wave fetch", y = "marine region") +
  scale_fill_manual(legend_title, values=c("#0073C2FF", "#EFC000FF")) + 
  coord_flip() + 
  theme(text = element_text(size=15))

grid.arrange(bxp1, bxp2, ncol = 1, nrow = 2) 
grid.arrange(bxp3, bxp4, ncol = 1, nrow = 2) 
grid.arrange(bxp5, bxp6, ncol = 1, nrow = 2) 
grid.arrange(bxp7, bxp8, ncol = 1, nrow = 2) 
grid.arrange(bxp9, bxp10, ncol = 1, nrow = 2) 

####PCA####
quantiles <- quantile(richness_new$richness, probs =c(0,0.25,0.5,0.75,1))
richness_new$richness <- cut(richness_new$richness, 
                       breaks = quantiles, 
                       include.lowest =TRUE, 
                       labels =c("Q1","Q2","Q3","Q4"))
head(richness_new)

pca <-prcomp(richness_new [,-c(1,11,12,13)], # remove richness, wf, zone and region etc
             center=TRUE,
             scale = TRUE) 

pca$scale
print(pca)
summary(pca) 

fviz_pca_biplot(pca,
                geom.ind = "point",
                pointshape = 21, pointsize = 2.5, 
                fill.ind = richness_new$richness,  col.ind = "#FFFFFF", 
                palette = "rickandmorty",
                alpha.var = 1, 
                col.var = "contrib", 
                gradient.cols = c("#999999","#0B0B0B","#9E1221"), 
                arrowsize = 1.0, labelsize = 6,
                legend.title=list(fill = "richness", color = "contribution"), 
                title = "Principal Component Analysis", subtitle = "Full data by richness",
                xlab="PC1 (39.5%)", ylab="PC2 (21.3%)", #check value in summary
) +
  theme(text = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14))


fviz_pca_biplot(pca,
                geom.ind = "point",
                pointshape = 21, pointsize = 2.5, 
                fill.ind = richness_new$mpa_region,  col.ind = "#FFFFFF", 
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

fviz_pca_biplot(pca,
                geom.ind = "point",
                pointshape = 21, pointsize = 2.5, 
                fill.ind = richness_new$zone,  col.ind = "#FFFFFF", 
                palette = "jco",
                alpha.var = 1, 
                col.var = "contrib", 
                gradient.cols = c("#999999","#0B0B0B","#9E1221"), 
                arrowsize = 1.0, labelsize = 6,
                legend.title=list(fill = "depth zone", color = "contribution"), 
                title = "Principal Component Analysis", subtitle = "Full data by depth zone",
                xlab="PC1 (39.5%)", ylab="PC2 (21.3%)", #check value in summary
                ) +
  theme(text = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14))
