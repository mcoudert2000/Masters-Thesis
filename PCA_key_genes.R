require("FactoMineR")
require("factoextra")


PCA_data <- read.csv('data/combined_normalized_counts_filtered.csv')
rownames(PCA_data) <- PCA_data$GENENAME
PCA_data <- PCA_data %>%
  dplyr::select(-c("X", "GENENAME")) %>% scale() %>% t()

pca_plot <- PCA(PCA_data, scale.unit = TRUE, graph = FALSE)

d <- as.data.frame(pca_plot$ind$coord)
d$label <- rownames(PCA_data)

fviz_pca_ind(pca_plot, col.ind = 'red', geom.ind = 'point', geom.var = F, repel = T) +
  geom_text(data = d[1:9,], aes(x = Dim.1, y = Dim.2, label = label), hjust = -.1, vjust =-.1)

#CLUSTER

#USE SIG-CLUST TO EVALUATE CLUSTER SIGNIFIGANCE
#Try sillouette as well
kplots2 <- fviz_cluster(kmeans(PCA_data, centers = 2), data = PCA_data, labelsize = 0) +
  geom_text(data = d[1:9,], aes(x = Dim.1, y = Dim.2, label = label))
kplots3 <- fviz_cluster(kmeans(PCA_data, centers = 3), data = PCA_data, labelsize = 0) +
  geom_text(data = d[1:9,], aes(x = Dim.1, y = Dim.2, label = label))
kplots4 <- fviz_cluster(kmeans(PCA_data, centers = 4), data = PCA_data,labelsize = 0) +
  geom_text(data = d[1:9,], aes(x = Dim.1, y = Dim.2, label = label))
kplots5 <- fviz_cluster(kmeans(PCA_data, centers = 5), data = PCA_data, labelsize = 0) +
  geom_text(data = d[1:9,], aes(x = Dim.1, y = Dim.2, label = label))

require(ggpubr)
kmeans_plot <- ggarrange(kplots2, kplots3, kplots4, kplots5, labels = c("2","3","4","5"))

#save plot
png(filename = "plots/KMeans Clustering on Results from Various Gene Signature Studies.png",
    width = 800, height = 400)
plot(kmeans_plot, main = "Kmeans Clustering on Results from Various Gene Signature Studies")
dev.off()

