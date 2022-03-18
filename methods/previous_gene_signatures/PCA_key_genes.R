require("FactoMineR")
require("factoextra")
require(dplyr)
require(cluster)

GR_Sig <- c('GRIA2', 'RYR3')
NIM_Sig <- c('NRG1', 'ITGA3', 'MAP1LC3A')
HOSS_Sig <- c('HOXC10', 'OSMR', 'SCARA3', 'SLC39A10')
PRGIT_Sig <- c('PTPRN', 'RGS14', 'G6PC3', 'IGFBP2', 'TIMP4')
BHLNSX_Sig <- c('BPIFB2', 'HOXA13', 'LRRC10', 'NELL1', 'SDR16C5', 'XIRP2')
genes_I_care_about <- c(GR_Sig, NIM_Sig, HOSS_Sig, PRGIT_Sig, BHLNSX_Sig)

load('data/combined_data.rdata')

PCA_key_genes_data <- combined_data_normalized
PCA_key_genes_data$E <- PCA_key_genes_data$E[rownames(combined_data_normalized$E) %in% genes_I_care_about,]
rownames(PCA_key_genes_data$E)

pca_plot <- PCA(t(PCA_key_genes_data$E), scale.unit = TRUE, graph = FALSE)

d <- as.data.frame(pca_plot$ind$coord)
d$label <- colnames(PCA_key_genes_data)

fviz_pca_ind(pca_plot, col.ind = F, geom.ind = F, geom.var = F, repel = T, axes = c(1,2)) +
  geom_point(aes(x = d$Dim.1, y = d$Dim.2, col = PCA_key_genes_data$targets$source)) +
  geom_text(data = d[1:9,], aes(x = Dim.1, y = Dim.2, label = label), hjust = -.1, vjust =-.1)


fviz_pca_ind(pca_plot, col.ind = F, geom.ind = F, geom.var = F, repel = T, axes = c(2,3)) +
  geom_point(aes(x = d$Dim.2, y = d$Dim.3, col = PCA_key_genes_data$targets$source)) +
  geom_text(data = d[1:9,], aes(x = Dim.2, y = Dim.3, label = label), hjust = -.1, vjust =-.1)


load('results/estimate_scores.rdata')

#Clustering on key genes colored by ESTIMATE tumor purity (dim1/2)
fviz_pca_ind(pca_plot, col.ind = F, geom.ind = F, geom.var = F, repel = T, axes = c(1,2)) +
  geom_point(aes(x = d$Dim.1, y = d$Dim.2, col = estimate_scores$tumor_purity)) +
  geom_text(data = d[1:9,], aes(x = Dim.1, y = Dim.2, label = label), hjust = -.1, vjust =-.1)

#Clustering on key genes colored by ESTIMATE tumor purity (dim2/3)
fviz_pca_ind(pca_plot, col.ind = F, geom.ind = F, geom.var = F, repel = T, axes = c(2,3)) +
  geom_point(aes(x = d$Dim.2, y = d$Dim.3, col = estimate_scores$tumor_purity)) +
  geom_text(data = d[1:9,], aes(x = Dim.2, y = Dim.3, label = label), hjust = -.1, vjust =-.1)


#CLUSTER
t(PCA_key_genes_data$E)
d

#USE SIG-CLUST TO EVALUATE CLUSTER SIGNIFIGANCE
#Try sillouette as well
kclust_pca_data <- scale(t(PCA_key_genes_data$E))

fviz_nbclust(kclust_pca_data, kmeans, method = "wss")

gap_stat <- clusGap(kclust_pca_data,
                    FUN = kmeans,
                    nstart = 25,
                    K.max = 10,
                    B = 50)

fviz_gap_stat(gap_stat)


#plot number of clusters vs. gap statistic
fviz_gap_stat(gap_stat)

kplots2 <- fviz_cluster(kmeans(kclust_pca_data, centers = 2), data = kclust_pca_data, labelsize = 10*as.numeric(PCA_key_genes_data$targets$source == "Astrid"), repel = F)
kplots3 <- fviz_cluster(kmeans(kclust_pca_data, centers = 3), data = kclust_pca_data, labelsize = 10*as.numeric(PCA_key_genes_data$targets$source == "Astrid"), repel = F)
kplots4 <- fviz_cluster(kmeans(kclust_pca_data, centers = 4), data = kclust_pca_data, labelsize = 10*as.numeric(PCA_key_genes_data$targets$source == "Astrid"), repel = F)

kplots9 <- fviz_cluster(kmeans(kclust_pca_data, centers = 9), data = kclust_pca_data, labelsize = 10*as.numeric(PCA_key_genes_data$targets$source == "Astrid"), repel = F) 
kplots9

require(ggpubr)
kmeans_plot <- ggarrange(kplots2, kplots3, kplots4, kplots9, labels = c("2","3","4","9"))
kmeans_plot
#save plot
png(filename = "plots/KMeans Clustering on Results from Various Gene Signature Studies.png",
    width = 800, height = 400)
plot(kmeans_plot, main = "Kmeans Clustering on Results from Various Gene Signature Studies")
dev.off()

#R1A          R1B          R2A          R2B          R2C 
#1            7            7            7            7 
#R1C           PA           PB           PC 
#7            6            1            1 

