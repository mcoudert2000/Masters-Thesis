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
load('results/estimate/isolated_tumor.rdata')

pure_tumor <- cbind(c_RNA1, c_RNA2)
pure_tumor <- pure_tumor[rownames(pure_tumor) %in% genes_I_care_about,]

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

#PCA with pure tumor
pca_plot_tum <- PCA(t(cbind(PCA_key_genes_data$E, pure_tumor)), scale.unit = T, graph = F)

d2 <- as.data.frame(pca_plot_tum$ind$coord)
d2$label <- c(colnames(PCA_key_genes_data), "C1", "C2")

#axes 2 and 3 w/ pure tumor #figure 15
fviz_pca_ind(pca_plot_tum, col.ind = F, geom.ind = F, geom.var = F, repel = T, axes = c(2,3)) + 
  geom_point(aes(x = d2$Dim.2, y = d2$Dim.3, col = c(PCA_key_genes_data$targets$source, "Wendler", "Wendler"))) + 
  geom_text_repel(data = d2[c(1:9, 399,400),], aes(x = Dim.2, y = Dim.3, label = label), hjust = -.1, vjust =-.1) +
  labs(col = "Source") +
  scale_color_manual(values = c(colors$CGGA, colors$TCGA, colors$Wendler)) +
  ggtitle("PCA plot of key genes with pure tumor included (2,3)")
  
#Axes 1 and 2 w/ pure tumor #figure 14
fviz_pca_ind(pca_plot_tum, col.ind = F, geom.ind = F, geom.var = F, repel = T, axes = c(1,2)) + 
  geom_point(aes(x = d2$Dim.1, y = d2$Dim.2, col = c(PCA_key_genes_data$targets$source, "Wendler", "Wendler"))) + 
  geom_text_repel(data = d2[c(1:9, 399,400),], aes(x = Dim.1, y = Dim.2, label = label), hjust = -.1, vjust =-.1) +
  labs(col = "Source") +
  scale_color_manual(values = c(colors$CGGA, colors$TCGA, colors$Wendler)) +
  ggtitle("PCA plot of key genes Dimensions: 1,2")



#### PCA with ALL genes and pure tumor
PCA_full_dataset <- combined_data_normalized

pca_full_plot <- PCA(t(cbind(PCA_full_dataset$E, (c_RNA1), (c_RNA2))), scale.unit = TRUE, graph = FALSE)

d3 <- as.data.frame(pca_full_plot$ind$coord)
d3$label <- c(colnames(PCA_key_genes_data), "C1", "C2")

#dim 2,3 figure 19
fviz_pca_ind(pca_full_plot, col.ind = F, geom.ind = F, geom.var = F, repel = T, axes = c(2,3)) + 
  geom_point(aes(x = d3$Dim.2, y = d3$Dim.3, col = c(PCA_key_genes_data$targets$source, "Wendler", "Wendler"))) + 
  geom_text_repel(data = d3[c(1:9, 399,400),], aes(x = Dim.2, y = Dim.3, label = label), hjust = -.1, vjust =-.1) +
  labs(col = "Source") +
  scale_color_manual(values = c(colors$CGGA, colors$TCGA, colors$Wendler)) +
  ggtitle("PCA Plot of All Genes with Pure Tumor, Dims: 2,3")

#dim 1,2 figure 18
require(ggrepel)
fviz_pca_ind(pca_full_plot, col.ind = F, geom.ind = F, geom.var = F, repel = T, axes = c(1,2)) + 
  geom_point(aes(x = d3$Dim.1, y = d3$Dim.2, col = c(PCA_key_genes_data$targets$source, "Wendler", "Wendler"))) + 
  #geom_text(data = d3[c(1:9, 399,400),], aes(x = Dim.1, y = Dim.2, label = label), hjust = -.1, vjust =-.1) +
  geom_text_repel(data = d3[c(1:9, 399,400),], aes(x = Dim.1, y = Dim.2, label = label))+
  labs(col = "Source") +
  scale_color_manual(values = c(colors$CGGA, colors$TCGA, colors$Wendler)) +
  ggtitle("PCA Plot of All Genes with Pure Tumor, Dims: 1,2")

#Investigating outliers on full plot
pca_outliers <- rownames(pca_full_plot$ind$coord[pca_full_plot$ind$coord[,'Dim.2'] > 100,])

#check if they are mostly wildtype or mutant (mostly mutant)
sum(combined_data_normalized$targets$IDH[rownames(combined_data_normalized$targets) %in% pca_outliers] == "Mutant") / length(pca_outliers)
sum(combined_data_normalized$targets$IDH[!rownames(combined_data_normalized$targets) %in% pca_outliers] == "Mutant") / (length(combined_data_normalized$targets[,1]) - length(pca_outliers))

t.test(combined_data_normalized$targets$IDH[rownames(combined_data_normalized$targets) %in% pca_outliers] == "Mutant",
       combined_data_normalized$targets$IDH[!rownames(combined_data_normalized$targets) %in% pca_outliers] == "Mutant")

#check if there's survival differences between groups
pca_outliers

require(survival)
require(survminer)
load('data/CGGA/CGGA_data.RDATA')
CGGA_data$samples$pca_outlier <- ifelse(rownames(CGGA_data$samples) %in% pca_outliers, 1, 0)

pca_outlier_survival_data <- data.frame(time = CGGA_data$samples$OS,
           event = CGGA_data$samples$Censor..alive.0..dead.1.,
           outlier = ifelse(rownames(CGGA_data$samples) %in% pca_outliers, T, F)) %>%
  dplyr::filter(!is.na(time))

ggsurvplot(survfit(Surv(time = pca_outlier_survival_data$time, 
                        event = pca_outlier_survival_data$event) ~ 
                     pca_outlier_survival_data$outlier, data =pca_outlier_survival_data), pval = T) +
  ggtitle("Survival Differences Between Outliers and Nonoutliers within the CGGA")



#CLUSTERING ON KEY GENES by ESTIMATE tumor purity

load('results/estimate_scores.rdata')

#Clustering on key genes colored by ESTIMATE tumor purity (dim1/2) 
#figure 16
fviz_pca_ind(pca_plot, col.ind = F, geom.ind = F, geom.var = F, repel = T, axes = c(1,2)) +
  geom_point(aes(x = d$Dim.1, y = d$Dim.2, col = estimate_scores$tumor_purity)) +
  geom_text_repel(data = d[1:9,], aes(x = Dim.1, y = Dim.2, label = label), hjust = -.1, vjust =-.1) +
  ggtitle("PCA Plot Colored by ESTIMATE Tumor Purity Dim: 1,2") +
  labs(col = "ESTIMATE tumor purity")

#Clustering on key genes colored by ESTIMATE tumor purity (dim2/3)
#figure 17
fviz_pca_ind(pca_plot, col.ind = F, geom.ind = F, geom.var = F, repel = T, axes = c(2,3)) +
  geom_point(aes(x = d$Dim.2, y = d$Dim.3, col = estimate_scores$tumor_purity)) +
  geom_text_repel(data = d[1:9,], aes(x = Dim.2, y = Dim.3, label = label), hjust = -.1, vjust =-.1) +
  ggtitle("PCA Plot Colored by ESTIMATE Tumor Purity Dim: 2,3") +
  labs(col = "ESTIMATE tumor purity")


#CLUSTER
t(PCA_key_genes_data$E)
d

#USE SIG-CLUST TO EVALUATE CLUSTER SIGNIFIGANCE
#Try sillouette as well
kclust_pca_data <- scale(t(PCA_key_genes_data$E))



gap_stat <- clusGap(kclust_pca_data,
                    FUN = kmeans,
                    nstart = 25,
                    K.max = 10,
                    B = 50)

#figure 20
set.seed(180018876)
annotate_figure(ggarrange(plotlist = list(silhoette = fviz_nbclust(kclust_pca_data, kmeans, method='silhouette') + labs(title = "Silhoette Width"), 
                          wss = fviz_nbclust(kclust_pca_data, kmeans, method='wss') + labs(title = "Elbow Selection"),
                          gap = fviz_nbclust(kclust_pca_data, kmeans, method='gap')), nrow = 1) + labs(title = "Gap Statistic"),
              top = "Results of Various Methods to Determine # Clusters for KMeans")


#plot number of clusters vs. gap statistic
fviz_gap_stat(gap_stat)
set.seed(2023)

kmeans2 <- kmeans(kclust_pca_data, centers = 2)
kmeans3 <- kmeans(kclust_pca_data, centers = 3)
kmeans4 <- kmeans(kclust_pca_data, centers = 4)
kplots2 <- fviz_cluster(kmeans2, data = kclust_pca_data, labelsize = 12*as.numeric(PCA_key_genes_data$targets$source == "Wendler"), repel = F)
kplots3 <- fviz_cluster(kmeans3, data = kclust_pca_data, labelsize = 12*as.numeric(PCA_key_genes_data$targets$source == "Wendler"), repel = F)
kplots4 <- fviz_cluster(kmeans4, data = kclust_pca_data, labelsize = 12*as.numeric(PCA_key_genes_data$targets$source == "Wendler"), repel = F)

sum(grepl("TCGA", rownames(dplyr::filter(data.frame(kmeans2$cluster), kmeans2.cluster == 1))))
sum(grepl("TCGA", rownames(dplyr::filter(data.frame(kmeans2$cluster), kmeans2.cluster == 2))))

sum(grepl("TCGA", rownames(dplyr::filter(data.frame(kmeans3$cluster), kmeans3.cluster == 3))))


require(ggpubr)
require(gridExtra)
kmeans_plot <- grid.arrange(grobs = list(kplots2, kplots3, kplots4), 
                         layout_matrix = matrix(c(1, 1, 2, 2, 4, 3, 3, 4), nrow = 2, byrow = TRUE),
                         top = "KMeans Clustering")
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



