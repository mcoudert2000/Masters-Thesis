#IDH mutant classifier
require("dplyr")
require("FactoMineR")
require("factoextra")
require("DESeq2")

mutant_list <- read.table('data/verhaak/IDH1-Mutated_Samples.txt', header = T)

transcription_data <- read.table("data/verhaak/unifiedScaled.txt", header = TRUE, row.names = 1, check.names = FALSE)



colnames(transcription_data) <- substr(colnames(transcription_data), start = 1, stop = 12)

#Check that all are covered
sum(colnames(transcription_data) %in% mutant_list$Tumor) / length(mutant_list$Tumor)

mutant_transcription_data <- transcription_data[colnames(transcription_data) %in% mutant_list$Tumor]
wild_type_transcription_data <- transcription_data[!colnames(transcription_data) %in% mutant_list$Tumor]

dim(wild_type_transcription_data)
dim(mutant_transcription_data) 
View(mutant_transcription_data)

pca_idh_data <-data.frame(t(transcription_data))
pca_idh_data$sample <- rownames(pca_idh_data)
pca_idh_data$idh_type <- ifelse(colnames(transcription_data) %in% mutant_list$Tumor, "Mutant", "Wild-Type")

pca_idh_data <- dplyr::select(pca_idh_data, c(sample, idh_type, everything()))

pca_plot <- PCA(pca_idh_data[,3:11863], scale.unit = TRUE, graph = FALSE)
pca_plot

fviz_pca_ind(pca_plot, col.ind = pca_idh_data$idh_type, label = F, axes = c(1,2))

require(ggplot2)
pca_idh_data %>% ggplot() +
  geom_point(aes(x = IDH1, y = ifelse(idh_type == 'Mutant', 0, 0.1), col = idh_type))

