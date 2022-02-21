#remotes::install_github("icbi-lab/immunedeconv")
require(immunedeconv)

set_cibersort_binary("/Users/matthewcoudert/Maths/2021:22/Masters-Thesis/data/Cibersort/CIBERSORT.R")
set_cibersort_mat("/Users/matthewcoudert/Maths/2021:22/Masters-Thesis/data/Cibersort/LM22.txt")

#Load full dataset
deconvolution_data <- as.matrix(read.csv('data/combined_normalized_counts.csv'))
rownames(deconvolution_data) <- deconvolution_data[,'GENENAME']
deconvolution_data <- deconvolution_data[,3:71]
#View(deconvolution_data)

deconv_res_cibersort <- deconvolute(deconvolution_data, method = 'cibersort')

deconv_res_xcell <- deconvolute(deconvolution_data, method = 'xcell')

immune_stroma_scores <- as.data.frame(deconv_res_xcell[37:38,])

rownames(immune_stroma_scores) <- immune_stroma_scores$cell_type
immune_stroma_scores <- dplyr::select(immune_stroma_scores, -cell_type)

meds <- c(median(as.numeric(immune_stroma_scores[1,])), median(as.numeric(immune_stroma_scores[2,])))
binary_results <- immune_stroma_scores
for(i in 1:2) {
  binary_results[i,] <- immune_stroma_scores[i,] > meds[i]
}
binary_results[1:9]


#Not finished yet
rforge <- "http://r-forge.r-project.org"
install.packages("estimate", repos=rforge, dependencies=TRUE)
require(estimate)
library(help = "estimate")

deconvolution_data
?filterCommonGenes
estimate_data <- read.csv('data/combined_normalized_counts.csv')


removed_duplicates <- merge_duplicate_rows(estimate_data[,2:71])
removed_duplicates$GeneSymbol <- rownames(removed_duplicates)

write.table(removed_duplicates, file = 'data/estimate_normalized_counts.txt', sep = '\t', quote = F)
eset <- ExpressionSet(assayData = as.matrix(removed_duplicates))
View(eset)

filterCommonGenes(input.f = 'data/estimate_normalized_counts.txt', 
                  output.f = 'results/estimate/GBM_10412genes.gct',
                  id = 'GeneSymbol')

estimateScore("results/estimate/GBM_10412genes.gct", "results/estimate/GBM_estimate_score.gct", platform="illumina") #check to make sure this is the correct platform

plotPurity(scores = "results/estimate/GBM_estimate_score.gct", platform = 'affymetrix')

estimate_scores <- read.delim(file="results/estimate/GBM_estimate_score.gct", skip=2, quote = " ")[,1:71]

estimate_scores <- t(estimate_scores)[2:71,]
colnames(estimate_scores) <- estimate_scores[1,]
estimate_scores <- as.data.frame(estimate_scores[2:70,])
estimate_scores[] <- lapply(estimate_scores, as.numeric)

require(gghighlight)

ggplot(estimate_scores[10:69,]) +
  geom_point(aes(x = StromalScore, y = ImmuneScore, col = 'TCGA')) +
  geom_point(data = estimate_scores[1:9,], aes(x = StromalScore, y = ImmuneScore, col = 'ASTRID\'s'))



