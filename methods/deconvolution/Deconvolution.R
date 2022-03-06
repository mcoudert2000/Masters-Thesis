#remotes::install_github("icbi-lab/immunedeconv")
require(immunedeconv)
require(ggplot2)
require(estimate)

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

library(help = "estimate")
estimate_data <- read.csv('data/combined_normalized_counts.csv')


removed_duplicates <- merge_duplicate_rows(estimate_data[,2:71])
removed_duplicates$GeneSymbol <- rownames(removed_duplicates)

#Normalize to log2
removed_duplicates[,1:69] <- log2(removed_duplicates[,1:69] + 1)

write.table(removed_duplicates, file = 'data/estimate_normalized_counts.txt', sep = '\t', quote = F)

filterCommonGenes(input.f = 'data/estimate_normalized_counts.txt', 
                  output.f = 'results/estimate/GBM_10412genes.gct',
                  id = 'GeneSymbol')

estimateScore("results/estimate/GBM_10412genes.gct", "results/estimate/GBM_estimate_score.gct", platform="illumina") #check to make sure this is the correct platform

#plotPurity(scores = "results/estimate/GBM_estimate_score.gct", platform = 'affymetrix')

estimate_scores <- read.delim(file="results/estimate/GBM_estimate_score.gct", skip=2, quote = " ")[,1:71]

estimate_scores <- t(estimate_scores)[2:71,]
colnames(estimate_scores) <- estimate_scores[1,]
estimate_scores <- as.data.frame(estimate_scores[2:70,])
estimate_scores[] <- lapply(estimate_scores, as.numeric)

estimate_scores$tumor_purity <- cos(estimate_scores$ESTIMATEScore *0.0001467884+0.6049872018) #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3826632/

#Plots
ggplot(estimate_scores[10:69,]) +
  geom_point(aes(x = StromalScore, y = ImmuneScore, col = 'TCGA')) +
  geom_point(data = estimate_scores[1:9,], aes(x = StromalScore, y = ImmuneScore, col = 'ASTRID\'s')) +
  geom_smooth(data = estimate_scores, aes(x = StromalScore, y = ImmuneScore), method = lm, se = T, linetype = 2, col = 'black') +
  ggtitle("Stromal Score vs Immune Score")

ggplot(estimate_scores[10:69,]) +
  geom_point(aes(x = 0, y = ESTIMATEScore, col = 'TCGA')) +
  geom_point(data = estimate_scores[1:9,], aes(x = 0, y = ESTIMATEScore, col = 'ASTRID\'s')) +
  geom_hline(yintercept = median(estimate_scores$ESTIMATEScore), linetype = 2)


ggplot(estimate_scores) +
  geom_boxplot(aes(x = tumor_purity)) +
  coord_flip()



### Linear Equation Solver

#Idea: RNA_tot = RNA_tumor * tumor_content + RNA_non_tumor * (1 - tumor_content)

lin_equation_data <- read.csv('data/Astrid_normalized_counts.csv')

rownames(lin_equation_data) <- lin_equation_data$X

lin_equation_data <- lin_equation_data[,3:11]
lin_equation_data$R2A

View(deconv_res_xcell)

lin_equation_data$R2C
estimate_scores[c(3,5),]$tumor_purity


tumor_rna <- t(solve(A,b))

colnames(tumor_rna) <- c("R2A", "R2C")
rownames(tumor_rna) <- rownames(lin_equation_data)

View(tumor_rna)

(lin_equation_data$R2A * estimate_scores[3,]$tumor_purity) / (2 * estimate_scores[3,]$tumor_purity - 1)
