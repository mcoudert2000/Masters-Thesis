#remotes::install_github("icbi-lab/immunedeconv")
require(immunedeconv)
require(ggplot2)
require(estimate)

set_cibersort_binary("/Users/matthewcoudert/Maths/2021:22/Masters-Thesis/data/Cibersort/CIBERSORT.R")
set_cibersort_mat("/Users/matthewcoudert/Maths/2021:22/Masters-Thesis/data/Cibersort/LM22.txt")

#Load full dataset
load('data/combined_data.rdata')

deconvolution_data <- as.matrix(combined_data_normalized$E)

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

#library(help = "estimate")
load('data/combined_data.rdata')
estimate_data <- combined_data_normalized$E

write.table(estimate_data, file = 'data/estimate_normalized_counts.txt', sep = '\t', quote = F)


filterCommonGenes(input.f = 'data/estimate_normalized_counts.txt', 
                  output.f = 'results/estimate/GBM_10412genes.gct',
                  id = 'GeneSymbol')

estimateScore("results/estimate/GBM_10412genes.gct", "results/estimate/GBM_estimate_score.gct", platform="illumina") #check to make sure this is the correct platform

#plotPurity(scores = "results/estimate/GBM_estimate_score.gct", platform = 'affymetrix')

estimate_scores <- read.delim(file="results/estimate/GBM_estimate_score.gct", skip=2, quote = " ")[,1:400]

estimate_scores <- t(estimate_scores)
colnames(estimate_scores) <- estimate_scores[1,]
estimate_scores <- as.data.frame(estimate_scores[3:400,])
estimate_scores[] <- lapply(estimate_scores, as.numeric)

estimate_scores$tumor_purity <- cos(estimate_scores$ESTIMATEScore *0.0001467884+0.6049872018) #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3826632/

dim(estimate_scores)
#Plots
ggplot(estimate_scores) +
  geom_point(data = estimate_scores[10:248,], aes(x = StromalScore, y = ImmuneScore, col = 'CGGA')) +
  geom_point(data = estimate_scores[1:9,], aes(x = StromalScore, y = ImmuneScore, col = 'ASTRID\'s')) +
  geom_point(data = estimate_scores[249:398,], aes(x = StromalScore, y = ImmuneScore, col = 'TCGA')) +
  geom_smooth(data = estimate_scores, aes(x = StromalScore, y = ImmuneScore), method = lm, se = T, linetype = 2, col = 'black') +
  ggtitle("Stromal Score vs Immune Score")

ggplot(estimate_scores[10:248,]) +
  geom_point(aes(x = -1, y = ESTIMATEScore, col = 'CGGA'), alpha =0.3) +
  geom_point(data = estimate_scores[1:9,], alpha = 0.3, aes(x = 0, y = ESTIMATEScore, col = 'ASTRID\'s')) +
  geom_point(data = estimate_scores[249:398,],alpha =0.3, aes(x = 1, y = ESTIMATEScore, col = 'TCGA')) +
  geom_hline(yintercept = median(estimate_scores$ESTIMATEScore), linetype = 2)


ggplot(estimate_scores) +
  geom_boxplot(aes(x = tumor_purity)) +
  coord_flip() 


### Linear Equation Solver

#Idea: RNA_tot = RNA_tumor * tumor_content + RNA_non_tumor * (1 - tumor_content)


R2A <- combined_data_normalized$E[,'R2A']
R2C <- combined_data_normalized$E[,'R2C']

R2A_tp <- estimate_scores[c(3),]$tumor_purity
R2C_tp <- estimate_scores[c(5),]$tumor_purity

R2A_tp - R2C_tp
c_RNA2 <- ( R2A - R2C + R2A_tp * R2C - R2C_tp * R2A ) / (R2A_tp - R2C_tp)

#This is (theoretically) isolated cancer RNA
c_RNA2

