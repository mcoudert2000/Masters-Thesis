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
