#MAIN

#Run this before running everything else to generate cleaned data sets 
#and create all variables that will be used

source("methods/Helper_functions.R")

#Data cleaning
source("methods/data_cleaning/Astrid_data_cleaning.R")
source("methods/data_cleaning/TCGA_data_cleaning.R")
source("methods/data_cleaning/CGGA_data_cleaning.R")
source("methods/data_cleaning/Combined_data_cleaning.R")

#Deconvolution
source("methods/deconvolution/Deconvolution.R")

#Gene signature based studies
source("methods/previous_gene_signatures/Prev_methods_calculation.R")
source("methods/previous_gene_signatures/Ensemble_prev_methods.R")
source("methods/previous_gene_signatures/PCA_key_genes.R")
source("methods/previous_gene_signatures/Dendrogram_prev_methods_results.R")

#Classifers
source("methods/classifiers/Verhaak_Classifier.R")
source("methods/classifiers/Differential_expression_CGGA_TCGA.R")
source("methods/classifiers/IDH_binary_classifier.R")

