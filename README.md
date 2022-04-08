# Masters-Thesis

The methods written are all within the "methods folder" with the following breakdown (listed in the order they appear within my dissertation).

I would suggest running source("methods/Main.R") before exploring the rest of the code, as it will create all the datasets and load objects needed to run the rest of the code. 

Data used from the following sources is available as follows:

TCGA: gbm.exp
CGGA: CGGA.mRNAseq_693.RSEM-genes.20200506.txt (expression data),
CGGA.mRNAseq_693.RSEM_clinical.20200506.txt (clinical information about patients)
Wendler: GeneCounts.RData

data_cleaning:
  - TCGA_data_cleaning2.R
  - Astrid_data_cleaning.R
  - CGGA_data_cleaning.R
  - Combined_data_cleaning2.R

deconvolution:
  - Deconvolution.R

previous_gene_signatures:
  - Prev_methods_calculation.R
  - Ensemble_prev_methods.R
  - PCA_key_genes.R
  - Dendrogram_prev_methods_results.R

classifiers:
  - Verhaak_classifier.R
  - Differential_expression_CGGA_TCGA.R
  - IDH_binary_classifier.R

