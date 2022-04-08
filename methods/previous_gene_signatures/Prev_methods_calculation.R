require(dplyr)
GR_Sig <- c('GRIA2', 'RYR3')
NIM_Sig <- c('NRG1', 'ITGA3', 'MAP1LC3A')
HOSS_Sig <- c('HOXC10', 'OSMR', 'SCARA3', 'SLC39A10')
PRGIT_Sig <- c('PTPRN', 'RGS14', 'G6PC3', 'IGFBP2', 'TIMP4')
BHLNSX_Sig <- c('BPIFB2', 'HOXA13', 'LRRC10', 'NELL1', 'SDR16C5', 'XIRP2')
genes_I_care_about <- c(GR_Sig, NIM_Sig, HOSS_Sig, PRGIT_Sig, BHLNSX_Sig)

#Load in data
load('data/combined_data.rdata')

#isolate dataset to only genes within the classifiers
key_genes_data <- data.frame(combined_data_normalized$E[rownames(combined_data_normalized$E) %in% genes_I_care_about,])

#get cpm data rather than log2 cpm to calculate surival for both as it is not specified in papers
key_genes_cpm_data <- data.frame(cpm(combined_data$counts[rownames(combined_data) %in% genes_I_care_about,]))

key_genes_data$GENENAME <- rownames(key_genes_data)
key_genes_cpm_data$GENENAME <- rownames(key_genes_cpm_data)

### Prev Methods -----
#caculates the summation for each of the samples
prev_method_calculator <- function(signature, coef, dat) {
  out <- 0
  for(g in signature) {
    out <- out + as.numeric(coef[g]) * dat %>% dplyr::filter(GENENAME == g) %>%
      dplyr::select(-GENENAME)
  }
  return(out)
}
#### GR-SIG ------
GR_Sig
GR_coef <- list(GRIA2 = -0.1, RYR3 = 0.1)
GR_sig_out <- prev_method_calculator(GR_Sig, GR_coef, key_genes_data)
GR_sig_out_cpm <- prev_method_calculator(GR_Sig, GR_coef, key_genes_cpm_data)


#### HOSS-Sig-----
HOSS_Sig
HOSS_coef <- list(HOXC10 = 0.089, OSMR = 0.238, SCARA3 = 0.238, SLC39A10 = -0.424) 
#This is from table 1 of paper
HOSS_sig_out <- prev_method_calculator(HOSS_Sig, HOSS_coef, key_genes_data)
HOSS_sig_out_cpm <- prev_method_calculator(HOSS_Sig, HOSS_coef, key_genes_cpm_data)

#### PRGIT-Sig------
PRGIT_Sig
PRGIT_coef <- list(PTPRN = 0.50894, RGS14 = 0.54671, 
                   G6PC3 = 1.20753, IGFBP2 = 0.25845, TIMP4 = -0.20684)
PRGIT_sig_out <- prev_method_calculator(PRGIT_Sig, PRGIT_coef, key_genes_data)
PRGIT_sig_out_cpm <- prev_method_calculator(PRGIT_Sig, PRGIT_coef, key_genes_cpm_data)


#### BHLNSX-Sig------
BHLNSX_Sig
BHLNSX_coef <- list(BPIFB2 = -0.0881, HOXA13 = -0.0854, LRRC10 = -0.0614, NELL1 = 0.151, SDR16C5 = 0.0945, XIRP2 = -0.1441)
BHLNSX_sig_out <- prev_method_calculator(BHLNSX_Sig, BHLNSX_coef, key_genes_data)
BHLNSX_sig_out_cpm <- prev_method_calculator(BHLNSX_Sig, BHLNSX_coef, key_genes_cpm_data)

#### NIM ------
#NRG1 × 0.132 + expression level of ITGA3 × 0.139 + expression level of MAP1LC3A × 0.269
NIM_Sig
NIM_coef <- list(NRG1 = 0.132, ITGA3 = 0.139, MAP1LC3A = 0.269)
NIM_sig_out <- prev_method_calculator(NIM_Sig, NIM_coef, key_genes_data)
NIM_sig_out_cpm <- prev_method_calculator(NIM_Sig, NIM_coef, key_genes_cpm_data)

#save results to dataframe
prev_methods_results <- data.frame(GR = t(GR_sig_out), NIM = t(NIM_sig_out), HOSS = t(HOSS_sig_out),
                                   PRGIT = t(PRGIT_sig_out), BHLNSX = t(BHLNSX_sig_out))

prev_methods_results_cpm <- data.frame(GR = t(GR_sig_out_cpm), NIM = t(NIM_sig_out_cpm), HOSS = t(HOSS_sig_out_cpm),
                                   PRGIT = t(PRGIT_sig_out_cpm), BHLNSX = t(BHLNSX_sig_out_cpm))


colnames(prev_methods_results) <- c("GR","NIM", "HOSS", "PRGIT", "BHLNSX")
colnames(prev_methods_results_cpm) <- c("GR","NIM", "HOSS", "PRGIT", "BHLNSX")

#running previous methods on deconvoluted samples
load('results/estimate/isolated_tumor.rdata')
key_genes_pure_tumor <- cbind(data.frame(c_RNA1), data.frame(c_RNA2)) %>%
  data.frame()

key_genes_pure_tumor <- key_genes_pure_tumor[rownames(key_genes_pure_tumor) %in% genes_I_care_about,]
key_genes_pure_tumor$GENENAME <- rownames(key_genes_pure_tumor)

GR_sig_tum <- prev_method_calculator(GR_Sig, GR_coef, key_genes_pure_tumor)
NIM_sig_tum <- prev_method_calculator(NIM_Sig, NIM_coef, key_genes_pure_tumor)
HOSS_sig_tum <- prev_method_calculator(HOSS_Sig, HOSS_coef, key_genes_pure_tumor)
PRGIT_sig_tum <- prev_method_calculator(PRGIT_Sig, PRGIT_coef, key_genes_pure_tumor)
BHLNSX_sig_tum <- prev_method_calculator(BHLNSX_Sig, BHLNSX_coef, key_genes_pure_tumor)


pure_tumor_prev_methods_results <- rbind(GR_sig_tum, NIM_sig_tum, HOSS_sig_tum, PRGIT_sig_tum, BHLNSX_sig_tum)
rownames(pure_tumor_prev_methods_results) <- c("GR", "NIM", "HOSS", "PRGIT", "BHLNSX")
pure_tumor_prev_methods_results <- pure_tumor_prev_methods_results %>% t()

pure_tumor_prev_methods_results
#save to disk
save(prev_methods_results, prev_methods_results_cpm, pure_tumor_prev_methods_results, file = 'results/prev_methods_results.rdata')



