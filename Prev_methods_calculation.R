require(dplyr)
GR_Sig <- c('GRIA2', 'RYR3')
HOSS_Sig <- c('HOXC10', 'OSMR', 'SCARA3', 'SLC39A10')
LMSZ_Sig <- c('LHX2', 'MEOX2', 'SNAI2', 'ZNF22')
PRGIT_Sig <- c('PTPRN', 'RGS14', 'G6PC3', 'IGFBP2', 'TIMP4')
DRCHP_Sig <- c('DES', 'RANBP17', 'CLEC5A', 'HOXC11', 'POSTN')
BHLNSX_Sig <- c('BPIFB2', 'HOXA13', 'LRRC10', 'NELL1', 'SDR16C5', 'XIRP2')
genes_I_care_about <- c(GR_Sig, HOSS_Sig, LMSZ_Sig, PRGIT_Sig, DRCHP_Sig, BHLNSX_Sig)

#Load in data
key_genes_data <- read.csv(file = 'data/combined_normalized_counts_filtered.csv') 

### Prev Methods -----
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
GR_coef <- list(GRIA2 = 0.1, RYR3 = -0.1)
GR_sig_out <- prev_method_calculator(GR_Sig, GR_coef, key_genes_data)


#### HOSS-Sig-----
HOSS_Sig
HOSS_coef <- list(HOXC10 = 0.089, OSMR = 0.238, SCARA3 = 0.238, SLC39A10 = -0.424) 
#This is from table 1 of paper
HOSS_sig_out <- prev_method_calculator(HOSS_Sig, HOSS_coef, key_genes_data)

#### PRGIT-Sig------
PRGIT_Sig
PRGIT_coef <- list(PTPRN = 0.50894, RGS14 = 0.54671, 
                   G6PC3 = 1.20753, IGFBP2 = 0.25845, TIMP4 = -0.20684)
PRGIT_sig_out <- prev_method_calculator(PRGIT_Sig, PRGIT_coef, key_genes_data)

#### DRCHP-Sig------
DRCHP_Sig
DRCHP_coef <- list(DES = 0.5536, RANBP17 = -0.7340, 
                   CLEC5A = 0.0995, HOXC11 = 0.2810, POSTN = 0.0566)
DRCHP_sig_out <- prev_method_calculator(DRCHP_Sig, DRCHP_coef, key_genes_data)

#### BHLNSX-Sig------
BHLNSX_Sig
BHLNSX_coef <- list(BPIFB2 = -0.0881, HOXA13 = -0.0854, LRRC10 = -0.0614, NELL1 = -0.151, SDR16C5 = -0.0945, XIRP2 = -0.1441)
BHLNSX_sig_out <- prev_method_calculator(DRCHP_Sig, DRCHP_coef, key_genes_data)

prev_methods_results <- data.frame(GR = t(GR_sig_out), HOSS = t(HOSS_sig_out), PRGIT = t(PRGIT_sig_out),
                      DRCHP = t(DRCHP_sig_out), BHLNSX = t(BHLNSX_sig_out))

prev_methods_results

write.csv(prev_methods_results, file = 'results/prev_methods_results.csv')

