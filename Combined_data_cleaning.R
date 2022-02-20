GR_Sig <- c('GRIA2', 'RYR3')
HOSS_Sig <- c('HOXC10', 'OSMR', 'SCARA3', 'SLC39A10')
LMSZ_Sig <- c('LHX2', 'MEOX2', 'SNAI2', 'ZNF22')
PRGIT_Sig <- c('PTPRN', 'RGS14', 'G6PC3', 'IGFBP2', 'TIMP4')
DRCHP_Sig <- c('DES', 'RANBP17', 'CLEC5A', 'HOXC11', 'POSTN')
BHLNSX_Sig <- c('BPIFB2', 'HOXA13', 'LRRC10', 'NELL1', 'SDR16C5', 'XIRP2')
genes_I_care_about <- c(GR_Sig, HOSS_Sig, LMSZ_Sig, PRGIT_Sig, DRCHP_Sig, BHLNSX_Sig)


TCGA_data <- read.csv("data/TCGA_normalized_counts.csv")
Astrid_data <- read.csv("data/Astrid_normalized_counts.csv")

require(dplyr)
Astrid_data <- Astrid_data %>%
  as.data.frame() %>%
  mutate(GENENAME = rownames(normalized_prev_data))  %>%
  dplyr::select(GENENAME, everything()) 


full_dataset <- inner_join(Astrid_data, TCGA_data, by = "GENENAME") %>%
  dplyr::select(-c("GENEID", "X.x", "X.y"))

key_genes_dataset <- full_dataset %>%
  dplyr::filter(GENENAME %in% genes_I_care_about)

write.csv(key_genes_dataset, 'data/combined_normalized_counts_filtered.csv')
write.csv(full_dataset, 'data/combined_normalized_counts.csv')
