#Combining all data
require(edgeR)
require(dplyr)

load('data/Astrid_data_counts.rdata')
load('data/TCGA_data_full.rdata')
print(load('data/CGGA/CGGA_data.RDATA'))


colnames(Astrid_data$samples)
colnames(CGGA_data$samples)
colnames(TCGA_data$samples)

#Columns to keep: lib.size, sample, IDH

A_samples <- Astrid_data$samples %>% dplyr::select(c('label', 'lib.size'))
A_samples$IDH <- rep("Unknown", 9)
colnames(A_samples) <- c("sample", "lib.size", "IDH")

C_samples <- CGGA_data$samples %>% 
  mutate(sample = rownames(CGGA_data$samples)) %>%
  dplyr::select(c("sample",'lib.size',"IDH_mutation_status"))

colnames(C_samples) <- c("sample", "lib.size", "IDH")

T_samples <- TCGA_data$samples %>%
  mutate(sample = rownames(TCGA_data$samples)) %>%
  dplyr::select(c("sample", "lib.size", "IDH"))

A_samples$source <- rep("Astrid",length(A_samples[,1]))
T_samples$source <- rep("TCGA",length(T_samples[,1]))
C_samples$source <- rep("CGGA",length(C_samples[,1]))


combined_data_samples <- rbind(A_samples, C_samples, T_samples)
combined_data_samples <- combined_data_samples %>% mutate(IDH = replace(IDH, IDH == "Wildtype", "WT"))

combined_data_samples$IDH <- factor(combined_data_samples$IDH, levels = c("Unknown", "WT", "Mutant"))

A_counts <- data.frame(t(Astrid_data$counts))
T_counts <- data.frame(t(TCGA_data$counts))
C_counts <- data.frame(t(CGGA_data$counts))

dim(A_counts)
dim(T_counts)
dim(C_counts)

genes_in_all <- intersect(intersect(colnames(A_counts), colnames(T_counts)), colnames(C_counts))

A_counts <- A_counts %>% dplyr::select(all_of(genes_in_all))
T_counts <- T_counts %>% dplyr::select(all_of(genes_in_all))
C_counts <- C_counts %>% dplyr::select(all_of(genes_in_all))

combined_data_counts <- rbind(A_counts, T_counts, C_counts)

combined_data <- DGEList(counts = t(combined_data_counts), samples = combined_data_samples)

combined_data_normalized <- voom(combined_data)
dim(combined_data_normalized)
table(combined_data_normalized$targets$source)

save(combined_data, combined_data_normalized, file = 'data/combined_data.rdata')
