require(dplyr)
require(edgeR)
#Load in our data
load('data/GeneCounts.RData')
time <- as.factor(c(2, 2, 3, 3, 3, 2, 1, 1, 1))
sample <- c("A002", "A004", "A005", "A007", "A012", "A013", "A014", "A016", "A018")
label <- c('R1A', 'R1B', 'R2A', 'R2B', 'R2C', 'R1C', 'PA', 'PB', 'PC')
tumor_content <- c(0.26, .93, .90, .71, .85, .87, .24, .26, .24)
sample_data <- data.frame(time = time, sample = sample, label = label, tumor_content = tumor_content)
rownames(sample_data) <- sample_data$label
#Swap gene names as before
prep_data <- function(txi, tx2gene) {
  txi_out <- txi
  #Filter out the water only sample
  txi_out[["counts"]] <- txi[['counts']][,1:9]
  
  #Swap geneID with symbol for comparison with online methods
  key <- tx2gene %>% dplyr::select(c('GENEID', 'symbol')) %>%
    distinct()
  
  #Get the current row names
  curr_rows <- rownames(txi_out$counts)
  new_rows <- rep('', length(curr_rows))
  
  for(i in seq_len(length(curr_rows))) { #Match keys
    new_rows[i] = key[which(curr_rows[i] == key), 2]
  }
  
  #MAKE SURE TO USE COUNTS HERE RATHER THAN ABUNDANCE
  rownames(txi_out$counts) <- new_rows
  return(txi_out$counts)
}

#Load and normalize data to CPM
prev_data <- prep_data(txi, tx2gene) 

Astrid_data <- DGEList(counts = prev_data, samples = sample_data)
colnames(Astrid_data) <- Astrid_data$samples$label
save(Astrid_data, file = 'data/Astrid_data_counts.rdata')

normalized_prev_data <- prev_data
for(i in 1:length(prev_data[1,])) {
  normalized_prev_data[, i] <- prev_data[, i] * (1e6 / sum(prev_data[, i]))
}

#Tidy data into the same format as TCGA
prev_data_frame <- normalized_prev_data %>%
  as.data.frame() %>%
  mutate(GENENAME = rownames(normalized_prev_data))  %>%
  dplyr::select(GENENAME, everything()) 

#Label our data for later analysis
colnames(prev_data_frame) <- c("GENENAME", sample_data$label) 

prev_data_frame
write.csv(prev_data_frame, file = 'data/Astrid_normalized_counts.csv')
key_genes_prev_data_frame <- prev_data_frame %>%
  dplyr::filter(GENENAME %in% genes_I_care_about)