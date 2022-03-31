require(dplyr)
require(edgeR)

#Load in our data
load('data/GeneCounts.RData') #provided by Dr. Andy Lynch
time <- as.factor(c(2, 2, 3, 3, 3, 2, 1, 1, 1))
sample <- c("A002", "A004", "A005", "A007", "A012", "A013", "A014", "A016", "A018")
label <- c('R1A', 'R1B', 'R2A', 'R2B', 'R2C', 'R1C', 'PA', 'PB', 'PC')
tumor_content <- c(0.26, .93, .90, .71, .85, .87, .24, .26, .24)
sample_data <- data.frame(time = time, sample = sample, label = label, tumor_content = tumor_content)
rownames(sample_data) <- sample_data$label


#Swap gene names from GENEID to HUGO symbol and remove water only sample
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

#Load data
prev_data <- prep_data(txi, tx2gene) 

Astrid_data <- DGEList(counts = prev_data, samples = sample_data)
colnames(Astrid_data) <- Astrid_data$samples$label

#Save to an rdata file for use in other analysis
save(Astrid_data, file = 'data/Astrid_data_counts.rdata')

