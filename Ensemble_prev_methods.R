ensemble_results <- read.csv('results/prev_methods_results.csv')

rownames(ensemble_results) <- ensemble_results$X
ensemble_results <- ensemble_results %>%
  dplyr::select(-'X')
#### Categorical Assignment -------
meds <- c()
binary_results <- ensemble_results
for(i in 1:length(ensemble_results[1,])) {
  meds <- c(meds, median(dendrogram_results[,i]))
  binary_results[,i] <- dendrogram_results[,i] > meds[i]
}
binary_results[1:9,]

#sorted in right order
binary_results[c(7:9,1,2,6,3:5),]

get_mode_per_row <- function(x) {
  out <- c()
  for(i in 1:length(x[,1])) {
    out <- c(out, ifelse(sum(x[i,]) > (length(x[1,]) / 2), yes = "HIGH", no = "LOW"))
  }
  return(out)
}

ensemble_binary_results <- get_mode_per_row(binary_results[c(7:9,1,2,6,3:5),])
ensemble_binary_results <- data.frame(sample = rownames(binary_results)[c(7:9,1,2,6,3:5)],
                                      risk_category = ensemble_binary_results)
write.csv(ensemble_binary_results, 'results/ensemble_prev_methods_results.csv')

stripchart(scale(ensemble_results$GR), ylim = c(0,5), xlim = c(-4,4))
stripchart(scale(ensemble_results$HOSS), at = 2, add = T)
stripchart(scale(ensemble_results$PRGIT), at = 3, add = T)
stripchart(scale(ensemble_results$DRCHP), at = 4, add = T)
stripchart(scale(ensemble_results$BHLNSX), at = 5, add = T)
          