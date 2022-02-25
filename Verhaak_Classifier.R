#Verhaak_classifier

m = read.table("data/verhaak/unifiedScaled.txt", header = TRUE, row.names = 1, check.names = FALSE)
m = as.matrix(m)

#This is the 'answer key' using other method to classify them (see original paper Verhaak et al.)
subtype = read.table("data/verhaak/TCGA_unified_CORE_ClaNC840.txt", sep = "\t", header = TRUE, 
                     check.names = FALSE, stringsAsFactors = FALSE)
subtype = structure(unlist(subtype[1, -(1:2)]), names = colnames(subtype)[-(1:2)])

m = m[, names(subtype)]


#Attempt to cluster on centroids provided by study
centroids <- dplyr::select(read.csv("data/verhaak/ClaNC840_centroids.csv"), -Composite_chr_coords)

verhaak_centroids <- centroids %>%
  as.data.frame() %>%
  dplyr::select(-Gene.Symbol) %>%
  t()

colnames(verhaak_centroids) <- gsub(pattern = '-', replacement = '.', centroids$Gene.Symbol)

#Reading in Astrid's data
astrid_data <- read.csv('data/Astrid_normalized_counts.csv') 

astrid_data <- astrid_data %>%
  dplyr::filter(GENENAME %in% centroids$Gene.Symbol)

astrid_data <- merge_duplicate_rows(astrid_data) %>%
  dplyr::select(-X) %>% t() %>% scale()


colnames(astrid_data) <- gsub("-", '.', colnames(astrid_data))

#Reading in data from verhaak et al
verhaak_data <- m %>%
  as.data.frame() %>%
  dplyr::filter(rownames(m) %in% centroids$Gene.Symbol) %>%
  t()

colnames(verhaak_data) <- gsub('-', '.', colnames(verhaak_data)) #Naming convention was slightly different


categories <- rep('NC', length(verhaak_data[,1]))
for(i in 1:length(verhaak_data[,1])) { #patient loop
  pro_dist <- 0
  class_dist <- 0
  mech_dist <- 0
  neur_dist <- 0
  for(g in colnames(verhaak_data)) {
    pro_dist = pro_dist + abs(verhaak_centroids['Proneural', g] - verhaak_data[i, g])
    class_dist = class_dist + abs(verhaak_centroids['Classical', g] - verhaak_data[i, g])
    mech_dist = mech_dist + abs(verhaak_centroids['Mesenchymal', g] - verhaak_data[i, g])
    neur_dist = neur_dist + abs(verhaak_centroids['Neural', g] - verhaak_data[i, g])
  }
  min_dist <- min(c(pro_dist, class_dist, mech_dist, neur_dist))
  if(pro_dist == min_dist) {
    categories[i] <- "Proneural"
  }
  if(class_dist == min_dist) {
    categories[i] <- "Classical"
  }
  if(mech_dist == min_dist) {
    categories[i] <- "Mesenchymal"
  }
  if(neur_dist == min_dist) {
    categories[i] <- "Neural"
  }
}
  
categories <- as.data.frame(categories)
  
rownames(categories) <- rownames(verhaak_data)

categories$patient <- rownames(categories)
subtype <- data.frame(subtype)
subtype$patient <- rownames(subtype)

verhaak_results <- dplyr::inner_join(categories, subtype, by = 'patient') %>%
  dplyr::select(c('patient', 'categories', 'subtype'))

#accuracy for prediction
sum(verhaak_results$categories == verhaak_results$subtype)/length(verhaak_results$categories)


### Applying same classifier to Astrid's data

astrid_categories <- rep('NC', length(astrid_data[,1]))
for(i in 1:length(astrid_data[,1])) { #patient loop
  pro_dist <- 0
  class_dist <- 0
  mech_dist <- 0
  neur_dist <- 0
  for(g in colnames(astrid_data)) {
    if(length(verhaak_centroids['Proneural', g]) == 0) {
      print(g)
    }
    pro_dist = pro_dist + abs(verhaak_centroids['Proneural', g] - astrid_data[i, g])
    class_dist = class_dist + abs(verhaak_centroids['Classical', g] - astrid_data[i, g])
    mech_dist = mech_dist + abs(verhaak_centroids['Mesenchymal', g] - astrid_data[i, g])
    neur_dist = neur_dist + abs(verhaak_centroids['Neural', g] - astrid_data[i, g])
    # print(pro_dist)
  }
  min_dist <- min(c(pro_dist, class_dist, mech_dist, neur_dist))
  if(pro_dist == min_dist) {
    astrid_categories[i] <- "Proneural"
  }
  if(class_dist == min_dist) {
    astrid_categories[i] <- "Classical"
  }
  if(mech_dist == min_dist) {
    astrid_categories[i] <- "Mesenchymal"
  }
  if(neur_dist == min_dist) {
    astrid_categories[i] <- "Neural"
  }
}

astrid_verhaak_results <- data.frame(sample = rownames(astrid_data), category = astrid_categories)

astrid_verhaak_results

