require(dendextend)
require(graphics)
require(dplyr)
#Try with these genes as well (table 1): https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5940130/

print(load('results/prev_methods_results.rdata'))
dendrogram_results <- scale(prev_methods_results_cpm)

#Calculates the distance between points
hc <- dendrogram_results %>% dist() %>% hclust()
dhc <- as.dendrogram(hc)

#This function just colors the leaves of our data and removes labels of others.
colLab <- function(n) {
  if(is.leaf(n)) {
    #I take the current attributes
    a = attributes(n)
    #print(a)
    #I deduce the line in the original data, and so the treatment and the specie.
    treatment = a$label
    if (grepl("CGGA", treatment)) {
      col_treatment = "#D55E00"
      lab.cex = 0.01
    } else {
    if (grepl("TCGA", treatment)) {
      col_treatment = "#56B4E9"
      lab.cex = 0.01
    } else{
      col_treatment = "#000000"
      lab.cex = .5
    }
    }
    #Modification of leaf attribute
    attr(n,"nodePar") <- c(a$nodePar,list(lab.cex=lab.cex,pch=10,col=col_treatment,lab.font=1))
  }
  return(n)
}
dL <- dendrapply(dhc, colLab)

png(filename = "plots/Clustering on Results from Various Gene Signature Studies.png",
    width = 1600, height = 1200)
plot(dL, main = "Clustering on Results from Various Gene Signature Studies")

dev.off()


