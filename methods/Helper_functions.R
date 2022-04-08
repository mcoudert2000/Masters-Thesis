#HELPER FUNCTION


#AVERAGE DUPLICATE ROWS
#If two rows have the same name then take the mean of the two and replace it with that
#Input should be have columns dat$GENENAME, samp1, samp2,..., sampn

merge_duplicate_rows <- function(dat) {
  library(data.table)
  library(dplyr)
  out <- aggregate(dat,by=list(dat$GENENAME),FUN=mean)
  out <- out %>% mutate(GENENAME = Group.1) %>%
    dplyr::select(-Group.1)
  print(paste(dim(dat)[1] - dim(out)[1], ' rows removed.'))
  rownames(out) <- out$GENENAME
  return(dplyr::select(out, -GENENAME))
}

colors = list(CGGA = "#D55E00", TCGA = "#56B4E9", 
  Wendler = "#000000")

colors_high_low = list(High = "#CC79A7", Low = "#E69F00")


#scales::show_col(c("#D55E00","#56B4E9","#000000"))

                   