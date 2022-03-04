require(dplyr)
require(AnnotationDbi)
require(EnsDb.Hsapiens.v86)
path <- "data/tcga_female_download/"
files <- list.files(path) 

#Load all the files into a dataframe
i = 0
for(f in files) {
  folder <- gsub(" ", "",paste(path,f))
  zipped = list.files(gsub(" ", "",paste(path,f))) 
  
  tryCatch(temp <- read.table((((gsub(" ", "",paste(folder,'/', zipped)))))),
           error = function(e) {
             print("Reading this file had an issue")
             temp = 0})
  if (i == 0) {
    genes <- data.frame(gene = temp[,1],counts = temp[,2])
    i = i + 1
  } else {
    temp <- data.frame(gene = temp[,1], counts = temp[,2])
    genes <- merge(genes, temp[0:60483,], by = 'gene')
    i = i + 1
  }
}
#Rename the columns
colnames(genes) <- c("GENEID",seq(1, 60))

#Remove the version number and decimal
genes$GENEID <- gsub("\\..*", "", genes$GENEID)
rows <- genes$GENEID

#Subset to only get genes that are relevant to the studies evaluated
GR_Sig <- c('GRIA2', 'RYR3')
HOSS_Sig <- c('HOXC10', 'OSMR', 'SCARA3', 'SLC39A10')
LMSZ_Sig <- c('LHX2', 'MEOX2', 'SNAI2', 'ZNF22')
PRGIT_Sig <- c('PTPRN', 'RGS14', 'G6PC3', 'IGFBP2', 'TIMP4')
DRCHP_Sig <- c('DES', 'RANBP17', 'CLEC5A', 'HOXC11', 'POSTN')
BHLNSX_Sig <- c('BPIFB2', 'HOXA13', 'LRRC10', 'NELL1', 'SDR16C5', 'XIRP2')
genes_I_care_about <- c(GR_Sig, HOSS_Sig, LMSZ_Sig, PRGIT_Sig, DRCHP_Sig, BHLNSX_Sig)

ourCols <- c("SYMBOL", "GENEID")
ourKeys <- genes$GENEID

annot <- AnnotationDbi::select(EnsDb.Hsapiens.v86, 
                               keys=ourKeys, 
                               columns=ourCols, 
                               keytype="GENEID") %>%
  dplyr::filter(nchar(GENEID) > 8) #This last line filters out non-ENSG gene names to avoid duplicate symbols

#Normalize to CPM (Need to figure out FPKM later)
genes_normalized <- genes
for(i in 2:length(genes)) {
  genes_normalized[, i] <- genes[, i] * (1e6 / sum(genes[, i]))
}

#Get the GENENAMES from the annotation key
counts_for_TCGA <- genes_normalized %>%
  dplyr::filter(GENEID %in% annot$GENEID) %>%
  mutate(GENENAME = annot$SYMBOL[match(GENEID, annot$GENEID)]) %>%
  dplyr::select(GENENAME, everything())
dim(counts_for_TCGA)

#Write to a file for access later
write.csv(counts_for_TCGA, "data/TCGA_normalized_counts.csv")

