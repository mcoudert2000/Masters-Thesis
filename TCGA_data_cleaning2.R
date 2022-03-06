#TCGA data cleaning

#Code provided by Dr. Andy Lynch to pull data from the TCGA database.
require(TCGAbiolinks)
require(edgeR)

query.exp <- GDCquery(
  
  project = "TCGA-GBM",
  
  legacy = TRUE,
  
  data.category = "Gene expression",
  
  data.type = "Gene expression quantification",
  
  platform = "Illumina HiSeq",
  
  file.type = "results",
  
  experimental.strategy = "RNA-Seq",
  
  sample.type = "Primary Tumor"
  
)

GDCdownload(query.exp)
gbm.exp <- GDCprepare(
  
  query = query.exp,
  
  save = TRUE,
  
  save.filename = "gbmExp.rda"
  
)


TCGA_samples <- data.frame(patient = gbm.exp$patient,
                           age_index = gbm.exp$age_at_index,
                           gender = gbm.exp$gender,
                           IDH = gbm.exp$paper_IDH.status,
                           subtype = gbm.exp@colData@listData[["paper_Transcriptome.Subtype"]])
TCGA_counts <- gbm.exp@assays@data@listData[["raw_count"]]

TCGA_counts <- TCGA_counts[, !is.na(TCGA_samples$IDH)]
TCGA_samples <- TCGA_samples[!is.na(TCGA_samples$IDH),]
  
rownames(TCGA_counts) <- gbm.exp@rowRanges@elementMetadata@listData[["gene_id"]]
colnames(TCGA_counts) <- TCGA_samples$patient

rownames(TCGA_counts) <- gsub('\\..*', '', rownames(TCGA_counts))

TCGA_counts <- TCGA_counts[,!duplicated(colnames(TCGA_counts))]
TCGA_samples <- TCGA_samples[!duplicated(TCGA_samples$patient),]

TCGA_counts <- TCGA_counts[!duplicated(rownames(TCGA_counts)),]

dim(TCGA_counts) #Check to see if dimensions make sense
dim(TCGA_samples)

TCGA_data <- DGEList(counts = TCGA_counts, samples = TCGA_samples,
                     group = TCGA_samples$IDH)


save(TCGA_data, file = 'data/TCGA_data_full.rdata')
