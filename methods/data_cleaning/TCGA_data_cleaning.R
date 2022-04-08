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

#save certain qualities
TCGA_samples <- data.frame(patient = gbm.exp$patient,
                           age_index = gbm.exp$age_at_index,
                           gender = gbm.exp$gender,
                           IDH = gbm.exp$paper_IDH.status,
                           subtype = gbm.exp@colData@listData[["paper_Transcriptome.Subtype"]],
                           time = gbm.exp$days_to_last_follow_up,
                           event = gbm.exp$vital_status)
TCGA_counts <- gbm.exp@assays@data@listData[["raw_count"]]

#Remove samples with no IDH-status recorded
TCGA_counts <- TCGA_counts[, !is.na(TCGA_samples$IDH)]
TCGA_samples <- TCGA_samples[!is.na(TCGA_samples$IDH),]
  
rownames(TCGA_counts) <- gbm.exp@rowRanges@elementMetadata@listData[["gene_id"]]
colnames(TCGA_counts) <- TCGA_samples$patient

#remove version number from gene names
rownames(TCGA_counts) <- gsub('\\..*', '', rownames(TCGA_counts))

#remove duplicate samples from the same patient
TCGA_counts <- TCGA_counts[,!duplicated(colnames(TCGA_counts))]
TCGA_samples <- TCGA_samples[!duplicated(TCGA_samples$patient),]

#remove duplicate genenames after the first
TCGA_counts <- TCGA_counts[!duplicated(rownames(TCGA_counts)),]

TCGA_data <- DGEList(counts = TCGA_counts, samples = TCGA_samples,
                     group = TCGA_samples$IDH)

###Get survival data

#Remove samples with no survival data
TCGA_survival_samples <- TCGA_samples[!is.na(TCGA_samples$time) & TCGA_samples$event != "Not Reported",]
TCGA_survival_counts <- TCGA_counts[,!is.na(TCGA_samples$time) & TCGA_samples$event != "Not Reported"]

TCGA_survival_data <- DGEList(counts = TCGA_survival_counts, samples = TCGA_survival_samples)

TCGA_survival_data$samples$time

#write both to disk
save(TCGA_data, file = 'data/TCGA_data_full.rdata')
save(TCGA_survival_data, file = 'data/TCGA_data_survival.rdata')
