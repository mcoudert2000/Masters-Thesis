#CGGA data processing
require(dplyr)

CGGA_clinical <- read.delim('data/CGGA/CGGA.mRNAseq_693_clinical.20200506.txt')
CGGA_expression <- read.delim('data/CGGA/CGGA.mRNAseq_693.RSEM-genes.20200506.txt')


CGGA_clinical <- dplyr::filter(CGGA_clinical, Grade == 'WHO IV')
dim(CGGA_clinical)
table(CGGA_clinical$IDH_mutation_status)

CGGA_expression$Gene_Name <- gsub('\\..*','',CGGA_expression$Gene_Name)

genes_CGGA <- CGGA_expression$Gene_Name

CGGA_expression <- CGGA_expression[,colnames(CGGA_expression) %in% CGGA_clinical$CGGA_ID] 

CGGA_expression$GENENAME <- genes_CGGA
CGGA_expression <- merge_duplicate_rows(CGGA_expression)

dim(CGGA_expression)
dim(CGGA_clinical)

#Remove NA values from IDH_mutation status

CGGA_expression <- CGGA_expression[, !is.na(CGGA_clinical$IDH_mutation_status)]
CGGA_clinical <- CGGA_clinical[!is.na(CGGA_clinical$IDH_mutation_status),]

CGGA_data <- DGEList(counts = CGGA_expression, samples = CGGA_clinical,
                     group = CGGA_clinical$IDH_mutation_status)

save(CGGA_expression, CGGA_clinical, file = 'data/CGGA/CGGA_data.RDATA')


