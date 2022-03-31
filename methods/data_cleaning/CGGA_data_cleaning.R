#CGGA data processing
require(dplyr)

#Read files downloaded from CGGA website
CGGA_clinical <- read.delim('data/CGGA/CGGA.mRNAseq_693_clinical.20200506.txt')
CGGA_expression <- read.delim('data/CGGA/CGGA.mRNAseq_693.RSEM-genes.20200506.txt')

#remove samples that include lower grade gliomas
CGGA_clinical <- dplyr::filter(CGGA_clinical, Grade == 'WHO IV')

table(CGGA_clinical$IDH_mutation_status)

#remove version number from gene name
CGGA_expression$Gene_Name <- gsub('\\..*','',CGGA_expression$Gene_Name)

genes_CGGA <- CGGA_expression$Gene_Name

#filter expression to exclude lower grade gliomas
CGGA_expression <- CGGA_expression[,colnames(CGGA_expression) %in% CGGA_clinical$CGGA_ID] 

#remove all but first duplicate of gene
CGGA_expression <- CGGA_expression[!duplicated(genes_CGGA),]

rownames(CGGA_expression) <- unique(genes_CGGA) #add rownames of genes to file

#Remove samples with no IDH status
CGGA_expression <- CGGA_expression[, !is.na(CGGA_clinical$IDH_mutation_status)]
CGGA_clinical <- CGGA_clinical[!is.na(CGGA_clinical$IDH_mutation_status),]

#save to DGEList and save to file
CGGA_data <- DGEList(counts = CGGA_expression, samples = CGGA_clinical,
                     group = CGGA_clinical$IDH_mutation_status)

save(CGGA_data, file = 'data/CGGA/CGGA_data.RDATA')


