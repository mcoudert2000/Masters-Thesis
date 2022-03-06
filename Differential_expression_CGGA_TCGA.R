#Differential Expression for TCGA and CGGA
require(edgeR)
require(limma)

#Following the workflow from this paper: https://www.bioconductor.org/packages/devel/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html#references

load('data/TCGA_data_full.rdata')
load('data/CGGA/CGGA_data.rdata')


TCGA_design <- model.matrix(~0+TCGA_data$samples$group)
CGGA_design <- model.matrix(~0+CGGA_data$samples$group)


colnames(TCGA_design) <- c("IDH_Mutant", "IDH_WT")
colnames(CGGA_design) <- c("IDH_Mutant", "IDH_WT")

TCGA_contrasts <- makeContrasts(MutantvsIDH = IDH_Mutant - IDH_WT, levels = colnames(TCGA_design))
CGGA_contrasts <- makeContrasts(MutantvsIDH = IDH_Mutant - IDH_WT, levels = colnames(CGGA_design))

TCGA_v <- voom(TCGA_data, TCGA_design) #This transforms the data to log2 cpm
CGGA_v <- voom(CGGA_data, CGGA_design)

TCGA_vfit <- lmFit(TCGA_v, TCGA_design)
TCGA_vfit <- contrasts.fit(TCGA_vfit, contrasts=TCGA_contrasts)
TCGA_efit <- eBayes(TCGA_vfit)


CGGA_vfit <- lmFit(CGGA_v, CGGA_design)
CGGA_vfit <- contrasts.fit(CGGA_vfit, contrasts=CGGA_contrasts)
CGGA_efit <- eBayes(CGGA_vfit)

par(mfrow = c(1,2))
plotSA(TCGA_efit, main="TCGA_Mean-Variance Trend")
plotSA(CGGA_efit, main="CGGA_Mean-Variance Trend")


TCGA_tfit <- treat(TCGA_vfit, lfc=1)
TCGA_dt <- decideTests(TCGA_tfit)
summary(TCGA_dt)

TCGA_common <- which(TCGA_dt[,1]!=0)
head(rownames(TCGA_tfit)[TCGA_common], n = 50)

CGGA_tfit <- treat(CGGA_vfit, lfc=1)
CGGA_dt <- decideTests(CGGA_tfit)
summary(CGGA_dt)

CGGA_common <- which(CGGA_dt[,1] != 0)
rownames(CGGA_tfit)[CGGA_common]

diff_expressed_genes <- list(TCGA = rownames(TCGA_tfit)[TCGA_common], CGGA = rownames(CGGA_tfit)[CGGA_common])

save(diff_expressed_genes, file = 'data/TCGA_CGGA_diff_expressed_genes.rdata')

rownames(CGGA_tfit)[CGGA_common] %in% rownames(TCGA_tfit)[TCGA_common]

diff_expressed_genes

