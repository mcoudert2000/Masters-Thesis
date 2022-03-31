require(immunedeconv)
require(ggplot2)
require(estimate)

#Load full dataset
load('data/combined_data.rdata')

deconvolution_data <- as.matrix(combined_data_normalized$E)

#deconvolute using xcell
deconv_res_xcell <- deconvolute(deconvolution_data, method = 'xcell')

#extract immune and stroma scores
immune_stroma_scores <- as.data.frame(deconv_res_xcell[37:38,])

rownames(immune_stroma_scores) <- immune_stroma_scores$cell_type
immune_stroma_scores <- dplyr::select(immune_stroma_scores, -cell_type)

#categorize samples as high or low immune/stroma scores
meds <- c(median(as.numeric(immune_stroma_scores[1,])), median(as.numeric(immune_stroma_scores[2,])))
binary_results <- immune_stroma_scores
for(i in 1:2) {
  binary_results[i,] <- immune_stroma_scores[i,] > meds[i]
}
binary_results[1:9]

### ESTIMATE
estimate_data <- combined_data_normalized$E

#THIS IS THE WORKFLOW FROM ESTIMATE VIGNETTE
write.table(estimate_data, file = 'data/estimate_normalized_counts.txt', sep = '\t', quote = F)
filterCommonGenes(input.f = 'data/estimate_normalized_counts.txt', 
                  output.f = 'results/estimate/GBM_10412genes.gct',
                  id = 'GeneSymbol')
estimateScore("results/estimate/GBM_10412genes.gct", "results/estimate/GBM_estimate_score.gct", platform="illumina") #check to make sure this is the correct platform
estimate_scores <- read.delim(file="results/estimate/GBM_estimate_score.gct", skip=2, quote = " ")[,1:400]

#format them nicely for later analysis
estimate_scores <- t(estimate_scores)
colnames(estimate_scores) <- estimate_scores[1,]
estimate_scores <- as.data.frame(estimate_scores[3:400,])
estimate_scores[] <- lapply(estimate_scores, as.numeric)

estimate_scores$tumor_purity <- cos(estimate_scores$ESTIMATEScore *0.0001467884+0.6049872018) #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3826632/

rownames(estimate_scores) <- gsub('\\.', '-', rownames(estimate_scores))
save(estimate_scores, file = 'results/estimate_scores.rdata')


#Checking if there is a difference in tumor between TCGA and CGGA samples.
t.test(estimate_scores[grep("TCGA", rownames(estimate_scores)),]$tumor_purity, estimate_scores[grep("CGGA", rownames(estimate_scores)),]$tumor_purity)

#seeing risk categories for Astrid Samples
data.frame(sample = rownames(estimate_scores)[1:9],
           stromal = ifelse(estimate_scores$StromalScore > median(estimate_scores$StromalScore), "HIGH", "LOW")[1:9],
           immune = ifelse(estimate_scores$StromalScore > median(estimate_scores$ImmuneScore), "HIGH", "LOW")[1:9])[c(7:9,1,2,6,3:5),]

#Estimate vs Astrid Tumor Purity
ggplot(data.frame(sample = rownames(estimate_scores)[1:9],
           tumor_purity_est = round(estimate_scores$tumor_purity[1:9],2),
           tumor_purity_astrid = round(Astrid_data$samples$tumor_content,2))[c(7:9,1,2,6,3:5),]) +
  geom_point(aes(x = tumor_purity_est, y = tumor_purity_astrid)) +
  geom_text(aes(x = tumor_purity_est, y = tumor_purity_astrid), label = rownames(estimate_scores)[c(7:9,1,2,6,3:5)], hjust = 0.1, vjust = 0.1) +
  geom_abline(slope = 1, intercept = 0) +
  xlab('Estimate Tumor Purity') + ylab('Astrid Tumor Purity') +
  ggtitle('Estimate vs Astrid Tumor Purity') +
  xlim(c(0.2, 1)) + ylim(c(0.2,1))

#Survival analysis of Stromal and Immune Scores
load('data/CGGA/CGGA_data.RDATA')

deconv_validation <- CGGA_data$samples
deconv_validation <- ensemble_validation[!is.na(ensemble_validation$OS), ]
#12 samples removed
#estimate scores for CGGA samples
CGGA_est_scores <- estimate_scores[grep('CGGA',rownames(estimate_scores)),]

CGGA_est_scores <- CGGA_est_scores[!is.na(CGGA_data$samples$OS),]
strom_bin <- ifelse(CGGA_est_scores$StromalScore > median(CGGA_est_scores$StromalScore), 'HIGH', "LOW")
immune_bin <- ifelse(CGGA_est_scores$ImmuneScore > median(CGGA_est_scores$ImmuneScore), 'HIGH', "LOW")
deconv_validation_data <- data.frame(stromalCat = strom_bin, immuneCat = immune_bin)
deconv_validation_data$sample <- rownames(deconv_validation)
deconv_validation_data$OS <- deconv_validation$OS
deconv_validation_data$event <- deconv_validation$Censor..alive.0..dead.1.
deconv_validation_data$estCat <- ifelse(CGGA_est_scores$ESTIMATEScore > median(CGGA_est_scores$ESTIMATEScore), 'HIGH', "LOW")

require(survival)
require(survminer)

#Validating on Stromal score
deconv_surv_obj <- Surv(time = deconv_validation_data$OS, event = deconv_validation_data$event)
deconv_fit <- survfit(deconv_surv_obj ~ stromalCat, data = deconv_validation_data)
stromal_plot <-ggsurvplot(deconv_fit, pval = T)

#Validating on Immune Score
deconv_surv_obj <- Surv(time = deconv_validation_data$OS, event = deconv_validation_data$event)
deconv_fit <- survfit(deconv_surv_obj ~ immuneCat, data = deconv_validation_data)
immune_plot <- ggsurvplot(deconv_fit, pval = T)

#Validating on Estimate Score
deconv_surv_obj <- Surv(time = deconv_validation_data$OS, event = deconv_validation_data$event)
deconv_fit <- survfit(deconv_surv_obj ~ estCat, data = deconv_validation_data)
estimate_plot <- ggsurvplot(deconv_fit, pval = T)

require(ggpubr)
#arrange the plots into 1
ggarrange(plotlist = list(stromal_plot$plot, immune_plot$plot, estimate_plot$plot),
          labels = c("Stromal", "Immune", "Estimate"), ncol = 3)

#Validating in TCGA data
load('data/TCGA_data_survival.rdata')

TCGA_est_scores <- estimate_scores[grep('TCGA',rownames(estimate_scores)),]

TCGA_est_scores <- TCGA_est_scores[gsub('\\.', '-',rownames(TCGA_est_scores)) %in% colnames(TCGA_survival_data),]
TCGA_stromCat <- ifelse(TCGA_est_scores$StromalScore > median(TCGA_est_scores$StromalScore), "HIGH", "LOW")
TCGA_immuneCat <- ifelse(TCGA_est_scores$ImmuneScore > median(TCGA_est_scores$ImmuneScore), "HIGH", "LOW")
TCGA_estCat <- ifelse(TCGA_est_scores$ESTIMATEScore > median(TCGA_est_scores$ESTIMATEScore), "HIGH", "LOW")

TCGA_time <- TCGA_survival_data$samples$time
TCGA_event <- ifelse(TCGA_survival_data$samples$event == "Dead", 1, 0)

TCGA_strom_p <- ggsurvplot(survfit(Surv(TCGA_time, TCGA_event) ~ TCGA_stromCat), data = TCGA_survival_samples, pval = T)
TCGA_imm_p <- ggsurvplot(survfit(Surv(TCGA_time, TCGA_event) ~ TCGA_immuneCat), data = TCGA_survival_samples, pval = T)
TCGA_est_p <- ggsurvplot(survfit(Surv(TCGA_time, TCGA_event) ~ TCGA_estCat), data = TCGA_survival_samples, pval = T)

ggarrange(plotlist = list(TCGA_strom_p$plot, TCGA_imm_p$plot, TCGA_est_p$plot),
          labels = c("Stromal", "Immune", "EstimateScore"), ncol = 3)



#Plots

#calculate correlation between immune and stromal score
cor(estimate_scores$ImmuneScore, estimate_scores$StromalScore)

#plot of immune vs stromal score
ggplot(estimate_scores) +
  geom_point(data = estimate_scores[10:248,], aes(x = StromalScore, y = ImmuneScore, col = 'CGGA')) +
  geom_point(data = estimate_scores[1:9,], aes(x = StromalScore, y = ImmuneScore, col = 'ASTRID\'s')) +
  geom_point(data = estimate_scores[249:398,], aes(x = StromalScore, y = ImmuneScore, col = 'TCGA')) +
  geom_smooth(data = estimate_scores, aes(x = StromalScore, y = ImmuneScore), method = lm, se = T, linetype = 2, col = 'black') +
  ggtitle("Stromal Score vs Immune Score") +
  labs(col = "Source")

ggplot(estimate_scores[10:248,]) +
  geom_point(aes(x = -1, y = ESTIMATEScore, col = 'CGGA'), alpha =0.3) +
  geom_point(data = estimate_scores[1:9,], alpha = 0.3, aes(x = 0, y = ESTIMATEScore, col = 'ASTRID\'s')) +
  geom_point(data = estimate_scores[249:398,],alpha =0.3, aes(x = 1, y = ESTIMATEScore, col = 'TCGA')) +
  geom_hline(yintercept = median(estimate_scores$ESTIMATEScore), linetype = 2)


ggplot(estimate_scores) +
  geom_boxplot(aes(x = tumor_purity)) +
  coord_flip() 


### Linear Equation Solver

#Idea: RNA_tot = RNA_tumor * tumor_content + RNA_non_tumor * (1 - tumor_content)
load('data/Astrid_data_counts.rdata')

R2A <- combined_data_normalized$E[,'R2A']
R2C <- combined_data_normalized$E[,'R2C']

R2A_tp <- estimate_scores[c(3),]$tumor_purity
R2C_tp <- estimate_scores[c(5),]$tumor_purity

R1B <- combined_data_normalized$E[,'R1B']
R1C <- combined_data_normalized$E[,'R1C']


R1C_tp <- estimate_scores[6,]$tumor_purity
R1B_tp <- estimate_scores[2,]$tumor_purity

#check for stability
R2A_tp - R2C_tp #very small difference in tumor content
R1C_tp - R1B_tp
c_RNA2 <- ( R2A - R2C + R2A_tp * R2C - R2C_tp * R2A ) / (R2A_tp - R2C_tp) #solve the linear equation
c_RNA1 <- ( R1B - R1C + R1B_tp * R1C - R1C_tp * R1B ) / (R1B_tp - R1C_tp)
#This is (theoretically) isolated cancer RNA

#write to file
save(c_RNA1, c_RNA2, file = 'results/estimate/isolated_tumor.rdata')

isolated_tumor <- data.frame(c_RNA1 = c_RNA1, c_RNA2 = c_RNA2)


#calculate isolated tumor, ESTIMATE tumor content 
filterCommonGenes(input.f = 'results/estimate/isolated_tumor.txt', 
                  output.f = 'results/estimate/isolated_tumor.gct',
                  id = 'GeneSymbol')
estimateScore("results/estimate/isolated_tumor.gct", "results/estimate/isolated_tumor_estimate_score.gct", platform="illumina") #check to make sure this is the correct platform
isolated_estimate_scores <- read.delim(file="results/estimate/isolated_tumor_estimate_score.gct", skip=2, quote = " ")[,2:4]

rownames(isolated_estimate_scores) <- isolated_estimate_scores$Description
isolated_estimate_scores <- isolated_estimate_scores %>%
  dplyr::select(-'Description') %>% t() %>% data.frame()



isolated_estimate_scores$tumor_purity <- cos(isolated_estimate_scores$ESTIMATEScore *0.0001467884+0.6049872018)
isolated_estimate_scores
