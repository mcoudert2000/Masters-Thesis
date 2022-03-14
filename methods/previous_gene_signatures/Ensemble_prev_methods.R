load('results/prev_methods_results.rdata')
require(dplyr)
require(ggplot2)
ensemble_results <- prev_methods_results

#### Categorical Assignment -------
meds <- c()
binary_results <- ensemble_results
binary_results_printable <- ensemble_results
for(i in 1:length(ensemble_results[1,])) {
  meds <- c(meds, median(ensemble_results[,i]))
  binary_results[,i] <- ensemble_results[,i] > meds[i]
  binary_results_printable[,i] <- ifelse(ensemble_results[,i] > meds[i], 'HIGH', 'LOW')
}
binary_results[1:9,]

#sorted in right order
binary_results_printable[c(7:9,1,2,6,3:5),]

#Table of outputs to see what percentage of them agree
table(as.factor(rowSums(binary_results)))

get_mode_per_row <- function(x) {
  out <- c()
  for(i in 1:length(x[,1])) {
    if (sum(x[i,]) == (length(x[1,]) / 2)) {
      out <- c(out, "INCONCLUSIVE")
    } else {
      out <- c(out, ifelse(sum(x[i,]) > (length(x[1,]) / 2), yes = "HIGH", no = "LOW"))
    }
  }
  return(out)
}

ensemble_binary_results <- get_mode_per_row(binary_results[c(7:9,1,2,6,3:5),])
ensemble_binary_results <- data.frame(sample = rownames(binary_results)[c(7:9,1,2,6,3:5)],
                                      risk_category = ensemble_binary_results)
ensemble_binary_results

write.csv(ensemble_binary_results, 'results/ensemble_prev_methods_results.csv')

#### Validating performance of ensemble results. 
print(load('data/CGGA/CGGA_data.RDATA'))


ensemble_validation <- CGGA_data$samples
ensemble_validation <- ensemble_validation[!is.na(ensemble_validation$OS), ]
#12 samples removed
ensemble_validation
ensemble_validation_data <- data.frame(RS = get_mode_per_row(binary_results[rownames(ensemble_results) %in% rownames(ensemble_validation),]))
  
ensemble_validation_data$sample <- rownames(ensemble_validation)
ensemble_validation_data$OS <- ensemble_validation$OS
ensemble_validation_data$event <- ensemble_validation$Censor..alive.0..dead.1.

require(survival)
require(survminer)

ens_surv_obj <- Surv(time = ensemble_validation_data$OS, event = ensemble_validation_data$event)
ens_fit <- survfit(ens_surv_obj ~ RS, data = ensemble_validation_data)
ggsurvplot(ens_fit, pval = T)

ggplot(ensemble_validation_data) +
  geom_boxplot(aes(x = OS, col = RS)) +
  ggtitle("Overall Survival Stratified by Ensemble Method Risk Score")


t.test(ensemble_validation_data[ensemble_validation_data$RS == "HIGH",]$OS,
       ensemble_validation_data[ensemble_validation_data$RS == "LOW",]$OS)

validation_samples <- rownames(ensemble_validation)
GR_bin <- binary_results[rownames(binary_results) %in% validation_samples,]$GR
NIM_bin <- binary_results[rownames(binary_results) %in% validation_samples,]$NIM
HOSS_bin <- binary_results[rownames(binary_results) %in% validation_samples,]$HOSS
PRGIT_bin <- binary_results[rownames(binary_results) %in% validation_samples,]$PRGIT
BHLNSX_bin <- binary_results[rownames(binary_results) %in% validation_samples,]$BHLNSX

validate_method <- function(binary_results, 
                            val_samples = validation_samples, 
                            OS = ensemble_validation_data$OS, 
                            event = ensemble_validation_data$event) {
  validate_data <- data.frame(samples = val_samples, 
                              RS = ifelse(binary_results, "HIGH", "LOW"),
                              OS = OS, 
                              event = event)
  ens_fit <- survfit(Surv(validate_data$OS, validate_data$event) ~ RS, data = validate_data)
  p <- ggsurvplot(ens_fit, pval = T)
  return(p)
}


ggsurvplot(survfit(Surv(ensemble_validation_data$OS,
                        ensemble_validation_data$event) ~
                     ifelse(GR_bin, "HIGH", "LOW"), data = ensemble_validation_data), pval = T)

ggsurvplot(survfit(Surv(ensemble_validation_data$OS,
                        ensemble_validation_data$event) ~
                     ifelse(NIM_bin, "HIGH", "LOW"), data = ensemble_validation_data), pval = T)

ggsurvplot(survfit(Surv(ensemble_validation_data$OS,
                        ensemble_validation_data$event) ~
                     ifelse(HOSS_bin, "HIGH", "LOW"), data = ensemble_validation_data), pval = T)

ggsurvplot(survfit(Surv(ensemble_validation_data$OS,
                        ensemble_validation_data$event) ~
                     ifelse(PRGIT_bin, "HIGH", "LOW"), data = ensemble_validation_data), pval = T)

ggsurvplot(survfit(Surv(ensemble_validation_data$OS,
                        ensemble_validation_data$event) ~
                     ifelse(BHLNSX_bin, "HIGH", "LOW"), data = ensemble_validation_data), pval = T)
#Checking if keeping all calculations within CGGA fixes the problem

CGGA_res <- prev_methods_results[grep('CGGA', rownames(prev_methods_results)),]
CGGA_res <- CGGA_res[!is.na(CGGA_data$samples$OS),]
dim(CGGA_res)
CGGA_binary_res <- CGGA_res
CGGA_meds <- c()

for(i in 1:length(CGGA_res[1,])) {
  CGGA_meds <- c(CGGA_meds, median(CGGA_res[,i]))
  CGGA_binary_res[,i] <- CGGA_res[,i] > CGGA_meds[i]
}

CGGA_bin_GR <- CGGA_binary_res$GR
CGGA_bin_NIM <- CGGA_binary_res$NIM
CGGA_bin_HOSS <- CGGA_binary_res$HOSS
CGGA_bin_PRGIT <- CGGA_binary_res$PRGIT
CGGA_bin_BHLNSX <- CGGA_binary_res$BHLNSX

ggsurvplot(survfit(Surv(ensemble_validation_data$OS,
                        ensemble_validation_data$event) ~
                     ifelse(CGGA_bin_GR, "HIGH", "LOW"), data = ensemble_validation_data), pval = T)


ggsurvplot(survfit(Surv(ensemble_validation_data$OS,
                        ensemble_validation_data$event) ~
                     ifelse(CGGA_bin_NIM, "HIGH", "LOW"), data = ensemble_validation_data), pval = T)


ggsurvplot(survfit(Surv(ensemble_validation_data$OS,
                        ensemble_validation_data$event) ~
                     ifelse(CGGA_bin_HOSS, "HIGH", "LOW"), data = ensemble_validation_data), pval = T)


ggsurvplot(survfit(Surv(ensemble_validation_data$OS,
                        ensemble_validation_data$event) ~
                     ifelse(CGGA_bin_PRGIT, "HIGH", "LOW"), data = ensemble_validation_data), pval = T)


ggsurvplot(survfit(Surv(ensemble_validation_data$OS,
                        ensemble_validation_data$event) ~
                     ifelse(CGGA_bin_BHLNSX, "HIGH", "LOW"), data = ensemble_validation_data), pval = T)


