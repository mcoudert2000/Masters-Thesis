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

pure_tum_bin_res <- pure_tumor_prev_methods_results
for(i in 1:length(ensemble_results[1,])) {
  pure_tum_bin_res[,i] <- ifelse(pure_tumor_prev_methods_results[,i] > meds[i], 'HIGH', 'LOW')
}
pure_tum_bin_res

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


{t.test(ensemble_validation_data[ensemble_validation_data$RS == "HIGH",]$OS,
       ensemble_validation_data[ensemble_validation_data$RS == "LOW",]$OS)

validation_samples <- rownames(ensemble_validation)
GR_bin <- binary_results[rownames(binary_results) %in% validation_samples,]$GR
NIM_bin <- binary_results[rownames(binary_results) %in% validation_samples,]$NIM
HOSS_bin <- binary_results[rownames(binary_results) %in% validation_samples,]$HOSS
PRGIT_bin <- binary_results[rownames(binary_results) %in% validation_samples,]$PRGIT
BHLNSX_bin <- binary_results[rownames(binary_results) %in% validation_samples,]$BHLNSX



CGGA_GR_p <- ggsurvplot(survfit(Surv(ensemble_validation_data$OS,
                        ensemble_validation_data$event) ~
                     ifelse(GR_bin, "HIGH", "LOW"), data = ensemble_validation_data), pval = T)

CGGA_NIM_p <- ggsurvplot(survfit(Surv(ensemble_validation_data$OS,
                        ensemble_validation_data$event) ~
                     ifelse(NIM_bin, "HIGH", "LOW"), data = ensemble_validation_data), pval = T)

CGGA_HOSS_p <- ggsurvplot(survfit(Surv(ensemble_validation_data$OS,
                        ensemble_validation_data$event) ~
                     ifelse(HOSS_bin, "HIGH", "LOW"), data = ensemble_validation_data), pval = T)

CGGA_PRGIT_p <- ggsurvplot(survfit(Surv(ensemble_validation_data$OS,
                        ensemble_validation_data$event) ~
                     ifelse(PRGIT_bin, "HIGH", "LOW"), data = ensemble_validation_data), pval = T)

CGGA_BHLNSX_p  <- ggsurvplot(survfit(Surv(ensemble_validation_data$OS,
                        ensemble_validation_data$event) ~
                     ifelse(BHLNSX_bin, "HIGH", "LOW"), data = ensemble_validation_data), pval = T)

CGGA_full_p <- ggarrange(plotlist = list(CGGA_GR_p$plot, CGGA_NIM_p$plot, CGGA_HOSS_p$plot, CGGA_PRGIT_p$plot, CGGA_BHLNSX_p$plot),
                           labels = c("GR", "NIM", "HOSS", "PRGIT", "BHLNSX"))
  }
CGGA_full_p
#Checking if keeping all calculations within CGGA fixes the problem

{CGGA_res <- prev_methods_results[grep('CGGA', rownames(prev_methods_results)),]
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

CGGA_GR_p <- ggsurvplot(survfit(Surv(ensemble_validation_data$OS,
                        ensemble_validation_data$event) ~
                     ifelse(CGGA_bin_GR, "HIGH", "LOW"), data = ensemble_validation_data), pval = T)


CGGA_NIM_p <- ggsurvplot(survfit(Surv(ensemble_validation_data$OS,
                        ensemble_validation_data$event) ~
                     ifelse(CGGA_bin_NIM, "HIGH", "LOW"), data = ensemble_validation_data), pval = T)


CGGA_HOSS_p <-ggsurvplot(survfit(Surv(ensemble_validation_data$OS,
                        ensemble_validation_data$event) ~
                     ifelse(CGGA_bin_HOSS, "HIGH", "LOW"), data = ensemble_validation_data), pval = T)


CGGA_PRGIT_p <-ggsurvplot(survfit(Surv(ensemble_validation_data$OS,
                        ensemble_validation_data$event) ~
                     ifelse(CGGA_bin_PRGIT, "HIGH", "LOW"), data = ensemble_validation_data), pval = T)


CGGA_BHLNSX_p <-ggsurvplot(survfit(Surv(ensemble_validation_data$OS,
                        ensemble_validation_data$event) ~
                     ifelse(CGGA_bin_BHLNSX, "HIGH", "LOW"), data = ensemble_validation_data), pval = T)

CGGA_within_p <- ggarrange(plotlist = list(CGGA_GR_p$plot, CGGA_NIM_p$plot, CGGA_HOSS_p$plot, CGGA_PRGIT_p$plot, CGGA_BHLNSX_p$plot),
          labels = c("GR", "NIM", "HOSS", "PRGIT", "BHLNSX"))
}
CGGA_within_p
#TCGA surv analysis
load('data/TCGA_data_survival.rdata')
{
TCGA_survival_data$samples$patient

TCGA_res <- prev_methods_results[gsub('\\.', '-', rownames(binary_results)) %in% TCGA_survival_data$samples$patient,]

TCGA_binary <- TCGA_res
TCGA_meds <- c()
for(i in 1:length(TCGA_res[1,])) {
  TCGA_meds <- c(TCGA_meds, median(TCGA_res[,i]))
  TCGA_binary[,i] <- TCGA_res[,i] > TCGA_meds[i]
}
colnames(TCGA_binary) <- c("GR", "NIM", "HOSS", "PRGIT", "BHLNSX")

TCGA_samples <- TCGA_survival_data$samples
TCGA_samples$event <- ifelse(TCGA_survival_data$samples$event == "Dead", 1, 0)

TCGA_event <- TCGA_samples$event
TCGA_time <- TCGA_samples$time

TCGA_GR <- ifelse(TCGA_binary$GR, "HIGH", "LOW")
TCGA_NIM <- ifelse(TCGA_binary$NIM, "HIGH", "LOW")
TCGA_HOSS <- ifelse(TCGA_binary$HOSS, "HIGH", "LOW")
TCGA_PRGIT <- ifelse(TCGA_binary$PRGIT, "HIGH", "LOW")
TCGA_BHLNSX <- ifelse(TCGA_binary$BHLNSX, "HIGH", "LOW")

TCGA_GR_p <- ggsurvplot(survfit(Surv(TCGA_time,
                        TCGA_event) ~
                     TCGA_GR, data = TCGA_samples), pval = T)

TCGA_NIM_p <- ggsurvplot(survfit(Surv(TCGA_time,
                        TCGA_event) ~
                     TCGA_NIM, data = TCGA_samples), pval = T)

TCGA_HOSS_p <-ggsurvplot(survfit(Surv(TCGA_time,
                        TCGA_event) ~
                     TCGA_HOSS, data = TCGA_samples), pval = T)

TCGA_PRGIT_p <-ggsurvplot(survfit(Surv(TCGA_time,
                        TCGA_event) ~
                     TCGA_PRGIT, data = TCGA_samples), pval = T)

TCGA_BHLNSX_p <-ggsurvplot(survfit(Surv(TCGA_time,
                        TCGA_event) ~
                     TCGA_BHLNSX, data = TCGA_samples), pval = T)

TCGA_p <- ggarrange(plotlist = list(TCGA_GR_p$plot, TCGA_NIM_p$plot,
                          TCGA_HOSS_p$plot, TCGA_PRGIT_p$plot,
                          TCGA_BHLNSX_p$plot), labels = c('GR', 'NIM', "HOSS", "PRGIT", "BHLNSX"), ncol = 5) +
  ggtitle("Within TCGA Log 2 CPM Survival Analysis")

ggsave('plots/prev_method_results/TCGA_log2cpm.png', plot = TCGA_p, width = 3000, height = 3000/1.846, units = 'px')
}

##REPEATED WITH CPM rather than log2 CPM
{TCGA_res <- prev_methods_results_cpm[gsub('\\.', '-', rownames(binary_results)) %in% TCGA_survival_data$samples$patient,]

TCGA_binary <- TCGA_res
TCGA_meds <- c()
for(i in 1:length(TCGA_res[1,])) {
  TCGA_meds <- c(TCGA_meds, median(TCGA_res[,i]))
  TCGA_binary[,i] <- TCGA_res[,i] > TCGA_meds[i]
}
TCGA_binary

TCGA_samples <- TCGA_survival_data$samples
TCGA_samples$event <- ifelse(TCGA_survival_data$samples$event == "Dead", 1, 0)

TCGA_event <- TCGA_samples$event
TCGA_time <- TCGA_samples$time

TCGA_GR <- ifelse(TCGA_binary$GR, "HIGH", "LOW")
TCGA_NIM <- ifelse(TCGA_binary$NIM, "HIGH", "LOW")
TCGA_HOSS <- ifelse(TCGA_binary$HOSS, "HIGH", "LOW")
TCGA_PRGIT <- ifelse(TCGA_binary$PRGIT, "HIGH", "LOW")
TCGA_BHLNSX <- ifelse(TCGA_binary$BHLNSX, "HIGH", "LOW")

TCGA_GR_cpm_p <- ggsurvplot(survfit(Surv(TCGA_time,
                        TCGA_event) ~
                     TCGA_GR, data = TCGA_samples), pval = T)

TCGA_NIM_cpm_p <- ggsurvplot(survfit(Surv(TCGA_time,
                        TCGA_event) ~
                     TCGA_NIM, data = TCGA_samples), pval = T)

TCGA_HOSS_cpm_p <- ggsurvplot(survfit(Surv(TCGA_time,
                        TCGA_event) ~
                     TCGA_HOSS, data = TCGA_samples), pval = T)

TCGA_PRGIT_cpm_p <- ggsurvplot(survfit(Surv(TCGA_time,
                        TCGA_event) ~
                     TCGA_PRGIT, data = TCGA_samples), pval = T)

TCGA_BHLNSX_cpm_p <- ggsurvplot(survfit(Surv(TCGA_time,
                        TCGA_event) ~
                     TCGA_BHLNSX, data = TCGA_samples), pval = T)

TCGA_cpm_p <- ggarrange(plotlist = list(TCGA_GR_cpm_p$plot, TCGA_NIM_cpm_p$plot,
                          TCGA_HOSS_cpm_p$plot, TCGA_PRGIT_cpm_p$plot,
                          TCGA_BHLNSX_cpm_p$plot), labels = c('GR', 'NIM', "HOSS", "PRGIT", "BHLNSX"), ncol = 5) +
  ggtitle("Within TCGA CPM Survival Analysis")

ggsave('plots/prev_method_results/TCGA_cpm.png', plot = TCGA_cpm_p, width = 3000, height = 3000/1.846, units = 'px')
}

#Applying ensemble method to only TCGA
ensemble_TCGA_data <- data.frame(RS = get_mode_per_row(binary_results[gsub('\\.', '-', rownames(binary_results)) %in% TCGA_survival_data$samples$patient,]))

ggsurvplot(survfit(Surv(TCGA_time, TCGA_event) ~ RS, data = ensemble_TCGA_data), pval = T) +
  ggtitle("Validating Ensemble Method with TCGA dataset")


#Survival between TCGA and CGGA
TCGA_v_CGGA_surv <- data.frame(sample = c(rownames(ensemble_validation), TCGA_samples$patient),
                               time = c(ensemble_validation$OS, TCGA_samples$time),
                               event = c(ensemble_validation$Censor..alive.0..dead.1., 
                                         ifelse(TCGA_samples$event == "Dead", 1, 0)),
                               source = c(rep("CGGA", length(ensemble_validation$OS)), 
                                          rep("TCGA", length(TCGA_samples$time))))

TCGA_v_CGGA_surv <-TCGA_v_CGGA_surv[!is.na(TCGA_v_CGGA_surv$time),] 

ggsurvplot(survfit(Surv(TCGA_v_CGGA_surv$time / 365, TCGA_v_CGGA_surv$event) ~
                     TCGA_v_CGGA_surv$source), pval = T, data = TCGA_v_CGGA_surv, risk.table = T) +
  ggtitle("Survival between CGGA and TCGA")
