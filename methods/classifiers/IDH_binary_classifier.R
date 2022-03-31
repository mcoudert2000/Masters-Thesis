#Binary Classifier IDH
library(tidyverse)
library(caret)
library(caretEnsemble)

library(rpart)
library(e1071)

set.seed(180018876)

load('data/TCGA_CGGA_diff_expressed_genes.rdata')
load('data/normalized_TCGA_CGGA_IDH.rdata')

IDH_train <- CGGA_v
IDH_test <- TCGA_v

#Split TCGA and CGGA into test and train set
IDH_train <- IDH_train[rownames(IDH_train) %in% diff_expressed_genes$CGGA,]
IDH_test <- IDH_test[rownames(IDH_test) %in% diff_expressed_genes$CGGA,]

outcome_train <- as.factor(IDH_train$design[,1])
outcome_test <- as.factor(IDH_test$design[,1])

IDH_train <- data.frame(t(IDH_train$E))
IDH_train$outcome <- outcome_train

IDH_test <- data.frame(t(IDH_test$E))
IDH_test$outcome <- outcome_test

table(IDH_train$outcome)
table(IDH_test$outcome)


indxTrain <- createDataPartition(y = IDH_train$outcome,p = 0.75,list = FALSE)
training <- IDH_train[indxTrain,]
testing <- IDH_train[-indxTrain,] #Check dimensions of the split > prop.table(table(data$Outcome)) * 100

prop.table(table(training$outcome)) * 100

x = training[,1:5]
y = training$outcome

#Train a decision tree with 5 cross validation folds
model = train(x,y,'rpart',trControl=trainControl(method='cv',number=5),
              parms = list(split = "gini"))
model
#confusiong matrix on CGGA test set
confusionMatrix(predict(model, newdata = testing), testing$outcome)

#Confusion matrix on TCGA test set
confusionMatrix(predict(model, newdata = IDH_test[,1:5]), IDH_test$outcome)


### Apply IDH classifier to ASTRID's data -----
load('data/Astrid_data_counts.rdata')
Astrid_v <- voom(Astrid_data)

Astrid_v <- Astrid_v[rownames(Astrid_v) %in% diff_expressed_genes$CGGA,]
Astrid_test <- data.frame(t(Astrid_v$E))

#Interesting it mostly predicts them all as IDH_WT except for 
astrid_IDH_results <- predict(model, newdata = data.frame(t(Astrid_v$E))) 

astrid_IDH_results <- data.frame(sample = colnames(Astrid_v), IDH_type = ifelse(astrid_IDH_results == 0, 'IDH_WT', 'IDH_Mutant'))
astrid_IDH_results


#Visualization ------
Astrid_test$source = 'Astrid'
IDH_train$source = 'CGGA'
IDH_test$source = 'TCGA'

Astrid_test$predict <- astrid_IDH_results$IDH_type
IDH_train$predict <- ifelse(predict(model, newdata = IDH_train[,c(1:5)]) == 0, "IDH_WT", "IDH_Mutant")
IDH_test$predict <- ifelse(predict(model, newdata = IDH_test[,c(1:5)]) == 0, "IDH_WT", "IDH_Mutant")

IDH_visualization_data <- rbind(Astrid_test, IDH_train[,c(1:5,7,8)], IDH_test[,c(1:5,7,8)])
IDH_visualization_data$source <- as.factor(IDH_visualization_data$source)

IDH_visualization_data$IDH <- c(rep("IDH_WT",9), ifelse(IDH_train$outcome == 1, "IDH_Mutant", "IDH_WT"), ifelse(IDH_test$outcome == 1, "IDH_Mutant", "IDH_WT"))

IDH_visualization_data
library(rpart.plot)
rpart.plot(model$finalModel)

#Checking % of samples that fit the condition in training and test set
sum(training$RBP1>=5.7) / length(training$RBP1)
sum(testing$RBP1>=5.7) / length(testing$RBP1)

require("FactoMineR")
require("factoextra")
IDH_PCA <- PCA(IDH_visualization_data[,1:5], scale.unit = T, graph = F)

IDH_PCA2 <- PCA(IDH_visualization_data[,c("RBP1", "HMX1")], scale.unit = T, graph = F)

IDH_PCA

Source <- IDH_visualization_data$source
Model_Prediction <- IDH_visualization_data$predict

d <- as.data.frame(IDH_PCA$ind$coord)
d$label <- rownames(IDH_PCA)

IDH_PCA_plot <- fviz_pca_ind(IDH_PCA, geom.var = T, geom.ind = F, repel = T) +
  geom_text(data = d[1:9,], aes(x = Dim.1, y = Dim.2, label = label), hjust = -.1, vjust =-.1) +
  geom_point(data = d, aes(x = Dim.1, y = Dim.2, col = IDH_visualization_data$IDH, pch = Source)) +
  ggtitle("PCA Plot of 5 Key Genes Used in IDH-Signature")

ggsave(IDH_PCA_plot, filename = 'plots/PCA Plot of 5 Key Genes Used in IDH-Signature.png')  


d2 <- as.data.frame(IDH_PCA2$ind$coord)
d2$label <- rownames(IDH_PCA2)

IDH_PCA2_plot <- fviz_pca_ind(IDH_PCA2, geom.var = T, geom.ind = F, repel = T) +
  geom_text(data = d2[1:9,], aes(x = Dim.1, y = Dim.2, label = label), hjust = -.1, vjust =-.1) +
  geom_point(data = d2, aes(x = Dim.1, y = Dim.2, col = Model_Prediction, pch = Source)) +
  ggtitle("PCA Plot of 2 Key Genes Used in IDH-Signature")

ggplot(IDH_visualization_data[,c("RBP1", "HMX1")]) +
  geom_point(aes(x = RBP1, y = HMX1, col = IDH_visualization_data$IDH, pch = Source)) +
  geom_segment(aes(x = 5.7, y = 2.2, xend = -5, yend = 2.2)) +
  geom_segment(aes(x = 5.7, y = 2.2, xend = 5.7, yend = 10)) +
  ggtitle("Plot of 2 Key Genes Used in IDH-Signature with Boundary Lines") +
  labs(col = "IDH-Type")


IDH_PCA2_plot
#Validation of survival differences between two groups

require(survival)
require(survminer)
load('data/TCGA_data_survival.rdata')

event <- ifelse(TCGA_survival_data$samples$event == "Dead", 1, 0)
time <- TCGA_survival_data$samples$time
IDH <- ifelse(TCGA_survival_data$samples$IDH == "WT", "WT", "Mutant")

survival_TCGA_IDH <- data.frame(event, time, IDH)

ggsurvplot(survfit(Surv(time, event) ~ IDH), data = survival_TCGA_IDH, pval = T) +
  ggtitle("Survival stratified by IDH_type within the TCGA")


load('data/CGGA/CGGA_data.RDATA')

CGGA_survival_samples <- CGGA_data$samples[!is.na(CGGA_data$samples$OS) &
                                             !is.na(CGGA_data$samples$Censor..alive.0..dead.1.),]


IDH_CGGA <- CGGA_survival_samples$IDH_mutation_status
event_CGGA <- CGGA_survival_samples$Censor..alive.0..dead.1.
time_CGGA <- CGGA_survival_samples$OS

survival_CGGA_IDH <- data.frame(event = event_CGGA,time = time_CGGA,IDH = IDH_CGGA)
ggsurvplot(survfit(Surv(time_CGGA, event_CGGA) ~ IDH_CGGA), data = survival_CGGA_IDH, pval = T) +
  ggtitle("Survival stratified by IDH_type within the CGGA")

survival_CGGA_IDH$source <- rep("CGGA", length(CGGA_survival_samples[,1]))
survival_TCGA_IDH$source <- rep("TCGA", length(TCGA_survival_samples[,1]))

survival_TCGA_IDH$IDH <- ifelse(survival_TCGA_IDH$IDH == "WT", "Wildtype", "Mutant")

survival_combined_IDH <- rbind(survival_TCGA_IDH, survival_CGGA_IDH)

#Survival for combined
ggsurvplot(survfit(Surv(survival_combined_IDH$time, survival_combined_IDH$event) ~ survival_combined_IDH$IDH), data =survival_combined_IDH, pval = T, risk.table = T) +
  ggtitle("Survival stratified by IDH type in full dataset")


#Seeing if misclassification is associated with tumor content

load('results/estimate_scores.rdata')

rownames(estimate_scores)
rownames(IDH_test) %in% rownames(estimate_scores)

IDH_misclassified <- rownames(IDH_visualization_data[IDH_visualization_data$predict != IDH_visualization_data$IDH,])[-(1:9)]
IDH_misclassified
IDH_classified <- rownames(IDH_visualization_data[IDH_visualization_data$predict == IDH_visualization_data$IDH,])[-(1:9)]

t.test(estimate_scores[IDH_misclassified,]$tumor_purity, estimate_scores[IDH_classified,]$tumor_purity)
