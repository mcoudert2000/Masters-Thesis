#Binary Classifier IDH
library(tidyverse)
library(caret)
library(caretEnsemble)
library(psych)
library(Amelia)
library(mice)
library(GGally)
library(rpart)
library(randomForest)
library(modelr)
library(broom)
library(ISLR)
library(e1071)

set.seed(180018876)

load('data/TCGA_CGGA_diff_expressed_genes.rdata')
load('data/normalized_TCGA_CGGA_IDH.rdata')

IDH_train <- CGGA_v
IDH_test <- TCGA_v

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

str(IDH_train)
str(IDH_test)

indxTrain <- createDataPartition(y = IDH_train$outcome,p = 0.75,list = FALSE)
training <- IDH_train[indxTrain,]
testing <- IDH_train[-indxTrain,] #Check dimensions of the split > prop.table(table(data$Outcome)) * 100

prop.table(table(training$outcome)) * 100

x = training[,1:5]
y = training$outcome


model = train(x,y,'rpart',trControl=trainControl(method='cv',number=5))
model
confusionMatrix(predict(model, newdata = testing), testing$outcome)

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

require("FactoMineR")
require("factoextra")
IDH_PCA <- PCA(IDH_visualization_data[,1:5], scale.unit = T, graph = F)

IDH_PCA

Source <- IDH_visualization_data$source
Model_Prediction <- IDH_visualization_data$predict

d <- as.data.frame(IDH_PCA$ind$coord)
d$label <- rownames(IDH_PCA)

IDH_PCA_plot <- fviz_pca_ind(IDH_PCA, geom.var = T, geom.ind = F, repel = T) +
  geom_text(data = d[1:9,], aes(x = Dim.1, y = Dim.2, label = label), hjust = -.1, vjust =-.1) +
  geom_point(data = d, aes(x = Dim.1, y = Dim.2, col = Model_Prediction, pch = Source)) +
  ggtitle("PCA Plot of 5 Key Genes Used in IDH-Signature")

ggsave(IDH_PCA_plot, filename = 'plots/PCA Plot of 5 Key Genes Used in IDH-Signature.png')  

