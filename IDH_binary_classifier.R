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

str(IDH_train)
str(IDH_test)

indxTrain <- createDataPartition(y = IDH_train$outcome,p = 0.75,list = FALSE)
training <- IDH_train[indxTrain,]
testing <- IDH_train[-indxTrain,] #Check dimensions of the split > prop.table(table(data$Outcome)) * 100

prop.table(table(training$outcome)) * 100

x = training[,1:5]
y = training$outcome


model = train(x,y,'nb',trControl=trainControl(method='cv',number=5), metric = 'Neg Pred Value')
model
confusionMatrix(predict(model, newdata = testing), testing$outcome)

confusionMatrix(predict(model, newdata = IDH_test[,1:5]), IDH_test$outcome)
