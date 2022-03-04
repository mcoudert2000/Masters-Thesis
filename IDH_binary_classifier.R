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

IDH_key_genes_expression <- scale(v$E[rownames(v$E) %in% rownames(dt)[de.common],])

outcome <- (ifelse(TCGA_sample_data_no_NA$IDH_status == 'WT', 0, 1))



idh_data <- data.frame(t((IDH_key_genes_expression)))
idh_data$outcome <- factor(outcome, levels = c(0,1))
#View(idh_data)

str(idh_data)


indxTrain <- createDataPartition(y = idh_data$outcome,p = 0.75,list = FALSE)
training <- idh_data[indxTrain,]
testing <- idh_data[-indxTrain,] #Check dimensions of the split > prop.table(table(data$Outcome)) * 100

prop.table(table(training$outcome)) * 100

x = training[,-9]
y = training$outcome


model = train(x,y,'nb',trControl=trainControl(method='cv',number=5), metric = 'Neg Pred Value')
model
confusionMatrix(predict(model, newdata = testing), testing$outcome)
