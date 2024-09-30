## Random Forest Reproductive Mode ##

library(randomForest)
library(scales)
library(caret)
library(rsample)
library(Metrics)
library(arrow)

#### Reproductive mode RF #### 


asex_sex_data_GLM2 <- read.csv2("/data/asex_sex_data_GLM_nocategory.csv")

# Split the data into training and testing sets
#initial split method prop 0.6
set.seed(42)
data_split <- initial_split(asex_sex_data_GLM2, prop = 0.6)
trainData <- training(data_split)
trainData$Reproductive.mode<-as.factor(trainData$Reproductive.mode)
testData <- testing(data_split)
testData$Reproductive.mode <- as.factor(testData$Reproductive.mode)
class(trainData$Reproductive.mode) #in this case it should be factor with two levels asex and sex (binary)
print(testData)

#For categories I shouldn't use RMSE but a different error approach, in this case OOB 

# Random Forest models
set.seed(42)
rfM0 <- randomForest(Reproductive.mode ~ dem_ele+MAP_chelsa+Latitude, data = trainData, importance = TRUE) #+veg_cover_100m+per_Phanerophytes

# Evaluate models
print(rfM0)


# Variable importance
importance(rfM0)
varImpPlot(rfM0)

#Predictions

x_pred = read_parquet("/data/x_pred.parquet")

mod_pred<-predict(rfM0, x_pred, type='prob')

write.csv(mod_pred[,2], file='../results/pred_rf_reproductivemode.csv')
