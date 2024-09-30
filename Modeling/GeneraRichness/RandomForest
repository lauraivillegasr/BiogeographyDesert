##Random Forest genera richness ##

library(randomForest)
library(ggplot2)
library(scales)
library(caret)
library(rsample)
library(Metrics)
library(rpart)
library(arrow)

#Load data
ecovar_data_GLM_nocategories <- read.csv2("/data/EcoVarData_generacurated_nocategories.csv")
#Normalize same variables that were normalized for the GLS
ecovar_data_GLM_nocategories$soil_thick_log<-(log((ecovar_data_GLM_nocategories$soil_thick)+1))

# Split data into training and test sets
#initial split method prop 0.6
set.seed(42)
data_split <- initial_split(ecovar_data_GLM_nocategories, prop = 0.6)
train_data <- training(data_split)
test_data <- testing(data_split)
class(train_data$composition_genera) #for regression should not be factor!
print(test_data)

# Set up the train control
train_control <- trainControl(method = "cv", number = 10)

# Define the grid of hyperparameters to search
tune_grid <- expand.grid(.mtry = c(1:6), 
                         .splitrule = "variance", 
                         .min.node.size = c(1, 3, 5))

# Train the Random Forest model
set.seed(42)
rf_model <- train(composition_genera ~Latitude+dem_ele+rock_type_int+soil_thick_log+MAP_chelsa+range_tmp_chelsa, data = train_data, 
                  method = "ranger",
                  trControl = train_control,
                  tuneGrid = tune_grid)

# Check the results
print(rf_model)
plot(rf_model)

# Make predictions and evaluate the model
predictions <- predict(rf_model, test_data)
predictions<- round(predictions)
predictions <- as.factor(predictions)

print(predictions)
actual<- as.factor(test_data$composition_genera)

#somehow th eocnfussion matrix doesn't work, it says I need to ensure both factors have same levels

# So...ensure both factors have the same levels
common_levels <- union(levels(predictions), levels(actual))
predictions <- factor(predictions, levels = common_levels)
actualnull <- factor(actual, levels = common_levels)

confmatrix<-confusionMatrix(predictions, actualnull)
print(confmatrix)

best_mtry <- rf_model$bestTune$mtry
print(paste("The optimal mtry value is:", best_mtry))
#the result was mtry=1! this validation was done for the null model



# Function to calculate RMSE
calculate_rmse <- function(model, test_data) {
  predictions <- predict(model, newdata = test_data)
  actuals <- test_data$composition_genera
  sqrt(mean((actuals - predictions)^2))
}


#Training all the models and testing them by using the RMSE function 


# Model 1: Null model
set.seed(42)
rf0_1 <- randomForest(composition_genera ~ Latitude+dem_ele+rock_type_int+soil_thick_log+MAP_chelsa+range_tmp_chelsa, data = train_data, importance = TRUE,mtry=1)
rmse_rf0_1 <- calculate_rmse(rf0_1, test_data)
variance_rf0_1 <- mean(rf0_1$rsq) * 100
print(rf0_1)
#plot(rf0_1)


#variable importance in NULL model

# Get variable importance from the model fit

importance(rf0_1)

#Predictions   


x_pred = read_parquet("/data/x_pred.parquet")

mod_pred<-predict(rf0_1, x_pred)

write.csv(mod_pred, file='../results/pred_rf_richness.csv')
