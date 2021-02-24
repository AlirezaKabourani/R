# Multiple Linear Regression

# Importing the dataset
dataset = read.csv('E://My Workspace//R Programming//R Anjoman//Part 2 - Regression//Section 5 - Multiple Linear Regression//R//50_Startups.csv')

# Encoding categorical data
dataset$State = factor(dataset$State,
                       levels = c('New York', 'California', 'Florida'),
                       labels = c(1, 2, 3))
unique(dataset$State)
# Splitting the dataset into the Training set and Test set
# install.packages('caTools')
# library(caTools)
set.seed(123)
dim(dataset)
split = sample(1:50, 50 * 0.8, replace = FALSE)
training_set = dataset[split, ]
test_set = dataset[-split, ]

# Feature Scaling
# training_set = scale(training_set)
# test_set = scale(test_set)

#-----------------------------------------------------------------------

# Fitting Multiple Linear Regression to the Training set
regressor = lm(formula = Profit ~ .,
               data = training_set)
summary(regressor)

# Predicting the Test set results
y_pred = predict(regressor, newdata = test_set)

data.frame(TrueProfit = test_set$Profit, PredictedProfit = y_pred)
#--------------------------------------------------------------------------

#  variable selection
base_model = lm(Profit ~ R.D.Spend, data = training_set)
step_model = step(base_model,scope = list(upper = regressor, lower = ~1), 
                  direction = "forward",trace = T)

forward_pred = predict(step_model, newdata = test_set)
data.frame(TrueProfit = test_set$Profit, PredictedProfit = forward_pred)

library(MLmetrics)
MSE(y_pred, test_set$Profit)
MSE(forward_pred, test_set$Profit)
AIC(regressor)
AIC(step_model)
summary(step_model)
base_model = lm(Profit ~ ., data = training_set)
step_model = step(base_model,scope = list(upper = regressor, lower = ~1), 
                  direction = "backward",trace = T)

base_model = lm(Profit ~ R.D.Spend, data = training_set)
step_model = step(base_model,scope = list(upper = regressor, lower = ~1), 
                  direction = "both",trace = T)
