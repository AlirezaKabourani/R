# Polynomial Regression

# Importing the dataset
dataset = read.csv('E:/My Workspace/R Programming/R Anjoman/Part 2 - Regression/Section 6 - Polynomial Regression/R/Position_Salaries.csv')
dataset = dataset[2:3]

plot(dataset$Salary)

# Splitting the dataset into the Training set and Test set

# Feature Scaling

# Fitting Linear Regression to the dataset
lin_reg = lm(formula = Salary ~ ., data = dataset)

# Fitting Polynomial Regression to the dataset
dataset$Level2 = dataset$Level^2
dataset$Level3 = dataset$Level^3
dataset$Level4 = dataset$Level^4
poly_reg = lm(formula = Salary ~ .,
              data = dataset)

# Visualising the Linear Regression results
library(ggplot2)
ggplot() +
  geom_point(aes(x = dataset$Level, y = dataset$Salary),
             colour = 'red') +
  geom_line(aes(x = dataset$Level, y = predict(lin_reg, newdata = dataset)),
            colour = 'blue') +
  ggtitle('Truth or Bluff (Linear Regression)') +
  xlab('Level') +
  ylab('Salary') + 
  theme_minimal()

ggplot(data = dataset, aes(x = Level, y = Salary)) +
  geom_point(colour = 'blue', size = 2) + 
  geom_smooth(method = lm, colour = 'red', linetype ='dashed') + 
  ggtitle('Truth or Bluff (Linear Regression)') +
  scale_y_continuous(breaks = seq(0,max(dataset$Salary),100000)) +
  scale_x_continuous(breaks = 1:10)+
  theme_minimal()

# Visualising the Polynomial Regression results
library(ggplot2)
ggplot() +
  geom_point(aes(x = dataset$Level, y = dataset$Salary),
             colour = 'red') +
  geom_line(aes(x = dataset$Level, y = predict(poly_reg, newdata = dataset)),
            colour = 'blue') +
  ggtitle('Truth or Bluff (Polynomial Regression)') +
  xlab('Level') +
  ylab('Salary')

# Visualising the Regression Model results (for higher resolution and smoother curve)
library(ggplot2)
x_grid = seq(min(dataset$Level), max(dataset$Level), 0.1)
ggplot() +
  geom_point(aes(x = dataset$Level, y = dataset$Salary),
             colour = 'red') +
  geom_line(aes(x = x_grid, y = predict(poly_reg,
                                        newdata = data.frame(Level = x_grid,
                                                             Level2 = x_grid^2,
                                                             Level3 = x_grid^3,
                                                             Level4 = x_grid^4))),
            colour = 'blue') +
  ggtitle('Truth or Bluff (Polynomial Regression)') +
  xlab('Level') +
  ylab('Salary') +
  theme_minimal()

# Predicting a new result with Linear Regression
predict(lin_reg, data.frame(Level = 6.5))

# Predicting a new result with Polynomial Regression
predict(poly_reg, data.frame(Level = 6.5,
                             Level2 = 6.5^2,
                             Level3 = 6.5^3,
                             Level4 = 6.5^4))

AIC(lin_reg)
AIC(poly_reg)
BIC(lin_reg)
BIC(poly_reg)
