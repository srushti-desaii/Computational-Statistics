## The goal is to fit the model using ridge
## regression and to find optimal lambda that 
## minimizes the testing error


rm(list=ls())

## Loading the Iris data 

data = iris

n = dim(data)[1]
n_training = n/2
n_testing = n - n_training
p = 3 #number of inputs


## Split into training and testing data.frame

training_data = data[seq(1, nrow(data), by = 2), ]
testing_data = data[seq(2, nrow(data), by = 2), ]

# Sepal.Length as output
trainingY = training_data$Sepal.Length
testingY = testing_data$Sepal.Length

# Sepal.Width, Petal.Length, and Petal.Width as inputs
trainingX = model.matrix(Sepal.Length ~ 0 + 
                           Sepal.Width + Petal.Length +
                           Petal.Width, data=training_data)
testingX = model.matrix(Sepal.Length ~ 0 + 
                          Sepal.Width + Petal.Length +
                          Petal.Width, data=testing_data)


## Scale the data

xBar = apply(trainingX, 2, mean)
s = apply(trainingX, 2, sd)

trainingZ = t((t(trainingX) - xBar) / s)
testingZ = t((t(testingX) - xBar) / s)

# Verify the scaled data
round(apply(trainingZ, 2, mean))
apply(trainingZ, 2, sd)


## Find the optimal lambda using testing data

lambda_values = seq(from = 0.6, to = 1.5, by = 0.001)
n_lambda_values = length(lambda_values)
mse = matrix(NA, nrow = n_lambda_values)

# loop to compute lambda values

for(i in 1:n_lambda_values){
  lambda = lambda_values[i]
  
  b0Hat = mean(trainingY)
  
  bHat = solve(t(trainingZ) %*% trainingZ + 
                 lambda * diag(p)) %*% t(trainingZ) %*%
    (trainingY - mean(trainingY))
  
  # Make predictions on testing data
  yHat = b0Hat + testingZ %*% bHat
  
  # Estimate mean-squared error
  mse[i] = mean((testingY - yHat)^2)
}


## Visual representation for lambda_values in order to 
## search for our range

plot(lambda_values, mse, bty = 'n', 
     lwd = 2, cex = 1.2)


## index that minimizes the vector

which.min(mse) #396


## Optimal lambda value
lambda_values[which.min(mse)] #0.995


### The optimal lambda value is 0.995
