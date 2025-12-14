### STAT 5P87 LAB 2 ###

### Implementing the  quadratic  discriminant ###
### analysis (QDA) classifier  as a nonlinear ###
### classification method  ###

# Srushti Desai #


rm(list = ls())
setwd("/Users/srushtidesai/Desktop/School/COURSES/STAT 5P87/data_files")


# load the training and testing data
training_data = read.csv('vowel_training.csv')
testing_data = read.csv('vowel_testing.csv')


# y as input and x1-x9 as outputs
Y = factor(training_data$y)
X = model.matrix(y ~ 0 + ., data=training_data)


# parameters 
p = dim(X)[2]
n = length(Y)
classes = levels(Y)
K = length(classes)


# Get number of observations in each class
# and estimate proportions
nK = matrix(0, nrow = K)
qK = matrix(0, nrow = K)

for(k in 1:K){
  nK[k] = sum(Y == classes[k])
  qK[k] = nK[k] / n
}


# Estimate a mean vector muHat_k for each class
muHat = matrix(NA, ncol = p, nrow = K)
for(k in 1:K){
  muHat[k,] = apply(X[Y == classes[k],], 2, mean)
}


## Estimate covariance matrix

# array to store our k covariance matrix
SigmaHat = array(0, dim = c(p, p, K))  

for(k in 1:K){
  muMatrix = matrix(rep(muHat[k,], nK[k]), 
                    nrow = nK[k], ncol = p,
                    byrow = TRUE)
  SigmaHat[,,k] = t(X[Y == classes[k],] - muMatrix) %*% (X[Y == classes[k], ] - muMatrix) / (nK[k] - 1)
}


## Compute discriminant for each observation / class

# Compute inverse covariance matrices
invSigma = array(0, dim = c(p, p, K))

for (k in 1:K) {
  invSigma[,,k] = solve(SigmaHat[,,k])
}

# I am using the log transformed formula for delta
delta = matrix(NA, nrow = n, ncol = K)

for (i in 1:n) {
  for (k in 1:K) {
    delta[i, k] = log(qK[k]) - 0.5 * log(det(SigmaHat[,,k])) + (X[i, ] %*% invSigma[,,k] %*% muHat[k, ])  - (0.5 * (t(X[i, ]) %*% invSigma[,,k] %*% X[i, ]))  - (0.5 * (muHat[k, ] %*% invSigma[,,k] %*% muHat[k, ]))
  }
}

# Predict all observation
yHat = apply(delta, 1, which.max)

# Compute training accuracy
mean(yHat == Y) #0.9791667



###--------------###



### Testing Accuracy



# Load testing data
testing_data = read.csv('vowel_testing.csv')

# y as input and x1-x9 as outputs
Y_testing = factor(testing_data$y)
X_testing = model.matrix(y ~ 0 + ., data=testing_data)

# parameter
n_testing = length(Y_testing)


## Compute discriminant for each observation / class

delta_testing = matrix(NA, nrow = n_testing, ncol = K)

for (i in 1:n_testing) {
  for (k in 1:K) {
    delta_testing[i, k] = log(qK[k]) - 0.5 * log(det(SigmaHat[,,k])) - (0.5 * (t(X_testing[i, ]) %*% invSigma[,,k] %*% X_testing[i, ])) + (X_testing[i, ] %*% invSigma[,,k] %*% muHat[k, ])-(0.5 * (muHat[k, ] %*% invSigma[,,k] %*% muHat[k, ]))
  }
}
 

# Find class that maximizes discriminant
yHat_testing = classes[apply(delta_testing, 1, which.max)]

# Compute testing accuracy
mean(yHat_testing == Y_testing)


## the testing accuracy is : 0.5043 that is 50.43 %



###--------------###



# Trying to Verify my solutions using MASS 
# package
 
library(MASS)

# Fitting QDA model for training data
model = qda(y ~ ., data = training_data)

# Predict on training data
yHat_training = predict(model, training_data)

# Compute training accuracy
mean(yHat_training$class == training_data$y)

# training accuracy : 0.9791667 (same as my ans)


# Predict on testing data
yHat_testing1 = predict(model, testing_data)

# Compute testing accuracy
mean(yHat_testing1$class == testing_data$y)

# testing accuracy : 0.504329
