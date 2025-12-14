# Kernel regression with adaptive bandwidth:
# handles repeated x values and zero-weight boundary cases,
# selects kernel and bandwidth via 10-fold cross-validation,
# and plots the fitted regression curve.



# Inputs
#   x0 - input to be predicted
#   X - vector of training inputs (n x 1)
#   Y - vector of training outputs (n x 1)
#   K - kernel function
#   h - bandwidth function

# Outputs
#   predicted y0 value

kernel_smoothing = function(x0, X, Y, K, h, lambda = 2){
  
  ## Situation 1: When bandwidth is 0
  
  # If there are at least lambda
  # observations with the same x0,
  # the bandwidth would be zero, so we
  # avoid calculating kernel weights.
  # Instead, we predict the value by
  # averaging the corresponding Y values
  
  if (length(which(x0 == X)) >= lambda){
    y0 = mean(Y[which(x0 == X)])  
  } else {
    
    ## Situation 2: If all k nearest
    ## neighbors are exactly on the
    ## boundary
    
    n = length(Y)
    distances = numeric(n)
    
    # Calculate the distance for all
    # observations
    for(i in 1:n){
      distances[i] = abs(X[i] - x0)
    }
    
    # Sort the distances and get the
    #indices of the k nearest neighbors
    sorted_distances = sort(distances)
    neighbours = sort(distances, index.return = TRUE)$ix[1:lambda]
    
    # If all nearest neighbors are on
    # the boundary of the kernel's support,
    # we avoid using weights and instead
    # predict based on the average of Y 
    # values where X is on the boundary.
    if (sum(abs(abs(X[neighbours] - x0) / h(lambda, x0, X))) == lambda){
      y0 = mean(Y[distances <= sorted_distances[lambda]])
    } 
    
    # If the nearest neighbors are not
    # on the boundary, we calculate weights
    else {
      W = K(abs(x0 - X) / h(lambda, x0, X))
      y0 = sum(W * Y) / sum(W)
    }
  }
  return(y0)
}


## Part B ##


# Using the modified kernel smoothing
# function to fit a kernel regression
# model and performing 10-fold CV to
# select a kernel (using Epanechnikov
# and Gaussian)


# load the concrete-data.csv
data = read.csv("concrete-data.csv")
head(data)

# strength as output and age as input
Y = data$strength
X = data$age


# adaptive bandwidth function #


# Input:
#   lambda - positive integer, number of nearest neighbors
#   x0 - scalar, point where we will compute bandwidth
#   X - vector of training observations

# Output: the distance from x0 to its lambda^th nearest neighbor in X

adaptive_bandwidth = function(lambda, x0, X){
  N = length(X)
  d = matrix(0, nrow = N)
  for(i in c(1:N)){
    d[i] = abs(x0 - X[i])
  } 
  d_sorted = sort(d)
  return(d_sorted[lambda])
}


# function to make folds # 


make_folds = function(Y, nFolds, stratified = FALSE, seed = 0){
  
  if(stratified & length(Y) == 1){
    stop('For stratified folds, Y must be a vector of outputs')
  }
  n = ifelse(length(Y) > 1, length(Y), Y)
  index = c(1:n)
  if(stratified){
    Y = factor(Y)
    classes = levels(Y)
    nClasses = length(classes)
    if(nClasses == 1){
      stop('stratified requires more than one class')
    }
    classfolds = list()
    for(class in 1:nClasses){
      classfolds[[class]] = list()
      classIndex = index[Y == classes[class]]
      n_class = sum(Y == classes[class])
      n_per_fold = floor(n_class / nFolds)
      shuffled_index = sample(classIndex)
      for(fold in c(1:(nFolds - 1))){
        classfolds[[class]][[fold]] = shuffled_index[c((1 + (fold - 1) * n_per_fold):(fold * n_per_fold))]
      }
      classfolds[[class]][[nFolds]] = shuffled_index[c(((nFolds - 1)*n_per_fold + 1):n_class)]
    }
    folds = list()
    for(fold in 1:nFolds){
      folds[[fold]] = classfolds[[1]][[fold]]
      for(class in 2:nClasses){
        folds[[fold]] = c(folds[[fold]], classfolds[[class]][[fold]])
      }
    }
  }else{
    folds = list()
    n_per_fold = floor(n / nFolds)
    shuffled_index = sample(index)
    for(fold in c(1:(nFolds - 1))){
      folds[[fold]] = shuffled_index[c((1 + (fold - 1) * n_per_fold):(fold * n_per_fold))]
    }  
    folds[[nFolds]] = shuffled_index[c(((nFolds - 1)*n_per_fold + 1):n)]
  }
  return(folds)
}

# generate folds
n = dim(data)[1]
nFolds = 10
folds = make_folds(n, nFolds)


## Using Gaussian Kernel ##

# Input: t - real number

# Output: D(t) - real number, where D()
# is the Gaussian kernel

gaussian_kernel = function(t){
  return((2*pi)^(-1/2) * exp(-t^2/2))
}


# Define set of lambda values
lambda_guassian_values = seq(from = 1, to = 50, by = 1)
n_lambda_guassian_values = length(lambda_guassian_values)
mse_gaussian = matrix(0, nrow = n_lambda_guassian_values, ncol = nFolds)


# Estimate MSE for each fold/lambda

for(fold in 1:nFolds){
  
  # Define training / testing dataframes based on fold
  training_data = data[-folds[[fold]],]
  testing_data = data[folds[[fold]],]
  
  n_training = dim(training_data)[1]
  n_testing = dim(testing_data)[1]
  
  # Define training and testing data
  trainingX = training_data$age
  trainingY = training_data$strength
  
  testingX = testing_data$age
  testingY = testing_data$strength
  
  # Scale inputs
  testingX = (testingX - mean(trainingX)) / sd(trainingX)
  trainingX = (trainingX - mean(trainingX)) / sd(trainingX)
  
  # Predict testing Y using the guassian kernel
  # smoother with adaptive bandwidth 
  for(i in 1:n_lambda_guassian_values){
    lambda = lambda_guassian_values[i]
    testingyHat = matrix(0, nrow = n_testing)
    for(j in 1:n_testing){
      testingyHat[j] = kernel_smoothing(testingX[j], trainingX, trainingY, gaussian_kernel, 
                                        adaptive_bandwidth, lambda = lambda)
    }
    # Compute MSE
    mse_gaussian[i, fold] = sum((testingY - testingyHat)^2) / n_testing
  }
}

# Average across folds
mse_gaussian = apply(mse_gaussian, 1, mean)

# visually see lambda values v/s MSE
plot(lambda_guassian_values, mse_gaussian)

# Get the optimal lambda value for
# Gaussian kernel 
lambda_guassian = lambda_guassian_values[which.min(mse_gaussian)]
gaussian_min_mse = min(mse_gaussian)

lambda_guassian
gaussian_min_mse


## Using Epanechnikov Kernel ##


# Input: t - real number

# Output: D(t) - real number, where D() is the Epanechnikov kernel

Epanechnikov_kernel = function(t){
  return(as.integer(abs(t) <= 1) * (3/4) * (1 - t^2))
}


# Define set of lambda values
lambda_Epanechnikov_values = seq(from = 1, to = 50, by = 1)
n_lambda_Epanechnikov_values = length(lambda_Epanechnikov_values)
mse_Epanechnikov = matrix(0, nrow = n_lambda_Epanechnikov_values, ncol = nFolds)


# Estimate MSE for each fold/lambda

for(fold in 1:nFolds){
  
  # Define training / testing dataframes based on fold
  training_data = data[-folds[[fold]],]
  testing_data = data[folds[[fold]],]
  
  n_training = dim(training_data)[1]
  n_testing = dim(testing_data)[1]
  
  # Define training and testing data
  trainingX = training_data$age
  trainingY = training_data$strength
  
  testingX = testing_data$age
  testingY = testing_data$strength
  
  # Scale inputs
  testingX = (testingX - mean(trainingX)) / sd(trainingX)
  trainingX = (trainingX - mean(trainingX)) / sd(trainingX)
  
  # Predict testing Y using the Epanechnikov kernel
  # smoother with adaptive bandwidth 
  for(i in 1:n_lambda_Epanechnikov_values){
    lambda = lambda_Epanechnikov_values[i]
    testingyHat = matrix(0, nrow = n_testing)
    for(j in 1:n_testing){
      testingyHat[j] = kernel_smoothing(testingX[j], trainingX, trainingY, Epanechnikov_kernel, 
                                        adaptive_bandwidth, lambda = lambda)
    }
    # Compute MSE
    mse_Epanechnikov[i, fold] = sum((testingY - testingyHat)^2) / n_testing
  }
}

# Average across folds
mse_Epanechnikov = apply(mse_Epanechnikov, 1, mean)

# visually see lambda values v/s MSE
plot(lambda_Epanechnikov_values, mse_Epanechnikov)

# Get the optimal lambda value for
# Epanechnikov kernel 
lambda_Epanechnikov = lambda_Epanechnikov_values[which.min(mse_Epanechnikov)]
Epanechnikov_min_mse = min(mse_Epanechnikov)

lambda_Epanechnikov
Epanechnikov_min_mse


# The MSE for both the kernals is close but
# Epanechnikov wins: 

# Gaussian : 168.1806 and Epanechnikov : 168.0878
# the respective optimal lambda values are 1 and 27!


## PART C ##


# Plot observed data along with the regression
# line for best model


# In order to plot regression line
x0 = seq(from = min(X), to = max(X), by = 0.01)
n0 = length(x0)

# predict values using Epanechnikov kernel
y0_Epanechnikov = matrix(0, nrow = n0)

for(j in c(1:n0)){
  y0_Epanechnikov[j] = kernel_smoothing(x0[j], X, Y, Epanechnikov_kernel, adaptive_bandwidth, lambda = lambda_Epanechnikov)
}

# Plot the observed value and add regression line
plot(X, Y)
lines(x0, y0_Epanechnikov, lwd = 2,
      col = 'darkseagreen3')