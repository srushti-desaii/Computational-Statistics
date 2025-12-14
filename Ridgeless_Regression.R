### STAT 5P87 LAB 5 ###

### Implementing Ridgeless Regression ###

# Srushti Desai #


## PART A ##
# Implementing ridgeless regression #


rm(list=ls())
set.seed(0)


n_training = 200
n_testing = 100

# range for gamma values
gamma_values = c(seq(0.1, 0.95, by = 0.1), seq(0.95, 1.05, by = 0.01), seq(1.05, 2, by = 0.1))

# maximum iteration
max_iteration = 6

# to converge the model
threshold = 1e-2  

# initialize matrix to store MSE for each
# iteration 
MSE_testing = matrix(0, nrow = length(gamma_values), ncol = max_iteration)


for (iter in 1:max_iteration) {
  for (i in 1:length(gamma_values)) {
    gamma = gamma_values[i]
    p = round(gamma * n_training)  
    
    # Generate training data and 
    # testing X 
    trainingX = matrix(rnorm(n_training * p), nrow = n_training, ncol = p)
    testingX = matrix(rnorm(n_testing * p), nrow = n_testing, ncol = p)
    
    # generate Beta values
    B = rnorm(p)
    B = B / sqrt(sum(B^2))  
    
    # Compute training and testing Y 
    trainingY = trainingX %*% B + rnorm(n_training)
    testingY = testingX %*% B + rnorm(n_testing)
    
    # Scale training data and testing
    # data
    
    trainingX_mean = colMeans(trainingX)
    trainingX_sd = apply(trainingX, 2, sd)
    
    trainingX = scale(trainingX, center = trainingX_mean, scale = trainingX_sd)
    testingX = scale(testingX, center = trainingX_mean, scale = trainingX_sd)
    
   
    # Ridgless regression 
    I = diag(p)  
    lambda = 1e-3  
    bHat_prev = rep(0, p)  
    b0 = mean(trainingY)
    
    repeat {
      bHat = solve(t(trainingX) %*% trainingX + lambda * I) %*% t(trainingX) %*% (trainingY - b0)
      
      # check for convergence 
      if (max(abs(bHat - bHat_prev)) < threshold) break
      
      bHat_prev = bHat
      lambda = lambda * 0.95 
    }
    
    # Predict yHat values
    yHat = b0 + testingX %*% bHat
    
    # Compute MSE
    MSE_testing[i, iter] = mean((testingY - yHat)^2)
  }
}

# Average MSE across iterations
MSE = rowMeans(MSE_testing)


## PART B ##
# Plot the Double Descent Curve #


plot(gamma_values, MSE, type = "o", 
     pch = 1, col = "navyblue", 
     xlab = "Gamma", 
     ylab = "MSE", 
     main = "Double Descent Curve", 
     lty = 1, 
     cex = 1) 

# i wanted my lines to be green :)
lines(gamma_values, MSE, col = "forestgreen", lwd = 2) 

