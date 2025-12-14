## The goal is to perform a simulation study to 
## explore the bias-variance trade-off by randomly 
## generate data, fit the linear regression model 
## for k = 1 to 15 and compute variance, squared 
## bias and MSE for each k and visually represent it


rm(list=ls())
set.seed(0)


## Simulate the data (given functions)

simulate_X = function(n, p) {
  X = matrix(2 * runif(n * p) - 1, nrow = n, ncol = p)
  X = cbind(matrix(1, nrow = n), X)  
  return(X)
}

simulate_Y = function(X, B, sigma) {
  n = dim(X)[1]
  E = matrix(rnorm(n, mean = 0, sd = sigma), nrow = n)
  Y = X %*% B + E
  return(Y)
}


p = 15 #no of inputs
n = 30 #size
n_iterations = 1000  #no of iteration 
sigma = 1 


## True model

B = matrix(0, nrow = p + 1)

# loop to compute beta values
for(j in c(1:p)) {
  B[j + 1] = 1 / j  #beta = 1/j
}


## x_new value to be predicted

x_new = matrix(1/2, ncol = p + 1)
x_new[1] = 1

# Initialize matrix to store results
hatY = matrix(0, nrow = n_iterations, ncol = p)


## Simulation loop

for (i in 1:n_iterations) {
  
  # Simulate training data (X and Y)
  X = simulate_X(n, p)
  y = simulate_Y(X, B, sigma)
  
  # Fit linear regression models with k predictors
  for (k in 1:p) {
    x = X[, 1:(k + 1)]  
    
    bHat = solve(t(x) %*% x)%*% t(x) %*% y  
    
    hatY[i, k] = x_new[1:(k + 1)] %*% bHat
  }
}


## Computing the variance, squared bias and MSE


# Compute the true value of y_new
Ey_new = x_new %*% B

# initialize vector to store the calculated value
sq_bias = rep(0, p)
variance = rep(0, p)
MSE = rep(0, p)

# computing all the values for k = 1 to 15 using loop
for (k in 1:p) {
  
  variance[k] = var(hatY[, k])
  
  sq_bias[k] = (mean(hatY[, k]) - Ey_new)^2
  
  MSE[k] = sigma^2 + variance[k] + sq_bias[k]
  
}


## Plot the variance, squared-bias and MSE as a
## function of k

plot(1:p, variance, type = "o", 
     col = "forestgreen", pch = 1, 
     xlab = "k",  ylab = "",
     ylim = c(0, max(c(variance, sq_bias, MSE))))

lines(1:p, sq_bias, type = "o" ,
      col = "navyblue", pch = 1)

lines(1:p, MSE, type = "o",
      col = "indianred", pch = 1)

legend("topright", 
       legend = c("Variance", "Squared Bias", "MSE"), 
       col = c("forestgreen", "navyblue", "indianred"),
       lty = 1 , cex = 0.8 )

