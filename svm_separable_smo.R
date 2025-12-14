### STAT 5P87 LAB 4 ###

### Implementing the sequential minimum ###
### optimization (SMO) algorithm for estimating ###
### the separable hyper plane model ###

# Srushti Desai #

rm(list=ls())
setwd("/Users/srushtidesai/Desktop/School/COURSES/STAT 5P87/data_files")

set.seed(0)

## PART A ##


# Applying SMO algorithm #


# load the separable data
data = read.csv("separable-data.csv")

# x1 and x2 as inputs 
X = model.matrix(y ~ 0 + . , data = data)

# y as input
Y = data$y

# defining parameters
n = dim(X)[1]
p = dim(X)[2]

# initialize alpha matrix
alpha = matrix(1/n, n , 1) 

# threshold where i want to converge 
threshold = 1e-6


max_iteration = 100

for (iter in 1:max_iteration) {
  
  # keep store new alpha values
  alpha_old = alpha  
  
  # iterate over different pairs
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      
      # compute bHat
      bHat = t(X) %*% (alpha * Y)
      
      # compute b0
      s = which(alpha > 0) 
      if (length(s) > 0) { 
        bHat_0 = (1 / Y[s]) - (X[s, ] %*% bHat)
        b0 = mean(bHat_0, na.rm = TRUE)
      }
      
      # Compute errors using the updated beta
      Ei = X[i, ] %*% bHat + b0 - Y[i]
      Ej = X[j, ] %*% bHat + b0 - Y[j]
      
      # Compute the eta 
      eta = t(X[i, ]) %*% X[i, ] + t(X[j, ]) %*% X[j, ] - 2 * t(X[i, ]) %*% X[j, ]
      
      # calculate alpha "j" value
      alpha_J = alpha[j] + (Y[j] * (Ei - Ej)) / eta
      
      # check on the constraints given 
      if (Y[i] == Y[j]) {
        if (alpha_J < 0) {
          alpha_J_new = 0
        } else if (alpha_J > (alpha[j] + alpha[i])) {
          alpha_J_new = alpha[j] + alpha[i]
        } else {
          alpha_J_new = alpha_J
        }
      } else {
        L = max(0, alpha[j] - alpha[i])
        if (alpha_J <= L) {
          alpha_J_new = L
        } else {
          alpha_J_new = alpha_J
        }
      }
      
      # calculate alpha "i" value
      # conditions fulfilled by checking on "j"
      alpha_I_new = alpha[i] + Y[i] * Y[j] * (alpha[j] - alpha_J_new)
      
      # update alpha values
      alpha[i] = alpha_I_new
      alpha[j] = alpha_J_new
      
    }
  }
  
  # check for convergence: if all alpha values remains
  # within the threshold we break the loop 
  if (max(abs(alpha - alpha_old)) < threshold) {
    break
  }
}


# compute final bHat
bHat = t(X) %*% (alpha * Y)
bHat 

# compute final b0
s = which(alpha > 0) 
if (length(s) > 0) { 
  bHat_0 = (1 / Y[s]) - (X[s, ] %*% bHat)
  b0 = mean(bHat_0, na.rm = TRUE)
}
b0


## PART B ##


# a.) Plot the resulting decision boundary #


x1range = c(min(data$x1), max(data$x1))
x2range = c(min(data$x2), max(data$x2))


plot(data$x1[data$y == -1], data$x2[data$y == -1], xlim = x1range, ylim = x2range, 
     bty = 'n', col = 'forestgreen', lwd = 2, xlab = 'x1', ylab = 'x2')

points(data$x1[data$y == 1], data$x2[data$y == 1], col = 'navyblue', lwd = 2)

points(X[alpha != 0,1], X[alpha != 0,2],col='indianred',cex=3, lwd = 2)

abline(a=-b0/bHat[2], b=-bHat[1]/bHat[2], col="black", lwd = 2)


# b.) Report final alpha values #
alpha
# [4,] 3.1898457
# [7,] 2.8103949
# [9,] 0.3794508
# rest are zero 


## Verify the solution ##

require(e1071)

# fit the SVM model
fit = svm(y ~ x1 + x2, data=data, scale = FALSE, 
          type = 'C', kernel = 'linear', cost = 1e100)

# get the index numbers for alpha
fit$index # 4,7,9

# get alpha values
fit$coefs / (-data$y[fit$index])
# [,1]
# [1,] 3.1909445
# [2,] 2.8117807
# [3,] 0.3791638

# compute bHat from SVM model 
t(fit$coefs) %*% fit$SV 
