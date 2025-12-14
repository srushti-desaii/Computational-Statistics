### STAT 5P87 LAB 1 ###

### Implementing Least Angle Regression (LAR) algorithm to solve LASSO paths ###

# Srushti Desai #


rm(list = ls())
setwd("/Users/srushtidesai/Desktop/School/COURSES/STAT 5P87/data_files")

# loading the prostate cancer data:
data = read.csv("prostate-data.csv")

# extracting the variables we are interested in:
data = data[,c(1, 2, 3, 9)]
head(data)

# lpsa as our output. 
y = data$lpsa

# lcavol , lweight  and age as our inputs and intercept as 0. 
X = model.matrix(lpsa ~ 0 + .  , data = data)
head(X)



### STEP 1 ###
#-- Scaling the data --# 

# Scaling the input data
xbar = apply(X , 2 , mean)
xbar
# lcavol   lweight       age 
# 1.350010  3.628943 63.865979 

# scaling the output
y = y - mean(y)

# scaling the inputs
Z = t((t(X) - xbar) / sqrt(colSums(t(t(X) - xbar)^2)))


# verify our scaling by scaling just first column
z = (X[,1] - xbar[1]) / sqrt(sum((X[,1]- xbar[1])^2)) 
identical(Z[,1],z) #TRUE

sd(Z[,1]) #0.1020621
round(mean(Z[,1])) #0 



### STEP 2 ###
#-- Applying LAR algorithm --# 

max_iterations = 4  
n = nrow(Z)         
p = ncol(Z)        


##  starting with bHat = 0

bHat = matrix(0, nrow = p, ncol = 1) 
new_bHat = matrix(0, nrow = p, ncol = 1) 


## in order to compute how much we want to increase b_j value

# creating matrices to store our penalty and beta values
t = rep(NA, max_iterations)                      
b = matrix(NA, nrow = max_iterations, ncol = p) 
b[1, ] = 0                                      
t[1] = 0                                         
J = c() #to store the current active set


## starting the loop 

for (iter in 2:max_iterations) {
  
  ## compute the correlations
  
  yHat = Z %*% bHat
  c = t(Z) %*% (y - yHat)  
  
  #maximum absolute correlation out of three
  C = max(abs(c))  
  
  #identify the inputs   that are maximally       correlated
  j = which(abs(c) == C) 
  J = sort(unique(c(j, J)))     
  
  ## create a reduced input matrix
  
  x_j = as.matrix(Z[, J])  
  
  # compute the adjusted matrix
  X_J = matrix(NA, nrow = n, ncol = length(J))
  for (i in 1:n) {
    for (j in 1:length(J)) {
      X_J[i, j] = x_j[i, j] * sign(c[J[j]]) 
    }
  }
  
  ## compute all the intermediate quantities required
  
  G_J = t(X_J) %*% X_J
  one = matrix(1, nrow = length(J), ncol = 1)
  A_J = sqrt(as.numeric(solve(t(one) %*% solve(G_J)   %*% one)))  
  w_J = as.numeric(A_J) * solve(G_J) %*% one 
  u_J = X_J %*% w_J           
  a = t(Z) %*% u_J                               
  
  ## identify the next most correlated variable
  
  gamma = numeric(0)
  for (j in 1:p) { 
    gamma1 = (C - c[j]) / (A_J - a[j])
    gamma2 = (C + c[j]) / (A_J + a[j])
    gamma1 = gamma1[gamma1 > 0]
    gamma2 = gamma2[gamma2 > 0]
    gamma = c(gamma, gamma1, gamma2)
  }
  gamma = gamma[is.finite(gamma) & gamma > 0]
  gammaHat = min(gamma)  #smallest positive gamma
  
  ## Update beta_j 
  
  for (j in 1:length(J)) { 
    new_bHat[J[j], 1] = bHat[J[j], 1] + sign(c[J[j]]) * gammaHat * w_J[j]
  }
  
  bHat = new_bHat
  b[iter, ] = t(new_bHat)
  t[iter] = sum(abs(new_bHat))  
}

## print the final outputs

b #all the beta values

#         [,1]     [,2]       [,3]
# [1,] 0.000000 0.000000  0.0000000
# [2,] 4.733879 0.000000  0.0000000
# [3,] 7.121754 2.387874  0.0000000
# [4,] 7.643797 3.057895 -0.8657997

t #all the penalty 

# [1]  0.000000  4.733879  9.509628
# [4] 11.567492



### STEP 3 ###
#-- Create LASSO path --# 

## plotting the sequence of bHat v/s the t values to create LASSO path

plot(
  t, b[, 1], type = "b", 
  xlab = "t", ylab = "Coefficient", 
  col = "indianred", lwd = 1, pch = 1, 
  ylim = range(b, na.rm = TRUE), 
  main = "LASSO Path"
)

lines(t, b[, 2], type = "b", col = "navyblue", lwd = 1, pch = 1)
lines(t, b[, 3], type = "b", col = "forestgreen", lwd = 1, pch = 1)

legend(
"topleft", legend = colnames(X), col = c("indianred", "navyblue", "forestgreen"), lty = 1, lwd = 1
)

