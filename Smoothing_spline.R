### STAT 5P87 LAB 3 ###

### Implementing an estimation method for  ###
### smoothing spline model ###

# Srushti Desai #


rm(list = ls())
setwd("/Users/srushtidesai/Desktop/School/COURSES/STAT 5P87/data_files")


# load piece-wise data
data = read.csv('piecewise-data.csv')

# x as input and y as output
X = data$X
Y = data$Y

# sort the data
o = order(X)
X = X[o]
Y = Y[o]
n = length(Y)


#-- Creating knots --#


xi = sort(X)
K = length(xi)


#-- Create basis expansion for knots at --#
#-- every point --#


N = matrix(0, nrow = n, ncol = K)

# Global terms 
N[,1] = 1
N[,2] = X

# "d" terms
d = matrix(0, nrow = n, ncol = K - 1)

for(k in c(1:(K-1))){
  d[,k] = (as.integer(X - xi[k] > 0) * (X - xi[k])^3 - as.integer(X - xi[K] > 0) * (X - xi[K])^3) / (xi[K] - xi[k])
}

for(k in c(1:(K-2))){
  N[, k + 2] = d[,k] - d[,K-1]
}

N


#-- Compute Omega using truncated polynomial --#
#-- basis for the natural spline --#


compute_omega = function(knots) {
  K = length(knots)
  Omega = matrix(0, K, K)
  for (i in 3:(K)) {
    for (j in 3:(K)) {
      if (i == j) {
        Omega[i,j] = 12 * (knots[K-1] - knots[i-2])^2 / (knots[K] - knots[i-2])
      } 
      else if(i > j) {
        Omega[i,j] = 6 * (knots[K-1] - knots[i-2]) * 
          ((3 * knots[j-2] - knots[K-1]) * knots[i-2] - knots[i-2]^2 + 
             2 * knots[K] * (knots[K-1] - knots[j-2]) - knots[j-2] * knots[K-1]) /
          ((knots[K] - knots[j-2]) * (knots[K] - knots[i-2]))
        Omega[j,i] = Omega[i,j]
      }     
    }
  }
  return(Omega)
}
 
Omega = compute_omega(xi)
Omega


# Part A #


#-- Estimate the smoothing spline model --##
#-- with lambda = 0.01 --##


lambda = 0.01  
bHat = solve(t(N) %*% N + lambda * Omega) %*% t(N) %*% Y
bHat


# Part B #


# Predict values 
x0 = seq(from=0, to=1, by=0.01)
n0 = length(x0)

# Construct basis for new points
N0 = matrix(0, nrow = length(x0), ncol = K)

N0[,1] = 1
N0[,2] = x0

d0 = matrix(0, nrow = length(x0), ncol = K - 1)

for(k in c(1:(K-1))){
  d0[,k] = (as.integer(x0 - xi[k] > 0) * (x0 - xi[k])^3 - 
              as.integer(x0 - xi[K] > 0) * (x0 - xi[K])^3) / (xi[K] - xi[k])
}

for(k in c(1:(K-2))){
  N0[, k + 2] = d0[,k] - d0[,K-1]
}

N0

# Compute yHat
yHat = N0 %*% bHat
yHat

#-- Plot for the problem with --#
#-- estimated regression line --#

plot(X, Y, 
     main="Smoothing Spline",
     col="forestgreen",
     xlab="X", ylab="Y", cex=0.8, pch=19, )
lines(x0, yHat, col="navyblue", lwd=2)


# Part C #


#-- Compute the corresponding df --#


df = sum(diag(N %*% solve(t(N) %*% N + lambda * Omega) %*% t(N))) 
df #3.3515


#-- Verify the solution --#


lambda = 0.01/(max(X) - min(X))^3
model = smooth.spline(X , Y , lambda = lambda )
model

df1 = model$df
df1 #3.3517
