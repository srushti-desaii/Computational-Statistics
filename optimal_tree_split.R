# Find optimal regression tree split (feature + threshold) by MSE

# Inputs
# X: n X p input matrix
# Y: output vector

# Output 
# pair (j, h) where j indicates which input (column) is used in the split 
# and h the threshold


optimal_split = function(X,Y){
  n = dim(X)[1]
  p = dim(X)[2]
  
  # initialize matrix to store MSE values 
  # and computed threshold
  MSE = matrix(0, ncol = p, nrow = n-1)
  threshold = matrix(0, ncol = p, nrow = n-1)
  
  # loop through all inputs and observation
  for (i in 1:p){
    rownames(X) = c(1:nrow(X))
    Xnew = X[order(X[,i]),]
    Ynew = Y[as.numeric(rownames(Xnew))]
    
    # loop through mid points (for threshold)
    for (j in 1:(n-1)){
      
      # compute the threshold (as midpoint
      # b/w obs as potential thresholds)
      threshold[j,i] = (Xnew[j,i] + Xnew[(j+1),i]) / 2
      h0 = threshold[j,i]
      
      # split obs based on threshold
      Y1 = Ynew[Xnew[,i] < h0]
      Y2 = Ynew[Xnew[,i] >= h0]
      
      # compute the MSE
      MSE1 = sum((Y1 - mean(Y1))^2)
      MSE2 = sum((Y2 - mean(Y2))^2)
      
      MSE[j, i] = MSE1 + MSE2
    }
  }
  # Get index of minimum MSE
  index = which(MSE == min(MSE), arr.ind = TRUE)
  
  # Get the optimal pair
  h = threshold[index[1,1],index[1,2]]
  j = index[1,2]
  pair = c(j,h)
  
  return(pair)
}

# test the function 
X = matrix(c(2, 4, 6, 8, 10,  
             1, 3, 5, 7, 9),  
           ncol = 2, byrow = TRUE)
Y = c(5, 7, 9, 11, 13) 

optimal_split(X, Y)


