# Fit polynomial-kernel SVM with stratified 5-fold CV

# Function to generate folds #
make_folds = function(Y, nFolds, stratified = FALSE, seed = 0){
  # K-Fold cross validation
  # Input:
  #   Y (either sample size, or vector of outputs)
  #   stratified (boolean): whether the folds should 
  #     be stratified. If TRUE then Y should be a vector of outputs
  # Output: list of vectors of fold indices
  set.seed(seed)
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


library(e1071)


# load the SA heart data 
data = read.csv('SAheart-data.csv')

# Convert chd = 0 to -1 
data$chd[data$chd == 0] = -1
n = nrow(data)

nFolds = 5
folds = make_folds(data$chd, nFolds, stratified = TRUE, seed = 0)

# Range of degrees to test
d_values = 1:10
n_d_values = length(d_values)

# Initialize accuracy matrix 
accuracy = matrix(0, nrow = n_d_values, ncol = nFolds)

# Cross-validation loop to get optimal degree
for (fold in 1:nFolds) {
  
  # Split the data into training and
  # testing
  training_data = data[-folds[[fold]], ]
  testing_data = data[folds[[fold]], ]
  
  for (i in 1:n_d_values) {
    d = d_values[i]
    
    # Fit SVM with polynomial kernel and 
    # gamma = 1 to predict chd using all 
    # inputs
    
    # Please note: cost is set to default
    # i.e 1
    fit = svm(chd ~ ., data = training_data,
              scale = TRUE, 
              type = "C-classification",
              kernel = "polynomial", 
              degree = d, gamma = 1)
    
    # Compute yHat
    yHat = predict(fit, newdata = testing_data)
    
    # Compute yhr accuracy
    accuracy[i, fold] = mean(yHat == testing_data$chd)
  }
}


accuracy = apply(accuracy, 1, mean)


# get which degree gives maximum accuracy
max(accuracy)
which.max(accuracy)


## Degree 1 gives the maximum accuracy with
## 72.7%


# fit the final model with degree 1
final_fit = svm(chd ~ ., data = data,
                scale = TRUE, 
                type = "C-classification",
                kernel = "polynomial",
                degree = 1, gamma = 1)
final_fit
