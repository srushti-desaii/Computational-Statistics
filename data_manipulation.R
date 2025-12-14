rm(list = ls())
data = read.csv('a1-q1-data.csv')


# rename the variable 'x' to 'x1'

data$variable[data$variable == 'x'] = 'x1'
data


# remove all the rows corresponding to 
# observation 2 

data = data[data$observation != 2,] 
row.names(data) = NULL
data


# add rows to the data frame for a new  
# observation 4 (x1 = 3, y = 2)

data = rbind(data, data.frame(observation = 4, variable = c( "y" , "x1"), value = c(3, 2)))
row.names(data) = NULL
data


# add rows to the data frame for a new variable

data = rbind(data, data.frame(observation = c(1, 3, 4), variable = "x2", value = c(3, 1, 5)))
row.names(data) = NULL
data


# Reorder rows so observations are grouped together

data = data[order(data$observation), ]
row.names(data) = NULL
data


# Create a new column named ‘value-squared'

data$value_squared = data$value^2
data

# the final data frame has 9 rows (observation)
# and 4 columns with obs, variable , value 
# and squared value
# dim(data)
# [1] 9 4





# Create unique student identifiers across schools

rm(list = ls())

data = read.csv("a1-q2-data.csv")
head(data)


## Initialize an empty id list

id = list()

## Loop to compute new id's

for (i in 1:nrow(data)) {
  
  # get the school of current row
  school = data$school[i]  
  
  # calculate and assign an id if it doesn't exist
  if (!school %in% names(id)) {
    id[[school]] = max(data$student, na.rm = TRUE)
  }
  
  data$student[i] = data$student[i] + id[[school]]
}

data

# the data set will have unique ids for each student
# starting from 1 and ending on 62