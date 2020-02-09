

#please watch some short R programming instructional videos and follow along in R while you watch. 

https://www.youtube.com/playlist?list=PLqzoL9-eJTNBDdKgJgJzaQcY6OXmsXAHU


#R data type

# Declare variables of different types
# Numeric
x <- 28
class(x)

# String
y <- "R is Fantastic"
class(y)

# Boolean
z <- TRUE
class(z)

# First way to declare a variable:  use the `<-`
name_of_variable <- value
# Second way to declare a variable:  use the `=`
name_of_variable = value

# Print variable x
x <- 42
x

y  <- 10
y

x-y

# Numerical vector
vec_num <- c(1, 10, 49)
vec_num

# Create the vectors
vect_1 <- c(1, 3, 5)
vect_2 <- c(2, 4, 6)
# Take the sum of A_vector and B_vector
sum_vect <- vect_1 + vect_2
# Print out total_vector
sum_vect
product_vect <- vect_1*vect_2
# Print out product_vector
product_vect

# Slice the first five rows of the vector
slice_vector <- c(1,2,3,4,5,6,7,8,9,10)
slice_vector[1:5]

# Create a vector from 1 to 10
logical_vector <- c(1:10)
logical_vector>5

# Print value strictly above 5
logical_vector[(logical_vector>5)]

# Print 5 and 6
logical_vector <- c(1:10)
logical_vector[(logical_vector>4) & (logical_vector<7)]


# Matrix operation:

# Construct a matrix with 5 rows that contain the numbers 1 up to 10 and byrow =  TRUE 

matrix_a <-matrix(1:10, byrow = TRUE, nrow = 5)

matrix_a

# Print dimension of the matrix with dim()
dim(matrix_a)

# Construct a matrix with 5 rows that contain the numbers 1 up to 10 and byrow =  FALSE
matrix_b <-matrix(1:10, byrow = FALSE, nrow = 5)
matrix_b

# concatenate c(1:5) to the matrix_a
matrix_a1 <- cbind(matrix_a, c(1:5))
# Check the dimension
dim(matrix_a1)

matrix_a1

matrix_a2 <-matrix(13:24, byrow = FALSE, ncol = 3)

matrix_c <-matrix(1:12, byrow = FALSE, ncol = 3)		
matrix_d <- cbind(matrix_a2, matrix_c)
dim(matrix_d)

matrix_c <-matrix(1:12, byrow = FALSE, ncol = 3)
# Create a vector of 3 columns
add_row <- c(1:3)
# Append to the matrix
matrix_c <- rbind(matrix_c, add_row)
# Check the dimension
dim(matrix_c)

#matrix_c[1,2] selects the element at the first row and second column.
#matrix_c[1:3,2:3] results in a matrix with the data on the rows 1, 2, 3 and columns 2, 3, 
#matrix_c[,1] selects all elements of the first column.
#matrix_c[1,] selects all elements of the first row.



### Generate Uniform random variables

x <- runif(1000, 0, 5)
hist(x)


### Generate Bernuolli random variables
install.packages('Rlab')
library(Rlab)
x <- rbern(1000, 0.3) # 1000 is sample size and 0.3 is the probability of taking value 1
hist(x)

### Generate Exponential random variables

x <- rexp(1000,5) # 1000 is sample size and 5 is the exponential parameter
hist(x)

### Generate Poisson random variables

x <- rpois(1000,5) # 1000 is sample size and 5 is the Poisson parameter
hist(x)

### Generate Bionomial random variables

x <- rbinom(1000,200, 0.5) #1000 is sample size, 200 is total number of trials and 0.5 is probability of taking value 1 for each trial.

hist(x)


### Generate Normal random variables

x<- rnorm(1000,0, 2) #1000 is sample size, 0 is mean, 2 is the standard deviation 
hist(x)

#### calculate Normal probability

#P(X < 2) 
#P(5>X>2) = P(X<5)- P(X<2)


pnorm(5, 0,2) - pnorm(2, 0, 2)
pnorm(2, 0, 2)


### Central Limit Theorem

n <-5000 
# your choice
CLT_results <- c() 
#initialize a vector
for(i in 1:n){ 
first_toss <- sample(1:6, size=1)
second_toss <- sample(1:6, size=1)
result_i <- first_toss + second_toss
CLT_results <- c(CLT_results, result_i) 
# append to vector 
}
hist(CLT_results, xlab="") 
#visualize 