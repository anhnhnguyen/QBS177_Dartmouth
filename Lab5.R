# Use the subset of the data to perform stepwise regression
data <- read.table("fulldat.txt", header=F)
data[1:3,1:10]
y <- data[,1]
x <- data[,3:50]
workingD <- data.frame(cbind(y,x))

# Random divide the data into training and testing
set.seed(1)
train <- sample(1:240,160,replace=FALSE) # randomly divide the data in to training and testing
test <- (1:240)[-train]

# Perform stepwise regression on training set
fit1 <- lm(y ~ ., workingD, subset=train)
fit2 <- lm(y ~ 1, workingD, subset=train)
step1 <- step(fit1,direction="both")  #both direction
step2 <- step(fit1,direction="backward") #backward 
step3 <- step(fit2,direction="forward",scope=list(upper=fit1,lower=fit2))  #forward

summary(step1)
summary(step2)
summary(step3)

# Check your training model on testing set
tempname <- names(step1$coef) 
colname <- match(tempname[-1], names(x)) #Model identification
predscore <- as.matrix(x[test,colname])%*%step1$coef[-1] #testing score
summary(lm(y[test]~predscore))

# Use the full data to run penalized regression
data <- read.table("fulldat.txt", header=F)
data<- as.matrix(data)
y <- data[,1]
x <- data[,3:7401]
workingD <- data.frame(cbind(y,x))

# Divide the data into training and testing
set.seed(1)
train <- sample(1:240,160,replace=FALSE) # randomly divide the data in to training and testing
test <- (1:240)[-train]
ytrain <- y[train]
ytest <- y[test]
xtrain <- x[train,]
xtest <- x[test,]

# Perform univariate regression on all genes on training set
pvalue <- 0
for (i in 1:7399){
  fit <- lm(ytrain~xtrain[,i])
  pvalue[i] <- summary(fit)$coef[2,4] 
}

# Perform penalized regression on training set
install.packages("glmnet", repos = "http://cran.us.r-project.org")
library(glmnet)
cvfit = cv.glmnet(xtrain, ytrain) # lasso regression
cvfit = cv.glmnet(xtrain, ytrain, alpha=0) # ridge regression
cvfit = cv.glmnet(xtrain, ytrain, alpha=0.5) #elastic net

plot(cvfit)
cvfit$lambda.min

# Check your training model on testing set
coef.min = coef(cvfit, s = "lambda.min")
active.min = which(coef.min != 0)
index.min = coef.min[active.min]

predscore <- xtest%*%coef.min[-1] #testing score

summary(lm(ytest~predscore))

















