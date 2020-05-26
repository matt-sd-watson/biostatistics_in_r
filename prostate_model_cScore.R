#prostate dataset
load("~/Documents/KU Leuven/Courses/Semester 2 2020/Statistical Methods for Bioinformatics/Assignment/prostate2.rdata")
# gain an impression of the predictor and response variables
head(prostate, 10)

# load some dependencies and useful libraries
library(plyr)
library(readr)
library(dplyr)
library(caret)
library(ggplot2)
library(repr)
library(glmnet)
library(corrplot)
library(dplyr)
library(mgcv)

# check if there are null values for any of the predictor variables and append all counts
# to a list
null_counts <- list()
for (i in colnames(prostate)) {
  null_counts[[i]] <- sum(is.na(prostate$i))
}

# we can see that there are no null values for any of the observations for the predictor variables
null_counts

# check to see if there are any impossible ages in the dataset for data semantics
impossible_age <- sum(prostate$age < 0)
impossible_age

# create a series of graphs to describe the data
names <- colnames(prostate)
names <- names[-1]
par(mfrow=c(2,4))
for (i in 1:7) {
  attach(prostate)
  x <- subset(prostate, select = names[[i]])
  final <- cbind(x, prostate$Cscore)
  colnames(final) <- c("x", "y")
  plot(final$x, final$y, xlab = names[[i]], ylab = "Cscore", col = "blue")
}
dev.off()

# create a correlation plot of the variables
res <- cor(prostate)
round(res, 2)
corrplot(res, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)

## creating the lasso model in r

# set the seed to establish partifioning for a training and validation set
# based on the descriptive visualizations, treat the svi variables as a binary factor

prostate$svi <- as.factor(prostate$svi)
set.seed(100)

# produce training and testing sets with 75% and 25% respectively
index = sample(1:nrow(prostate), 0.75*nrow(prostate))
training_set <- prostate[index,]
validation_set <- prostate[-index,]

dim(training_set)
dim(validation_set)

# create the dummy variables for the lasso regression

dummies <- dummyVars(Cscore ~ ., data = prostate)
train_dummies = predict(dummies, newdata = training_set)
test_dummies = predict(dummies, newdata = validation_set)

print(dim(train_dummies)); print(dim(test_dummies))

# establish a function to create a dataframe with the performance metrics used at testing
eval_results <- function(true, predicted, df) {
  SSE <- sum((predicted - true)^2)
  SST <- sum((true - mean(true))^2)
  R_square <- 1 - SSE / SST
  RMSE = sqrt(SSE/nrow(df))

  data.frame(
    RMSE = RMSE,
    Rsquare = R_square
  )
  
}

lambdas <- 10^seq(2, -3, by = -.1)

#establish the training objects to build the LASSO model
y_train = training_set %>%
  select(Cscore) %>%
  unlist() %>%
  as.numeric()
x = as.matrix(train_dummies)

# train the lasso model using the training set and find the optimal laambda value
lasso_reg <- cv.glmnet(x, y_train, alpha = 1, lambda = lambdas, standardize = TRUE, nfolds = 5)

# find the minimum of the lasso lambda values to apply as a penalty
lambda_best <- lasso_reg$lambda.min 
lambda_best

# evaluate the model on the test dataset
x_test <- as.matrix(test_dummies)
y_test <- validation_set %>%
  select(Cscore) %>%
  unlist() %>%
  as.numeric()
lasso_model <- glmnet(x_test, y_test, alpha = 1, lambda = lambda_best, standardize = TRUE)
predictions_test <- predict(lasso_model, s = lambda_best, newx = x_test)

# extract the test coefficients from the model building
test_coef <- predict(lasso_model, type = "coefficients", s = lambda_best)
test_coef

# observe the results of the model test using RMSE and R2 value
eval_results(y_test, predictions_test, validation_set)

## creating the optimal non linear model

# fit a gam model for each of the predictor variables in the dataset
# use a smooth spline for each of the variables except for the binary factor
# a lambda value of 1 is used based on experimentation to optimize the smoothing fit

gam_model <- gam(Cscore ~ s(lpsa, bs='cr', sp=1) + s(lcavol, bs='cr', sp=1) + s(lcp, bs='cr', sp=1) + s(lweight, bs='cr', sp=1) + svi, data = prostate)
summary(gam_model)
dev.off()
par(mfrow=c(2,3))
plot(gam_model)
