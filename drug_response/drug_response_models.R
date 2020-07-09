library(ggplot2)
library(pls)
library(ggbiplot)
library(pROC)
library(xgboost)
library(caret)
library(nnet)
library(Matrix)

# read the data into variables
responses <- read.table(file.choose())
mutations <- read.table(file.choose())

# look at the header to gain an impression of the datasets
head(responses, 10)
head(mutations, 10)

# perform some intial data cleaning and plot the 
# frequency of the lnIC50 values before establishing thresholds
responses_plot <- as.data.frame(responses[2:101,])
colnames(responses_plot) <- c("sample", "lnic50")
responses_plot$lnic50 <- as.numeric(as.character(responses_plot$lnic50))
head(responses_plot, 10)
hist <- hist(responses_plot$lnic50,col = "blue")
text(hist$mids,hist$counts,labels=hist$counts, adj=c(0.5, -0.5))

# observe the min and max values of the response variable
min(responses_plot$lnic50)
max(responses_plot$lnic50)

# create a plot of lnIC50 values for each sample
# include a geom line of the mean for visual comparison
plot <- ggplot(responses_plot, aes(x = sample, y = lnic50)) + geom_point() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
plot + geom_hline(yintercept = mean(responses_plot$lnic50), color = "blue", size = 2)

# create the PCA object for the mutations
mutations_for_pca <- t(as.data.frame(mutations))
mutations_pca<- prcomp(mutations_for_pca, center = TRUE,scale. = TRUE)
mutations_pca

# prepare the merged dataset to link mutations to the reponse
binary_mutations <- mutations_for_pca[ order(row.names(mutations_for_pca)), ]
merged_data <- as.data.frame(cbind(binary_mutations, responses_plot$lnic50))
multi_class_data <- merged_data
multi_class_data$V51 <- as.numeric(as.character(multi_class_data$V51))

# create the threholds for the classifications
# below 0: sensitive
# above 4.5: resistant
# between: neither
multi_class_data$class <- ifelse(multi_class_data$V51 >= 4.5, 2, ifelse(multi_class_data$V51 <= 0, 0, 1)) 
multi_class_data
multi_class_data <- multi_class_data[ -c(51) ]

# create letter labels for the PCA plot
multi_class_data$groups <- ifelse(multi_class_data$class == 0, "sensitive", ifelse(multi_class_data$class == 2, "resistant", "neither"))

# create PCA plot with cell classes
pca_model <- ggbiplot(mutations_pca, groups = multi_class_data$groups, ellipse = TRUE, circle = TRUE) + ggtitle("PCA, resistant and sensitive")
pca_model

# using the first 2 PCA components, identify
# a preliminary set of genes that may be involved in cell responses
components <- as.data.frame(mutations_pca$rotation)
genes_for_sensitivity <- subset(components, PC1 > 0.25)
genes_for_sensitivity <- row.names(genes_for_sensitivity)
genes_for_sensitivity

genes_for_resistance <- subset(components, PC2 < -0.15)
genes_for_resistance <- row.names(genes_for_resistance)
genes_for_resistance

multi_class_data <- multi_class_data[ -c(52) ]

# create the training and test datasets
set.seed(99)
training_size = sample(1:nrow(multi_class_data), 0.8*nrow(multi_class_data))
training_set <- multi_class_data[training_size,]
test_set <- multi_class_data[-training_size,]
nrow(training_set)
nrow(test_set)

# view the test dataset to observe the balance of cell line classes
test_set

# create the matrices required for the xg boost 

# xg boost matrices for training set
train_labels <- as.vector(training_set$class)
train_sparse <- sparse.model.matrix(class ~ ., data = training_set)
multi_train_matrix <- xgb.DMatrix(data=train_sparse, label = train_labels)

# xg boost matrices for the test dataset
test_labels <- as.vector(test_set$class)
test_sparse <- sparse.model.matrix(class ~ ., data = test_set)
multi_test_matrix <- xgb.DMatrix(data=test_sparse, label = test_labels)

# make a variable for the number of classes for cell lines
num_class = length(levels(as.factor(multi_class_data$class)))

# set the parameters for boosting model
params <- list(
  booster="gbtree",
  eta=0.001,
  max_depth=4,
  gamma=3,
  subsample=0.75,
  colsample_bytree=1,
  objective="multi:softprob",
  eval_metric="mlogloss",
  num_class=num_class
)

# create the boosting model
multi_boost <- xgb.train(
  params=params,
  data=multi_train_matrix,
  nrounds=10000,
  nthreads=1,
  early_stopping_rounds=10,
  watchlist=list(val1=multi_train_matrix, val2=multi_test_matrix),
  verbose=0
)

# view the boosting model
multi_boost

#make predictions for the test dataset 
multi_boost_predict <- predict(multi_boost, multi_test_matrix, reshape=T)
multi_boost_predict <- as.data.frame(multi_boost_predict)
colnames(multi_boost_predict) = levels(multi_class_data$class)

# create the confusion matrix to assess the model accuracy
confusion <- confusionMatrix(as.factor(multi_boost_predict), as.factor(test_set$class))
confusion

# retrieve the important genes from xg boost
mat <- xgb.importance (feature_names = colnames(multi_train_matrix), model = multi_boost)
xgb.plot.importance (importance_matrix = mat[1:10])

# we apply some of the genes retrieved from the multinomial boost to 
# another multinomial classification algorithm
# Gene5 is included because it improves performance
multi_model <- nnet::multinom(class ~ Gene3 + Gene5 + Gene12 + Gene9 + Gene36 + Gene37, data = training_set)
summary(multi_model)

# create predictions and confuson matrices for
# the second multinomial classifier
predicted.classes <- multi_model %>% predict(test_set)
confusion <- confusionMatrix(as.factor(predicted.classes), as.factor(test_set$class))
confusion

# view the ROC curves for the training and test datasets
multi_predict_train <- predict(multi_model, newdata = training_set)
test_roc <- multiclass.roc(as.factor(training_set$class) ~ as.numeric(multi_predict_train), plot = TRUE, print.auc = TRUE)

multi_predict  <- predict(multi_model, newdata = test_set)
test_roc <- multiclass.roc(as.factor(test_set$class) ~ as.numeric(multi_predict), plot = TRUE, print.auc = TRUE)
