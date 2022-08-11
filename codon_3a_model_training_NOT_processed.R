# Clear memory
rm(list=ls())
packages <- c("ggplot2",
              "nnet",
              "ggplot2",
              "caret",
              "randomForest",
              "C50"
              )

package.check <- lapply(
  packages,
  FUN = function(x)
  {
    if (!require(x, character.only = TRUE))
    {
      install.packages(x, dependencies = TRUE,
                       repos = "http://cran.us.r-project.org")
      library(x, character.only = TRUE)
    }
  }
)

#  NO PROCESSING rpart (CART MODEL)
b <- read.csv("codon_train.csv", header = TRUE)
names(b)
b$DNAtype <- factor(b$DNAtype)
# Drop "X" column
# Get only predictors
bf<- subset(b, select = -c(1))
# Observe the data type
sapply(bf, class)

# include only the predictors
xbase<- subset(bf, select = -c(1))
ybase <- bf$DNAtype

names(xbase)
class(ybase)

te <- read.csv("codon_test.csv", header = TRUE)
names(te)
te$DNAtype <- factor(te$DNAtype)

# Drop "X" column
dp<- subset(te, select = -c(1))

# Observe the data type
sapply(dp, class)

# include only the predictors
tex<- subset(dp, select = -c(1))
actual_y <- dp$DNAtype
names(tex)
class(actual_y)

TC <- trainControl(method="cv", 
                   number=10)

#---------------------r-part--------------------
# define training control using k-fold Cross Validation
#train the regression tree model
set.seed(123)
cartFit <- train(xbase,ybase,
                 trControl=TC,
                 method="rpart",
                 tuneLength = 5,
                 parms = list(split="information"))

print(cartFit)
# Accuracy was used to select the optimal model using the largest value.
# The final value used for the model was cp = 0.007063572.
cartFit$resample

cartFit$bestTune
plot(cartFit, metric = "Accuracy")
plot(varImp(cartFit))
# Apply the model to the test data
model_cart <- predict(cartFit, tex)
table(actual_y,model_cart)
save(model_cart, file = "CART_N0_process_model.RData")
save(cartFit, file = "CART_NO_process_model_CV.RData")


# actual_y    0    1    2
# 0 1842    5    6
# 1   27  545    7
# 2   49    0  114
library(rattle)
fancyRpartPlot(cartFit$finalModel)
#---------------------treebag--------------------
# define training control using k-fold Cross Validation
set.seed(123)
treebagFit <- train(DNAtype ~ .,
                      method     = "treebag",
                      trControl  = TC,
                      metric     = "Accuracy",
                      data       = bf)
treebagFit$resample
# 0.9932628 - 0.9826756
# 1%

# Apply the model to the test data
treebag_model <- predict(treebagFit, tex)
table(actual_y, treebag_model)
# 0    1    2
# 0 1847    1    5
# 1    9  569    1
# 2   13    3  147
save(treebagFit, file = "treebag_N0_process_model.RData")

#---------------------C5.0--------------------
# define training control using k-fold Cross Validation
fitC5.0 <- train(DNAtype ~ .,
                 data = bf,
                 method = "C5.0",
                 trControl = TC)
# Examine the variability among our 10 folds.
fitC5.0$resample
# 0.9961501 - 0.9894129
plot(fitC5.0, metric = "Kappa")
plot(fitC5.0, metric = "Accuracy")
fitC5.0$bestTune
fitC5.0$finalModel
# Number of boosting iterations: 20 
# Average number of rules: 30.4 

# Apply the model to the test data
C5.0_model <- predict(fitC5.0, tex)
table(actual_y,C5.0_model)

#    0    1    2
# 0 1852    0    1
# 1    5  573    1
# 2    4    3  156
save(fitC5.0, file = "model_C5_N0_process.RData")
#save(rfFit, file = "model_RF_NO_process_CV.RData")

#---------------------KNN--------------------
# define training control using k-fold Cross Validation
set.seed(123)

myknnmod = train(DNAtype ~ .,
                 data = bf,
                 method = "knn",
                 trControl = TC,
                 tuneGrid = expand.grid(k = seq(1, 13, by = 2))
)

ggplot(myknnmod, highlight = TRUE)
myknnmod$bestTune
# k = 1 is the best
myknnmod$finalModel
# Apply the model to the test data
KNN_model <- predict(myknnmod, tex)
table(actual_y,KNN_model)
# 0    1    2
# 0 1845    2    6
# 1    1  577    1
# 2    1    1  161
save(myknnmod, file = "model_knn_N0_process.RData")
#---------------------XGboost--------------------
# define training control using k-fold Cross Validation
set.seed(123)
parametersGrid <-  expand.grid(eta = 0.1, 
                               colsample_bytree=c(0.5,0.7),
                               max_depth=c(3,6,9, 13),
                               nrounds=100,
                               gamma=1,
                               min_child_weight=2, 
                               subsample=0.5
)

xgbFit <- train(DNAtype ~ .,
               method     = "xgbTree",
               trControl  = TC,
               metric     = "Accuracy",
               data       = bf, 
               tuneGrid   = parametersGrid
)

print(xgbFit)
xgbFit$resample
save(xgbFit, file = "xgb13_model.RData")

ggplot(xgbFt, highlight = TRUE)
xgbFit$bestTune
# k = 1 is the best
xgbFit$finalModel
# Apply the model to the test data
xgb_model <- predict(xgbFit, tex)
table(actual_y,xgb_model)
#---------------------plsDA--------------------
# define training control using k-fold Cross Validation
set.seed(123)
plsda_fit <- train(DNAtype ~ .,
                   data = bf,
                   method = "pls",
                   tuneLength = 10, # number of components
                   trControl = TC,
                   metric = "Accuracy")
# Examine the variability among our 10 folds.
plsda_fit$resample
save(plsda_fit, file = "model_plsda.RData")
plot(plsda_fit, metric = "Kappa")
plot(plsda_fit, metric = "Accuracy")
plsda_fit$finalModel
# Partial least squares classification, fitted with the 
# orthogonal scores algorithm.
# The softmax function was used to compute class probabilities.
plot(plsda_fit$finalModel$scores)
# Apply the model to the test data
plsda_model <- predict(plsda_fit, tex)
table(actual_y,plsda_model)
#       0    1    2
# 0 1852    0    1
# 1   14  561    4
# 2   56    1  106
#---------------------Multinomial Regression--------------------
# define training control using k-fold Cross Validation
set.seed(123)
myLoRfit <- train(
  DNAtype~., 
  data = bf,
  method = "multinom",
  trControl = TC,
  trace = FALSE
)
# Examine the variability among our 10 folds.
myLoRfit$resample
summary(myLoRfit)
ggplot(myLoRfit, highlight = TRUE)
myLoRfit$bestTune
# k = 5 is the best
myLoRfit$finalModel
save(myLoRfit, file = "multinomLR_model.RData")
# Apply the model to the test data
LoR_model <- predict(myLoRfit, tex)
table(actual_y,LoR_model)
######################################################################
#---------------------RF--------------------
# define training control using k-fold Cross Validation
parametersGrid <-  expand.grid(mtry = 1:64)
set.seed(123)
npRFfit <- train(DNAtype ~ .,
             method     = "rf",
             trControl  = TC,
             metric     = "Accuracy",
             tuneGrid   = parametersGrid,
             data       = bf)
summary(npRFfit)
# Examine the variability among our 10 folds.
npRFfit$resample
ggplot(npRFfit, highlight = TRUE)
plot(varImp(npRFfit))
npRFfit$bestTune
# mtry
# 22   22
npRFfit$finalModel
# randomForest(x = x, y = y, mtry = min(param$mtry, ncol(x))) 
# Type of random forest: classification
# Number of trees: 500
# No. of variables tried at each split: 22
# 
# OOB estimate of  error rate: 0.91%
# Confusion matrix:
#   0    1   2 class.error
# 0 7404    7   3  0.00134880
# 1   35 2279   6  0.01767241
# 2   37    7 609  0.06738132
save(npRFfit, file = "random_forests_model.RData")
# Apply the model to the test data
RF_model <- predict(npRFfit, tex)
table(actual_y,RF_model) 

#       0    1    2
# 0 1851    1    1
# 1    9  569    1
# 2    5    2  156
save(RF_model, file = "random_forests_mtry22_model.RData")
nprf.model <- randomForest(DNAtype ~ ., 
                         data = bf, 
                         mtry = 10)

# Apply the model to the test data
nprf <- predict(nprf.model, tex)
table(actual_y,nprf)

#    0    1    2
# 0 1852    0    1
# 1   10  569    0
# 2    6    2  155
###########################################################################
set.seed(123)
parametersGrid <-  expand.grid(size = 1:10,
                               decay = seq(from = 0.1, to = 0.5, by=0.1)
)

nnfit <- train(DNAtype ~ .,
             method     = "nnet",
             trControl  = TC,
             metric     = "Accuracy",
             data       = bf,
             tuneGrid   = parametersGrid
)

print(nnfit)
summary(nnfit)
ggplot(nnfit, highlight = TRUE)
nnfit$bestTune
# k = 5 is the best
nnfit$finalModel
plot(nnfit, metric = "Accuracy")
save(nnfit, file = "nn_not_processed_model.RData")
# Apply the model to the test data
nnmodel <- predict(nnfit, tex)
table(actual_y,nnmodel)
# actual_y    0    1      2
# 0         1851  1       1
# 1         7     567     5
# 2         3    0      160