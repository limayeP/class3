# Clear memory
rm(list=ls())
packages <- c("ggplot2",
              "nnet",
              "ggplot2",
              "caret",
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

b <- read.csv("codon_train_transformed.csv", header = TRUE)
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

te <- read.csv("codon_test_transformed.csv", header = TRUE)
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

# define training control using k-fold Cross Validation
TC <- trainControl(method="cv", 
                   number=10)
#summaryFunction = multiClassSummary)

#---------------------r-part--------------------
set.seed(123)
#train the regression tree model
pcartFit <- train(xbase,ybase,
                 trControl=TC,
                 method="rpart",
                 tuneLength = 10,
                 parms = list(split="information"))

print(pcartFit)
pcartFit$resample
# The variability in Accuracy is 4% acceptable
pcartFit$bestTune
plot(pcartFit, metric = "Accuracy")

plot(varImp(pcartFit))

# Apply the model to the test data
pmodel_cart <- predict(pcartFit, tex)
table(actual_y,pmodel_cart)

save(pmodel_cart, file = "CART_process_model.RData")
save(pcartFit, file = "CART_process_model_CV.RData")

#     0    1    2
# 0 1840    3   10
# 1   16  559    4
# 2   19    1  143
library(rattle)
fancyRpartPlot(pcartFit$finalModel)

#---------------------treebag--------------------
# define training control using k-fold Cross Validation
set.seed(123)
ptreebagFit <- train(DNAtype ~ .,
                      method     = "treebag",
                      trControl  = TC,
                      metric     = "Accuracy",
                      data       = bf)
ptreebagFit$resample
# 4%
plot(varImp(ptreebagFit))
# Apply the model to the test data
ptreebag_model <- predict(ptreebagFit, tex)
table(actual_y, ptreebag_model)
# actual_y    0    1    2
ptreebag_model

# ptreebag_model
# actual_y    0    1    2
# 0 1847    1    5
# 1    9  569    1
# 2   14    3  146
save(ptreebagFit, file = "treebag_process_model.RData")
#---------------------C5.0--------------------
# define training control using k-fold Cross Validation
pfitC5.0 <- train(DNAtype ~ .,
                 data = bf,
                 method = "C5.0",
                 trControl = TC)
# Examine the variability among our 10 folds.
pfitC5.0$resample
# 0.9961501 - 0.9894129
plot(pfitC5.0, metric = "Kappa")
plot(pfitC5.0, metric = "Accuracy")
pfitC5.0$bestTune
#  trials model winnow
#9     20 rules  FALSE
pfitC5.0$finalModel
# Number of boosting iterations: 20 
# Average number of rules: 62.2 


# Apply the model to the test data
pC5.0_model <- predict(pfitC5.0, tex)
table(actual_y,pC5.0_model)

# actual_y    0    1    2
# 0 1851    1    1
# 1    7  570    2
# 2    6    3  154

save(pfitC5.0, file = "model_C5_process.RData")
#save(rfFit, file = "model_RF_NO_process_CV.RData")

#---------------------KNN--------------------
# define training control using k-fold Cross Validation
set.seed(123)

Pmyknnmod = train(DNAtype ~ .,
                 data = bf,
                 method = "knn",
                 trControl = TC,
                 tuneGrid = expand.grid(k = seq(1, 13, by = 2))
)

# Examine the variability among our 10 folds.
Pmyknnmod$resample
ggplot(Pmyknnmod, highlight = TRUE)
Pmyknnmod$bestTune
# k = 1 is the best
Pmyknnmod$finalModel
# Apply the model to the test data
PKNN_model <- predict(Pmyknnmod, tex)
table(actual_y,PKNN_model)
# actual_y    0    1    2
# 0 1849    1    3
# 1    0  577    2
# 2    0    0  163
save(Pmyknnmod, file = "model_knn_process.RData")
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
Pplsda_fit <- train(DNAtype ~ .,
                   data = bf,
                   method = "pls",
                   tuneLength = 10, # number of components
                   trControl = TC,
                   metric = "Accuracy")
# Examine the variability among our 10 folds.
Pplsda_fit$resample
save(Pplsda_fit, file = "model_plsda_PROCESSED.RData")
plot(Pplsda_fit, metric = "Kappa")
plot(Pplsda_fit, metric = "Accuracy")
Pplsda_fit$finalModel
# Partial least squares classification, fitted with the 
# orthogonal scores algorithm.
# The softmax function was used to compute class probabilities.
plot(Pplsda_fit$finalModel$scores)
# Apply the model to the test data
Pplsda_model <- predict(Pplsda_fit, tex)
table(actual_y,Pplsda_model)
# actual_y    0    1    2
# 0 1851    1    1
# 1   18  556    5
# 2   32    1  130

#---------------------Multinomial Regression--------------------
# define training control using k-fold Cross Validation
set.seed(123)
pmyLoRfit <- train(
  DNAtype~., 
  data = bf,
  method = "multinom",
  trControl = TC,
  trace = FALSE
)
# Examine the variability among our 10 folds.
pmyLoRfit$resample
summary(pmyLoRfit)
ggplot(pmyLoRfit, highlight = TRUE)
pmyLoRfit$bestTune
# k = 5 is the best
pmyLoRfit$finalModel
save(pmyLoRfit, file = "multinomLR_model_porcessed.RData")
# Apply the model to the test data
pLoR_model <- predict(pmyLoRfit, tex)
table(actual_y,pLoR_model)
# actual_y    0    1    2
# 0 1851    1    1
# 1    1  577    1
# 2    2    1  160

#---------------------RF--------------------
# define training control using k-fold Cross Validation
parametersGrid <-  expand.grid(mtry = 1:64)
set.seed(123)
pRFfit <- train(DNAtype ~ .,
                 method     = "rf",
                 trControl  = TC,
                 metric     = "Accuracy",
                 tuneGrid   = parametersGrid,
                 data       = bf)
summary(pRFfit)
# Examine the variability among our 10 folds.
pRFfit$resample
ggplot(pRFfit, highlight = TRUE)
plot(varImp(pRFfit))
varImp(pRFfit)
pRFfit$bestTune
# mtry
# 18
pRFfit$finalModel
# Call:
#   randomForest(x = x, y = y, mtry = min(param$mtry, ncol(x))) 
# Type of random forest: classification
# Number of trees: 500
# No. of variables tried at each split: 18
# 
# OOB estimate of  error rate: 0.97%
# Confusion matrix:
#   0    1   2 class.error
# 0 7406    6   2  0.00107904
# 1   36 2276   8  0.01896552
# 2   42    7 604  0.07503828
save(pRFfit, file = "random_forests_mtry22_after_transformationmodel.RData")
# Apply the model to the test data
RF_model <- predict(pRFfit, tex)
table(actual_y,RF_model) 
# actual_y    0    1    2
# 0 1851    1    1
# 1    9  570    0
# 2    5    2  156

library(randomForest)
prf.model <- randomForest(DNAtype ~ ., 
                           data = bf, 
                           mtry = 10)


# Apply the model to the test data
prf <- predict(prf.model, tex)
table(actual_y,prf)
#---------------------NN--------------------
# define training control using k-fold Cross Validation
parametersGrid <-  expand.grid(size = 1:10,
                               decay = seq(from = 0.1, to = 0.5, by=0.1)
)
set.seed(123)
pnnfit <- train(DNAtype ~ .,
               method     = "nnet",
               trControl  = TC,
               metric     = "Accuracy",
               data       = bf,
               tuneGrid   = parametersGrid
)

print(pnnfit)
# Accuracy was used to select the optimal model using the largest value.
# The final values used for the model were size = 3 and decay = 0.1
# Examine the variability among our 10 folds.
pnnfit$resample
summary(pnnfit)
ggplot(pnnfit, highlight = TRUE)
pnnfit$bestTune
# k = 5 is the best
pnnfit$finalModel
# a 64-3-3 network with 207 weights
# inputs: UUU UUC UUA UUG CUU CUC CUA CUG AUU AUC AUA AUG GUU GUC GUA GUG GCU GCC GCA GCG CCU CCC CCA CCG UGG GGU GGC GGA GGG UCU UCC UCA UCG AGU AGC ACU ACC ACA ACG UAU UAC CAA CAG AAU AAC UGU UGC CAU CAC AAA AAG CGU CGC CGA CGG AGA AGG GAU GAC GAA GAG UAA UAG UGA 
# output(s): .outcome 
# options were - softmax modelling  decay=0.1
plot(pnnfit, metric = "Accuracy")
save(pnnfit, file = "nn_processed_model.RData")

# Apply the model to the test data
pnnmodel <- predict(pnnfit, tex)
table(actual_y,pnnmodel)