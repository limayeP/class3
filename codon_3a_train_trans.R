# https://topepo.github.io/caret/pre-processing.html
packages <- c("ggplot2",
              "rattle",
              "Hmisc",
              "psych",
              "caret",
              "e1071",
              "corrplot",
              "cowplot",
              "ggpubr",
              "DataExplorer")

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

d <- read.csv("~/Documents/CCSU/DATA_599/project_working/codon_train.csv")
View(d)
d$X <- NULL
# View type of each column
sapply(d, class)
# Convert DNAtype to a factor
d$DNAtype <- as.factor(d$DNAtype)
dim(d)
summary(d)
# ANy NAs
sum(is.na(d))
View(d)

#############################################################
# With DataExplorer understand the data
# https://boxuancui.github.io/DataExplorer/

create_report(d, y = "DNAtype")

# View basic description
introduce(d)
plot_intro(d)
plot_missing(d)

# Removing the target: DNAtype to get predictors only
adt = subset(d, select = -c(1))
View(adt)
names(adt)
# skewness values for codon frequencies determined with caret package
skewValues <- apply(adt, 2, skewness)
# This statistic can be used to prioritize variables for visualizing.
skewValues
plot_histogram(adt)
# All the codon frequencies were right skewed.

# correlations between predictors
cor(adt)
rdc <- cor.mtest(adt, conf.level = 0.95)

plot_qq(adt)
plot_qq(adt$AAU, title = "AAU before squareroot transformation")
plot_qq(adt$ACA, title = "ACA before squareroot transformation")
###########################################################################
TRANSFORMATION
# Apply square root transformation since the data have negative values
fun <- function(x) {
  sqrt(x)
}
rt <- apply(adt,2,fun)
rs <- data.frame(rt)
dna <- subset(d, select = c(1))
det <-cbind(dna, rs)
View(det)
dim(det)
# checking for symmetrical distribution after transformation
plot_histogram(rs)
plot_qq(rs)
plot_qq(rs$AAU, title = "AAU after squareroot transformation")
plot_qq(rs$ACA, title = "ACA after squareroot transformation")

##################################################################################
# The values are more or less uniformly distributed
skewValues_rs <- apply(rs, 2, skewness)
# This statistic can be used to prioritize variables for visualizing.
skewValues_rs
############################################################
# find the standard deviations of all the codons
sd <- apply(rs,2,sd)
##############################################################3
# Save train and test sets as separate csv files
write.csv(det, "codon_train_transformed.csv", row.names = TRUE)

###########################################################################
# Plot correlations between predictors in the training set
df5 <- read.csv("codon_train_transformed.csv", header = TRUE)
names(df5)
# Drop "X" column
# Get only predictors
df6<- subset(df5, select = -c(1))
names(df6)
# Observe the data type
sapply(df6, class)
sapply(df6, typeof)
# [1] 12982    65
df6$DNAtype <- as.factor(df6$DNAtype)
df7<- subset(df6, select = -c(1))
# ANy NAs
sum(is.na(df7))
summary(df7)
View(df7)
######################################################################
##############################################################
# TRANFORMING THE TEST DATA AND SCALING AND CENTERING IT
#######################################################################
a <- read.csv("~/Documents/CCSU/DATA_599/project_working/codon_test.csv")

View(a)
a$X <- NULL
# View type of each column
sapply(a, class)
# Convert DNAtype to a factor
a$DNAtype <- as.factor(a$DNAtype)
dim(a)
summary(a)
# ANy NAs
sum(is.na(a))
View(a)
# Removing the target: DNAtype to get predictors only
dt<- subset(a, select = -c(1))
# TRANSFORMATION
# Apply square root transformation since the data have negative values
fun <- function(x) {
  sqrt(x)
}
t <- apply(dt,2,fun)
s <- data.frame(t)
dna <- subset(a, select = c(1))
et <-cbind(dna, s)

###############################
# Save train and test sets as separate csv files
write.csv(et, "codon_test_transformed.csv", row.names = TRUE)

##########################################################################3
# Check if the dataset is imbalanced.
t1 <- table(d$DNAtype)
#   0    1    2 
# 9267 2899  816 

length(d$DNAtype)
# [1] 10387

# choose the value with which every row value is divided
someVar = 10387
# apply for each row, where row_val equals the values in rows
round(apply(t1, 1, function(row_val) row_val / someVar), 3)
#    0    1    2 
#   0.714 0.223 0.063 

# Check for class imbalance
TC <- trainControl(
  method = "CV",
  number = 10
)

fit.ch <- train(DNAtype ~.,
                data = d,
                method = "rpart",
                trControl = TC)

# Check for overfitting
fit.ch$resample
# Check on the test set
a <- read.csv("~/Documents/CCSU/DATA_599/project_working/codon_test.csv")
View(a)
a$X <- NULL
# View type of each column
sapply(a, class)
# Convert DNAtype to a factor
a$DNAtype <- as.factor(a$DNAtype)
dim(a)
summary(a)
testsetpreds.ch <- predict(fit.ch, a)
table(a$DNAtype, testsetpreds.ch)
testsetpreds.ch
#      0    1    2
# 0 1841    5    7
# 1   30  545    4
# 2   60    0  103


ads <- adt
########################################################################
#Transformations that did not work
#########################################################################
# Codon frequecies are right skewed. 
# Transformation applied. Did not work
FUN <- function(x){
  if(any(!(x > 0))){ 
    x + 10^-6
  }
  else x
}

dat<- sapply(ads, FUN)
dat <- data.frame(dat)
View(dat)

boxcox_skewness = function(x){
  BCT = BoxCoxTrans(x)
  x_bc = predict( BCT, x )
  skewness(x_bc) 
}

datty <- apply(dat, 2, boxcox_skewness)
class(datty)

trans<- preProcess(dat, method = c("BoxCox", "center", "scale"))
df <- predict(trans, dat)
# Plotting to see if uniformity has been induced.
## View histogram of continuous variables
plot_histogram(df)
plot_qq(df)
# Worsks well
###########******************************######################
###########*################################8########################3
# Codon frequecies are right skewed. 
# Transformation applied. Did not work
by<- preProcess(df, method = c("center", "YeoJohnson"))
vy <- predict(by, df)
plot_histogram(vy)
# Did not see uniformity
# Transformation applied. Did not work
ey<- preProcess(df, method = c("center","expoTrans" ))
gy <- predict(ey, df)
plot_histogram(gy)
# Did not see uniformity
##############################################################
#******************************************************************#



