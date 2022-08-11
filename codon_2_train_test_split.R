# Train test split and validation of the split
# Load packages
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

d <- read.csv("codon_usage_c_imputed_nuc_mito_chloro.csv", header = TRUE)
names(d)
# Drop "X" column
# Get only predictors
df<- subset(d, select = -c(1))
names(df)
# Observe the data type
sapply(df, class)
# [1] 12982    65
df$DNAtype <- as.factor(df$DNAtype)
# ANy NAs
sum(is.na(df))
summary(df)
View(df)
##################################################################################################
# Create the train test sets and save them
set.seed(107)

inTrain <- createDataPartition(
  y = df$DNAtype, p = 0.80,
  list = FALSE)

# vector containing row indices of the records
# destined to be the training set
df.train <- df[ inTrain,]
df.test <- df[-inTrain,]

summary(df.train$DNAtype)
#   0    1    2 
# 7414 2320  653
# 0
# 7414/10387
# 0.7137768

# 1
# 2320/10387
# 0.2233561

# 2
# 653/10387
# 0.06286705

summary(df.test$DNAtype)
#    0    1    2 
# 1853  579  163 
# 2595

# 0
# 1853/2595
# 0.7140655

# 1
# 579/2595
# 0.2231214

# 2
# 163/2595
# 0.0628131

# validating class 0 partition
# H0: Ptraining = Ptest
# H1: Ptraining not equal to Ptest
# Where the P's represent the population proportions of
# records with class0= True
pt0 <- matrix(c(7414, 2973,1853, 742), nrow=2)
colnames(pt0) <- c("False", "True")
rownames(pt0) <- c("Training", "Test")
pt0
round(prop.table(pt0, 1), 4)
chisq.test(pt0, correct=FALSE)
# Since the value is 0.9768 is not less than 0.05
# do not reject H0

# validating class 1 partition
# H0: Ptraining = Ptest
# H1: Ptraining not equal to Ptest
# Where the P's represent the population proportions of
# records with class1= True
pt1 <- matrix(c(2320, 8067, 579, 2016), nrow=2)
colnames(pt1) <- c("False", "True")
rownames(pt1) <- c("Training", "Test")
pt1
round(prop.table(pt1, 1), 4)
chisq.test(pt1, correct=FALSE)
# Since the value is 0.9795 is not less than 0.05
# do not reject H0

# validating class 2 partition
# H0: Ptraining = Ptest
# H1: Ptraining not equal to Ptest
# Where the P's represent the population proportions of
# records with class2= True
pt2 <- matrix(c(653, 9734, 163, 2432), nrow=2)
colnames(pt2) <- c("False", "True")
rownames(pt2) <- c("Training", "Test")
pt2
round(prop.table(pt2, 1), 4)
chisq.test(pt2, correct=FALSE)
# Since the value is 0.9919 is not less than 0.05
# do not reject H0
#############################################
# The Kruskal-Wallis test is a non-parametric
# alternative to the analysis of variance or 2 sample T test
# It has less stringent assumptions as follows:
# 1) The groups are independent
# 2) The distributions are similar
# The assumption that the groups are independent is good 
# since the records were randomly assigned to each group
# check if the distributions are similar
# with boxplot by group
# NEED TO MERGE TRAIN AND TEST SET
df.train$part <- 
  rep("train", nrow(df.train))

df.test$part <- 
  rep("test", nrow(df.test))

# Merge the two data sets using rbind
df.all <- rbind(df.train, df.test)
view(df.all)

# Boxplot by group
boxplot(UUU ~as.factor(part), data = df.all)
kruskal.test(UUU ~as.factor(part), 
             data = df.all)
# Since the p-value is larger than 0.05
# the partition is validated with respect to 
# UUU variable

# Boxplot by group
a2 = boxplot(GUG ~as.factor(part), data = df.all)
kruskal.test(GUG ~as.factor(part), 
             data = df.all)
# Since the p-value is larger than 0.05
# the partition is validated with respect to 
# GUG variable

# Boxplot by group
a3 = boxplot(UAG ~as.factor(part), data = df.all)
kruskal.test(UAG ~as.factor(part), 
             data = df.all)
# Since the p-value is larger than 0.05
# the partition is validated with respect to 
# UAG variable

# Boxplot by group
a4 = boxplot(UGA ~as.factor(part), data = df.all)
kruskal.test(UGA ~as.factor(part), 
             data = df.all)
# Since the p-value is larger than 0.05
# the partition is validated with respect to 
# UAG variable

a5 = boxplot(AAA ~as.factor(part), data = df.all)
kruskal.test(AAA ~as.factor(part), 
             data = df.all)
# Since the p-value is larger than 0.05
# the partition is validated with respect to 
# AAA variable

a6 = boxplot(AUG ~as.factor(part), data = df.all)
kruskal.test(AUG ~as.factor(part), 
             data = df.all)
# Since the p-value is larger than 0.05
# the partition is validated with respect to 
# AAA variable
# Remove the variable part from df
df.train$part <- NULL
df.test$part <- NULL
View(df)

# Save train and test sets as separate csv files
write.csv(df.train, "codon_train.csv", row.names = TRUE)
write.csv(df.test, "codon_test.csv", row.names = TRUE)
