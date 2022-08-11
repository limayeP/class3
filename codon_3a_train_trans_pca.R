packages <- c("ggplot2",
              "rattle",
              "Hmisc",
              "psych",
              "caret",
              "e1071",
              "corrplot",
              "cowplot",
              "ggpubr",
              "DataExplorer",
              "MASS",
              "car")

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

bt <- read.csv("codon_train_transformed.csv")
View(bt)
bt$X <- NULL
View(bt)
# View type of each column
sapply(bt, class)
# Convert DNAtype to a factor
# d$DNAtype <- as.factor(d$DNAtype)
dim(bt)
names(bt)
summary(bt)
View(bt)
# Remove the target DNAtype
tpca <- subset(bt, select = -c(1))
names(tpca)

#########################################################################3
# PCA with varimax rotation all 64 predictors (nfactors = 64)
pca64<- principal(tpca, nfactors = 64,
                   rotate = "varimax",
                   scores = TRUE)
pca64
# Eigenvalues
pca64$values

# Loadings matrix, variance explained
pca64$loadings

print(pca64$loadings,
      cutoff = 0.5)

# Scree plot
plot(pca64$values,
     type = "b",
     main = "screeplot with pca64");
abline(h = 1, lty = 2)

round(pca64$communality, 2)

##################################################################################################################
#Extract information from the pca pbject predictors as rows and principal components as columns
loadtable_64a <-data.frame(pca64$loadings[1:62, 1:62])
write.csv(loadtable_64a,"codon_train_trans_pca_64_loadings.csv", row.names = TRUE)
#pca loadings as rows and principal components as columns
output_64 <-data.frame(pca64$Vaccounted)
write.csv(output_64,"codon_train_trans_pca_64_Vaccounted_only.csv", row.names = TRUE)
#Extract information from the pca pbject predictors as rows and principal components as columns
# COMBINE COMPONENT SOCRES WITH THE ORIGINAL DATA SET
c_s_64 <- cbind(tpca, pca64$scores)
write.csv(c_s_64,"codon_train_trans_64scores_dataset.csv", row.names = TRUE)
# correlations to understand the relationship between components and DNAtype
dna <- bt["DNAtype"]
# For building models the scores will be used
c_s_dna_64 <- cbind(dna,pca64$scores)
write.csv(c_s_dna_64 ,"codon_train_trans_64pcascores_for_model_building.csv", row.names = TRUE)
rcorr(as.matrix(c_s_dna_64))


##########################################################
########################################################################################################################test
#######################################################################################################3
# validate with test set
########################################################################################################3
#################################################################################################################3

testbt <- read.csv("codon_test_transformed.csv")
View(testbt)
testbt$X <- NULL
View(testbt)
# Convert DNAtype to a factor
# d$DNAtype <- as.factor(d$DNAtype)
dim(testbt)
names(testbt)
summary(testbt)
##############################################################################
# Remove the target variables "DNAtype"
testdpca <- subset(testbt, select = -c(1))
names(testdpca)
head(testdpca)
# PCA with varimax rotation 59 predictors (nfactors = 64)
testpca64<- principal(testdpca, nfactors = 64,
                      rotate = "varimax",
                      scores = TRUE)
testpca64

print(testpca64$loadings,
      cutoff = 0.5)
DNAtype <- testbt$DNAtype
# For building models the scores will be used
c_s_dna_64 <- cbind(DNAtype,testpca64$scores)
write.csv(c_s_dna_64 ,"codon_test_trans_64pcascores.csv", row.names = TRUE)
rcorr(as.matrix(c_s_dna_64))


#####################################################################################

# wdbc.pr <- prcomp(tpca, center = TRUE)
# summary(wdbc.pr)
# 
# screeplot(wdbc.pr, type = "l", npcs = 15, main = "Screeplot of the first 10 PCs")
# abline(h = 1, col="red", lty=5)
# legend("topright", legend=c("Eigenvalue = 1"),
#        col=c("red"), lty=5, cex=0.6)
# 
# cumpro <- cumsum(wdbc.pr$sdev^2 / sum(wdbc.pr$sdev^2))
# plot(cumpro[0:15], xlab = "PC #", ylab = "Amount of explained variance", main = "Cumulative variance plot")
# abline(v = 6, col="blue", lty=5)
# abline(h = 0.88759, col="blue", lty=5)
# legend("topleft", legend=c("Cut-off @ PC6"),
#        col=c("blue"), lty=5, cex=0.6)


#################################################################################




#