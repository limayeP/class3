# Trinary classification of DNA using frequency of codon usage
# Pallavi Limaye
# DATA 599 Special Project
# 11/8/2021
##############################################################################
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
              "skimr")

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

b <- read.csv("codon_usage_a_original_with_blanks.csv")
View(b)
names(b)
b$X <- NULL
View(b)
skimmed <- skim(b)
skimmed

#convert column 'UUC' from factor to numeric
b$UUC <- as.numeric(as.character(b$UUC))

# extracting positions of NA values
print ("Row and Col positions of NA values")
which(is.na(b), arr.ind=TRUE)

# Locate missing values
apply(is.na(b), 2, which)
# $UUU
# [1]  487 5064

# $UUC
# [1] 5064

# $UUA
# [1] 5064

# To impute values using KNN
# Develop the model that will be doing the imputing:
imputation_model <- preProcess(b,
                               method = c("knnImpute"))
# Use this model to predict the values for the missing data, 
# that is, to impute missing values.
b_imputed <- predict(imputation_model, b)
View(b_imputed)

# b_imputed contains the standardized values 
# to find the actual values de-standardize with caret
# na.rm for removing any missing values before calculating the mean
uuu_mean = mean(b$UUU,
               na.rm = TRUE)
uuu_sd = sd(b$UUU, na.rm = TRUE)

uuc_mean = mean(b$UUC,
                na.rm = TRUE)
uuc_sd = sd(b$UUC, na.rm = TRUE)

uua_mean = mean(b$UUA,
                na.rm = TRUE)
uua_sd = sd(b$UUA, na.rm = TRUE)

# What are the mean values
uuu_mean
# [1] 0.0248177
uuu_sd
# [1] 0.0176278

uuc_mean
# [1] 0.02344158

uuc_sd
# [1] 0.01159921

uua_mean
# [1] 0.02063674
uua_sd
# [1] 0.02070969

# In new columns To replace the NA values with imputed values
b$uuui <- ifelse(
  test = is.na(b$UUU) == TRUE, yes = b_imputed$UUU
  * uuu_sd + uuu_mean,
  no = b$UUU)

b$uuci <- ifelse(
  test = is.na(b$UUC) == TRUE, yes = b_imputed$UUC
  * uuc_sd + uuc_mean,
  no = b$UUC)

b$uuai <- ifelse(
  test = is.na(b$UUA) == TRUE, yes = b_imputed$UUA
  * uua_sd + uua_mean,
  no = b$UUA)

View(b)
# What are the imputed values
b$uuui[487]
# [1] 0.021128
b$uuui[5064]
# [1] 0.02737
b$uuci[5064]
# [1] 0.019106
b$uuai[5064]
#[1] 0.025466

# Replace NA in b with imputed values
b$UUU[487]
b$UUU[487] <- b$uuui[487]

b$UUU[5064]
b$UUU[5064] <- b$uuui[5064]

b$UUC[5064]
b$UUC[5064] <- b$uuci[5064]

b$UUA[5064]
b$UUA[5064] <- b$uuai[5064] 

View(b)
ab = subset(b, select = -c(uuui,uuci, uuai))
View(ab)
write.csv(ab,"codon_usage_b_imputed_without_reorganization.csv", row.names = TRUE)

##############################################################################
# For Trinary classification of DNAtypes, keep only 3 classes of DNA
##############################################################################
# Open the imputed file
ab <- read.csv("codon_usage_b_imputed_without_reorganization.csv")
summary(ab)
ab$X <- NULL
View(ab)

# Check for any NA values
sum(is.na(ab))
# View type of each column
sapply(ab, class)
# Convert DNAtype to a factor
ab$DNAtype <- as.factor(ab$DNAtype)
# View type of each column
sapply(ab, class)
summary(ab$DNAtype)
###############################################################################
# plot all DNA types
abc= ggplot(data=ab, aes(x=DNAtype, fill=DNAtype))+
  geom_bar(color = "Black")+
  scale_fill_brewer(palette="Paired")+
  ggtitle("Distribution of All DNAtypes")+
  theme(plot.title = element_text(size = 20, face = "bold"))

e = abc + theme(panel.border = element_rect(color = "black",
                                            fill = NA,
                                            size = 5))
cb = e+theme(axis.text=element_text(size=18),
             axis.title=element_text(size=18,face="bold"))

cb+theme(legend.key.height= unit(1, 'cm'),
         legend.key.width= unit(1, 'cm'))+
  theme(legend.title = element_text(size=30, face="bold")+
          theme(legend.text = element_text(size=30, face="bold")))
###########################################################################
# Keep only rows with DNAtype <= 2 (Nuclear, Mitochondrial and chloroplast DNA)
ef <- read.csv("codon_usage_b_imputed_without_reorganization.csv")
View(ef)
ef$X <- NULL
View(ef)
df <- subset(ef, DNAtype <= 2) 
# Check for DNAType more than 2 exist
any(df$DNAtype>=3)
# FALSE
##******************##
df$DNAtype <- as.factor(df$DNAtype)
# Plot distribution of the 3 DNAtypes
a = ggplot(data=df, aes(x=DNAtype, fill=DNAtype))+
  geom_bar(color = "Black")+
  scale_fill_brewer(palette="Paired")+
  ggtitle("Distribution of DNAtype: Nuclear, Mitochondrial and Chloroplast")+
  theme(plot.title = element_text(size = 20, face = "bold"))

b = a + theme(panel.border = element_rect(color = "black",
                                          fill = NA,
                                          size = 5))
c = b+theme(axis.text=element_text(size=18),
            axis.title=element_text(size=18,face="bold"))

c+theme(legend.key.height= unit(1, 'cm'),
        legend.key.width= unit(1, 'cm'))+
  theme(legend.title = element_text(size=30, face="bold")+
          theme(legend.text = element_text(size=30, face="bold")))
#################################################################################3
# Plot distribution of the 3 DNAtypes by Kingdom
d = ggplot(data=df, aes(x=DNAtype, fill=Kingdom))+
  geom_bar(color = "Black")+
  scale_fill_brewer(palette="Paired")+
  ggtitle("Distribution of DNAtype by kingdom")+
  theme(plot.title = element_text(size = 20, face = "bold"))

e = d + theme(panel.border = element_rect(color = "black",
                                          fill = NA,
                                          size = 5))
f = e+theme(axis.text=element_text(size=18),
            axis.title=element_text(size=18,face="bold"))

f+theme(legend.key.height= unit(1, 'cm'),
        legend.key.width= unit(1, 'cm'))+
  theme(legend.title = element_text(size=30, face="bold")+
          theme(legend.text = element_text(size=30, face="bold")))

# Normalized distribution of 3 DNAtypes by Kingdom
h = ggplot(data = df,aes(x=DNAtype, fill=Kingdom))+
  geom_bar(color = "black",
                 position = "fill")+
  scale_fill_brewer(palette="Paired")+
  ggtitle("Normalized Distribution of DNAtype by Kingdom")+
  theme(plot.title = element_text(size = 20, face = "bold"))

i = h + theme(panel.border = element_rect(color = "black",
                                          fill = NA,
                                          size = 5))
j = i+theme(axis.text=element_text(size=18),
            axis.title=element_text(size=18,face="bold"))

j+theme(legend.key.height= unit(1, 'cm'),
        legend.key.width= unit(1, 'cm'))+
  theme(legend.title = element_text(size=30, face="bold")+
          theme(legend.text = element_text(size=30, face="bold")))

#######################################################################
# Find the porportion of each DNAtype 
summary(df$DNAtype)

##############################################################################
names(df)
#  Remove “Kingdom”, “SpeciesName”, “Ncodons” and “SpeciesID
dt <- subset(df, select = -c(1,3,4,5))
names(dt)
View(dt)
dim(dt)
# [1] 12982    65

# Save the dataset with impute data, only DNAtypes:0, 1 and 2 and 64 codons
write.csv(dt,"codon_usage_c_imputed_nuc_mito_chloro.csv", row.names = TRUE)


