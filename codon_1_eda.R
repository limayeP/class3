packages <- c("ggplot2",
              "rattle",
              "Hmisc",
              "psych",
              "caret",
              "e1071",
              "corrplot",
              "cowplot",
              "ggpubr",
              "skimr",
              "DataExplorer",
              "purrr",
              "dplyr"
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

ad <- read.csv("codon_usage_c_imputed_nuc_mito_chloro.csv", header = TRUE)
names(ad)
# Drop "X" column
# Get only predictors
adf<- subset(ad, select = -c(1))
names(adf)
# Observe the data type
sapply(adf, class)
# [1] 12982    65
adf$DNAtype <- as.factor(adf$DNAtype)
# ANy NAs
sum(is.na(adf))
summary(adf)
View(adf)
skimmed <- skim(adf)
skimmed
adt <- adf
#############################################################################
# Find the porportion of each DNAtype 
summary(adt$DNAtype)
t1 <- table(adt$DNAtype)
#   0    1    2 
# 9267 2899  816 

length(adt$DNAtype)
# [1] 12982

# choose the value with which every row value is divided
someVar = 12982
# apply for each row, where row_val equals the values in your rows
round(apply(t1, 1, function(row_val) row_val / someVar), 3)
# 0       1       2
# 0.714 0.223  0.063
t2 = round(apply(t1, 1, function(row_val) row_val / someVar)*100, 2)
DNAtype = c('Nuclear','Mitochondrial', 'Chloroplast')
t3 = data.frame(rbind(DNAtype, t2))
colnames(t3) = c('DNAtype','Percent distribution')
library(data.table)
t4 = transpose(t3)
colnames(t4) = c("DNAtype", "dist")

t4$dist = as.numeric(as.character(t4$dist))
t4$DNAtype = as.factor(as.character(t4$DNAtype))

# plot of categorical target
ggplot(t4,aes(x=DNAtype, y=dist, fill=DNAtype))+
  geom_bar(stat='identity', color="Black")+
  coord_flip()+
  xlab("DNA class")+
  ylab("Percent Distribution")+
  scale_fill_grey(start = 0.25, end = 0.75) +
  theme(legend.position="none")+
  theme(text = element_text(size = 14, color="black"))+
  theme(axis.text = element_text(color="black"))+
  theme(panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 2))

# plot of categorical target
ggplot(adt,aes(DNAtype))+
  geom_bar(aes(fill=DNAtype))+
  xlab("DNAtype")+
  ylab("Frequency")+
  ggtitle("Bar Chart of DNAtype")
##########################################################################
# With DataExplorer understand the data
create_report(adt, y = "DNAtype")
# View basic description
introduce(adt)
plot_intro(adt)
plot_missing(adt)
## View histogram of continuous variables
plot_histogram(adt)
#############################################################################
# All the codon frequencies were right skewed.

# skewness values for codon frequencies determined with caret package
# Get the predictors only
ads = subset(adt, select = -c(1))
View(ads)
names(ads)
skewValues <- apply(ads, 2, skewness)
skewValues

## View estimated density distribution of all continuous variables
plot_density(adt)

## View quantile-quantile plot of all continuous variables
plot_qq(adt)

# ## View quantile-quantile plot of all continuous variables by feature `DNAtype`
plot_qq(adt, by = "DNAtype")

### View overall correlation heatmap of just the predictors
plot_correlation(ads)

#View bivariate continuous distribution based on `DNAtype`
plot_boxplot(adt, by = "DNAtype")
#############################################################################
# To filter for near-zero variance predictors
nearZeroVar(ads)

describeBy(adt$UAG, adt$DNAtype)
kruskal.test(UAG~ DNAtype,
             data = adt)
# UAG not to be removed.
#######################################################################
# Plot correlations between predictors in the transformed training set
correlations  <- cor(ads)
corrplot(correlations, order= "hclust")
rdc <- cor.mtest(ads, conf.level = 0.95)
rdc$p
write.csv(rdc$p,"p_values_correlations.csv", row.names = TRUE)

corrplot(correlations, p.mat = rdc$p, order = "hclust", sig.level = 0.05, insig = "blank", addrect=1)
corrplot(correlations, p.mat = rdc$p, order = "hclust", sig.level = 0.05, insig = "pch", addrect=1)

# To eliminate the minimum # of predictors to achieve all pairwise
# correlations more than 0.75.
highCorr <- findCorrelation(correlations, cutoff = 0.75)
length(highCorr)
highCorr
print(names(ads[ , highCorr]))
# Find the correlation between these variables
pairs_correlation <- ads[ , highCorr]
pairs(pairs_correlation)

describeBy(adt$GGC, adt$DNAtype)
kruskal.test(GGC~ DNAtype,
             data = adt)

describeBy(adt$CGC, adt$DNAtype)
kruskal.test(CGC~ DNAtype,
             data = adt)

describeBy(adt$AUU, adt$DNAtype)
kruskal.test(AUU~ DNAtype,
             data = adt)

describeBy(adt$CCG, adt$DNAtype)
kruskal.test(CCG~ DNAtype,
             data = adt)

describeBy(adt$AAU, adt$DNAtype)
kruskal.test(AAU~ DNAtype,
             data = adt)

describeBy(adt$UUA, adt$DNAtype)
kruskal.test(UUA~ DNAtype,
             data = adt)

describeBy(adt$GCG, adt$DNAtype)
kruskal.test(GCG~ DNAtype,
             data = adt)

describeBy(adt$CUA, adt$DNAtype)
kruskal.test(CUA~ DNAtype,
             data = adt)
# Scatterplots of correlated variables
pairs_correlation <- ads[ , highCorr]
pairs.panels(pairs_correlation,
             method = "pearson",
             hist.col = "#00AFBB",
             density = FALSE,
             ellipses = FALSE,
             stars = TRUE)
#######################################################################
pca1 <- principal(ads, nfactors = 8,
                  rotate = "varimax",
                  scores = TRUE)
pca1
project.pca <- prcomp(t(ads), scale = F)
summary(project.pca)   

#Determine the proportion of variance of each component
#Proportion of variance equals (PC stdev^2) / (sum all PCs stdev^2)
project.pca.proportionvariances <- ((project.pca$sdev^2) / (sum(project.pca$sdev^2)))*100
df_pca_data=data.frame(PC1 = project.pca$x[,1], PC2 = project.pca$x[,2], sample = colnames(ads))
# df_pca_data=data.frame(PC1 = project.pca$x[,1], PC2 = project.pca$x[,2], sample = 1:19)

png('pcaplot_1.png')
ggplot(df_pca_data, aes(PC1,PC2, color = sample))+
  geom_point(size=8)+
  labs(x=paste0("PC1, ", round(project.pca.proportionvariances[1], 2), "%"),y=paste0("PC2, ", round(project.pca.proportionvariances[2], 2), "%"))
dev.off()
#Scree plot
barplot(project.pca.proportionvariances, cex.names=1, xlab=paste("Principal component (PC), 1-", length(project.pca$sdev)), ylab="Proportion of variation (%)", main="Scree plot", ylim=c(0,100))

######################################################################
ggplot(adt, aes(x=CCG,y=GCC, color = DNAtype, shape = DNAtype)) +geom_point()
ggplot(adt, aes(x=CGC,y=AUU, color = DNAtype, shape = DNAtype)) +geom_point()

ggplot(adt, aes(x=AUU,y=CCG, color = DNAtype, shape = DNAtype)) +geom_point()

ggplot(adt, aes(x=CCG,y=AAU, color = DNAtype, shape = DNAtype)) +geom_point()

ggplot(adt, aes(x=AAU,y=UUA, color = DNAtype, shape = DNAtype)) +geom_point()

ggplot(adt, aes(x=UUA,y=GCG, color = DNAtype, shape = DNAtype)) +geom_point()

ggplot(adt, aes(x=GCG,y=CUA, color = DNAtype, shape = DNAtype)) +geom_point()
ggplot(adt, aes(x=GGC,y=CGC, color = DNAtype, shape = DNAtype)) +geom_point()
ggplot(adt, aes(x=GGC,y=CUA, color = DNAtype, shape = DNAtype)) +geom_point()
ggplot(adt, aes(x=UUA,y=AUU, color = DNAtype, shape = DNAtype)) +geom_point()
ggplot(adt, aes(x=GCG,y=CCG, color = DNAtype, shape = DNAtype)) +geom_point()
ggplot(adt, aes(x=GGC,y=CCG, color = DNAtype, shape = DNAtype)) +geom_point()
ggplot(adt, aes(x=GGC,y=GCG, color = DNAtype, shape = DNAtype)) +geom_point()

# Find linear combinations (linear dependencies)
comboInfo <- findLinearCombos(ads)
comboInfo
d <- adt
# scatterplots
plot(d$UUU,d$AAA)
plot(d$UUU,d$UUA)
plot(d$UUU,d$UUG)
# Scatterplot with Overlay of Churn
par(mfrow=c(1,3))
ggplot(d, aes(x=UUU,y=AAA, color = DNAtype, shape = DNAtype)) +geom_point()+
  geom_smooth(se=FALSE)
ggplot(d, aes(x=UUU,y=UUA, color = DNAtype, shape = DNAtype)) +geom_point()
ggplot(d, aes(x=UUU,y=UUG, color = DNAtype, shape = DNAtype)) +geom_point()
###################################################################
# Having visualised the relationships between X and Y, 
# We can only say which variables are likely to be important to predict Y.
# It may not be wise to conclude which variables are NOT important.
# Because sometimes, variables with uninteresting pattern can help 
# explain certain aspects of Y that the visually important variables may not.
# So to be safe, let’s not arrive at conclusions 
# about excluding variables prematurely.
######################################################################
# 1) Glycine
gly<- adt %>% 
  group_by(DNAtype) %>% 
  summarise(
    CAA_mean = mean(CAA),
    CAG_mean = mean(CAG),
    CAA_sd = sd(CAA),
    CAG_sd = sd(CAG),
    CAA_se = CAA_sd / length(CAA),
    CAG_se = CAA_sd / length(CAG),
    )
gly <- data.frame(gly) 

affp<- ggplot(gly, aes(x = DNAtype, y = CAA_mean, fill=DNAtype)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=CAA_mean - CAA_se, ymax= CAA_mean + CAA_se),
                width=.2,
                position=position_dodge(.9))

affq = affp + xlab("CAA")+
  ylab("Mean")+
  ggtitle("Mean CAA (Glycine) by DNAtype ")+
  theme(plot.title = element_text(size = 20, face = "bold"))
affq
###################################################################
# Histogram of predictors with target overlay
# Shortened plot functions
# Names on the plot have to be changed manually
overlap_plot <- function(data, n) {
  g = ggplot(data,aes(n))+
    geom_histogram(aes(fill=DNAtype),
                   color = "black", bins=30)+
    xlab("UGA")+
    ylab("Count")+
    ggtitle("UGA with DNAtype")+
    theme(plot.title = element_text(size = 16, face = "bold"))
  
  y = g + theme(panel.border = element_rect(color = "black",
                                            fill = NA,
                                            size = 2))
  f = y+theme(axis.text=element_text(size=14),
              axis.title=element_text(size=15,face="bold"))
  
  f+theme(legend.key.height= unit(1, 'cm'),
          legend.key.width= unit(1, 'cm'))+
    theme(legend.title = element_text(size=30, face="bold")+
            theme(legend.text = element_text(size=30, face="bold")))
}

norm_overlap_plot <- function(data, n){
  h = ggplot(data,aes(n))+
    geom_histogram(aes(fill = DNAtype), color = "black", 
                   bins= 30,
                   position = "fill")+
    xlab("UGA")+
    ylab("Count")+
    ggtitle("Normalized UGA with DNAtype")+
    theme(plot.title = element_text(size = 16, face = "bold"))
  
  
  j = h+theme(panel.border = element_rect(color = "black",
                                          fill = NA,
                                          size = 2))#change legend title font size
  
  i = j+theme(axis.text=element_text(size=14),
              axis.title=element_text(size=15,face="bold"))
  i +theme(legend.key.height= unit(1, 'cm'),
           legend.key.width= unit(1, 'cm'))+
    theme(legend.title = element_text(size=30, face="bold")+
            theme(legend.text = element_text(size=30, face="bold")))
}

boxplot_fig <- function(data, n){
  ggplot(data, aes(x=DNAtype, y = UGA, color=DNAtype), color = "black") +
    geom_boxplot()+
    ggtitle("boxplot of UGA with DNAtype")+
    theme(panel.border = element_rect(color = "black",
                                      fill = NA,
                                      size = 2))+#change legend title font size
    theme(plot.title = element_text(size = 16, face = "bold"))+
    theme(axis.text=element_text(size=14),
          axis.title=element_text(size=15,face="bold"))+
    theme(legend.key.height= unit(1, 'cm'),
          legend.key.width= unit(1, 'cm'))+
    theme(legend.title = element_text(size=30, face="bold")+
            theme(legend.text = element_text(size=30, face="bold")))
}

########################################################################
bo1 = overlap_plot(adt, adt$GAG)
bo2 = norm_overlap_plot(adt, adt$GAG)
bo3 = boxplot_fig(adt, adt$GAG)
plot_grid(bo1,bo2,bo3,nrow = 1, rel_widths = c(.5,.5,.5), rel_heights = c(.5, .5,.5))
# Check significance of difference in boxplots
box_sig <- function(data, n){
describeBy(adt$AAA, adt$DNAtype)
kruskal.test(AAA~ DNAtype,
             data = adt)
}

# To check for Glutamine (Q) CAA (most preferred), CAG (least preferred)
describeBy(adt$CAA, adt$DNAtype)
kruskal.test(CAA~ DNAtype,
             data = adt)

describeBy(adt$CAG, adt$DNAtype)
kruskal.test(CAG~ DNAtype,
             data = adt)
##############################################################
# To check for Leucine (L) CUC and CUG (least preferred)
describeBy(adt$CUC, adt$DNAtype)
kruskal.test(CUC~ DNAtype,
             data = adt)

describeBy(adt$CUG, adt$DNAtype)
kruskal.test(CUG~ DNAtype,
             data = adt)

describeBy(adt$CUA, adt$DNAtype)
kruskal.test(CUA~ DNAtype,
             data = adt)

describeBy(adt$UUG, adt$DNAtype)
kruskal.test(UUG~ DNAtype,
             data = adt)

describeBy(adt$CUU, adt$DNAtype)
kruskal.test(CUU~ DNAtype,
             data = adt)

describeBy(adt$CUA, adt$DNAtype)
kruskal.test(CUA~ DNAtype,
             data = adt)

describeBy(adt$UUA, adt$DNAtype)
kruskal.test(UUA~ DNAtype,
             data = adt)
#############################################################

# To check for Asparagine(N)AAC, AAU
describeBy(adt$AAC, adt$DNAtype)
kruskal.test(AAC~ DNAtype,
             data = adt)

describeBy(adt$AAU, adt$DNAtype)
kruskal.test(AAU~ DNAtype,
             data = adt)
##############################################################
# To check for Glycine (G) GGA, GGU, GGG, GCC
describeBy(adt$GGA, adt$DNAtype)
kruskal.test(GGA~ DNAtype,
             data = adt)

describeBy(adt$GGU, adt$DNAtype)
kruskal.test(GGU~ DNAtype,
             data = adt)

describeBy(adt$GGG, adt$DNAtype)
kruskal.test(GGG~ DNAtype,
             data = adt)

describeBy(adt$GGC, adt$DNAtype)
kruskal.test(GGC~ DNAtype,
             data = adt)
##################################################################
# To check foR threonine (T) ACU, ACC, ACA, ACG
describeBy(adt$ACU, adt$DNAtype)
kruskal.test(ACU~ DNAtype,
             data = adt)

describeBy(adt$ACC, adt$DNAtype)
kruskal.test(ACC~ DNAtype,
             data = adt)

describeBy(adt$ACA, adt$DNAtype)
kruskal.test(ACA~ DNAtype,
             data = adt)

describeBy(adt$ACG, adt$DNAtype)
kruskal.test(ACG~ DNAtype,
             data = adt)
#################################################################
# # To check for tyrosine (Y) UAU, UAC
describeBy(adt$UAU, adt$DNAtype)
kruskal.test(UAU~ DNAtype,
             data = adt)

describeBy(adt$UAC, adt$DNAtype)
kruskal.test(UAC~ DNAtype,
             data = adt)

#################################################################
# # To check foR serine () UAU, UAC
describeBy(adt$AGU, adt$DNAtype)
kruskal.test(AGU~ DNAtype,
             data = adt)

describeBy(adt$AGC, adt$DNAtype)
kruskal.test(AGC~ DNAtype,
             data = adt)
###############################################################
# # To check for STOP codons UAG, UGA, UAA
describeBy(adt$UAG, adt$DNAtype)
kruskal.test(UAG~ DNAtype,
             data = adt)

describeBy(adt$UGA, adt$DNAtype)
kruskal.test(UGA~ DNAtype,
             data = adt)

describeBy(adt$UAA, adt$DNAtype)
kruskal.test(UAA~ DNAtype,
             data = adt)
##################################################################
# The p-value is quite close to significance!
# We can pretty much reject the hypothesis that the two
# groups have the same distribution.
# In other words, we’ve got some evidence here that
# AAA, AAU,  is a significant predictive contributer