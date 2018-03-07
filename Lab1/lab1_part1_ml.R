####################################################################
# Machine Learning - MIRI Master
# Lluís A. Belanche
# Lab instructor: Marta Arias <marias@cs.upc.edu>

####################################################################
# LAB 1: Data pre-processing
# version of February 2018
# Cedric Bhihe, <cedric.bhihe@gmail.com>

# The goal is to predict whether someone should be granted a credit
# The target variable is called "Assessment"
rm(list=ls(all=TRUE))
   
####################################################################
# SECTION 1: READING THE FILE CREDSCO.TXT (loan data: credit scoring)
####################################################################

# Reading properly a data set is non-trivial because you need to know
# data format: decimal separator, column separator, is there a
# header?, how are strings quoted, how (if any) are missing values
# coded?, should character vectors be converted to factors?, should
# white spaces be stripped ?, ...

# Consult ?read.csv and play with useful control parameters:
# header=TRUE
# na.strings="?"
# dec = "."
# sep = ";"
# quote = "\"
# .. and some others

# After opening the file "credsco.csv" and inspecting it, we decide the following settings:
setwd("~/Documents/Study/UPC/MIRI-subjects/ML_machine-learning/Labs/Lab1/")
Credit <- read.csv("credsco.csv", header = TRUE, quote = "\"", dec = ".", check.names=TRUE)

# the dimensions of the data set are 

dim(Credit)
# ...which means 4,455 examples described by 14 variables

# Basic questions:
# Which is the target variable? where is it? how many different values? is it a classification problem or a regression problem?
###### Answer: the target variable is located in column 1 and is called 'Assessment'; it has two possible values
# What are the other variables?

names(Credit)

# inspect the first 4 examples

Credit[1:4,]

# inspect predictive variables 4, 5, 6 and 7 for the first example

Credit[1,5:8]

####################################################################
# SECTION 2: BASIC INSPECTION OF THE DATASET
####################################################################

# Perform a basic inspection of the dataset. Have a look at the minimum and maximum values for 
# each variable; find possible errors and abnormal values (outliers); find possible missing values 
# (these are marked by '?', '99999', 'NA' and sometimes symply by '0'); decide which variables should
# be continuous and which categorical; if there are mixed types, we have three options: recode continuous
# to categorical, recode categorical to continuous or leave them as they are. The thing is that some ML 
# methods feel more comfortable with one of these three options ...

# In essence, try to get a first idea of the best decision for each of these issues

summary(Credit)

# Assessment,Housing,MaritalStatus,Records,TypeOfJob are categorical and need to be treated properly
# In particular, Assessment is the target variable; we need to identify correct values

# Capital, ChargesOnCapital and Income present abnormally high maximums (99999999)
# There are also suspicious zeros, in both types of variables, which we identify as missing values

####################################################################
# SECTION 3: DEALING WITH MISSING VALUES
####################################################################

# Make a decision on a sensible treatment for the missing values and apply it; it is wise to 
# write down the possible consequences of this decision and the alternatives that could be
# considered in case the final results are not satisfactory

# the easiest way is of course to eliminate the involved rows or
# columns; this can be done only partially. For example, we could decide to
# eliminate the variables with the highest proportion of missing values.

# Deleting instances and/or variables containing missing values results
# in loss of relevant data and is also frustrating because of the effort
# in collecting the sacrificed information.

# CAREFUL! R does not know magically which entries are missing values

# this does not do anything, because there are no explicitly declared NA's:

Credit.complete <- na.omit(Credit)
dim(Credit.complete)

# the previous code does nothing! because there (still) are no explicitly declared NA's ...
# In the present case we decide to perform a step-by-step treatment, separate for the 
# categorical and continuous information

# We first decide to remove those rows with with missing values in the categorical variables 
# (there are a few)

table(Credit[,1]==0)
table(Credit[,3]==0)
table(Credit[,6]==0)
table(Credit[,8]==0)

Credit <- Credit[Credit[,1] != 0 & Credit[,3] != 0 & Credit[,6] != 0 & Credit[,8] != 0,]

dim(Credit)

# Process rows with missing values in the continuous variables (code 99999999)

attach(Credit)


# look at that:

hist(Income)
hist(Income[Income != 99999999])
hist(Income[Income != 99999999 & Income != 0], breaks=15)

# these are then clearly incorrect

table(Income == 99999999)

table(Income == 0)

table(Capital == 99999999)

table(ChargesOnCapital == 99999999)

# what do we do with this one? let's assume it is correct

table(YearsInJob == 0)


# Continuous variables have too many missing values, we can not eliminate them just like that: we must devise a treatment for these missing values

# first we mark them to 'NA' (so R knows!), including those from no 'Income'

Income[Income == 99999999 | Income == 0] <- NA
Capital[Capital == 99999999] <- NA
ChargesOnCapital[ChargesOnCapital == 99999999] <- NA

# see the difference

summary(Credit[,10])   # summary of raw Income column in Credit data table
summary(Income)      # summary of cleaned Income column

# imputation by 1st-NN (first nearest-neighbour): for every data example with a missing 'Income', we look for the most similar one (according to the remaining variables) and copy its 'Income' value
library(class) # load module "class" 

# Imputation of 'Income'
aux <- Credit[,-10]
dim(aux)
aux1 <- aux[!is.na(Income),]
dim(aux1)
aux2 <- aux[is.na(Income),]
dim(aux2)
# Neither of aux1, aux2 can contain NAs  
knn.inc <- knn(aux1,aux2,Income[!is.na(Income)])  # 1st-NN imputation
Income[is.na(Income)] <- as.numeric(as.character(knn.inc))  # translate factores in numerical values
 
# Imputation of 'Capital'
aux <- Credit[,-11]
aux1 <- aux[!is.na(Capital),]
aux2 <- aux[is.na(Capital),]
knn.cap <- knn(aux1,aux2,Capital[!is.na(Capital)])
Capital[is.na(Capital)] <- as.numeric(as.character(knn.cap))

# Imputation of 'ChargesOnCapital'
aux <- Credit[,-12]
aux1 <- aux[!is.na(ChargesOnCapital),]
aux2 <- aux[is.na(ChargesOnCapital),]
knn.cac <- knn(aux1,aux2, ChargesOnCapital[!is.na(ChargesOnCapital)])
ChargesOnCapital[is.na(ChargesOnCapital)] <- as.numeric(as.character(knn.cac))

ChargesOnCapital[Capital==0] <- 0

# assign back to the dataframe

Credit[,10] <- Income
Credit[,11] <- Capital
Credit[,12] <- ChargesOnCapital

# inspect again the result, especially the new statistics

dim(Credit)
summary(Credit)

####################################################################
# SECTION 4: TREATMENT OF MIXED DATA TYPES
####################################################################

# In this case we decide to keep the original type and leave the decision for later, 
# depending on the specific analysis
# First we explicitly declare some as categorical ('factor' in R)

Assessment    <- as.factor(Assessment)
Housing       <- as.factor(Housing)
MaritalStatus <- as.factor(MaritalStatus)
Records       <- as.factor(Records)
TypeOfJob     <- as.factor(TypeOfJob)

levels(Assessment)
levels(Housing)
levels(MaritalStatus)
levels(Records)
levels(TypeOfJob)

# not very nice, let's recode

levels(Assessment) <- c("positive","negative")
levels(Housing) <- c("rent","owner","private","ignore","parents","other")
levels(MaritalStatus) <- c("single","married","widower","split","divorced")
levels(Records) <- c("no","yes")
levels(TypeOfJob) <- c("indefinite","temporal","self-employed","other")

# WARNING! Some R programmers do not like 'attach', with some reason.
# We'll tackle that issue later.

is.factor(Assessment)
is.factor(Credit[,1]) # Original untransformed column "Assessment"

####################################################################
# SECTION 5: DERIVATION OF NEW VARIABLES: FEATURE EXTRACTION
####################################################################

# Decide whether it can be sensible to derive new variables; we decide to extract two new 
#   continuous and one new categorical variable (for the sake of illustration):

# Financing ratio (continuous)

FinancingRatio <- 100*AmountRequested/MarketPrice

hist(FinancingRatio, breaks=20)  # default is 'breaks=20'

# Saving capacity (continuous)

SavingCapacity <- (Income-Expenses-(ChargesOnCapital/100))/(AmountRequested/Deadline)

hist(SavingCapacity, breaks=16)

# Amount Requested greater than 1.2 times the median by people 1.5 times older than the mean (categorical):

Dubious <- rep("No", nrow(Credit))
Dubious[AmountRequested > 1.2*median(AmountRequested, na.rm = TRUE) & Age > 1.5*mean (Age, na.rm = TRUE)] <- "Yes"
Dubious <- as.factor(Dubious)

table(Dubious,Assessment)



####################################################################
# SECTION 6: WHAT WE HAVE DONE SO FAR
####################################################################

# create a new dataframe that gathers everything and inspect it again

Credit.new <- data.frame(Assessment,YearsInJob,Housing,Deadline,Age,MaritalStatus,
                         Records,TypeOfJob,Expenses,Income,Capital,ChargesOnCapital,
                         AmountRequested,MarketPrice,FinancingRatio,SavingCapacity,Dubious)
                   
summary(Credit.new)
dim(Credit.new)

detach(Credit)
attach(Credit.new)
is.factor(Credit.new[,1])

####################################################################
# SECTION 7: GAUSSIANITY AND TRANSFORMATIONS
####################################################################

# Goal: perform a graphical summary of some of the variables (both categorical and continuous)

# For continuous data:
# histograms and boxplots

hist(Income)
hist(Income,col=5)

hist(Capital)
hist(log10(Capital), breaks=20)

boxplot (Deadline)
title ("These are the deadlines")

boxplot (Age, col = "lightgray")
title ("and these are the ages")

boxplot(Credit.new[,9:16], outline=TRUE) 
boxplot(Credit.new[,9:16], outline=FALSE) # much better, but would be nicer one by one

# the previous plots suggest to take logs on some variables: Capital and ChargesOnCapital (we'll do it later)

# For categorical data:
# Frequency tables, Contingency tables, Bar charts, Pie charts

# should we treat Age as categorical? probably not

table(Age)                            

min(Age)                                                          
max(Age)                                                          
Age.cat <- cut(Age, breaks = seq(30, 90, 10)) # WARNING! we are generating NAs               
Age.cat

Age.cat <- cut(Age, breaks = seq(15, 75, 10))   
Age.tab <- table(Age.cat)                              
Age.tab
barplot(Age.tab)                                     # bar chart
pie(Age.tab)                                         # pie chart

# incidentally, this is how we could generate another new variable based on Age:

Age2.cat <- factor(as.integer(Age < 55))        
levels(Age2.cat) <- c("over55","under55")

TypeOfJob.Age <- table(TypeOfJob, Age.cat)            # contingency table
TypeOfJob.Age
margin.table(TypeOfJob.Age, 1)                        # row sums
margin.table(TypeOfJob.Age, 2)                        # column sums

prop.table(TypeOfJob.Age)                             # relative freqencies
round(prop.table(TypeOfJob.Age), digits=3)            # idem, rounded to 3 digits
round(prop.table(TypeOfJob.Age) * 100, digits=3)      # total percentages

round(prop.table(TypeOfJob.Age, 1), digits=3)         # table of relative frequencies (row-wise)
TypeOfJob.Age.rel <- round(prop.table(TypeOfJob.Age, 2), digits=3)      # table of relative frequencies (column-wise)

barplot(TypeOfJob.Age)                                # basic stacked bar chart

barplot(TypeOfJob.Age.rel, yaxt="n", xlab="Age", ylab="proportion", 
        col = c("white", "grey80", "grey40", "black"), 
        main = "TypeOfJob by Age", xlim=c(0,9))       # stacked bar chart

legend("bottomleft", legend=rownames(TypeOfJob.Age.rel), col="black", 
      fill = c("white", "grey80", "grey40", "black"), cex=0.65)

axis(2, at=seq(0, 1, 0.2))

barplot(TypeOfJob.Age.rel, beside = TRUE)             # grouped bar chart

# we can perform graphical comparisons between some pairs of variables (both categorical 
# and continuous), using the plot(),  pairs() and identify() procedures

par(mfrow=c(1,2))
plot (AmountRequested, Capital, main = "Amount req. vs. Market price", 
      cex = .5, col = "dark red")

plot (log10(AmountRequested), log10(Capital+1), 
      main = "Amount req. vs. Market price (better)", 
      cex = .5, col = "dark red")

plot(AmountRequested, Capital+1, 
     main = "Amount req. vs. Market price (even better)", 
     cex = .5, col = "dark red", log="xy")

# adding a center (dashed) and a regression line (blue)

abline(v = mean(log10(AmountRequested)), lty = 2, col = "grey40")
abline(h = mean(log10(Capital+1)), lty = 2, col = "grey40")
abline(lm(log10(Capital+1) ~ log10(AmountRequested)), col = "blue")

# (note that log10(x+1)=0 for x=0, so our transformation keeps the zeros)
# On the other hand, these same zeros spoil the regression ...


par(mfrow=c(1,1))   # go back to plotting 1 plot on one row

plot (TypeOfJob, AmountRequested)

plot (Age,Assessment, col="red") # WARNING!
plot (Assessment, Age, col="red") # better

plot (Assessment, TypeOfJob,col="blue", xlab="Assessment",ylab="TypeOfJob")

pairs(~ AmountRequested + Capital + Age)   # Each term in the formula ~ x + y + z
                                           # will give a separate variable in the
                                           # pairs plot

# Plotting a variable against the normal (Gaussian) pdf in red
hist.with.normal <- function (x, main, xlabel=deparse(substitute(x)), ...)
{
  h <- hist(x,plot=F, ...)
  s <- sd(x)
  m <- mean(x)
  ylim <- range(0,h$density,dnorm(0,sd=s))
  hist(x,freq=FALSE,ylim=ylim,xlab=xlabel, main=main, ...)
  curve(dnorm(x,m,s),add=T,col="red")
}

hist.with.normal(Expenses, "Expenses")

par(mfrow=c(2,4))   # plot 4 charts per row on 2 rows
for (i in 0:7) 
  { plot(density(Credit.new[,i+9]), xlab="", main=names(Credit.new)[i+9]) }

# do any of the continuous variables "look" Gaussian? 
# features to look for in comparing to a Gaussian: outliers, asymmetries, long tails

# A useful tool for "Gaussianization" is the Box-Cox power transformation

library(MASS)

par(mfrow=c(1,3))

hist(Capital, main="Look at that ...")

bx <- boxcox(I(Capital+1) ~ . - Assessment, data = Credit.new,
             lambda = seq(-0.25, 0.25, length = 10))

lambda <- bx$x[which.max(bx$y)]

Capital.BC <- (Capital^lambda - 1)/lambda

hist(Capital.BC, main="Look at that now!")

par (mfrow=c(1,1))


####################################################################
# SECTION 8: DEALING WITH MIXED TYPES OF VARIABLES (REVISITED)
####################################################################

# RECODING CATEGORICAL VARIABLES (DUMMY BINARY CODING)

# Most R functions for modelling accept factors mixed with numeric data; the factors are silently recoded using a dummy binary code, so there is no need to handle this explicitly

# RECODING THE CONTINUOUS VARIABLES (TO ASSESS NON-LINEARITIES)

YearsInJob.CAT         <- cut(YearsInJob, breaks=c(-1,1,3,8,14,99))
Deadline.CAT           <- cut(Deadline, breaks=c(0,12,24,36,48,99))
Age.CAT                <- cut(Age, breaks=c(0,25,30,40,50,99))
Expenses.CAT           <- cut(Expenses, breaks=c(0,40,50,60,80,9999))
Income.CAT             <- cut(Income, breaks=c(0,80,110,140,190,9999))
Capital.CAT            <- cut(Capital, breaks=c(-1,0,3000,5000,8000,999999))
ChargesOnCapital.CAT   <- cut(ChargesOnCapital, breaks=c(-1,0,500,1500,2500,999999))
AmountRequested.CAT    <- cut(AmountRequested, breaks=c(0,600,900,1100,1400,99999))
MarketPrice.CAT        <- cut(MarketPrice, breaks=c(0,1000,1300,1500,1800,99999))
FinancingRatio.CAT     <- cut(FinancingRatio, breaks=c(0,50,70,80,90,100))
SavingCapacity.CAT     <- cut(SavingCapacity, breaks=c(-99,0,2,4,6,99))

levels(YearsInJob.CAT)         <- paste("YearsInJob",levels(YearsInJob.CAT))
levels(Deadline.CAT)           <- paste("Deadline",levels(Deadline.CAT))
levels(Age.CAT)                <- paste("Age",levels(Age.CAT))
levels(Expenses.CAT)           <- paste("Expenses",levels(Expenses.CAT))
levels(Income.CAT)             <- paste("Income",levels(Income.CAT))
levels(Capital.CAT)            <- paste("Capital",levels(Capital.CAT))
levels(ChargesOnCapital.CAT)   <- paste("ChargesOnCapital",levels(ChargesOnCapital.CAT))
levels(AmountRequested.CAT)    <- paste("AmountRequested",levels(AmountRequested.CAT))
levels(MarketPrice.CAT)        <- paste("MarketPrice",levels(MarketPrice.CAT))
levels(FinancingRatio.CAT)     <- paste("FinancingRatio",levels(FinancingRatio.CAT))
levels(SavingCapacity.CAT)     <- paste("SavingCapacity",levels(SavingCapacity.CAT))

Credit.new.CAT <- data.frame(Assessment,YearsInJob.CAT,Housing,Deadline.CAT,Age.CAT,MaritalStatus,Records,TypeOfJob,Expenses.CAT,Income.CAT,Capital.CAT,ChargesOnCapital.CAT,AmountRequested.CAT,MarketPrice.CAT,FinancingRatio.CAT,SavingCapacity.CAT)


# notice the difference:

Credit.new[1,]
dim(Credit.new)

Credit.new.CAT[1,]
dim(Credit.new.CAT)


##############################################
# SECTION 9: STATISTICAL ANALYSIS
##############################################

summary(Credit.new) # yet again

# basic summary statistics by groups

library(psych) # describeBy

range.cont <- c(2,4,5,9:16)

describeBy (Credit.new[,range.cont], Credit.new$Assessment)

# Plot of all pairs of some continuous variables according to the class (Assessment variable)

pairs(Credit.new[,c(2,5,9,16)], main = "Credit Scoring DataBase", col = (1:length(levels(Credit.new$Assessment)))[unclass(Credit.new$Assessment)])

# (this rather shows how difficult this problem is)
                        
#### Feature selection analysis

# Feature selection for continuous variables: Fisher's F

names(Credit.new)[range.cont]

pvalcon <- NULL

varc <- list(YearsInJob, Deadline, Age, Expenses, Income, Capital, ChargesOnCapital, AmountRequested, MarketPrice, FinancingRatio, SavingCapacity)

for (i in 1:length(varc)) 
  pvalcon[i] <- (oneway.test (varc[[i]]~Credit.new$Assessment))$p.value

pvalcon <- matrix(pvalcon)
row.names(pvalcon) <- colnames(Credit.new)[range.cont]

# Ordered list of continuous variables according to their association to Assessment (greater to lesser)

sort(pvalcon[,1])

# Graphical representation of Assessment

ncon <- nrow(pvalcon)

par (mfrow=c(3,4))

for (i in 1:ncon) 
{
  barplot (tapply(varc[[i]], Credit.new$Assessment, mean),main=paste("Means by",row.names(pvalcon)[i]), las=2, cex.names=0.75)
  abline (h=mean(varc[[i]]))
  legend (0,mean(varc[[i]]),"Global mean",bty="n") 
}


mosthighlycorrelated <- function(mydataframe,numtoreport)
{
  # find the correlations
  cormatrix <- cor(mydataframe)
  # set the correlations on the diagonal or lower triangle to zero,
  # so they will not be reported as the highest ones:
  diag(cormatrix) <- 0
  cormatrix[lower.tri(cormatrix)] <- 0
  # flatten the matrix into a dataframe for easy sorting
  fm <- as.data.frame(as.table(cormatrix))
  # assign human-friendly names
  names(fm) <- c("First.Variable", "Second.Variable","Correlation")
  # sort and print the top n correlations
  head(fm[order(abs(fm$Correlation),decreasing=TRUE),],n=numtoreport)
}

mosthighlycorrelated(Credit.new[,range.cont], 10)


####################################################################
# SECTION 9: MISCELLANEOUS ISSUES
####################################################################
			    
# Shuffling the data (to avoid possible biases)

set.seed (104)
Credit.new <- Credit.new[sample.int(nrow(Credit.new)),]
 

####################################################################
# SECTION 10: SAVE THE PREPROCESSED DATA INTO A FILE FOR FUTURE USE
####################################################################

# WARNING! This is a binary file, unless 'ascii' is TRUE

save(Credit.new, file = "Credsco-processed.Rdata")

# remove (almost) everything in the working environment
rm(list = ls())

Credit.new

load("Credsco-processed.Rdata")
objects()
summary(Credit.new)