####################################################################
# Machine Learning - MIRI Master
# Lluís A. Belanche

# LAB 1: Data pre-processing (practise)
# version of February 2018
####################################################################

# This exercise involves the use of the 'Auto' data set, which can be found in the file 'Auto.data'. 
# The file contains a number of variables for cars of different times.

graphics.off()      # reset/close all graphical devices 


##### Reading the file 'Auto.data'

Auto <- read.table("Auto.data", header=TRUE, na.strings="?")

# put proper country of origin
Auto[,"origin"] <- factor(c("USA","EU","Japan")[Auto[,"origin"]])

# convert "miles per gallon" to "liters per km"
Auto[,"mpg"] <- 235.4/Auto[,"mpg"]
colnames(Auto)[1] <- "L100km"

# car name is not useful for modelling, but it may be handy to keep it as the row name

# WARNING! surprisingly, car names are not unique, so we first prefix them by their row number

Auto$name <- paste (1:nrow(Auto), Auto$name)
rownames(Auto) <- Auto$name
Auto <- subset (Auto, select=-name)

# so this is your departing data set
summary(Auto)
attach(Auto)

# maybe you remember a plot from Lecture 1 ... (it was mpg against HP, though)

Auto.lm <- lm(L100km ~ horsepower, Auto)

plot(Auto[,"horsepower"],Auto[,"L100km"],
     pch=20,
     xlab="horsepower",ylab="fuel consumption (l/100km)",
     main="Linear regression")

# here 'x' is "horsepower"

a <- Auto.lm$coefficients["(Intercept)"]
b <- Auto.lm$coefficients["horsepower"]
abline(a=a,b=b,col="blue")
text(50,25,sprintf("y(x)=%.3fx+%.2f",b,a),col="red",pos=4)

# Now if you want to create a pdf with the previous plot, just do

# pdf("horsepower.pdf")

# the previous code would be embedded here
  
# dev.off()
  
# In order to crate quick LaTeX code, try this:
  
install.packages("xtable")
library(xtable)

xtable(Auto[1:4,])

# Was that nice? 
# this is a list of R objects that can be embedded into a LaTeX table code:

methods(xtable)

# Wait ... with all these graphics and facilities, we may have gotten a bit obfuscated ...

# Probably we went too fast ... 

table(Auto$cylinders)

# that's strange, some cars have an odd number of cylinders (are these errors?)

subset(Auto,cylinders==3)

# These Mazdas wear a Wankel engine, so this is correct

subset(Auto,cylinders==5)

# Yes, these Audis displayed five-cylinder engines, so the data is correct

# but, from summary(Auto) above we see that horsepower has 5 NA's that we'll need to take care of ...


###################################################################################
# Assignment for the lab session
###################################################################################

# 1. print the dimensions of the data set 
# 2. identify possible target variables according to classification or regression problems
# 3. inspect the first 4 examples and the predictive variables 6 and 7 for the tenth example
# 4. perform a basic inspection of the dataset. Have a look at the minimum and maximum values for each variable; find possible errors and abnormal values (outliers); find possible missing values; decide which variables are continuous and which are categorical
# 5. make a decision on a sensible treatment for the missing values and apply it; 
#    WARNING: 'origin' is categorical and cannot be used for knn imputation, unless you make it binary temporarily
# 6. derive one new continuous variable: weight/horsepower; derive one new categorical variable: sports_car, satisfying horsepower > 1.2*mean(horsepower) AND acceleration < median(acceleration); do you think this new variable is helpful in predicting 'origin' ?
# 7. create a new dataframe that gathers everything and inspect it again
# 8. perform a graphical summary of some of the variables (both categorical and continuous)
# 9. perform a graphical comparison between some pairs of variables (both categorical and continuous)
# 10. do any of the continuous variables "look" Gaussian? can you transform some variable so that it looks more so?
# 11. choose one of the continuous variables and recode it as categorical; choose one of the categorical variables and recode it as continuous
# 12. create a new dataframe that gathers everything and inspect it again; consider 'origin' as the target variable; perform a basic statistical analysis as indicated in SECTION 9
# 13. shuffle the final dataset and save it into a file for future use

# Your code starts here ...

