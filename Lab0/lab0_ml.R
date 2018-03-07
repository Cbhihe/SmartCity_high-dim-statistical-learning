####################################################################
# Machine Learning - MIRI Master
# Lluís A. Belanche

# LAB 0: Introduction to R (and a bit of statistics)
# version of February 2017
####################################################################

# Welcome to R!

##########################
## The R help system
##########################

# One of the first things you would like to do is getting help

# Given an R command or function x, entering ?x will bring up its help file

?exp

# The next command searches local matches for the string 'exponential' and provides a list of help files in which this term can be found in the concept or title (note that in R single '' and double "" quotes are the same)

help.search("exponential")

help.start()        # hypertext (HTML) version of R's online doc

RSiteSearch("mean") # broad search for 'mean' on http://search.r-project.org (will open in browser)

help.search("mean") 
??mean              # `??' is a shortcut for `help.search()'

help(mean)     # help on function mean in local search path
?mean          # shortcut for `help(mean)'

example(mean)  # run all the examples proposed in `help(mean)'

apropos("str") # everythng containing 'str' in its name


##########################
## The R environment
##########################

3-    1 # R will ignore these spaces

3 -     # the '+' prompt signals an incomplete command
1

# R handles special situations explicitly

-1/0; 0/0; 1/0; sqrt(-2+0i) 

# notice we use ';' to evaluate distinct expressions in a single command line

1:5         # The vector of numbers 1, 2, 3, 4, 5 
mean(1:5) 
sum(1:5)    

(2:5)^10    # 2 to the power of 10, 3 to the power of 10, ... a vector

# use concatenate 'c()' to glue objects

log2(c(0.5, 1, 2, 4, 8))

# R is case-sensitive, so foo, Foo and FOO are all considered different names

x <- 3
X
x

## Handling vectors/sequences

a <- c(1,3,2,4,5)

# the () forces to print the evaluation result

(a <- c(1,3,2,4,5))

c(0,2^2,mean(a),pi,"hello")

1:5                 # forward sequence with operator ':'
5:1                 # backward sequence with operator ':'
2*5:1               # operator ':' dominates other operators

seq(from=0,to=10,length=5)   # function 'seq()' generates 
seq(from=0,to=10,by=5)       # sorted sequences of given length/step

## Repeated sequences

rep(x=1:3,times=2)     # function 'rep()' generates repeated sequences
rep(x=1:3,times=1:3)
rep(c("ETH","ZURICH"),c(3,2))

paste("X",1:5,sep="")              # function 'paste()' and recycling of values
paste(c("Jan","Feb"),2010,sep="-")

gl(n=3,k=2,label=c("L","M","H"))   # function 'gl()' generates levels
                                   # (see 'factor()' below)

## Random sequences

set.seed(1234) # set seed of the random number generator

rnorm(n=5,mean=0,sd=2) # 100 draws from N(0,2) with 'rnorm()'
hist(rnorm(n=100,mean=0,sd=2))  # plot histogram of previous draws' distribution 

sample(x=1:10,size=5) # 'sample()' draws samples of specified size without replacing drawn items
sample(x=1:10,size=5, replace=TRUE)  # same but replacing drawn items.

sort(sample(x=1:10,size=10, replace=TRUE)) # 'sort()' sorts in ascending/descending order


##########################
## directory/save/delete
##########################

ls()                # list stored objects in current R session
ls(pattern="^x")    # list objects starting with 'x'

getwd()              # display current working directory
setwd(dir="C:/")     # set current working directory at 'C:\'
                     # use '\\' or '/' instead of '\' in paths

save(m,file="M.Rdata")         # save R object m in file M.Rdata
save.image("myfile.Rdata")     # save current workspace in myfile.RData
savehistory("myhist.Rhistory") # save command history in myhist.Rhistory

load("myfile.Rdata");ls()      # load myfile.RData's content into workspace
loadhistory("myhist.Rhistory") # load command history from hist

write("print('Hello!')","hello.txt") # write 'print(Hello!)' in file "hello.txt"
source("hello.txt")  # execute all R commands in file "hello.txt"

##########################
## Vectors
##########################

vector(mode = "logical",  # 'vector()' generates a vector of given mode
length = 3)               #  and length with default components

numeric(10);logical(3)   # 'numeric()', logical(), etc...
                         # generate logical/logical vectors of specified length

x <- as.vector(c(2,4,8)); x; mode(x)  # 'as.vector()' transforms
                                      # concatenated object into vector
                                      # 'mode()' displays the vector's mode

names(x) <- c("banana", "apple","nut") # 'names()' sets components' name
names(x); is.vector(x); length(x)        # 'names()' also displays the
                                         # name/labels of the vector's elements
                                         # 'is.vector()' tests vector's identity,
                                         # 'length()' displays vector's length
 
x.old <- x

names(x) <- NULL          # erase the names of x's components

length(x) <- 5            # set the length of x to 5; R fills in the holes with NAs
                          # NAs represent "missing values"

x

mode(x)                      # NAs do not alter vector's mode

mode(x) <- "character"    # change numeric vector to
                             # character (note that NAs are not converted)

x

x <- as.numeric(x);length(x) <- 3  # change character vector to
                                  # numeric and shorten vector to length 3

x

x[2]; x[-2]; x[c(1,3)]; x[rep(1:3,2)]   # indexing using integers

x[c(TRUE,TRUE,FALSE)]      # index using logical values TRUE/FALSE
x[(x>=2)&(x<=6)]

x.old[c("nut","banana")]    # index vector of character strings
                           
x.old[c("nut","banana")] <- c(9,0); # modifies these components

x.old

##########################
## Factors
##########################

# factors are used to implement categorical (nominal/ordinal) values

gender <- rep(c("Female","Male"),c(3,2)); gender; mode(gender) # generate a character vector

levels(gender)  # 'levels()' displays factor levels
                # vectors have no levels

gender <- factor(gender); gender # 'factor()' converts vectors into factors

as.character(gender)             # back to the character vector

cont <- cut(1:9,breaks=3); cont # 'cut()' codes the components of the
                                # sequence according to specified breaks
                                # and then creates a factor
 
# much better
(cont <- cut(1:9,breaks=3,labels=c("low","medium","high")))


####################################################
# Matrices
####################################################

A <- matrix(1:16,4,4)
A
A[2,3]
A[c(1,3),c(2,4)]
A[1:3,2:4]
A[1:2,]
A[,1:2]
A[1,]
A[-c(1,3),]
A[-c(1,3),-c(1,3,4)]
dim(A)

(B <- matrix(1:6,nrow=2,ncol=3,byrow=TRUE))
(B <- matrix(1:6,nrow=2,ncol=3,byrow=FALSE))

B[-(2:3), ]

B[-(2:3), 2]

rbind(1:3,4:6)        # function 'rbind()' for binding rows

cbind(1:2,3,4)        # function 'cbind()' with recycling

array(1:6,dim=c(2,3)) # function 'array()'

V <- 1:6; dim(V) <- c(2,3); V   # function 'dim()' sets/retrieves dimensions

V2 <- V     # copy 'V' in 'V2'

dimnames(V)   # function 'dimnames()' sets/rerieves row/col names

dimnames(V2) <- list(NULL, c("C1","C2","C3")); V2

V[-2,c(1,3)]; V2[-2,c("C1","C3")]

X <- matrix(1:4,ncol=2)   # define X, a square matrix
X
Y <- matrix(4:1,ncol=2)   # define Y, a square matrix
Y 
X + Y       # element-wise addition (the dimensions must agree) 
X * Y       # element-wise multiplication (the dimensions must agree) 

X %*% X                     # standard matrix multiplication operator

apply(X,MARGIN=1,FUN=sum) # function 'apply()' applies function FUN
                          # to specified MARGIN ('1' means rows)

t(X)                      # function 't()' transposes a matrix

solve(X) # solve linear system of the form X %*% z = b for z
         # if b is not specified, compute the inverse of X

solve(X,matrix(c(1,0,0,1),ncol=2)) %*% X # check that 1st term calculate X^(-1)

# Solve X %*% z = Y for z (X is n-by-k, z k-by-m, Y n-by-m)
solve(X,Y)

eigen(X)            # eigenvalues and eigenvectors of X

svd(X)      # Singular Value Decomposition of X, 
            # where d = singular values,
            # u = left singular vector, v = right singular vector 
svd(X)$u %*% matrix(c(svd(X)$d[1],0,0,svd(X)$d[2]),ncol=2) %*% t(svd(X)$v) 

######################################
## Data frames
######################################

# Data frames are matrices, but they allow named columns and non-numerical data

mystuff <- data.frame(name=c(
    "John",
    "Robert",
    "William",
    "Charles",
    "David",
    "James",
    "Richard",
    "Johann",
    "George",
    "Paul"),
           age=c(11,10,18,6,5,6,7,8,9,2),
           country=c("US","UK","USA","USA","NZ","USA","USA","IRE","USA","CAN"))

mystuff

summary(mystuff)

names(mystuff) # names of the variables in the data frame

dim(mystuff)  # dimensions of the data frame 

str(mystuff)  # 'str()' gives the structure of R objects

mystuff[,2]; mystuff[,2][1:3]   # index vectors on the row and col dimensions

mystuff$age; mystuff$age[1:3]   # extract components with operator `$' and names

mystuff[[2]]; mystuff[[2]][1:3] # operator `[[i]]' extracts the ith component
 
mystuff[[2]]                   # '[[.]]' drops names from 2nd component
 
mystuff[2]                     # '[.]' keep names of the 2nd component

age

attach(mystuff)  # make columns visible at top level

age; age[1:3]

age <- age/10

mystuff$age            # WARNING! changes in attached objects do not modify
                       # the source object

mystuff$age <- mystuff$age + 1 # this works

mystuff

detach(mystuff)                # detach dataframe (undoes the 'attach')

# rows from 'mystuff' displayed in decreasing order with respect to 'age'
mystuff[order(-mystuff$age),]    

# Data frames can be subset

(USAnames <- subset(mystuff, country == "USA", select = name))

# be careful: USAnames is still a column with the same possible values
summary(USAnames)

# if you would like to 'fix' this, we must re-factor
USAnames$name <- factor(USAnames$name)

summary(USAnames)


## Managing missing values: example 1
(Z <- data.frame(x = c(1, 2, 3), y = c(0, 10, NA)))

na.omit(Z)   # function 'na.omit()' removes rows with NAs
             # see also functions 'na.fail()', 'na.exclude()', 'na.pass()'.

(z <- c(1,-99,-999,"."))

(z[is.element(z,c(-99,-999,"."))] <- NA)   # converts '-99' and  '-999'
                                           # and '.' into NAs
(z <- as.numeric(z))                       # converts 'z' into numeric

z[!is.na(z)]                               # z vector without NAs

## Managing missing values: example 2

library(MASS) # there you can find the Pima Indians data set, full of missing values

table(is.na(Pima.tr2$bp))   # 13 NAs out of 300 rows
sapply(Pima.tr2, function(x) sum(is.na(x))) # guess what this does?


#####################################
## Lists
#####################################

mylist <- list(Day=c("Mo","Tu","Su"),   # a list of R objects
               Place="Sabadell",eigen(X),obj4=mystuff)

names(mylist)               # names of the objects in 'mylist'

mylist$Day                  # display the object Day in 'mylist'
mylist$Day[2]               # display the second element from
                            # 'Day' in 'mylist'
mylist$obj4$age             # shows element 'age' from 'obj4' in 'mylist'

mylist[[3]]$values          # the i-th component can also be accessed
                            # with the double squared brackets `[[i]]' 

# check the difference between `[[.]]' and `[.]'
mylist[[1]]

mylist[1]



#####################################
# Control structures
#####################################

# Let's do some iteration and conditioning

for (i in 1:10) 
  if ((g<-rnorm(1))>0.25) 
  { cat(paste(i, round(g,4),"\n")) } else { cat("no!\n") }


#####################################
# Function definitions
#####################################

# Example 1

my.sqr <- function(x) x^2

my.sqr(1:5)

# Example 2

drawCards = function (number)
{
  faceValues <- c(1:10, "Jack", "Queen", "King")
  suits <- c("Hearts", "Diamonds", "Spades", "Clubs")
  f <- sample(faceValues, number, replace=TRUE)
  s <- sample(suits, number, replace=TRUE)
  paste(f, "of", s)
}

drawCards(5)

# Note that the function returns as value the result of its last evaluation

# Example 3

# Plotting a variable against the normal pdf

hist.with.normal <- function (x, xlabel=deparse(substitute(x)), ...)
{
  h <- hist(x,plot=F, ...)
  s <- sd(x)
  m <- mean(x)
  ylim <- range(0,h$density,dnorm(0,sd=s))
  hist(x,freq=FALSE,ylim=ylim,xlab=xlabel, main="", ...)
  curve(dnorm(x,m,s),add=TRUE, col = "violet")
}

hist.with.normal (mtcars$mpg)

# Note the advanced ways of passing parameters:

# deparse(substitute(x)) is a trick to get the parameter's name
# parameters are optional, and can be given a default value
# ... means additional (unspecified) parameters

#####################################
## Importing/exporting .csv files
#####################################

    #create a data frame
    dates <- c("3/27/1995", "4/3/1995",
               "4/10/1995", "4/18/1995")
    prices <- c(11.1, 7.9, 1.9, 7.3)
    d <- data.frame(dates=dates, prices=prices)

    #create the .csv file
    dir.create ("labo")
    filename <- "labo/temp.csv"
    write.table(d, file = filename, sep = ",", row.names = FALSE)

    #read the .csv file
    read.table(file = filename, sep = ",", header = TRUE)
    read.csv(file = filename) #same thing


#####################################
# Processing a data file
#####################################

# Consider the 'trees' data set, providing measurements of the girth, height and volume of timber in 31 felled black cherry trees.

 trees

# R can provide very quick summary information

 summary(trees)
 names(trees)

 boxplot(trees)
 pairs(trees)
 hist(trees$Height)

 hist.with.normal (trees$Volume)

 plot(trees$Height, trees$Volume)
 
 par(mfrow=c(3,1))    # # sets a 1x1 display

 for (i in 1:3) 
 { plot(density(trees[,i]), xlab="Value", main=names(trees)[i]) }

# Can we model Volume from Girth and Height?

# 2D plot

par(mfrow=c(1,1))  # sets old display

# elegant way of accessing variables in a dataframe
with(trees, plot(Girth,Volume,main = "Black cherry trees growth", pch = 17))

# add centered axes

with(trees, abline(v = mean(Girth), col = "blue"))
with(trees, abline(h = mean(Volume), col = "blue"))

# add regression line

with(trees, abline( lm(Volume ~ Girth), col = "red" ))

# Individual variables can be directly named if we 'attach' the data frame to the current workspace

 attach(trees)
 Girth

# if we attach() we do not need to specify the dataframe all the time 
# (handle this with care, it may mask errors when updating attached information)

# now 4 plots in one 2x2 device

par(mfrow = c(2,2))              
plot(Girth,Height)                 

# these are specially useful for spatial data, but can be applied:

image(kde2d(Girth,Height))         #colored image
persp(kde2d(Girth,Height))         #perspective plot
contour(kde2d(Girth,Height))       #contour plot

library(rgl)

# and a 3d scatterplot (that can be rotated with the mouse)
plot3d(Girth,Height,Volume, size = 5, col='green')

# We can make more complicated (though still linear) models

 (model2 <- lm(log(Volume) ~ Girth + I(Girth^2)))

 (model3 <- lm(log(Volume) ~ poly(Girth, degree = 3)))

# We can get (or set) the names of an object

 names(model2)

# and its attributes

 attributes(model2)
 model2$fitted.values

# or display the internal structure of an object

 str(trees)
 str(model2)
 unclass(model2)

# though in general summary() will be enough

 summary(model3)
 detach(trees)


#########################
# Creating random samples
#########################

# Let's create different iid random samples: Normal (Gaussian), uniform, Chi-square and Weibull; we show a histogram (left) and a density estimation (right) against the normal pdf

N <- 10000

par(mfrow=c(2,3))

# Normal (Gaussian) distribution

nd1 <- rnorm(N, mean = 0, sd = 1)         
nd2 <- rnorm(N, mean = 2, sd = 5)        
nd3 <- rnorm(N, mean = -10, sd = 0.01)       

hist.with.normal (nd1, "Standard Normal")
hist.with.normal (nd2, "Normal N(2,5)")
hist.with.normal (nd3, "Normal N(-10,0.01)")

plot(density(nd1), col='red', main="estimated density")
plot(density(nd2), col='red', main="estimated density")
plot(density(nd2), col='red', main="estimated density")

# Chi-square distribution

cq1 <- rchisq(N, 1)                  
cq2 <- rchisq(N, 2)                  
cq3 <- rchisq(N, 10)

hist.with.normal (cq1, "Chi-square (df=1)")
hist.with.normal (cq2, "Chi-square (df=2)")
hist.with.normal (cq3, "Chi-square (df=10)")

plot(density(cq1), col='blue', main="estimated density")
plot(density(cq2), col='blue', main="estimated density")
plot(density(cq2), col='blue', main="estimated density")

# Weibull distribution

wd1 <- rweibull(N, shape = 2, scale = 1)
wd2 <- rweibull(N, shape = 1, scale = pi)   
wd3 <- rweibull(N, shape = 5, scale = 0.5)   

hist.with.normal (wd1, "Weibull (2,1)")
hist.with.normal (wd2, "Weibull (1,pi)")
hist.with.normal (wd3, "Weibull (5,0.5)")

plot(density(wd1), col='darkcyan', main="estimated density")
plot(density(wd2), col='darkcyan', main="estimated density")
plot(density(wd3), col='darkcyan', main="estimated density")

# the boxplot is a very useful 1D tool

boxplot(nd1)                          
boxplot(nd3)                          
boxplot(cq1)                          
boxplot(cq3)                          
boxplot(wd1)                          
boxplot(wd3)                          

# various 2D plots

par(mfrow=c(2,2))

plot(nd1, nd1^2)                          
plot(nd2^2, cq2)                          

plot(nd2*wd2, wd1/wd3, xlab = "HORIZONTAL Title", ylab = "VERTICAL Title")
plot(nd2,sin(nd2), main = "nice sine plot", pch = 19, col="tomato2")

#####################
## Plots and graphics
#####################

par(mfrow=c(1,1))

x <- rnorm(100)
y <- rnorm(100)

plot(x,y,xlab="this is the x-axis",ylab="this is the y-axis",main="Plot of X vs Y")

pdf("Figure.pdf")

plot(x,y,col="green")

dev.off()

x <- seq(-pi,pi,length=50)
y <- x

f <- outer(x,y,function(x,y)cos(y)/(1+x^2))
contour(x,y,f)
contour(x,y,f,nlevels=45,add=T)
fa <- (f-t(f))/2
contour(x,y,fa,nlevels=15)
image(x,y,fa)
persp(x,y,fa)
persp(x,y,fa,theta=30)
persp(x,y,fa,theta=30,phi=20)
persp(x,y,fa,theta=30,phi=70)
persp(x,y,fa,theta=30,phi=40)


## Still more plots ...

par(mfrow=c(2,3)) # create a 2x3 plotting device
set.seed(12); x <- rnorm(10)    # random draw from N(0,1)

for (i in c("p","l","b","c","h","s")) # loop through plot types
plot(x=x,type=i, main=paste("plot(type=",i,")",sep=""))

## Simple 2D plots: hist(), boxplot()

par(mfrow=c(1,3))       # parametrize a 1x3 plotting device

plot(x=-5:5,y=(-5:5)^2, # simple scatter plot of y=x^2
main=expression(y==x^2))  # with math expression in the title

hist(x,freq=FALSE)          # histogram with prob. densities
lines(density(x),col="red") # add red kernel density to histogram
rug(x)                      # add rug representation to histogram

boxplot(x,rt(10,3),rchisq(10,3), # boxplot
names=c("Normal","Student","ChiSquare"))  # series labels
title("Box-plot")                # add title to boxplot

graphics.off()      # reset/close all graphical devices 

## more 3D plots: persp(), contour()

x <- seq(-10,10, length=30); y <- x       # x-y dimensions
f <- function(x, y)                       # z dimension
  {r <- sqrt(x^2+y^2)
   10 * sin(r)/r + 3}

z <- outer(x, y, f)   # apply f to "outer product" of x and y

persp(x, y, z, zlim=range(0,max(z)),      # 3D plot
theta=30, phi=30, expand=0.5, col="lightblue")

contour(x, y, z) # contour plot


############################
## Packages/Libraries
############################

library()              # list of all installed packages
 
install.packages("np") # install package called np
install.packages(pkgs="D:/DAAG_0.91.zip", repos=NULL) 
update.packages()      # update installed packages


# library(package) and require(package) both load the package with name 'package' in the local library file

library(np)            # load package np in the current session
require(np)            # require returns FALSE if the package does not exist and goes on
                       # these functions do not reload a package which is already loaded.
 
                       # once the desired has been downloaded
?"np"                  # provide help on package np
vignette("np")         # provide short intro to package np

# For a listing of all routines in the np package type: 

library(help="np")

############################
## Computations with big datasets
############################

M <- matrix(rnorm(5e6),ncol=50) 
dim(M) 

## Times are: user cpu, system cpu, elapsed

system.time(M+1) 

# see the difference!

M.df <- data.frame(M) 

system.time(M.df+1) 

# see the difference!

M <- matrix(rnorm(2e6), nrow=2000)
dim(M)

system.time(apply(M,1,sum))             # row sums
system.time(M %*% rep(1,1000))           # row sums
system.time(rowSums(M))                 # row sums
