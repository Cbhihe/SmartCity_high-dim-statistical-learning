####################################################################
# Machine Learning - MIRI Master
# Lluís A. Belanche

# LAB 2: Bias/Variance and VC-dimension
# version of February 2017
####################################################################

####################################################################
#  Underfitting/Overfitting and Bias/Variance
####################################################################

set.seed(2222)
par(mfrow=c(1, 1))

## First we are going to do a simple modeling in one dimension with polynomials
## of different degrees: we are going to observe the underfitting/overfitting phenomenon

## Suppose we have a process which generates data in the form of a sine + some noise \varepsilon:

## t = sin(2·pi·x) + varepsilon, where varepsilon ~ N(0, 0.04)

## the function f(x) = sin(2·pi·x) is the best possible solution (the regression function)
## we therefore aim at getting close to this function using only a limited sample of N data pairs (x,t)

## If we can assume that our target function is continuous, then the polynomials are a good choice
## for our function class:

## We start with a simple linear model: y(x) = beta_0 + beta_1·x

## and progressively add more higher terms:

##      y(x) = beta_0 + beta_1·x + beta_2·x^2 + ... + beta_M·x^M 

## These additional terms will improve the fit to the training data (reduce the empirical error), 
## but at the expense of an steady increase in complexity

library(PolynomF)

N <- 10

f <- function(x) sin(2 * pi * x)  # this is the (unknown) target function and the best solution to our problem

x <- runif(N)                     # generate the x according to a uniform distribution p(x)

t <- f(x) + rnorm(N, sd=0.2)      # generate the t according to a gaussian conditional distribution p(t|x)
                                  # note stdev (sd) = sqrt(var) = sqrt(0.04) = 0.2

## Plot the available data sample (what we observe)
plot(data.frame(x, t),ylim=c(-1.25,1.25))

## Plot the deterministic part (unknown in practice, and the best possible solution = the regression function)
curve(f, type="l", col="green", add=TRUE)

## Fitting polynomials is easy in R with the poly() function
## this is a 1-degree polynomial (a line)
polyfit <- lm(t ~ poly(x, 1, raw=TRUE))

## We will be analysing this summary output next session ... though probably you can understand many things
summary(polyfit)

## Today we are not so interested in how to optimize the model *coefficients*, but in
## optimizing the *degree* and, particularly, in the effect of changing the degree on the modeling

## Let us see how our fit is:

p <- polynom(coef(polyfit))                           # this is the polynomial of degree 1 that we just fit
p
points(data.frame(x, predict(polyfit)), col="red")    # and these are its predictions on same data
curve(p, col="red", add=TRUE)                         # and the polynomial as a function (i.e., a model)

## Now let us fit a bunch of polynomials of different degrees
par(mfrow=c(2, 3))

for (i in c(1, 2, 3, 4, 6, 9)) 
{
  plot(data.frame(x, t), xlab=paste("Polynomial fit of degree", i), ylab="f(x)",ylim=c(-1.25,1.25))
  curve(f, type="l", col="green", add=TRUE)
  polyfit <- lm(t ~ poly(x, i, raw=TRUE))
  p <- polynom(coef(polyfit))
  curve(p, col="red", add=TRUE)
}

## We can see that adding more polynomial terms improves the fit. The limit case is for the 9th polynomial
## (degree 9, having 10 coefficients) which is able to interpolate the data (it passes through every data point)

## So which one is better?

## In terms of training error, the 9th polynomial (because it shows 0 error)
## In terms of complexity, the 1th polynomial (because it has less coefficients)

## The best is the one with a better trade-off between training error and complexity

## We have seen in class that there are at least two ways of characterizing this trade-off: 
## Vapnik's theorems (using VC-dimension) and Bias/Variance analysis (due to Geman & Geman)

## Let's do a Bias/Variance analysis first (notice that, without the green line, we still do not know which polynomial is underfitting/overfitting the data)

## To show Bias, let us concentrate on the 1th polynomial:
## we generate 10 independent datasets and model as before keeping the degree constant to 1

par(mfrow=c(2, 5))

for (i in 1:10)
{
  x <- runif(N)                 
  t <- f(x) + rnorm(N, sd=0.2)
  
  plot(data.frame(x, t), xlab=paste("Polynomial fit of degree 1, data sample", i), ylab="f(x)",xlim=c(0,1),ylim=c(-1.25,1.25))
  curve(f, type="l", col="green", add=TRUE)
  polyfit <- lm(t ~ poly(x, 1, raw=TRUE))
  p <- polynom(coef(polyfit))
  curve(p, col="red", add=TRUE)
}

## We can see two things:
## 1) the model is quite stable (therefore it has low Variance)
## 2) the model is quite bad on average (therefore it has high Bias)

## To show Variance, let us concentrate on the 6th polynomial:
## we generate 10 independent datasets and model as before keeping the degree constant to 6

par(mfrow=c(2, 5))

for (i in 1:10)
{
  x <- runif(N)                 
  t <- f(x) + rnorm(N, sd=0.2)
  
  plot(data.frame(x, t), xlab=paste("Polynomial fit of degree 6, data sample", i), ylab="f(x)",xlim=c(0,1),ylim=c(-1.25,1.25))
  curve(f, type="l", col="green", add=TRUE)
  polyfit <- lm(t ~ poly(x, 6, raw=TRUE))
  p <- polynom(coef(polyfit))
  curve(p, col="red", add=TRUE)
}

## We can see two things:
## 1) the model is quite unstable (therefore it has high Variance)
## 2) the model is quite good *on average* (therefore it has low Bias)

## In practice, we cannot generate many different data samples ... so we cannot compute Bias and Variance.

## What are we going to do? How can we avoid overfitting/underfitting? (note that most often the real danger is in overfitting; this is because many ML methods tend to be very flexible, i.e., they are able to represent complex models)

## There are several ways to control this:
  
## 1) Get more training data (typically out of our grasp)
## 2) Use (that is, sacrifice!) part of the data for validation
## 3) Use an explicit complexity control (aka regularization)

## Let's see the effect of more training data first
## Let us fit a 9th-polynomial for different data sizes:

par(mfrow=c(2, 2))

for (N in c(10, 25, 50, 100)) 
{
  x <- runif(N)                 
  t <- f(x) + rnorm(N, sd=0.2)
  
  plot(data.frame(x, t), xlab=paste("Training data size:", N), ylab="f(x)",xlim=c(0,1),ylim=c(-1.25,1.25))
  curve(f, type="l", col="green", add=TRUE)
  polyfit <- lm(t ~ poly(x, 6, raw=TRUE))
  p <- polynom(coef(polyfit))
  curve(p, col="red", add=TRUE)
}

## So getting more data increasingly reduces the chances of overfitting for the same model
## More technically, the chance of overfitting is "moved ahead"; this means that
## as we get more data, only very complex models will have a possibility to overfit

## Now let us use part of the data for validation (this is what cross-validation does, although in a more systematic way: we will see this in a future lab session)

## We are going to generate the same plot as before, but this time we will be computing the prediction error on the validation data, which we create as an independent (and larger) sample

## generate training data again
N <- 10
x <- runif(N)                     # generate the x according to a uniform distribution p(x)
t <- f(x) + rnorm(N, sd=0.3)      # generate the t according to a gaussian conditional distribution p(t|x)

par(mfrow=c(2, 3))

## generate validation data (note generation mechanism MUST be the same)
N.val <- 1000
x.val <- runif(N.val)                     # generate the x according to a uniform distribution p(x)
t.val <- f(x.val) + rnorm(N.val, sd=0.3)      # generate the t according to a gaussian conditional distribution p(t|x)
val <- data.frame(x=x.val, t=t.val)

errors <- matrix (nrow=6, ncol=3)
colnames(errors) <- c("Degree","TR.NRMSE","VA.NRMSE")
degrees <- c(1, 2, 3, 4, 6, 9)

## Same loop as before, but this time we keep track of training *and* validation error 
## notice we plot true f (in green), along with validation data

## This time we switch to a much better assessment of error: the NRMSE, which is a root normalized MSE

for (i in 1:length(degrees))
{
  plot(data.frame(x.val, t.val), xlab=paste("Polynomial prediction of degree", degrees[i]), ylab="f(x)",ylim=c(-1.25,1.25), col="black")
  curve(f, type="l", col="green", add=TRUE)
  polyfit <- lm(t ~ poly(x, degrees[i], raw=TRUE))
  
  # fill in degree, training error and validation error (both are NRMSEs)
  errors[i, "Degree"] <- degrees[i]
  errors[i, "TR.NRMSE"] <- sqrt( sum(polyfit$residuals^2) / ((N-1)*var(t)) )
                                
  predictions <- predict(polyfit, newdata=val)
  errors[i, "VA.NRMSE"] <- sqrt( sum((t.val - predictions)^2) / ((N.val-1)*var(t.val)) )
  
  points(data.frame(x.val, predict(polyfit, newdata=val)), col="red") # these are the predictions on validation data
}

errors

## Inspect the results. Which degree yields the lowest error on validation data? There is a close difference between degrees 3 and 4; typically that for degree 3 is a bit smaller

## Notice how the error on training data keeps decreasing
## and at the same time complexity (indicated by the degree) keeps increasing

## Sadly, this method (known as the 'holdout') does not always work so well
## Moreover, both training and validation data samples should be large at the same time ...
## On top of that, we depend on fortunate/unfortunate data splits
## We will return to this topic ('resampling methods') in a future lab session 

## The third method (using explicit complexity control, a.k.a. regularization) will actually be the topic of the next lab session, so we defer this analysis for then


####################################################################
# The VC dimension: example 1
####################################################################

# A function of a single parameter having infinite VC dimension

N <- 10 # this could be any natural number
x <- 10^-seq(1,N,1) # this is our choice of inputs

t <- c(-1, +1, -1, +1, -1, -1, +1, -1, +1, -1) # this could be any of the 2^N choices of labels

(alpha <- pi*(1+sum((1-t)*10^seq(1,N,1))/2))


y <- sign(sin(x*alpha))

# Indeed it works ...
N == length(t==y)

# I have not found a satisfactory way of plotting the result ...
# Moreover, the numerical precision required for larger values of N is huge (so this is really a theoretical example)

####################################################################
# The VC dimension: example 2
####################################################################

## Let's use VC dimension for Structural risk minimization (SRM); we play again with polynomial fitting

par(mfrow=c(1, 1))
set.seed(1234)

N <- 25

f <- function(x) sin(2 * pi * x)  # this is the (unknown) target function and the best solution to our problem

x <- runif(N)                     # generate the x according to a uniform distribution p(x)

t <- f(x) + rnorm(N, sd=0.2)      # generate the t according to a gaussian conditional distribution p(t|x)
# note stdev (sd) = sqrt(var) = sqrt(0.04) = 0.2

## Plot the available data sample (what we observe)
plot(data.frame(x, t),ylim=c(-1.25,1.25))

## Plot the deterministic part (unknown in practice, and the best possible solution = the regression function)
curve(f, type="l", col="green", add=TRUE)

VCdim <- function (d,k) { choose(d+k,k) }   # as seen in class

# as seen in class
H <- function (N,h,eta) { sqrt((h*(log(2*N/h)+1) - log(eta/4))/N) }

eta <- 0.01

## Now let us fit a bunch of polynomials of different degrees; for each of them we compute the NRMSE in training, the VC-dimension and the bound

results <- matrix (nrow=9, ncol=5)
colnames(results) <- c("Degree","TR error","VC dim", "complexity", "bound")
degrees <- seq(9)

for (i in 1:length(degrees)) 
{
  polyfit <- lm(t ~ poly(x, degrees[i], raw=TRUE))

  # fill in degree, training error and VC-dimension
  results[i, "Degree"] <- degrees[i]
  results[i, "TR error"] <- sqrt( sum(polyfit$residuals^2) / ((N-1)*var(t)) )
  results[i, "VC dim"] <- VCdim (1,degrees[i])
  results[i, "complexity"] <- H(N,results[i, "VC dim"],eta)
  results[i, "bound"] <- results[i, "TR error"] + H(N,results[i, "VC dim"],eta)
}

results

## It sort of works ... the SRM principle says that we should choose a polynomial of degree 3 because it minimizes the right hand side of the risk bound

## So we have used another decomposition of the prediction error to choose the optimal model