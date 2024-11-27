library(rstan); library(quantreg)
options(mc.cores = parallel::detectCores())

n <- 1000;
betax1.true <- 0.5; # association between X1 and Y
betax2.true <- 0.5; # association between X2 and Y
sigmay.true <- 1; # variance of errors of Y

# generate z from a beta distribution for correlated binary predictor variables
z <-rbeta(n, 0.25, 0.25)
# generate exposure variable x
x1 <- rbinom(n, 1, p=z); 
x2 <- rbinom(n, 1, p=z); 

## Generate Y from a rounded Gaussian random variable response
y <- round(rnorm(n, mean=betax1.true * x1 + betax2.true * x2, sd=sigmay.true * (1+ 0.5 * x1)), digits=0) 

## Alternatively, generate generate Y from a zero-inflated Poisson random variable response
y <- rpois(n, lambda=exp(betax1.true * x1 + betax2.true * x2)) 
y[sample(1:n, 0.25*n)] <- 0

## Conventional quantle regression
rq(y ~ x1 + x2, tau=0.5, method="br")

## Bayesian quantile regression using stan
# Prepare the dataset 
tau <- 0.5
data.stan <- list(
  n=length(y), 
  p=3, 
  y=y, 
  x=cbind(1,x1,x2), 
  tau=tau
) 

# MCMC sampling from posterior distribution
results <- stan(
  file="bqr.stan",
  data=data.stan, 
  chains=5, 
  iter=4000, 
  warmup=1000, 
  control=list(max_treedepth=5), 
  verbose=TRUE
) 
  
summary(results, pars=c("beta"))

## MCMC diagnostics
plot(tmp, plotfun="trace")
plot(tmp, plotfun="hist")
  
