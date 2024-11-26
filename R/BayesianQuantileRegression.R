library(rstan); library(quantreg)
options(mc.cores = parallel::detectCores())

n <- 1000;
betax1.true <- 0.5; # association between X1 and Y
betax2.true <- 0.5; # association between X2 and Y
sigmay.true <- 1; # variance of errors of Y

tau <- 0.5

# generate z from a beta distribution for correlated binary predictor variables
z <-rbeta(n, beta.cor, beta.cor)
# generate exposure variable x
x1 <- rbinom(n, 1, p=z); 
x2 <- rbinom(n, 1, p=z); 

## Generate Y from a rounded Gaussian random variable response
y <- round(rnorm(n, mean=betax1.true * x1 + betax2.true * x2, sd=sigmay.true * (1+ 0.5 * x1)), digits=0) ## Model 1

## Alternatively, generate generate Y from a rounded Poisson random variable response
y <- rpois(n, lambda=exp(betax1.true * x1 + betax2.true * x2)) ## Model 2

## Conventional quantle regression
rq(y ~ x1 + x2, tau=tau, method="br")

## Bayesian quantile regression using statn
# Prepare teh dataset 
data.stan <- list(
  n=length(y), 
  p=3, 
  y=y, 
  x=cbind(1,x1,x2), 
  tau=tau
) 

# MCMC sampling from posterior distribution
results <- stan(
  file="bqr.stan"
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
  
## Sandwich adjusted posterior variance  
# Calculate estimate of scale paraemter
tmp.rq <- rq(data.stan$y ~ data.stan$x[,-1], tau=data.stan$tau, method="br") 
ae <- data.stan$tau*tmp.rq$residuals; 
ae[tmp.rq$residuals<0] <- (data.stan$tau-1)*(tmp.rq$residuals[tmp.rq$residuals<0])
sum.ae <- sum(ae)
sig.hat <- sum.ae/data.stan$n

# Calcuate corrected posterior covariance matrix
tmp3 <- extract(tmp, permuted=TRUE, inc_warmup=FALSE);
D1.m1 <- var(tmp3$beta)/1
D0 <- (t(data.stan$x) %*% data.stan$x)
