# Robust posterior credible intervals for model parameters using the sandwich adjustment to the 
# posterior covariance proposed by Yang, Wang and He (2016) International Statistical Review, 84, 327-344

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

## Alternatively, generate generate Y from a rounded Poisson random variable response
y <- rpois(n, lambda=exp(betax1.true * x1 + betax2.true * x2)) 
y[sample(1:n, 0.25*n)] <- 0

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

results <- summary(results, pars=c("beta"))
   
# Calculate maximum likelihood estimate of scale parameter
tmp.rq <- rq(data.stan$y ~ data.stan$x[,-1], tau=data.stan$tau, method="br") 
ae <- data.stan$tau*tmp.rq$residuals; 
ae[tmp.rq$residuals<0] <- (data.stan$tau-1)*(tmp.rq$residuals[tmp.rq$residuals<0])
sum.ae <- sum(ae)
sig.hat <- sum.ae/data.stan$n

tmp <- extract(results, permuted=TRUE, inc_warmup=FALSE);
D1.m1 <- var(tmp$beta)
D0 <- (t(data.stan$x) %*% data.stan$x)
yang.variance <- ((data.stan$tau*(1-data.stan$tau))/sig.hat^2) * ((D1.m1 %*% D0) %*% D1.m1)

alpha <- 0.05
summary(results)$summary[1:3,1] - qnorm(1-alpha/2)*sqrt(diag(yang.variance))
summary(results)$summary[1:3,1] + qnorm(1-alpha/2)*sqrt(diag(yang.variance))
