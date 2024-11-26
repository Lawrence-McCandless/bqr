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
x1 <- rbinom(n, 1, p=z); 
# generate confounding variable vector c
c <- matrix(NA, nrow=n, ncol=p); for (l in 1:p){c[,l] <- rbinom(n, 1, p=z)}

## Generate Y from a rounded Gaussian random variable response
y <- round(rnorm(n, mean=beta0.true + betax.true * x + c %*% rep(betac1.true, p), sd=sigmay.true * (1+gammx.true * x)), digits=0) ## rounded linear heteroscdasticity

## Alternatively, generate generate Y from a rounded Poisson random variable response
y <- rpois(n, lambda=exp(beta0.true + betax.true * x + c %*% rep(betac1.true, p))) ## rounded linear heteroscdasticity

## Conventional quantle regression
rq(y ~ x + c, tau=tau, method="br"); tp <- coef(summary(tmp)); print(tp)

## Bayesian quantile regression using statn
# Prepare teh dataset 
data.stan <- list(n=length(y), p=p+2, y=y, x=cbind(1,x,c), tau=tau) 

# MCMC sampling from posterior distribution
tmp <- sampling(b.aldsig, data=data.stan, chains=5, iter=4000, warmup=1000, thin=1, control=list(max_treedepth=7), verbose=TRUE) 
  
summary(tmp, pars=c("beta"))

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
