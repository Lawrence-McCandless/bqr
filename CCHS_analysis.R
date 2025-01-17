library(rstan); library(quantreg)
options(mc.cores = parallel::detectCores())

load('cchs.RData'); 

rq(PHQ9 ~ ccc140 + male + income + canadian + smk + alc + cannabis + hlthprovider + belonging + stressed, tau=0.5, data=cchs)


attach(cchs)
x <- cbind(1, ccc015, male, income, canadian, smk, alc, cannabis, hlthprovider, belonging, stressed)
detach(cchs)

data.stan <- list(n=dim(x)[1], p=dim(x)[2], y=cchs$PHQ9, x=x, tau=0.5) 


## Bayesian quantile regression using stan
# Prepare the dataset 
tau <- 0.5
data.stan <- list(
  n=dim(x)[1], 
  p=dim(x)[2], 
  y=cchs$PHQ9, 
  x=x, 
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
  
