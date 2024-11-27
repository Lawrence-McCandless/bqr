# Bayesian method to assess the presence of heteroscedastic errors of the outcome variable

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
data.stan <- list(
  n=length(y), 
  p=3, 
  y=y, 
  x=cbind(1,x1,x2)
) 

tau1 <- 0.5; tau2 <- 0.75 ## Two values of tau for assessing heteroscedasticity

# MCMC sampling from posterior distribution
data.stan$tau <- tau1
results1 <- stan(
  file="bqr.stan",
  data=data.stan, 
  chains=5, 
  iter=4000, 
  warmup=1000, 
  control=list(max_treedepth=5), 
  verbose=TRUE
) 

data.stan$tau <- tau2
results2 <- stan(
  file="bqr.stan",
  data=data.stan, 
  chains=5, 
  iter=4000, 
  warmup=1000, 
  control=list(max_treedepth=5), 
  verbose=TRUE
) 

# Calculate maximum likelihood estimate of scale parameter for tau1
tmp.rq <- rq(data.stan$y ~ data.stan$x[,-1], tau=tau1, method="br") 
ae <- tau1*tmp.rq$residuals; 
ae[tmp.rq$residuals<0] <- (data.stan$tau-1)*(tmp.rq$residuals[tmp.rq$residuals<0])
sum.ae <- sum(ae)
sig.hat1 <- sum.ae/data.stan$n

# Calculate maximum likelihood estimate of scale parameter for tau2
tmp.rq <- rq(data.stan$y ~ data.stan$x[,-1], tau=tau2, method="br") 
ae <- tau2*tmp.rq$residuals; 
ae[tmp.rq$residuals<0] <- (data.stan$tau-1)*(tmp.rq$residuals[tmp.rq$residuals<0])
sum.ae <- sum(ae)
sig.hat2 <- sum.ae/data.stan$n

D0 <- (t(data.stan$x) %*% data.stan$x)   
 
sumr1 <- summary(results1, pars=c("beta"))
tmp1 <- extract(results1, permuted=TRUE, inc_warmup=FALSE)
D1.m1.1 <- var(tmp1$beta)
  
sumr2 <- summary(results2, pars=c("beta"))
tmp2 <- extract(results2, permuted=TRUE, inc_warmup=FALSE)
D1.m1.2 <- var(tmp2$beta)

yang.variance1 <- ((tau1*(1-tau1))/sig.hat1^2) * ((D1.m1.1 %*% D0) %*% D1.m1.1)
yang.variance2 <- ((tau2*(1-tau2))/sig.hat2^2) * ((D1.m1.2 %*% D0) %*% D1.m1.2)
yang.covariance <- ((min(tau1, tau2)-tau1*tau2)/(sig.hat1 * sig.hat2)) * ((D1.m1.1 %*% D0) %*% D1.m1.2)
  
covar <- rbind(cbind(yang.variance1,yang.covariance), cbind(yang.covariance,yang.variance2)) 
R <- matrix(0, ncol=1, nrow=2*dim(D0)[1]); R[c(2,2+dim(D0)[1])] <- c(-1,1); 
se <- sqrt((t(R) %*% covar) %*% R) 
pnorm(abs((summary(results1)$summary[2,1] - summary(results2)$summary[2,1])/se)) 


