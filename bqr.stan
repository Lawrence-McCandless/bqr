functions{
real ald_lpdf(real x, real u, real sig, real tau){ // user-defined log pdf of asymetric Laplace distribution;
  real y;
  if (x<u)
    y = log((tau*(1-tau)/sig)*exp(-(tau-1)*(x-u)/sig));
  if (x>=u)
    y = log((tau*(1-tau)/sig)*exp(-tau*(x-u)/sig)); 
  return y;
}
}
data {
  int<lower=0> n;
  int<lower=0> p;
  vector[n] y; // outcome variable;
  array[n] row_vector[p] x; // predictor variables;
  real<lower=0> tau; // quantile;
  
}
parameters {
  vector[p] beta; // regression coefficients;
  real sig; // scale parameter; 
}
model {
  for (i in 1:n){
  real eta; 
  eta =  x[i] * beta;
  y[i] ~ ald(eta, sig, tau); 
  }
 }
