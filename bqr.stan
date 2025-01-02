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
