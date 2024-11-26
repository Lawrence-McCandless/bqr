functions{
real ald_lpdf(real x, real u, real sig, real tau){
  real y;
  if (x<u)
    y = log((tau*(1-tau)/sig)*exp(-(tau-1)*(x-u)/sig));// for formula see page 4 of OCR notes;
  if (x>=u)
    y = log((tau*(1-tau)/sig)*exp(-tau*(x-u)/sig)); 
  return y;
}
}
data {
  int<lower=0> n;
  int<lower=0> p;
  vector[n] y;
  array[n] row_vector[p] x;
  real<lower=0> tau;
  
}
parameters {
  vector[p] beta;
  real sig;
}
model {
  for (i in 1:n){
  real eta; 
  eta =  x[i] * beta;
  y[i] ~ ald(eta, sig, tau); 
  }
 }
