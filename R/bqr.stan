data {
  int<lower=0> N;         // Number of observations
  vector[N] x;            // Predictor
  vector[N] y;            // Response
}

parameters {
  real beta0;             // Intercept
  real beta1;             // Slope
  real<lower=0> sigma;    // Noise standard deviation
}

model {
  // Priors
  beta0 ~ normal(0, 10);
  beta1 ~ normal(0, 10);
  sigma ~ cauchy(0, 2);
  
  // Likelihood
  y ~ normal(beta0 + beta1 * x, sigma);
}