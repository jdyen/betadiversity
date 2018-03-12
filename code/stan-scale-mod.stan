data {
  int<lower=0> n;
  int<lower=0> ns;
  int<lower=0> nr;
  vector[n] y;
  int<lower=0> scale[n];
  int<lower=0> region[n];
  int<lower=0> npred;
  int<lower=0> scale_pred[npred];
  int<lower=0> region_pred[npred];
}

parameters {
  real alpha;
  vector[ns] beta_scale;
  vector[nr] beta_region;
  real<lower=0> sigma_scale;
  real<lower=0> sigma_region;
  real<lower=0> phi;
}

model {
  for (i in 1:n)
    y[i] ~ beta(inv_logit(alpha + beta_scale[scale[i]] + beta_region[region[i]]) * phi,
             (1 - inv_logit(alpha + beta_scale[scale[i]] + beta_region[region[i]])) * phi);

  beta_scale ~ normal(0, sigma_scale);
  beta_region ~ normal(0, sigma_region);

}

generated quantities {
  vector[n] log_lik;
  vector[n] y_fitted;
  vector[npred] y_pred;

  for (i in 1:n)
    y_fitted[i] = beta_rng(inv_logit(alpha + beta_scale[scale[i]] + beta_region[region[i]]) * phi,
             (1 - inv_logit(alpha + beta_scale[scale[i]] + beta_region[region[i]])) * phi);

  for (i in 1:npred)
    y_pred[i] = beta_rng(inv_logit(alpha + beta_scale[scale_pred[i]] + beta_region[region_pred[i]]) * phi, (1 - inv_logit(alpha + beta_scale[scale_pred[i]] + beta_region[region_pred[i]])) * phi);

  for (i in 1:n)
    log_lik[i] = beta_lpdf(y[i] | inv_logit(alpha + beta_scale[scale[i]] + beta_region[region[i]]) * phi, (1 - inv_logit(alpha + beta_scale[scale[i]] + beta_region[region[i]])) * phi);
}

