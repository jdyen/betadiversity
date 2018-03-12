data {
  int<lower=0> n;
  int<lower=0> ng;
  int<lower=0> nreg;
  vector[n] y;
  int<lower=0> guild[n];
  int<lower=0> region[n];
  int<lower=0> npred;
  int<lower=0> guild_pred[npred];
}

parameters {
  real alpha;
  vector[ng] beta_guild;
  vector[nreg] gamma_region;
  real<lower=0> sigma_guild;
  real<lower=0> sigma_region;
  real<lower=0> phi;
}

model {
  for (i in 1:n)
    y[i] ~ beta(inv_logit(alpha + beta_guild[guild[i]] + gamma_region[region[i]]) * phi,
             (1 - inv_logit(alpha + beta_guild[guild[i]] + gamma_region[region[i]])) * phi);

  beta_guild ~ normal(0, sigma_guild);
  gamma_region ~ normal(0, 0.00001);
  sigma_region ~ normal(0, 1);
  sigma_guild ~ normal(0, 1);
  
  phi ~ normal(0, 1);

}

generated quantities {
  vector[n] log_lik;
  vector[n] y_fitted;
  vector[npred] y_pred;

  for (i in 1:n)
    y_fitted[i] = beta_rng(inv_logit(alpha + beta_guild[guild[i]] + gamma_region[region[i]]) * phi,
             (1 - inv_logit(alpha + beta_guild[guild[i]] + gamma_region[region[i]])) * phi);

  for (i in 1:npred)
    y_pred[i] = beta_rng(inv_logit(alpha + beta_guild[guild_pred[i]]) * phi, (1 - inv_logit(alpha + beta_guild[guild_pred[i]])) * phi);

  for (i in 1:n)
    log_lik[i] = beta_lpdf(y[i] | inv_logit(alpha + beta_guild[guild[i]] + gamma_region[region[i]]) * phi, (1 - inv_logit(alpha + beta_guild[guild[i]]) + gamma_region[region[i]]) * phi);
}

