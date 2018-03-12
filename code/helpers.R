# helper function
remove_na_fun <- function(x) {
  if (any(is.na(x))) {
    x[which(is.na(x))] <- mean(x, na.rm = TRUE)
  }
  x
}

# full mod function
stan_beta_mod <- function(y, x,
                          region, range, year = NULL,
                          n.its = 5000, n.ch = 4,
                          par.run = TRUE,
                          int_set = TRUE,
                          quad_set = TRUE,
                          region_set = FALSE,
                          ...)
{
  if (is.null(year)) {
    mod.file <- stan_gen_beta_mod(spatial = FALSE, cross_val = FALSE)
  } else {
    mod.file <- stan_gen_beta_mod(spatial = TRUE, cross_val = FALSE)
  }
  print(mod.file)
  if (region_set) {
    x <- cbind(x, region)
    colnames(x)[ncol(x)] <- 'region'
  }
  # prepare data
  x_tmp <- apply(x, 2, remove_na_fun)
  form_set <- paste(' ~ ', prep_formula(x_tmp, int = int_set, quad = quad_set, region_int = region_set), sep = '')
  data_tmp <- model.frame(form_set, data = as.data.frame(x_tmp))
  x_use <- model.matrix(terms(data_tmp), data = data_tmp)[, -1]

  # remove missing observations
  if (any(is.na(y))) {
    rows.to.rm <- which(is.na(y))
    y <- y[-rows.to.rm]
    x_use <- x_use[-rows.to.rm, ]
    region <- as.integer(as.factor(region[-rows.to.rm]))
    range <- as.integer(as.factor(range[-rows.to.rm]))
    if (!is.null(year)) {
      year <- as.integer(as.factor(year[-rows.to.rm]))
    }
  } else {
    region <- as.integer(as.factor(region))
    range <- as.integer(as.factor(range))
    if (!is.null(year)) {
      year <- as.integer(as.factor(year))
    }
  }
  y <- ifelse(y == 1, 0.99, y)
  y <- ifelse(y == 0, 0.01, y)
  if (any(is.na(region))) {
    region[which(is.na(region))] <- 1
  }
  if (any(is.na(range))) {
    range[which(is.na(range))] <- max(range, na.rm = TRUE) + 1
  }
  if (!is.null(year)) {
    if (any(is.na(year))) {
      year[which(is.na(year))] <- max(year, na.rm = TRUE) + 1
    }
  }
  x_use <- sweep(x_use, 2, apply(x_use, 2, mean), '-')
  x_use <- sweep(x_use, 2, apply(x_use, 2, sd), '/')
  
  # collate data set for stan model
  data.set <- list(y = y,
                   x = x_use,
                   region = as.integer(region),
                   range = as.integer(range),
                   n = length(y),
                   nk = ncol(x_use),
                   nrange = length(unique(range)),
                   nregion = length(unique(region)))
  if (!is.null(year)) {
    data.set$year <- as.integer(year)
    data.set$nyear <- length(unique(year))
  }
  
  mod <- stan(file = mod.file,
              data = data.set,
              iter = n.its,
              pars = c('mu', 'log_lik', 'beta_shrink', 'alpha', 'cov_eff', 'y_fitted'),
              chains = n.ch,
              init = 'random',
              seed = 102241451,
              cores = ifelse(par.run, parallel::detectCores(), 1),
              ...)
  
  return(list(mod = mod, data = data.set))
}

# full mod function
stan_beta_mod_cv <- function(y, x,
                             region, range, year = NULL,
                             n.its = 5000, n.ch = 4,
                             n.cv = 10,
                             par.run = TRUE,
                             int_set = TRUE,
                             quad_set = TRUE,
                             region_set = FALSE,
                             ...)
{
  if (region_set) {
    x <- cbind(x, region)
    colnames(x)[ncol(x)] <- 'region'
  }
  # prepare data
  x_tmp <- apply(x, 2, remove_na_fun)
  form_set <- paste(' ~ ', prep_formula(x_tmp, int = int_set, quad = quad_set, region_int = region_set), sep = '')
  data_tmp <- model.frame(form_set, data = as.data.frame(x_tmp))
  x_use <- model.matrix(terms(data_tmp), data = data_tmp)[, -1]
  
  # remove missing observations
  if (any(is.na(y))) {
    rows.to.rm <- which(is.na(y))
    y <- y[-rows.to.rm]
    x_use <- x_use[-rows.to.rm, ]
    region <- as.integer(as.factor(region[-rows.to.rm]))
    range <- as.integer(as.factor(range[-rows.to.rm]))
    if (!is.null(year)) {
      year <- as.integer(as.factor(year[-rows.to.rm]))
    }
  } else {
    region <- as.integer(as.factor(region))
    range <- as.integer(as.factor(range))
    if (!is.null(year)) {
      year <- as.integer(as.factor(year))
    }
  }
  y <- ifelse(y == 1, 0.99, y)
  y <- ifelse(y == 0, 0.01, y)
  if (any(is.na(region))) {
    region[which(is.na(region))] <- 1
  }
  if (any(is.na(range))) {
    range[which(is.na(range))] <- max(range, na.rm = TRUE) + 1
  }
  if (!is.null(year)) {
    if (any(is.na(year))) {
      year[which(is.na(year))] <- max(year, na.rm = TRUE) + 1
    }
  }
  x_use <- sweep(x_use, 2, apply(x_use, 2, mean), '-')
  x_use <- sweep(x_use, 2, apply(x_use, 2, sd), '/')
  
  # collate data set for stan model
  data.set <- list(y = y,
                   x = x_use,
                   region = as.integer(region),
                   range = as.integer(range),
                   n = length(y),
                   nk = ncol(x_use),
                   nregion = length(unique(region)),
                   nrange = length(unique(range)))
  if (!is.null(year)) {
    data.set$year <- year
    data.set$nyear <- length(unique(year))
  }
  
  # set stan model file
  if (is.null(year)) {
    mod.file <- stan_gen_beta_mod(spatial = FALSE, cross_val = TRUE)
  } else {
    mod.file <- stan_gen_beta_mod(spatial = TRUE, cross_val = TRUE)
  }

  mod.def <- stan_model(file = mod.file)
  
  if (par.run & (Sys.info()['sysname'] != 'Windows')) {
    mod <- mclapply(1:n.cv, stan_cv_int, mod.def,
                    data.set, n.cv, n.its, n.ch, par.run = FALSE,
                    mc.cores = parallel::detectCores(),
                    ...)
  } else {
    mod <- lapply(1:n.cv, stan_cv_int, mod.def,
                  data.set, n.cv, n.its, n.ch,
                  par.run = par.run,
                  ...)
  }
  
  out <- do.call('rbind', mod)
  
  r2 <- cor(out[, 1], out[, 2]) ** 2

  return(list(r2_cv = r2, pred = out[, 1], obs = out[, 2]))
}

stan_cv_int <- function(i, mod.def, data, n.cv, n.its, n.ch, par.run, ...) {
  # set up rows to holdout
  n.sample <- floor(data$n / n.cv)
  if (i < n.cv) {
    rows.to.cv <- ((i - 1) * n.sample + 1):(i * n.sample)
  } else {
    rows.to.cv <- ((i - 1) * n.sample + 1):(data$n)
  }
  
  # setup data
  new.y <- data$y[-rows.to.cv]
  new.x <- data$x[-rows.to.cv, ]
  new.reg <- as.integer(as.factor(data$region[-rows.to.cv]))
  new.ran <- as.integer(as.factor(data$range[-rows.to.cv]))
  x.holdout <- data$x[rows.to.cv, ]
  y.holdout <- data$y[rows.to.cv]
  data.cv <- list(y = new.y,
                  x = new.x,
                  region = new.reg,
                  range = new.ran,
                  n = length(new.y),
                  nk = data$nk,
                  nregion = length(unique(new.reg)),
                  nrange = length(unique(new.ran)),
                  n_holdout = length(rows.to.cv),
                  x_holdout = x.holdout)

  if (data.cv$nregion == 1) {
    data.cv$region <- sample(c(1, 2), size = data.cv$n, replace = TRUE)
  }
  data.cv$nregion = length(unique(data.cv$region))
  if (data.cv$nrange == 1) {
    data.cv$range <- sample(c(1, 2), size = data.cv$n, replace = TRUE)
  }
  data.cv$nrange = length(unique(data.cv$range))
  
  if (!is.null(data$year)) {
    new.yr <- as.integer(as.factor(data$year[-rows.to.cv]))
    data.cv$year = new.yr
    data.cv$nyear = length(unique(new.yr))
  }
  
  # fit model
  mod = sampling(object = mod.def,
                 data = data.cv,
                 iter = n.its,
                 pars = c('y_pred'),
                 chains = n.ch,
                 init = 'random',
                 seed = 102241451,
                 cores = ifelse(par.run, parallel::detectCores(), 1),
                 ...)

  pred_vals <- summary(mod, pars = 'y_pred')$summary[, 'mean']
  real_vals <- y.holdout
  if (length(pred_vals) == length(real_vals)) {
    out <- data.frame(pred = pred_vals,
                      obs = real_vals)
  }
  out
}

stan_beta_mod_summary <- function(mod) {
  mean_fitted <- summary(mod$mod, pars = 'y_fitted')$summary[, 'mean']
  y_obs <- mod$data$y
  r2 <- cor(mean_fitted, y_obs) ** 2
  rmsd <- sqrt(mean((mean_fitted - y_obs) ** 2))
  md <- mean((mean_fitted - y_obs))
  mod_ll <- extract_log_lik(mod$mod)
  waic_elpd <- waic(mod_ll)
  loo_elpd <- loo(mod_ll)
  return(list(r2 = r2, rmsd = rmsd, md = md,
              waic_elpd = waic_elpd, loo_elpd = loo_elpd))
}

prep_formula <- function(x, quad = TRUE, int = TRUE, region_int = FALSE) {
  out <- colnames(x)
  if (quad) {
    out2 <- c(out, paste('c(', out, '^2)', sep = ''))
  } else {
    out2 <- out
  }
  if (int) {
    for (i in 1:(length(out) - 1)) {
      for (j in (i + 1):length(out)) {
        out2 <- c(out2, paste(out[i], out[j], sep = ':'))
      }
    }
  }
  if (region_int) {
    out2 <- c(out2, paste(out2, '*region', sep = ''))
  }
  out <- paste(out2, collapse = '+')
  out
}

stan_gen_beta_mod <- function(spatial = FALSE, cross_val = FALSE) {
  cat(
    'data {
    int<lower=0> n;
    int<lower=0> nk;
    int<lower=0> nregion;
    int<lower=0> nrange;
    vector[n] y;
    matrix[n, nk] x;
    int<lower=0> region[n];
    int<lower=0> range[n];',
    if (spatial) {
      '\n  int<lower=0> year[n];
      int<lower=0> nyear;'
    },
    if (cross_val) {
      '\n  int<lower=0> n_holdout;
      matrix[n_holdout, nk] x_holdout;'
    },
    '\n}
    
    parameters {
    real alpha;
    real<lower=0> phi;
    vector[nregion] gamma_region;
    vector[nrange] gamma_range;
    vector[nk] cov_eff;
    real<lower=0> sigma_region;
    real<lower=0> sigma_beta;
    real<lower=0> sigma_range;
    vector<lower=0>[nk] lambda;',
    if (spatial) {
      '\n  vector[nyear] gamma_year;
      real<lower=0> sigma_year;'
    },
    '\n}
    
    transformed parameters {
    vector[n] mu;
    vector[n] a;
    vector[n] b;
    for (i in 1:n)\n',
    if (spatial) {
      '   mu[i] = inv_logit(alpha + x[i, ] * cov_eff + gamma_region[region[i]] + gamma_range[range[i]] + gamma_year[year[i]]);\n'
    } else {
      '   mu[i] = inv_logit(alpha + x[i, ] * cov_eff + gamma_region[region[i]] + gamma_range[range[i]]);\n'
    },
    ' a = mu * phi;
    b = (1 - mu) * phi;
    }
    
    model {
    y ~ beta(a, b);
    cov_eff ~ normal(0, (sigma_beta * lambda));
    lambda ~ cauchy(0, 5);
    sigma_beta ~ normal(0, 2);
    alpha ~ normal(0, 2);
    gamma_region ~ normal(0, sigma_region);
    sigma_region ~ normal(0, 2);
    gamma_range ~ normal(0, sigma_range);
    sigma_range ~ normal(0, 2);',
    if (spatial) {
      '\n  gamma_year ~ normal(0, sigma_year);
      sigma_year ~ normal(0, 2);'
    },
    '\n}
    
    generated quantities {
    vector[n] log_lik;
    vector[nk] beta_shrink;
    vector[n] y_fitted;',
  if (cross_val) {
    '\n  vector[n_holdout] mu_pred;
    vector[n_holdout] y_pred;
    vector[n_holdout] a_holdout;
    vector[n_holdout] b_holdout;'
  },
  '\n  for (i in 1:n)
  y_fitted[i] = beta_rng(a[i], b[i]);',
  if (cross_val) {
    '\n  for (i in 1:n_holdout)
    mu_pred[i] = inv_logit(alpha + x_holdout[i, ] * cov_eff);
    for (i in 1:n_holdout)
    a_holdout[i] = mu_pred[i] * phi;
    for (i in 1:n_holdout)
    b_holdout[i] = (1 - mu_pred[i]) * phi;
    for (i in 1:n_holdout)
    y_pred[i] = beta_rng(a_holdout[i], b_holdout[i]);'
  },
  '\n  for (i in 1:n)
  log_lik[i] = beta_lpdf(y[i] | a[i], b[i]);
  for (k in 1:nk)
  beta_shrink[k] = 1 / (1 + (lambda[k] * lambda[k]));
    }\n\n',
  file = (mod.file <- paste0(tempfile(), ".stan")))
  
  mod.file
}

# full mod function for guild models
stan_guild_mod <- function(y, guild, region, year = NULL,
                           n.its = 5000, n.ch = 4,
                           ...)
{
  # remove missing observations
  if (any(is.na(y))) {
    rows.to.rm <- which(is.na(y))
    y <- y[-rows.to.rm]
    guild <- as.integer(as.factor(guild[-rows.to.rm]))
    region <- as.integer(as.factor(region[-rows.to.rm]))
    if (!is.null(year)) {
      year <- as.integer(as.factor(year[-rows.to.rm]))
    }
  } else {
    guild <- as.integer(as.factor(guild))
    region <- as.integer(as.factor(region))
    if (!is.null(year)) {
      year <- as.integer(as.factor(year))
    }
  }
  y <- ifelse(y == 1, 0.99, y)
  y <- ifelse(y == 0, 0.01, y)
  if (any(is.na(guild))) {
    guild[which(is.na(guild))] <- 1
  }
  if (any(is.na(region))) {
    region[which(is.na(region))] <- 1
  }
  
  # set up predictor matrix
  guild_pred <- unique(guild)
  npred <- length(guild_pred)
  
  # collate data set for stan model
  data.set <- list(y = y,
                   guild = as.integer(guild),
                   region = as.integer(region),
                   n = length(y),
                   nreg = length(unique(region)),
                   ng = length(unique(guild)),
                   npred = npred,
                   guild_pred = guild_pred)
  
  if (!is.null(year)) {
    if (any(is.na(year))) {
      year[which(is.na(year))] <- max(year, na.rm = TRUE) + 1
    }
    data.set$year <- year
    data.set$ny <- length(unique(year))
    mod.file <- './r-code/stan-guild-mod-spat.stan'
  } else {
    mod.file = './r-code/stan-guild-mod.stan'  
  }
  
  print(tapply(data.set$y, data.set$guild, mean))
  print(data.set$guild_pred)
  
  mod.def <- stan_model(file = mod.file)
  
  mod <- sampling(object = mod.def,
                  data = data.set,
                  iter = n.its,
                  chains = n.ch,
                  init = '0',
                  cores = n.ch,
                  ...)
  
  fitted <- summary(mod, pars = 'y_fitted')$summary[, 'mean']
  r2 <- cor(fitted, data.set$y) ** 2
  y_pred <- summary(mod, pars = 'y_pred',
                    p = c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975))$summary[, c('mean', 'sd', '2.5%', '10%', '25%', '50%', '75%', '90%', '97.5%')]
  loglik <- extract_log_lik(mod)
  loo_est <- loo(loglik)
  waic_est <- waic(loglik)
  rhats <- summary(mod)$summary[, 'Rhat']
  neff <- summary(mod)$summary[, 'n_eff']
  return(list(r2 = r2,
              fitted = fitted, obs = data.set$y,
              y_pred = y_pred,
              loo = loo_est, waic = waic_est,
              rhat = rhats, neff = neff,
              mod = mod,
              data_set = data.set))
}

# full mod function for scale models
stan_scale_mod <- function(y, region, scale, year = NULL,
                           n.its = 5000, n.ch = 4,
                           ...)
{
  # remove missing observations
  if (any(is.na(y))) {
    rows.to.rm <- which(is.na(y))
    y <- y[-rows.to.rm]
    region <- as.integer(as.factor(region[-rows.to.rm]))
    scale <- as.integer(as.factor(scale[-rows.to.rm]))
    if (!is.null(year)) {
      year <- as.integer(as.factor(year[-rows.to.rm]))
    }
  } else {
    region <- as.integer(as.factor(region))
    scale <- as.integer(as.factor(scale))
    if (!is.null(year)) {
      year <- as.integer(as.factor(year))
    }
  }
  y <- ifelse(y == 1, 0.99, y)
  y <- ifelse(y == 0, 0.01, y)
  if (any(is.na(region))) {
    region[which(is.na(region))] <- 1
  }
  if (any(is.na(scale))) {
    scale[which(is.na(scale))] <- 1
  }
  
  # set up predictor matrix
  region_pred <- rep(c(1, 2), each = 2)
  scale_pred <- rep(c(1, 2), times = 2)
  npred = length(region_pred)
  
  # collate data set for stan model
  data.set <- list(y = y,
                   region = as.integer(region),
                   scale = as.integer(scale),
                   n = length(y),
                   ns = length(unique(scale)),
                   nr = length(unique(region)),
                   npred = npred,
                   scale_pred = scale_pred,
                   region_pred = region_pred)
  
  if (!is.null(year)) {
    if (any(is.na(year))) {
      year[which(is.na(year))] <- max(year, na.rm = TRUE) + 1
    }
    data.set$year <- year
    data.set$ny <- length(unique(year))
    mod.file <- './r-code/stan-scale-mod-spat.stan'
  } else {
    mod.file = './r-code/stan-scale-mod.stan'  
  }
  
  mod.def <- stan_model(file = mod.file)
  
  mod <- sampling(object = mod.def,
                  data = data.set,
                  iter = n.its,
                  chains = n.ch,
                  init = '0',
                  cores = n.ch,
                  ...) 
  
  fitted <- summary(mod, pars = 'y_fitted')$summary[, 'mean']
  r2 <- cor(fitted, data.set$y) ** 2
  y_pred <- summary(mod, pars = 'y_pred')$summary[, c('mean', 'sd', '2.5%', '50%', '97.5%')]
  loglik <- extract_log_lik(mod)
  loo_est <- loo(loglik)
  waic_est <- waic(loglik)
  rhats <- summary(mod)$summary[, 'Rhat']
  neff <- summary(mod)$summary[, 'n_eff']
  return(list(r2 = r2,
              fitted = fitted, obs = data.set$y,
              y_pred = y_pred,
              loo = loo_est, waic = waic_est,
              rhat = rhats, neff = neff,
              mod = mod,
              data_set = data.set))
}
