# R code to fit beta regression models of beta diversity against
#  spatial extents

# load libraries
source("./code/install_packages.R")

# optional: set working directory
# setwd("PATH/TO/DIR")

# load_data
load("./data/pre_loaded_data.RData")

# prepare data
source("./code/prepare-lm-data.R")
source("./code/helpers.R")

# set MCMC settings
n.its <- 5000
n.chains <- 4

bird.temp.beta <- rbind(bird.temp.pt$beta, bird.temp.can$beta)
bird.temp.region <- c(bird.temp.pt$region, bird.temp.can$region)
bird.temp.scale <- c(rep(1, nrow(bird.temp.pt$beta)), rep(2, nrow(bird.temp.can$beta)))

mod_bird_temp1 <- stan_scale_mod(y = bird.temp.beta[, 1],
                                 region = bird.temp.region,
                                 scale = bird.temp.scale,
                                 n.its = n.its,
                                 n.ch = n.chains,
                                 control = list(adapt_delta = 0.95))
save(mod_bird_temp1, file = "./outputs/mod-bird-temp-turnover.R")

mod_bird_temp2 <- stan_scale_mod(y = bird.temp.beta[, 2],
                                 region = bird.temp.region,
                                 scale = bird.temp.scale,
                                 n.its = n.its,
                                 n.ch = n.chains,
                                 control = list(adapt_delta = 0.95))
save(mod_bird_temp2, file = "./outputs/mod-bird-temp-nestedness.R")

bird.spat.beta <- rbind(bird.spat.pt$beta, bird.spat.can$beta)
bird.spat.region <- c(bird.spat.pt$region, bird.spat.can$region)
bird.spat.scale <- c(rep(1, nrow(bird.spat.pt$beta)), rep(2, nrow(bird.spat.can$beta)))
bird.spat.year <- c(bird.spat.pt$year, bird.spat.can$year)

mod_bird_spat1 <- stan_scale_mod(y = bird.spat.beta[, 1],
                                 region = bird.spat.region,
                                 scale = bird.spat.scale,
                                 year = bird.spat.year,
                                 n.its = n.its,
                                 n.ch = n.chains,
                                 control = list(adapt_delta = 0.95))
save(mod_bird_spat1, file = "./outputs/mod-bird-spat-turnover.R")

mod_bird_spat2 <- stan_scale_mod(y = bird.spat.beta[, 2],
                                 region = bird.spat.region,
                                 scale = bird.spat.scale,
                                 year = bird.spat.year,
                                 n.its = n.its,
                                 n.ch = n.chains,
                                 control = list(adapt_delta = 0.95))
save(mod_bird_spat2, file = "./outputs/mod-bird-spat-nestedness.R")

bfs.temp.beta <- rbind(bfs.temp.tran$beta, bfs.temp.can$beta)
bfs.temp.region <- c(bfs.temp.tran$region, bfs.temp.can$region)
bfs.temp.scale <- c(rep(1, nrow(bfs.temp.tran$beta)), rep(2, nrow(bfs.temp.can$beta)))

mod_bfs_temp1 <- stan_scale_mod(y = bfs.temp.beta[, 1],
                                 region = bfs.temp.region,
                                 scale = bfs.temp.scale,
                                 n.its = n.its,
                                 n.ch = n.chains,
                                 control = list(adapt_delta = 0.95))
save(mod_bfs_temp1, file = "./outputs/mod-bfs-temp-turnover.R")

mod_bfs_temp2 <- stan_scale_mod(y = bfs.temp.beta[, 2],
                                region = bfs.temp.region,
                                scale = bfs.temp.scale,
                                n.its = n.its,
                                n.ch = n.chains,
                                control = list(adapt_delta = 0.95))
save(mod_bfs_temp2, file = "./outputs/mod-bfs-temp-nestedness.R")

bfs.spat.beta <- rbind(bfs.spat.tran$beta, bfs.spat.can$beta)
bfs.spat.region <- c(bfs.spat.tran$region, bfs.spat.can$region)
bfs.spat.scale <- c(rep(1, nrow(bfs.spat.tran$beta)), rep(2, nrow(bfs.spat.can$beta)))
bfs.spat.year <- c(bfs.spat.tran$year, bfs.spat.can$year)

mod_bfs_spat1 <- stan_scale_mod(y = bfs.spat.beta[, 1],
                                region = bfs.spat.region,
                                scale = bfs.spat.scale,
                                year = bfs.spat.year,
                                n.its = n.its,
                                n.ch = n.chains,
                                control = list(adapt_delta = 0.95))
save(mod_bfs_spat1, file = "./outputs/mod-bfs-spat-turnover.R")

mod_bfs_spat2 <- stan_scale_mod(y = bfs.spat.beta[, 2],
                                region = bfs.spat.region,
                                scale = bfs.spat.scale,
                                year = bfs.spat.year,
                                n.its = n.its,
                                n.ch = n.chains,
                                control = list(adapt_delta = 0.95))
save(mod_bfs_spat2, file = "./outputs/mod-bfs-spat-nestedness.R")
