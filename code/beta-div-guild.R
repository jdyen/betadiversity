# R code to fit beta regression models of beta diversity against
#  functional groupings of species

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
n.chains <- 3

# bird nesting guilds - point
beta_temp <- do.call('rbind', beta.nest.pt[c(1, 2, 3)])
guild_temp <- rep(names(beta.nest.pt[c(1, 2, 3)]),
                  times = sapply(beta.nest.pt[c(1, 2, 3)], nrow))
region_temp <- rep(bird.temp.pt$region, times = 3)

mod_temp <- stan_guild_mod(y = beta_temp[, 1],
                           guild = guild_temp,
                           region = region_temp,
                           n.its = n.its,
                           n.ch = n.chains,
                           control = list(adapt_delta = 0.98))
save(mod_temp, file = './outputs/mod-bird-nest-turnover-pt.R')

mod_temp <- stan_guild_mod(y = beta_temp[, 2],
                           guild = guild_temp,
                           region = region_temp,
                           n.its = n.its,
                           n.ch = n.chains,
                           control = list(adapt_delta = 0.98))
save(mod_temp, file = './outputs/mod-bird-nest-nestedness-pt.R')

# bird nesting guilds - canyon
beta_temp <- do.call('rbind', beta.nest.can[c(1, 2, 3)])
guild_temp <- rep(names(beta.nest.can[c(1, 2, 3)]),
                  times = sapply(beta.nest.can[c(1, 2, 3)], nrow))
region_temp <- rep(bird.temp.can$region, times = 3)

mod_temp <- stan_guild_mod(y = beta_temp[, 1],
                           guild = guild_temp,
                           region = region_temp,
                           n.its = n.its,
                           n.ch = n.chains,
                           control = list(adapt_delta = 0.98))
save(mod_temp, file = './outputs/mod-bird-nest-turnover-can.R')

mod_temp <- stan_guild_mod(y = beta_temp[, 2],
                           guild = guild_temp,
                           region = region_temp,
                           n.its = n.its,
                           n.ch = n.chains,
                           control = list(adapt_delta = 0.98))
save(mod_temp, file = './outputs/mod-bird-nest-nestedness-can.R')

# bird riparian guilds - point
beta_temp <- do.call('rbind', beta.ripa.pt[1:3])
guild_temp <- rep(names(beta.ripa.pt[1:3]),
                  times = sapply(beta.ripa.pt[1:3], nrow))
region_temp <- rep(bird.temp.pt$region, times = 3)

mod_temp <- stan_guild_mod(y = beta_temp[, 1],
                           guild = guild_temp,
                           region = region_temp,
                           n.its = n.its,
                           n.ch = n.chains,
                           control = list(adapt_delta = 0.95))
save(mod_temp, file = './outputs/mod-bird-ripa-turnover-pt.R')

mod_temp <- stan_guild_mod(y = beta_temp[, 2],
                           guild = guild_temp,
                           region = region_temp,
                           n.its = n.its,
                           n.ch = n.chains,
                           control = list(adapt_delta = 0.95))
save(mod_temp, file = './outputs/mod-bird-ripa-nestedness-pt.R')

# bird riparian guilds - canyon
beta_temp <- do.call('rbind', beta.ripa.can[1:3])
guild_temp <- rep(names(beta.ripa.can[1:3]),
                  times = sapply(beta.ripa.can[1:3], nrow))
region_temp <- rep(bird.temp.can$region, times = 3)

mod_temp <- stan_guild_mod(y = beta_temp[, 1],
                           guild = guild_temp,
                           region = region_temp,
                           n.its = n.its,
                           n.ch = n.chains,
                           control = list(adapt_delta = 0.95))
save(mod_temp, file = './outputs/mod-bird-ripa-turnover-can.R')

mod_temp <- stan_guild_mod(y = beta_temp[, 2],
                           guild = guild_temp,
                           region = region_temp,
                           n.its = n.its,
                           n.ch = n.chains,
                           control = list(adapt_delta = 0.95))
save(mod_temp, file = './outputs/mod-bird-ripa-nestedness-can.R')

# butterfly overwintering guilds - transect
beta_temp <- do.call('rbind', beta.owint.tran)
guild_temp <- rep(names(beta.owint.tran),
                  times = sapply(beta.owint.tran, nrow))
region_temp <- rep(bfs.temp.tran$region, times = 4)

mod_temp <- stan_guild_mod(y = beta_temp[, 1],
                           guild = guild_temp,
                           region = region_temp,
                           n.its = n.its,
                           n.ch = n.chains,
                           control = list(adapt_delta = 0.95))
save(mod_temp, file = './outputs/mod-bfs-owint-turnover-tran.R')

mod_temp <- stan_guild_mod(y = beta_temp[, 2],
                           guild = guild_temp,
                           region = region_temp,
                           n.its = n.its,
                           n.ch = n.chains,
                           control = list(adapt_delta = 0.95))
save(mod_temp, file = './outputs/mod-bfs-owint-nestedness-tran.R')

# butterfly overwintering guilds - canyon
beta_temp <- do.call('rbind', beta.owint.can)
guild_temp <- rep(names(beta.owint.can),
                  times = sapply(beta.owint.can, nrow))
region_temp <- rep(bfs.temp.can$region, times = 4)

mod_temp <- stan_guild_mod(y = beta_temp[, 1],
                           guild = guild_temp,
                           region = region_temp,
                           n.its = n.its,
                           n.ch = n.chains,
                           control = list(adapt_delta = 0.95))
save(mod_temp, file = './outputs/mod-bfs-owint-turnover-can.R')

mod_temp <- stan_guild_mod(y = beta_temp[, 2],
                           guild = guild_temp,
                           region = region_temp,
                           n.its = n.its,
                           n.ch = n.chains,
                           control = list(adapt_delta = 0.95))
save(mod_temp, file = './outputs/mod-bfs-owint-nestedness-can.R')

# butterfly vagility guilds - transect
beta_temp <- do.call('rbind', beta.vgty.tran)
guild_temp <- rep(names(beta.vgty.tran),
                  times = sapply(beta.vgty.tran, nrow))
region_temp <- rep(bfs.temp.tran$region, times = 4)

mod_temp <- stan_guild_mod(y = beta_temp[, 1],
                           guild = guild_temp,
                           region = region_temp,
                           n.its = n.its,
                           n.ch = n.chains,
                           control = list(adapt_delta = 0.95))
save(mod_temp, file = './outputs/mod-bfs-vgty-turnover-tran.R')

mod_temp <- stan_guild_mod(y = beta_temp[, 2],
                           guild = guild_temp,
                           region = region_temp,
                           n.its = n.its,
                           n.ch = n.chains,
                           control = list(adapt_delta = 0.95))
save(mod_temp, file = './outputs/mod-bfs-vgty-nestedness-tran.R')

# butterfly overwintering guilds - canyon
beta_temp <- do.call('rbind', beta.vgty.can)
guild_temp <- rep(names(beta.vgty.can),
                  times = sapply(beta.vgty.can, nrow))
region_temp <- rep(bfs.temp.can$region, times = 4)

mod_temp <- stan_guild_mod(y = beta_temp[, 1],
                           guild = guild_temp,
                           region = region_temp,
                           n.its = n.its,
                           n.ch = n.chains,
                           control = list(adapt_delta = 0.95))
save(mod_temp, file = './outputs/mod-bfs-vgty-turnover-can.R')

mod_temp <- stan_guild_mod(y = beta_temp[, 2],
                           guild = guild_temp,
                           region = region_temp,
                           n.its = n.its,
                           n.ch = n.chains,
                           control = list(adapt_delta = 0.95))
save(mod_temp, file = './outputs/mod-bfs-vgty-nestedness-can.R')

# bird nesting guilds - point (spatial)
beta_temp <- do.call('rbind', beta.spat.nest.pt[c(1, 2, 3)])
guild_temp <- rep(names(beta.nest.pt[c(1, 2, 3)]),
                  times = sapply(beta.spat.nest.pt[c(1, 2, 3)], nrow))
year_temp <- rep(bird.spat.pt$year, times = length(beta.spat.nest.pt[c(1, 2, 3)]))
region_temp <- rep(bird.spat.pt$region, times = 3)

mod_temp <- stan_guild_mod(y = beta_temp[, 1],
                           guild = guild_temp,
                           year = year_temp,
                           region = region_temp,
                           n.its = n.its,
                           n.ch = n.chains,
                           control = list(adapt_delta = 0.95))
save(mod_temp, file = './outputs/mod-bird-spat-nest-turnover-pt.R')

mod_temp <- stan_guild_mod(y = beta_temp[, 2],
                           guild = guild_temp,
                           year = year_temp,
                           region = region_temp,
                           n.its = n.its,
                           n.ch = n.chains,
                           control = list(adapt_delta = 0.95))
save(mod_temp, file = './outputs/mod-bird-spat-nest-nestedness-pt.R')

# bird nesting guilds - canyon (spatial)
beta_temp <- do.call('rbind', beta.spat.nest.can[c(1, 2, 3)])
guild_temp <- rep(names(beta.nest.can[c(1, 2, 3)]),
                  times = sapply(beta.spat.nest.can[c(1, 2, 3)], nrow))
year_temp <- rep(bird.spat.can$year, times = length(beta.spat.nest.can[c(1, 2, 3)]))
region_temp <- rep(bird.spat.can$region, times = 3)

mod_temp <- stan_guild_mod(y = beta_temp[, 1],
                           guild = guild_temp,
                           year = year_temp,
                           region = region_temp,
                           n.its = n.its,
                           n.ch = n.chains,
                           control = list(adapt_delta = 0.95))
save(mod_temp, file = './outputs/mod-bird-spat-nest-turnover-can.R')

mod_temp <- stan_guild_mod(y = beta_temp[, 2],
                           guild = guild_temp,
                           year = year_temp,
                           region = region_temp,
                           n.its = n.its,
                           n.ch = n.chains,
                           control = list(adapt_delta = 0.95))
save(mod_temp, file = './outputs/mod-bird-spat-nest-nestedness-can.R')

# bird riparian guilds - point (spatial)
beta_temp <- do.call('rbind', beta.spat.ripa.pt[1:3])
guild_temp <- rep(names(beta.ripa.pt[1:3]),
                  times = sapply(beta.spat.ripa.pt[1:3], nrow))
year_temp <- rep(bird.spat.pt$year, times = length(beta.spat.ripa.pt[1:3]))
region_temp <- rep(bird.spat.pt$region, times = 3)

mod_temp <- stan_guild_mod(y = beta_temp[, 1],
                           guild = guild_temp,
                           year = year_temp,
                           region = region_temp,
                           n.its = n.its,
                           n.ch = n.chains,
                           control = list(adapt_delta = 0.95))
save(mod_temp, file = './outputs/mod-bird-spat-ripa-turnover-pt.R')

mod_temp <- stan_guild_mod(y = beta_temp[, 2],
                           guild = guild_temp,
                           year = year_temp,
                           region = region_temp,
                           n.its = n.its,
                           n.ch = n.chains,
                           control = list(adapt_delta = 0.95))
save(mod_temp, file = './outputs/mod-bird-spat-ripa-nestedness-pt.R')

# bird riparian guilds - canyon (spatial)
beta_temp <- do.call('rbind', beta.spat.ripa.can[1:3])
guild_temp <- rep(names(beta.ripa.can[1:3]),
                  times = sapply(beta.spat.ripa.can[1:3], nrow))
year_temp <- rep(bird.spat.can$year, times = length(beta.spat.ripa.can[1:3]))
region_temp <- rep(bird.spat.can$region, times = 3)

mod_temp <- stan_guild_mod(y = beta_temp[, 1],
                           guild = guild_temp,
                           year = year_temp,
                           region = region_temp,
                           n.its = n.its,
                           n.ch = n.chains,
                           control = list(adapt_delta = 0.95))
save(mod_temp, file = './outputs/mod-bird-spat-ripa-turnover-can.R')

mod_temp <- stan_guild_mod(y = beta_temp[, 2],
                           guild = guild_temp,
                           year = year_temp,
                           region = region_temp,
                           n.its = n.its,
                           n.ch = n.chains,
                           control = list(adapt_delta = 0.95))
save(mod_temp, file = './outputs/mod-bird-spat-ripa-nestedness-can.R')

# butterfly overwintering guilds - transect (spatial)
beta_temp <- do.call('rbind', beta.spat.owint.tran)
guild_temp <- rep(names(beta.spat.owint.tran),
                  times = sapply(beta.spat.owint.tran, nrow))
year_temp <- rep(bfs.spat.tran$year, times = length(beta.spat.owint.tran))
region_temp <- rep(bfs.spat.tran$region, times = 4)

mod_temp <- stan_guild_mod(y = beta_temp[, 1],
                           guild = guild_temp,
                           year = year_temp,
                           region = region_temp,
                           n.its = n.its,
                           n.ch = n.chains,
                           control = list(adapt_delta = 0.95))
save(mod_temp, file = './outputs/mod-bfs-spat-owint-turnover-tran.R')

mod_temp <- stan_guild_mod(y = beta_temp[, 2],
                           guild = guild_temp,
                           year = year_temp,
                           region = region_temp,
                           n.its = n.its,
                           n.ch = n.chains,
                           control = list(adapt_delta = 0.95))
save(mod_temp, file = './outputs/mod-bfs-spat-owint-nestedness-tran.R')

# butterfly overwintering guilds - canyon (spatial)
beta_temp <- do.call('rbind', beta.spat.owint.can)
guild_temp <- rep(names(beta.spat.owint.can),
                  times = sapply(beta.spat.owint.can, nrow))
year_temp <- rep(bfs.spat.can$year, times = length(beta.spat.owint.can))
region_temp <- rep(bfs.spat.can$region, times = 4)

mod_temp <- stan_guild_mod(y = beta_temp[, 1],
                           guild = guild_temp,
                           year = year_temp,
                           region = region_temp,
                           n.its = n.its,
                           n.ch = n.chains,
                           control = list(adapt_delta = 0.95))
save(mod_temp, file = './outputs/mod-bfs-spat-owint-turnover-can.R')

mod_temp <- stan_guild_mod(y = beta_temp[, 2],
                           guild = guild_temp,
                           year = year_temp,
                           region = region_temp,
                           n.its = n.its,
                           n.ch = n.chains,
                           control = list(adapt_delta = 0.95))
save(mod_temp, file = './outputs/mod-bfs-spat-owint-nestedness-can.R')

# butterfly vagility guilds - transect (spatial)
beta_temp <- do.call('rbind', beta.spat.vgty.tran)
guild_temp <- rep(names(beta.spat.vgty.tran),
                  times = sapply(beta.spat.vgty.tran, nrow))
year_temp <- rep(bfs.spat.tran$year, times = length(beta.spat.vgty.tran))
region_temp <- rep(bfs.spat.tran$region, times = 4)

mod_temp <- stan_guild_mod(y = beta_temp[, 1],
                           guild = guild_temp,
                           year = year_temp,
                           region = region_temp,
                           n.its = n.its,
                           n.ch = n.chains,
                           control = list(adapt_delta = 0.95))
save(mod_temp, file = './outputs/mod-bfs-spat-vgty-turnover-tran.R')

mod_temp <- stan_guild_mod(y = beta_temp[, 2],
                           guild = guild_temp,
                           year = year_temp,
                           region = region_temp,
                           n.its = n.its,
                           n.ch = n.chains,
                           control = list(adapt_delta = 0.95))
save(mod_temp, file = './outputs/mod-bfs-spat-vgty-nestedness-tran.R')

# butterfly overwintering guilds - canyon (spatial)
beta_temp <- do.call('rbind', beta.spat.vgty.can)
guild_temp <- rep(names(beta.spat.vgty.can),
                  times = sapply(beta.spat.vgty.can, nrow))
year_temp <- rep(bfs.spat.can$year, times = length(beta.spat.vgty.can))
region_temp <- rep(bfs.spat.can$region, times = 4)

mod_temp <- stan_guild_mod(y = beta_temp[, 1],
                           guild = guild_temp,
                           year = year_temp,
                           region = region_temp,
                           n.its = n.its,
                           n.ch = n.chains,
                           control = list(adapt_delta = 0.95))
save(mod_temp, file = './outputs/mod-bfs-spat-vgty-turnover-can.R')

mod_temp <- stan_guild_mod(y = beta_temp[, 2],
                           guild = guild_temp,
                           year = year_temp,
                           region = region_temp,
                           n.its = n.its,
                           n.ch = n.chains,
                           control = list(adapt_delta = 0.95))
save(mod_temp, file = './outputs/mod-bfs-spat-vgty-nestedness-can.R')


# clear workspace
rm(list = ls())

# load libraries
library(rstan)
library(loo)
library(parallel)

# set working directory
setwd('~/Dropbox/Erica_GB_indicators/turnover-analysis/')

# load_data
load('./data/full-data-load.RData')

# prepare data
source('./r-code/prepare-lm-data.R')
source('./r-code/r-code-scale-model.R')

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
save(mod_bird_temp1, file = './outputs/mod-bird-temp-turnover.R')

mod_bird_temp2 <- stan_scale_mod(y = bird.temp.beta[, 2],
                                 region = bird.temp.region,
                                 scale = bird.temp.scale,
                                 n.its = n.its,
                                 n.ch = n.chains,
                                 control = list(adapt_delta = 0.95))
save(mod_bird_temp2, file = './outputs/mod-bird-temp-nestedness.R')

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
save(mod_bird_spat1, file = './outputs/mod-bird-spat-turnover.R')

mod_bird_spat2 <- stan_scale_mod(y = bird.spat.beta[, 2],
                                 region = bird.spat.region,
                                 scale = bird.spat.scale,
                                 year = bird.spat.year,
                                 n.its = n.its,
                                 n.ch = n.chains,
                                 control = list(adapt_delta = 0.95))
save(mod_bird_spat2, file = './outputs/mod-bird-spat-nestedness.R')

bfs.temp.beta <- rbind(bfs.temp.tran$beta, bfs.temp.can$beta)
bfs.temp.region <- c(bfs.temp.tran$region, bfs.temp.can$region)
bfs.temp.scale <- c(rep(1, nrow(bfs.temp.tran$beta)), rep(2, nrow(bfs.temp.can$beta)))

mod_bfs_temp1 <- stan_scale_mod(y = bfs.temp.beta[, 1],
                                region = bfs.temp.region,
                                scale = bfs.temp.scale,
                                n.its = n.its,
                                n.ch = n.chains,
                                control = list(adapt_delta = 0.95))
save(mod_bfs_temp1, file = './outputs/mod-bfs-temp-turnover.R')

mod_bfs_temp2 <- stan_scale_mod(y = bfs.temp.beta[, 2],
                                region = bfs.temp.region,
                                scale = bfs.temp.scale,
                                n.its = n.its,
                                n.ch = n.chains,
                                control = list(adapt_delta = 0.95))
save(mod_bfs_temp2, file = './outputs/mod-bfs-temp-nestedness.R')

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
save(mod_bfs_spat1, file = './outputs/mod-bfs-spat-turnover.R')

mod_bfs_spat2 <- stan_scale_mod(y = bfs.spat.beta[, 2],
                                region = bfs.spat.region,
                                scale = bfs.spat.scale,
                                year = bfs.spat.year,
                                n.its = n.its,
                                n.ch = n.chains,
                                control = list(adapt_delta = 0.95))
save(mod_bfs_spat2, file = './outputs/mod-bfs-spat-nestedness.R')

