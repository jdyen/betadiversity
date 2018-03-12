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
beta_temp <- do.call("rbind", beta.nest.pt[c(1, 2, 3)])
guild_temp <- rep(names(beta.nest.pt[c(1, 2, 3)]),
                  times = sapply(beta.nest.pt[c(1, 2, 3)], nrow))
region_temp <- rep(bird.temp.pt$region, times = 3)

mod_temp <- stan_guild_mod(y = beta_temp[, 1],
                           guild = guild_temp,
                           region = region_temp,
                           n.its = n.its,
                           n.ch = n.chains,
                           control = list(adapt_delta = 0.98))
save(mod_temp, file = "./outputs/mod-bird-nest-turnover-pt.R")

mod_temp <- stan_guild_mod(y = beta_temp[, 2],
                           guild = guild_temp,
                           region = region_temp,
                           n.its = n.its,
                           n.ch = n.chains,
                           control = list(adapt_delta = 0.98))
save(mod_temp, file = "./outputs/mod-bird-nest-nestedness-pt.R")

# bird nesting guilds - canyon
beta_temp <- do.call("rbind", beta.nest.can[c(1, 2, 3)])
guild_temp <- rep(names(beta.nest.can[c(1, 2, 3)]),
                  times = sapply(beta.nest.can[c(1, 2, 3)], nrow))
region_temp <- rep(bird.temp.can$region, times = 3)

mod_temp <- stan_guild_mod(y = beta_temp[, 1],
                           guild = guild_temp,
                           region = region_temp,
                           n.its = n.its,
                           n.ch = n.chains,
                           control = list(adapt_delta = 0.98))
save(mod_temp, file = "./outputs/mod-bird-nest-turnover-can.R")

mod_temp <- stan_guild_mod(y = beta_temp[, 2],
                           guild = guild_temp,
                           region = region_temp,
                           n.its = n.its,
                           n.ch = n.chains,
                           control = list(adapt_delta = 0.98))
save(mod_temp, file = "./outputs/mod-bird-nest-nestedness-can.R")

# bird riparian guilds - point
beta_temp <- do.call("rbind", beta.ripa.pt[1:3])
guild_temp <- rep(names(beta.ripa.pt[1:3]),
                  times = sapply(beta.ripa.pt[1:3], nrow))
region_temp <- rep(bird.temp.pt$region, times = 3)

mod_temp <- stan_guild_mod(y = beta_temp[, 1],
                           guild = guild_temp,
                           region = region_temp,
                           n.its = n.its,
                           n.ch = n.chains,
                           control = list(adapt_delta = 0.95))
save(mod_temp, file = "./outputs/mod-bird-ripa-turnover-pt.R")

mod_temp <- stan_guild_mod(y = beta_temp[, 2],
                           guild = guild_temp,
                           region = region_temp,
                           n.its = n.its,
                           n.ch = n.chains,
                           control = list(adapt_delta = 0.95))
save(mod_temp, file = "./outputs/mod-bird-ripa-nestedness-pt.R")

# bird riparian guilds - canyon
beta_temp <- do.call("rbind", beta.ripa.can[1:3])
guild_temp <- rep(names(beta.ripa.can[1:3]),
                  times = sapply(beta.ripa.can[1:3], nrow))
region_temp <- rep(bird.temp.can$region, times = 3)

mod_temp <- stan_guild_mod(y = beta_temp[, 1],
                           guild = guild_temp,
                           region = region_temp,
                           n.its = n.its,
                           n.ch = n.chains,
                           control = list(adapt_delta = 0.95))
save(mod_temp, file = "./outputs/mod-bird-ripa-turnover-can.R")

mod_temp <- stan_guild_mod(y = beta_temp[, 2],
                           guild = guild_temp,
                           region = region_temp,
                           n.its = n.its,
                           n.ch = n.chains,
                           control = list(adapt_delta = 0.95))
save(mod_temp, file = "./outputs/mod-bird-ripa-nestedness-can.R")

# butterfly overwintering guilds - transect
beta_temp <- do.call("rbind", beta.owint.tran)
guild_temp <- rep(names(beta.owint.tran),
                  times = sapply(beta.owint.tran, nrow))
region_temp <- rep(bfs.temp.tran$region, times = 4)

mod_temp <- stan_guild_mod(y = beta_temp[, 1],
                           guild = guild_temp,
                           region = region_temp,
                           n.its = n.its,
                           n.ch = n.chains,
                           control = list(adapt_delta = 0.95))
save(mod_temp, file = "./outputs/mod-bfs-owint-turnover-tran.R")

mod_temp <- stan_guild_mod(y = beta_temp[, 2],
                           guild = guild_temp,
                           region = region_temp,
                           n.its = n.its,
                           n.ch = n.chains,
                           control = list(adapt_delta = 0.95))
save(mod_temp, file = "./outputs/mod-bfs-owint-nestedness-tran.R")

# butterfly overwintering guilds - canyon
beta_temp <- do.call("rbind", beta.owint.can)
guild_temp <- rep(names(beta.owint.can),
                  times = sapply(beta.owint.can, nrow))
region_temp <- rep(bfs.temp.can$region, times = 4)

mod_temp <- stan_guild_mod(y = beta_temp[, 1],
                           guild = guild_temp,
                           region = region_temp,
                           n.its = n.its,
                           n.ch = n.chains,
                           control = list(adapt_delta = 0.95))
save(mod_temp, file = "./outputs/mod-bfs-owint-turnover-can.R")

mod_temp <- stan_guild_mod(y = beta_temp[, 2],
                           guild = guild_temp,
                           region = region_temp,
                           n.its = n.its,
                           n.ch = n.chains,
                           control = list(adapt_delta = 0.95))
save(mod_temp, file = "./outputs/mod-bfs-owint-nestedness-can.R")

# butterfly vagility guilds - transect
beta_temp <- do.call("rbind", beta.vgty.tran)
guild_temp <- rep(names(beta.vgty.tran),
                  times = sapply(beta.vgty.tran, nrow))
region_temp <- rep(bfs.temp.tran$region, times = 4)

mod_temp <- stan_guild_mod(y = beta_temp[, 1],
                           guild = guild_temp,
                           region = region_temp,
                           n.its = n.its,
                           n.ch = n.chains,
                           control = list(adapt_delta = 0.95))
save(mod_temp, file = "./outputs/mod-bfs-vgty-turnover-tran.R")

mod_temp <- stan_guild_mod(y = beta_temp[, 2],
                           guild = guild_temp,
                           region = region_temp,
                           n.its = n.its,
                           n.ch = n.chains,
                           control = list(adapt_delta = 0.95))
save(mod_temp, file = "./outputs/mod-bfs-vgty-nestedness-tran.R")

# butterfly overwintering guilds - canyon
beta_temp <- do.call("rbind", beta.vgty.can)
guild_temp <- rep(names(beta.vgty.can),
                  times = sapply(beta.vgty.can, nrow))
region_temp <- rep(bfs.temp.can$region, times = 4)

mod_temp <- stan_guild_mod(y = beta_temp[, 1],
                           guild = guild_temp,
                           region = region_temp,
                           n.its = n.its,
                           n.ch = n.chains,
                           control = list(adapt_delta = 0.95))
save(mod_temp, file = "./outputs/mod-bfs-vgty-turnover-can.R")

mod_temp <- stan_guild_mod(y = beta_temp[, 2],
                           guild = guild_temp,
                           region = region_temp,
                           n.its = n.its,
                           n.ch = n.chains,
                           control = list(adapt_delta = 0.95))
save(mod_temp, file = "./outputs/mod-bfs-vgty-nestedness-can.R")

# bird nesting guilds - point (spatial)
beta_temp <- do.call("rbind", beta.spat.nest.pt[c(1, 2, 3)])
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
save(mod_temp, file = "./outputs/mod-bird-spat-nest-turnover-pt.R")

mod_temp <- stan_guild_mod(y = beta_temp[, 2],
                           guild = guild_temp,
                           year = year_temp,
                           region = region_temp,
                           n.its = n.its,
                           n.ch = n.chains,
                           control = list(adapt_delta = 0.95))
save(mod_temp, file = "./outputs/mod-bird-spat-nest-nestedness-pt.R")

# bird nesting guilds - canyon (spatial)
beta_temp <- do.call("rbind", beta.spat.nest.can[c(1, 2, 3)])
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
save(mod_temp, file = "./outputs/mod-bird-spat-nest-turnover-can.R")

mod_temp <- stan_guild_mod(y = beta_temp[, 2],
                           guild = guild_temp,
                           year = year_temp,
                           region = region_temp,
                           n.its = n.its,
                           n.ch = n.chains,
                           control = list(adapt_delta = 0.95))
save(mod_temp, file = "./outputs/mod-bird-spat-nest-nestedness-can.R")

# bird riparian guilds - point (spatial)
beta_temp <- do.call("rbind", beta.spat.ripa.pt[1:3])
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
save(mod_temp, file = "./outputs/mod-bird-spat-ripa-turnover-pt.R")

mod_temp <- stan_guild_mod(y = beta_temp[, 2],
                           guild = guild_temp,
                           year = year_temp,
                           region = region_temp,
                           n.its = n.its,
                           n.ch = n.chains,
                           control = list(adapt_delta = 0.95))
save(mod_temp, file = "./outputs/mod-bird-spat-ripa-nestedness-pt.R")

# bird riparian guilds - canyon (spatial)
beta_temp <- do.call("rbind", beta.spat.ripa.can[1:3])
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
save(mod_temp, file = "./outputs/mod-bird-spat-ripa-turnover-can.R")

mod_temp <- stan_guild_mod(y = beta_temp[, 2],
                           guild = guild_temp,
                           year = year_temp,
                           region = region_temp,
                           n.its = n.its,
                           n.ch = n.chains,
                           control = list(adapt_delta = 0.95))
save(mod_temp, file = "./outputs/mod-bird-spat-ripa-nestedness-can.R")

# butterfly overwintering guilds - transect (spatial)
beta_temp <- do.call("rbind", beta.spat.owint.tran)
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
save(mod_temp, file = "./outputs/mod-bfs-spat-owint-turnover-tran.R")

mod_temp <- stan_guild_mod(y = beta_temp[, 2],
                           guild = guild_temp,
                           year = year_temp,
                           region = region_temp,
                           n.its = n.its,
                           n.ch = n.chains,
                           control = list(adapt_delta = 0.95))
save(mod_temp, file = "./outputs/mod-bfs-spat-owint-nestedness-tran.R")

# butterfly overwintering guilds - canyon (spatial)
beta_temp <- do.call("rbind", beta.spat.owint.can)
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
save(mod_temp, file = "./outputs/mod-bfs-spat-owint-turnover-can.R")

mod_temp <- stan_guild_mod(y = beta_temp[, 2],
                           guild = guild_temp,
                           year = year_temp,
                           region = region_temp,
                           n.its = n.its,
                           n.ch = n.chains,
                           control = list(adapt_delta = 0.95))
save(mod_temp, file = "./outputs/mod-bfs-spat-owint-nestedness-can.R")

# butterfly vagility guilds - transect (spatial)
beta_temp <- do.call("rbind", beta.spat.vgty.tran)
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
save(mod_temp, file = "./outputs/mod-bfs-spat-vgty-turnover-tran.R")

mod_temp <- stan_guild_mod(y = beta_temp[, 2],
                           guild = guild_temp,
                           year = year_temp,
                           region = region_temp,
                           n.its = n.its,
                           n.ch = n.chains,
                           control = list(adapt_delta = 0.95))
save(mod_temp, file = "./outputs/mod-bfs-spat-vgty-nestedness-tran.R")

# butterfly overwintering guilds - canyon (spatial)
beta_temp <- do.call("rbind", beta.spat.vgty.can)
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
save(mod_temp, file = "./outputs/mod-bfs-spat-vgty-turnover-can.R")

mod_temp <- stan_guild_mod(y = beta_temp[, 2],
                           guild = guild_temp,
                           year = year_temp,
                           region = region_temp,
                           n.its = n.its,
                           n.ch = n.chains,
                           control = list(adapt_delta = 0.95))
save(mod_temp, file = "./outputs/mod-bfs-spat-vgty-nestedness-can.R")
