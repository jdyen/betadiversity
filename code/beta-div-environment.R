# R code to fit beta regression models of beta diversity against
#  environmental predictors

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
n.its.cv <- 4000
n.chains <- 4
n.ch.cv <- 4

# define model to be fitted
quad_mod <- quad_mod_cv <- TRUE       # include quadratic terms?
int_mod <- int_mod_cv <- FALSE        # include pairwise interactions?
region_mod <- region_mod_cv <- TRUE   # allow differences in effects between regions?

file_suffix <- paste0(ifelse(quad_mod, "quad_", "lin_"),
                      ifelse (region_mod, "reg_", "one_"),
                      ifelse(int_mod, "int", "idp"))

mod_bird <- stan_beta_mod(y = bird.temp.pt$beta[, 1],
                          x = bird.temp.pt$pred,
                          region = bird.temp.pt$region,
                          range = bird.temp.pt$range,
                          year = NULL,
                          n.its = n.its,
                          n.ch = n.chains,
                          par.run = TRUE,
                          quad_set = quad_mod,
                          int_set = int_mod,
                          region_set = region_mod,
                          control = list(adapt_delta = 0.98, max_treedepth = 20))
save(mod_bird, file = paste0("./outputs/mod-bird-pt1_", file_suffix, ".R"))
rm("mod_bird")
mod_bird_cv <- stan_beta_mod_cv(y = bird.temp.pt$beta[, 1],
                                x = bird.temp.pt$pred,
                                region = bird.temp.pt$region,
                                range = bird.temp.pt$range,
                                year = NULL,
                                n.its = n.its.cv,
                                n.ch = n.ch.cv,
                                n.cv = 10,
                                par.run = TRUE,
                                quad_set = quad_mod_cv,
                                int_set = int_mod_cv,
                                region_set = region_mod_cv,
                                control = list(adapt_delta = 0.95))
save(mod_bird_cv, file = paste0("./outputs/mod-bird-pt1_", file_suffix, "_cv.R"))
rm("mod_bird_cv")
mod_bird <- stan_beta_mod(y = bird.temp.pt$beta[, 2],
                          x = bird.temp.pt$pred,
                          region = bird.temp.pt$region,
                          range = bird.temp.pt$range,
                          year = NULL,
                          n.its = n.its,
                          n.ch = n.chains,
                          par.run = TRUE,
                          quad_set = quad_mod,
                          int_set = int_mod,
                          region_set = region_mod,
                          control = list(adapt_delta = 0.98, max_treedepth = 20))
save(mod_bird, file = paste0("./outputs/mod-bird-pt2_", file_suffix, ".R"))
rm("mod_bird")
mod_bird_cv <- stan_beta_mod_cv(y = bird.temp.pt$beta[, 2],
                                x = bird.temp.pt$pred,
                                region = bird.temp.pt$region,
                                range = bird.temp.pt$range,
                                year = NULL,
                                n.its = n.its.cv,
                                n.ch = n.ch.cv,
                                n.cv = 10,
                                par.run = TRUE,
                                quad_set = quad_mod_cv,
                                int_set = int_mod_cv,
                                region_set = region_mod_cv)
save(mod_bird_cv, file = paste0("./outputs/mod-bird-pt2_", file_suffix, "_cv.R"))
rm("mod_bird_cv")
mod_bird <- stan_beta_mod(y = bird.temp.can$beta[, 1],
                          x = bird.temp.can$pred,
                          region = bird.temp.can$region,
                          range = bird.temp.can$range,
                          year = NULL,
                          n.its = n.its,
                          n.ch = n.chains,
                          par.run = TRUE,
                          quad_set = quad_mod,
                          int_set = int_mod,
                          region_set = region_mod,
                          control = list(adapt_delta = 0.98, max_treedepth = 20))
save(mod_bird, file = paste0("./outputs/mod-bird-can1_", file_suffix, ".R"))
rm("mod_bird")
mod_bird_cv <- stan_beta_mod_cv(y = bird.temp.can$beta[, 1],
                                x = bird.temp.can$pred,
                                region = bird.temp.can$region,
                                range = bird.temp.can$range,
                                year = NULL,
                                n.its = n.its.cv,
                                n.ch = n.ch.cv,
                                n.cv = 10,
                                par.run = TRUE,
                                quad_set = quad_mod_cv,
                                int_set = int_mod_cv,
                                region_set = region_mod_cv)
save(mod_bird_cv, file = paste0("./outputs/mod-bird-can1_", file_suffix, "_cv.R"))
rm("mod_bird_cv")
mod_bird <- stan_beta_mod(y = bird.temp.can$beta[, 2],
                          x = bird.temp.can$pred,
                          region = bird.temp.can$region,
                          range = bird.temp.can$range,
                          year = NULL,
                          n.its = n.its,
                          n.ch = n.chains,
                          par.run = TRUE,
                          quad_set = quad_mod,
                          int_set = int_mod,
                          region_set = region_mod,
                          control = list(adapt_delta = 0.98, max_treedepth = 20))
save(mod_bird, file = paste0("./outputs/mod-bird-can2_", file_suffix, ".R"))
rm("mod_bird")
mod_bird_cv <- stan_beta_mod_cv(y = bird.temp.can$beta[, 2],
                                x = bird.temp.can$pred,
                                region = bird.temp.can$region,
                                range = bird.temp.can$range,
                                year = NULL,
                                n.its = n.its.cv,
                                n.ch = n.ch.cv,
                                n.cv = 10,
                                par.run = TRUE,
                                quad_set = quad_mod_cv,
                                int_set = int_mod_cv,
                                region_set = region_mod_cv)
save(mod_bird_cv, file = paste0("./outputs/mod-bird-can2_", file_suffix, "_cv.R"))
rm("mod_bird_cv")
mod_bird <- stan_beta_mod(y = bird.spat.pt$beta[, 1],
                          x = bird.spat.pt$pred,
                          region = bird.spat.pt$region,
                          range = bird.spat.pt$range,
                          year = bird.spat.pt$year,
                          n.its = n.its,
                          n.ch = n.chains,
                          par.run = TRUE,
                          quad_set = quad_mod,
                          int_set = int_mod,
                          region_set = region_mod,
                          control = list(adapt_delta = 0.98, max_treedepth = 20))
save(mod_bird, file = paste0("./outputs/mod-bird-spat1_", file_suffix, ".R"))
rm("mod_bird")
mod_bird_cv <- stan_beta_mod_cv(y = bird.spat.pt$beta[, 1],
                                x = bird.spat.pt$pred,
                                region = bird.spat.pt$region,
                                range = bird.spat.pt$range,
                                year = bird.spat.pt$year,
                                n.its = n.its.cv,
                                n.ch = n.ch.cv,
                                n.cv = 10,
                                par.run = TRUE,
                                quad_set = quad_mod_cv,
                                int_set = int_mod_cv,
                                region_set = region_mod_cv)
save(mod_bird_cv, file = paste0("./outputs/mod-bird-spat1_", file_suffix, "_cv.R"))
rm("mod_bird_cv")
mod_bird <- stan_beta_mod(y = bird.spat.pt$beta[, 2],
                          x = bird.spat.pt$pred,
                          region = bird.spat.pt$region,
                          range = bird.spat.pt$range,
                          year = bird.spat.pt$year,
                          n.its = n.its,
                          n.ch = n.chains,
                          par.run = TRUE,
                          quad_set = quad_mod,
                          int_set = int_mod,
                          region_set = region_mod,
                          control = list(adapt_delta = 0.98, max_treedepth = 20))
save(mod_bird, file = paste0("./outputs/mod-bird-spat2_", file_suffix, ".R"))
rm("mod_bird")
mod_bird_cv <- stan_beta_mod_cv(y = bird.spat.pt$beta[, 2],
                                x = bird.spat.pt$pred,
                                region = bird.spat.pt$region,
                                range = bird.spat.pt$range,
                                year = bird.spat.pt$year,
                                n.its = n.its.cv,
                                n.ch = n.ch.cv,
                                n.cv = 10,
                                par.run = TRUE,
                                quad_set = quad_mod_cv,
                                int_set = int_mod_cv,
                                region_set = region_mod_cv)
save(mod_bird_cv, file = paste0("./outputs/mod-bird-spat2_", file_suffix, "_cv.R"))
rm("mod_bird_cv")

mod_bfs <- stan_beta_mod(y = bfs.temp.tran$beta[, 1],
                         x = bfs.temp.tran$pred,
                         region = bfs.temp.tran$region,
                         range = bfs.temp.tran$range,
                         year = NULL,
                         n.its = n.its,
                         n.ch = n.chains,
                         par.run = TRUE,
                         quad_set = quad_mod,
                         int_set = int_mod,
                         region_set = region_mod,
                         control = list(adapt_delta = 0.98, max_treedepth = 20))
save(mod_bfs, file = paste0("./outputs/mod-bfs-tran1_", file_suffix, ".R"))
rm("mod_bfs")
mod_bfs_cv <- stan_beta_mod_cv(y = bfs.temp.tran$beta[, 1],
                               x = bfs.temp.tran$pred,
                               region = bfs.temp.tran$region,
                               range = bfs.temp.tran$range,
                               year = NULL,
                               n.its = n.its.cv,
                               n.ch = n.ch.cv,
                               n.cv = 10,
                               par.run = TRUE,
                               quad_set = quad_mod_cv,
                               int_set = int_mod_cv,
                               region_set = region_mod_cv)
save(mod_bfs_cv, file = paste0("./outputs/mod-bfs-tran1_", file_suffix, "_cv.R"))
rm("mod_bfs_cv")
mod_bfs <- stan_beta_mod(y = bfs.temp.tran$beta[, 2],
                         x = bfs.temp.tran$pred,
                         region = bfs.temp.tran$region,
                         range = bfs.temp.tran$range,
                         year = NULL,
                         n.its = n.its,
                         n.ch = n.chains,
                         par.run = TRUE,
                         quad_set = quad_mod,
                         int_set = int_mod,
                         region_set = region_mod,
                         control = list(adapt_delta = 0.98, max_treedepth = 20))
save(mod_bfs, file = paste0("./outputs/mod-bfs-tran2_", file_suffix, ".R"))
rm("mod_bfs")
mod_bfs_cv <- stan_beta_mod_cv(y = bfs.temp.tran$beta[, 2],
                               x = bfs.temp.tran$pred,
                               region = bfs.temp.tran$region,
                               range = bfs.temp.tran$range,
                               year = NULL,
                               n.its = n.its.cv,
                               n.ch = n.ch.cv,
                               n.cv = 10,
                               par.run = TRUE,
                               quad_set = quad_mod_cv,
                               int_set = int_mod_cv,
                               region_set = region_mod_cv)
save(mod_bfs_cv, file = paste0("./outputs/mod-bfs-tran2_", file_suffix, "_cv.R"))
rm("mod_bfs_cv")
mod_bfs <- stan_beta_mod(y = bfs.temp.can$beta[, 1],
                         x = bfs.temp.can$pred,
                         region = bfs.temp.can$region,
                         range = bfs.temp.can$range,
                         year = NULL,
                         n.its = n.its,
                         n.ch = n.chains,
                         par.run = TRUE,
                         quad_set = quad_mod,
                         int_set = int_mod,
                         region_set = region_mod,
                         control = list(adapt_delta = 0.98, max_treedepth = 20))
save(mod_bfs, file = paste0("./outputs/mod-bfs-can1_", file_suffix, ".R"))
rm("mod_bfs")
mod_bfs_cv <- stan_beta_mod_cv(y = bfs.temp.can$beta[, 1],
                               x = bfs.temp.can$pred,
                               region = bfs.temp.can$region,
                               range = bfs.temp.can$range,
                               year = NULL,
                               n.its = n.its.cv,
                               n.ch = n.ch.cv,
                               n.cv = 10,
                               par.run = TRUE,
                               quad_set = quad_mod_cv,
                               int_set = int_mod_cv,
                               region_set = region_mod_cv)
save(mod_bfs_cv, file = paste0("./outputs/mod-bfs-can1_", file_suffix, "_cv.R"))
rm("mod_bfs_cv")
mod_bfs <- stan_beta_mod(y = bfs.temp.can$beta[, 2],
                         x = bfs.temp.can$pred,
                         region = bfs.temp.can$region,
                         range = bfs.temp.can$range,
                         year = NULL,
                         n.its = n.its,
                         n.ch = n.chains,
                         par.run = TRUE,
                         quad_set = quad_mod,
                         int_set = int_mod,
                         region_set = region_mod,
                         control = list(adapt_delta = 0.98, max_treedepth = 20))
save(mod_bfs, file = paste0("./outputs/mod-bfs-can2_", file_suffix, ".R"))
rm("mod_bfs")
mod_bfs_cv <- stan_beta_mod_cv(y = bfs.temp.can$beta[, 2],
                               x = bfs.temp.can$pred,
                               region = bfs.temp.can$region,
                               range = bfs.temp.can$range,
                               year = NULL,
                               n.its = n.its.cv,
                               n.ch = n.ch.cv,
                               n.cv = 10,
                               par.run = TRUE,
                               quad_set = quad_mod_cv,
                               int_set = int_mod_cv,
                               region_set = region_mod_cv)
save(mod_bfs_cv, file = paste0("./outputs/mod-bfs-can2_", file_suffix, "_cv.R"))
rm("mod_bfs_cv")
mod_bfs <- stan_beta_mod(y = bfs.spat.tran$beta[, 1],
                         x = bfs.spat.tran$pred,
                         region = bfs.spat.tran$region,
                         range = bfs.spat.tran$range,
                         year = bfs.spat.tran$year,
                         n.its = n.its,
                         n.ch = n.chains,
                         par.run = TRUE,
                         quad_set = quad_mod,
                         int_set = int_mod,
                         region_set = region_mod,
                         control = list(adapt_delta = 0.98, max_treedepth = 20))
save(mod_bfs, file = paste0("./outputs/mod-bfs-spat1_", file_suffix, ".R"))
rm("mod_bfs")
mod_bfs_cv <- stan_beta_mod_cv(y = bfs.spat.tran$beta[, 1],
                               x = bfs.spat.tran$pred,
                               region = bfs.spat.tran$region,
                               range = bfs.spat.tran$range,
                               year = bfs.spat.tran$year,
                               n.its = n.its.cv,
                               n.ch = n.ch.cv,
                               n.cv = 10,
                               par.run = TRUE,
                               quad_set = quad_mod_cv,
                               int_set = int_mod_cv,
                               region_set = region_mod_cv)
save(mod_bfs_cv, file = "./outputs/mod-bfs-spat1_lin_reg_int-cv.R")
save(mod_bfs_cv, file = paste0("./outputs/mod-bfs-spat1_", file_suffix, "_cv.R"))
rm("mod_bfs_cv")
mod_bfs <- stan_beta_mod(y = bfs.spat.tran$beta[, 2],
                         x = bfs.spat.tran$pred,
                         region = bfs.spat.tran$region,
                         range = bfs.spat.tran$range,
                         year = bfs.spat.tran$year,
                         n.its = n.its,
                         n.ch = n.chains,
                         par.run = TRUE,
                         quad_set = quad_mod,
                         int_set = int_mod,
                         region_set = region_mod,
                         control = list(adapt_delta = 0.98, max_treedepth = 20))
save(mod_bfs, file = paste0("./outputs/mod-bfs-spat2_", file_suffix, ".R"))
rm("mod_bfs")
mod_bfs_cv <- stan_beta_mod_cv(y = bfs.spat.tran$beta[, 2],
                               x = bfs.spat.tran$pred,
                               region = bfs.spat.tran$region,
                               range = bfs.spat.tran$range,
                               year = bfs.spat.tran$year,
                               n.its = n.its.cv,
                               n.ch = n.ch.cv,
                               n.cv = 10,
                               par.run = TRUE,
                               quad_set = quad_mod_cv,
                               int_set = int_mod_cv,
                               region_set = region_mod_cv)
save(mod_bfs_cv, file = paste0("./outputs/mod-bfs-spat2_", file_suffix, "_cv.R"))
rm("mod_bfs_cv")
