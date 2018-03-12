if (!require(rstan)) {
  install.packages("rstan")
  library(rstan)
}
if (!require(loo)) {
  install.packages("loo")
  library(loo)
}
if (!require(parallel)) {
  install.packages("parallel")
  library(parallel)
}
