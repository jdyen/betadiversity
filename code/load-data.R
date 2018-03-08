library(betapart)

source('./r-code/beta-div-helper-funs.R')

if (load.from.scratch) {
  # load compositional data and env data
  source('./r-code/load-data-gb.R')
  
  # calculate beta diversity indices
  source('./r-code/calc-beta-div.R')
  
  # calculate beta diversity for guilds
  source('./r-code/calc-beta-div-guilds.R')
  
} else {
  load('./data/full-data-load.RData')
}

