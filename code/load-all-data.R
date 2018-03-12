# script to load full data
# Requires:
#  data/cent-birds-survey-data.csv
#  data/west-birds-survey-data.csv
#  data/cent-bfs-survey-data.csv
#  data/west-bfs-survey-data.csv
#  data/spp_exclusions.csv
#  data/cent-birds-cov-pt.csv
#  data/cent-birds-cov-can.csv
#  data/west-birds-cov-pt.csv
#  data/west-birds-cov-can.csv
#  data/cent-bfs-cov-tran.csv
#  data/cent-bfs-cov-can.csv
#  data/west-bfs-cov-tran.csv
#  data/west-bfs-cov-can.csv
#  data/bird-guilds.csv
#  data/bfs-guilds.csv
#  data/cent-pt-coords.csv
#  data/west-pt-coords.csv
#  data/west-bird-can-coords.csv
#  data/cent-bfs-tran-coords.csv
#  data/west-bfs-tran-coords.csv

# load betapart package
if (!require(betapart)) {
  install.packages("betapart")
  library(betapart)
}

source("./code/load-data-helpers.R")

# load compositional data and env data
source("./r-code/load-pa-data.R")

# calculate beta diversity indices
source("./r-code/calculate-beta-diversity.R")
