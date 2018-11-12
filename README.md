# betadiversity 
Code to regress estimates of beta diversity against three sets of predictor variables: spatial scale, functional groupings of species, and environmental covariates. Regressions are based on beta regression models, implemented in Stan through the rstan R package.

The analyses outlined here are presented in detail in Yen, Fleishman, Fogarty, and Dobkin (in press) Relating beta diversity of birds and butterflies in the Great Basin to spatial resolution, environmental variables and trait-based groups. Global Ecology and Biogeography (accepted 1 October 2018).

Copyright &copy; 2018, Jian Yen

*****

## License details
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

*****

## Overview
Beta regression models to relate estimates of beta diversity to three sets of predictor variables.

Maintainer: Jian Yen (jdl.yen@gmail.com)

Updated: 12 November 2018

## Usage
Several scripts are provided. The main model runs are provided in beta-div-[predictor].R files, where [predictor] is the predictor variable of interest. These scripts require a pre-loaded data set (available in the data directory). Several helper scripts are provided, including helpers.R, install_packages.R, and prepare-lm-data.R. The load-all-data.R script loads all data from scratch, calling on many individual data files and two helper scripts (load-data-helpers.R and load-pa-data.R). Raw data files are available in the data directory. A pre-loaded data set is provided because loading data from raw files takes approximately half an hour on a MacBook Air.

The [filename].stan files are Stan models with spatial scale or functional groupings (guilds) as predictor variables. Models with environmental predictor variables are created dynamically using the stan_gen_beta_mod() function in helpers.R.

Data files are included here for reproducibility purposes. Users interested in using these data for alternative analyses should refer to the original data archives:
https://www.fs.usda.gov/rds/archive/Product/RDS-2011-0002.3, https://www.fs.usda.gov/rds/archive/Product/RDS-2011-0003.3, https://www.fs.usda.gov/rds/archive/Product/RDS-2013-0007-2, https://www.fs.usda.gov/rds/archive/Product/RDS-2015-0031, https://www.fs.usda.gov/rds/archive/Product/RDS-2015-0030, and https://www.fs.usda.gov/rds/archive/Product/RDS-2015-0032.

