# betadiversity 
Code to regress estimates of beta diversity against three sets of predictor variables: spatial scale, functional groupings of species, and environmental covariates. Regressions are based on beta regression models, implemented in Stan through the rstan R package.

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

Updated: 13 March 2018

## Usage
Several scripts are provided. The main model runs are provided in beta-div-[predictor].R files, where [predictor] is the predictor variable of interest. These scripts require a pre-loaded data set (to be made available on publication). Several helper scripts are provided, including helpers.R, install_packages.R, and prepare-lm-data.R. The load-all-data.R script loads all data from scratch, calling on many individual data files and two helper scripts (load-data-helpers.R and load-pa-data.R). Raw data files will be made available on publication, alongside a pre-loaded data set (loading data from raw files takes approximately half an hour on a MacBook Air).

The [filename].stan files are Stan models with spatial scale or functional groupings (guilds) as predictor variables. Models with environmental predictor variables are created dynamically using the stan_gen_beta_mod() function in helpers.R.


