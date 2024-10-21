
# Registry data simulations for the analysis of GLP1 use and dementia in a cohort of diabetes patients
Reproducible analysis pipeline for a longitudinal simulation study of the effect of 2nd line diabetes drugs from EHR data similar to the Danish registry. Simulation results found in the arXiv paper [here](https://arxiv.org/abs/2310.03235).


## Overview

This repository contains a pipeline for longitudinal simulation studies from EHR data similar to the Danish registry. The pipeline is adapted from the applied analysis functions used to analyze the Danish registry, which use the `targets` package framework and reproducible pipeline. The simulations study, however, does not use a `targets` pipeline itself due to computation time and for increased flexibility in running parts of the simulation. It uses a series of scripts that can be run in order to reproduce the simulation study using the same functions used in the applied analysis. The pipeline is designed to be reproducible and to be run or adapted by any reader of the manuscript who has access to a high-performance computing cluster, as the whole simulation is computationally expensive and prohibitively slow without parallelization.


## Table of Contents

1. [Installation](#installation)
2. [Functions](#functions)
3. [Usage](#usage)
4. [License](#license)


## Installation

1. Install required packages by running:

```R
   install.packages(c("tidyverse", "targets", "fst","lava","ltmle","data.table","tidyverse","glmnet","Matrix","Publish","matrixStats","speedglm","parallel","caret","foreach","clustermq"))
```

## Functions

The simulation study primarily runs off functions sourced from the Ltmle folder, which are the `ltmle` package functions for longitudinal maximum likelihood estimation, adapted to run without error on Statistics Denmark, the remote server location of the Danish Registry and to allow different estimator specifications outside of GLM and SuperLearner. Other functions in the /functions/ folder are wrapper functions to run the simulation and applied analysis without error on Statistics Denmark as well as the remote server used to run the simulation.


## Usage

The scripts to reproduce the simulation are in the /src/ folder labeled in order, and the script `0_run_simulations.R` runs each script in order. The file contents, in brief, are as follows:

* `1_simulate_data.R`: using the lava package and exported coefficients using the Danish registry, we simulate the datasets to be used for the simulation and save this in the /data/ subfolder. 
* `2_run_simulation_for_point_estimators.R`: use function `calc_realistic_truth` to intervene on the DGP to set desired causal contrasts and simulate from them with large N, approximating the truth.
* `3_run_null_simulations_for_point_estimates.R`:
* `4_calc_performance.R`:
* `5_calc_performance_null.R`:
* `6_run_bootstrap_variance.R`:
* `7_calc_bootstrap_performance.R`:
* `8_create_latex_tables.R`:
* `9_example_DAG.R`: creates the example DAG used in the publication. 

``` R
source(here("1_simulate_data.R"))
source(here("2_run_simulation_for_point_estimators.R"))
source(here("3_run_null_simulations_for_point_estimates.R"))
source(here("4_calc_performance.R"))
source(here("5_calc_performance_null.R"))
source(here("6_run_bootstrap_variance.R"))
source(here("7_calc_bootstrap_performance.R"))
source(here("8_create_latex_tables.R"))
source(here("9_example_DAG.R"))
```

## Notes
- The deterministic Q function (which encodes the setting of the probability of the outcome dementia deterministically to 1 once death, the competing event, has occurred) is built into the updates to the ltmle package, found in the repo. The exact code used for this can be found here: Ltmle/Augmentation/event_node_manipulator.R

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details.

---

For any questions or issues, please open an issue on the GitHub repository or contact the maintainer at amertens@berkeley.edu.
