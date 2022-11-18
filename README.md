
<!-- README.md is generated from README.Rmd. Please edit that file -->

# cellssm

<!-- badges: start -->
<!-- badges: end -->

## Overview

cellssm provides an easy way to analyze the time-series of distances of
cells or organelles from an external stimulus. Main functions are based
on the Bayesian inference of parameters in the state-space model.

## Installation

You can install the development version of cellssm from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("hnishio/cellssm")
```

## Example

A real data example of chloroplast accumulation responses to a blue
microbeam:

``` r
# Load packages
library(cellssm)

# Set the path to which CmdStan was installed
cmdstanr::set_cmdstan_path("~/cmdstan/")

# Load data of chloroplast movements
data("cell1", "cell2", "cell3", "cell4")
cell_list <- list(cell1[,1:3])
# This is for a chloroplast in a cell as an example.
# If you want to run the modeling for multiple data frames, execute:
# cell_list <- list(cell1, cell2, cell3, cell4)

# Execution of state-space modeling
ssm_individual(cell_list = cell_list, out = "02_ssm_individual",
               res_name = "chloroplast", ex_name = "microbeam",
               unit1 = "micrometer", unit2 = "min")
#> Running MCMC with 4 parallel chains...
#> 
#> Chain 1 Iteration:    1 / 6000 [  0%]  (Warmup) 
#> Chain 2 Iteration:    1 / 6000 [  0%]  (Warmup) 
#> Chain 3 Iteration:    1 / 6000 [  0%]  (Warmup) 
#> Chain 4 Iteration:    1 / 6000 [  0%]  (Warmup) 
#> Chain 4 Iteration: 1200 / 6000 [ 20%]  (Warmup) 
#> Chain 3 Iteration: 1200 / 6000 [ 20%]  (Warmup) 
#> Chain 1 Iteration: 1200 / 6000 [ 20%]  (Warmup) 
#> Chain 2 Iteration: 1200 / 6000 [ 20%]  (Warmup) 
#> Chain 3 Iteration: 2400 / 6000 [ 40%]  (Warmup) 
#> Chain 3 Iteration: 3001 / 6000 [ 50%]  (Sampling) 
#> Chain 1 Iteration: 2400 / 6000 [ 40%]  (Warmup) 
#> Chain 1 Iteration: 3001 / 6000 [ 50%]  (Sampling) 
#> Chain 2 Iteration: 2400 / 6000 [ 40%]  (Warmup) 
#> Chain 3 Iteration: 4200 / 6000 [ 70%]  (Sampling) 
#> Chain 1 Iteration: 4200 / 6000 [ 70%]  (Sampling) 
#> Chain 3 Iteration: 5400 / 6000 [ 90%]  (Sampling) 
#> Chain 3 Iteration: 6000 / 6000 [100%]  (Sampling) 
#> Chain 3 finished in 5.7 seconds.
#> Chain 2 Iteration: 3001 / 6000 [ 50%]  (Sampling) 
#> Chain 1 Iteration: 5400 / 6000 [ 90%]  (Sampling) 
#> Chain 1 Iteration: 6000 / 6000 [100%]  (Sampling) 
#> Chain 1 finished in 6.7 seconds.
#> Chain 2 Iteration: 4200 / 6000 [ 70%]  (Sampling) 
#> Chain 4 Iteration: 2400 / 6000 [ 40%]  (Warmup) 
#> Chain 4 Iteration: 3001 / 6000 [ 50%]  (Sampling) 
#> Chain 2 Iteration: 5400 / 6000 [ 90%]  (Sampling) 
#> Chain 2 Iteration: 6000 / 6000 [100%]  (Sampling) 
#> Chain 2 finished in 10.3 seconds.
#> Chain 4 Iteration: 4200 / 6000 [ 70%]  (Sampling) 
#> Chain 4 Iteration: 5400 / 6000 [ 90%]  (Sampling) 
#> Chain 4 Iteration: 6000 / 6000 [100%]  (Sampling) 
#> Chain 4 finished in 13.2 seconds.
#> 
#> All 4 chains finished successfully.
#> Mean chain execution time: 8.9 seconds.
#> Total execution time: 13.2 seconds.
```

<img src="man/figures/ssm_individual_cell1_chloroplast1.jpg" style="width:50.0%" />

This figure is an example of the output files. Dots, solid lines and
shaded regions are the observed values, median and 95% credible
intervals of the Bayesian inference, respectively. Orange solid lines
represent the start time estimated by the model. Shaded and light
regions represent the period without and with the explanatory variable,
respectively.
