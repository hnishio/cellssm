
<!-- README.md is generated from README.Rmd. Please edit that file -->

# cellssm: State-Space Modeling for the Directional Movement of Cells

<!-- badges: start -->
<!-- badges: end -->

## Overview

An easy way for molecular biologists to analyse the time series of
distances of cells or organelles from an external stimulus. Using this
package, one can estimate the true dynamics from noisy movement data,
extract the common dynamics among multiple cells or organelles, and
estimate the start time of the movements.

## Details

The package only requires the preparation of csv files containing the
columns time (min, sec, etc.), presence of stimulus (0 or 1), and
distances of cells or organelles from the stimulus (millimeter,
micrometer, etc.). The Main functions are based on the Bayesian
inference of the parameters in the state-space model, including the
time-varying coefficient of regression. This package consists of main
functions to implement state-space modelling, and sub-functions to
perform minor tasks.

#### Main functions:

-   ssm_individual : Bayesian inference of the state-space model for
    individual dynamics
-   ssm_common : Bayesian inference of the state-space model to extract
    common dynamics
-   ssm_KFAS : Inference of the state-space model for individual
    dynamics using Kalman filter
-   nomodel : Estimation of movement without the state-space model

## Installation

You can install the development version of cellssm from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("hnishio/cellssm")
```

## Example

Real data example of chloroplast accumulation responses to a blue
microbeam:

``` r
# Load package
library(cellssm)

# Set the path to which CmdStan was installed
cmdstanr::set_cmdstan_path("~/cmdstan/")

# Load data of chloroplast movements
data("cell1", "cell2", "cell3", "cell4")
cell_list <- list(cell1, cell2, cell3, cell4)

# Execution of state-space modeling
ssm_individual(cell_list = cell_list, out = "02_ssm_individual",
               res_name = "chloroplast", ex_name = "microbeam",
               unit1 = "micrometer", unit2 = "min")
```

<img src="man/figures/ssm_individual_cell1_chloroplast1.jpg" width="40.0%" />

This figure is an example of the output files. Dots, solid lines and
shaded regions are the observed values, medians, and 95% credible
intervals of the Bayesian inference, respectively. Orange solid lines
represent the start time estimated by the model. Shaded and light
regions represent the period without and with the explanatory variable,
respectively.
