% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ssm_individual.R
\name{ssm_individual}
\alias{ssm_individual}
\title{Bayesian inference of the state-space model (individual model)}
\usage{
ssm_individual(
  cell_list,
  visual = NULL,
  out,
  seed = 123,
  warmup = 1000,
  sampling = 1000,
  thin = 3,
  stepwise = c(95, 90),
  start_sensitivity = 5,
  ex_sign = "negative",
  df_name = "cell",
  res_name = "organelle",
  ex_name,
  df_idx = NULL,
  res_idx = NULL,
  graph = TRUE,
  graph_title = "Cell",
  unit1,
  unit2,
  shade = TRUE,
  start_line = TRUE,
  ps = 7,
  theme_plot = "bw",
  diagnosis = TRUE
)
}
\arguments{
\item{cell_list}{(list of data frame) The input time-series data. First column, time (column name "time");
second column, an explanatory variable (0 or 1, column name "ex"); third to the last columns,
distances of cells or organelles from the explanatory variable (
any column names are accepted). The velocities of the third to the last columns are
calculated from the distance and time, and used as response variables in the
modelling. See the following \strong{Examples} for further details.}

\item{visual}{(data frame) The optional data of visual estimation of the start time of
the influence of an explanatory variable. First column, cells (column name "cell");
second column, index (column name "index"); third column, start time. The default is \code{NULL}.}

\item{out}{(character string) The path of the output directory.}

\item{seed}{(positive integer) A seed for random number generation to pass to CmdStan.}

\item{warmup}{(positive integer) The number of warmup iterations of MCMC sampling after thinning.
The default is 1000.}

\item{sampling}{(positive integer) The number of post-warmup iterations of MCMC sampling after thinning.
The default is 1000.}

\item{thin}{(positive integer) Intervals of MCMC samples. This is useful for
reducing the autocorrelation of the MCMC samples and improving the convergence of
MCMC. The default is 3.}

\item{stepwise}{(integer vector of length two chosen from (95, 90, 80, 70, 60, 50))
The confidence intervals for stepwise estimation of the start time of the movement.
For example, when the first value is set to 95, the start of the movement is
estimated from the 95 \% confidence interval of the time-varying coefficient of
the explanatory variable. When the second value is set to the lower values and
the start time can not be estimated with the 95 \% confidence interval, the
threshold is sequentially lowered to the second value until the start time is
determined. If the start time is still not able to be determined, it is set to
infinity. The default is c(95, 90).}

\item{start_sensitivity}{(positive integer) The sensitivity to estimate the start time
of the movement. Larger values indicate a higher sensitivity and adopt the earlier
start time estimated with lower confidence intervals. In detail, this parameter
relates to how much earlier the start time estimated by the lower threshold of
confidence intervals is than that estimated by the 95 \% confidence interval.
When this difference is larger than the value obtained by dividing the period of
the explanatory variable by \code{start_sensitivity}, the start time estimated with
lower confidence intervals is adopted. When \code{start_sensitivity} is set to the
value lower than 1, no threshold lower than the 95 \% confidence interval is
considerred. The default is 5.}

\item{ex_sign}{(character string) "positive" or "negative". This is used to
estimate the start time of the positive or negative influence of the explanatory
variable on the distances of cells or organelles.
If cells or organelles are moving away from the explanatory variable, \code{ex_sign}
should be "positive". If cells or organelles are approaching to the explanatory
variable, \code{ex_sign} should be "negative".}

\item{df_name}{(character string) The name of the data frame. This is used for
file names. The default is "cell".}

\item{res_name}{(character string) The name of the response variable. This is used
for file names and graph labels. The default is "organelle".}

\item{ex_name}{(character string) The name of the explanatory variable. This is used
for graph labels.}

\item{df_idx}{(integer vector) Indexes of the data frame. This should be set
only when you want to set the indexes manually. This is used for
file names and graph titles. For example, setting "\code{df_idx} = c(1,3) and \code{res_idx}
= c(2,8)" with the default \code{df_name} and \code{res_name} results in the file names
including "cell1_organelle2" and "cell3_organelle8". The default is \code{NULL} and
the indexes are automatically set.}

\item{res_idx}{(integer vector) Indexes of the response variable. This should be set
only when you want to set the indexes manually. This is used for
file names and graph titles. For example, setting "\code{df_idx} = c(1,3) and \code{res_idx}
= c(2,8)" with the default \code{df_name} and \code{res_name} results in the file names
including "cell1_organelle2" and "cell3_organelle8". The default is \code{NULL} and
the indexes are automatically set.}

\item{graph}{(logical) Whether to output the graphs of the estimation.
See \strong{Value} for more details. The default is \code{TRUE}.}

\item{graph_title}{(character string) The name of the data frame. This is used for
graph titles. The default is "Cell".}

\item{unit1}{(character string) The unit of the response variable. One of "meter",
"centimeter", "millimeter", "micrometer", "nanometer". If another character
string is given, it is used as it is. This is used for graph labels.}

\item{unit2}{(character string) The unit of time. This is used for graph labels.}

\item{shade}{(logical) Whether to draw shade in graphs during the period without
the explanatory variable. The default is \code{TRUE}.}

\item{start_line}{(logical) Whether to draw a line at the start time of the influence of
the explanatory variable in graphs. The default is \code{TRUE}.}

\item{ps}{(positive integer) Font size of graphs specified in pt. The default is 7 pt.
Plot sizes are automatically adjusted according to the font size.}

\item{theme_plot}{(character string) A plot theme of the \link{ggplot2} package. One of "bw", "light",
"classic", "gray", "dark", "test", "minimal" and "void". The default is "bw".}

\item{diagnosis}{(logical) Whether to output the visualised diagnoses of
MCMC sampling. See \strong{Value} for more details. The default is \code{TRUE}.}
}
\value{
A directory named after the \code{out} parameter is created, which has three subdirectories.
\itemize{
\item A subdirectory "csv" includes "ssm_individual_cell \code{i} _ \code{res_name} \code{j} .csv",
"ssm_individual_cell \code{i} _ \code{res_name} \code{j} _sd.csv" and "ssm_individual_mvtime.csv"
with \code{i} the indexes of cells, \code{j} the indexes of the response variables.
"ssm_individual_cell \code{i} _ \code{res_name} \code{j} .csv" contains the Bayesian credible
intervals of time-varying parameters, where "Y" is the observed velocity, "w" is
the white noise, "alpha" is the true state of velocity, and "b_ex" is
the time-varying coefficient of the explanatory variable. "ssm_individual_cell \code{i} _ \code{res_name} \code{j} _sd.csv"
contains the Bayesian credible intervals of non-time-varying parameters, where
"s_w" and "s_b_ex" are the white noise and system noise, respectively.
"ssm_individual_mvtime.csv" contains the estimated start time, end time, and period of
the directional movement.
\item A subdirectory "pdf" includes "ssm_individual_cell \code{i} _ \code{res_name} \code{j} .pdf".
This is the visualised results of the model. The figure consists of four panels:
(1) observed distance of \code{res_name} from \code{ex_name}, (2) velocity of \code{res_name},
(3) regression coefficient of \code{ex_name}, and (4) random fluctuations
in the velocity. In these panels,
dots, solid lines, and shaded regions are the observed values, medians, and 95\%
credible intervals of the Bayesian inference, respectively. When the optional visual estimation of
the start time is given to the \code{visual} parameter, orange solid lines and
green dashed lines represent the start time estimated using the model and
the visual observation, respectively. When the \code{shade} parameter is \code{TRUE},
the shaded and light regions represent the periods without and with the explanatory
variable, respectively.
\item A subdirectory "diagnosis" includes "ssm_individual_combo_cell \code{i} _ \code{res_name} \code{j} .pdf" and
"ssm_individual_rhat_cell \code{i} _ \code{res_name} \code{j} .pdf". These are the visualised diagnoses of
MCMC sampling. "ssm_individual_combo_cell \code{i} _ \code{res_name} \code{j} .pdf" shows the
density plots of the posterior distributions and the trace plots for "b_ex" and
"alpha" at the start and end of the time series and the parameter with the worst Rhat value.
"ssm_individual_rhat_cell \code{i} _ \code{res_name} \code{j} .pdf" indicates the Rhat values of the parameters.
These are drawn using the \link{bayesplot} package.
}
}
\description{
\code{ssm_individual} performs the Bayesian inference of the parameters of the
state-space model to analyse the velocity of individual cells or organelles.
It estimates the Bayesian credible intervals of the parameters and the start time of
the influence of an explanatory variable in the regression model. Before using
this function, the cmdstanr package (>= 0.5.2)
(\url{https://mc-stan.org/cmdstanr/index.html}) and CmdStan
(\url{https://mc-stan.org/docs/cmdstan-guide/cmdstan-installation.html})
must be installed. Installation of CmdStan from GitHub would be easier for beginners.
Install Git and execute a few commands in the Terminal app on Mac or
the Git Bash app on Windows.
}
\examples{
# For the first-time usage, install the cmdstanr package (>= 0.5.2)
# (https://mc-stan.org/cmdstanr/index.html) and CmdStan
# (https://mc-stan.org/docs/cmdstan-guide/cmdstan-installation.html)


### Real data example of chloroplast accumulation responses to a blue microbeam ###

# Load packages
library(cellssm)

# Load data of chloroplast movements
data("cell1", "cell2", "cell3", "cell4", "visual")
cell_list <- list(cell1, cell2, cell3, cell4)

# Check the format of input data
cell_list
visual

\dontrun{
# Execution of state-space modelling

  # Set the path where CmdStan was installed
  cmdstanr::set_cmdstan_path("~/cmdstan/")

  # When you do not want to compare the statistical and visual estimations of the start time
  ssm_individual(cell_list = cell_list, out = "02_ssm_individual",
                 res_name = "chloroplast", ex_name = "microbeam",
                 unit1 = "micrometer", unit2 = "min")

  # When you do want to compare the statistical and visual estimations of the start time
  ssm_individual(cell_list = cell_list, visual = visual, out = "02_ssm_individual",
                 res_name = "chloroplast", ex_name = "microbeam",
                 unit1 = "micrometer", unit2 = "min")
}



### Simulated data example of Paramecium escape responses from laser heating ###

# Load packages
library(cellssm)

# Load data
data("Paramecium")
cell_list <- list(Paramecium)

# Check the format of input data
cell_list

\dontrun{
# Execution of state-space modelling

  # Set the path where CmdStan was installed
  cmdstanr::set_cmdstan_path("~/cmdstan/")

  ssm_individual(cell_list = cell_list, out = "12_ssm_individual",
                 warmup=1000, sampling=1000, thin=6,
                 start_sensitivity = 3, ex_sign = "positive", df_name = "experiment",
                 res_name = "Paramecium", ex_name = "heat",
                 graph_title = "Experiment", unit1 = "millimeter", unit2 = "sec")
}

}
