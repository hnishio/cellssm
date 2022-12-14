% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nomodel.R
\name{nomodel}
\alias{nomodel}
\title{Estimation of movement without the state-space model}
\usage{
nomodel(
  cell_list,
  visual = NULL,
  out,
  ex_sign = "negative",
  sma_period = 10,
  fold = 2,
  df_name = "cell",
  res_name,
  ex_name,
  unit1,
  unit2,
  shade = TRUE,
  start_line = TRUE,
  ps = 7,
  theme_plot = "bw"
)
}
\arguments{
\item{cell_list}{(list of data frame) The input time-series data. First column, time (column name "time");
second column, an explanatory variable (0 or 1, column name "ex"); third to the last columns,
distances of cells or organelles from the explanatory variable (
any column names are accepted). See the following \strong{Examples} for further details.}

\item{visual}{(data frame) The optional data of visual estimation of the start time of
the influence of an explanatory variable. First column, cells (column name "cell");
second column, index (column name "index"); third column, the start time. The default is \code{NULL}.}

\item{out}{(character string) The path of the output directory.}

\item{ex_sign}{(character string) "positive" or "negative". This is used to
estimate the start time of the positive or negative influence of the explanatory
variable on the distances of cells or organelles.
If cells or organelles are moving away from the explanatory variable, \code{ex_sign}
should be "positive". If cells or organelles are approaching to the explanatory
variable, \code{ex_sign} should be "negative".}

\item{sma_period}{(positive integer) One of the definitions of the start time:
the simple moving averages of the first-order differences for \code{sma_period} are
\code{ex_sign} (positive or negative). The default is 10.}

\item{fold}{(positive real) One of the definitions of the start time:
difference from the distance at five time points ahead is \code{fold} times larger
than the average change rate. The default is 2.}

\item{df_name}{(character string) The name of the data frame. This is used for
file names and graph titles. The default is "cell".}

\item{res_name}{(character string) The name of the response variable. This is used
for file names and graph labels.}

\item{ex_name}{(character string) The name of the explanatory variable. This is
used for graph labels.}

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
}
\value{
A directory named after the \code{out} parameter is created.

"nomodel_cell \code{i} _ \code{res_name} \code{j} .pdf" (with \code{i} the indexes of cells,
\code{j} the indexes of the response variables) is the visualised result of the estimation.
The observed distance of \code{res_name} from \code{ex_name} is shown as a solid line.
When the optional visual estimation of the start time is given to the \code{visual}
parameter, orange solid lines and green dashed lines represent the start time
estimated by the model and the visual observation, respectively. When the \code{shade}
parameter is \code{TRUE}, the shaded and light regions represent the period without
and with the explanatory variable, respectively.

"nomodel_mvtime.csv" contains the estimated start time, end time, and period of
the directional movement.
}
\description{
\code{nomodel} estimates the start time of the movement as the time
when the time-varying distances of cells or organelles from
the explanatory variable is increasing or decreasing both in the short term and
in the long term.
}
\examples{
### Real data example of chloroplast accumulation responses to a blue microbeam ###

# Load package
library(cellssm)

# Load data of chloroplast movements
data("cell1", "cell2", "cell3", "cell4", "visual")
cell_list <- list(cell1, cell2, cell3, cell4)

# Check the format of input data
cell_list
visual

\dontrun{
# Predict movement

# When you do not want to compare the statistical and visual estimations of the start time
nomodel(cell_list = cell_list, out = "04_nomodel",
        res_name = "chloroplast", ex_name = "microbeam",
        unit1 = "micrometer", unit2 = "min")

# When you do want to compare the statistical and visual estimations of the start time
nomodel(cell_list = cell_list, visual = visual, out = "04_nomodel",
        res_name = "chloroplast", ex_name = "microbeam",
        unit1 = "micrometer", unit2 = "min")
}



### Simulated data example of Paramecium escape responses from laser heating ###

# Load package
library(cellssm)

# Load data
data("Paramecium")
cell_list <- list(Paramecium)

# Check the format of input data
cell_list

\dontrun{
# Predict movement
nomodel(cell_list = cell_list, out = "14_nomodel",
        ex_sign = "positive", fold = 1, df_name = "experiment",
        res_name = "Paramecium", ex_name = "heat",
        unit1 = "millimeter", unit2 = "sec")
}

}
