% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lm_dist_start.R
\name{lm_dist_start}
\alias{lm_dist_start}
\title{Linear regression (x: distance at time 0, y: start time)}
\usage{
lm_dist_start(
  cell_list,
  mvtime,
  robust = FALSE,
  graph_title = "Cell",
  ex_name,
  unit1,
  unit2,
  eq_pos = "topleft",
  ps = 7,
  theme_plot = "bw"
)
}
\arguments{
\item{cell_list}{(list of data frame) The input time-series data. First column, time (column name "time");
second column, an explanatory variable (0 or 1, column name "ex"); third to the last columns,
distances of cells or organelles from the explanatory variable (
any column names are accepted). See the following \strong{Examples} for further details.}

\item{mvtime}{(data frame) The movement time estimated by \link{ssm_individual}
("ssm_individual_mvtime.csv") or \link{ssm_KFAS} ("ssm_KFAS_mvtime.csv").}

\item{robust}{(logical) When \code{TRUE}, performs robust regression using repeated medians
(Siegel regression). When \code{FALSE}, performs linear regression. The default is \code{FALSE}.}

\item{graph_title}{(character string) The name of the data frame. This is used for
graph titles. The default is "Cell".}

\item{ex_name}{(character string) The name of the explanatory variable. This is used for graph labels.}

\item{unit1}{(character string) The unit of the response variable. One of "meter",
"centimeter", "millimeter", "micrometer", "nanometer". If another character
string is given, it is used as it is. This is used for graph labels.}

\item{unit2}{(character string) The unit of time. This is used for graph labels.}

\item{eq_pos}{(character string) The position of equations. One of "topleft", "topright",
"bottomleft", "bottomright". The default is "topleft".}

\item{ps}{(positive integer) Font size of graphs specified in pt. The default is 7 pt.
Plot sizes are automatically adjusted according to the font size.}

\item{theme_plot}{(character string) A plot theme of the \link{ggplot2} package. One of "bw", "light",
"classic", "gray", "dark", "test", "minimal" and "void". The default is "bw".}
}
\value{
A list of ggplot objects is returned. In each plot, dots, solid lines
and shaded regions are the observed values, regression lines, and 95\% confidence
intervals, respectively.
}
\description{
\code{lm_dist_start} performs linear regression to analyse the
relationship between the distance at time 0 and the start time of the
influence of an explanatory variable obtained by the state-space model.
}
\examples{
### Real data example of chloroplast accumulation responses to a blue microbeam ###

# Load package
library(cellssm)

# Load data
data("cell1", "cell2", "cell3", "cell4", "chloroplast_mvtime")
cell_list <- list(cell1, cell2, cell3, cell4)

# dist vs. start
glist <- lm_dist_start(cell_list = cell_list, mvtime = chloroplast_mvtime,
                       ex_name = "microbeam", unit1 = "micrometer", unit2 = "min")

g <- glist[[1]] + glist[[2]] + glist[[3]] + glist[[4]] + glist[[5]] +
  patchwork::plot_layout(ncol = 3)

\dontrun{
# Create an output directory
out <- "06_lm_dist_start"
if(file.exists(out)==FALSE){
  dir.create(out, recursive=TRUE)
}

# Save output
suppressWarnings(ggplot2::ggsave(paste0(out, "/individual_chloroplast_lm_dist_start.pdf"),
                                 g, height = 110, width = 50*3, units = "mm"))
}



### Simulated data example of Paramecium escape responses from a laser heating ###

# Load package
library(cellssm)

# Load data of chloroplast movements
data("Paramecium", "Paramecium_mvtime")
cell_list <- list(Paramecium)

# dist vs. start
glist <- lm_dist_start(cell_list = cell_list, mvtime = Paramecium_mvtime,
                       graph_title = "Experiment", ex_name = "heat",
                       unit1 = "millimeter", unit2 = "sec")

\dontrun{
# Create an output directory
out <- "16_lm_dist_start"
if(file.exists(out)==FALSE){
  dir.create(out, recursive=TRUE)
}

# Save output
g <- glist[[1]]
suppressWarnings(ggplot2::ggsave(paste0(out, "/individual_Paramecium_lm_dist_start.pdf"),
                                 g, height = 50, width = 50, units = "mm"))
}

}
