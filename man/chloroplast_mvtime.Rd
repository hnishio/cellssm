% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/chloroplast_mvtime.R
\docType{data}
\name{chloroplast_mvtime}
\alias{chloroplast_mvtime}
\title{Model estimation of the start and end times of chloroplast accumulation responses in
four plant cells}
\usage{
data("chloroplast_mvtime")
}
\description{
A real dataset containing the start and end times of chloroplast
movements in response to a blue microbeam, estimated by \link{ssm_individual}.
}
\details{
A data frame with 50 rows and 5 columns.
}
\section{Columns}{

\itemize{
\item cell : index of cells (1--4)
\item chloroplast : index of chloroplasts in each cell
\item start_time : start time of chloroplast movements in minutes
\item end_time : end time of chloroplast movements in minutes
\item move_time : duration of chloroplast movements in minutes
}
}

\keyword{datasets}
