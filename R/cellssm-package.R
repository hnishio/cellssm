#' cellssm: State-Space Modelling for the Directional Movement of Cells
#'
#' @description An easy way for molecular biologists to analyse the time series of
#' distances of cells or organelles from an external stimulus. Using this package,
#' one can estimate the true dynamics from noisy movement data, extract the common
#' dynamics among multiple cells or organelles, and estimate the start time of the
#' movements.
#'
#' @details The package only requires the preparation of csv files containing the
#' columns time (min, sec, etc.), presence of stimulus (0 or 1), and distances of
#' cells or organelles from the stimulus (mm, \eqn{\mu}m, etc.). The Main functions
#' are based on the Bayesian inference of the parameters in the state-space model,
#' including the time-varying coefficient of regression. This package consists of
#' main functions to implement state-space modelling, and sub-functions to perform
#' minor tasks.
#'
#' @section Main functions:
#' * [ssm_individual] : Bayesian inference of the state-space model for individual dynamics
#' * [ssm_common] : Bayesian inference of the state-space model to extract common dynamics
#' * [ssm_KFAS] : Inference of the state-space model for individual dynamics using Kalman filter
#' * [nomodel] : Estimation of movement without the state-space model
#'
#' @section Sub-functions:
#' * [dist_vis] : Visualisation of the distance from an explanatory variable
#' * [lm_dist_beta] : Robust linear regression (x: distance at time 0, y: coefficient of explanatory variable)
#' * [lm_dist_start] : Robust linear regression (x: distance at time 0, y: start time)
#' * [lm_signal] : Estimation of signal transfer speed by robust linear regression (x: start time, y: distance at time 0)
#'
#' @section Datasets:
#' * [cell1] : Time series of chloroplast accumulation responses in cell 1
#' * [cell2] : Time series of chloroplast accumulation responses in cell 2
#' * [cell3] : Time series of chloroplast accumulation responses in cell 3
#' * [cell4] : Time series of chloroplast accumulation responses in cell 4
#' * [visual] : Visual estimation of the start time of chloroplast accumulation responses in
#' four plant cells
#' * [chloroplast_mvtime] : Model estimation of the start and end times of chloroplast accumulation responses in
#' the four plant cells
#' * [Paramecium] : Time series of \emph{Paramecium} escape responses
#' * [Paramecium_mvtime] : Model estimation of the start and end times of \emph{Paramecium} escape responses
#'
#' @author Haruki Nishio, \email{harukin218@@gmail.com}, \url{https://orcid.org/0000-0002-6124-782X}
#' @author Satoyuki Hirano, \email{mc226592@@s.utsunomiya-u.ac.jp}
#' @author Yutaka Kodama, \email{kodama@@cc.utsunomiya-u.ac.jp}
#' @docType package
#' @name cellssm-package
#' @aliases cellssm
#' @import ggplot2
#' @import patchwork
#'
NULL
