
### Quantile Function
quantile99 <- function(x){
  stats::quantile(x, probs = c(0.005, 0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.975, 0.995), names = TRUE)
}



#' Bayesian inference of the state-space model (individual model)
#'
#' \code{ssm_individual} performs the Bayesian inference of the parameters of the
#' state-space model to analyze the velocity of individual cell or organelle.
#' It estimates the Bayesian credible intervals of parameters and the start time of
#' the influence of an explanatory variable in a regression model. Before using
#' this function, you need to install the cmdstanr package (>= 0.5.2)
#' (\url{https://mc-stan.org/cmdstanr/index.html}) and CmdStan
#' (\url{https://mc-stan.org/docs/cmdstan-guide/cmdstan-installation.html}).
#' To install CmdStan, installation from GitHub would be easier for beginners:
#' you need to install Git and execute a few commands in the Terminal app on Mac or
#' the Git Bash app on Windows.
#'
#' @param cell_list (list of data frame) The input time-series data. First column, time (column name "time");
#' second column, an explanatory variable (0 or 1, column name "ex"); third to the last columns,
#' distances of cells or organelles from the explanatory variable (
#' any column names are accepted). The velocity of third to the last columns is
#' calculated from the distance and time, and used as response variables in the
#' modeling. See examples below for more details.
#' @param visual (data frame) The optional data of visual estimation of the start time of
#' the influence of an explanatory variable. First column, cells (column name "cell");
#' second column, index (column name "index"); third column, the start time. The default is `NULL`.
#' @param out (character string) The path of the output directory.
#' @param warmup (positive integer) The number of warmup iterations of MCMC sampling after thinning.
#' The default is 1000.
#' @param sampling (positive integer) The number of post-warmup iterations of MCMC sampling after thinning.
#' The default is 1000.
#' @param thin (positive integer) Interval of MCMC sampling. This is useful to reduce autocorrelation of
#' MCMC samples and improve convergence of MCMC. The default is 3.
#' @param start_sensitivity (positive integer) The sensitivity to detect the start time
#' of movements. Larger values indicate higher sensitivities. The default is 5.
#' @param ex_sign (character string) "positive" or "negative". This is used to
#' estimate the start time of the positive or negative influence of the explanatory
#' variable on the distances of cells or organelles.
#' If cells or organelles are moving away from the explanatory variable, `ex_sign`
#' should be "positive". If cells or organelles are approaching to the explanatory
#' variable, `ex_sign` should be "negative".
#' @param df_name (character string) The name of the data frame. This is used for
#' file names and graph titles. The default is "cell".
#' @param res_name (character string) The name of the response variable. This is used
#' for file names and graph labels.
#' @param ex_name (character string) The name of the explanatory variable. This is used
#' for graph labels.
#' @param unit1 (character string) The unit of a response variable. One of "meter",
#' "centimeter", "millimeter", "micrometer", "nanometer". If another character
#' string is given, it is used as it is. This is used for graph labels.
#' @param unit2 (character string) The unit of time. This is used for graph labels.
#' @param shade (logical) Whether to draw shade in graphs during the absence of
#' the explanatory variable. The default is `TRUE`.
#' @param start_line (logical) Whether to draw a line at the start time of the influence of
#' the explanatory variable in graphs. The default is `TRUE`.
#' @param ps (positive integer) Font size of graphs specified in pt. The default is 7 pt.
#' Plot sizes are automatically adjusted according to the font size.
#' @param theme_plot (character string) A ggplot theme. One of "bw", "light",
#' "classic", "gray", "dark", "test", "minimal" and "void". The default is "bw".
#' @returns A directory named after the `out` parameter is created, which have three subdirectories.
#' * A subdirectory "csv" includes "ssm_individual_cell `i` _ `res_name` `j` .csv",
#' "ssm_individual_cell `i` _ `res_name` `j` _sd.csv" and "ssm_individual_mvtime.csv"
#' with `i` the indexes of cells, `j` the indexes of the response variables.
#' "ssm_individual_cell `i` _ `res_name` `j` .csv" contains the Bayesian credible
#' intervals of time-varying parameters, where "Y" is the observed velocity, "w" is
#' the white noise, "alpha" is the true state of velocity and "b_ex" is
#' the time-varying coefficient of the explanatory variable. "ssm_individual_cell `i` _ `res_name` `j` _sd.csv"
#' contains the Bayesian credible intervals of non-time-varying parameters, where
#' "s_w", "s_b_ex" and "s_Y" are the white noise, system noise and observation error
#' as standard deviations, respectively.
#' "ssm_individual_mvtime.csv" contains the estimated start time, end time and period of
#' the directional movement.
#' * A subdirectory "pdf" includes "ssm_individual_cell `i` _ `res_name` `j` .pdf".
#' This is the visualized results of the model. The figure consists of four panels,
#' (1) the observed distance of `res_name` from `ex_name`, (2) the velocity of `res_name`,
#' (3) the regression coefficient of `ex_name`, (4) the random fluctuations
#' of the velocity. In these panels,
#' dots, solid lines and shaded regions are the observed values, median and 95%
#' credible intervals of the Bayesian inference, respectively. When the optional visual estimation of
#' the start time is given to the `visual` parameter, orange solid lines and
#' green dashed lines represent the start time estimated by the model and
#' the visual observation, respectively. When the `shade` parameter is `TRUE`,
#' the shaded and light regions represent the period without and with the explanatory
#' variable, respectively.
#' * A subdirectory "diagnosis" includes "ssm_individual_combo_cell `i` _ `res_name` `j` .pdf" and
#' "ssm_individual_rhat_cell `i` _ `res_name` `j` .pdf". These are the visualized diagnosis of
#' MCMC sampling. "ssm_individual_combo_cell `i` _ `res_name` `j` .pdf" shows the
#' density plots of the posterior distributions and the trace plots for "b_ex" and
#' "alpha" at the start and end of time series and the parameter with the worst Rhat value.
#' "ssm_individual_rhat_cell `i` _ `res_name` `j` .pdf" shows the Rhat values of parameters.
#' These are drawn using the [bayesplot] package.
#' @examples
#' # For the first time usage, install the cmdstanr package (>= 0.5.2)
#' # (https://mc-stan.org/cmdstanr/index.html) and CmdStan
#' # (https://mc-stan.org/docs/cmdstan-guide/cmdstan-installation.html)
#'
#'
#' ### A real data example of chloroplast accumulation responses to a blue microbeam ###
#'
#' # Load packages
#' library(cellssm)
#'
#' # Set the path to which CmdStan was installed
#' cmdstanr::set_cmdstan_path("~/cmdstan/")
#'
#' # Load data of chloroplast movements
#' data("cell1", "cell2", "cell3", "cell4", "visual")
#' cell_list <- list(cell1)
#' # If you want to run the modeling for multiple data frames, execute:
#' # cell_list <- list(cell1, cell2, cell3, cell4)
#'
#' # Check the format of the input data
#' cell_list
#' visual
#'
#' # Execution of state-space modeling
#'
#' \dontrun{
#' # When you don't want to compare the statistical and visual estimation of the start time
#' ssm_individual(cell_list = cell_list, out = "02_ssm_individual",
#'                res_name = "chloroplast", ex_name = "microbeam",
#'                unit1 = "micrometer", unit2 = "min")
#'
#' # When you want to compare the statistical and visual estimation of the start time
#' ssm_individual(cell_list = cell_list, visual = visual, out = "02_ssm_individual",
#'                res_name = "chloroplast", ex_name = "microbeam",
#'                unit1 = "micrometer", unit2 = "min")
#' }
#'
#'
#' ### A simulated data example of Paramecium escape responses from a laser heating ###
#'
#' # Load packages
#' library(cellssm)
#'
#' # Set the path to which CmdStan was installed
#' cmdstanr::set_cmdstan_path("~/cmdstan/")
#'
#' # Load data
#' data("Paramecium")
#' cell_list <- list(Paramecium)
#'
#' # Check the format of the input data
#' cell_list
#'
#' # Execution of state-space modeling
#'
#' \dontrun{
#' ssm_individual(cell_list = cell_list, out = "12_ssm_individual",
#'                warmup=1000, sampling=1000, thin=6,
#'                start_sensitivity = 3, ex_sign = "positive", df_name = "experiment",
#'                res_name = "Paramecium", ex_name = "heat",
#'                unit1 = "millimeter", unit2 = "sec")
#' }
#'
#' @export
#'
ssm_individual <- function(cell_list, visual=NULL, out, warmup=1000, sampling=1000, thin=3,
                           start_sensitivity = 5, ex_sign = "negative", df_name = "cell",
                           res_name, ex_name, unit1, unit2,
                           shade = TRUE, start_line = TRUE, ps = 7, theme_plot = "bw"){


  ## Binding variables locally to the function
  index <- time <- `alpha_2.5%` <- `alpha_97.5%` <- `alpha_50%` <- Y <-
    `b_ex_2.5%` <- `b_ex_97.5%` <- `b_ex_50%` <-
    `w_2.5%` <- `w_97.5%` <- `w_50%` <- NULL


  ## Set up

  # Create output directories
  if(file.exists(paste0(out, "/csv"))==F){
    dir.create(paste0(out, "/csv"), recursive=T)
  }
  if(file.exists(paste0(out, "/pdf"))==F){
    dir.create(paste0(out, "/pdf"), recursive=T)
  }
  if(file.exists(paste0(out, "/diagnosis"))==F){
    dir.create(paste0(out, "/diagnosis"), recursive=T)
  }
  if(file.exists("tmp")==F){
    dir.create("tmp", recursive=T)
    output_dir <- "tmp"
  }else{
    dir.create("tmp99", recursive=T)
    output_dir <- "tmp99"
  }

  # Prepare a container for movement time
  df_mv <- data.frame(NULL)

  # Compile stan file
  stan_file <- system.file("extdata", "individual_model.stan", package = "cellssm")
  model <- cmdstanr::cmdstan_model(stan_file)


  ## Execution of the Bayesian inference

  # for loop of cells
  for(i in 1:length(cell_list)){

    # for loop of response variable
    for (j in 1:(ncol(cell_list[[i]])-2)){

      # Prepare data_list
      data_list <- list(
        N = nrow(cell_list[[i]])-1,
        N_ex = length(which(cell_list[[i]]$ex==1)),
        ex = cell_list[[i]]$ex[-1],
        Y = diff(cell_list[[i]][,j+2])
      )

      # Get the boundary indexes of ex
      zero_idx <- which(data_list$ex == 0)
      boundary1 <- zero_idx[which(diff(zero_idx) != 1)]
      boundary2 <- zero_idx[which(diff(zero_idx) != 1) + 1]

      if(length(zero_idx[which(diff(zero_idx) != 1)]) < 1){
        boundary1 <- max(zero_idx)
        boundary2 <- data_list$N + 1
      }

      # Modify data_list
      data_list <- c(data_list, list(boundary1=boundary1, boundary2=boundary2))

      # Execute MCMC
      fit <- model$sample(
        data = data_list,
        seed = 1,
        iter_warmup = warmup*thin,
        iter_sampling = sampling*thin,
        chains = 4,
        parallel_chains = 4,
        refresh = floor(warmup/2.5*thin),
        #show_messages = F,
        #sig_figs = 4,
        output_dir = output_dir,
        output_basename = paste0("cell", i, "_", res_name, j),
        adapt_delta = 0.95,
        thin = thin
      )

      # 99% Bayesian credible intervals
      outcsv_name <- list.files(output_dir)
      outcsv_name <- outcsv_name[grep(paste0("cell", i, "_", res_name, j, "-"), outcsv_name)]
      tmp_csv_w <- NULL
      tmp_csv_b_ex <- NULL
      tmp_csv_alpha <- NULL
      tmp_csv_s_w <- NULL
      tmp_csv_s_b_ex <- NULL
      tmp_csv_s_Y <- NULL

      for(k in 1:length(outcsv_name)){
        tmp_csv <- as.data.frame(data.table::fread(cmd = paste0("grep -v '^#' ", output_dir, "/", outcsv_name[k])))
        tmp_csv_w <- rbind(tmp_csv_w, tmp_csv[,stringr::str_starts(names(tmp_csv), "w")])
        tmp_csv_b_ex <- rbind(tmp_csv_b_ex, tmp_csv[,stringr::str_starts(names(tmp_csv), "b_ex")])
        tmp_csv_alpha <- rbind(tmp_csv_alpha, tmp_csv[,stringr::str_starts(names(tmp_csv), "alpha")])
        tmp_csv_s_w <- c(tmp_csv_s_w, tmp_csv[,stringr::str_starts(names(tmp_csv), "s_w")])
        tmp_csv_s_b_ex <- c(tmp_csv_s_b_ex, tmp_csv[,stringr::str_starts(names(tmp_csv), "s_b_ex")])
        tmp_csv_s_Y <- c(tmp_csv_s_Y, tmp_csv[,stringr::str_starts(names(tmp_csv), "s_Y")])
      }

      df_w <- as.data.frame(t(apply(tmp_csv_w, 2, quantile99)))
      df_b_ex <- as.data.frame(t(apply(tmp_csv_b_ex, 2, quantile99)))
      df_alpha <- as.data.frame(t(apply(tmp_csv_alpha, 2, quantile99)))
      df_s <- t(data.frame(s_w = quantile99(tmp_csv_s_w),
                           s_b_ex = quantile99(tmp_csv_s_b_ex),
                           s_Y = quantile99(tmp_csv_s_Y)))
      df_s <- cbind(data.frame(s_name = row.names(df_s)), df_s)

      # Save output
      colnames(df_w) <- paste0("w_", colnames(df_w))
      colnames(df_b_ex) <- paste0("b_ex_", colnames(df_b_ex))
      colnames(df_alpha) <- paste0("alpha_", colnames(df_alpha))
      df <- cbind(data.frame(time = cell_list[[i]]$time[-1],
                             ex = cell_list[[i]]$ex[-1],
                             Y = diff(cell_list[[i]][,j+2])),
                  df_w, df_alpha,
                  as.data.frame(rbind(matrix(NA, nrow = data_list$boundary1, ncol = 21), as.matrix(df_b_ex), matrix(NA, nrow = data_list$N - data_list$boundary2 + 1, ncol = 21))))
      data.table::fwrite(df, file = paste0(out, "/csv/ssm_individual_", df_name, i, "_", res_name, j, ".csv"))
      data.table::fwrite(df_s, file = paste0(out, "/csv/ssm_individual_", df_name, i, "_", res_name, j, "_sd.csv"))


      ## Movement time
      ex_period <- max(cell_list[[i]]$time[cell_list[[i]]$ex == 1], na.rm = T) - min(cell_list[[i]]$time[cell_list[[i]]$ex == 1], na.rm = T) + 1

      if(ex_sign == "positive"){  # ex_sign == "positive"

        if(length(which(df$`b_ex_2.5%` > 0)) > 0){
          start_time <- df$time[which(df$`b_ex_2.5%` > 0)][1]
          end_time <- df$time[which(df$`b_ex_2.5%` > 0)][length(which(df$`b_ex_2.5%` > 0))] + 1
          move_time <- end_time - start_time
          if((start_time - df$time[which(df$`b_ex_5%` > 0)][1]) > ex_period/start_sensitivity){
            start_time <- df$time[which(df$`b_ex_5%` > 0)][1]
            end_time <- df$time[which(df$`b_ex_5%` > 0)][length(which(df$`b_ex_5%` > 0))] + 1
            move_time <- end_time - start_time
          }
          if((start_time - df$time[which(df$`b_ex_10%` > 0)][1]) > ex_period/start_sensitivity){
            start_time <- df$time[which(df$`b_ex_10%` > 0)][1]
            end_time <- df$time[which(df$`b_ex_10%` > 0)][length(which(df$`b_ex_10%` > 0))] + 1
            move_time <- end_time - start_time
          }
          if((start_time - df$time[which(df$`b_ex_15%` > 0)][1]) > ex_period/start_sensitivity){
            start_time <- df$time[which(df$`b_ex_15%` > 0)][1]
            end_time <- df$time[which(df$`b_ex_15%` > 0)][length(which(df$`b_ex_15%` > 0))] + 1
            move_time <- end_time - start_time
          }
          if((start_time - df$time[which(df$`b_ex_20%` > 0)][1]) > ex_period/start_sensitivity){
            start_time <- df$time[which(df$`b_ex_20%` > 0)][1]
            end_time <- df$time[which(df$`b_ex_20%` > 0)][length(which(df$`b_ex_20%` > 0))] + 1
            move_time <- end_time - start_time
          }
          if((start_time - df$time[which(df$`b_ex_25%` > 0)][1]) > ex_period/start_sensitivity){
            start_time <- df$time[which(df$`b_ex_25%` > 0)][1]
            end_time <- df$time[which(df$`b_ex_25%` > 0)][length(which(df$`b_ex_25%` > 0))] + 1
            move_time <- end_time - start_time
          }

        }else if(length(which(df$`b_ex_5%` > 0)) > 0){
          start_time <- df$time[which(df$`b_ex_5%` > 0)][1]
          end_time <- df$time[which(df$`b_ex_5%` > 0)][length(which(df$`b_ex_5%` > 0))] + 1
          move_time <- end_time - start_time
          if((start_time - df$time[which(df$`b_ex_10%` > 0)][1]) > ex_period/start_sensitivity){
            start_time <- df$time[which(df$`b_ex_10%` > 0)][1]
            end_time <- df$time[which(df$`b_ex_10%` > 0)][length(which(df$`b_ex_10%` > 0))] + 1
            move_time <- end_time - start_time
          }
          if((start_time - df$time[which(df$`b_ex_15%` > 0)][1]) > ex_period/start_sensitivity){
            start_time <- df$time[which(df$`b_ex_15%` > 0)][1]
            end_time <- df$time[which(df$`b_ex_15%` > 0)][length(which(df$`b_ex_15%` > 0))] + 1
            move_time <- end_time - start_time
          }
          if((start_time - df$time[which(df$`b_ex_20%` > 0)][1]) > ex_period/start_sensitivity){
            start_time <- df$time[which(df$`b_ex_20%` > 0)][1]
            end_time <- df$time[which(df$`b_ex_20%` > 0)][length(which(df$`b_ex_20%` > 0))] + 1
            move_time <- end_time - start_time
          }
          if((start_time - df$time[which(df$`b_ex_25%` > 0)][1]) > ex_period/start_sensitivity){
            start_time <- df$time[which(df$`b_ex_25%` > 0)][1]
            end_time <- df$time[which(df$`b_ex_25%` > 0)][length(which(df$`b_ex_25%` > 0))] + 1
            move_time <- end_time - start_time
          }

        }else if(length(which(df$`b_ex_10%` > 0)) > 0){
          start_time <- df$time[which(df$`b_ex_10%` > 0)][1]
          end_time <- df$time[which(df$`b_ex_10%` > 0)][length(which(df$`b_ex_10%` > 0))] + 1
          move_time <- end_time - start_time
          if((start_time - df$time[which(df$`b_ex_15%` > 0)][1]) > ex_period/start_sensitivity){
            start_time <- df$time[which(df$`b_ex_15%` > 0)][1]
            end_time <- df$time[which(df$`b_ex_15%` > 0)][length(which(df$`b_ex_15%` > 0))] + 1
            move_time <- end_time - start_time
          }
          if((start_time - df$time[which(df$`b_ex_20%` > 0)][1]) > ex_period/start_sensitivity){
            start_time <- df$time[which(df$`b_ex_20%` > 0)][1]
            end_time <- df$time[which(df$`b_ex_20%` > 0)][length(which(df$`b_ex_20%` > 0))] + 1
            move_time <- end_time - start_time
          }
          if((start_time - df$time[which(df$`b_ex_25%` > 0)][1]) > ex_period/start_sensitivity){
            start_time <- df$time[which(df$`b_ex_25%` > 0)][1]
            end_time <- df$time[which(df$`b_ex_25%` > 0)][length(which(df$`b_ex_25%` > 0))] + 1
            move_time <- end_time - start_time
          }

        }else if(length(which(df$`b_ex_15%` > 0)) > 0){
          start_time <- df$time[which(df$`b_ex_15%` > 0)][1]
          end_time <- df$time[which(df$`b_ex_15%` > 0)][length(which(df$`b_ex_15%` > 0))] + 1
          move_time <- end_time - start_time
          if((start_time - df$time[which(df$`b_ex_20%` > 0)][1]) > ex_period/start_sensitivity){
            start_time <- df$time[which(df$`b_ex_20%` > 0)][1]
            end_time <- df$time[which(df$`b_ex_20%` > 0)][length(which(df$`b_ex_20%` > 0))] + 1
            move_time <- end_time - start_time
          }
          if((start_time - df$time[which(df$`b_ex_25%` > 0)][1]) > ex_period/start_sensitivity){
            start_time <- df$time[which(df$`b_ex_25%` > 0)][1]
            end_time <- df$time[which(df$`b_ex_25%` > 0)][length(which(df$`b_ex_25%` > 0))] + 1
            move_time <- end_time - start_time
          }

        }else if(length(which(df$`b_ex_20%` > 0)) > 0){
          start_time <- df$time[which(df$`b_ex_20%` > 0)][1]
          end_time <- df$time[which(df$`b_ex_20%` > 0)][length(which(df$`b_ex_20%` > 0))] + 1
          move_time <- end_time - start_time
          if((start_time - df$time[which(df$`b_ex_25%` > 0)][1]) > ex_period/start_sensitivity){
            start_time <- df$time[which(df$`b_ex_25%` > 0)][1]
            end_time <- df$time[which(df$`b_ex_25%` > 0)][length(which(df$`b_ex_25%` > 0))] + 1
            move_time <- end_time - start_time
          }

        }else if(length(which(df$`b_ex_25%` > 0)) > 0){
          start_time <- df$time[which(df$`b_ex_25%` > 0)][1]
          end_time <- df$time[which(df$`b_ex_25%` > 0)][length(which(df$`b_ex_25%` > 0))] + 1
          move_time <- end_time - start_time

        }else{
          start_time <- Inf
          end_time <- Inf
          move_time <- Inf
        }

      }else{  # ex_sign == "negative"

        if(length(which(df$`b_ex_97.5%` < 0)) > 0){
          start_time <- df$time[which(df$`b_ex_97.5%` < 0)][1]
          end_time <- df$time[which(df$`b_ex_97.5%` < 0)][length(which(df$`b_ex_97.5%` < 0))] + 1
          move_time <- end_time - start_time
          if((start_time - df$time[which(df$`b_ex_95%` < 0)][1]) > ex_period/start_sensitivity){
            start_time <- df$time[which(df$`b_ex_95%` < 0)][1]
            end_time <- df$time[which(df$`b_ex_95%` < 0)][length(which(df$`b_ex_95%` < 0))] + 1
            move_time <- end_time - start_time
          }
          if((start_time - df$time[which(df$`b_ex_90%` < 0)][1]) > ex_period/start_sensitivity){
            start_time <- df$time[which(df$`b_ex_90%` < 0)][1]
            end_time <- df$time[which(df$`b_ex_90%` < 0)][length(which(df$`b_ex_90%` < 0))] + 1
            move_time <- end_time - start_time
          }
          if((start_time - df$time[which(df$`b_ex_85%` < 0)][1]) > ex_period/start_sensitivity){
            start_time <- df$time[which(df$`b_ex_85%` < 0)][1]
            end_time <- df$time[which(df$`b_ex_85%` < 0)][length(which(df$`b_ex_85%` < 0))] + 1
            move_time <- end_time - start_time
          }
          if((start_time - df$time[which(df$`b_ex_80%` < 0)][1]) > ex_period/start_sensitivity){
            start_time <- df$time[which(df$`b_ex_80%` < 0)][1]
            end_time <- df$time[which(df$`b_ex_80%` < 0)][length(which(df$`b_ex_80%` < 0))] + 1
            move_time <- end_time - start_time
          }
          if((start_time - df$time[which(df$`b_ex_75%` < 0)][1]) > ex_period/start_sensitivity){
            start_time <- df$time[which(df$`b_ex_75%` < 0)][1]
            end_time <- df$time[which(df$`b_ex_75%` < 0)][length(which(df$`b_ex_75%` < 0))] + 1
            move_time <- end_time - start_time
          }

        }else if(length(which(df$`b_ex_95%` < 0)) > 0){
          start_time <- df$time[which(df$`b_ex_95%` < 0)][1]
          end_time <- df$time[which(df$`b_ex_95%` < 0)][length(which(df$`b_ex_95%` < 0))] + 1
          move_time <- end_time - start_time
          if((start_time - df$time[which(df$`b_ex_90%` < 0)][1]) > ex_period/start_sensitivity){
            start_time <- df$time[which(df$`b_ex_90%` < 0)][1]
            end_time <- df$time[which(df$`b_ex_90%` < 0)][length(which(df$`b_ex_90%` < 0))] + 1
            move_time <- end_time - start_time
          }
          if((start_time - df$time[which(df$`b_ex_85%` < 0)][1]) > ex_period/start_sensitivity){
            start_time <- df$time[which(df$`b_ex_85%` < 0)][1]
            end_time <- df$time[which(df$`b_ex_85%` < 0)][length(which(df$`b_ex_85%` < 0))] + 1
            move_time <- end_time - start_time
          }
          if((start_time - df$time[which(df$`b_ex_80%` < 0)][1]) > ex_period/start_sensitivity){
            start_time <- df$time[which(df$`b_ex_80%` < 0)][1]
            end_time <- df$time[which(df$`b_ex_80%` < 0)][length(which(df$`b_ex_80%` < 0))] + 1
            move_time <- end_time - start_time
          }
          if((start_time - df$time[which(df$`b_ex_75%` < 0)][1]) > ex_period/start_sensitivity){
            start_time <- df$time[which(df$`b_ex_75%` < 0)][1]
            end_time <- df$time[which(df$`b_ex_75%` < 0)][length(which(df$`b_ex_75%` < 0))] + 1
            move_time <- end_time - start_time
          }

        }else if(length(which(df$`b_ex_90%` < 0)) > 0){
          start_time <- df$time[which(df$`b_ex_90%` < 0)][1]
          end_time <- df$time[which(df$`b_ex_90%` < 0)][length(which(df$`b_ex_90%` < 0))] + 1
          move_time <- end_time - start_time
          if((start_time - df$time[which(df$`b_ex_85%` < 0)][1]) > ex_period/start_sensitivity){
            start_time <- df$time[which(df$`b_ex_85%` < 0)][1]
            end_time <- df$time[which(df$`b_ex_85%` < 0)][length(which(df$`b_ex_85%` < 0))] + 1
            move_time <- end_time - start_time
          }
          if((start_time - df$time[which(df$`b_ex_80%` < 0)][1]) > ex_period/start_sensitivity){
            start_time <- df$time[which(df$`b_ex_80%` < 0)][1]
            end_time <- df$time[which(df$`b_ex_80%` < 0)][length(which(df$`b_ex_80%` < 0))] + 1
            move_time <- end_time - start_time
          }
          if((start_time - df$time[which(df$`b_ex_75%` < 0)][1]) > ex_period/start_sensitivity){
            start_time <- df$time[which(df$`b_ex_75%` < 0)][1]
            end_time <- df$time[which(df$`b_ex_75%` < 0)][length(which(df$`b_ex_75%` < 0))] + 1
            move_time <- end_time - start_time
          }

        }else if(length(which(df$`b_ex_85%` < 0)) > 0){
          start_time <- df$time[which(df$`b_ex_85%` < 0)][1]
          end_time <- df$time[which(df$`b_ex_85%` < 0)][length(which(df$`b_ex_85%` < 0))] + 1
          move_time <- end_time - start_time
          if((start_time - df$time[which(df$`b_ex_80%` < 0)][1]) > ex_period/start_sensitivity){
            start_time <- df$time[which(df$`b_ex_80%` < 0)][1]
            end_time <- df$time[which(df$`b_ex_80%` < 0)][length(which(df$`b_ex_80%` < 0))] + 1
            move_time <- end_time - start_time
          }
          if((start_time - df$time[which(df$`b_ex_75%` < 0)][1]) > ex_period/start_sensitivity){
            start_time <- df$time[which(df$`b_ex_75%` < 0)][1]
            end_time <- df$time[which(df$`b_ex_75%` < 0)][length(which(df$`b_ex_75%` < 0))] + 1
            move_time <- end_time - start_time
          }

        }else if(length(which(df$`b_ex_80%` < 0)) > 0){
          start_time <- df$time[which(df$`b_ex_80%` < 0)][1]
          end_time <- df$time[which(df$`b_ex_80%` < 0)][length(which(df$`b_ex_80%` < 0))] + 1
          move_time <- end_time - start_time
          if((start_time - df$time[which(df$`b_ex_75%` < 0)][1]) > ex_period/start_sensitivity){
            start_time <- df$time[which(df$`b_ex_75%` < 0)][1]
            end_time <- df$time[which(df$`b_ex_75%` < 0)][length(which(df$`b_ex_75%` < 0))] + 1
            move_time <- end_time - start_time
          }

        }else if(length(which(df$`b_ex_75%` < 0)) > 0){
          start_time <- df$time[which(df$`b_ex_75%` < 0)][1]
          end_time <- df$time[which(df$`b_ex_75%` < 0)][length(which(df$`b_ex_75%` < 0))] + 1
          move_time <- end_time - start_time

        }else{
          start_time <- Inf
          end_time <- Inf
          move_time <- Inf
        }

      }

      mv_time <- data.frame(start_time = start_time, end_time = end_time, move_time = move_time)
      df_mv <- rbind(df_mv, mv_time)


      ## Plotting

      # Visual
      if(!is.null(visual)){
        vis <- dplyr::filter(visual, cell == i & index == j)$time
        label_statistical <- "Statistical"
        label_visual <- "Visual"
      }else{
        vis <- NULL
        label_statistical <- NULL
        label_visual <- NULL
      }

      # Shade
      if(shade == T){
        alpha = 0.3
      }else{
        alpha = 0
      }

      # Start time
      if(start_line == T){
        col1 = "orange"
        col2 = "aquamarine3"
      }else{
        col1 = "transparent"
        col2 = "transparent"
      }

      # Theme
      if(theme_plot == "bw"){
        theme_plot2 <- theme_bw(base_size = ps)
      }else if(theme_plot == "light"){
        theme_plot2 <- theme_light(base_size = ps)
      }else if(theme_plot == "classic"){
        theme_plot2 <- theme_classic(base_size = ps)
      }else if(theme_plot == "gray"){
        theme_plot2 <- theme_gray(base_size = ps)
      }else if(theme_plot == "dark"){
        theme_plot2 <- theme_dark(base_size = ps)
      }else if(theme_plot == "test"){
        theme_plot2 <- theme_test(base_size = ps)
      }else if(theme_plot == "minimal"){
        theme_plot2 <- theme_minimal(base_size = ps)
      }else if(theme_plot == "void"){
        theme_plot2 <- theme_void(base_size = ps)
      }

      # label
      if(unit1=="meter"){
        label_y <- bquote(atop(paste("Distance of ", .(res_name)), paste("from ", .(ex_name), " ", (m))))
        label_alpha <- bquote(atop("Velocity of movement", (m/.(unit2))))
        label_beta <- bquote(atop(paste("Coefficient of ", .(ex_name)), (m/.(unit2))))
        label_random <- bquote(atop("Random fluctuation", (m/.(unit2))))
      }else if(unit1=="centimeter"){
        label_y <- bquote(atop(paste("Distance of ", .(res_name)), paste("from ", .(ex_name), " ", (cm))))
        label_alpha <- bquote(atop("Velocity of movement", (cm/.(unit2))))
        label_beta <- bquote(atop(paste("Coefficient of ", .(ex_name)), (cm/.(unit2))))
        label_random <- bquote(atop("Random fluctuation", (cm/.(unit2))))
      }else if(unit1=="millimeter"){
        label_y <- bquote(atop(paste("Distance of ", .(res_name)), paste("from ", .(ex_name), " ", (mm))))
        label_alpha <- bquote(atop("Velocity of movement", (mm/.(unit2))))
        label_beta <- bquote(atop(paste("Coefficient of ", .(ex_name)), (mm/.(unit2))))
        label_random <- bquote(atop("Random fluctuation", (mm/.(unit2))))
      }else if(unit1=="micrometer"){
        label_y <- bquote(atop(paste("Distance of ", .(res_name)), paste("from ", .(ex_name), " ", (mu*m))))
        label_alpha <- bquote(atop("Velocity of movement", (mu*m/.(unit2))))
        label_beta <- bquote(atop(paste("Coefficient of ", .(ex_name)), (mu*m/.(unit2))))
        label_random <- bquote(atop("Random fluctuation", (mu*m/.(unit2))))
      }else if(unit1=="nanometer"){
        label_y <- bquote(atop(paste("Distance of ", .(res_name)), paste("from ", .(ex_name), " ", (nm))))
        label_alpha <- bquote(atop("Velocity of movement", (nm/.(unit2))))
        label_beta <- bquote(atop(paste("Coefficient of ", .(ex_name)), (nm/.(unit2))))
        label_random <- bquote(atop("Random fluctuation", (nm/.(unit2))))
      }else{
        label_y <- bquote(atop(paste("Distance of ", .(res_name)), paste("from ", .(ex_name), " (",  .(unit1), ")")))
        label_alpha <- bquote(atop("Velocity of movement", (.(unit1)/.(unit2))))
        label_beta <- bquote(atop(paste("Coefficient of ", .(ex_name)), paste("(", .(unit1), " / ", .(unit2), ")")))
        label_random <- bquote(atop("Random fluctuation", paste("(", .(unit1), " / ", .(unit2), ")")))
      }

      # Title of the plots
      if(length(cell_list) == 1){
        titles <- paste(stringr::str_to_title(res_name), " ", j, sep="")
      }else{
        titles <- paste(stringr::str_to_title(df_name), " ", i, ", ", res_name, " ", j, sep="")
      }

      # X-axis min and max of shade
      shade_xmin <- min(cell_list[[i]]$time[cell_list[[i]]$ex == 0])
      shade_xmax <- max(cell_list[[i]]$time[cell_list[[i]]$ex == 0])
      zero_time <- (cell_list[[i]]$time[cell_list[[i]]$ex == 0])
      boundary1 <- zero_time[which(diff(zero_time) != 1)]
      boundary2 <- zero_time[which(diff(zero_time) != 1)+1]
      shade_xmin <- c(shade_xmin, boundary2)
      shade_xmax <- c(boundary1, shade_xmax)

      # Location of text
      text_x1 <- max(cell_list[[i]]$time) - (max(cell_list[[i]]$time) - min(cell_list[[i]]$time)) * 0.13
      text_x2 <- max(cell_list[[i]]$time) - (max(cell_list[[i]]$time) - min(cell_list[[i]]$time)) * 0.16

      # Distance
      ymax <- max(cell_list[[i]][,j+2])
      ymin <- min(cell_list[[i]][,j+2])
      yrange <- (ymax - ymin)
      yceiling <-  ymax + yrange * 0.05
      yfloor <- ymin - yrange * 0.05

      g_dist <- ggplot(data = cell_list[[i]]) +
        annotate("rect", xmin = shade_xmin, xmax = shade_xmax,
                 ymin = yfloor, ymax = yceiling, alpha = alpha, fill = "gray50") +
        geom_line(aes(x = time, y = cell_list[[i]][,j+2]), size=0.5) +
        geom_vline(xintercept = mv_time$start_time, linetype="solid", col = col1) +
        geom_vline(xintercept = vis, linetype="dashed", col = col2) +
        annotate("text", x=text_x1, y=yceiling-yrange*0.08, label=label_statistical, col=col1, size = ps/ggplot2::.pt) +
        annotate("text", x=text_x2, y=yceiling-yrange*0.2, label=label_visual, col=col2, size = ps/ggplot2::.pt) +
        scale_x_continuous(expand = c(0,0)) +
        scale_y_continuous(expand = c(0,0)) +
        theme_plot2 +
        theme(legend.position = "none",
              axis.title=element_text(size = ps),
              axis.title.x=element_blank(),
              axis.text = element_text(size = ps),
              plot.title = element_text(size = ps, face = "bold")) +
        labs(title = titles,
             y = label_y)

      # alpha
      ymax <- max(df$Y)
      ymin <- min(df$Y)
      yrange <- (ymax - ymin)
      yceiling <-  ymax + yrange * 0.05
      yfloor <- ymin - yrange * 0.05

      g_alpha <- ggplot(data = df, aes(x = time)) +
        annotate("rect", xmin = shade_xmin, xmax = shade_xmax,
                 ymin = yfloor, ymax = yceiling, alpha = alpha, fill = "gray50") +
        geom_ribbon(aes(ymin = `alpha_2.5%`, ymax = `alpha_97.5%`), alpha = 0.5) +
        geom_line(aes(y = `alpha_50%`), size = 0.5) +
        geom_point(aes(y = Y), alpha = 0.5, size=0.5) +
        geom_vline(xintercept = mv_time$start_time, linetype="solid", col = col1) +
        geom_vline(xintercept = vis, linetype="dashed", col = col2) +
        geom_hline(yintercept = 0, linetype="dashed") +
        scale_x_continuous(expand = c(0,0)) +
        scale_y_continuous(expand = c(0,0)) +
        theme_plot2 +
        theme(legend.position = "none",
              axis.title=element_text(size = ps),
              axis.title.x=element_blank(),
              axis.text = element_text(size = ps),
              plot.title = element_blank()) +
        labs(y = label_alpha)

      # beta_ex
      ymax <- max(df$`b_ex_97.5%`, na.rm = T)
      ymin <- min(df$`b_ex_2.5%`, na.rm = T)
      yrange <- (ymax - ymin)
      yceiling <-  ymax + yrange * 0.05
      yfloor <- ymin - yrange * 0.05

      g_b_ex <- ggplot(data = df, aes(x = time)) +
        annotate("rect", xmin = shade_xmin, xmax = shade_xmax,
                 ymin = yfloor, ymax = yceiling, alpha = alpha, fill = "gray50") +
        geom_ribbon(aes(ymin = `b_ex_2.5%`, ymax = `b_ex_97.5%`), alpha = 0.5) +
        geom_line(aes(y = `b_ex_50%`), size = 0.5) +
        geom_vline(xintercept = mv_time$start_time, linetype="solid", col = col1) +
        geom_vline(xintercept = vis, linetype="dashed", col = col2) +
        geom_hline(yintercept = 0, linetype="dashed") +
        scale_x_continuous(expand = c(0,0)) +
        scale_y_continuous(expand = c(0,0)) +
        theme_plot2 +
        theme(legend.position = "none",
              axis.title=element_text(size = ps),
              axis.title.x=element_blank(),
              axis.text = element_text(size = ps),
              plot.title = element_blank()) +
        labs(y = label_beta)

      # Random movement
      ymax <- max(df$`w_97.5%`)
      ymin <- min(df$`w_2.5%`)
      yrange <- (ymax - ymin)
      yceiling <-  ymax + yrange * 0.05
      yfloor <- ymin - yrange * 0.05

      g_w <- ggplot(data = df, aes(x = time)) +
        annotate("rect", xmin = shade_xmin, xmax = shade_xmax,
                 ymin = yfloor, ymax = yceiling, alpha = alpha, fill = "gray50") +
        geom_ribbon(aes(ymin = `w_2.5%`, ymax = `w_97.5%`), alpha = 0.5) +
        geom_line(aes(y = `w_50%`), size = 0.5) +
        geom_vline(xintercept = mv_time$start_time, linetype="solid", col = col1) +
        geom_vline(xintercept = vis, linetype="dashed", col = col2) +
        geom_hline(yintercept = 0, linetype="dashed") +
        scale_x_continuous(expand = c(0,0)) +
        scale_y_continuous(expand = c(0,0)) +
        theme_plot2 +
        theme(legend.position = "none",
              axis.title=element_text(size = ps),
              axis.text = element_text(size = ps),
              plot.title = element_blank()) +
        labs(x = paste("Time (", unit2, ")", sep=""),
             y = label_random)

      # Integrate plots
      g <- g_dist + g_alpha + g_b_ex + g_w +
        plot_layout(ncol = 1, heights = c(1, 1, 1, 1))
      suppressWarnings(
        ggsave(paste0(out, "/pdf/ssm_individual_", df_name, i, "_", res_name, j, ".pdf"),
               g, height = ps*20*4/4, width = ps*10*1.2, units = "mm")
      )


      ## Diagnosis of MCMC

      # Check of Rhat
      bayesplot::color_scheme_set("viridisC")
      bayesplot::bayesplot_theme_set(bayesplot::theme_default(base_size = ps+2, base_family = "sans"))
      g <- bayesplot::mcmc_rhat(bayesplot::rhat(fit))
      ggsave(paste0(out, "/diagnosis/ssm_individual_rhat_", df_name, i, "_", res_name, j, ".pdf"),
             g, height = ps*20, width = ps*20, units = "mm")
      max_rhat <- names(which.max(bayesplot::rhat(fit)))

      # Confirmation of convergence
      g <- bayesplot::mcmc_combo(
        fit$draws(),
        combo = c("dens_overlay", "trace"),
        widths = c(1, 1),
        pars = c("b_ex[1]", paste0("b_ex[", data_list$N_ex, "]"), "alpha[1]", paste0("alpha[", data_list$N, "]"), max_rhat),
        gg_theme = theme_classic(base_size = ps+2)
      )
      ggsave(paste0(out, "/diagnosis/ssm_individual_combo_", df_name, i, "_", res_name, j, ".pdf"),
             g, height = ps*20, width = ps*20, units = "mm")


      ## Remove temporary files
      file.remove(paste0(output_dir, "/", outcsv_name))


      ## Release memory
      rm(fit, tmp_csv, tmp_csv_w, tmp_csv_b_ex, tmp_csv_alpha,
         tmp_csv_s_w, tmp_csv_s_b_ex, tmp_csv_s_Y,
         df)
      gc(reset = T);gc(reset = T)
    }
  }


  ## Remove the temporary directory
  unlink(output_dir, recursive = T)


  ## Save movement time
  if(!is.null(visual)){
    df_visual_mv <- cbind(visual, df_mv)
    names(df_visual_mv)[1:3] <- c(df_name, res_name, "visual_start_time")
  }else{
    res <- NULL
    cell <- NULL
    for(i in 1:length(cell_list)){
      res <- c(res, 1:(ncol(cell_list[[i]]) - 2))
      cell <- c(cell, rep(i, length = (ncol(cell_list[[i]]) - 2)))
    }
    df_res_name <- data.frame(cell = cell, res = res)
    names(df_res_name) <- c(df_name, res_name)
    df_visual_mv <- cbind(df_res_name, df_mv)
  }

  data.table::fwrite(df_visual_mv, file = paste0(out, "/csv/ssm_individual_", res_name, "_mvtime.csv"))
}
