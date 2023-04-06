
### Quantile Function
quantile99 <- function(x){
  stats::quantile(x, probs = c(0.005, 0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.975, 0.995), names = TRUE)
}



#' Bayesian inference of the state-space model (common model)
#'
#' \code{ssm_common} performs the Bayesian inference of the parameters of the
#' state-space model to extract the common dynamics among the replicated input data.
#' It estimates the Bayesian credible intervals of the parameters and the start time of
#' the influence of an explanatory variable in the regression model. Before using
#' this function, CmdStan (\url{https://mc-stan.org/docs/cmdstan-guide/cmdstan-installation.html})
#' must be installed. Installation of CmdStan from GitHub would be easier for beginners.
#' Install Git and execute a few commands in the Terminal app on Mac or
#' the Git Bash app on Windows.
#'
#' @param cell_list (list of data frame) The input time-series data. First column, time (column name "time");
#' second column, an explanatory variable (0 or 1, column name "ex"); third to the last columns,
#' distances of cells or organelles from the explanatory variable (
#' any column names are accepted). The velocities of the third to the last columns are
#' calculated from the distance and time, and used as response variables in the
#' modelling. See the following \strong{Examples} for further details.
#' @param mvtime (data frame) The movement time estimated by [ssm_individual]
#' ("ssm_individual_mvtime.csv"). This parameter is optional but should be set if the start time of
#' movement is affected by the distance from the explanatory variable.
#'
#' When a data frame is given for this parameter,
#' the model assumes the virtual response variable the distance of which from the explanatory variable is zero.
#' The observed dynamics of all response variables in a data frame are then assumed to be produced from
#' the dynamics of the virtual response variable (the common dynamics) with modifications of
#' (1) the time lag of the start time depending on the distance from the explanatory variable, and
#' (2) system noise and observation error.
#'
#' When nothing is given for this parameter, the time lag of the start time is not assumed, and
#' the observed dynamics of all response variables in a data frame are assumed to be produced from
#' the common dynamics.
#' @param out (character string) The path of the output directory.
#' @param seed (positive integer) A seed for random number generation to pass to CmdStan.
#' @param warmup (positive integer) The number of warmup iterations of MCMC sampling after thinning.
#' The default is 1000.
#' @param sampling (positive integer) The number of post-warmup iterations of MCMC sampling after thinning.
#' The default is 1000.
#' @param thin (positive integer) Intervals of MCMC samples. This is useful to reduce the autocorrelation of
#' MCMC samples and improve the convergence. The default is 3.
#' @param df_name (character string) The name of the data frame. This is used for
#' file names and graph titles. The default is "cell".
#' @param res_name (character string) The name of the response variable. This is used
#' for file names and graph labels.
#' @param ex_name (character string) The name of the explanatory variable. This is
#' used for graph labels.
#' @param df_idx (integer vector) Indexes of the data frame. This should be set
#' only when you want to set the indexes manually. This is used for
#' file names and graph titles. The default is `NULL` and the indexes are automatically set.
#' @param unit1 (character string) The unit of the response variable. One of "meter",
#' "centimeter", "millimeter", "micrometer", "nanometer". If another character
#' string is given, it is used as it is. This is used for graph labels.
#' @param unit2 (character string) The unit of time. This is used for graph labels.
#' @param shade (logical) Whether to draw shade in graphs during the period without
#' the explanatory variable. The default is `TRUE`.
#' @param ps (positive integer) Font size of graphs specified in pt. The default is 7 pt.
#' Plot sizes are automatically adjusted according to the font size.
#' @param theme_plot (character string) A plot theme of the [ggplot2] package. One of "bw", "light",
#' "classic", "gray", "dark", "test", "minimal" and "void". The default is "bw".
#' @returns A directory named after the `out` parameter is created, which has three subdirectories.
#' * A subdirectory "csv" includes "ssm_common_cell `j` .csv" and
#' "ssm_common_cell `j` _sd.csv" with `j` the indexes of data frames.
#' "ssm_common_cell `j` .csv" contains the Bayesian credible
#' intervals of time-varying parameters, where "Y" is the observed velocity, "w" is
#' the white noise, "alpha" is the common velocity, "dist" is the common distance, "b_ex" is
#' the common time-varying coefficient of an explanatory variable, "alpha_each" is
#' the true state of each velocity, and "b_ex_each" is the true state of each
#' time-varying coefficient. "ssm_common_cell `j` _sd.csv"
#' contains the Bayesian credible intervals of non-time-varying parameters, where
#' "s_w" and "s_b_ex" are the white noise and system noise, respectively.
#' * A subdirectory "pdf" includes "ssm_common_cell `j` .pdf".
#' This is the visualised results of the model. The figure consists of four panels:
#' (1) observed distance of `res_name` from `ex_name`, (2) velocity of `res_name`,
#' (3) regression coefficient of `ex_name`, and (4) random fluctuations
#' in the velocity. In these panels,
#' solid lines and shaded regions are the medians and 95%
#' credible intervals of the Bayesian inference, respectively. When the `shade` parameter is `TRUE`,
#' the shaded and light regions represent the periods without and with the explanatory
#' variable, respectively.
#' * A subdirectory "diagnosis" includes "ssm_common_combo_cell `j` .pdf" and
#' "ssm_common_rhat_cell `j` .pdf". These are the visualised diagnoses of
#' MCMC sampling. "ssm_common_combo_cell `j` .pdf" shows the
#' density plots of the posterior distributions and the trace plots for "b_ex" and
#' "alpha" at the start and end of the time series and the parameter with the worst Rhat value.
#' "ssm_common_rhat_cell `j` .pdf" indicates the Rhat values of parameters.
#' These are drawn using the [bayesplot] package.
#' @examples
#' # For the first-time usage, install the cmdstanr package (>= 0.5.2)
#' # (https://mc-stan.org/cmdstanr/index.html) and CmdStan
#' # (https://mc-stan.org/docs/cmdstan-guide/cmdstan-installation.html)
#'
#'
#' ### Real data example of chloroplast accumulation responses to a blue microbeam ###
#'
#' # Load packages
#' library(cellssm)
#'
#' # Load data
#' data("cell1", "cell2", "cell3", "cell4", "chloroplast_mvtime")
#' cell_list <- list(cell1, cell2, cell3, cell4)
#'
#' # Check the format of the input data
#' cell_list
#' chloroplast_mvtime
#'
#' \dontrun{
#' # Execution of state-space modeling
#'
#'   # Set the path to which CmdStan was installed
#'   cmdstanr::set_cmdstan_path("~/cmdstan/")
#'
#'   # With the data frame of the movement time. This is recommended.
#'   ssm_common(cell_list = cell_list, mvtime = chloroplast_mvtime, out = "08_ssm_common",
#'              res_name = "chloroplast", ex_name = "microbeam",
#'              unit1 = "micrometer", unit2 = "min")
#'
#'   # Without the data frame of the movement time.
#'   ssm_common(cell_list = cell_list, out = "08_ssm_common",
#'              res_name = "chloroplast", ex_name = "microbeam",
#'              unit1 = "micrometer", unit2 = "min")
#' }
#'
#'
#'
#' ### Simulated data example of Paramecium escape responses from laser heating ###
#'
#' # Load package
#' library(cellssm)
#'
#' # Load data
#' data("Paramecium", "Paramecium_mvtime")
#' cell_list <- list(Paramecium)
#'
#' # Check the format of input data
#' cell_list
#' Paramecium_mvtime
#'
#' \dontrun{
#' # Execution of state-space modelling
#'
#'   # Set the path where CmdStan was installed
#'   cmdstanr::set_cmdstan_path("~/cmdstan/")
#'
#'   # With the data frame of the movement time. This is recommended
#'   ssm_common(cell_list = cell_list, mvtime = Paramecium_mvtime, out = "18_ssm_common",
#'              df_name = "experiment", res_name = "Paramecium", ex_name = "heat",
#'              unit1 = "millimeter", unit2 = "sec")
#'
#'   # Without the data frame of the movement time
#'   ssm_common(cell_list = cell_list, out = "18_ssm_common",
#'              df_name = "experiment", res_name = "Paramecium", ex_name = "heat",
#'             unit1 = "millimeter", unit2 = "sec")
#' }
#'
#' @export
#'
ssm_common <- function(cell_list, mvtime=NULL, out, seed = 123, warmup=1000, sampling=1000, thin=3,
                       df_name = "cell", res_name, ex_name, df_idx = NULL, unit1, unit2,
                       shade = TRUE, ps = 7, theme_plot = "bw"){

  ## Dependency on cmdstanr
  if(!requireNamespace("cmdstanr", quietly = TRUE)){
    stop("Package \"cmdstanr\" must be installed to use this function.",
         call. = FALSE)
  }


  ## Binding variables locally to the function
  time <- `dist_2.5%` <- `dist_97.5%` <- `dist_50%` <-
    `alpha_2.5%` <- `alpha_97.5%` <- `alpha_50%` <-
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

  # Adjust data.frame
  if(!is.null(mvtime)){
    if(ncol(mvtime) > 5){
      mvtime <- mvtime[,1:4]
      names(mvtime)[1:4] <- c("cell", "each", "visual", "predicted")
    }else{
      mvtime <- mvtime[,1:3]
      names(mvtime)[1:3] <- c("cell", "each", "predicted")
    }
    distance <- NULL
    for(i in 1:length(cell_list)){
      distance <- c(distance, as.numeric(cell_list[[i]][min(which(cell_list[[i]]$ex == 1))-1,-(1:2)]))
    }
    mvtime$distance <- distance

    # Remove infinity from cell_list and mvtime
    null_cell <- mvtime[is.infinite(rowSums(mvtime)),]$cell
    null_each <- mvtime[is.infinite(rowSums(mvtime)),]$each
    if(length(null_cell) > 0){
      for(i in 1:length(null_cell)){
        cell_list[[null_cell[i]]] <- cell_list[[null_cell[i]]][,-(null_each[i]+2)]
        mvtime <- mvtime[!is.infinite(rowSums(mvtime)),]
      }
    }
  }

  # Compile stan file
  if(is.null(mvtime)){
    stan_file <- system.file("extdata", "common_model_mvtimenull.stan", package = "cellssm")
  }else{
    stan_file <- system.file("extdata", "common_model.stan", package = "cellssm")
  }
  model <- cmdstanr::cmdstan_model(stan_file)


  ## Execution of the Bayesian inference

  # for loop of cells
  for(i in 1:length(cell_list)){

    # File name
    if(!is.null(df_idx)){
      file_name <- paste0(df_name, df_idx[i])
      start_time <- mvtime$predicted[mvtime$cell==df_idx[i]]
    }else{
      file_name <- paste0(df_name, i)
      start_time <- mvtime$predicted[mvtime$cell==i]
    }

    # Estimation of observation error
    sd_vel_all <- NULL
    for (j in 1:(ncol(cell_list[[i]])-2)){
      vel <- diff(cell_list[[i]][,j+2])
      sp_vel <- stats::smooth.spline(1:length(vel), vel)
      pred_vel <- stats::predict(sp_vel, 1:length(vel))
      sd_vel <- stats::sd(vel - pred_vel$y)
      sd_vel_all <- c(sd_vel_all, sd_vel)
    }
    obs <- stats::median(sd_vel_all)

    # Prepare data_list
    if(is.null(mvtime)){
      data_list <- list(
        N = nrow(cell_list[[i]])-1,
        N_ex = length(which(cell_list[[i]]$ex==1)),
        N_each = ncol(cell_list[[i]])-2,
        ex = cell_list[[i]]$ex[-1],
        Y = apply(cell_list[[i]][,-(1:2)], 2, diff),
        obs = obs
      )
    }else{
      data_list <- list(
        N = nrow(cell_list[[i]])-1,
        N_ex = length(which(cell_list[[i]]$ex==1)),
        N_each = ncol(cell_list[[i]])-2,
        ex = cell_list[[i]]$ex[-1],
        Y = apply(cell_list[[i]][,-(1:2)], 2, diff),
        obs = obs,
        start = start_time
      )
    }

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
      seed = seed,
      iter_warmup = warmup*thin,
      iter_sampling = sampling*thin,
      chains = 4,
      parallel_chains = 4,
      refresh = floor(warmup/2.5*thin),
      #show_messages = F,
      #sig_figs = 4,
      output_dir = output_dir,
      output_basename = file_name,
      adapt_delta = 0.95,
      thin = thin
    )

    # 99% Bayesian credible intervals
    outcsv_name <- list.files(output_dir)
    outcsv_name <- outcsv_name[grep(paste0(file_name, "-"), outcsv_name)]
    tmp_csv_w <- NULL
    tmp_csv_b_ex <- NULL
    tmp_csv_alpha <- NULL
    tmp_csv_b_ex_each <- NULL
    tmp_csv_alpha_each <- NULL
    tmp_csv_dist <- NULL
    tmp_csv_s_w <- NULL
    tmp_csv_s_b_ex <- NULL

    for(k in 1:length(outcsv_name)){
      tmp_csv <- as.data.frame(data.table::fread(cmd = paste0("grep -v '^#' ", output_dir, "/", outcsv_name[k])))
      tmp_csv_w <- rbind(tmp_csv_w, tmp_csv[,stringr::str_starts(names(tmp_csv), "w\\.")])
      tmp_csv_b_ex <- rbind(tmp_csv_b_ex, tmp_csv[,stringr::str_starts(names(tmp_csv), "b_ex\\.")])
      tmp_csv_alpha <- rbind(tmp_csv_alpha, tmp_csv[,stringr::str_starts(names(tmp_csv), "alpha\\.")])
      tmp_csv_b_ex_each <- rbind(tmp_csv_b_ex_each, tmp_csv[,stringr::str_starts(names(tmp_csv), "b_ex_each")])
      tmp_csv_alpha_each <- rbind(tmp_csv_alpha_each, tmp_csv[,stringr::str_starts(names(tmp_csv), "alpha_each")])
      tmp_csv_dist <- rbind(tmp_csv_dist, tmp_csv[,stringr::str_starts(names(tmp_csv), "dist")])
      tmp_csv_s_w <- c(tmp_csv_s_w, tmp_csv[,stringr::str_starts(names(tmp_csv), "s_w")])
      tmp_csv_s_b_ex <- c(tmp_csv_s_b_ex, tmp_csv[,stringr::str_starts(names(tmp_csv), "s_b_ex")])
    }

    df_w <- as.data.frame(t(apply(tmp_csv_w, 2, quantile99)))
    df_b_ex <- as.data.frame(t(apply(tmp_csv_b_ex, 2, quantile99)))
    df_alpha <- as.data.frame(t(apply(tmp_csv_alpha, 2, quantile99)))
    df_b_ex_each <- as.data.frame(t(apply(tmp_csv_b_ex_each, 2, quantile99)))
    df_alpha_each <- as.data.frame(t(apply(tmp_csv_alpha_each, 2, quantile99)))
    df_dist <- as.data.frame(t(apply(tmp_csv_dist, 2, quantile99)))
    df_s <- t(data.frame(s_w = quantile99(tmp_csv_s_w),
                         s_b_ex = quantile99(tmp_csv_s_b_ex)))
    df_s <- cbind(data.frame(s_name = row.names(df_s)), df_s)

    # Save output
    colnames(df_w) <- paste0("w_", colnames(df_w))
    colnames(df_b_ex) <- paste0("b_ex_", colnames(df_b_ex))
    colnames(df_alpha) <- paste0("alpha_", colnames(df_alpha))
    colnames(df_b_ex_each) <- paste0("b_ex_each_", colnames(df_b_ex_each))
    colnames(df_alpha_each) <- paste0("alpha_each_", colnames(df_alpha_each))
    colnames(df_dist) <- paste0("dist_", colnames(df_dist))
    df_b_ex <- as.data.frame(rbind(matrix(NA, nrow = data_list$boundary1, ncol = 21), as.matrix(df_b_ex), matrix(NA, nrow = data_list$N - data_list$boundary2 + 1, ncol = 21)))

    for(l in 1: data_list$N_each){
      eval(parse(text = paste0(
        "df_alpha_each", l, " <- df_alpha_each[stringr::str_ends(row.names(df_alpha_each), '\\\\.", l, "'),]
           df_b_ex_each", l, " <- df_b_ex_each[stringr::str_ends(row.names(df_b_ex_each), '\\\\.", l, "'),]
           colnames(df_alpha_each", l, ") <- paste0(colnames(df_alpha_each", l, "), '_", l, "')
           colnames(df_b_ex_each", l, ") <- paste0(colnames(df_b_ex_each", l, "), '_", l, "')
           df_b_ex_each", l, " <- as.data.frame(rbind(matrix(NA, nrow = data_list$boundary1, ncol = 21), as.matrix(df_b_ex_each", l, "), matrix(NA, nrow = data_list$N - data_list$boundary2 + 1, ncol = 21)))"
      )))
    }

    colnames(data_list$Y) <- stringr::str_replace(colnames(data_list$Y), "distance", "Y")
    dfs <- cbind(data.frame(time = cell_list[[i]]$time[-1]),
                 ex = cell_list[[i]]$ex[-1],
                 data.frame(data_list$Y),
                 df_w, df_alpha, df_dist, df_b_ex)

    for(l in 1: data_list$N_each){
      eval(parse(text = paste0(
        "dfs <- cbind(dfs, df_alpha_each", l, ")
           dfs <- cbind(dfs, df_b_ex_each", l, ")"
      )))
    }

    data.table::fwrite(dfs, file = paste0(out, "/csv/ssm_common_", file_name, ".csv"))
    data.table::fwrite(df_s, file = paste0(out, "/csv/ssm_common_", file_name, "_sd.csv"))


    ## Diagnosis of MCMC

    # Check of Rhat
    bayesplot::color_scheme_set("viridisC")
    bayesplot::bayesplot_theme_set(bayesplot::theme_default(base_size = ps+2, base_family = "sans"))
    suppressWarnings(g <- bayesplot::mcmc_rhat(bayesplot::rhat(fit)))
    suppressWarnings(ggsave(paste0(out, "/diagnosis/ssm_common_rhat_", file_name, ".pdf"),
                            g, height = ps*20, width = ps*20, units = "mm"))
    max_rhat <- names(which.max(bayesplot::rhat(fit)))

    # Confirmation of convergence
    g <- bayesplot::mcmc_combo(
      fit$draws(),
      combo = c("dens_overlay", "trace"),
      widths = c(1, 1),
      pars = c("b_ex[1]", paste0("b_ex[", data_list$N_ex, "]"), "alpha[1]", paste0("alpha[", data_list$N, "]"), max_rhat),
      gg_theme = theme_classic(base_size = ps+2)
    )
    ggsave(paste0(out, "/diagnosis/ssm_common_combo_", file_name, ".pdf"),
           g, height = ps*20, width = ps*20, units = "mm")


    # Remove temporary files
    file.remove(paste0(output_dir, "/", outcsv_name))

    rm(fit, tmp_csv, tmp_csv_w, tmp_csv_b_ex, tmp_csv_alpha, tmp_csv_dist,
       tmp_csv_b_ex_each, tmp_csv_alpha_each,
       tmp_csv_s_w, tmp_csv_s_b_ex,
       dfs)
  }


  ## Remove the temporary directory
  unlink(output_dir, recursive = T)


  ## Y range of plots
  ymax_dist_all <- NULL
  ymin_dist_all <- NULL
  ymax_alpha_all <- NULL
  ymin_alpha_all <- NULL
  ymax_b_ex_all <- NULL
  ymin_b_ex_all <- NULL
  ymax_w_all <- NULL
  ymin_w_all <- NULL

  # for loop of cells
  for(i in 1:length(cell_list)){

    # File name
    if(!is.null(df_idx)){
      file_name <- paste0(df_name, df_idx[i])
    }else{
      file_name <- paste0(df_name, i)
    }

    # Load data
    dfs <- data.table::fread(file = paste0(out, "/csv/ssm_common_", file_name, ".csv"))

    ymax_dist <- max(dfs$`dist_97.5%`)
    ymax_dist_all <- cbind(ymax_dist_all, ymax_dist)
    ymin_dist <- min(dfs$`dist_2.5%`)
    ymin_dist_all <- cbind(ymin_dist_all, ymin_dist)

    ymax_alpha <- max(dfs$`alpha_97.5%`)
    ymax_alpha_all <- cbind(ymax_alpha_all, ymax_alpha)
    ymin_alpha <- min(dfs$`alpha_2.5%`)
    ymin_alpha_all <- cbind(ymin_alpha_all, ymin_alpha)

    ymax_b_ex <- max(dfs$`b_ex_97.5%`, na.rm = T)
    ymax_b_ex_all <- cbind(ymax_b_ex_all, ymax_b_ex)
    ymin_b_ex <- min(dfs$`b_ex_2.5%`, na.rm = T)
    ymin_b_ex_all <- cbind(ymin_b_ex_all, ymin_b_ex)

    ymax_w <- max(dfs$`w_97.5%`)
    ymax_w_all <- cbind(ymax_w_all, ymax_w)
    ymin_w <- min(dfs$`w_2.5%`)
    ymin_w_all <- cbind(ymin_w_all, ymin_w)
  }


  ## Plotting

  # for loop of cells
  for(i in 1:length(cell_list)){

    # File name
    if(!is.null(df_idx)){
      file_name <- paste0(df_name, df_idx[i])
    }else{
      file_name <- paste0(df_name, i)
    }

    # Load data
    dfs <- data.table::fread(file = paste0(out, "/csv/ssm_common_", file_name, ".csv"))

    # Shade
    if(shade == T){
      alpha = 0.3
    }else{
      alpha = 0
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
    if(!is.null(df_idx)){
      titles <- paste(stringr::str_to_title(df_name), " ", df_idx[i], sep="")
    }else{
      titles <- paste(stringr::str_to_title(df_name), " ", i, sep="")
    }

    # X-axis min and max of shade
    shade_xmin <- min(cell_list[[i]]$time[cell_list[[i]]$ex == 0])
    shade_xmax <- max(cell_list[[i]]$time[cell_list[[i]]$ex == 0])
    zero_time <- (cell_list[[i]]$time[cell_list[[i]]$ex == 0])
    boundary1 <- zero_time[which(diff(zero_time) != 1)]
    boundary2 <- zero_time[which(diff(zero_time) != 1)+1]
    shade_xmin <- c(shade_xmin, boundary2)
    shade_xmax <- c(boundary1, shade_xmax)

    # distance
    ymax_dist <- max(ymax_dist_all)
    ymin_dist <- min(ymin_dist_all)
    yrange_dist <- (ymax_dist - ymin_dist)
    yceiling_dist <-  ymax_dist + yrange_dist * 0.05
    yfloor_dist <- ymin_dist - yrange_dist * 0.05

    g_dist <- ggplot(data = dfs, aes(x = time)) +
      annotate("rect", xmin = shade_xmin, xmax = shade_xmax,
               ymin = yfloor_dist, ymax = yceiling_dist, alpha = alpha, fill = "gray50") +
      geom_ribbon(aes(ymin = `dist_2.5%`, ymax = `dist_97.5%`), alpha = 0.5) +
      geom_line(aes(y = `dist_50%`), linewidth=0.5) +
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

    # Velocity
    ymax_alpha <- max(ymax_alpha_all)
    ymin_alpha <- min(ymin_alpha_all)
    yrange_alpha <- (ymax_alpha - ymin_alpha)
    yceiling_alpha <-  ymax_alpha + yrange_alpha * 0.05
    yfloor_alpha <- ymin_alpha - yrange_alpha * 0.05

    g_velocity <- ggplot(data = dfs, aes(x = time)) +
      annotate("rect", xmin = shade_xmin, xmax = shade_xmax,
               ymin = yfloor_alpha, ymax = yceiling_alpha, alpha = alpha, fill = "gray50") +
      geom_ribbon(aes(ymin = `alpha_2.5%`, ymax = `alpha_97.5%`), alpha = 0.5) +
      geom_line(aes(y = `alpha_50%`), linewidth = 0.5) +
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
    ymax_b <- max(ymax_b_ex_all)
    ymin_b <- min(ymin_b_ex_all)
    yrange_b <- (ymax_b - ymin_b)
    yceiling_b <-  ymax_b + yrange_b * 0.05
    yfloor_b <- ymin_b - yrange_b * 0.05

    g_b_ex <- ggplot(data = dfs, aes(x = time)) +
      annotate("rect", xmin = shade_xmin, xmax = shade_xmax,
               ymin = yfloor_b, ymax = yceiling_b, alpha = alpha, fill = "gray50") +
      geom_ribbon(aes(ymin = `b_ex_2.5%`, ymax = `b_ex_97.5%`), alpha = 0.5) +
      geom_line(aes(y = `b_ex_50%`), linewidth = 0.5) +
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
    ymax_w <- max(ymax_w_all)
    ymin_w <- min(ymin_w_all)
    yrange_w <- (ymax_w - ymin_w)
    yceiling_w <-  ymax_w + yrange_w * 0.05
    yfloor_w <- ymin_w - yrange_w * 0.05

    g_w <- ggplot(data = dfs, aes(x = time)) +
      annotate("rect", xmin = shade_xmin, xmax = shade_xmax,
               ymin = yfloor_w, ymax = yceiling_w, alpha = alpha, fill = "gray50") +
      geom_ribbon(aes(ymin = `w_2.5%`, ymax = `w_97.5%`), alpha = 0.5) +
      geom_line(aes(y = `w_50%`), linewidth = 0.5) +
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
    g <- g_dist + g_velocity + g_b_ex + g_w +
      plot_layout(ncol = 1, heights = c(1, 1, 1, 1))
    suppressWarnings(
      ggsave(paste0(out, "/pdf/ssm_common_", file_name, ".pdf"),
             g, height = ps*20*4/4, width = ps*10*1.2, units = "mm")
    )

  }

}
