
#' Linear regression (x: distance at time 0, y: start time)
#'
#' \code{lm_dist_start} performs linear regression to analyse the
#' relationship between the distance at time 0 and the start time of the
#' influence of an explanatory variable obtained by the state-space model.
#'
#' @param cell_list (list of data frame) The input time-series data. First column, time (column name "time");
#' second column, an explanatory variable (0 or 1, column name "ex"); third to the last columns,
#' distances of cells or organelles from the explanatory variable (
#' any column names are accepted). See the following \strong{Examples} for further details.
#' @param mvtime (data frame) The movement time estimated by [ssm_individual]
#' ("ssm_individual_mvtime.csv") or [ssm_KFAS] ("ssm_KFAS_mvtime.csv").
#' @param robust (logical) When `TRUE`, performs robust regression using repeated medians
#' (Siegel regression). When `FALSE`, performs linear regression. The default is `FALSE`.
#' @param graph_title (character string) The name of the data frame. This is used for
#' graph titles. The default is "Cell".
#' @param ex_name (character string) The name of the explanatory variable. This is used for graph labels.
#' @param unit1 (character string) The unit of the response variable. One of "meter",
#' "centimeter", "millimeter", "micrometer", "nanometer". If another character
#' string is given, it is used as it is. This is used for graph labels.
#' @param unit2 (character string) The unit of time. This is used for graph labels.
#' @param eq_pos (character string) The position of equations. One of "topleft", "topright",
#' "bottomleft", "bottomright". The default is "topleft".
#' @param ps (positive integer) Font size of graphs specified in pt. The default is 7 pt.
#' Plot sizes are automatically adjusted according to the font size.
#' @param theme_plot (character string) A plot theme of the [ggplot2] package. One of "bw", "light",
#' "classic", "gray", "dark", "test", "minimal" and "void". The default is "bw".
#' @returns A list of ggplot objects is returned. In each plot, dots, solid lines
#' and shaded regions are the observed values, regression lines, and 95% confidence
#' intervals, respectively.
#' @examples
#' ### Real data example of chloroplast accumulation responses to a blue microbeam ###
#'
#' # Load package
#' library(cellssm)
#'
#' # Load data
#' data("cell1", "cell2", "cell3", "cell4", "chloroplast_mvtime")
#' cell_list <- list(cell1, cell2, cell3, cell4)
#'
#' # dist vs. start
#' glist <- lm_dist_start(cell_list = cell_list, mvtime = chloroplast_mvtime,
#'                        ex_name = "microbeam", unit1 = "micrometer", unit2 = "min")
#'
#' g <- glist[[1]] + glist[[2]] + glist[[3]] + glist[[4]] + glist[[5]] +
#'   patchwork::plot_layout(ncol = 3)
#'
#' \dontrun{
#' # Create an output directory
#' out <- "06_lm_dist_start"
#' if(file.exists(out)==FALSE){
#'   dir.create(out, recursive=TRUE)
#' }
#'
#' # Save output
#' suppressWarnings(ggplot2::ggsave(paste0(out, "/individual_chloroplast_lm_dist_start.pdf"),
#'                                  g, height = 110, width = 50*3, units = "mm"))
#' }
#'
#'
#'
#' ### Simulated data example of Paramecium escape responses from a laser heating ###
#'
#' # Load package
#' library(cellssm)
#'
#' # Load data of chloroplast movements
#' data("Paramecium", "Paramecium_mvtime")
#' cell_list <- list(Paramecium)
#'
#' # dist vs. start
#' glist <- lm_dist_start(cell_list = cell_list, mvtime = Paramecium_mvtime,
#'                        graph_title = "Experiment", ex_name = "heat",
#'                        unit1 = "millimeter", unit2 = "sec")
#'
#' \dontrun{
#' # Create an output directory
#' out <- "16_lm_dist_start"
#' if(file.exists(out)==FALSE){
#'   dir.create(out, recursive=TRUE)
#' }
#'
#' # Save output
#' g <- glist[[1]]
#' suppressWarnings(ggplot2::ggsave(paste0(out, "/individual_Paramecium_lm_dist_start.pdf"),
#'                                  g, height = 50, width = 50, units = "mm"))
#' }
#'
#' @export
#'
lm_dist_start <- function(cell_list, mvtime, robust = FALSE,
                          graph_title = "Cell", ex_name,
                          unit1, unit2, eq_pos = "topleft", ps = 7,
                          theme_plot = "bw"){


  ## Binding variables locally to the function
  fit <- lwr <- upr <- predicted <- ..rr.label.. <- NULL


  ## Error message
  if(length(cell_list) != length(unique(mvtime[,1]))){
    stop(paste("The length of 'cell_list' should be the same as the unique identifier at 'mvtime$cell' !!"))
  }


  ## Adjust data.frame
    if(ncol(mvtime) > 5){
      df <- mvtime[,1:4]
      names(df)[1:4] <- c("cell", "each", "visual", "predicted")
    }else{
      df <- mvtime[,1:3]
      names(df)[1:3] <- c("cell", "each", "predicted")
    }
    distance <- NULL
    for(i in 1:length(cell_list)){
      distance <- c(distance, as.numeric(cell_list[[i]][min(which(cell_list[[i]]$ex == 1))-1,-(1:2)]))
    }
    df$distance <- distance
    df <- df[!is.infinite(rowSums(df)),]
    df <- df[!is.na(rowSums(df)),]


  ## label
  if(unit1=="meter"){
    label_x <- bquote(paste("Distance from ", .(ex_name), " ", (m)))
    label_y <- bquote(paste("Start time ", (.(unit2))))
  }else if(unit1=="centimeter"){
    label_x <- bquote(paste("Distance from ", .(ex_name), " ", (cm)))
    label_y <- bquote(paste("Start time ", (.(unit2))))
  }else if(unit1=="millimeter"){
    label_x <- bquote(paste("Distance from ", .(ex_name), " ", (mm)))
    label_y <- bquote(paste("Start time ", (.(unit2))))
  }else if(unit1=="micrometer"){
    label_x <- bquote(paste("Distance from ", .(ex_name), " ", (mu*m)))
    label_y <- bquote(paste("Start time ", (.(unit2))))
  }else if(unit1=="nanometer"){
    label_x <- bquote(paste("Distance from ", .(ex_name), " ", (nm)))
    label_y <- bquote(paste("Start time ", (.(unit2))))
  }else{
    label_x <- bquote(paste("Distance from ", .(ex_name), " (", .(unit1), ")"))
    label_y <- bquote(paste("Start time ", (.(unit2))))
  }


  ## Theme
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


  ## Plotting for each cell
  glist <- list()
  for(i in 1:length(cell_list)){

    # Title of the plots
    if(length(cell_list) == 1){
      titles <- NULL
    }else{
      titles <- paste(graph_title, " ", i, sep="")
    }

    data <- df[df$cell==i,]

    range_x <- diff(range(df$distance))
    min_axis_x <- min(df$distance)
    max_axis_x <- max(df$distance)
    range_y <- diff(range(df$predicted))
    min_axis_y <- min(df$predicted) - range_y*0.1
    max_axis_y <- max(df$predicted) + range_y*0.1

    ##### Linear regression (x: distance, y: predicted) #####
    if(robust == T){
      model <- RobustLinearReg::siegel_regression(predicted ~ distance, data = data)
    }else{
      model <- stats::lm(predicted ~ distance, data = data)
    }
    suppressWarnings(conf_interval <- stats::predict(model, interval="confidence", level = 0.95))
    conf_interval2 <- as.data.frame(cbind(data$distance, conf_interval)[order(data$distance, decreasing = F),])
    names(conf_interval2)[1] <- "distance"

    r2 <- format(summary(model)$r.squared, digits=2, nsmall = 2)
    r2lab <- bquote(paste(italic(R^2), " = ", .(r2), sep=""))

    intercept <- format(summary(model)$coefficients[1,1], digits=2, nsmall = 2)
    coefficient <- format(summary(model)$coefficients[2,1], digits=2, nsmall = 2)
    equationlab <- bquote(paste(italic(y), " = ", .(intercept), " + ", .(coefficient), " ", italic(x), sep=""))
    if(as.numeric(coefficient) < 0){
      coefficient <- format(-summary(model)$coefficients[2,1], digits=2, nsmall = 2)
      equationlab <- bquote(paste(italic(y), " = ", .(intercept), " - ", .(coefficient), " ", italic(x), sep=""))
    }

    # Position of equations
    if(eq_pos == "topleft"){
      eq_x=min_axis_x+range_x*0.02; eq_y1=max_axis_y-range_y*0.05
      eq_x=min_axis_x+range_x*0.02; eq_y2=max_axis_y-range_y*0.17
    }else if(eq_pos == "topright"){
      eq_x=max_axis_x-range_x*0.63; eq_y1=max_axis_y-range_y*0.05
      eq_x=max_axis_x-range_x*0.63; eq_y2=max_axis_y-range_y*0.17
    }else if(eq_pos == "bottomleft"){
      eq_x=min_axis_x+range_x*0.02; eq_y1=min_axis_y+range_y*0.17
      eq_x=min_axis_x+range_x*0.02; eq_y2=min_axis_y+range_y*0.05
    }else if(eq_pos == "bottomright"){
      eq_x=max_axis_x-range_x*0.63; eq_y1=min_axis_y+range_y*0.17
      eq_x=max_axis_x-range_x*0.63; eq_y2=min_axis_y+range_y*0.05
    }

    glist[[i]] <- ggplot() +
      geom_line(data = conf_interval2, aes(x=distance, y=fit), color = "steelblue", linewidth = 1) +
      geom_ribbon(data = conf_interval2, aes(x=distance, ymin = lwr, ymax = upr), alpha = 0.4, fill = "steelblue") +
      geom_point(data = data, aes(x=distance, y=predicted), size=0.8, alpha=0.5) +
      annotate("text", x=eq_x, y=eq_y1,
               label=equationlab, size=ps/ggplot2::.pt, hjust = 0) +
      annotate("text", x=eq_x, y=eq_y2,
               label=r2lab, size=ps/ggplot2::.pt, hjust = 0) +
      coord_cartesian(xlim=c(min_axis_x, max_axis_x), ylim=c(min_axis_y, max_axis_y), clip='on') +
      theme_plot2 +
      theme(plot.title = element_text(size=ps, face = "bold"),
            axis.title=element_text(size=ps),
            axis.text=element_text(size=ps),
            plot.tag = element_text(size = 12, face = "bold"),
            plot.margin=unit(c(3,3,3,3), 'mm')) +
      labs(title = titles,
           x=label_x,
           y = label_y)
  }


  ## Plotting for all cells
  # Title of the plots
  if(length(cell_list) != 1){
    titles <- c(paste0("All ", graph_title, "s"))

    data <- df

    range_x <- diff(range(df$distance))
    min_axis_x <- min(df$distance)
    max_axis_x <- max(df$distance)
    range_y <- diff(range(df$predicted))
    min_axis_y <- min(df$predicted) - range_y*0.1
    max_axis_y <- max(df$predicted) + range_y*0.1

    ##### Linear regression (x: distance, y: predicted) #####
    if(robust == T){
      model <- RobustLinearReg::siegel_regression(predicted ~ distance, data = data)
    }else{
      model <- stats::lm(predicted ~ distance, data = data)
    }
    suppressWarnings(conf_interval <- stats::predict(model, interval="confidence", level = 0.95))
    conf_interval2 <- as.data.frame(cbind(data$distance, conf_interval)[order(data$distance, decreasing = F),])
    names(conf_interval2)[1] <- "distance"

    r2 <- format(summary(model)$r.squared, digits=2, nsmall = 2)
    r2lab <- bquote(paste(italic(R^2), " = ", .(r2), sep=""))

    intercept <- format(summary(model)$coefficients[1,1], digits=2, nsmall = 2)
    coefficient <- format(summary(model)$coefficients[2,1], digits=2, nsmall = 2)
    equationlab <- bquote(paste(italic(y), " = ", .(intercept), " + ", .(coefficient), " ", italic(x), sep=""))
    if(as.numeric(coefficient) < 0){
      coefficient <- format(-summary(model)$coefficients[2,1], digits=2, nsmall = 2)
      equationlab <- bquote(paste(italic(y), " = ", .(intercept), " - ", .(coefficient), " ", italic(x), sep=""))
    }

    glist[[length(cell_list) + 1]] <- ggplot() +
      geom_line(data = conf_interval2, aes(x=distance, y=fit), color = "steelblue", linewidth = 1) +
      geom_ribbon(data = conf_interval2, aes(x=distance, ymin = lwr, ymax = upr), alpha = 0.4, fill = "steelblue") +
      geom_point(data = data, aes(x=distance, y=predicted), size=0.8, alpha=0.5) +
      annotate("text", x=eq_x, y=eq_y1,
               label=equationlab, size=ps/ggplot2::.pt, hjust = 0) +
      annotate("text", x=eq_x, y=eq_y2,
               label=r2lab, size=ps/ggplot2::.pt, hjust = 0) +
      coord_cartesian(xlim=c(min_axis_x, max_axis_x), ylim=c(min_axis_y, max_axis_y), clip='on') +
      theme_plot2 +
      theme(plot.title = element_text(size=ps, face = "bold"),
            axis.title=element_text(size=ps),
            axis.text=element_text(size=ps),
            plot.tag = element_text(size = 12, face = "bold"),
            plot.margin=unit(c(3,3,3,3), 'mm')) +
      labs(title = titles,
           x=label_x,
           y = label_y)
  }

  return(glist)
}
