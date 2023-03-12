
#' Visualisation of the distance from an explanatory variable
#'
#' \code{dist_vis} visualises the distance of a response variable from an
#' explanatory variable
#'
#' @param cell_list (list of data frame) The input time-series data. First column, time (column name "time");
#' second column, an explanatory variable (0 or 1, column name "ex"); third to the last columns,
#' distances of cells or organelles from the explanatory variable (
#' any column names are accepted). See the following \strong{Examples} for further details.
#' @param df_name (character string) The name of the data frame. This is used for
#' graph titles. The default is "cell".
#' @param res_name (character string) The name of the response variable. This is used
#' for graph labels.
#' @param ex_name (character string) The name of the explanatory variable. This is used
#' for graph labels.
#' @param unit1 (character string) The unit of the response variable. One of "meter",
#' "centimeter", "millimeter", "micrometer", "nanometer". If another character
#' string is given, it is used as it is. This is used for graph labels.
#' @param unit2 (character string) The unit of time. This is used for graph labels.
#' @param col (character string) A colour map of [viridisLite::viridis] to draw the lines.
#' One of "viridis", "magma", "plasma", "inferno", "cividis", "mako", "rocket", and "turbo".
#' The default is "viridis".
#' @param shade (logical) Whether to draw shade in graphs during the period without
#' the explanatory variable. The default is `TRUE`.
#' @param ps (positive integer) Font size of graphs specified in pt. The default is 7 pt.
#' Plot sizes are automatically adjusted according to the font size.
#' @param theme_plot (character string) A plot theme of the [ggplot2] package. One of "bw", "light",
#' "classic", "gray", "dark", "test", "minimal" and "void". The default is "bw".
#' @returns A list of ggplot objects is returned. In each plot, coloured lines represent
#' different `res_name` in each data frame. When the `shade` parameter is `TRUE`,
#' the shaded and light regions represent the periods without and with the explanatory
#' variable, respectively.
#' @examples
#' ### Real data example of chloroplast accumulation responses to a blue microbeam ###
#'
#' # Load package
#' library(cellssm)
#'
#' # Load real data of chloroplast movements
#' data("cell1", "cell2", "cell3", "cell4")
#' cell_list <- list(cell1, cell2, cell3, cell4)
#'
#' # Check the format of input data
#' cell_list
#'
#' # Plotting
#' glist <- dist_vis(cell_list = cell_list,
#'                   res_name = "chloroplast", ex_name = "microbeam",
#'                   unit1 = "micrometer", unit2 = "min")
#'
#' g <- glist[[1]] + glist[[2]] + glist[[3]] + glist[[4]] +
#'  patchwork::plot_layout(ncol = 2, heights = c(1, 1))
#'
#' \dontrun{
#' # Create an output directory
#' out <- "01_dist_vis"
#' if(file.exists(out)==FALSE){
#'   dir.create(out, recursive=TRUE)
#' }
#'
#' # Save output
#' ggplot2::ggsave(paste0(out, "/distance_chloroplast.pdf"),
#'                 g, height = 120, width = 180, units = "mm")
#' }
#'
#'
#'
#' ### Simulated data example of Paramecium escape responses from a laser heating ###
#'
#' # Load package
#' library(cellssm)
#'
#' # Load simulated data of Paramecium movement
#' data("Paramecium")
#' cell_list <- list(Paramecium)
#'
#' # Check the format of input data
#' cell_list
#'
#' # Plotting
#' glist <- dist_vis(cell_list = cell_list,
#'                   res_name = "Paramecium", ex_name = "heat",
#'                   unit1 = "millimeter", unit2 = "sec")
#'
#' g <- glist[[1]]
#'
#' \dontrun{
#' # Create an output directory
#' out <- "11_dist_vis"
#' if(file.exists(out)==FALSE){
#'   dir.create(out, recursive=TRUE)
#' }
#'
#' # Save output
#' ggplot2::ggsave(paste0(out, "/distance_paramecium.pdf"),
#'                 g, height = 60, width = 90, units = "mm")
#' }
#'
#' @export
#'
dist_vis <- function(cell_list, df_name = "cell", res_name, ex_name, unit1, unit2,
                     col = "viridis", shade = TRUE, ps = 7, theme_plot = "bw"){

  # Binding variables locally to the function
  time <- NULL

  # Set axis limits
  xmaxs <- NULL
  xmins <- NULL
  ymaxs <- NULL
  ymins <- NULL
  for(i in 1:length(cell_list)){
    xmax <- max(cell_list[[i]]$time)
    xmin <- min(cell_list[[i]]$time)
    ymax <- max(cell_list[[i]][,3:ncol(cell_list[[i]])])
    ymin <- min(cell_list[[i]][,3:ncol(cell_list[[i]])])
    xmaxs <- cbind(xmaxs, xmax)
    xmins <- cbind(xmins, xmin)
    ymaxs <- cbind(ymaxs, ymax)
    ymins <- cbind(ymins, ymin)
  }

  xmax <- max(xmaxs)
  xmin <- min(xmins)
  xrange <- xmax - xmin
  ymax <- max(ymaxs)
  ymin <- min(ymins)
  yrange <- ymax - ymin

  max_yaxis <- ymax + yrange*0.05
  min_yaxis <- ymin - yrange*0.05
  max_xaxis <- xmax
  min_xaxis <- xmin


  ## Plotting

  # List of results
  glist <- list()

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
  }else if(unit1=="centimeter"){
    label_y <- bquote(atop(paste("Distance of ", .(res_name)), paste("from ", .(ex_name), " ", (cm))))
  }else if(unit1=="millimeter"){
    label_y <- bquote(atop(paste("Distance of ", .(res_name)), paste("from ", .(ex_name), " ", (mm))))
  }else if(unit1=="micrometer"){
    label_y <- bquote(atop(paste("Distance of ", .(res_name)), paste("from ", .(ex_name), " ", (mu*m))))
  }else if(unit1=="nanometer"){
    label_y <- bquote(atop(paste("Distance of ", .(res_name)), paste("from ", .(ex_name), " ", (nm))))
  }else{
    label_y <- bquote(atop(paste("Distance of ", .(res_name)), paste("from ", .(ex_name), " (",  .(unit1), ")")))
  }

  # for loop of cells
  for(i in 1:length(cell_list)){

    # Title of the plots
    if(length(cell_list) == 1){
      titles <- NULL
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

    g <- ggplot(data = cell_list[[i]], aes(x = time)) +
      annotate("rect", xmin = shade_xmin, xmax = shade_xmax,
               ymin = min_yaxis, ymax = max_yaxis, alpha = alpha, fill = "gray50") +
      coord_cartesian(xlim=c(min_xaxis, max_xaxis), ylim=c(min_yaxis, max_yaxis), clip='on') +
      scale_x_continuous(expand = c(0,0)) +
      scale_y_continuous(expand = c(0,0)) +
      theme_plot2 +
      theme(legend.position = "none",
            axis.title = element_text(size = ps),
            axis.text = element_text(size = ps),
            plot.title = element_text(size = ps, face = "bold"),
            plot.margin=unit(c(3,3,3,3), 'mm')) +
      labs(title = titles,
           x = paste("Time (", unit2, ")", sep=""),
           y = label_y)

    # for loop of response variable
    viridis_test <- viridisLite::plasma(10)

    for(j in 1:(ncol(cell_list[[i]])-2)){
      if(col == "viridis"){
        g <- g + eval(parse(text = paste("geom_line(aes(y = cell_list[[", i, "]][,", j+2, "]), col = as.vector(viridisLite::viridis(", (ncol(cell_list[[i]])-2), "))[", j, "])", sep="")))
      }else if(col == "magma"){
        g <- g + eval(parse(text = paste("geom_line(aes(y = cell_list[[", i, "]][,", j+2, "]), col = as.vector(viridisLite::magma(", (ncol(cell_list[[i]])-2), "))[", j, "])", sep="")))
      }else if(col == "plasma"){
        g <- g + eval(parse(text = paste("geom_line(aes(y = cell_list[[", i, "]][,", j+2, "]), col = as.vector(viridisLite::plasma(", (ncol(cell_list[[i]])-2), "))[", j, "])", sep="")))
      }else if(col == "inferno"){
        g <- g + eval(parse(text = paste("geom_line(aes(y = cell_list[[", i, "]][,", j+2, "]), col = as.vector(viridisLite::inferno(", (ncol(cell_list[[i]])-2), "))[", j, "])", sep="")))
      }else if(col == "cividis"){
        g <- g + eval(parse(text = paste("geom_line(aes(y = cell_list[[", i, "]][,", j+2, "]), col = as.vector(viridisLite::cividis(", (ncol(cell_list[[i]])-2), "))[", j, "])", sep="")))
      }else if(col == "mako"){
        g <- g + eval(parse(text = paste("geom_line(aes(y = cell_list[[", i, "]][,", j+2, "]), col = as.vector(viridisLite::mako(", (ncol(cell_list[[i]])-2), "))[", j, "])", sep="")))
      }else if(col == "rocket"){
        g <- g + eval(parse(text = paste("geom_line(aes(y = cell_list[[", i, "]][,", j+2, "]), col = as.vector(viridisLite::rocket(", (ncol(cell_list[[i]])-2), "))[", j, "])", sep="")))
      }else if(col == "turbo"){
        g <- g + eval(parse(text = paste("geom_line(aes(y = cell_list[[", i, "]][,", j+2, "]), col = as.vector(viridisLite::turbo(", (ncol(cell_list[[i]])-2), "))[", j, "])", sep="")))
      }
    }

    glist[[i]] <- g
  }

  return(glist)

}

