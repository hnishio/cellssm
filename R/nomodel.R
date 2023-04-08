
#' Estimation of movement without the state-space model
#'
#' \code{nomodel} estimates the start time of the movement as the time
#' when the time-varying distances of cells or organelles from
#' the explanatory variable is increasing or decreasing both in the short term and
#' in the long term.
#'
#' @param cell_list (list of data frame) The input time-series data. First column, time (column name "time");
#' second column, an explanatory variable (0 or 1, column name "ex"); third to the last columns,
#' distances of cells or organelles from the explanatory variable (
#' any column names are accepted). See the following \strong{Examples} for further details.
#' @param visual (data frame) The optional data of visual estimation of the start time of
#' the influence of an explanatory variable. First column, cells (column name "cell");
#' second column, index (column name "index"); third column, the start time. The default is `NULL`.
#' @param out (character string) The path of the output directory.
#' @param ex_sign (character string) "positive" or "negative". This is used to
#' estimate the start time of the positive or negative influence of the explanatory
#' variable on the distances of cells or organelles.
#' If cells or organelles are moving away from the explanatory variable, `ex_sign`
#' should be "positive". If cells or organelles are approaching to the explanatory
#' variable, `ex_sign` should be "negative".
#' @param period,fold (positive integer, positive real) One of the definitions of the start time:
#' difference from the distance at `period` time points ahead is `fold` times larger
#' than the average change rate. The default is 5 for `period` and 2 for `fold`.
#' @param df_name (character string) The name of the data frame. This is used for
#' file names and graph titles. The default is "cell".
#' @param res_name (character string) The name of the response variable. This is used
#' for file names and graph labels. The default is "organelle".
#' @param ex_name (character string) The name of the explanatory variable. This is
#' used for graph labels.
#' @param df_idx (integer vector) Indexes of the data frame. This should be set
#' only when you want to set the indexes manually. This is used for
#' file names and graph titles. For example, setting "`df_idx` = c(1,3) and `res_idx`
#' = c(2,8)" with the default `df_name` and `res_name` results in the file names
#' including "cell1_organelle2" and "cell3_organelle8". The default is `NULL` and
#' the indexes are automatically set.
#' @param res_idx (integer vector) Indexes of the response variable. This should be set
#' only when you want to set the indexes manually. This is used for
#' file names and graph titles. For example, setting "`df_idx` = c(1,3) and `res_idx`
#' = c(2,8)" with the default `df_name` and `res_name` results in the file names
#' including "cell1_organelle2" and "cell3_organelle8". The default is `NULL` and
#' the indexes are automatically set.
#' @param unit1 (character string) The unit of the response variable. One of "meter",
#' "centimeter", "millimeter", "micrometer", "nanometer". If another character
#' string is given, it is used as it is. This is used for graph labels.
#' @param unit2 (character string) The unit of time. This is used for graph labels.
#' @param shade (logical) Whether to draw shade in graphs during the period without
#' the explanatory variable. The default is `TRUE`.
#' @param start_line (logical) Whether to draw a line at the start time of the influence of
#' the explanatory variable in graphs. The default is `TRUE`.
#' @param ps (positive integer) Font size of graphs specified in pt. The default is 7 pt.
#' Plot sizes are automatically adjusted according to the font size.
#' @param theme_plot (character string) A plot theme of the [ggplot2] package. One of "bw", "light",
#' "classic", "gray", "dark", "test", "minimal" and "void". The default is "bw".
#' @returns A directory named after the `out` parameter is created.
#'
#' "nomodel_cell `i` _ `res_name` `j` .pdf" (with `i` the indexes of cells,
#' `j` the indexes of the response variables) is the visualised result of the estimation.
#' The observed distance of `res_name` from `ex_name` is shown as a solid line.
#' When the optional visual estimation of the start time is given to the `visual`
#' parameter, orange solid lines and green dashed lines represent the start time
#' estimated by the model and the visual observation, respectively. When the `shade`
#' parameter is `TRUE`, the shaded and light regions represent the period without
#' and with the explanatory variable, respectively.
#'
#' "nomodel_mvtime.csv" contains the estimated start time, end time, and period of
#' the directional movement.
#' @examples
#' ### Real data example of chloroplast accumulation responses to a blue microbeam ###
#'
#' # Load package
#' library(cellssm)
#'
#' # Load data of chloroplast movements
#' data("cell1", "cell2", "cell3", "cell4", "visual")
#' cell_list <- list(cell1, cell2, cell3, cell4)
#'
#' # Check the format of input data
#' cell_list
#' visual
#'
#' \dontrun{
#' # Predict movement
#'
#' # When you do not want to compare the statistical and visual estimations of the start time
#' nomodel(cell_list = cell_list, out = "04_nomodel",
#'         res_name = "chloroplast", ex_name = "microbeam",
#'         unit1 = "micrometer", unit2 = "min")
#'
#' # When you do want to compare the statistical and visual estimations of the start time
#' nomodel(cell_list = cell_list, visual = visual, out = "04_nomodel",
#'         res_name = "chloroplast", ex_name = "microbeam",
#'         unit1 = "micrometer", unit2 = "min")
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
#' data("Paramecium")
#' cell_list <- list(Paramecium)
#'
#' # Check the format of input data
#' cell_list
#'
#' \dontrun{
#' # Predict movement
#' nomodel(cell_list = cell_list, out = "14_nomodel",
#'         ex_sign = "positive", fold = 1, df_name = "experiment",
#'         res_name = "Paramecium", ex_name = "heat",
#'         unit1 = "millimeter", unit2 = "sec")
#' }
#'
#' @export
#'
nomodel <- function(cell_list, visual = NULL, out,
                    ex_sign = "negative", period = 5, fold = 2,
                    df_name = "cell", res_name = "organelle", ex_name,
                    df_idx = NULL, res_idx = NULL, unit1, unit2,
                    shade = TRUE, start_line = TRUE, ps = 7, theme_plot = "bw"){


  ## Binding variables locally to the function
  index <- x <- y <- NULL


  ## Set up

  # Create an output directory
  if(file.exists(out)==F){
    dir.create(out, recursive=T)
  }

  # Prepare a container for movement time
  df_mv <- data.frame(NULL)


  ## Estimation of the start time

  # for loop of cells
  for(i in 1:length(cell_list)){

    # for loop of the response variable
    for(j in 1:(ncol(cell_list[[i]])-2)){

      # File name
      if(!is.null(df_idx) & !is.null(res_idx)){
        file_name <- paste0(df_name, df_idx[j], "_", res_name, res_idx[j])
      }else{
        file_name <- paste0(df_name, i, "_", res_name, j)
      }

      # Distance
      Y = cell_list[[i]][,j+2]

      if(ex_sign == "positive"){  # positive

        ## Estimate the start of movement
        ## definition1: approaching to light
        st_index <- min(which(cell_list[[i]]$ex == 1)) : max(which(cell_list[[i]]$ex == 1))
        st_index1 <- st_index[Y[st_index] < Y[st_index+1]] # increase at the next point
        st_index1 <- st_index1[-length(st_index1)]

        ## definition2: moving average of differences (sma_period points) are positive
        N_ex <- sum(cell_list[[i]]$ex == 1)
        diff <- diff(Y[st_index]) #difference between adjacent data
        sma <- NULL
        sma_period <- round(N_ex / 10, 0)
        for(k in 1:(length(diff)-(sma_period-1))){
          sma[k] <- mean(diff[k:(k+(sma_period-1))]) # sma of next sma_period differences
        }
        st_index2 <- st_index[which(sma > 0)]

        ## definition3: difference from the data after period tp is larger than fold*period/N_ex times of max - min
        st_index3 <- st_index[(Y[st_index+period] - Y[st_index]) >
                                (fold*period/N_ex)*(max(Y[st_index]) - min(Y[st_index]))]
        st_index3 <- st_index3[-((length(st_index3)-(period-1)):length(st_index3))]


        ## Overlap of def1, def2 and def3
        st_index1_2 <- intersect(st_index1, st_index2)
        st_index1_2_3 <- intersect(st_index1_2, st_index3)
        st_index_first <- min(st_index1_2_3, na.rm = T)
        st_index_last <- max(st_index1_2_3, na.rm = T)

        start_time <- cell_list[[i]]$time[st_index_first] #time at st_index_first
        end_time <- cell_list[[i]]$time[st_index_last] #time at st_index_last
        move_time <- end_time - start_time

        mv_time <- data.frame(start_time = start_time, end_time = end_time, move_time = move_time)
        df_mv <- rbind(df_mv, mv_time)

      }else{  # negative

        ## Estimate the start of movement
        ## definition1: approaching to light
        st_index <- min(which(cell_list[[i]]$ex == 1)) : max(which(cell_list[[i]]$ex == 1))
        st_index1 <- st_index[Y[st_index] > Y[st_index+1]] # decrease at the next point
        st_index1 <- st_index1[-length(st_index1)]

        ## definition2: moving average of differences (sma_period points) are negative
        N_ex <- sum(cell_list[[i]]$ex == 1)
        diff <- diff(Y[st_index]) #difference between adjacent data
        sma <- NULL
        sma_period <- round(N_ex / 10, 0)
        for(k in 1:(length(diff)-(sma_period-1))){
          sma[k] <- mean(diff[k:(k+(sma_period-1))]) # sma of next sma_period differences
        }
        st_index2 <- st_index[which(sma < 0)]

        ## definition3: difference from the data after period tp is larger than fold*period/N_ex times of max - min
        st_index3 <- st_index[(Y[st_index] - Y[st_index+period]) >
                                (fold*period/N_ex)*(max(Y[st_index]) - min(Y[st_index]))]
        st_index3 <- st_index3[-((length(st_index3)-(period-1)):length(st_index3))]

        ## Overlap of def1, def2 and def3
        st_index1_2 <- intersect(st_index1, st_index2)
        st_index1_2_3 <- intersect(st_index1_2, st_index3)
        st_index_first <- min(st_index1_2_3, na.rm = T)
        st_index_last <- max(st_index1_2_3, na.rm = T)

        start_time <- cell_list[[i]]$time[st_index_first] #time at st_index_first
        end_time <- cell_list[[i]]$time[st_index_last] #time at st_index_last
        move_time <- end_time - start_time

        mv_time <- data.frame(start_time = start_time, end_time = end_time, move_time = move_time)
        df_mv <- rbind(df_mv, mv_time)

      }

      ## Plotting

      # Visual
      if(!is.null(visual)){
        if(!is.null(df_idx) & !is.null(res_idx)){
          vis <- dplyr::filter(visual, cell == df_idx[j] & index == res_idx[j])$time
        }else{
          vis <- dplyr::filter(visual, cell == i & index == j)$time
        }
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

      # Title of the plots
      if(!is.null(df_idx) & !is.null(res_idx)){
        titles <- paste(stringr::str_to_title(df_name), " ", df_idx[j], ", ", res_name, " ", res_idx[j], sep="")
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
      text_x2 <- max(cell_list[[i]]$time) - (max(cell_list[[i]]$time) - min(cell_list[[i]]$time)) * 0.158

      # Distance
      df_g <- data.frame(x = cell_list[[i]]$time, y = Y)
      ymax <- max(df_g$y)
      ymin <- min(df_g$y)
      yrange <- (ymax - ymin)
      yceiling <-  ymax + yrange * 0.05
      yfloor <- ymin - yrange * 0.05

      g_dist <- ggplot(data = df_g) +
        annotate("rect", xmin = shade_xmin, xmax = shade_xmax,
                 ymin = yfloor, ymax = yceiling, alpha = alpha, fill = "gray50") +
        geom_line(aes(x = x, y = y), linewidth=0.5) +
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
             x = paste("Time (", unit2, ")", sep=""),
             y = label_y)

      # Integrate plots
      suppressWarnings(
        ggsave(paste0(out, "/nomodel_", file_name, ".pdf"),
               g_dist, height = ps*20*1/4, width = ps*10*1.2, units = "mm")
      )

    }
  }

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

  data.table::fwrite(df_visual_mv, file = paste0(out, "/nomodel_", res_name, "_mvtime.csv"))

}

