# @title stack_densities
#
# Takes density lines of equal x range and stack them on top of each other symmetrically aournd zero.
#
# @param data a data frame
# @param bw bandwidth of kernal density
# @param n_grid number of x points that should be calculated. The higher the more smooth plot.
# @param center_fun a function that returns the y center for each possible x in range of x.
# @param method which estimation should be applied. Default is LOESS.
#
# @return a data frame
#
#
# @export
stack_densities <- function(data, bw = bw, n_grid = n_grid, center_fun = center_fun, method = method){

  if (is.null(center_fun)) {
    center_fun <- function(x) {
      return(0)
    }
  }

  fun <- switch(method,
                density = make_smooth_density,
                loess = make_smooth_loess,
                raw = make_connect_dots)

  x_range <- range(data$x)

  list <- lapply(split(data, data$group), fun, bw = bw, n_grid = n_grid, min_x = x_range[1], max_x = x_range[2])

  data <- do.call(rbind, list)

  data$group_tmp <- as.numeric(factor(data$group))

  data <- data[order(data$x, data$group_tmp),]

  data_list <- lapply(split(data, data$x), calc_y_offsets)

  data <- do.call(rbind, data_list)

  data$ymin <- data$ymin + center_fun(data$x)

  data$ymax <- data$ymax + center_fun(data$x)

  data <- lapply(split(data, data$group), extend_data)

  do.call(rbind, data)

}



compute_stacks <- function(df, method = 'themeRiver') {

  range_x <- range(df$x)

  full_values <- data.frame(x = seq(range_x[1], range_x[2], 1/10^max(sapply(df$x, decimal_places))))

  others <- df[!names(df) %in% c("x", "y")]

  sub_df <- df[names(df) %in% c("x", "y", "group")]

  list <- split(sub_df, df$group)

  list <- lapply(list, function(x) merge(full_values, x, by = "x", all.x = TRUE))

  list <- lapply(list,  replace_values)

  list <- lapply(list, function(x) setNames(x, c("x", unique(as.character(x$group)), "group")))

  list <- lapply(list, function(x) x[names(x) != "group"])

  values <- as.matrix(Reduce(function(x, y) merge(x, y, by = "x"), list))

  col_names <- colnames(values)

  xval <- values[, 1]

  values <- as.matrix(values[, -1])

  colnames(values) <- col_names[-1]

  dims <- dim(values)

  if (is.null(dims[1])) {

    timePoints <- length(values)
    nStreams <- 1
    values <- as.matrix(values)

  } else {

    timePoints <- dims[1]
    nStreams <- dim(values)[2]

  }

  out <- switch(method,
                base = base(values, timePoints, nStreams),
                themeRiver = themeRiver(values, timePoints, nStreams),
                newWiggle = newWiggle(values, timePoints, nStreams),
                minimizedWiggle = minimizedWiggle(values, timePoints, nStreams))

  out <- merge(out, unique(others))

  out$x <- rep(c(xval, rev(xval)), length(unique(out$group)))

  out


}
