base <- function(values, timePoints, nStreams) {

  yy <- matrix(0, timePoints, (nStreams * 2))

  for (iStream in 1 : nStreams) {

    tmpVals <- predict(smooth.spline(values[, iStream]))$y

    tmpVals[tmpVals < 0] <- 0

    if (iStream > 1) {

      yy[, iStream * 2 - 1] <- yy[, (iStream - 1) * 2]

      yy[, iStream * 2] <- yy[, iStream * 2 - 1] + tmpVals

    } else {

        yy[, 2] <- tmpVals

    }
  }

  groups <- colnames(values)

  if (is.null(groups)) {

    groups <- "-1"

  }

  y <- vector()

  for (iStream in 1:nStreams)
  {
    y <- cbind(y, c(yy[, iStream * 2], rev(yy[, iStream * 2 - 1])))

  }

  y_groups <- stack(setNames(as.data.frame(y), groups))

  data.frame(x = rep(c(1:timePoints, timePoints:1), length(groups)),
             y = y_groups$values,
             group = y_groups$ind)

}


newWiggle <- function(values, timePoints, nStreams) {

  thin2large <- sort(apply(values, 2, FUN = var),
                     decreasing = FALSE, index.return = TRUE)$ix

  idxStreams <- c(thin2large[seq(1, length(thin2large), 2)],
                  thin2large[seq(2, length(thin2large), 2)])

  yy <- matrix(0, timePoints, (nStreams * 2))

  for (iStream in 1 : nStreams) {

    tmpVals <- predict(smooth.spline(values[, idxStreams[iStream]]))$y

    if (iStream > 1) {

      yy[, iStream * 2 - 1] <- yy[, (iStream - 1) * 2]

      yy[, iStream * 2] <- yy[, iStream * 2 - 1] + tmpVals

    } else {

      baseline <- rowSums(matrix((nStreams - 1 : nStreams - .5),
                                 nrow = timePoints, ncol = nStreams, byrow = TRUE) * values)

      yy[, 1] <- -predict(smooth.spline(baseline / nStreams))$y

      yy[, 2] <- yy[, iStream * 2 - 1] + tmpVals

    }

  }

  groups <- colnames(values)[idxStreams]

  if (is.null(groups)) {

    groups <- "-1"

  }

  y <- vector()

  for (iStream in 1:nStreams)
  {
    y <- cbind(y, c(yy[, iStream * 2], rev(yy[, iStream * 2 - 1])))

  }

  y_groups <- stack(setNames(as.data.frame(y), groups))

  data.frame(x = rep(c(1:timePoints, timePoints:1), length(groups)),
             y = y_groups$values,
             group = y_groups$ind)

}

themeRiver <- function(values, timePoints, nStreams) {

  yy <- matrix(0, timePoints, (nStreams * 2))

  for (iStream in 1 : nStreams) {

    tmpVals <- predict(smooth.spline(values[, iStream]))$y

    if (iStream > 1) {

      yy[, iStream * 2 - 1] <- yy[, (iStream - 1) * 2]

      yy[, iStream * 2] <- yy[, iStream * 2 - 1] + tmpVals

    } else {

      yy[, 1] <- -(1/2) * predict(smooth.spline(rowSums(values)))$y

      yy[, 2] <- yy[, iStream * 2 - 1] + tmpVals
    }
  }

  groups <- colnames(values)

  if (is.null(groups)) {

    groups <- "-1"

  }

  y <- vector()

  for (iStream in 1:nStreams)
  {
    y <- cbind(y, c(yy[, iStream * 2], rev(yy[, iStream * 2 - 1])))

  }

  y_groups <- stack(setNames(as.data.frame(y), groups))

  data.frame(x = rep(c(1:timePoints, timePoints:1), length(groups)),
         y = y_groups$values,
         group = y_groups$ind)

}

minimizedWiggle <- function(values, timePoints, nStreams) {

  yy <- matrix(0, timePoints, (nStreams * 2))

  for (iStream in 1 : nStreams) {

    tmpVals <- predict(smooth.spline(values[, iStream]))$y

    if (iStream > 1) {

      yy[, iStream * 2 - 1] <- yy[, (iStream - 1) * 2]

      yy[, iStream * 2] <- yy[, iStream * 2 - 1] + tmpVals

    } else {

      baseline <- array(0, timePoints)

      for (ipoint in 1 : timePoints) {

        for (jStream in 1 : nStreams) {

          baseline[ipoint] = baseline[ipoint] +
            + (nStreams - jStream - .5) * values[ipoint, jStream]

        }

        baseline[ipoint] = baseline[ipoint] / nStreams
      }

      yy[, 1] <- -predict(smooth.spline(baseline))$y

      yy[, 2] <- yy[, iStream * 2 - 1] + tmpVals

    }
  }

  groups <- colnames(values)

  y <- vector()

  for (iStream in 1:nStreams)
  {
    y <- cbind(y, c(yy[, iStream * 2], rev(yy[, iStream * 2 - 1])))

  }

  y_groups <- stack(setNames(as.data.frame(y), groups))

  data.frame(x = rep(c(1:timePoints, timePoints:1), length(groups)),
             y = y_groups$values,
             group = y_groups$ind)
}


# @title make_smooth_density
#
# Takes points and turns them into a density line.
#
# @param .df a data frame that must contain x and y
# @param bw bandwidth of kernal density
# @param n_grid number of x points that should be calculated. The higher the more smooth plot.
# @param min_x minimum x value of all groups
# @param max_x maximum x value of all groups
#
# @return a data frame
make_smooth_density <- function(df, bw = bw, n_grid = n_grid, min_x, max_x) {

  group <- df$group[[1]]

  group_min_x <- min(df$x, na.rm = T)

  group_max_x <- max(df$x, na.rm = T)

  range_dist <- max_x - min_x

  group_average_y <- mean(df$y, na.rm = TRUE)

  df <-  df[stats::complete.cases(df), ]

  bwidth <- bw

  w <- df$y / sum(df$y)

  m <- stats::density(df$x, weights = w, from = min_x - range_dist, to = max_x + range_dist, n = n_grid, bw = bwidth)

  df <- data.frame(x = m$x,
                   y = m$y)

  df <- df[(df$x <= max_x & df$x >= min_x) | df$y > 1/10000 * max(df$y), ]

  # Un-normalize density so that height matches true data relative size

  mulitplier <- abs(group_max_x - group_min_x) * group_average_y

  df$y <- df$y * mulitplier

  data.frame(x = df$x,
             y = df$y,
             group = group)
}


# @title make_connect_dots
#
# Returns n number of points from data that perfectly fits data.
#
# @param df a data frame that must contain x and y
# @param n_grid number of x points that should be calculated. The higher the more smooth plot.
# @param min_x minimum x value of all groups
# @param max_x maximum x value of all groups
#
# @return a data frame
#
# @export
make_connect_dots <- function(df, n_grid = n_grid, min_x, max_x, ...){

  new_y <- sapply(split(df, list(df$x, df$group)), function(x) mean(x$y))

  df <- unique(df[c("x", "group")])

  df <- df[order(df$x),]

  df$y <- new_y

  group <- df$group[[1]]

  df <-  df[stats::complete.cases(df), ]

  range_x <- seq(min_x, max_x, length.out = n_grid)

  steps <- df$x[-1]

  outside_local_range <- range_x[range_x < min(df$x) | range_x > max(df$x)]

  inside_local_range <- range_x[range_x >= min(df$x) & range_x <= max(df$x)]

  df <- data.frame(xlag = c(NA, utils::head(df$x, -1)),
                   ylag = c(NA, utils::head(df$y, -1)),
                   x = df$x,
                   y = df$y)

  df <-  df[stats::complete.cases(df), ]

  inrange <- apply(df, 1, build_range, inside_range = inside_local_range)

  inrange_out <- do.call(rbind, inrange)

  out <- rbind(inrange_out,
               data.frame(x = outside_local_range,
                          y = rep(0, length(outside_local_range))))

  out$group <- group

  out[order(out$x),]

}
