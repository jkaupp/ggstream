
newWiggle <- function(values, timePoints, nStreams) {

  thin2large <- sort(apply(values, 2, FUN = var),
                     decreasing = FALSE, index.return = TRUE)$ix

  idxStreams <- c(thin2large[seq(1, length(thin2large), 2)],
                  thin2large[seq(length(thin2large), 2, -2)])

  yy <- matrix(0, timePoints, (nStreams * 2))

  for (iStream in 1 : nStreams) {

    tmpVals <- values[, idxStreams[iStream]]

    if (iStream > 1) {

      yy[, iStream * 2 - 1] <- yy[, (iStream - 1) * 2]

      yy[, iStream * 2] <- yy[, iStream * 2 - 1] + tmpVals

    } else {

      baseline <- rowSums(matrix((nStreams - 1 : nStreams - .5),
                                 nrow = timePoints, ncol = nStreams, byrow = TRUE) * values)

      yy[, 1] <- - (baseline / nStreams)

      yy[, 2] <- yy[, iStream * 2 - 1] + tmpVals

    }


  }

  return(yy)

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

  y <- vector()

  for (iStream in 1:nStreams)
  {
    y <- cbind(y, c(yy[, iStream * 2], rev(yy[, iStream * 2 - 1])))

  }

  y_groups <- stack(setNames(as.data.frame(y), groups))

  data.frame(x = rep(c(1:timePoints, timePoints:1), length(groups)),
         y = y_groups$values,
         fill = y_groups$ind)

}

minimizedWiggle <- function(values, timePoints, nStreams) {

  yy <- matrix(0, timePoints, (nStreams * 2))

  for (iStream in 1 : nStreams) {

    tmpVals <- values[, iStream]

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

      yy[, 1] <- - baseline

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


compute_stacks <- function(df, method = 'themeRiver') {

  range_x <- range(df$x)

  full_values <- data.frame(x = seq(range_x[1], range_x[2], 1/10^max(sapply(df$x, decimal_places))))

  others <- df[!names(df) %in% c("x", "y")]

  sub_df <- df[names(df) %in% c("x", "y", "fill")]

  list <- split(sub_df, df$fill)

  list <- lapply(list, function(x) merge(full_values, x, by = "x", all.x = TRUE))

  list <- lapply(list,  replace_values)

  list <- lapply(list, function(x) setNames(x, c("x", unique(as.character(x$fill)), "fill")))

  list <- lapply(list, function(x) x[names(x) != "fill"])

  values <- as.matrix(Reduce(function(x, y) merge(x, y, by = "x"), list))

  xval <- values[, 1]

  values <- values[, -1]

  timePoints <- dim(values)[1]

  nStreams <- dim(values)[2]

  out <- switch(method,
                themeRiver = themeRiver(values, timePoints, nStreams),
                newWiggle = newWiggle(values, timePoints, nStreams),
                minimizedWiggle = minimizedWiggle(values, timePoints, nStreams))

  out <- merge(out, unique(others))

  out$x <- rep(c(xval, rev(xval)), length(unique(out$group)))

  out <- out[names(out) %in% c("x", "y", "group")]

  out


}
