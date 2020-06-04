
computeStacks <- function(df, method = 'ThemeRiver'){

  list <- split(df, df$group)

  list <- lapply(list, function(x) setNames(x, c("x", unique(x[, 3]), "group")))

  list <- lapply(list, function(x) x[names(x) != "group"])

  values <- as.matrix(Reduce(function(x, y) merge(x, y, by = "x"), list))

  values <- values[, -1]

  timePoints <- dim(values)[1]

  nStreams <- dim(values)[2]

  if (method == "newWiggle"){

    thin2large <- sort(apply(values, 2, FUN=var),
                       decreasing = FALSE, index.return = TRUE)$ix

    idxStreams <- c(thin2large[seq(1, length(thin2large), 2)],
                    thin2large[seq(length(thin2large), 2, -2)])
  }

  yy <- matrix(0, timePoints, (nStreams * 2))

  for (iStream in 1 : nStreams) {

    tmpVals <- values[, iStream]

    if (method == "newWiggle")
      tmpVals <- values[, idxStreams[iStream]]

    if (iStream > 1){
      yy[, iStream * 2 - 1] <- yy[, (iStream - 1) * 2]
      yy[, iStream * 2] <- yy[, iStream * 2 - 1] + tmpVals
    } else {
      switch(method,
             ThemeRiver = {
               yy[, 1] <- -(1/2) * rowSums(values)
               yy[, 2] <- yy[, iStream * 2 - 1] + tmpVals},
             zero = {
               yy[, 2] <- tmpVals},
             minimizedWiggle = {
               baseline <- array(0, timePoints)
               for (ipoint in 1 : timePoints) {
                 for (jStream in 1 : nStreams) {
                   baseline[ipoint] = baseline[ipoint] +
                     + (nStreams - jStream - .5) * values[ipoint, jStream]}
                 baseline[ipoint] = baseline[ipoint] / nStreams}
               yy[, 1] <- - baseline
               yy[, 2] <- yy[, iStream * 2 - 1] + tmpVals},
             newWiggle = {
               baseline <- rowSums(matrix((nStreams - 1 : nStreams - .5),
                                          nrow = timePoints, ncol = nStreams, byrow = TRUE) * values)
               yy[, 1] <- - (baseline / nStreams)
               yy[, 2] <- yy[, iStream * 2 - 1] + tmpVals},
             { # default
               print(paste0(baseline, 'not recognized'))
               print('baseline value can be zero, ThemeRiver, minimizedWiggle or newWiggle')}
      )
    }# end: if (iStream > 1){
  }# end:
  return(yy)
}
