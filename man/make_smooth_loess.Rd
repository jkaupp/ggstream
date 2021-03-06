% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geom_stream.R
\name{make_smooth_loess}
\alias{make_smooth_loess}
\title{make_smooth_loess

Takes points and turns them into a LOESS-estimated line.}
\usage{
make_smooth_loess(.df, bw = bw, n_grid = n_grid, min_x, max_x)
}
\arguments{
\item{.df}{a data frame that must contain x and y}

\item{bw}{bandwidth of kernal density}

\item{n_grid}{number of x points that should be calculated. The higher the more smooth plot.}

\item{min_x}{minimum x value of all groups}

\item{max_x}{maximum x value of all groups}
}
\value{
a data frame
}
\description{
make_smooth_loess

Takes points and turns them into a LOESS-estimated line.
}
