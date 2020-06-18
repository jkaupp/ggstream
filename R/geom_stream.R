StatStreamDensity <- ggplot2::ggproto(
  "StatStreamDensity",
  ggplot2::Stat,
  required_aes = c("x", "y"),
  extra_params = c("bw", "n_grid", "na.rm", "center_fun", "method"),

  setup_data = function(data, params) {

    .panels <- unique(data$PANEL)

    per_panel <- lapply(
          split(data, .panels),
          stack_streams,
          bw = params$bw,
          n_grid = params$n_grid,
          center_fun = params$center_fun,
          method = params$method)

    per_panel <- lapply(seq_along(per_panel), function(x) {
        data.frame(per_panel[[x]], PANEL = .panels[x])
      })

    per_panel <- do.call(rbind, per_panel)

    per_panel$PANEL <- factor(per_panel$PANEL)

    chars <- unique(data[!names(data) %in% c("x", "y")])

    chars$id  <- 1:nrow(chars)

    per_panel$p_id <- 1:nrow(per_panel)

    out <- merge(chars, per_panel, all.x = FALSE)

    out <- out[order(out$id, out$p_id), ]

    out <- out[!names(out) %in% c("id", "p_id")]

    rownames(out) <- NULL

    out

  },

  compute_group = function(data, scales) {
    data
  }
)


#' @title geom_stream
#'
#' stat to compute `geom_stream`
#'
#' @param mapping provide you own mapping. both x and y need to be numeric.
#' @param data provide you own data
#' @param geom change geom
#' @param position change position
#' @param na.rm remove missing values
#' @param show.legend show legend in plot
#' @param bw bandwidth of kernal density estimation
#' @param n_grid number of x points that should be calculated. The higher the more smooth plot.
#' @param center_fun a function that returns the y center for each possible x in range of x.
#' @param method Which method of estimation should be used. Default is LOESS, similar to `geom_smooth` but sets negative values to zero.
#' @param inherit.aes should the geom inherits aesthetics
#' @param ... other arguments to be passed to the geom
#'
#' @return a data frame
#' @export
geom_stream <- function(mapping = NULL, data = NULL, geom = "polygon",
                        position = "identity", show.legend = NA,
                        inherit.aes = TRUE, na.rm = T, bw = 0.75, n_grid = 3000, method = c("loess", "density", "raw", "base", "themeRiver", "newWiggle", "minimizedWiggle"), center_fun = NULL, ...) {
  method <- match.arg(method)
  ggplot2::layer(
    stat = StatStreamDensity, data = data, mapping = mapping, geom = geom,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, bw = bw, center_fun = center_fun, n_grid = n_grid, method = method, ...)
  )
}
