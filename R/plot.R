#' Plot partial covariate effects
#'
#' @param x A `statnnet` object.
#' @param variable One or more predictor names. The default plots all supported
#'   predictors.
#' @param uncertainty One of `"delta"`, `"simulation"`, or `"none"`.
#' @param level Confidence level.
#' @param ... Additional arguments passed to [pce()].
#'
#' @return The plotted PCE data, invisibly.
#' @export
plot.statnnet <- function(x, variable = NULL,
                          uncertainty = c("delta", "simulation", "none"),
                          level = 0.95, ...) {
  uncertainty <- match.arg(uncertainty)
  if (is.null(variable)) {
    variable <- x$variables
  }
  if (!is.character(variable) || !length(variable) ||
      any(!variable %in% x$variables)) {
    stop("`variable` must contain predictor names used in the fitted formula.", call. = FALSE)
  }

  plotted <- lapply(variable, function(name) {
    pce(
      x,
      variable = name,
      uncertainty = uncertainty,
      level = level,
      ...
    )
  })
  names(plotted) <- variable

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  if (length(plotted) > 1L) {
    columns <- ceiling(sqrt(length(plotted)))
    rows <- ceiling(length(plotted) / columns)
    graphics::par(mfrow = c(rows, columns))
  }

  for (name in variable) {
    values <- plotted[[name]]
    categorical <- all(is.na(values$value)) || any(!is.na(values$contrast))
    if (categorical) {
      positions <- seq_len(nrow(values))
      limits <- range(c(values$conf_low, values$conf_high, values$estimate), na.rm = TRUE)
      if (!all(is.finite(limits))) limits <- range(values$estimate, na.rm = TRUE)
      graphics::plot(
        positions,
        values$estimate,
        xaxt = "n",
        xlab = name,
        ylab = "Partial covariate effect",
        main = name,
        ylim = limits,
        pch = 19
      )
      graphics::axis(1, at = positions, labels = values$contrast)
      if (any(is.finite(values$conf_low))) {
        graphics::arrows(
          positions,
          values$conf_low,
          positions,
          values$conf_high,
          angle = 90,
          code = 3,
          length = 0.05
        )
      }
    } else {
      limits <- range(c(values$conf_low, values$conf_high, values$estimate), na.rm = TRUE)
      if (!all(is.finite(limits))) limits <- range(values$estimate, na.rm = TRUE)
      graphics::plot(
        values$value,
        values$estimate,
        type = "l",
        xlab = name,
        ylab = "Partial covariate effect",
        main = name,
        ylim = limits
      )
      if (any(is.finite(values$conf_low))) {
        graphics::lines(values$value, values$conf_low, lty = 2)
        graphics::lines(values$value, values$conf_high, lty = 2)
      }
    }
    graphics::abline(h = 0, lty = 3, col = "grey60")
  }
  invisible(if (length(plotted) == 1L) plotted[[1L]] else plotted)
}
